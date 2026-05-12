"""
NextGen-Amplicon -- PICRUSt2 Functional Prediction Wrapper
Runs PICRUSt2 inside its conda environment against an ASV FASTA + abundance table.

Usage (called from main.py):
    from picrust2_pipeline import run_picrust2
    result = run_picrust2(job_id, asv_fasta, asv_table_csv, output_dir)
"""

import subprocess
import os
import sys
import json
import shutil
import tempfile
import csv
from pathlib import Path


# ── Locate PICRUSt2 conda environment ──────────────────────────────────────

def _find_picrust2_env() -> str | None:
    """Return path to the picrust2 conda env, or None if not found."""
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True, text=True, timeout=20
        )
        for line in result.stdout.splitlines():
            if "picrust2" in line.lower():
                parts = line.split()
                for p in parts:
                    if os.path.isdir(p):
                        return p
    except Exception:
        pass
    # Common fallback locations
    home = Path.home()
    for candidate in [
        home / "miniconda3" / "envs" / "picrust2",
        home / "anaconda3"  / "envs" / "picrust2",
        home / "conda"      / "envs" / "picrust2",
        Path("/opt/conda/envs/picrust2"),
        Path("/usr/local/envs/picrust2"),
    ]:
        if candidate.exists():
            return str(candidate)
    return None


def _picrust2_bin(env_path: str, cmd: str) -> str:
    """Return full path to a picrust2 binary inside the conda env."""
    return str(Path(env_path) / "bin" / cmd)


# ── Convert DADA2 ASV table CSV → BIOM-compatible TSV ──────────────────────

def _asv_csv_to_tsv(asv_csv: str, out_tsv: str):
    """
    Convert DADA2-style ASV table CSV (samples x ASVs, ASV_ID column or seq column)
    to a tab-separated OTU table with ASVs as rows and samples as columns.
    """
    with open(asv_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        raise ValueError("ASV table CSV is empty")

    fieldnames = reader.fieldnames or []

    # Find sample ID column
    id_col = None
    for name in ["sample_id", "SampleID", "#SampleID", "sample"]:
        if name in fieldnames:
            id_col = name
            break

    # ASV columns = everything except the id column
    asv_cols = [c for c in fieldnames if c != id_col and c != "sequence"]

    if not asv_cols:
        raise ValueError("No ASV columns found in ASV table")

    # Build a {asv_id: {sample: count}} dict
    samples = [r[id_col] if id_col else f"sample_{i}" for i, r in enumerate(rows)]
    asv_data: dict[str, dict] = {asv: {} for asv in asv_cols}

    for row in rows:
        samp = row[id_col] if id_col else f"sample_{rows.index(row)}"
        for asv in asv_cols:
            asv_data[asv][samp] = row.get(asv, "0") or "0"

    with open(out_tsv, "w") as f:
        f.write("# OTU ID\t" + "\t".join(samples) + "\n")
        for asv in asv_cols:
            counts = "\t".join(str(asv_data[asv].get(s, "0")) for s in samples)
            f.write(f"{asv}\t{counts}\n")


# ── Main entry point ────────────────────────────────────────────────────────

def run_picrust2(
    job_id: str,
    asv_fasta: str,
    asv_table_csv: str,
    output_dir: str,
    threads: int = 4,
    checkpoint_file: str | None = None,
) -> dict:
    """
    Run the full PICRUSt2 pipeline:
        1. place_seqs   — phylogenetic placement of ASVs
        2. hsp           — predict marker/function copy numbers
        3. metagenome_pipeline — predict metagenome
        4. pathway_pipeline   — infer MetaCyc pathway abundances

    Args:
        job_id:          Job identifier (used for logging)
        asv_fasta:       Path to ASV FASTA file
        asv_table_csv:   Path to ASV table CSV (DADA2 format: samples x ASVs)
        output_dir:      Directory to write PICRUSt2 outputs
        threads:         CPU threads
        checkpoint_file: Optional path to write JSON progress updates

    Returns:
        dict with keys: success (bool), message (str), output_dir (str), files (list)
    """
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "picrust2.log")

    def _cp(step: str, pct: int, msg: str = ""):
        info = {"step": step, "pct": pct, "msg": msg}
        if checkpoint_file:
            try:
                with open(checkpoint_file, "w") as cf:
                    json.dump(info, cf)
            except Exception:
                pass
        with open(log_path, "a") as lf:
            lf.write(f"[CHECKPOINT] {step} ({pct}%) {msg}\n")

    _cp("init", 2, "PICRUSt2 starting")

    # Verify input files
    if not os.path.exists(asv_fasta):
        return {"success": False, "message": f"ASV FASTA not found: {asv_fasta}"}
    if not os.path.exists(asv_table_csv):
        return {"success": False, "message": f"ASV table CSV not found: {asv_table_csv}"}

    # Find PICRUSt2 conda env
    _cp("find_env", 4, "Locating PICRUSt2 conda environment")
    env_path = _find_picrust2_env()
    if env_path is None:
        return {
            "success": False,
            "message": (
                "PICRUSt2 conda environment not found. "
                "Install with: conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y"
            )
        }

    with open(log_path, "a") as lf:
        lf.write(f"PICRUSt2 env: {env_path}\n")

    # Convert ASV table CSV → TSV
    _cp("convert_table", 8, "Converting ASV table")
    asv_tsv = os.path.join(output_dir, "asv_table.tsv")
    try:
        _asv_csv_to_tsv(asv_table_csv, asv_tsv)
    except Exception as e:
        return {"success": False, "message": f"ASV table conversion failed: {e}"}

    # ── 1. place_seqs ───────────────────────────────────────────────
    _cp("place_seqs", 15, "Placing ASV sequences onto reference tree")
    place_out = os.path.join(output_dir, "placed_seqs.tre")

    try:
        p = subprocess.run(
            [
                _picrust2_bin(env_path, "place_seqs.py"),
                "--study_fasta", asv_fasta,
                "--out_tree",    place_out,
                "--processes",   str(threads),
            ],
            capture_output=True, text=True, timeout=3600
        )
        with open(log_path, "a") as lf:
            lf.write(p.stdout + p.stderr)
        if p.returncode != 0:
            return {
                "success": False,
                "message": f"place_seqs failed (exit {p.returncode}): {p.stderr[-500:]}"
            }
    except subprocess.TimeoutExpired:
        return {"success": False, "message": "place_seqs timed out (>1 hour)"}
    except FileNotFoundError:
        return {
            "success": False,
            "message": f"place_seqs.py not found in env {env_path}. Is PICRUSt2 installed?"
        }

    # ── 2. hsp (hidden state prediction) ────────────────────────────
    _cp("hsp_marker", 30, "Predicting 16S copy numbers (EC + KO)")

    marker_out  = os.path.join(output_dir, "hsp_marker.tsv")
    func_ko_out = os.path.join(output_dir, "hsp_ko.tsv")
    func_ec_out = os.path.join(output_dir, "hsp_ec.tsv")
    func_cog_out= os.path.join(output_dir, "hsp_cog.tsv")

    for trait, out_f in [("16S", marker_out), ("KO", func_ko_out), ("EC", func_ec_out), ("COG", func_cog_out)]:
        try:
            p = subprocess.run(
                [
                    _picrust2_bin(env_path, "hsp.py"),
                    "-i", place_out,
                    "-t", trait,
                    "-o", out_f,
                    "-p", str(threads),
                    "--calculate_NSTI",
                ],
                capture_output=True, text=True, timeout=1800
            )
            with open(log_path, "a") as lf:
                lf.write(p.stdout + p.stderr)
            if p.returncode != 0:
                with open(log_path, "a") as lf:
                    lf.write(f"WARNING: hsp.py {trait} failed: {p.stderr[-300:]}\n")
        except subprocess.TimeoutExpired:
            with open(log_path, "a") as lf:
                lf.write(f"WARNING: hsp.py {trait} timed out\n")

    # ── 3. metagenome_pipeline ──────────────────────────────────────
    _cp("metagenome", 55, "Predicting metagenome functional abundances")

    mg_out = os.path.join(output_dir, "metagenome")
    os.makedirs(mg_out, exist_ok=True)

    try:
        p = subprocess.run(
            [
                _picrust2_bin(env_path, "metagenome_pipeline.py"),
                "-s",            asv_tsv,
                "-f",            func_ko_out,
                "-m",            marker_out,
                "--out_dir",     mg_out,
                "--max_nsti",    "2.0",
            ],
            capture_output=True, text=True, timeout=3600
        )
        with open(log_path, "a") as lf:
            lf.write(p.stdout + p.stderr)
        if p.returncode != 0:
            return {
                "success": False,
                "message": f"metagenome_pipeline failed (exit {p.returncode}): {p.stderr[-500:]}"
            }
    except subprocess.TimeoutExpired:
        return {"success": False, "message": "metagenome_pipeline timed out"}

    # ── 4. pathway_pipeline ─────────────────────────────────────────
    _cp("pathway", 75, "Inferring MetaCyc pathway abundances")

    pw_out = os.path.join(output_dir, "pathways")
    os.makedirs(pw_out, exist_ok=True)

    ko_metagenome = os.path.join(mg_out, "pred_metagenome_unstrat.tsv.gz")
    if not os.path.exists(ko_metagenome):
        ko_metagenome = os.path.join(mg_out, "pred_metagenome_unstrat.tsv")

    if os.path.exists(ko_metagenome):
        try:
            p = subprocess.run(
                [
                    _picrust2_bin(env_path, "pathway_pipeline.py"),
                    "-i", ko_metagenome,
                    "-o", pw_out,
                    "-p", str(threads),
                ],
                capture_output=True, text=True, timeout=3600
            )
            with open(log_path, "a") as lf:
                lf.write(p.stdout + p.stderr)
            if p.returncode != 0:
                with open(log_path, "a") as lf:
                    lf.write(f"WARNING: pathway_pipeline failed: {p.stderr[-300:]}\n")
        except subprocess.TimeoutExpired:
            with open(log_path, "a") as lf:
                lf.write("WARNING: pathway_pipeline timed out\n")
    else:
        with open(log_path, "a") as lf:
            lf.write("WARNING: KO metagenome output not found, skipping pathway_pipeline\n")

    # ── 5. add_descriptions ─────────────────────────────────────────
    _cp("descriptions", 88, "Adding KEGG/MetaCyc descriptions")

    pred_file = os.path.join(mg_out, "pred_metagenome_unstrat.tsv.gz")
    if not os.path.exists(pred_file):
        pred_file = os.path.join(mg_out, "pred_metagenome_unstrat.tsv")

    if os.path.exists(pred_file):
        desc_out = os.path.join(output_dir, "KO_predicted_with_descriptions.tsv")
        try:
            p = subprocess.run(
                [
                    _picrust2_bin(env_path, "add_descriptions.py"),
                    "-i", pred_file,
                    "-m", "KO",
                    "-o", desc_out,
                ],
                capture_output=True, text=True, timeout=300
            )
            with open(log_path, "a") as lf:
                lf.write(p.stdout + p.stderr)
        except Exception as e:
            with open(log_path, "a") as lf:
                lf.write(f"add_descriptions skipped: {e}\n")

    pw_abund = os.path.join(pw_out, "path_abun_unstrat.tsv.gz")
    if not os.path.exists(pw_abund):
        pw_abund = os.path.join(pw_out, "path_abun_unstrat.tsv")
    if os.path.exists(pw_abund):
        pw_desc_out = os.path.join(output_dir, "pathways_with_descriptions.tsv")
        try:
            p = subprocess.run(
                [
                    _picrust2_bin(env_path, "add_descriptions.py"),
                    "-i", pw_abund,
                    "-m", "METACYC",
                    "-o", pw_desc_out,
                ],
                capture_output=True, text=True, timeout=300
            )
            with open(log_path, "a") as lf:
                lf.write(p.stdout + p.stderr)
        except Exception as e:
            with open(log_path, "a") as lf:
                lf.write(f"pathway add_descriptions skipped: {e}\n")

    # ── COG metagenome + descriptions ──────────────────────────────
    _cp("cog", 92, "Predicting COG functional categories")

    if os.path.exists(func_cog_out):
        cog_mg_out = os.path.join(output_dir, "COG")
        os.makedirs(cog_mg_out, exist_ok=True)
        try:
            p = subprocess.run(
                [
                    _picrust2_bin(env_path, "metagenome_pipeline.py"),
                    "-s",        asv_tsv,
                    "-f",        func_cog_out,
                    "-m",        marker_out,
                    "--out_dir", cog_mg_out,
                    "--max_nsti", "2.0",
                ],
                capture_output=True, text=True, timeout=3600
            )
            with open(log_path, "a") as lf:
                lf.write(p.stdout + p.stderr)

            # add COG descriptions
            cog_pred = os.path.join(cog_mg_out, "pred_metagenome_unstrat.tsv.gz")
            if not os.path.exists(cog_pred):
                cog_pred = os.path.join(cog_mg_out, "pred_metagenome_unstrat.tsv")
            if os.path.exists(cog_pred):
                cog_desc_out = os.path.join(output_dir, "COG_predicted_with_descriptions.tsv")
                subprocess.run(
                    [
                        _picrust2_bin(env_path, "add_descriptions.py"),
                        "-i", cog_pred,
                        "-m", "COG",
                        "-o", cog_desc_out,
                    ],
                    capture_output=True, text=True, timeout=300
                )
                cat_path = os.path.join(output_dir, "COG_categories.tsv")
                subprocess.run(
                    [
                        _picrust2_bin(env_path, "add_descriptions.py"),
                        "-i", cog_pred,
                        "-m", "COG_category",
                        "-o", cat_path,
                    ],
                    capture_output=True, text=True, timeout=300
                )
                with open(log_path, "a") as lf:
                    lf.write("COG predictions + descriptions written.\n")
                cat("  ✓ COG functional categories predicted\n") if False else None
        except Exception as e:
            with open(log_path, "a") as lf:
                lf.write(f"COG metagenome skipped: {e}\n")
    else:
        with open(log_path, "a") as lf:
            lf.write("WARNING: COG hsp output not found — skipping COG metagenome\n")

    # ── Done ────────────────────────────────────────────────────────
    _cp("done", 100, "PICRUSt2 complete")

    output_files = []
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            rel = os.path.relpath(os.path.join(root, fname), output_dir)
            output_files.append(rel)

    summary = {
        "success":    True,
        "message":    "PICRUSt2 functional prediction complete",
        "output_dir": output_dir,
        "files":      output_files,
    }
    with open(os.path.join(output_dir, "picrust2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    return summary


# ── CLI entry point for standalone use ─────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="PICRUSt2 pipeline wrapper")
    parser.add_argument("--job_id",       required=True)
    parser.add_argument("--asv_fasta",    required=True)
    parser.add_argument("--asv_table",    required=True, help="ASV table CSV (DADA2 format)")
    parser.add_argument("--output_dir",   required=True)
    parser.add_argument("--threads",      type=int, default=4)
    parser.add_argument("--checkpoint",   default=None, help="Path for checkpoint JSON")
    args = parser.parse_args()

    result = run_picrust2(
        job_id          = args.job_id,
        asv_fasta       = args.asv_fasta,
        asv_table_csv   = args.asv_table,
        output_dir      = args.output_dir,
        threads         = args.threads,
        checkpoint_file = args.checkpoint,
    )
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["success"] else 1)
