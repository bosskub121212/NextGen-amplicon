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

_DNA = set("ACGTURYKMSWBDHVN")


def _is_seq(x: str) -> bool:
    x = (x or "").strip().upper()
    return len(x) >= 50 and set(x) <= _DNA


def _asv_csv_to_tsv(asv_csv: str, out_tsv: str, out_fasta: str | None = None) -> str | None:
    """
    Write a PICRUSt2 feature table (ASVs as rows, samples as columns, header
    "# OTU ID") and, when the sequences are available, a FASTA whose IDs are
    the SAME IDs. Returns the FASTA path written, or None.

    Two layouts arrive here:
      * DADA2 asv_table.csv: one row per ASV, one column per sample, plus a
        "sequence" column (and sometimes an unnamed index column).
      * samples x ASVs (COX1, older outputs): one row per sample, one column
        per ASV (IDs or sequences), plus an optional sample-ID column.

    The earlier version assumed the second layout for everything, so on DADA2
    output it wrote sample names as feature IDs; none matched the ASV1..ASV500
    IDs in asvs.fasta (which also covered only the top 500 ASVs) and PICRUSt2
    had nothing to predict from.
    """
    with open(asv_csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = [c for c in (reader.fieldnames or [])]
    if not rows:
        raise ValueError("ASV table CSV is empty")

    seq_col = next((c for c in fieldnames if c.lower() == "sequence"), None)
    if seq_col is None and all(_is_seq(r.get(fieldnames[-1], "")) for r in rows[:20]):
        seq_col = fieldnames[-1]
    id_col = next((c for c in fieldnames
                   if c in ("sample_id", "SampleID", "#SampleID", "sample", "")), None)

    feats: list[tuple[str, str | None, dict]] = []   # (id, sequence, {sample: count})
    if seq_col is not None:
        # ASVs as rows
        samples = [c for c in fieldnames if c not in (seq_col, id_col) and c != ""]
        for i, r in enumerate(rows, start=1):
            counts = {s_: (r.get(s_) or "0") for s_ in samples}
            feats.append((f"ASV{i}", (r.get(seq_col) or "").strip(), counts))
    else:
        # samples as rows
        asv_cols = [c for c in fieldnames if c != id_col and c != ""]
        samples = [(r.get(id_col) if id_col else None) or f"sample_{i}" for i, r in enumerate(rows)]
        for j, a in enumerate(asv_cols, start=1):
            counts = {samples[i]: (r.get(a) or "0") for i, r in enumerate(rows)}
            seq = a.strip() if _is_seq(a) else None
            feats.append((f"ASV{j}" if seq else a, seq, counts))

    # drop features with no reads anywhere
    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return 0.0
    feats = [f for f in feats if sum(_num(v) for v in f[2].values()) > 0]
    if not feats:
        raise ValueError("ASV table has no non-zero features")

    with open(out_tsv, "w") as f:
        f.write("# OTU ID\t" + "\t".join(samples) + "\n")
        for fid, _seq, counts in feats:
            f.write(fid + "\t" + "\t".join(str(counts.get(s_, "0")) for s_ in samples) + "\n")

    if out_fasta and all(seq for _id, seq, _c in feats):
        with open(out_fasta, "w") as f:
            for fid, seq, _c in feats:
                f.write(f">{fid}\n{seq}\n")
        return out_fasta
    return None


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

    # PICRUSt2's scripts call their helpers (hmmalign, epa-ng, gappa, castor via
    # Rscript) by name, so the env's bin/ must be on PATH for every child — the
    # backend's own PATH does not include it, which stopped place_seqs with
    # "No such file or directory: 'hmmalign'". Passed per call, not exported, so
    # the backend process itself is unchanged.
    _penv = dict(os.environ)
    _penv["PATH"] = str(Path(env_path) / "bin") + os.pathsep + _penv.get("PATH", "")
    _penv["CONDA_PREFIX"] = str(env_path)

    def _run(*a, **kw):
        kw.setdefault("env", _penv)
        return subprocess.run(*a, **kw)

    # Convert ASV table CSV → TSV
    _cp("convert_table", 8, "Converting ASV table")
    asv_tsv = os.path.join(output_dir, "asv_table.tsv")
    try:
        own_fasta = _asv_csv_to_tsv(asv_table_csv, asv_tsv,
                                    os.path.join(output_dir, "study_seqs.fna"))
    except Exception as e:
        return {"success": False, "message": f"ASV table conversion failed: {e}"}
    # The table's own sequences, under the table's own IDs, whenever the table
    # carries them — asvs.fasta holds only the top 500 ASVs under IDs that the
    # table does not use.
    if own_fasta:
        asv_fasta = own_fasta
    if not os.path.exists(asv_fasta):
        return {"success": False, "message": f"ASV FASTA not found: {asv_fasta}"}

    # ── 1. place_seqs ───────────────────────────────────────────────
    _cp("place_seqs", 15, "Placing ASV sequences onto reference tree")
    place_out = os.path.join(output_dir, "placed_seqs.tre")

    try:
        p = _run(
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
            p = _run(
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
        p = _run(
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
            p = _run(
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
            p = _run(
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
            p = _run(
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
            p = _run(
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
                _run(
                    [
                        _picrust2_bin(env_path, "add_descriptions.py"),
                        "-i", cog_pred,
                        "-m", "COG",
                        "-o", cog_desc_out,
                    ],
                    capture_output=True, text=True, timeout=300
                )
                cat_path = os.path.join(output_dir, "COG_categories.tsv")
                _run(
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
