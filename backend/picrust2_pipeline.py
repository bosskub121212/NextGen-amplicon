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

    def _log(msg: str):
        with open(log_path, "a") as lf:
            lf.write(msg if msg.endswith("\n") else msg + "\n")

    def _step(cmd: list, timeout: int, label: str):
        """Run one PICRUSt2 script; log its output; return (ok, stderr tail)."""
        try:
            p = _run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            _log(f"WARNING: {label} timed out")
            return False, f"{label} timed out"
        _log(p.stdout + p.stderr)
        if p.returncode != 0:
            _log(f"WARNING: {label} failed (exit {p.returncode}): {p.stderr[-300:]}")
            return False, p.stderr[-500:]
        return True, ""

    def _first(*paths):
        for x in paths:
            if os.path.exists(x):
                return x
        return None

    # ── 2. hsp (hidden state prediction) ────────────────────────────
    # PICRUSt2 2.6 CLI: hsp.py -i <trait table> -t <placed tree> -o <out> -p N
    # (-i takes the trait name, -t the tree; -n adds the NSTI column, which
    # metagenome_pipeline needs for --max_nsti, so it goes on the marker run).
    _cp("hsp_marker", 30, "Predicting 16S copy numbers")

    marker_out   = os.path.join(output_dir, "marker_predicted_and_nsti.tsv.gz")
    func_ko_out  = os.path.join(output_dir, "KO_predicted.tsv.gz")
    func_ec_out  = os.path.join(output_dir, "EC_predicted.tsv.gz")
    func_cog_out = os.path.join(output_dir, "COG_predicted.tsv.gz")

    ok, err = _step([_picrust2_bin(env_path, "hsp.py"),
                     "-i", "16S", "-t", place_out, "-o", marker_out,
                     "-p", str(threads), "-n"], 1800, "hsp.py 16S")
    if not ok or not os.path.exists(marker_out):
        return {"success": False, "message": f"hsp.py 16S failed: {err}"}

    hsp_pct = {"KO": 38, "EC": 45, "COG": 50}
    func_tabs = {}
    for trait, out_f in [("EC", func_ec_out), ("KO", func_ko_out), ("COG", func_cog_out)]:
        _cp(f"hsp_{trait.lower()}", hsp_pct[trait], f"Predicting {trait} copy numbers")
        ok, _ = _step([_picrust2_bin(env_path, "hsp.py"),
                       "-i", trait, "-t", place_out, "-o", out_f,
                       "-p", str(threads)], 3600, f"hsp.py {trait}")
        if ok and os.path.exists(out_f):
            func_tabs[trait] = out_f
    if not func_tabs:
        return {"success": False, "message": "hsp.py failed for every function table (EC, KO, COG)"}

    # ── 3. metagenome_pipeline ──────────────────────────────────────
    # 2.6 CLI: metagenome_pipeline.py -i <table> -m <marker> -f <func> -o <dir>
    mg_dirs = {"EC": "EC_metagenome_out", "KO": "KO_metagenome_out", "COG": "COG_metagenome_out"}
    mg_pct = {"EC": 58, "KO": 64, "COG": 70}
    mg_done = {}
    first_err = ""
    for trait, func in func_tabs.items():
        _cp(f"metagenome_{trait.lower()}", mg_pct[trait], f"Predicting {trait} metagenome")
        outd = os.path.join(output_dir, mg_dirs[trait])
        ok, err = _step([_picrust2_bin(env_path, "metagenome_pipeline.py"),
                         "-i", asv_tsv, "-m", marker_out, "-f", func,
                         "-o", outd, "--max_nsti", "2.0"], 3600,
                        f"metagenome_pipeline {trait}")
        unstrat = _first(os.path.join(outd, "pred_metagenome_unstrat.tsv.gz"),
                         os.path.join(outd, "pred_metagenome_unstrat.tsv"))
        if ok and unstrat:
            mg_done[trait] = unstrat
        elif not first_err:
            first_err = err
    if not mg_done:
        return {"success": False, "message": f"metagenome_pipeline failed: {first_err}"}

    # ── 4. pathway_pipeline ─────────────────────────────────────────
    # MetaCyc pathways are inferred from the EC metagenome, not KO.
    _cp("pathway", 78, "Inferring MetaCyc pathway abundances")
    pw_out = os.path.join(output_dir, "pathways_out")
    pw_abund = None
    if "EC" in mg_done:
        ok, _ = _step([_picrust2_bin(env_path, "pathway_pipeline.py"),
                       "-i", mg_done["EC"], "-o", pw_out, "-p", str(threads)],
                      3600, "pathway_pipeline")
        pw_abund = _first(os.path.join(pw_out, "path_abun_unstrat.tsv.gz"),
                          os.path.join(pw_out, "path_abun_unstrat.tsv"))
    else:
        _log("WARNING: EC metagenome not available — skipping pathway_pipeline")

    # ── 5. add_descriptions ─────────────────────────────────────────
    _cp("descriptions", 88, "Adding KEGG / EC / MetaCyc / COG descriptions")
    desc_jobs = [(mg_done.get(t), t, os.path.join(output_dir, mg_dirs[t],
                  "pred_metagenome_unstrat_descrip.tsv.gz")) for t in ("KO", "EC", "COG")]
    desc_jobs.append((pw_abund, "METACYC",
                      os.path.join(pw_out, "path_abun_unstrat_descrip.tsv.gz")))
    for src_f, mapping, out_f in desc_jobs:
        if src_f:
            _step([_picrust2_bin(env_path, "add_descriptions.py"),
                   "-i", src_f, "-m", mapping, "-o", out_f], 600,
                  f"add_descriptions {mapping}")

    missing = [t for t in ("EC", "KO", "COG") if t not in mg_done]
    if pw_abund is None:
        missing.append("pathways")

    # ── Done ────────────────────────────────────────────────────────
    _cp("done", 100, "PICRUSt2 complete")

    output_files = []
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            rel = os.path.relpath(os.path.join(root, fname), output_dir)
            output_files.append(rel)

    summary = {
        "success":    True,
        "message":    "PICRUSt2 functional prediction complete"
                      + (f" (missing: {', '.join(missing)})" if missing else ""),
        "missing":    missing,
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
