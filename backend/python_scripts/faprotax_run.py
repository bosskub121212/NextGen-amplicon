#!/usr/bin/env python3
"""
NextGen-Amplicon — FAPROTAX Functional Annotation
Predicts functional groups from 16S taxonomy (bacteria/archaea only).

Usage:
  python3 faprotax_run.py \
    --otu_table  /path/to/job/r_tables/otu_relative.csv \
    --taxonomy   /path/to/job/exported/taxonomy/taxonomy.tsv \
    --output_dir /path/to/job/faprotax \
    --group_col  treatment

Requirements:
  pip install faprotax
  OR: download from http://www.loucalab.com/archive/FAPROTAX/
"""

import argparse, os, sys, subprocess, json
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="FAPROTAX functional annotation")
    parser.add_argument("--otu_table",  required=True, help="OTU/ASV table CSV (rows=taxa, cols=samples)")
    parser.add_argument("--taxonomy",   required=True, help="QIIME2 taxonomy TSV")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--group_col",  default="treatment", help="Group column in metadata")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"[FAPROTAX] Output: {args.output_dir}")

    # ── Load OTU table ──────────────────────────────────────────────────────
    otu_df = pd.read_csv(args.otu_table, index_col=0)
    print(f"  OTU table: {otu_df.shape[0]} taxa × {otu_df.shape[1]} samples")

    # ── Load taxonomy ───────────────────────────────────────────────────────
    tax_df = pd.read_csv(args.taxonomy, sep="\t", comment="#")
    tax_col = "Taxon" if "Taxon" in tax_df.columns else tax_df.columns[1]
    tax_df = tax_df.set_index(tax_df.columns[0])

    # ── Build FAPROTAX-compatible OTU table ─────────────────────────────────
    # Format: #OTU_ID, sample1, sample2, ..., taxonomy
    common = otu_df.index.intersection(tax_df.index)
    if len(common) == 0:
        print("[ERROR] No matching taxa between OTU table and taxonomy")
        sys.exit(1)

    combined = otu_df.loc[common].copy()
    combined["taxonomy"] = tax_df.loc[common, tax_col].values

    faprotax_input = os.path.join(args.output_dir, "otu_for_faprotax.tsv")
    combined.index.name = "#OTU ID"
    combined.to_csv(faprotax_input, sep="\t")
    print(f"  Written: {faprotax_input}")

    # ── Try running FAPROTAX ────────────────────────────────────────────────
    # Try installed faprotax first, then local download
    faprotax_script = None
    for candidate in [
        "collapse_table.py",                        # if installed via pip
        os.path.expanduser("~/faprotax/collapse_table.py"),
        "/usr/local/lib/faprotax/collapse_table.py"
    ]:
        if os.path.exists(candidate) or _cmd_exists(candidate.split()[0]):
            faprotax_script = candidate
            break

    if faprotax_script is None:
        print("[WARN] FAPROTAX collapse_table.py not found")
        print("  Install: pip install faprotax")
        print("  Or download from: http://www.loucalab.com/archive/FAPROTAX/")
        _write_instructions(args.output_dir, faprotax_input)
        return

    func_table = os.path.join(args.output_dir, "functional_table.tsv")
    log_file   = os.path.join(args.output_dir, "faprotax.log")

    cmd = [
        "python3", faprotax_script,
        "-i", faprotax_input,
        "-o", func_table,
        "--collapse_columns", "taxonomy",
        "--force"
    ]

    print(f"  Running FAPROTAX...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        with open(log_file, "w") as f:
            f.write(result.stdout + result.stderr)

        if result.returncode == 0 and os.path.exists(func_table):
            func_df = pd.read_csv(func_table, sep="\t", index_col=0)
            # Remove zero rows
            func_df = func_df.loc[func_df.sum(axis=1) > 0]
            func_df.to_csv(os.path.join(args.output_dir, "faprotax_results.csv"))

            summary = {
                "n_functions": len(func_df),
                "top_functions": func_df.sum(axis=1).nlargest(10).to_dict()
            }
            with open(os.path.join(args.output_dir, "faprotax_summary.json"), "w") as f:
                json.dump(summary, f, indent=2)

            print(f"  ✓ FAPROTAX complete: {len(func_df)} functional groups")
            print(f"  Top functions: {', '.join(list(summary['top_functions'].keys())[:5])}")
        else:
            print(f"[WARN] FAPROTAX failed — check {log_file}")
    except subprocess.TimeoutExpired:
        print("[WARN] FAPROTAX timed out (>5 min)")
    except Exception as e:
        print(f"[ERROR] FAPROTAX: {e}")


def _cmd_exists(cmd):
    import shutil
    return shutil.which(cmd) is not None


def _write_instructions(output_dir, input_file):
    instr = f"""FAPROTAX not installed. To run manually:

1. Install:
   pip install faprotax

2. Run:
   python3 -m faprotax.collapse_table \\
     -i {input_file} \\
     -o {output_dir}/functional_table.tsv \\
     --collapse_columns taxonomy --force

3. Results will be in: {output_dir}/functional_table.tsv
"""
    with open(os.path.join(output_dir, "faprotax_instructions.txt"), "w") as f:
        f.write(instr)
    print(f"  Instructions written to {output_dir}/faprotax_instructions.txt")


if __name__ == "__main__":
    main()
