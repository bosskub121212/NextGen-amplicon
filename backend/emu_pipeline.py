#!/usr/bin/env python3
"""
NextGen-Amplicon — ONT 16S pipeline via Emu
Handles: V7-V8 sub-region and V1-V9 full-length ONT reads

Usage (called by main.py):
    python emu_pipeline.py --input <fastq_dir> --output <out_dir>
                           --db_path <emu_db_dir>
                           [--primer_f SEQ --primer_r SEQ]
                           [--min_abundance 0.0001]
                           [--threads 4]
                           [--metadata metadata.csv]
                           [--topN 30]
                           [--region V7-V8|V1-V9]
                           [--job_name NAME]
"""

import argparse, csv, json, os, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

# ── Progress reporter (mirrors dada2_pipeline.R format for main.py SSE) ───────
def progress(pct: int, label: str):
    print(f"PROGRESS:{pct}|{label}", flush=True)

def log(msg: str):
    print(msg, flush=True)

# ── Argument parsing ────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",         required=True)
    p.add_argument("--output",        required=True)
    p.add_argument("--db_path",       required=True,  help="Path to Emu database directory")
    p.add_argument("--primer_f",      default="",     help="Forward primer sequence")
    p.add_argument("--primer_r",      default="",     help="Reverse primer sequence")
    p.add_argument("--min_abundance", type=float, default=0.0001)
    p.add_argument("--threads",       type=int,   default=4)
    p.add_argument("--metadata",      default="")
    p.add_argument("--topN",          type=int,   default=30)
    p.add_argument("--region",        default="V1-V9")
    p.add_argument("--job_name",      default="")
    p.add_argument("--marker",        default="ONT-16S")
    return p.parse_args()

# ── Helpers ─────────────────────────────────────────────────────────────────────
def find_fastq(input_dir: Path) -> dict[str, Path]:
    """Find FASTQ files and map sample_name → path."""
    samples = {}
    for ext in ("*.fastq.gz", "*.fastq", "*.fq.gz", "*.fq"):
        for f in sorted(input_dir.glob(ext)):
            name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", f.name)
            samples[name] = f
    return samples

def count_reads(fastq: Path) -> int:
    """Count reads in a fastq(.gz) file."""
    try:
        if str(fastq).endswith(".gz"):
            result = subprocess.run(["zcat", str(fastq)], capture_output=True)
            lines = result.stdout.count(b"\n")
        else:
            lines = fastq.read_bytes().count(b"\n")
        return lines // 4
    except Exception:
        return 0

def run_cmd(cmd: list, desc: str = "") -> bool:
    """Run a subprocess, stream stdout/stderr, return success."""
    if desc:
        log(f"  → {desc}")
    log(f"  CMD: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        log(f"  [WARN] Command returned {result.returncode}: {desc}")
        return False
    return True

# ── Step 1: Trim primers with cutadapt ─────────────────────────────────────────
def trim_primers(samples: dict, out_dir: Path, primer_f: str, primer_r: str,
                 threads: int) -> dict[str, Path]:
    """Trim primers from each sample. Returns trimmed sample map."""
    if not (primer_f or primer_r):
        log("  Primers not specified — skipping cutadapt trimming")
        return samples

    cutadapt = shutil.which("cutadapt") or shutil.which(
        str(Path.home() / ".local" / "bin" / "cutadapt"))
    if not cutadapt:
        log("  [WARN] cutadapt not found — skipping primer trimming")
        return samples

    trim_dir = out_dir / "trimmed"
    trim_dir.mkdir(exist_ok=True)
    trimmed = {}

    for name, fq in samples.items():
        out_fq = trim_dir / fq.name
        cmd = [cutadapt, "-j", str(threads)]
        if primer_f:
            cmd += ["-g", primer_f]
        if primer_r:
            cmd += ["-a", f"rc_{primer_r}", "-G", primer_r, "-A", f"rc_{primer_f}"]
        cmd += ["--discard-untrimmed", "-m", "50",
                "-o", str(out_fq), str(fq)]
        run_cmd(cmd, f"Trim {name}")
        trimmed[name] = out_fq if out_fq.exists() else fq

    return trimmed

# ── Step 2: Run Emu on each sample ─────────────────────────────────────────────
def run_emu(samples: dict, out_dir: Path, db_path: str,
            min_abundance: float, threads: int) -> dict[str, Path]:
    """Run emu abundance on each sample. Returns map of sample → TSV output."""
    # Build the emu command prefix.
    # Prefer 'conda run -n emu' to avoid venv/PATH conflicts when backend runs
    # in a Python virtualenv — emu uses #!/usr/bin/env python3 which would
    # otherwise pick up the venv Python instead of the conda emu env.
    conda = shutil.which("conda")
    emu_conda_bin = None
    for candidate in [
        Path.home() / "miniconda3" / "envs" / "emu" / "bin" / "emu",
        Path.home() / "anaconda3" / "envs" / "emu" / "bin" / "emu",
    ]:
        if candidate.exists():
            emu_conda_bin = str(candidate)
            break

    if emu_conda_bin:
        # Prefer calling the emu env's Python directly — avoids issues where
        # 'conda run --no-capture-output' doesn't fully activate the env and
        # fails to find packages like pysam even when they're installed.
        emu_python = str(Path(emu_conda_bin).parent / "python3")
        if Path(emu_python).exists():
            emu_prefix = [emu_python, emu_conda_bin]
        elif conda:
            # Fallback: conda run (some systems need this)
            emu_prefix = [conda, "run", "--no-capture-output", "-n", "emu", "emu"]
        else:
            emu_prefix = [emu_conda_bin]
    else:
        emu_sys = shutil.which("emu")
        if not emu_sys:
            raise RuntimeError(
                "Emu not found. Install with: conda create -n emu -c bioconda -c conda-forge emu -y")
        emu_prefix = [emu_sys]

    emu_out_dir = out_dir / "emu_output"
    emu_out_dir.mkdir(exist_ok=True)
    results = {}

    for name, fq in samples.items():
        sample_out = emu_out_dir / name
        sample_out.mkdir(exist_ok=True)
        cmd = emu_prefix + [
            "abundance",
            "--type", "map-ont",
            "--db", db_path,
            "--threads", str(threads),
            "--min-abundance", str(min_abundance),
            "--keep-files",
            "--output-dir", str(sample_out),
            str(fq),
        ]
        success = run_cmd(cmd, f"Emu: {name}")
        # Emu names output as <input_basename>_rel-abundance.tsv
        tsv_candidates = list(sample_out.glob("*rel-abundance.tsv"))
        if tsv_candidates:
            results[name] = tsv_candidates[0]
        elif not success:
            log(f"  [WARN] Emu failed for sample: {name}")

    return results

# ── Step 3: Combine Emu TSV outputs → abundance matrix ──────────────────────────
def combine_emu_outputs(emu_results: dict[str, Path],
                        out_dir: Path) -> tuple[list, dict]:
    """
    Read per-sample Emu TSVs, combine into abundance matrix.
    Returns (sample_names, taxa_dict) where
    taxa_dict[tax_id] = {taxonomy info + {sample: count}}
    """
    sample_names = sorted(emu_results.keys())
    taxa_dict = {}   # tax_id → {species, genus, family, order, class, phylum, kingdom, sample_counts}

    for sample, tsv_path in emu_results.items():
        if not tsv_path.exists():
            continue
        with open(tsv_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                tax_id = str(row.get("tax_id", "")).strip()
                if not tax_id:
                    continue
                if tax_id not in taxa_dict:
                    taxa_dict[tax_id] = {
                        "species":  row.get("species", ""),
                        "genus":    row.get("genus", ""),
                        "family":   row.get("family", ""),
                        "order":    row.get("order", ""),
                        "class":    row.get("class", ""),
                        "phylum":   row.get("phylum", ""),
                        "kingdom":  row.get("superkingdom", "Bacteria"),
                        "counts":   {s: 0 for s in sample_names},
                    }
                # Emu gives "estimated counts" as float
                est_counts = float(row.get("estimated counts", row.get("count", 0)))
                taxa_dict[tax_id]["counts"][sample] = round(est_counts)

    return sample_names, taxa_dict

# ── Step 4: Write ASV table + taxonomy table ────────────────────────────────────
def write_tables(sample_names: list, taxa_dict: dict, out_dir: Path, top_n: int):
    """Write asv_table.csv and taxonomy.csv in DADA2-compatible format."""

    # Sort taxa by total abundance desc, take topN for viz
    sorted_taxa = sorted(
        taxa_dict.items(),
        key=lambda x: sum(x[1]["counts"].values()),
        reverse=True
    )

    # asv_table.csv: rows = ASVs (tax_id), cols = samples
    asv_path = out_dir / "asv_table.csv"
    with open(asv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + sample_names)
        for tax_id, info in sorted_taxa:
            row = [tax_id] + [info["counts"][s] for s in sample_names]
            writer.writerow(row)
    log(f"  Saved: asv_table.csv  ({len(sorted_taxa)} taxa × {len(sample_names)} samples)")

    # taxonomy.csv: Kingdom, Phylum, Class, Order, Family, Genus, Species
    tax_path = out_dir / "taxonomy.csv"
    with open(tax_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ASV", "Kingdom", "Phylum", "Class", "Order",
                         "Family", "Genus", "Species"])
        for tax_id, info in sorted_taxa:
            sp = info["species"] or ""
            genus = info["genus"] or ""
            # Species field: "Genus species" if both present
            species_full = sp if sp else (f"{genus} sp." if genus else "")
            writer.writerow([
                tax_id,
                info["kingdom"] or "Bacteria",
                info["phylum"] or "",
                info["class"] or "",
                info["order"] or "",
                info["family"] or "",
                genus,
                species_full,
            ])
    log(f"  Saved: taxonomy.csv")

    return sorted_taxa

# ── Step 5: Write read_tracking.csv (simplified for ONT) ───────────────────────
def write_read_tracking(samples_raw: dict, samples_trimmed: dict,
                        out_dir: Path):
    """Write simplified read tracking: input and (optionally) trimmed counts."""
    rows = []
    for name in sorted(samples_raw.keys()):
        raw_fq   = samples_raw.get(name)
        trim_fq  = samples_trimmed.get(name)
        n_input  = count_reads(raw_fq) if raw_fq else 0
        n_final  = count_reads(trim_fq) if (trim_fq and trim_fq != raw_fq) else n_input
        rows.append({"sample": name, "input": n_input, "final": n_final})

    path = out_dir / "read_tracking.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "input", "final"])
        writer.writeheader()
        writer.writerows(rows)
    log("  Saved: read_tracking.csv")
    return rows

# ── Step 6: Write summary.json ──────────────────────────────────────────────────
def write_summary(out_dir: Path, sample_names: list, n_taxa: int,
                  read_rows: list, region: str, args):
    total_reads = sum(r["input"] for r in read_rows)
    summary = {
        "job_name":       args.job_name or out_dir.name,
        "marker":         args.marker,
        "platform":       "ONT",
        "region":         region,
        "n_samples":      len(sample_names),
        "n_taxa":         n_taxa,
        "total_reads":    total_reads,
        "has_taxonomy":   True,
        "taxonomy_level": "species",
        "db_path":        args.db_path,
        "timestamp":      datetime.now().isoformat(),
        "pipeline":       "emu",
        "emu_min_abundance": args.min_abundance,
    }
    path = out_dir / "summary.json"
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    log("  Saved: summary.json")

# ── Step 7: Run viz_pipeline.R ──────────────────────────────────────────────────
def run_viz(out_dir: Path, metadata_path: str, top_n: int, threads: int,
            marker: str):
    """Run the shared R visualization pipeline."""
    script = Path(__file__).parent / "r_scripts" / "viz_pipeline.R"
    if not script.exists():
        log("  [WARN] viz_pipeline.R not found — skipping visualizations")
        return

    cmd = [
        "Rscript", str(script),
        "--input_dir",  str(out_dir),
        "--output_dir", str(out_dir),
        "--marker",     marker,
        "--topN",       str(top_n),
        "--threads",    str(threads),
    ]
    if metadata_path and Path(metadata_path).exists():
        cmd += ["--metadata", metadata_path]

    log("  Running viz_pipeline.R ...")
    run_cmd(cmd, "viz_pipeline.R")

# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    input_dir  = Path(args.input)
    out_dir    = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log("=" * 50)
    log(f"  NextGen-Amplicon — ONT 16S Pipeline (Emu)")
    log(f"  Region: {args.region}  |  DB: {args.db_path}")
    log("=" * 50)

    # ── 1. Find samples ────────────────────────────────────────────────────────
    progress(5, "Step 1/6 — Scanning input FASTQ files")
    samples_raw = find_fastq(input_dir)
    if not samples_raw:
        log(f"[ERROR] No FASTQ files found in: {input_dir}")
        sys.exit(1)
    log(f"  Found {len(samples_raw)} samples: {', '.join(sorted(samples_raw))}")

    # ── 2. Trim primers ────────────────────────────────────────────────────────
    progress(15, "Step 2/6 — Trimming primers (cutadapt)")
    samples_trimmed = trim_primers(
        samples_raw, out_dir, args.primer_f, args.primer_r, args.threads)

    # ── 3. Run Emu ─────────────────────────────────────────────────────────────
    progress(25, f"Step 3/6 — Running Emu on {len(samples_trimmed)} samples")
    emu_results = run_emu(
        samples_trimmed, out_dir, args.db_path, args.min_abundance, args.threads)

    if not emu_results:
        log("[ERROR] Emu produced no output. Check db_path and input FASTQ files.")
        sys.exit(1)
    log(f"  Emu completed for {len(emu_results)}/{len(samples_trimmed)} samples")

    # ── 4. Combine outputs ─────────────────────────────────────────────────────
    progress(60, "Step 4/6 — Combining abundance tables")
    sample_names, taxa_dict = combine_emu_outputs(emu_results, out_dir)

    # ── 5. Write tables ────────────────────────────────────────────────────────
    progress(70, "Step 5/6 — Writing ASV & taxonomy tables")
    sorted_taxa = write_tables(sample_names, taxa_dict, out_dir, args.topN)

    # Write read tracking
    read_rows = write_read_tracking(samples_raw, samples_trimmed, out_dir)

    # Write summary
    write_summary(out_dir, sample_names, len(sorted_taxa), read_rows,
                  args.region, args)

    # ── 6. Visualizations ──────────────────────────────────────────────────────
    progress(78, "Step 6/6 — Generating plots (viz_pipeline.R)")
    run_viz(out_dir, args.metadata, args.topN, args.threads, args.marker)

    progress(100, "Complete")
    log("=" * 50)
    log("  emu_pipeline.py completed successfully")
    log("=" * 50)

if __name__ == "__main__":
    main()
