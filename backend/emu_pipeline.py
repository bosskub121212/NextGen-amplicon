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

import argparse, csv, gzip, json, os, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

# ── Region-appropriate length bounds (post primer-trim) ────────────────────────
# Used to give cutadapt a realistic --maximum-length so mis-primed / chimeric-length
# reads get discarded instead of flowing through to Emu. Falls back to a loose
# sanity bound for regions not in this table.
REGION_LENGTH_BOUNDS = {
    "V7-V8": (250, 450),
    "V1-V9": (1000, 1700),
}
DEFAULT_TRIM_BOUNDS = (50, 2000)

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
    # QC / cleanup steps (added to match a manual QIIME2 ONT pipeline used for
    # accuracy comparison — see gap_analysis_partner_vs_our_app.md 2026-08-16).
    # All have sensible defaults so no UI/main.py changes are required; each
    # step is skipped gracefully (with a warning) if its tool isn't installed.
    p.add_argument("--qc_min_qual",   type=int, default=10,   help="Chopper min mean read quality (Q)")
    p.add_argument("--qc_min_len",    type=int, default=100,  help="Chopper min read length (raw, pre-primer-trim)")
    p.add_argument("--qc_max_len",    type=int, default=2000, help="Chopper max read length (raw, pre-primer-trim)")
    p.add_argument("--skip_qc",       action="store_true",    help="Skip Chopper QC filtering step")
    p.add_argument("--skip_chimera",  action="store_true",    help="Skip vsearch chimera-removal step")
    p.add_argument("--cutadapt_error_rate", type=float, default=0.20, help="cutadapt -e (max error rate)")
    p.add_argument("--cutadapt_overlap",    type=int,   default=10,   help="cutadapt -O (min primer overlap)")
    return p.parse_args()

# ── Helpers ─────────────────────────────────────────────────────────────────────
def find_fastq(input_dir: Path) -> dict[str, Path]:
    """Find FASTQ files and map sample_name → path.

    If sample_manifest.json exists in input_dir (written by the frontend's
    manual sample/file pairing UI), use it as the source of truth for sample
    names instead of raw filenames. Each entry: {"sample": "...", "file1": "..."}.
    (file2 is ignored here — Emu/ONT reads are inherently single, long reads.)
    """
    manifest_path = input_dir / "sample_manifest.json"
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            samples = {}
            for row in manifest:
                sample = row.get("sample", "").strip()
                file1  = row.get("file1", "").strip()
                if not sample or not file1:
                    continue
                # Sanitize: sample name becomes a directory name (emu_output/<sample>/),
                # so strip anything that could be interpreted as a path separator.
                safe_sample = re.sub(r"[\\/]+", "_", sample).strip()
                if not safe_sample:
                    log(f"  [WARN] manifest sample name invalid after sanitizing: {sample!r} — skipping")
                    continue
                if safe_sample != sample:
                    log(f"  [WARN] sample name '{sample}' contained path separators — sanitized to '{safe_sample}'")
                fpath = input_dir / file1
                if not fpath.exists():
                    log(f"  [WARN] manifest file not found, skipping: {file1}")
                    continue
                samples[safe_sample] = fpath
            if samples:
                log(f"  Using sample_manifest.json ({len(samples)} sample(s))")
                return samples
            log("  [WARN] sample_manifest.json present but yielded no valid samples — falling back to filename scan")
        except Exception as e:
            log(f"  [WARN] Failed to read sample_manifest.json: {e} — falling back to filename scan")

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

def run_cmd(cmd: list, desc: str = "", env=None) -> bool:
    """Run a subprocess, stream stdout/stderr, return success."""
    if desc:
        log(f"  → {desc}")
    log(f"  CMD: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=False, env=env)
    if result.returncode != 0:
        log(f"  [WARN] Command returned {result.returncode}: {desc}")
        return False
    return True

def _rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (IUPAC aware)."""
    table = str.maketrans(
        "ACGTacgtMKRYWSBVHDNmkrywsbvhdn",
        "TGCAtgcaKMYRWSVBDHNkmyrwsvbdhn"
    )
    return seq.translate(table)[::-1]

def _fastq_iter(path: Path):
    """Yield (header, seq, plus, qual) 4-line records from a fastq(.gz) file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not q:
                break
            yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

# ── Step 1: QC filter raw reads with Chopper ────────────────────────────────────
def qc_filter_chopper(samples: dict, out_dir: Path, min_qual: int, min_len: int,
                      max_len: int, threads: int) -> dict[str, Path]:
    """Quality + length filter raw ONT reads before primer trimming (matches a
    typical Chopper pre-processing step used ahead of manual QIIME2 ONT pipelines).
    Skips gracefully (keeps raw reads) if chopper isn't installed."""
    chopper = shutil.which("chopper")
    if not chopper:
        log("  [WARN] chopper not found — skipping QC filtering "
            "(install with: conda install -c bioconda chopper)")
        return samples

    qc_dir = out_dir / "qc_filtered"
    qc_dir.mkdir(exist_ok=True)
    filtered = {}

    for name, fq in samples.items():
        out_fq = qc_dir / f"{name}.fastq.gz"
        cat_cmd = ["zcat", str(fq)] if str(fq).endswith(".gz") else ["cat", str(fq)]
        chopper_cmd = [chopper, "-q", str(min_qual), "-l", str(min_len),
                       "--maxlength", str(max_len), "-t", str(threads)]
        log(f"  → QC filter (chopper): {name}")
        log(f"  CMD: {' '.join(cat_cmd)} | {' '.join(str(c) for c in chopper_cmd)} | gzip")
        try:
            with open(out_fq, "wb") as out_f:
                p1 = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(chopper_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
                p1.stdout.close()
                p3 = subprocess.Popen(["gzip"], stdin=p2.stdout, stdout=out_f)
                p2.stdout.close()
                p3.communicate()
            if out_fq.exists() and out_fq.stat().st_size > 0:
                filtered[name] = out_fq
            else:
                log(f"  [WARN] chopper produced empty output for {name} — using unfiltered reads")
                filtered[name] = fq
        except Exception as e:
            log(f"  [WARN] chopper failed for {name}: {e} — using unfiltered reads")
            filtered[name] = fq

    return filtered

# ── Step 2: Trim primers with cutadapt ─────────────────────────────────────────
def trim_primers(samples: dict, out_dir: Path, primer_f: str, primer_r: str,
                 threads: int, error_rate: float = 0.20, overlap: int = 10,
                 min_len: int = 50, max_len: int = 0) -> dict[str, Path]:
    """Trim primers from ONT single-end reads with cutadapt.
    Uses single-end mode only (-g / -a) — never -G/-A which trigger paired-end.

    error_rate/overlap default to a tolerant 20% / 10bp overlap (vs cutadapt's
    own defaults of 10% / 3bp) — ONT reads have much higher per-base error than
    Illumina, so an untuned cutadapt call silently --discard-untrimmed's a large
    fraction of genuinely-good reads whose primer site has a couple of errors.
    --revcomp checks both orientations, since ONT reads aren't orientation-fixed.
    min_len/max_len bound the POST-trim amplicon length (region-appropriate,
    e.g. ~250-450bp for V7-V8) to catch mis-primed / chimeric-length artifacts."""
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
        # ONT = single-end: forward primer at 5' end, RC(reverse primer) at 3' end
        # NEVER use -G/-A — those flags trigger paired-end mode and require two files.
        cmd = [cutadapt, "-j", str(threads), "-e", str(error_rate), "-O", str(overlap)]
        if primer_f:
            cmd += ["-g", primer_f]          # 5' forward primer
        if primer_r:
            cmd += ["-a", _rev_comp(primer_r)]  # 3' reverse primer (as RC)
        if primer_f or primer_r:
            cmd += ["--revcomp"]             # also check the RC strand, keep whichever matches
        cmd += ["--discard-untrimmed", "-m", str(min_len)]
        if max_len:
            cmd += ["-M", str(max_len)]
        cmd += ["-o", str(out_fq), str(fq)]
        run_cmd(cmd, f"Trim {name}")
        trimmed[name] = out_fq if out_fq.exists() else fq

    return trimmed

# ── Step 3: Remove chimeras with vsearch (reference-based) ─────────────────────
def remove_chimeras(samples: dict, out_dir: Path, db_path: str,
                    threads: int) -> dict[str, Path]:
    """Filter out chimeric reads with vsearch --uchime_ref against the Emu
    database's own sequences.fasta. Reference-based (not --uchime_denovo)
    because raw/trimmed ONT reads are ~100% unique — de novo chimera detection
    relies on abundance ratios between near-identical sequences that don't
    meaningfully exist here. Skips gracefully if vsearch or the reference
    fasta isn't available."""
    vsearch = shutil.which("vsearch")
    if not vsearch:
        log("  [WARN] vsearch not found — skipping chimera removal "
            "(install with: conda install -c bioconda vsearch)")
        return samples

    ref_fasta = Path(db_path) / "sequences.fasta"
    if not ref_fasta.exists():
        log(f"  [WARN] {ref_fasta} not found — skipping chimera removal "
            "(only works with Emu DBs built by our build_emu_db*.py scripts)")
        return samples

    chim_dir = out_dir / "chimera_filtered"
    chim_dir.mkdir(exist_ok=True)
    cleaned = {}

    for name, fq in samples.items():
        fa_in  = chim_dir / f"{name}.in.fasta"
        fa_ok  = chim_dir / f"{name}.nonchimeric.fasta"
        fa_bad = chim_dir / f"{name}.chimeric.fasta"
        out_fq = chim_dir / f"{name}.fastq.gz"

        # cutadapt/chopper output is always plain fastq or fastq.gz — convert to
        # fasta since --uchime_ref requires it.
        n_reads = 0
        with open(fa_in, "w") as out_f:
            for h, seq, _, _ in _fastq_iter(fq):
                rid = h[1:].split()[0] if h.startswith("@") else h.split()[0]
                out_f.write(f">{rid}\n{seq}\n")
                n_reads += 1

        if n_reads == 0:
            log(f"  [WARN] {name}: no reads to check for chimeras — skipping")
            cleaned[name] = fq
            continue

        cmd = [vsearch, "--uchime_ref", str(fa_in), "--db", str(ref_fasta),
               "--nonchimeras", str(fa_ok), "--chimeras", str(fa_bad),
               "--threads", str(threads), "--fasta_width", "0"]
        ok = run_cmd(cmd, f"Chimera check: {name}")

        if not ok or not fa_ok.exists():
            log(f"  [WARN] vsearch chimera check failed for {name} — keeping all reads")
            cleaned[name] = fq
            for p in (fa_in, fa_ok, fa_bad):
                try: p.unlink()
                except Exception: pass
            continue

        keep_ids = set()
        with open(fa_ok) as f:
            for line in f:
                if line.startswith(">"):
                    keep_ids.add(line[1:].strip().split()[0])

        n_kept = 0
        with gzip.open(out_fq, "wt") as out_f:
            for h, seq, plus, qual in _fastq_iter(fq):
                rid = h[1:].split()[0] if h.startswith("@") else h.split()[0]
                if rid in keep_ids:
                    out_f.write(f"{h}\n{seq}\n{plus}\n{qual}\n")
                    n_kept += 1

        log(f"  {name}: {n_reads} reads → {n_kept} non-chimeric ({n_reads - n_kept} chimeras removed)")
        cleaned[name] = out_fq if n_kept > 0 else fq

        for p in (fa_in, fa_ok, fa_bad):
            try: p.unlink()
            except Exception: pass

    return cleaned

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

    # Build environment: inject emu conda bin dir into PATH so minimap2 is found
    emu_env = os.environ.copy()
    if emu_conda_bin:
        # Prefer calling the emu env's Python directly — avoids issues where
        # 'conda run --no-capture-output' doesn't fully activate the env and
        # fails to find packages like pysam even when they're installed.
        emu_python = str(Path(emu_conda_bin).parent / "python3")
        emu_bin_dir = str(Path(emu_conda_bin).parent)
        # Always add emu env bin to PATH so minimap2 and other tools are found
        emu_env["PATH"] = emu_bin_dir + os.pathsep + emu_env.get("PATH", "")
        log(f"  emu PATH: {emu_bin_dir} prepended")
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
        sample_out.mkdir(parents=True, exist_ok=True)
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
        success = run_cmd(cmd, f"Emu: {name}", env=emu_env)
        # Prefer the thresholded file (already filtered by --min-abundance).
        # Emu names: <basename>_rel-abundance-threshold-<N>.tsv (current)
        #            <basename>_rel-abundance.tsv              (fallback / older)
        tsv_thresh = list(sample_out.glob("*rel-abundance-threshold*.tsv"))
        tsv_plain  = list(sample_out.glob("*rel-abundance.tsv"))
        # Exclude threshold files from plain list
        tsv_plain  = [p for p in tsv_plain if "threshold" not in p.name]
        if tsv_thresh:
            results[name] = tsv_thresh[0]
        elif tsv_plain:
            results[name] = tsv_plain[0]
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
                # Emu outputs relative abundance (0–1).
                # Try "estimated counts" (old Emu ≤2.x), then "abundance" (current).
                raw = row.get("estimated counts") or row.get("estimated_counts") or row.get("abundance", "0")
                val = float(raw) if raw else 0.0
                # If value looks like relative abundance (≤ 1.0), scale to ppm
                # so phyloseq and diversity tools get integer-like counts.
                if val <= 1.0:
                    est_counts = round(val * 1_000_000)
                else:
                    est_counts = round(val)
                taxa_dict[tax_id]["counts"][sample] = est_counts

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
def write_read_tracking(samples_raw: dict, samples_qc: dict, samples_trimmed: dict,
                        samples_final: dict, out_dir: Path):
    """Write read tracking through each pipeline stage: input → QC-filtered →
    primer-trimmed → chimera-filtered (final). Any stage that was skipped
    (tool not installed) will simply show the same count as the stage before it."""
    rows = []
    for name in sorted(samples_raw.keys()):
        raw_fq   = samples_raw.get(name)
        qc_fq    = samples_qc.get(name)
        trim_fq  = samples_trimmed.get(name)
        final_fq = samples_final.get(name)
        n_input  = count_reads(raw_fq) if raw_fq else 0
        n_qc     = count_reads(qc_fq)    if (qc_fq    and qc_fq    != raw_fq)  else n_input
        n_trim   = count_reads(trim_fq)  if (trim_fq  and trim_fq  != qc_fq)   else n_qc
        n_final  = count_reads(final_fq) if (final_fq and final_fq != trim_fq) else n_trim
        rows.append({"sample": name, "input": n_input, "qc_filtered": n_qc,
                     "primer_trimmed": n_trim, "final": n_final})

    path = out_dir / "read_tracking.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "input", "qc_filtered", "primer_trimmed", "final"])
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
    progress(5, "Step 1/8 — Scanning input FASTQ files")
    samples_raw = find_fastq(input_dir)
    if not samples_raw:
        log(f"[ERROR] No FASTQ files found in: {input_dir}")
        sys.exit(1)
    log(f"  Found {len(samples_raw)} samples: {', '.join(sorted(samples_raw))}")

    # ── 2. QC filter (Chopper) ─────────────────────────────────────────────────
    if args.skip_qc:
        log("  Step 2/8 — QC filtering skipped (--skip_qc)")
        samples_qc = samples_raw
    else:
        progress(12, "Step 2/8 — QC filtering (chopper)")
        samples_qc = qc_filter_chopper(
            samples_raw, out_dir, args.qc_min_qual, args.qc_min_len,
            args.qc_max_len, args.threads)

    # ── 3. Trim primers ────────────────────────────────────────────────────────
    progress(22, "Step 3/8 — Trimming primers (cutadapt)")
    trim_min_len, trim_max_len = REGION_LENGTH_BOUNDS.get(args.region, DEFAULT_TRIM_BOUNDS)
    samples_trimmed = trim_primers(
        samples_qc, out_dir, args.primer_f, args.primer_r, args.threads,
        error_rate=args.cutadapt_error_rate, overlap=args.cutadapt_overlap,
        min_len=trim_min_len, max_len=trim_max_len)

    # ── 4. Remove chimeras (vsearch, reference-based) ──────────────────────────
    if args.skip_chimera:
        log("  Step 4/8 — Chimera removal skipped (--skip_chimera)")
        samples_clean = samples_trimmed
    else:
        progress(32, "Step 4/8 — Chimera removal (vsearch)")
        samples_clean = remove_chimeras(
            samples_trimmed, out_dir, args.db_path, args.threads)

    # ── 5. Run Emu ─────────────────────────────────────────────────────────────
    progress(42, f"Step 5/8 — Running Emu on {len(samples_clean)} samples")
    emu_results = run_emu(
        samples_clean, out_dir, args.db_path, args.min_abundance, args.threads)

    if not emu_results:
        log("[ERROR] Emu produced no output. Check db_path and input FASTQ files.")
        sys.exit(1)
    log(f"  Emu completed for {len(emu_results)}/{len(samples_clean)} samples")

    # ── 6. Combine outputs ─────────────────────────────────────────────────────
    progress(62, "Step 6/8 — Combining abundance tables")
    sample_names, taxa_dict = combine_emu_outputs(emu_results, out_dir)

    # ── 7. Write tables ────────────────────────────────────────────────────────
    progress(72, "Step 7/8 — Writing ASV & taxonomy tables")
    sorted_taxa = write_tables(sample_names, taxa_dict, out_dir, args.topN)

    # Write read tracking
    read_rows = write_read_tracking(samples_raw, samples_qc, samples_trimmed,
                                    samples_clean, out_dir)

    # Write summary
    write_summary(out_dir, sample_names, len(sorted_taxa), read_rows,
                  args.region, args)

    # ── 8. Visualizations ──────────────────────────────────────────────────────
    progress(82, "Step 8/8 — Generating plots (viz_pipeline.R)")
    run_viz(out_dir, args.metadata, args.topN, args.threads, args.marker)

    progress(100, "Complete")
    log("=" * 50)
    log("  emu_pipeline.py completed successfully")
    log("=" * 50)

if __name__ == "__main__":
    main()
