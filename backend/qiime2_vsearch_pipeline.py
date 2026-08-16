#!/usr/bin/env python3
"""
NextGen-Amplicon — QIIME2/VSEARCH OTU-picking pipeline for 16S-family markers
==============================================================================
This is the "QIIME2 (VSEARCH OTU)" Sequencer Type option in the 16S rRNA
pipeline card — an alternative to DADA2 that mirrors a manually-run QIIME2/
VSEARCH OTU-clustering workflow (Chopper QC -> tolerant cutadapt -> VSEARCH
dereplicate -> reference-based chimera removal -> OTU clustering @ X% identity
-> map reads back to OTU centroids -> taxonomy). Built to let results be
compared directly against a partner lab's manual VSEARCH-based pipeline,
and to process high per-read-error long reads (e.g. ONT) that DADA2's
exact-dereplication approach cannot handle regardless of amplicon length.

Chimera removal uses vsearch --uchime_ref (against the same taxonomy
reference database) rather than the partner script's literal --uchime_denovo:
de novo detection needs abundance differences between near-identical sequences
to work, which high-error long reads don't reliably have after dereplication
(virtually every read is unique), and vsearch's --uchime_denovo is single-
threaded regardless of --threads — a severe bottleneck on ONT-scale unique-
sequence counts. Same reasoning already applied to emu_pipeline.py's ONT-16S
chimera step. See resolve_reference_fasta()/remove_chimeras_ref() below.

Calls the `vsearch` / `cutadapt` / `chopper` binaries directly (NOT the
`qiime` CLI) — same tools already required for the Emu ONT-16S pipeline's
QC/chimera steps, so no extra multi-GB QIIME2 conda environment is needed
on deployment machines. Taxonomy reference accepts EITHER an Emu-format
directory (sequences.fasta + seq2taxid.tsv + taxonomy.tsv) OR a plain
DADA2-trainset-style SILVA FASTA file (e.g. silva_nr99_v138.2_toGenus_
trainset.fa.gz) — the same files already used by the DADA2 pipeline's own
Taxonomy Database picker, auto-detected from whether db_path is a
directory or a file. See assign_taxonomy() below.

Usage (called by main.py):
    python qiime2_vsearch_pipeline.py --input <fastq_dir> --output <out_dir>
                                      --db_path <reference_db_dir>
                                      [--primer_f SEQ --primer_r SEQ]
                                      [--min_len 250 --max_len 450]
                                      [--otu_similarity 0.97]
                                      [--threads 4]
                                      [--metadata metadata.csv]
                                      [--topN 30]
                                      [--job_name NAME]
"""

import argparse, csv, gzip, json, os, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

# ── Progress reporter (mirrors dada2_pipeline.R / emu_pipeline.py format) ──────
def progress(pct: int, label: str):
    print(f"PROGRESS:{pct}|{label}", flush=True)

def log(msg: str):
    print(msg, flush=True)

# ── Argument parsing ────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",         required=True)
    p.add_argument("--output",        required=True)
    p.add_argument("--db_path",       required=True, help="Path to reference database directory (Emu format)")
    p.add_argument("--primer_f",      default="",    help="Forward primer sequence")
    p.add_argument("--primer_r",      default="",    help="Reverse primer sequence")
    p.add_argument("--min_len",       type=int, default=50,   help="Min read length post primer-trim")
    p.add_argument("--max_len",       type=int, default=2000, help="Max read length post primer-trim")
    p.add_argument("--otu_similarity", type=float, default=0.97, help="VSEARCH OTU clustering identity (0-1)")
    p.add_argument("--tax_identity",   type=float, default=0.80, help="VSEARCH taxonomy-assignment identity (0-1)")
    p.add_argument("--threads",       type=int,   default=4)
    p.add_argument("--metadata",      default="")
    p.add_argument("--topN",          type=int,   default=30)
    p.add_argument("--job_name",      default="")
    p.add_argument("--marker",        default="16S")
    p.add_argument("--cutadapt_error_rate", type=float, default=0.20)
    p.add_argument("--cutadapt_overlap",    type=int,   default=10)
    p.add_argument("--qc_min_qual",   type=int, default=10)
    p.add_argument("--skip_qc",       action="store_true")
    return p.parse_args()

# ── Helpers (self-contained — duplicated from emu_pipeline.py by design, so ────
#    each pipeline script stays independently deployable/editable) ─────────────
def find_fastq(input_dir: Path) -> dict[str, Path]:
    """Find FASTQ files and map sample_name -> path (manifest-aware, single-end)."""
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
                safe_sample = re.sub(r"[\\/]+", "_", sample).strip()
                if not safe_sample:
                    log(f"  [WARN] manifest sample name invalid after sanitizing: {sample!r} — skipping")
                    continue
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
    if desc:
        log(f"  → {desc}")
    log(f"  CMD: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=False, env=env)
    if result.returncode != 0:
        log(f"  [WARN] Command returned {result.returncode}: {desc}")
        return False
    return True

def _rev_comp(seq: str) -> str:
    table = str.maketrans(
        "ACGTacgtMKRYWSBVHDNmkrywsbvhdn",
        "TGCAtgcaKMYRWSVBDHNkmyrwsvbdhn"
    )
    return seq.translate(table)[::-1]

def _fastq_iter(path: Path):
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
def qc_filter_chopper(samples: dict, out_dir: Path, min_qual: int,
                      threads: int) -> dict[str, Path]:
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
        chopper_cmd = [chopper, "-q", str(min_qual), "-t", str(threads)]
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
                 threads: int, error_rate: float, overlap: int,
                 min_len: int, max_len: int) -> dict[str, Path]:
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
        cmd = [cutadapt, "-j", str(threads), "-e", str(error_rate), "-O", str(overlap)]
        if primer_f:
            cmd += ["-g", primer_f]
        if primer_r:
            cmd += ["-a", _rev_comp(primer_r)]
        if primer_f or primer_r:
            cmd += ["--revcomp"]
        cmd += ["--discard-untrimmed", "-m", str(min_len)]
        if max_len:
            cmd += ["-M", str(max_len)]
        cmd += ["-o", str(out_fq), str(fq)]
        run_cmd(cmd, f"Trim {name}")
        trimmed[name] = out_fq if out_fq.exists() else fq
    return trimmed

# ── Step 3: Pool all samples -> FASTA (for OTU discovery only) ─────────────────
def pool_samples_to_fasta(samples: dict, out_path: Path) -> int:
    n = 0
    with open(out_path, "w") as out_f:
        for name, fq in samples.items():
            for i, (h, seq, _, _) in enumerate(_fastq_iter(fq)):
                out_f.write(f">{name}_{i}\n{seq}\n")
                n += 1
    return n

# ── Step 4: Dereplicate pooled reads ────────────────────────────────────────────
def dereplicate(vsearch: str, pooled_fasta: Path, out_dir: Path, threads: int) -> Path | None:
    uniques = out_dir / "uniques.fasta"
    cmd = [vsearch, "--derep_fulllength", str(pooled_fasta),
           "--sizeout", "--relabel", "Uniq", "--output", str(uniques),
           "--threads", str(threads)]
    ok = run_cmd(cmd, "Dereplicate pooled reads")
    return uniques if ok and uniques.exists() else None

# ── Step 5: Reference-based chimera removal (on abundance-sorted uniques) ──────
def resolve_reference_fasta(db_path: str) -> Path | None:
    """Locate a reference FASTA for --uchime_ref (and reused by assign_taxonomy's
    Emu-format branch). Accepts a plain FASTA/.fa.gz FILE directly, or an Emu-format
    DIRECTORY containing sequences.fasta — checked both at db_path itself and at
    db_path's parent, since `emu build-database <name>` (the step our
    build_emu_db*.py scripts print as the "next command") creates the compiled
    index in a subdirectory named <name>, one level below where sequences.fasta
    actually lives (e.g. .../silva_qiime2_v7v8/silva_qiime2_v7v8/). If ont_db_path
    was registered as that inner compiled-index folder, sequences.fasta is one
    directory up."""
    dbp = Path(db_path)
    if dbp.is_file():
        return dbp
    if dbp.is_dir():
        direct = dbp / "sequences.fasta"
        if direct.exists():
            return direct
        parent = dbp.parent / "sequences.fasta"
        if parent.exists():
            return parent
    return None

def remove_chimeras_ref(vsearch: str, uniques_fasta: Path, out_dir: Path,
                        ref_fasta: Path, threads: int) -> Path:
    """Reference-based chimera detection (vsearch --uchime_ref) against the same
    taxonomy reference database, checking each unique sequence individually.
    Used instead of --uchime_denovo: de novo detection relies on abundance ratios
    between near-identical sequences to spot a chimera as a rare recombinant of
    two more-abundant parents, but ONT reads are ~100% unique after dereplication
    (no real signal for it to work with) — and critically, vsearch's own
    --uchime_denovo implementation is single-threaded regardless of --threads
    ("the uchime_denovo command does not support multithreading"), making it a
    severe bottleneck on the large unique-sequence sets ONT data produces.
    --uchime_ref does not have either problem and fully supports --threads."""
    sorted_fa = out_dir / "uniques.sorted.fasta"
    run_cmd([vsearch, "--sortbysize", str(uniques_fasta), "--output", str(sorted_fa),
              "--threads", str(threads)], "Sort uniques by abundance")
    src = sorted_fa if sorted_fa.exists() else uniques_fasta
    nonchim = out_dir / "uniques.nonchimeras.fasta"
    ok = run_cmd([vsearch, "--uchime_ref", str(src), "--db", str(ref_fasta),
                  "--nonchimeras", str(nonchim), "--threads", str(threads)],
                 "Reference-based chimera detection (uchime_ref)")
    return nonchim if ok and nonchim.exists() else src

# ── Step 6: OTU clustering ──────────────────────────────────────────────────────
def cluster_otus(vsearch: str, nonchim_fasta: Path, out_dir: Path,
                 similarity: float, threads: int) -> Path | None:
    otus = out_dir / "otus.fasta"
    cmd = [vsearch, "--cluster_size", str(nonchim_fasta), "--id", str(similarity),
           "--centroids", str(otus), "--relabel", "OTU_", "--sizein", "--sizeout",
           "--threads", str(threads)]
    ok = run_cmd(cmd, f"Cluster OTUs @ {similarity * 100:.1f}% identity")
    return otus if ok and otus.exists() else None

# ── Step 7: Map each sample's reads back to OTU centroids ──────────────────────
def map_reads_to_otus(vsearch: str, samples: dict, otus_fasta: Path, out_dir: Path,
                      similarity: float, threads: int):
    """Returns (sample_names, otu_counts) where otu_counts[otu_id] = {sample: count}."""
    map_dir = out_dir / "otu_mapping"
    map_dir.mkdir(exist_ok=True)
    sample_names = sorted(samples.keys())
    otu_counts: dict[str, dict[str, int]] = {}
    mapped_totals = {}

    for name in sample_names:
        fq = samples[name]
        fa_path = map_dir / f"{name}.reads.fasta"
        n_reads = 0
        with open(fa_path, "w") as f:
            for i, (h, seq, _, _) in enumerate(_fastq_iter(fq)):
                f.write(f">{name}_{i}\n{seq}\n")
                n_reads += 1

        hits_path = map_dir / f"{name}.hits.tsv"
        cmd = [vsearch, "--usearch_global", str(fa_path), "--db", str(otus_fasta),
               "--id", str(similarity), "--top_hits_only",
               "--userout", str(hits_path), "--userfields", "query+target",
               "--threads", str(threads)]
        run_cmd(cmd, f"Map reads to OTUs: {name}")

        n_mapped = 0
        if hits_path.exists():
            with open(hits_path) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    otu_id = parts[1]
                    otu_counts.setdefault(otu_id, {s: 0 for s in sample_names})
                    otu_counts[otu_id][name] += 1
                    n_mapped += 1
        mapped_totals[name] = (n_reads, n_mapped)
        log(f"  {name}: {n_mapped}/{n_reads} reads mapped to an OTU")

        for p in (fa_path, hits_path):
            try: p.unlink()
            except Exception: pass

    return sample_names, otu_counts, mapped_totals

# ── Step 8: Assign taxonomy to OTU centroids ────────────────────────────────────
# Supports TWO reference formats, auto-detected from db_path:
#   (a) a DIRECTORY in Emu format: sequences.fasta + seq2taxid.tsv + taxonomy.tsv
#       (e.g. emu_qiime2_v7v8 built from a partner lab's QIIME2 .qza artifacts)
#   (b) a single FILE — a standard DADA2-trainset-style SILVA FASTA, e.g.
#       silva_nr99_v138.2_toGenus_trainset.fa.gz — the same files already used
#       by the DADA2 pipeline's own Taxonomy Database picker, no rebuild needed.
def assign_taxonomy(vsearch: str, otus_fasta: Path, db_path: str, out_dir: Path,
                    identity: float, threads: int) -> dict[str, dict]:
    dbp = Path(db_path)
    if dbp.is_dir():
        return _assign_taxonomy_emu_format(vsearch, otus_fasta, dbp, out_dir, identity, threads)
    if dbp.is_file():
        return _assign_taxonomy_trainset_format(vsearch, otus_fasta, dbp, out_dir, identity, threads)
    log(f"  [WARN] reference database not found at {db_path} — taxonomy will be blank")
    return {}

def _assign_taxonomy_emu_format(vsearch: str, otus_fasta: Path, db_dir: Path, out_dir: Path,
                                identity: float, threads: int) -> dict[str, dict]:
    ref_fasta     = db_dir / "sequences.fasta"
    seq2tax_path  = db_dir / "seq2taxid.tsv"
    taxonomy_path = db_dir / "taxonomy.tsv"
    if not (ref_fasta.exists() and seq2tax_path.exists() and taxonomy_path.exists()):
        log(f"  [WARN] Emu-format reference incomplete at {db_dir} "
            "(need sequences.fasta + seq2taxid.tsv + taxonomy.tsv) — taxonomy will be blank")
        return {}

    seq2tax = {}
    with open(seq2tax_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                seq2tax[parts[0]] = parts[1]

    tax_lookup: dict[str, dict] = {}
    with open(taxonomy_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            tax_lookup[parts[0]] = dict(zip(header, parts))

    hits_path = out_dir / "otu_taxonomy_hits.tsv"
    cmd = [vsearch, "--usearch_global", str(otus_fasta), "--db", str(ref_fasta),
           "--id", str(identity), "--top_hits_only",
           "--userout", str(hits_path), "--userfields", "query+target+id",
           "--threads", str(threads)]
    run_cmd(cmd, f"Assign OTU taxonomy (vsearch @ {identity * 100:.0f}% identity vs. Emu-format reference)")

    otu_tax: dict[str, dict] = {}
    if hits_path.exists():
        with open(hits_path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                otu_id, target = parts[0], parts[1]
                if otu_id in otu_tax:
                    continue
                taxid = seq2tax.get(target)
                otu_tax[otu_id] = tax_lookup.get(taxid, {}) if taxid else {}
    n_assigned = sum(1 for v in otu_tax.values() if v)
    log(f"  {n_assigned}/{len(otu_tax)} OTUs got a taxonomy hit")
    return otu_tax

def _assign_taxonomy_trainset_format(vsearch: str, otus_fasta: Path, ref_fasta: Path, out_dir: Path,
                                     identity: float, threads: int) -> dict[str, dict]:
    """DADA2-trainset-style FASTA (e.g. silva_nr99_v138.2_toGenus_trainset.fa.gz). Headers ARE
    the semicolon-delimited lineage (Kingdom;Phylum;Class;Order;Family;Genus;) — DADA2's
    assignTaxonomy() format, capped at 6 ranks (no species — matches how the DADA2 pipeline
    already treats these same files). vsearch reads .gz input natively, no decompression needed,
    and reports the matched header verbatim as `target`, so no separate id->lineage lookup file
    is required at all."""
    hits_path = out_dir / "otu_taxonomy_hits.tsv"
    cmd = [vsearch, "--usearch_global", str(otus_fasta), "--db", str(ref_fasta),
           "--id", str(identity), "--top_hits_only",
           "--userout", str(hits_path), "--userfields", "query+target",
           "--threads", str(threads)]
    run_cmd(cmd, f"Assign OTU taxonomy (vsearch @ {identity * 100:.0f}% identity vs. SILVA trainset FASTA)")

    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus"]
    otu_tax: dict[str, dict] = {}
    if hits_path.exists():
        with open(hits_path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                otu_id, lineage = parts[0], parts[1]
                if otu_id in otu_tax:
                    continue
                pieces = [p.strip() for p in lineage.split(";") if p.strip()]
                otu_tax[otu_id] = dict(zip(ranks, pieces))
    n_assigned = sum(1 for v in otu_tax.values() if v)
    log(f"  {n_assigned}/{len(otu_tax)} OTUs got a taxonomy hit (species rank unavailable — "
        "trainset format caps at genus, same limitation as the DADA2 pipeline)")
    return otu_tax

# ── Step 9: Write ASV/taxonomy tables (same format as emu_pipeline.py) ─────────
def write_tables(sample_names: list, otu_counts: dict, otu_tax: dict, out_dir: Path):
    sorted_otus = sorted(otu_counts.items(), key=lambda x: sum(x[1].values()), reverse=True)

    asv_path = out_dir / "asv_table.csv"
    with open(asv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + sample_names)
        for otu_id, counts in sorted_otus:
            writer.writerow([otu_id] + [counts[s] for s in sample_names])
    log(f"  Saved: asv_table.csv  ({len(sorted_otus)} OTUs x {len(sample_names)} samples)")

    tax_path = out_dir / "taxonomy.csv"
    with open(tax_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ASV", "Kingdom", "Phylum", "Class", "Order",
                         "Family", "Genus", "Species"])
        for otu_id, _ in sorted_otus:
            row = otu_tax.get(otu_id, {})
            genus = row.get("genus", "") or ""
            sp    = row.get("species", "") or ""
            species_full = sp if sp else (f"{genus} sp." if genus else "")
            writer.writerow([
                otu_id,
                row.get("superkingdom", "") or "Bacteria",
                row.get("phylum", "") or "",
                row.get("class", "") or "",
                row.get("order", "") or "",
                row.get("family", "") or "",
                genus,
                species_full,
            ])
    log("  Saved: taxonomy.csv")
    return sorted_otus

# ── Step 10: read_tracking.csv (multi-stage) ────────────────────────────────────
def write_read_tracking(samples_raw: dict, samples_qc: dict, samples_trimmed: dict,
                        mapped_totals: dict, out_dir: Path):
    rows = []
    for name in sorted(samples_raw.keys()):
        raw_fq  = samples_raw.get(name)
        qc_fq   = samples_qc.get(name)
        trim_fq = samples_trimmed.get(name)
        n_input = count_reads(raw_fq) if raw_fq else 0
        n_qc    = count_reads(qc_fq)   if (qc_fq   and qc_fq   != raw_fq) else n_input
        n_trim  = count_reads(trim_fq) if (trim_fq and trim_fq != qc_fq)  else n_qc
        n_final = mapped_totals.get(name, (n_trim, n_trim))[1]
        rows.append({"sample": name, "input": n_input, "qc_filtered": n_qc,
                     "primer_trimmed": n_trim, "final": n_final})
    path = out_dir / "read_tracking.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "input", "qc_filtered", "primer_trimmed", "final"])
        writer.writeheader()
        writer.writerows(rows)
    log("  Saved: read_tracking.csv")
    return rows

# ── Step 11: summary.json ───────────────────────────────────────────────────────
def write_summary(out_dir: Path, sample_names: list, n_otus: int,
                  read_rows: list, args):
    total_reads = sum(r["input"] for r in read_rows)
    summary = {
        "job_name":       args.job_name or out_dir.name,
        "marker":         args.marker,
        "platform":       "QIIME2-VSEARCH",
        "n_samples":      len(sample_names),
        "n_taxa":         n_otus,
        "total_reads":    total_reads,
        "has_taxonomy":   True,
        "taxonomy_level": "species",
        "db_path":        args.db_path,
        "timestamp":      datetime.now().isoformat(),
        "pipeline":       "qiime2_vsearch",
        "otu_similarity": args.otu_similarity,
    }
    path = out_dir / "summary.json"
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    log("  Saved: summary.json")

# ── Step 12: Run viz_pipeline.R ─────────────────────────────────────────────────
def run_viz(out_dir: Path, metadata_path: str, top_n: int, threads: int, marker: str):
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
    input_dir = Path(args.input)
    out_dir   = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log("=" * 50)
    log("  NextGen-Amplicon — QIIME2/VSEARCH OTU Pipeline")
    log(f"  OTU identity: {args.otu_similarity*100:.1f}%  |  DB: {args.db_path}")
    log("=" * 50)

    vsearch = shutil.which("vsearch")
    if not vsearch:
        log("[ERROR] vsearch not found — this pipeline requires it (unlike Emu's optional "
            "chimera step, VSEARCH OTU clustering IS this pipeline's core algorithm).")
        log("        Install with: conda install -c bioconda vsearch")
        sys.exit(1)

    # ── 1. Find samples ────────────────────────────────────────────────────────
    progress(5, "Step 1/9 — Scanning input FASTQ files")
    samples_raw = find_fastq(input_dir)
    if not samples_raw:
        log(f"[ERROR] No FASTQ files found in: {input_dir}")
        sys.exit(1)
    log(f"  Found {len(samples_raw)} samples: {', '.join(sorted(samples_raw))}")

    # ── 2. QC filter ────────────────────────────────────────────────────────────
    if args.skip_qc:
        log("  Step 2/9 — QC filtering skipped (--skip_qc)")
        samples_qc = samples_raw
    else:
        progress(12, "Step 2/9 — QC filtering (chopper)")
        samples_qc = qc_filter_chopper(samples_raw, out_dir, args.qc_min_qual, args.threads)

    # ── 3. Trim primers ─────────────────────────────────────────────────────────
    progress(20, "Step 3/9 — Trimming primers (cutadapt)")
    samples_trimmed = trim_primers(
        samples_qc, out_dir, args.primer_f, args.primer_r, args.threads,
        error_rate=args.cutadapt_error_rate, overlap=args.cutadapt_overlap,
        min_len=args.min_len, max_len=args.max_len)

    # ── 4. Pool + dereplicate ───────────────────────────────────────────────────
    progress(30, "Step 4/9 — Pooling samples & dereplicating (vsearch)")
    pooled_fasta = out_dir / "pooled.fasta"
    n_pooled = pool_samples_to_fasta(samples_trimmed, pooled_fasta)
    log(f"  Pooled {n_pooled} reads across {len(samples_trimmed)} samples")
    uniques = dereplicate(vsearch, pooled_fasta, out_dir, args.threads)
    if not uniques:
        log("[ERROR] Dereplication failed — check vsearch install and pooled reads.")
        sys.exit(1)

    # ── 5. Chimera removal (reference-based, on abundance-sorted uniques) ──────
    progress(42, "Step 5/9 — Chimera removal vs. reference (vsearch)")
    ref_fasta_chim = resolve_reference_fasta(args.db_path)
    if ref_fasta_chim:
        nonchim = remove_chimeras_ref(vsearch, uniques, out_dir, ref_fasta_chim, args.threads)
    else:
        log(f"  [WARN] No reference FASTA found at {args.db_path} for chimera check "
            "— skipping chimera removal")
        nonchim = uniques

    # ── 6. OTU clustering ────────────────────────────────────────────────────────
    progress(52, f"Step 6/9 — Clustering OTUs @ {args.otu_similarity*100:.1f}% identity")
    otus_fasta = cluster_otus(vsearch, nonchim, out_dir, args.otu_similarity, args.threads)
    if not otus_fasta:
        log("[ERROR] OTU clustering failed.")
        sys.exit(1)

    # ── 7. Map reads back to OTU centroids ──────────────────────────────────────
    progress(64, f"Step 7/9 — Mapping reads to {sum(1 for _ in open(otus_fasta) if _.startswith('>'))} OTUs")
    sample_names, otu_counts, mapped_totals = map_reads_to_otus(
        vsearch, samples_trimmed, otus_fasta, out_dir, args.otu_similarity, args.threads)

    if not otu_counts:
        log("[ERROR] No reads mapped to any OTU. Check db_path and input FASTQ files.")
        sys.exit(1)

    # ── 8. Taxonomy assignment ──────────────────────────────────────────────────
    progress(76, "Step 8/9 — Assigning taxonomy (vsearch vs. reference)")
    otu_tax = assign_taxonomy(vsearch, otus_fasta, args.db_path, out_dir,
                              args.tax_identity, args.threads)

    # ── 9. Write tables + viz ───────────────────────────────────────────────────
    progress(86, "Step 9/9 — Writing tables & generating plots")
    sorted_otus = write_tables(sample_names, otu_counts, otu_tax, out_dir)
    read_rows = write_read_tracking(samples_raw, samples_qc, samples_trimmed, mapped_totals, out_dir)
    write_summary(out_dir, sample_names, len(sorted_otus), read_rows, args)
    run_viz(out_dir, args.metadata, args.topN, args.threads, args.marker)

    # Cleanup large intermediates (keep otus.fasta + tables; drop pooled/uniques)
    for p in (pooled_fasta, uniques, out_dir / "uniques.sorted.fasta",
              out_dir / "uniques.nonchimeras.fasta"):
        try:
            if p.exists(): p.unlink()
        except Exception:
            pass

    progress(100, "Complete")
    log("=" * 50)
    log("  qiime2_vsearch_pipeline.py completed successfully")
    log("=" * 50)

if __name__ == "__main__":
    main()
