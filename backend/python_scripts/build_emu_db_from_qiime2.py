#!/usr/bin/env python3
"""
build_emu_db_from_qiime2.py — Build an Emu-compatible database from QIIME2 artifacts.

WHY THIS SCRIPT EXISTS (vs build_emu_db.py):
  build_emu_db.py sources sequences from the SILVA *DADA2-formatted* trainset
  (silva_nr99_v138.1_train_set.fa.gz). That file's FASTA headers are capped at
  6 taxonomy ranks (Kingdom→Genus) by design — DADA2's assignTaxonomy() only
  ever resolves to Genus, with a *separate* addSpecies() step for species. As a
  result, any Emu database built from it has an always-empty "species" column,
  because the source data never had species-level names to begin with.

  build_emu_db.py also extracts the target region (e.g. V7-V8) using a plain
  regex primer match with ZERO mismatch tolerance beyond IUPAC degeneracy. Real
  "universal" 16S primers have known binding-site mismatches against some
  lineages (Cyanobacteria is a classic, well-documented example) — a 0-mismatch
  regex silently drops every reference sequence with even a single mismatch,
  systematically under-representing those lineages in the database.

  This script sidesteps BOTH problems by building from QIIME2 artifacts that
  already contain the full SILVA taxonomy (all 7 ranks, real species names)
  and an already region-extracted sequence set (extracted via QIIME2's own
  primer-trimming, done once by whoever ran the QIIME2 pipeline — no re-extraction
  needed here).

INPUT (already unzipped from .qza — .qza files are just zip archives):
  --taxonomy-tsv   path to <uuid>/data/taxonomy.tsv from a FeatureData[Taxonomy] .qza
                   (columns: "Feature ID", "Taxon", ...)
  --sequences-fasta path to <uuid>/data/dna-sequences.fasta from a
                   FeatureData[Sequence] .qza (already region-extracted)

To unzip a .qza:
  unzip -o -q your_taxonomy.qza -d tax_extracted
  unzip -o -q your_ref_seqs.qza -d seqs_extracted
  # then point this script at:
  #   tax_extracted/*/data/taxonomy.tsv
  #   seqs_extracted/*/data/dna-sequences.fasta

OUTPUT (in --output-dir):
  sequences.fasta   — seq_N identifiers
  seq2taxid.tsv     — seq_id TAB taxid
  taxonomy.tsv      — Emu format: tax_id, superkingdom, phylum, class, order,
                       family, genus, species (wide, header'd — NOT the 2-column
                       alias format)

After this script, run (with emu conda env active), from INSIDE output-dir:
  cd <output_dir>
  conda run -n emu emu build-database \\
    --sequences  sequences.fasta \\
    --seq2tax    seq2taxid.tsv \\
    --taxonomy-list taxonomy.tsv \\
    <db_name>
"""

import sys
import pathlib
import argparse
import textwrap


def parse_fasta(path):
    header, seq_parts = None, []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]   # Feature ID = first whitespace-delimited token
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def load_taxonomy(path) -> dict[str, str]:
    """Read a QIIME2 taxonomy.tsv → {Feature ID: Taxon lineage string}."""
    tax = {}
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # QIIME2 exports "Feature ID" and "Taxon" as the first two columns
        fid_idx = 0
        taxon_idx = header.index("Taxon") if "Taxon" in header else 1
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(fid_idx, taxon_idx):
                continue
            fid = parts[fid_idx].strip()
            taxon = parts[taxon_idx].strip()
            if fid and taxon:
                tax[fid] = taxon
    return tax


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxonomy-tsv", required=True, help="Path to extracted taxonomy.tsv from a QIIME2 taxonomy .qza")
    ap.add_argument("--sequences-fasta", required=True, help="Path to extracted dna-sequences.fasta from a QIIME2 sequence .qza")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading taxonomy from {args.taxonomy_tsv} …", flush=True)
    tax_map = load_taxonomy(args.taxonomy_tsv)
    print(f"  {len(tax_map):,} Feature ID → Taxon entries loaded", flush=True)

    print(f"Parsing {args.sequences_fasta} …", flush=True)
    TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    lineage_to_taxid: dict[str, int] = {}
    next_taxid = 1
    seq_records = []
    n_total = 0
    n_no_taxonomy = 0

    for i, (fid, seq) in enumerate(parse_fasta(args.sequences_fasta)):
        if i % 50_000 == 0 and i > 0:
            print(f"  {i:,} sequences scanned, {len(seq_records):,} kept…", flush=True)
        n_total += 1
        lineage = tax_map.get(fid)
        if not lineage:
            n_no_taxonomy += 1
            continue
        if lineage not in lineage_to_taxid:
            lineage_to_taxid[lineage] = next_taxid
            next_taxid += 1
        taxid = lineage_to_taxid[lineage]
        seq_id = f"seq_{i + 1}"
        seq_records.append((seq_id, taxid, seq.upper().replace("U", "T")))

    n_seqs = len(seq_records)
    n_taxa = len(lineage_to_taxid)
    print(f"  Scanned: {n_total:,}  |  Kept: {n_seqs:,}  |  No taxonomy match: {n_no_taxonomy:,}", flush=True)
    print(f"  Unique taxa: {n_taxa:,}", flush=True)

    if n_seqs == 0:
        raise RuntimeError("No sequences matched a taxonomy entry — check Feature ID formats line up "
                            "between --sequences-fasta and --taxonomy-tsv.")

    # ── sequences.fasta ──
    seqs_path = out_dir / "sequences.fasta"
    print(f"Writing {seqs_path} …", flush=True)
    with open(seqs_path, "w") as fh:
        for seq_id, taxid, seq in seq_records:
            fh.write(f">{seq_id}\n")
            for chunk in textwrap.wrap(seq, 80):
                fh.write(chunk + "\n")

    # ── seq2taxid.tsv ──
    s2t_path = out_dir / "seq2taxid.tsv"
    print(f"Writing {s2t_path} …", flush=True)
    with open(s2t_path, "w") as fh:
        for seq_id, taxid, _ in seq_records:
            fh.write(f"{seq_id}\t{taxid}\n")

    # ── taxonomy.tsv — Emu format (wide, header'd, one column per rank) ──
    tax_path = out_dir / "taxonomy.tsv"
    print(f"Writing {tax_path} …", flush=True)
    n_with_species = 0
    with open(tax_path, "w") as fh:
        fh.write("tax_id\t" + "\t".join(TAX_LEVELS) + "\n")
        for lineage, taxid in sorted(lineage_to_taxid.items(), key=lambda x: x[1]):
            parts = [p.strip() for p in lineage.split(";")]
            parts += [""] * (len(TAX_LEVELS) - len(parts))
            parts = parts[:len(TAX_LEVELS)]
            if parts[6]:
                n_with_species += 1
            fh.write(f"{taxid}\t" + "\t".join(parts) + "\n")

    print(f"  {n_with_species:,} / {n_taxa:,} unique taxa have a species-level name "
          f"({n_with_species/n_taxa*100:.1f}%)", flush=True)

    print("\nIntermediate files ready. Now run:\n")
    print(
        f"  cd {out_dir} && \\\n"
        f"  conda run --no-capture-output -n emu \\\n"
        f"    emu build-database \\\n"
        f"    --sequences  sequences.fasta \\\n"
        f"    --seq2tax    seq2taxid.tsv \\\n"
        f"    --taxonomy-list taxonomy.tsv \\\n"
        f"    silva_qiime2_v7v8\n"
    )
    print(f"Expected output: {out_dir}/silva_qiime2_v7v8/")
    print(f"  {n_seqs:,} sequences, {n_taxa:,} taxa")


if __name__ == "__main__":
    main()
