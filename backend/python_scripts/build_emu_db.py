#!/usr/bin/env python3
"""
build_emu_db.py — Build an Emu-compatible database from a SILVA DADA2 trainset FASTA.

Works with SILVA DADA2 format where FASTA headers ARE the full taxonomy string, e.g.:
  >Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio;

Tested with:
  silva_nr99_v138.1_wSpecies_train_set.fa.gz
  silva_nr99_v138.2_toGenus_trainset.fa.gz

Usage:
  python3 build_emu_db.py <silva_fasta.fa.gz> <output_dir>

Output <output_dir>/silva_db/ will contain:
  species_taxid.fasta   — sequences with seq_N IDs (Emu format)
  taxonomy.tsv          — taxid→lineage map (Emu format)

After this script, run (with emu conda env active):
  emu build-database \\
    --sequences  <output_dir>/sequences.fasta \\
    --seq2tax    <output_dir>/seq2taxid.tsv \\
    --taxonomy-list <output_dir>/taxonomy_list.tsv \\
    --db-name    silva_db \\
    --output-dir <output_dir>
"""

import sys
import gzip
import pathlib
import textwrap


def open_fasta(path):
    """Open plain or gzipped FASTA."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_fasta(path):
    """Yield (header, seq) pairs."""
    header, seq_parts = None, []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def normalise_lineage(header):
    """
    Convert SILVA DADA2 header → clean semicolon-delimited lineage.
    Strips trailing semicolons / whitespace; replaces spaces with underscores.
    """
    lin = header.strip().rstrip(";").strip()
    # Replace spaces inside taxon names
    parts = [p.strip().replace(" ", "_") for p in lin.split(";")]
    return ";".join(p for p in parts if p)


def build_intermediate_files(silva_fasta, out_dir):
    """
    Parse SILVA DADA2 FASTA and write:
      sequences.fasta   — seq_N identifiers
      seq2taxid.tsv     — seq_id TAB taxid
      taxonomy_list.tsv — taxid TAB lineage
    Returns (n_seqs, n_taxa).
    """
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lineage_to_taxid = {}          # lineage_str -> int taxid
    next_taxid = 1

    seq_records = []               # list of (seq_id, taxid, sequence)

    print(f"Parsing {silva_fasta} …", flush=True)
    for i, (header, seq) in enumerate(parse_fasta(silva_fasta)):
        if i % 50_000 == 0 and i > 0:
            print(f"  {i:,} sequences processed…", flush=True)

        lineage = normalise_lineage(header)
        if not lineage:
            continue  # skip blank headers

        if lineage not in lineage_to_taxid:
            lineage_to_taxid[lineage] = next_taxid
            next_taxid += 1

        taxid = lineage_to_taxid[lineage]
        seq_id = f"seq_{i + 1}"
        seq_records.append((seq_id, taxid, seq))

    n_seqs = len(seq_records)
    n_taxa = len(lineage_to_taxid)
    print(f"  Total: {n_seqs:,} sequences, {n_taxa:,} unique taxa", flush=True)

    # --- sequences.fasta ---
    seqs_path = out_dir / "sequences.fasta"
    print(f"Writing {seqs_path} …", flush=True)
    with open(seqs_path, "w") as fh:
        for seq_id, taxid, seq in seq_records:
            fh.write(f">{seq_id}\n")
            for chunk in textwrap.wrap(seq, 80):
                fh.write(chunk + "\n")

    # --- seq2taxid.tsv --- (seq_id TAB taxid — Emu expects this order)
    s2t_path = out_dir / "seq2taxid.tsv"
    print(f"Writing {s2t_path} …", flush=True)
    with open(s2t_path, "w") as fh:
        for seq_id, taxid, _ in seq_records:
            fh.write(f"{seq_id}\t{taxid}\n")

    # --- taxonomy_list.tsv --- (taxid TAB lineage)
    tax_path = out_dir / "taxonomy_list.tsv"
    print(f"Writing {tax_path} …", flush=True)
    with open(tax_path, "w") as fh:
        for lineage, taxid in sorted(lineage_to_taxid.items(), key=lambda x: x[1]):
            fh.write(f"{taxid}\t{lineage}\n")

    return n_seqs, n_taxa


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    silva_fasta = sys.argv[1]
    out_dir = pathlib.Path(sys.argv[2])

    if not pathlib.Path(silva_fasta).exists():
        print(f"ERROR: FASTA not found: {silva_fasta}")
        sys.exit(1)

    n_seqs, n_taxa = build_intermediate_files(silva_fasta, out_dir)

    db_dir = out_dir / "silva_db"
    print("\nIntermediate files ready. Now run:\n")
    print(
        f"  conda run --no-capture-output -n emu \\\n"
        f"    emu build-database \\\n"
        f"    --sequences  {out_dir}/sequences.fasta \\\n"
        f"    --seq2tax    {out_dir}/seq2taxid.tsv \\\n"
        f"    --taxonomy-list {out_dir}/taxonomy_list.tsv \\\n"
        f"    --db-name    silva_db \\\n"
        f"    --output-dir {out_dir}\n"
    )
    print(f"Expected output: {db_dir}/species_taxid.fasta  +  {db_dir}/taxonomy.tsv")
    print(f"  ({n_seqs:,} sequences, {n_taxa:,} taxa)")


if __name__ == "__main__":
    main()
