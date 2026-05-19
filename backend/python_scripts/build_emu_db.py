#!/usr/bin/env python3
"""
build_emu_db.py — Build an Emu-compatible database from a SILVA DADA2 trainset FASTA.

Supports full-length AND region-specific databases (e.g. V7–V8 for ONT partial reads).
When --primer-f / --primer-r are provided, the script extracts the amplicon region
from each SILVA sequence before building the database — matching what Biostrings
matchPattern() does, but in pure Python without R.

Works with SILVA DADA2 format where FASTA headers ARE the full taxonomy string, e.g.:
  >Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio;

Tested with:
  silva_nr99_v138.1_train_set.fa.gz
  silva_nr99_v138.1_wSpecies_train_set.fa.gz

Usage (full-length, default):
  python3 build_emu_db.py --silva-db silva.fa.gz --output-dir ./emu_silva

Usage (region-specific, e.g. V7–V8):
  python3 build_emu_db.py \\
    --silva-db silva.fa.gz \\
    --output-dir ./emu_v7v8 \\
    --primer-f AACMGGATTAGATACCCKG \\
    --primer-r ACGTCATCCCCACCTTCC \\
    --min-len 250 --max-len 450

After this script, run (with emu conda env active):
  emu build-database \\
    --sequences  <output_dir>/sequences.fasta \\
    --seq2tax    <output_dir>/seq2taxid.tsv \\
    --taxonomy-list <output_dir>/taxonomy_list.tsv \\
    --db-name    silva_db \\
    --output-dir <output_dir>
"""

import sys
import re
import gzip
import pathlib
import argparse
import textwrap

# ── IUPAC degenerate base expansion ──────────────────────────────────────────
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]",  "Y": "[CT]",  "S": "[GC]",  "W": "[AT]",
    "K": "[GT]",  "M": "[AC]",  "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                            "TGCAARYSWKMBDHVNtgcaaryswkmbdhvn")


def iupac_to_regex(primer: str) -> re.Pattern:
    """Convert an IUPAC primer string to a compiled regex (case-insensitive)."""
    pattern = "".join(IUPAC.get(b.upper(), b.upper()) for b in primer)
    return re.compile(pattern, re.IGNORECASE)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def find_primer(seq: str, primer_re: re.Pattern):
    """Return the first match object, or None."""
    return primer_re.search(seq)


# ── FASTA helpers ─────────────────────────────────────────────────────────────
def open_fasta(path):
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


def normalise_lineage(header: str) -> str:
    """Convert SILVA DADA2 header → clean semicolon-delimited lineage."""
    lin = header.strip().rstrip(";").strip()
    parts = [p.strip().replace(" ", "_") for p in lin.split(";")]
    return ";".join(p for p in parts if p)


# ── Region extraction ─────────────────────────────────────────────────────────
def extract_amplicon(seq: str, fwd_re: re.Pattern, rev_rc_re: re.Pattern,
                     min_len: int, max_len: int) -> str | None:
    """
    Find amplicon between forward primer and reverse-complement of reverse primer.
    Returns the region BETWEEN primers (primers trimmed), or None if not found.
    Tries both orientations (handles reads on either strand).
    """
    for s in (seq, reverse_complement(seq)):
        fwd_m = find_primer(s, fwd_re)
        rev_m = find_primer(s, rev_rc_re)
        if fwd_m and rev_m:
            amp_start = fwd_m.end()          # after forward primer
            amp_end   = rev_m.start()        # before rev-comp of reverse primer
            if amp_end > amp_start:
                amplicon = s[amp_start:amp_end]
                if min_len <= len(amplicon) <= max_len:
                    return amplicon
    return None


# ── Core builder ──────────────────────────────────────────────────────────────
def build_intermediate_files(silva_fasta: str, out_dir: pathlib.Path,
                              fwd_re=None, rev_rc_re=None,
                              min_len: int = 0, max_len: int = 99999):
    """
    Parse SILVA DADA2 FASTA and write Emu intermediate files:
      sequences.fasta   — seq_N identifiers
      seq2taxid.tsv     — seq_id TAB taxid
      taxonomy_list.tsv — taxid TAB lineage

    If fwd_re / rev_rc_re are provided, extract amplicon region first.
    Returns (n_seqs, n_taxa).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    region_mode = fwd_re is not None and rev_rc_re is not None
    if region_mode:
        print(f"Mode: REGION-SPECIFIC extraction (min={min_len} bp, max={max_len} bp)", flush=True)
    else:
        print("Mode: FULL-LENGTH (no primer extraction)", flush=True)

    lineage_to_taxid: dict[str, int] = {}
    next_taxid = 1
    seq_records = []

    n_total = 0
    n_skipped = 0

    print(f"Parsing {silva_fasta} …", flush=True)
    for i, (header, seq) in enumerate(parse_fasta(silva_fasta)):
        if i % 50_000 == 0 and i > 0:
            print(f"  {i:,} sequences scanned, {len(seq_records):,} kept…", flush=True)
        n_total += 1

        lineage = normalise_lineage(header)
        if not lineage:
            n_skipped += 1
            continue

        # ── Region extraction ──
        if region_mode:
            final_seq = extract_amplicon(seq, fwd_re, rev_rc_re, min_len, max_len)
            if final_seq is None:
                n_skipped += 1
                continue
        else:
            final_seq = seq.upper().replace("U", "T")

        if lineage not in lineage_to_taxid:
            lineage_to_taxid[lineage] = next_taxid
            next_taxid += 1

        taxid  = lineage_to_taxid[lineage]
        seq_id = f"seq_{i + 1}"
        seq_records.append((seq_id, taxid, final_seq))

    n_seqs = len(seq_records)
    n_taxa = len(lineage_to_taxid)
    print(f"  Scanned: {n_total:,}  |  Kept: {n_seqs:,}  |  Skipped: {n_skipped:,}", flush=True)
    print(f"  Unique taxa: {n_taxa:,}", flush=True)

    if n_seqs == 0:
        raise RuntimeError(
            "No sequences passed the primer filter!\n"
            "Check that your primers match SILVA sequences and adjust --min-len / --max-len."
        )

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

    # ── taxonomy.tsv — Emu format ─────────────────────────────
    # Emu reads this with pandas and expects a 'tax_id' column header,
    # followed by separate columns for each taxonomy level.
    # Levels: superkingdom, phylum, class, order, family, genus, species
    TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    tax_path = out_dir / "taxonomy.tsv"
    print(f"Writing {tax_path} …", flush=True)
    with open(tax_path, "w") as fh:
        fh.write("tax_id\t" + "\t".join(TAX_LEVELS) + "\n")
        for lineage, taxid in sorted(lineage_to_taxid.items(), key=lambda x: x[1]):
            parts = lineage.split(";")
            parts += [""] * (len(TAX_LEVELS) - len(parts))   # pad if short
            parts = parts[:len(TAX_LEVELS)]                   # trim if too long
            fh.write(f"{taxid}\t" + "\t".join(parts) + "\n")
    # Keep taxonomy_list.tsv (original format) for reference
    alias_path = out_dir / "taxonomy_list.tsv"
    if not alias_path.exists():
        with open(alias_path, "w") as fh:
            for lineage, taxid in sorted(lineage_to_taxid.items(), key=lambda x: x[1]):
                fh.write(f"{taxid}\t{lineage}\n")

    return n_seqs, n_taxa


# ── CLI ───────────────────────────────────────────────────────────────────────

# Common primer presets for quick reference
PRIMER_PRESETS = {
    "V7-V8":  {"f": "AACMGGATTAGATACCCKG",   "r": "ACGTCATCCCCACCTTCC",   "min": 250, "max": 450},
    "V1-V9":  {"f": "AGRGTTYGATYMTGGCTCAG",   "r": "RGYTACCTTGTTACGACTT",  "min": 1000,"max": 1800},
    "V3-V4":  {"f": "CCTACGGGNGGCWGCAG",      "r": "GACTACHVGGGTATCTAATCC","min": 350, "max": 500},
    "V4":     {"f": "GTGYCAGCMGCCGCGGTAA",    "r": "GGACTACNVGGGTWTCTAAT", "min": 200, "max": 280},
    "V4-V5":  {"f": "GTGYCAGCMGCCGCGGTAA",    "r": "CCGYCAATTYMTTTRAGTTT", "min": 300, "max": 480},
}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--silva-db",   required=True,
                        help="SILVA DADA2 trainset FASTA (plain or .gz)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for intermediate files")

    # Region extraction options
    grp = parser.add_argument_group("Region extraction (optional)")
    grp.add_argument("--primer-f",  default="",
                     help="Forward primer sequence (IUPAC). E.g. AACMGGATTAGATACCCKG")
    grp.add_argument("--primer-r",  default="",
                     help="Reverse primer sequence (IUPAC). E.g. ACGTCATCCCCACCTTCC")
    grp.add_argument("--region",    choices=list(PRIMER_PRESETS.keys()),
                     help="Preset region (sets primers + lengths). Overrides --primer-f/r/min/max.")
    grp.add_argument("--min-len",   type=int, default=None,
                     help="Min amplicon length after primer trim (default: 250 for region mode)")
    grp.add_argument("--max-len",   type=int, default=None,
                     help="Max amplicon length after primer trim (default: 450 for region mode)")

    # Legacy positional args (backward compat with old 2-arg form)
    parser.add_argument("_silva_pos", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("_outdir_pos", nargs="?", help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Backward compatibility: if called as  build_emu_db.py <silva> <outdir>
    if args._silva_pos and not args.silva_db:
        args.silva_db = args._silva_pos
    if args._outdir_pos and not args.output_dir:
        args.output_dir = args._outdir_pos

    # Apply preset if requested
    if args.region:
        preset = PRIMER_PRESETS[args.region]
        if not args.primer_f: args.primer_f = preset["f"]
        if not args.primer_r: args.primer_r = preset["r"]
        if args.min_len is None: args.min_len = preset["min"]
        if args.max_len is None: args.max_len = preset["max"]
        print(f"Using preset {args.region}: F={args.primer_f}  R={args.primer_r}")

    # Defaults for min/max
    if args.min_len is None:
        args.min_len = 250 if args.primer_f else 0
    if args.max_len is None:
        args.max_len = 450 if args.primer_f else 99999

    silva_path = pathlib.Path(args.silva_db)
    out_dir    = pathlib.Path(args.output_dir)

    if not silva_path.exists():
        print(f"ERROR: SILVA FASTA not found: {silva_path}")
        sys.exit(1)

    # Build regex objects
    fwd_re = rev_rc_re = None
    if args.primer_f and args.primer_r:
        fwd_re    = iupac_to_regex(args.primer_f)
        rev_rc_re = iupac_to_regex(reverse_complement(args.primer_r))
        print(f"Forward primer : {args.primer_f}")
        print(f"Reverse primer : {args.primer_r}  (rev-comp for search: {reverse_complement(args.primer_r)})")
        print(f"Length filter  : {args.min_len}–{args.max_len} bp")
    elif args.primer_f or args.primer_r:
        print("WARNING: Both --primer-f and --primer-r must be provided to enable extraction. Running full-length.")

    n_seqs, n_taxa = build_intermediate_files(
        silva_fasta = str(silva_path),
        out_dir     = out_dir,
        fwd_re      = fwd_re,
        rev_rc_re   = rev_rc_re,
        min_len     = args.min_len,
        max_len     = args.max_len,
    )

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
    region_tag = f"  Region: {args.region or 'custom primers'}" if fwd_re else "  Region: full-length"
    print(f"Expected output: {out_dir}/silva_db/")
    print(f"  {n_seqs:,} sequences, {n_taxa:,} taxa")
    print(region_tag)


if __name__ == "__main__":
    main()
