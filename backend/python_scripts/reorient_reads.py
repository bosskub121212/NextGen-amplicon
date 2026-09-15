#!/usr/bin/env python3
"""
reorient_reads.py — detect and fix mixed R1/R2 orientation in paired-end
amplicon FASTQ files.

WHY THIS EXISTS
---------------
In a normal paired-end amplicon library every read pair is sequenced in the
same orientation: R1 starts at the forward primer, R2 starts at the reverse
primer. Some library preps (and some sequencing runs) produce a library where
a substantial fraction of the molecules are loaded the other way round, so for
those pairs R1 starts at the REVERSE primer and R2 at the FORWARD primer.

DADA2 cannot merge those "flipped" pairs. mergePairs() expects R1 and
revcomp(R2) to overlap in the middle of the amplicon; for a flipped pair the
two mates are biologically misaligned, so the pair is silently dropped at the
merge step no matter what truncLen is used. cutadapt's paired-end mode also
searches a FIXED orientation (-g primer_f on R1, -G primer_r on R2), so it
leaves flipped pairs untrimmed rather than fixing them.

The symptom is a merge rate far below expectation (e.g. 25% instead of 70%+)
with no error message anywhere.

WHAT THIS DOES
--------------
Streams both mates in lockstep and, for every pair, decides its orientation by
looking for the forward/reverse primer near the start of each mate (IUPAC
degenerate codes supported, mismatches tolerated):

  R1 has primer_f, R2 has primer_r  -> already correct, copy through
  R1 has primer_r, R2 has primer_f  -> FLIPPED, swap the two mates
  only one mate matched             -> decided by that mate alone
  neither mate matched              -> ambiguous (see --ambiguous)

Swapping the mates is enough: after the swap the pair is in the standard
orientation, so the normal cutadapt call and mergePairs() both work.

USAGE
-----
  # report only, no files written (fast — samples the first N pairs)
  python3 reorient_reads.py --check \\
      --r1 S1_1.fq.gz --r2 S1_2.fq.gz \\
      --primer-f CCTACGGGNGGCWGCAG --primer-r GGACTACNVGGGTWTCTAAT

  # write reoriented files
  python3 reorient_reads.py \\
      --r1 S1_1.fq.gz --r2 S1_2.fq.gz --outdir ./reoriented \\
      --primer-f CCTACGGGNGGCWGCAG --primer-r GGACTACNVGGGTWTCTAAT

Both modes print a JSON report on stdout.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import re
import sys
from pathlib import Path

# ── IUPAC degenerate nucleotide codes ────────────────────────────────────────
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}

_COMP = str.maketrans("ACGTUacgtuRYSWKMBDHVNrysWkmbdhvn",
                      "TGCAAtgcaaYRSWMKVHDBNyrsWmkvhdbn")


def revcomp(seq: str) -> str:
    """Reverse complement, IUPAC aware."""
    return seq.translate(_COMP)[::-1]


def iupac_regex(primer: str) -> re.Pattern:
    """Compile a primer with IUPAC degenerate codes into a regex."""
    parts = []
    for ch in primer.upper():
        allowed = IUPAC.get(ch)
        if allowed is None:
            # Unknown character — match anything rather than failing outright
            parts.append(".")
        elif len(allowed) == 1:
            parts.append(allowed)
        else:
            parts.append("[" + allowed + "]")
    return re.compile("".join(parts))


def _expand(primer: str) -> list[str]:
    """Per-position allowed-base sets, for the mismatch-tolerant scan."""
    return [IUPAC.get(ch, "ACGT") for ch in primer.upper()]


class PrimerMatcher:
    """
    Finds a primer near the start of a read.

    Fast path is a compiled IUPAC regex (C-level, catches the large majority of
    reads). Only reads that fail the fast path go through the slower
    mismatch-tolerant scan, which keeps the whole thing fast enough to stream
    millions of reads.
    """

    def __init__(self, primer: str, window: int = 30, max_mismatch: int = 2):
        self.primer = primer.upper()
        self.plen = len(self.primer)
        self.window = max(window, self.plen + 6)
        self.max_mismatch = max_mismatch
        self.rx = iupac_regex(self.primer)
        self.sets = _expand(self.primer)

    def match(self, seq: str) -> bool:
        head = seq[: self.window]
        if self.rx.search(head):
            return True
        if self.max_mismatch <= 0:
            return False
        last_start = len(head) - self.plen
        if last_start < 0:
            return False
        sets = self.sets
        limit = self.max_mismatch
        for start in range(last_start + 1):
            mism = 0
            ok = True
            for i in range(self.plen):
                if head[start + i] not in sets[i]:
                    mism += 1
                    if mism > limit:
                        ok = False
                        break
            if ok:
                return True
        return False


# ── FASTQ streaming ──────────────────────────────────────────────────────────
def open_fastq(path: str | Path, mode: str = "rt"):
    """Open a FASTQ file, transparently handling gzip."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fastq(handle):
    """Yield (header, seq, plus, qual) 4-line records."""
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not q:
            return
        yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")


def _out_name(src: Path, outdir: Path, suffix: str) -> Path:
    """Build the output filename, preserving .fq.gz / .fastq.gz correctly."""
    name = src.name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(ext):
            stem = name[: -len(ext)]
            # always write gzipped
            out_ext = ext if ext.endswith(".gz") else ext + ".gz"
            return outdir / f"{stem}{suffix}{out_ext}"
    return outdir / f"{name}{suffix}.fastq.gz"


# ── Core ─────────────────────────────────────────────────────────────────────
def analyze_and_fix(
    r1_path: str,
    r2_path: str,
    primer_f: str,
    primer_r: str,
    outdir: str | None = None,
    suffix: str = "_reoriented",
    check_only: bool = False,
    sample_size: int = 20000,
    window: int = 30,
    max_mismatch: int = 2,
    ambiguous: str = "keep",
) -> dict:
    """
    Analyze (and optionally rewrite) a read pair.

    ambiguous: what to do with pairs where neither mate matched a primer
        "keep" — write through unchanged (default; never loses data)
        "drop" — leave them out of the output entirely
    """
    mf = PrimerMatcher(primer_f, window, max_mismatch)
    mr = PrimerMatcher(primer_r, window, max_mismatch)

    stats = {
        "total_pairs": 0,
        "already_correct": 0,
        "flipped": 0,
        "ambiguous": 0,
        "written_pairs": 0,
        "dropped_pairs": 0,
    }

    out1 = out2 = None
    h_out1 = h_out2 = None
    if not check_only:
        if not outdir:
            raise ValueError("outdir is required unless check_only=True")
        outp = Path(outdir)
        outp.mkdir(parents=True, exist_ok=True)
        out1 = _out_name(Path(r1_path), outp, suffix)
        out2 = _out_name(Path(r2_path), outp, suffix)
        h_out1 = gzip.open(out1, "wt", compresslevel=6)
        h_out2 = gzip.open(out2, "wt", compresslevel=6)

    try:
        with open_fastq(r1_path) as h1, open_fastq(r2_path) as h2:
            it1, it2 = read_fastq(h1), read_fastq(h2)
            for rec1, rec2 in zip(it1, it2):
                stats["total_pairs"] += 1

                s1, s2 = rec1[1], rec2[1]
                f_in_1 = mf.match(s1)
                r_in_2 = mr.match(s2)
                r_in_1 = mr.match(s1)
                f_in_2 = mf.match(s2)

                # score each hypothesis; 2 = both mates agree, 1 = one mate
                fwd_score = int(f_in_1) + int(r_in_2)
                rev_score = int(r_in_1) + int(f_in_2)

                if fwd_score > rev_score:
                    orientation = "forward"
                elif rev_score > fwd_score:
                    orientation = "flipped"
                else:
                    orientation = "ambiguous"

                if orientation == "forward":
                    stats["already_correct"] += 1
                    a, b = rec1, rec2
                elif orientation == "flipped":
                    stats["flipped"] += 1
                    a, b = rec2, rec1          # swap the mates
                else:
                    stats["ambiguous"] += 1
                    a, b = rec1, rec2

                if check_only:
                    if stats["total_pairs"] >= sample_size > 0:
                        break
                    continue

                if orientation == "ambiguous" and ambiguous == "drop":
                    stats["dropped_pairs"] += 1
                    continue

                h_out1.write(f"{a[0]}\n{a[1]}\n{a[2]}\n{a[3]}\n")
                h_out2.write(f"{b[0]}\n{b[1]}\n{b[2]}\n{b[3]}\n")
                stats["written_pairs"] += 1
    finally:
        if h_out1:
            h_out1.close()
        if h_out2:
            h_out2.close()

    total = max(stats["total_pairs"], 1)
    stats["pct_already_correct"] = round(stats["already_correct"] / total * 100, 2)
    stats["pct_flipped"] = round(stats["flipped"] / total * 100, 2)
    stats["pct_ambiguous"] = round(stats["ambiguous"] / total * 100, 2)
    stats["sampled"] = bool(check_only and sample_size > 0)
    stats["needs_reorientation"] = stats["pct_flipped"] >= 5.0

    if not check_only:
        stats["out_r1"] = str(out1)
        stats["out_r2"] = str(out2)

    # Human-readable verdict
    if stats["pct_flipped"] >= 5.0:
        stats["verdict"] = (
            f"{stats['pct_flipped']:.1f}% of read pairs are in the opposite "
            f"orientation. These cannot be merged by DADA2 as-is — reorienting "
            f"is strongly recommended."
        )
    elif stats["pct_ambiguous"] >= 50.0:
        stats["verdict"] = (
            f"{stats['pct_ambiguous']:.1f}% of read pairs matched neither primer. "
            f"Check that the primer sequences are correct for this data."
        )
    else:
        stats["verdict"] = "Read orientation looks consistent — no reorientation needed."

    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Detect and fix mixed R1/R2 orientation in paired-end FASTQ files.")
    ap.add_argument("--r1", required=True, help="R1 FASTQ (.fq/.fastq, optionally .gz)")
    ap.add_argument("--r2", required=True, help="R2 FASTQ (.fq/.fastq, optionally .gz)")
    ap.add_argument("--primer-f", required=True, help="Forward primer (IUPAC codes OK)")
    ap.add_argument("--primer-r", required=True, help="Reverse primer (IUPAC codes OK)")
    ap.add_argument("--outdir", default="", help="Output directory (required unless --check)")
    ap.add_argument("--suffix", default="_reoriented", help="Suffix for output filenames")
    ap.add_argument("--check", action="store_true",
                    help="Report only — scan a subsample and write nothing")
    ap.add_argument("--sample-size", type=int, default=20000,
                    help="Pairs to scan in --check mode (0 = all)")
    ap.add_argument("--window", type=int, default=30,
                    help="How many leading bases to search for the primer")
    ap.add_argument("--max-mismatch", type=int, default=2,
                    help="Mismatches tolerated when matching a primer")
    ap.add_argument("--ambiguous", choices=["keep", "drop"], default="keep",
                    help="What to do with pairs matching neither primer")
    args = ap.parse_args(argv)

    if not args.check and not args.outdir:
        ap.error("--outdir is required unless --check is given")

    stats = analyze_and_fix(
        r1_path=args.r1,
        r2_path=args.r2,
        primer_f=args.primer_f,
        primer_r=args.primer_r,
        outdir=args.outdir or None,
        suffix=args.suffix,
        check_only=args.check,
        sample_size=args.sample_size,
        window=args.window,
        max_mismatch=args.max_mismatch,
        ambiguous=args.ambiguous,
    )
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
