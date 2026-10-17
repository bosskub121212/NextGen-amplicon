#!/usr/bin/env python3
"""Simulated amplicon reads from a mock community with a known composition.

The reads are made from real reference sequences (refs_16s.fasta, refs_its.fasta,
refs_coi.fasta — NCBI RefSeq / GenBank, see each header), cut out with the same
primers the pipeline trims, sampled in known proportions and given realistic
substitution errors with matching quality scores. Because the truth is known
exactly, a pipeline that shifts ranks, swaps samples, drops a taxon or invents one
fails the check in check_mock.py — which is how the rank-shift, PICRUSt2 and PacBio
orientation bugs of v2.14 would have been caught.

    python3 simulate.py --refs refs_16s.fasta --composition zymo_16s \\
        --primer_f CCTACGGGNGGCWGCAG --primer_r GACTACHVGGGTATCTAATCC \\
        --platform illumina --out /tmp/mock16s

Platforms: illumina (paired 2x250, quality falling along the read),
           illumina_single (R1 only), pacbio (HiFi, full amplicon, half the reads
           reverse-complemented), ont (single reads, ~1% errors, both strands).
Standard library only; deterministic for a given --seed.
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Theoretical composition of the ZymoBIOMICS Microbial Community Standard (D6300),
# "16S only" column (copy-number weighted), Zymo protocol Table 1.
COMPOSITIONS = {
    "zymo_16s": {"Pseudomonas_aeruginosa": 4.2, "Escherichia_coli": 10.1, "Salmonella_enterica": 10.4,
                 "Limosilactobacillus_fermentum": 18.4, "Enterococcus_faecalis": 9.9,
                 "Staphylococcus_aureus": 15.5, "Listeria_monocytogenes": 14.1, "Bacillus_subtilis": 17.4},
    "fungi_its": {"Saccharomyces_cerevisiae": 30, "Cryptococcus_neoformans": 10, "Aspergillus_niger": 25,
                  "Candida_albicans": 20, "Penicillium_chrysogenum": 15},
    "insects_coi": {"Drosophila_melanogaster": 30, "Apis_mellifera": 25, "Musca_domestica": 20,
                    "Bombyx_mori": 15, "Aedes_aegypti": 10},
}

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", "I": "ACGT"}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def rc(s: str) -> str:
    return s.upper().translate(COMP)[::-1]


def read_fasta(p: Path) -> dict[str, str]:
    out, name = {}, None
    for ln in p.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += ln.strip().upper()
    return out


def find_primer(seq: str, primer: str, max_mm: int = 5) -> int:
    """Start of the best degenerate match (fewest mismatches)."""
    best, pos = max_mm + 1, -1
    L = len(primer)
    for i in range(len(seq) - L + 1):
        mm = 0
        for j, c in enumerate(primer):
            if seq[i + j] not in IUPAC.get(c, c):
                mm += 1
                if mm > max_mm:
                    break
        if mm < best:
            best, pos = mm, i
            if mm == 0:
                break
    return pos


def resolve(primer: str, rng: random.Random) -> str:
    """One concrete primer sequence out of a degenerate one."""
    return "".join(rng.choice(IUPAC.get(c, c)) for c in primer.upper())


def amplicon(seq: str, pf: str, pr: str, rng: random.Random) -> str:
    """Primer-to-primer product. A primer site missing from a partial reference is
    added, so the read still starts with a primer the pipeline can trim."""
    i = find_primer(seq, pf) if pf else 0
    j = find_primer(seq, rc(pr)) if pr else -1
    start = i + len(pf) if i >= 0 else 0
    end = j if j > start else len(seq)
    insert = seq[start:end]
    return (resolve(pf, rng) if pf else "") + insert + (rc(resolve(pr, rng)) if pr else "")


def qual_profile(n: int, platform: str, rng: random.Random) -> list[int]:
    if platform == "pacbio":
        return [min(93, max(20, int(rng.gauss(40, 3)))) for _ in range(n)]
    if platform == "ont":
        return [max(5, int(rng.gauss(20, 4))) for _ in range(n)]
    # Illumina: high at the start, falling towards the 3' end
    out = []
    for k in range(n):
        mu = 37 - 12 * (k / max(n - 1, 1)) ** 2.2
        out.append(max(2, min(41, int(rng.gauss(mu, 2.5)))))
    return out


def mutate(seq: str, quals: list[int], rng: random.Random) -> str:
    s = list(seq)
    for k, q in enumerate(quals):
        if rng.random() < 10 ** (-q / 10):
            s[k] = rng.choice([b for b in "ACGT" if b != s[k]])
    return "".join(s)


def fq(name, seq, quals):
    return f"@{name}\n{seq}\n+\n{''.join(chr(q + 33) for q in quals)}\n"


def multinomial(n: int, probs: list[float], rng: random.Random) -> list[int]:
    counts = [0] * len(probs)
    cum = []
    acc = 0.0
    for p in probs:
        acc += p
        cum.append(acc)
    for _ in range(n):
        r = rng.random() * acc
        for k, c in enumerate(cum):
            if r <= c:
                counts[k] += 1
                break
    return counts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refs", required=True)
    ap.add_argument("--composition", required=True,
                    help=f"one of {', '.join(COMPOSITIONS)} or a JSON file {{name: percent}}")
    ap.add_argument("--primer_f", default="")
    ap.add_argument("--primer_r", default="")
    ap.add_argument("--platform", default="illumina",
                    choices=["illumina", "illumina_single", "pacbio", "ont"])
    ap.add_argument("--read_len", type=int, default=250)
    ap.add_argument("--reads", type=int, default=4000, help="reads (pairs) per sample")
    ap.add_argument("--samples", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    rng = random.Random(a.seed)
    refs = read_fasta(Path(a.refs) if Path(a.refs).exists() else HERE / a.refs)
    comp = (COMPOSITIONS[a.composition] if a.composition in COMPOSITIONS
            else json.loads(Path(a.composition).read_text()))
    missing = [k for k in comp if k not in refs]
    if missing:
        raise SystemExit(f"no reference for {missing}")
    names = list(comp)
    tot = sum(comp.values())
    probs = [comp[k] / tot for k in names]
    # Reference records can carry ambiguity codes (N, Y, R ...). A real organism
    # has one base there, and filterAndTrim(maxN = 0) drops any read with an N —
    # which silently removed E. coli and Listeria from the first test run.
    for k in names:
        refs[k] = "".join(rng.choice(IUPAC.get(c, "ACGT")) for c in refs[k])
    amps = {k: amplicon(refs[k], a.primer_f, a.primer_r, rng) for k in names}

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    truth = {"composition": {k: round(100 * p, 3) for k, p in zip(names, probs)},
             "amplicon_bp": {k: len(v) for k, v in amps.items()},
             "platform": a.platform, "primer_f": a.primer_f, "primer_r": a.primer_r,
             "samples": {}}
    for si in range(a.samples):
        sname = f"mock{si + 1}"
        counts = multinomial(a.reads, probs, rng)
        truth["samples"][sname] = dict(zip(names, counts))
        r1 = gzip.open(out / (f"{sname}_R1.fastq.gz" if a.platform.startswith("illumina")
                              else f"{sname}.fastq.gz"), "wt")
        r2 = gzip.open(out / f"{sname}_R2.fastq.gz", "wt") if a.platform == "illumina" else None
        n = 0
        for k, c in zip(names, counts):
            amp = amps[k]
            for _ in range(c):
                n += 1
                rid = f"{sname}.{n} {k}"
                if a.platform.startswith("illumina"):
                    L = min(a.read_len, len(amp))
                    q1 = qual_profile(L, "illumina", rng)
                    r1.write(fq(rid + " 1:N:0", mutate(amp[:L], q1, rng), q1))
                    if r2:
                        q2 = qual_profile(L, "illumina", rng)
                        q2 = [max(2, q - 2) for q in q2]          # R2 is the worse read
                        r2.write(fq(rid + " 2:N:0", mutate(rc(amp)[:L], q2, rng), q2))
                else:
                    s = amp if rng.random() < 0.5 else rc(amp)  # both orientations
                    q = qual_profile(len(s), a.platform, rng)
                    r1.write(fq(rid, mutate(s, q, rng), q))
        r1.close()
        if r2:
            r2.close()
    (out / "truth.json").write_text(json.dumps(truth, indent=2))
    print(f"wrote {a.samples} samples x {a.reads} reads to {out}")
    for k in names:
        print(f"  {k:32s} {truth['composition'][k]:6.2f}%  amplicon {len(amps[k])} bp")


if __name__ == "__main__":
    main()
