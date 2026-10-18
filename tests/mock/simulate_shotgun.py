#!/usr/bin/env python3
"""Simulated shotgun metagenome reads from genomes in known proportions.

    python3 simulate_shotgun.py --genomes spec.json --platform illumina|ont --reads 200000
        --out DIR --sample NAME [--host host.fa --host_frac 0.2] [--seed 1]

spec.json: [{"name": "Escherichia coli", "taxid": 562, "fasta": "GCF_x.fna.gz", "cells": 10}, ...]
"cells" are relative cell numbers; reads are drawn in proportion to cells x
genome length (a bigger genome gives more reads per cell), so the truth file
gives both: cell share and read share.

Illumina: 2x150 pairs, insert ~350 bp (sd 50), substitutions rising from 0.1 %
to 1 % along the read. ONT: log-normal lengths (median ~4 kb), ~2 % errors
(subs / ins / del), both strands. Standard library only, deterministic.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import random
from pathlib import Path

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def load(fa):
    op = gzip.open if str(fa).endswith(".gz") else open
    seqs, cur = [], []
    with op(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).upper())
    return seqs


def pick(contigs, cum, rng, L):
    """A random window of length L from a genome, contig chosen by length."""
    for _ in range(50):
        x = rng.random() * cum[-1]
        i = next(k for k, c in enumerate(cum) if c >= x)
        c = contigs[i]
        if len(c) > L:
            s = rng.randrange(0, len(c) - L)
            return c[s:s + L]
    return None


def err_illumina(s, rng):
    n = len(s)
    out = []
    for i, ch in enumerate(s):
        p = 0.001 + 0.009 * i / max(1, n - 1)
        out.append(rng.choice([b for b in "ACGT" if b != ch]) if rng.random() < p else ch)
    q = "".join(chr(33 + max(2, int(38 - 18 * i / n))) for i in range(n))
    return "".join(out), q


def err_ont(s, rng, rate=0.02):
    out = []
    for ch in s:
        x = rng.random()
        if x < rate * 0.35:
            continue
        if x < rate * 0.8:
            out.append(rng.choice("ACGT"))
        else:
            out.append(ch)
        if rng.random() < rate * 0.2:
            out.append(rng.choice("ACGT"))
    s = "".join(out)
    return s, "5" * len(s)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--platform", choices=["illumina", "ont"], default="illumina")
    ap.add_argument("--reads", type=int, default=200000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="mock1")
    ap.add_argument("--host")
    ap.add_argument("--host_frac", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(a.genomes).read_text())
    G = []
    for g in spec:
        contigs = load(g["fasta"])
        L = sum(map(len, contigs))
        cum, t = [], 0
        for c in contigs:
            t += len(c)
            cum.append(t)
        G.append(dict(g, contigs=contigs, cum=cum, length=L))
    host = None
    if a.host and a.host_frac > 0:
        hc = load(a.host)
        hcum, t = [], 0
        for c in hc:
            t += len(c)
            hcum.append(t)
        host = (hc, hcum)
    w = [g["cells"] * g["length"] for g in G]
    W = sum(w)
    n_host = int(a.reads * a.host_frac) if host else 0
    n_mic = a.reads - n_host
    plan = [(G[i], round(n_mic * w[i] / W)) for i in range(len(G))]
    truth = {"cells": {}, "reads": {}, "host_reads": n_host, "platform": a.platform}
    totc = sum(g["cells"] for g in G)
    for g, n in plan:
        truth["cells"][str(g["taxid"])] = g["cells"] / totc
        truth["reads"][str(g["taxid"])] = n / max(1, n_mic)
    jobs = [(g["contigs"], g["cum"], n) for g, n in plan]
    if host:
        jobs.append((host[0], host[1], n_host))
    recs = []
    for contigs, cum, n in jobs:
        for _ in range(n):
            if a.platform == "illumina":
                ins = max(200, int(rng.gauss(350, 50)))
                frag = pick(contigs, cum, rng, ins)
                if not frag:
                    continue
                if rng.random() < 0.5:
                    frag = rc(frag)
                r1, q1 = err_illumina(frag[:150], rng)
                r2, q2 = err_illumina(rc(frag)[:150], rng)
                recs.append((r1, q1, r2, q2))
            else:
                L = int(min(30000, max(300, rng.lognormvariate(math.log(4000), 0.6))))
                frag = pick(contigs, cum, rng, L)
                if not frag:
                    continue
                if rng.random() < 0.5:
                    frag = rc(frag)
                s, q = err_ont(frag, rng)
                recs.append((s, q))
    rng.shuffle(recs)
    if a.platform == "illumina":
        with gzip.open(out / f"{a.sample}_R1.fastq.gz", "wt", compresslevel=3) as f1, \
                gzip.open(out / f"{a.sample}_R2.fastq.gz", "wt", compresslevel=3) as f2:
            for i, (r1, q1, r2, q2) in enumerate(recs):
                f1.write(f"@{a.sample}.{i}/1\n{r1}\n+\n{q1}\n")
                f2.write(f"@{a.sample}.{i}/2\n{r2}\n+\n{q2}\n")
    else:
        with gzip.open(out / f"{a.sample}.fastq.gz", "wt", compresslevel=3) as f1:
            for i, (s, q) in enumerate(recs):
                f1.write(f"@{a.sample}.{i}\n{s}\n+\n{q}\n")
    (out / f"{a.sample}_truth.json").write_text(json.dumps(truth, indent=1))
    print(f"{a.sample}: {len(recs):,} {'pairs' if a.platform == 'illumina' else 'reads'} "
          f"({n_host:,} host)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
