#!/usr/bin/env python3
"""GTDB 16S (SSU) sequences -> DADA2 reference files.

GTDB publishes the 16S genes of its species representatives with the genome
taxonomy in the header:

  >RS_GCF_000005845.2 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...;
   s__Escherichia coli [locus_tag=...] [location=...] [ssu_len=1542] [contig_len=...]

This writes the two files the DADA2 pipeline uses:

  gtdb_r<N>_ssu_toGenus_trainset.fa.gz   >Bacteria;Pseudomonadota;...;Escherichia;
  gtdb_r<N>_ssu_assignSpecies.fa.gz      >RS_GCF_000005845.2 Escherichia coli

Genus names keep GTDB's alphabetic suffixes (Bacillus_A) — they are distinct
genera in GTDB. Placeholder species ("sp030638685") are left out of the species
file, which only names described species. Sequences shorter than 1,200 bp are
dropped (partial genes give weak classifier training).

  python3 build_gtdb_dada2.py --release 232 --out backend/databases/GTDB
  python3 build_gtdb_dada2.py --bac bac120_ssu_reps_r232.fna.gz --ar ar53_ssu_reps_r232.fna.gz \\
                              --release 232 --out .
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import sys
import urllib.request
from pathlib import Path

BASE = "https://data.gtdb.ecogenomic.org/releases/release{r}/{r}.0/genomic_files_reps/"
MD5 = {  # from the release's MD5SUM.txt
    "232": {"bac120_ssu_reps_r232.fna.gz": "3f0a53ff7fd8e3c56fc7246fe3b04342",
            "ar53_ssu_reps_r232.fna.gz": "96dddfbdb5bebafd8e459c02785d5a02"},
}


def fetch(url: str, dest: Path, md5: str | None) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  have {dest.name}")
    else:
        print(f"  downloading {url}")
        tmp = dest.with_suffix(dest.suffix + ".part")
        req = urllib.request.Request(url, headers={"User-Agent": "NextGen-Amplicon"})
        with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as fh:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                fh.write(b)
        tmp.replace(dest)
    if md5:
        h = hashlib.md5(dest.read_bytes()).hexdigest()
        if h != md5:
            dest.unlink()
            raise SystemExit(f"checksum mismatch for {dest.name}: {h} != {md5}")
        print(f"  md5 ok {dest.name}")
    return dest


def records(path: Path):
    name, seq = None, []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = ln[1:].rstrip("\n"), []
            else:
                seq.append(ln.strip())
    if name:
        yield name, "".join(seq)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", default="232")
    ap.add_argument("--bac")
    ap.add_argument("--ar")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_len", type=int, default=1200)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    r = a.release
    srcs = []
    for kind, given in (("bac120", a.bac), ("ar53", a.ar)):
        fn = f"{kind}_ssu_reps_r{r}.fna.gz"
        srcs.append(Path(given) if given else
                    fetch(BASE.format(r=r) + fn, out / fn, MD5.get(r, {}).get(fn)))

    tg = out / f"gtdb_r{r}_ssu_toGenus_trainset.fa.gz"
    sp = out / f"gtdb_r{r}_ssu_assignSpecies.fa.gz"
    n = n_sp = short = 0
    with gzip.open(tg, "wt") as ftg, gzip.open(sp, "wt") as fsp:
        for src in srcs:
            for hdr, seq in records(src):
                if len(seq) < a.min_len:
                    short += 1
                    continue
                acc, rest = hdr.split(" ", 1)
                lineage = rest.split(" [", 1)[0]
                parts = [re.sub(r"^[a-z]__", "", x) for x in lineage.split(";")]
                if len(parts) < 7:
                    continue
                ftg.write(">" + ";".join(parts[:6]) + ";\n" + seq + "\n")
                n += 1
                species = parts[6]
                bits = species.split(" ")
                if len(bits) == 2 and not re.match(r"^sp\d", bits[1]):
                    fsp.write(f">{acc} {species}\n{seq}\n")
                    n_sp += 1
    print(f"wrote {tg.name} ({n} sequences, {short} shorter than {a.min_len} bp dropped)")
    print(f"wrote {sp.name} ({n_sp} sequences with a described species)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
