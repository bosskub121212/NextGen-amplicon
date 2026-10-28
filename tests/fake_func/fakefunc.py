#!/usr/bin/env python3
"""Stand-ins for metaphlan, humann, humann_* utilities and seqkit (tests only).

Each wrapper calls:  python3 fakefunc.py <tool> <args...>
They write files in the real tools' formats, small and fixed, so the pipeline,
tables, panels, stats and report can be tested without the 50 GB databases.
FAKE_HUMANN_FAIL=s2,s3 makes humann fail for those sample basenames.
FAKE_SEQKIT_LOG=path records seqkit argument lists (one line per call).
"""
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path


def opt(args, *names, default=None):
    for n in names:
        if n in args:
            return args[args.index(n) + 1]
    return default


def opener(p):
    if str(p).endswith(".bz2"):
        import bz2
        return bz2.open(p, "rt")
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read(p):
    with open(p) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    return lines[0], lines[1:]


def write(p, head, rows):
    with open(p, "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def metaphlan(a):
    out = opt(a, "-o", "--output_file")
    rows = [("k__Bacteria", "2", "100.0"),
            ("k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia", "561", "60.0"),
            ("k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli", "562", "60.0"),
            ("k__Bacteria|p__Bacillota|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus", "1386", "40.0"),
            ("k__Bacteria|p__Bacillota|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_subtilis", "1423", "40.0")]
    with open(out, "w") as fh:
        fh.write("#mpa_vJun23_CHOCOPhlAnSGB_202307\n#/fake metaphlan\n#1000 reads processed\n")
        fh.write("#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n")
        for r in rows:
            fh.write("\t".join(r) + "\t\n")


def humann(a):
    inp, out, base = opt(a, "--input", "-i"), opt(a, "--output", "-o"), opt(a, "--output-basename")
    base = base or Path(inp).name.split(".")[0]
    if base in (os.environ.get("FAKE_HUMANN_FAIL") or "").split(","):
        print("CRITICAL ERROR: fake failure", file=sys.stderr)
        sys.exit(1)
    o = Path(out)
    o.mkdir(parents=True, exist_ok=True)
    ec, bs = "g__Escherichia.s__Escherichia_coli", "g__Bacillus.s__Bacillus_subtilis"
    write(o / f"{base}_genefamilies.tsv", ["# Gene Family", f"{base}_Abundance-RPKs"], [
        ("UNMAPPED", 400), ("UniRef90_K1", 100), (f"UniRef90_K1|{ec}", 70), (f"UniRef90_K1|{bs}", 30),
        ("UniRef90_B1", 50), (f"UniRef90_B1|{bs}", 50), ("UniRef90_N1", 450), ("UniRef90_N1|unclassified", 450)])
    write(o / f"{base}_pathabundance.tsv", ["# Pathway", f"{base}_Abundance"], [
        ("UNMAPPED", 300), ("UNINTEGRATED", 500), ("PWY-1: fake pathway one", 150),
        (f"PWY-1: fake pathway one|{ec}", 150), ("PWY-2: fake pathway two", 50)])
    write(o / f"{base}_pathcoverage.tsv", ["# Pathway", f"{base}_Coverage"], [
        ("PWY-1: fake pathway one", 1.0), (f"PWY-1: fake pathway one|{ec}", 1.0),
        ("PWY-2: fake pathway two", 0.3)])
    log = opt(a, "--o-log")
    if log:
        Path(log).write_text("01/01/2026 INFO: Unaligned reads after nucleotide alignment: 61.2000000000 %\n"
                             "01/01/2026 INFO: Unaligned reads after translated alignment: 40.1000000000 %\n")


def join(a):
    d, out, name = Path(opt(a, "-i", "--input")), opt(a, "-o", "--output"), opt(a, "--file_name")
    files = sorted(p for p in (d.rglob("*") if "--search-subdirectories" in a else d.iterdir())
                   if p.is_file() and (not name or name in p.name))
    cols, data, first = [], defaultdict(dict), None
    for f in files:
        h, rows = read(f)
        first = first or h[0]
        cols.append(h[1])
        for r in rows:
            data[r[0]][h[1]] = r[1]
    write(out, [first] + cols, [[k] + [data[k].get(c, 0) for c in cols] for k in data])


def renorm(a):
    h, rows = read(opt(a, "-i", "--input"))
    tot = [sum(float(r[i + 1]) for r in rows if "|" not in r[0]) or 1 for i in range(len(h) - 1)]
    if "--update-snames" in a:
        h = [h[0]] + [c.replace("-RPKs", "-CPM") for c in h[1:]]
    write(opt(a, "-o", "--output"), h,
          [[r[0]] + [round(float(r[i + 1]) / tot[i] * 1e6, 4) for i in range(len(h) - 1)] for r in rows])


def regroup(a):
    m = {}
    with opener(opt(a, "-c", "--custom")) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            for u in f[1:]:
                m[u] = f[0]
    h, rows = read(opt(a, "-i", "--input"))
    acc = defaultdict(lambda: [0.0] * (len(h) - 1))
    for r in rows:
        feat, _, strat = r[0].partition("|")
        g = feat if feat == "UNMAPPED" else m.get(feat, "UNGROUPED")
        key = g + ("|" + strat if strat else "")
        for i, v in enumerate(r[1:]):
            acc[key][i] += float(v)
    write(opt(a, "-o", "--output"), h, [[k] + v for k, v in acc.items()])


def rename(a):
    names = {}
    with opener(opt(a, "-c", "--custom")) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            names[f[0]] = f[1]
    h, rows = read(opt(a, "-i", "--input"))
    out = []
    for r in rows:
        feat, sep, strat = r[0].partition("|")
        if feat in names:
            feat = f"{feat}: {names[feat]}"
        out.append([feat + sep + strat] + r[1:])
    write(opt(a, "-o", "--output"), h, out)


def split(a):
    inp, od = Path(opt(a, "-i", "--input")), Path(opt(a, "-o", "--output"))
    h, rows = read(inp)
    stem = inp.name.rsplit(".", 1)[0]
    write(od / f"{stem}_stratified.tsv", h, [r for r in rows if "|" in r[0]])
    write(od / f"{stem}_unstratified.tsv", h, [r for r in rows if "|" not in r[0]])


def seqkit(a):
    lg = os.environ.get("FAKE_SEQKIT_LOG")
    if lg:
        with open(lg, "a") as fh:
            fh.write(" ".join(a) + "\n")
    out = opt(a, "-o", "--out-file")
    src = [x for x in a[1:] if Path(x).exists() and x != out]
    if out and src:
        with opener(src[0]) as fi, (gzip.open(out, "wt") if out.endswith(".gz") else open(out, "w")) as fo:
            fo.write(fi.read())


def kraken2(a):
    rep = opt(a, "--report")
    Path(rep).write_text(
        "100.00\t10\t0\tR\t1\troot\n100.00\t10\t0\tD\t2\t  Bacteria\n"
        "60.00\t6\t0\tG\t561\t    Escherichia\n60.00\t6\t6\tS\t562\t      Escherichia coli\n"
        "40.00\t4\t0\tG\t1386\t    Bacillus\n40.00\t4\t4\tS\t1423\t      Bacillus subtilis\n")


TOOLS = {"kraken2": kraken2, "metaphlan": metaphlan, "humann": humann, "humann_join_tables": join,
         "humann_renorm_table": renorm, "humann_regroup_table": regroup,
         "humann_rename_table": rename, "humann_split_stratified_table": split, "seqkit": seqkit}

if __name__ == "__main__":
    tool, args = sys.argv[1], sys.argv[2:]
    if "--version" in args:
        print({"humann": "humann v3.9", "metaphlan": "MetaPhlAn version 4.1.1 (11 Mar 2024)"}.get(tool, "1.0"))
        sys.exit(0)
    TOOLS[tool](args)
