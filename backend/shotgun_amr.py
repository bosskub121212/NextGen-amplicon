#!/usr/bin/env python3
"""shotgun_amr.py — resistome from shotgun reads (beta), Illumina or ONT.

Called by shotgun_pipeline.py for every sample, after QC / host removal and
before Kraken 2:

  1. minimap2 aligns the reads to the nucleotide sequences of the NCBI
     AMRFinderPlus reference catalogue (AMR_CDS; the database the WGS pipeline
     installs with AMRFinderPlus — nothing extra to download)
  2. alignments below the identity threshold or too short are dropped; reads
     are pooled per gene family (blaTEM-1, -116 … → blaTEM) and a family is
     called present when its reads cover enough of the gene (breadth)
  3. the read IDs carrying each gene are handed back so that Kraken 2's
     per-read output names the organism the gene sits in — reliable for
     Nanopore reads (kilobases of flanking DNA), tentative for short reads
     and for mobile genes found in many species

Acquired genes only: point-mutation resistance (gyrA, parC, 23S …) cannot be
called from unassembled metagenome reads and is not attempted.
"""
from __future__ import annotations

import csv
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

_PAIR_SUFFIX = re.compile(r"(/[12])+$")


# ── database ─────────────────────────────────────────────────────────────────
def _cds_file(d: Path) -> Path | None:
    for n in ("AMR_CDS", "AMR_CDS.fa", "AMR_CDS.fasta"):
        if (d / n).is_file():
            return d / n
    return None


def find_db(dbp: dict, conda_bases: list[Path], amrfinder_bin: str | None = None) -> Path | None:
    """The AMRFinderPlus data folder holding AMR_CDS and fam.tab."""
    cands = []
    if dbp.get("amrfinder_db"):
        cands.append(Path(dbp["amrfinder_db"]))
    if amrfinder_bin:
        cands.append(Path(amrfinder_bin).resolve().parent.parent / "share" / "amrfinderplus" / "data")
    for b in conda_bases:
        for env in ("ngamp-wgs", "wgs", "amrfinder", "amrfinderplus"):
            cands.append(b / "envs" / env / "share" / "amrfinderplus" / "data")
        envs = b / "envs"
        if envs.is_dir():
            cands += sorted(envs.glob("*/share/amrfinderplus/data"))
    for c in cands:
        for d in (c, c / "latest"):
            if d.is_dir() and _cds_file(d) and (d / "fam.tab").exists():
                return d
        if c.is_dir():   # data/<version>/ folders without a "latest" link
            vs = sorted((x for x in c.iterdir() if x.is_dir() and _cds_file(x)), reverse=True)
            if vs and (vs[0] / "fam.tab").exists():
                return vs[0]
    return None


def db_version(d: Path) -> str:
    v = d / "version.txt"
    return v.read_text().strip() if v.exists() else d.name


def load_catalog(d: Path) -> dict:
    """{sequence id: {allele, family, type, subtype, class, subclass, product, length}}."""
    fam = {}
    with open(d / "fam.tab") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            v = line.rstrip("\n").split("\t")
            if len(v) >= 18:
                fam[v[0]] = {"parent": v[1], "symbol": v[2], "type": v[13], "subtype": v[14], "class": v[15],
                             "subclass": v[16], "family_name": v[17]}
    cat, cur, n = {}, None, 0
    with open(_cds_file(d)) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    cat[cur]["length"] = n
                hdr = line[1:].strip()
                sid = hdr.split()[0]
                f = sid.split("|")
                allele = f[5] if len(f) > 6 else sid
                node = f[6] if len(f) > 6 else allele
                product = (f[7] if len(f) > 7 else "").replace("_", " ")
                # the family's own row; a few alleles point at a node without type —
                # walk up until one has it
                info, k, hops = fam.get(node, {}), node, 0
                while info and not info.get("type") and hops < 5:
                    k = info.get("parent", "")
                    info, hops = fam.get(k, {}), hops + 1
                sym = fam.get(node, {}).get("symbol") or node
                cat[sid] = {"allele": allele, "family": node, "symbol": re.sub(r"_gen$", "", sym),
                            "type": info.get("type", ""),
                            "subtype": info.get("subtype", ""), "class": info.get("class", ""),
                            "subclass": info.get("subclass", ""),
                            "product": product or fam.get(node, {}).get("family_name", "")}
                cur, n = sid, 0
            else:
                n += len(line.strip())
        if cur:
            cat[cur]["length"] = n
    return cat


# ── alignment ────────────────────────────────────────────────────────────────
def align(minimap2, ref: Path, reads: list[Path], platform: str, threads: int, paf: Path,
          logf: Path, run) -> bool:
    preset = "sr" if platform == "illumina" else "map-ont"
    args = ["-c", "-x", preset, "--secondary=no", "-t", threads, ref] + list(reads)
    return run(minimap2, args, logf, stdout_to=paf) == 0 and paf.exists()


def _union(iv: list[tuple[int, int]]) -> int:
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        tot += ce - cs
    return tot


def summarise(paf: Path, cat: dict, platform: str, total_reads: int,
              min_id: float, min_breadth: float) -> tuple[list[dict], dict]:
    """Gene families present in one sample, and {family: set(read ids)}."""
    hits = defaultdict(list)       # family -> [(allele id, tstart, tend, ident, qname)]
    for line in open(paf):
        v = line.rstrip("\n").split("\t")
        if len(v) < 12:
            continue
        qname, qlen, qs, qe = v[0], int(v[1]), int(v[2]), int(v[3])
        tid, tlen, ts, te = v[5], int(v[6]), int(v[7]), int(v[8])
        nm, alen = int(v[9]), int(v[10])
        if tid not in cat or alen == 0:
            continue
        de = next((float(x[5:]) for x in v[12:] if x.startswith("de:f:")), None)
        ident = 1.0 - de if de is not None else nm / alen
        if ident < min_id:
            continue
        span = te - ts
        at_end = ts <= 10 or te >= tlen - 10
        if platform == "illumina":
            # most of the read aligned, unless the read runs off either end of the gene
            if span < 50 or ((qe - qs) < 0.8 * qlen and not at_end):
                continue
        elif span < min(200, 0.5 * tlen):
            continue
        hits[cat[tid]["family"]].append((tid, ts, te, ident, _PAIR_SUFFIX.sub("", qname)))

    rows, reads = [], {}
    for famname, hs in hits.items():
        per_allele = Counter(h[0] for h in hs)
        rep = per_allele.most_common(1)[0][0]
        rlen = cat[rep]["length"] or 1
        iv, bases = [], 0
        for tid, ts, te, ident, q in hs:
            L = cat[tid]["length"] or rlen
            if abs(L - rlen) <= 0.05 * rlen:      # same-length alleles: same coordinates
                iv.append((min(ts, rlen), min(te, rlen)))
            bases += te - ts
        breadth = 100.0 * _union(iv) / rlen
        if breadth < min_breadth:
            continue
        names = {h[4] for h in hs}
        n = len(names)
        c = cat[rep]
        rows.append({"gene": c["symbol"], "family_id": famname,
                     "allele": re.sub(r"_gen$", "", c["allele"]), "type": c["type"],
                     "subtype": c["subtype"], "class": c["class"], "subclass": c["subclass"],
                     "product": c["product"], "length": rlen, "reads": n,
                     "breadth": round(breadth, 1), "depth": round(bases / rlen, 2),
                     "rpkm": round(n / (rlen / 1e3) / max(1e-9, total_reads / 1e6), 2),
                     "identity": round(100 * sum(h[3] for h in hs) / len(hs), 2)})
        reads[famname] = names
    rows.sort(key=lambda r: (r["type"] != "AMR", -r["rpkm"]))
    return rows, reads


# ── Kraken 2 per-read output → organism carrying each gene ───────────────────
def kraken_keep(proc_stdout, wanted: set) -> dict:
    """Stream Kraken 2 --output lines, keep {read id: taxid} for wanted reads."""
    out = {}
    for line in proc_stdout:
        if not wanted:
            continue
        p = line.split("\t", 3)
        if len(p) < 3:
            continue
        rid = _PAIR_SUFFIX.sub("", p[1])
        if rid in wanted:
            try:
                out[rid] = int(p[2])
            except ValueError:
                pass
    return out


def report_tree(report: Path) -> dict:
    """{taxid: (rank code, name, parent taxid)} from a Kraken 2 report."""
    tree, stack = {}, []
    for line in open(report):
        v = line.rstrip("\n").split("\t")
        if len(v) < 6:
            continue
        rank, tid, name = v[3].strip(), int(v[4]), v[5]
        depth = (len(name) - len(name.lstrip(" "))) // 2
        while stack and stack[-1][0] >= depth:
            stack.pop()
        tree[tid] = (rank, name.strip(), stack[-1][1] if stack else 0)
        stack.append((depth, tid))
    return tree


def _up_to(tree: dict, tid: int, want: str) -> str | None:
    for _ in range(60):
        if tid not in tree:
            return None
        rank, name, parent = tree[tid]
        if rank == want:
            return name
        if not rank or rank[0] in "RDKPCOU" or rank == "F":
            return None   # above the wanted rank already
        tid = parent
    return None


def assign_hosts(rows: list[dict], reads: dict, read_tax: dict, tree: dict) -> None:
    """host_taxon / host_share / host_reads on each row (majority of classified reads)."""
    for r in rows:
        ids = reads.get(r.get("family_id", r["gene"]), set())
        gen, spc = Counter(), Counter()
        for rid in ids:
            t = read_tax.get(rid, 0)
            if not t:
                continue
            g = _up_to(tree, t, "G")
            if g:
                gen[g] += 1
                s = _up_to(tree, t, "S")
                if s:
                    spc[s] += 1
        n = sum(gen.values())
        r["host_reads"] = n
        if not n:
            r["host_taxon"], r["host_share"] = "", ""
            continue
        g, gc = gen.most_common(1)[0]
        name = g
        sp = [(s, c) for s, c in spc.items() if s.split()[0] == g]
        if sp:
            s, sc = max(sp, key=lambda x: x[1])
            if sc >= 0.7 * gc:
                name = s
        r["host_taxon"] = name
        r["host_share"] = round(100 * gc / n, 1)


COLS = ["sample", "gene", "family_id", "allele", "type", "subtype", "class", "subclass", "product", "length",
        "reads", "breadth", "depth", "rpkm", "identity", "host_taxon", "host_share", "host_reads"]


def write_tables(out: Path, per_sample: dict[str, list[dict]]) -> int:
    """amr_genes.csv (long) and amr_rpkm.csv (gene × sample). Returns families found."""
    with open(out / "amr_genes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for s, rows in per_sample.items():
            for r in rows:
                w.writerow({"sample": s, **r})
    names = list(per_sample)
    genes = {}
    for s, rows in per_sample.items():
        for r in rows:
            genes.setdefault(r["family_id"], r)
    order = sorted(genes, key=lambda g: (genes[g]["type"] != "AMR", genes[g]["class"], genes[g]["gene"]))
    with open(out / "amr_rpkm.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "family_id", "type", "class"] + names)
        for g in order:
            val = {s: next((r["rpkm"] for r in per_sample[s] if r["family_id"] == g), 0) for s in names}
            w.writerow([genes[g]["gene"], g, genes[g]["type"], genes[g]["class"]] + [val[s] for s in names])
    return len(genes)
