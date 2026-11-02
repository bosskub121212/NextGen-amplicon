#!/usr/bin/env python3
"""fish_mito_ref.py — a fish mitochondrial BLAST reference from two peer-reviewed sources.

  MIDORI2   Leray et al. 2022, Environmental DNA 4:894-907. GenBank mitochondrial genes,
            quality-controlled, re-released with every GenBank release (newest found on
            Zenodo when this runs).
  MitoFish  Zhu et al. 2023, Molecular Biology and Evolution 40:msad035. Fish mitogenomes
            and partial sequences re-annotated by MitoAnnotator; monthly snapshots. Records
            MitoFish judges heterospecific in the gene region (mislabelled or chimeric) are
            left out.

    python3 fish_mito_ref.py latest
    python3 fish_mito_ref.py build --gene Cytb|12S --midori F.fasta.gz --mitofish mitofishdb.fa.gz
            --tables DIR --out DIR [--makeblastdb PATH] [--midori_release GB273 --mitofish_month 2026/09]

pandas + pyarrow are needed only for the MitoFish tables (env ngamp-amp).
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
MIDORI_FALLBACK = {"release": "GB273", "record": "22685701",
                   "url": "https://zenodo.org/api/records/22685701/files/MIDORI2_GB273_DADA2.zip/content"}
MIDORI_LATEST = "https://zenodo.org/api/records/22685701/versions/latest"
MITOFISH = "https://mitofish.aori.u-tokyo.ac.jp"
GENES = {"Cytb": {"midori": "Cytb", "mitofish": "Cytb"}, "12S": {"midori": "srRNA", "mitofish": "12S_rRNA"}}
UA = {"User-Agent": "NextGen-Amplicon (reference builder; python urllib)"}
_PLACEHOLDER = re.compile(r"^(phylum|class|order|family|genus|species)_", re.I)
_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def _get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return r.read()


def _head_type(url: str) -> str:
    req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-1000"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.headers.get("Content-Type", "")


def latest_midori2(fetch=_get) -> dict:
    """Newest MIDORI2 DADA2 bundle on Zenodo; GB273 when Zenodo cannot be asked."""
    try:
        d = json.loads(fetch(MIDORI_LATEST))
        rel = d["metadata"]["version"]
        key = f"MIDORI2_{rel}_DADA2.zip"
        if not any(f.get("key") == key for f in d.get("files", [])):
            raise ValueError(f"{key} not in the latest record")
        return {"release": rel, "record": str(d["id"]),
                "url": f"https://zenodo.org/api/records/{d['id']}/files/{key}/content"}
    except Exception:
        return dict(MIDORI_FALLBACK)


def latest_mitofish(head=_head_type, today: dt.date | None = None, months: int = 6) -> str | None:
    """Newest monthly MitoFish snapshot ("YYYY/MM") that is published."""
    d = today or dt.date.today()
    y, m = d.year, d.month
    for _ in range(months + 1):
        ym = f"{y}/{m:02d}"
        try:
            if "gzip" in head(f"{MITOFISH}/download/fullseq/{ym}/mitofishdb.fa.gz"):
                return ym
        except Exception:
            pass
        y, m = (y, m - 1) if m > 1 else (y - 1, 12)
    return None


def parse_midori_header(h: str) -> dict:
    """DADA2 'Chordata_7711;...;Cyprinus carpio_7962;' or BLAST 'ACC.1.s.e###root_1;Eukaryota_2759;...;
    Cyprinus_carpio_7962' -> Phylum..Species; placeholder ranks ("order_X") left empty."""
    h = h.lstrip(">").strip().split("###", 1)[-1]
    parts = [p for p in h.split(";") if p != ""]
    out = {}
    for rank, p in zip(RANKS, parts[-6:] if len(parts) >= 6 else parts):
        name = re.sub(r"_\d+$", "", p).strip()
        out[rank] = "" if _PLACEHOLDER.match(name) else name.replace("_", " ")
    for r in RANKS:
        out.setdefault(r, "")
    return out


def _fasta(path: Path):
    if str(path).endswith(".zip"):                     # MIDORI2 BLAST fasta ships zipped
        import io
        import zipfile
        z = zipfile.ZipFile(path)
        fh0 = io.TextIOWrapper(z.open(z.namelist()[0]))
        op = lambda *_: fh0
    else:
        op = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def midori_accession(h: str) -> str:
    """'MH910097.1.16215.17432###...' -> 'MH910097' ('' for DADA2-format headers)."""
    h = h.lstrip(">")
    return h.split("###", 1)[0].split(".")[0] if "###" in h else ""


def midori2_records(path: Path):
    """(accession or MIDORI2_<n>, sequence, lineage); the BLAST-format files carry accessions."""
    for i, (h, s) in enumerate(_fasta(path)):
        yield midori_accession(h) or f"MIDORI2_{i}", s.upper(), parse_midori_header(h)


def _mitofish_tables(tables: Path, gene: str):
    import pandas as pd
    ann = pd.read_parquet(tables / "seq_annotation.parquet")
    ann = ann[ann["gene"] == gene]
    het = pd.read_parquet(tables / "seq_heterospecific.parquet", columns=["gene", "accession", "judge"])
    bad = set(het[(het["gene"] == gene) & (het["judge"] == "Heterospecific")]["accession"])
    st = pd.read_parquet(tables / "seq_taxonid.parquet", columns=["accession", "taxon_id"])
    ts = pd.read_parquet(tables / "taxonid_speciesid.parquet")
    sl = pd.read_parquet(tables / "speciesid_lineageid.parquet")
    tn = pd.read_parquet(tables / "taxonid_name.parquet")
    names = dict(zip(tn["lineage_id"].astype(str), tn["lineage_name"]))
    want = {r.lower(): r for r in RANKS}
    lin = {}
    for sid, rank, lid in zip(sl["species_id"].astype(str), sl["rank_name"], sl["lineage_id"].astype(str)):
        if rank in want:
            lin.setdefault(sid, {})[want[rank]] = names.get(lid, "")
    sp_of = dict(zip(ts["taxon_id"].astype(str), ts["species_id"].astype(str)))
    acc_tax = dict(zip(st["accession"], st["taxon_id"].astype(str)))
    regions = {}
    for a, s, e, strand in zip(ann["accession"], ann["start"], ann["end"], ann["strand"]):
        regions.setdefault(a, (int(s), int(e), int(strand)))
    return regions, bad, acc_tax, sp_of, lin


def mitofish_regions(fasta: Path, tables: Path, gene: str, stats: dict | None = None):
    """(accession, region sequence, lineage) for every MitoFish record with that gene annotated."""
    regions, bad, acc_tax, sp_of, lin = _mitofish_tables(Path(tables), gene)
    n_bad = 0
    for head, seq in _fasta(Path(fasta)):
        for acc in head.split()[0].split(";"):             # identical sequences are collapsed "A;B;C"
            if acc not in regions:
                continue
            if acc in bad:
                n_bad += 1
                break
            s, e, strand = regions[acc]
            sub = seq[s - 1:e].upper()
            if strand < 0:
                sub = sub.translate(_COMP)[::-1]
            lineage = {r: "" for r in RANKS}
            lineage.update(lin.get(sp_of.get(acc_tax.get(acc, ""), ""), {}))
            if sub and lineage.get("Species"):
                yield acc, sub, lineage
            break
    if stats is not None:
        stats["n_heterospecific_removed"] = n_bad


def _join(vals) -> str:
    out = []
    for v in vals:
        for x in v.split(";"):
            if x and x not in out:
                out.append(x)
    return ";".join(out)


def shared_lineage(group: list[dict]) -> dict:
    """Lowest rank every record names alike; higher ranks from the agreed value or, where
    sources file the taxon under different higher names (a renamed order), the first record.
    Records that disagree at Class get the plain top-down walk."""
    for i in range(len(RANKS) - 1, -1, -1):
        vals = {g.get(RANKS[i], "") for g in group}
        if len(vals) == 1 and "" not in vals and len({g.get("Class", "") for g in group} - {""}) <= 1:
            out = {k: "" for k in RANKS}
            for k in RANKS[:i + 1]:
                v = {g.get(k, "") for g in group} - {""}
                out[k] = v.pop() if len(v) == 1 else group[0].get(k, "")
            return out
    out = {k: "" for k in RANKS}
    for k in RANKS:
        v = {g.get(k, "") for g in group} - {""}
        if len(v) > 1:
            break
        if v:
            out[k] = v.pop()
    return out


def merge_shared(recs: list[dict]) -> tuple[list[dict], int]:
    """One record per sequence. A sequence filed under several species (identical
    haplotypes of sister species, or a mislabelled record) becomes one record named to
    their lowest common rank, with every species listed in species_all — so a BLAST hit
    to it can never pick one species by chance, and the hit list cap cannot hide one."""
    by = {}
    for r in recs:
        by.setdefault(r["seq"], []).append(r)
    out, n = [], 0
    for seq, group in by.items():
        if len(group) == 1:
            out.append(group[0])
            continue
        n += 1
        lin = shared_lineage(group)
        sp = sorted({g.get("Species", "") for g in group} - {""})
        out.append({"seq": seq, "source": _join(g["source"] for g in group),
                    "accession": _join(g["accession"] for g in group), **lin,
                    "species_all": ";".join(sp) if len(sp) > 1 else ""})
    return out, n


def build(out_dir: Path, gene: str, midori_fa: Path, mitofish_fa: Path | None, tables: Path | None,
          makeblastdb: str | None, versions: dict | None = None) -> dict:
    """Merged reference fishmito_<gene>.fasta + ids.tsv (+ BLAST db when makeblastdb is given)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    recs, seen = [], {}
    n_mid = n_mf = 0
    stats = {}

    n_amb = [0]

    def add(src, acc, seq, lineage):
        if seq and sum(1 for c in seq if c not in "ACGT") > 0.01 * len(seq):
            n_amb[0] += 1                                # N runs make 60 % hits of a 100 % match
            return
        key = (seq, lineage.get("Species", ""))
        if key in seen:
            r = recs[seen[key]]
            if src not in r["source"].split(";"):
                r["source"] += ";" + src
                r["accession"] += ";" + acc
            return
        seen[key] = len(recs)
        recs.append({"seq": seq, "source": src, "accession": acc, **lineage})

    bad = set()
    if tables and Path(tables, "seq_heterospecific.parquet").exists():
        import pandas as pd
        het = pd.read_parquet(Path(tables) / "seq_heterospecific.parquet", columns=["gene", "accession", "judge"])
        g = GENES[gene]["mitofish"] if gene in GENES else gene
        bad = set(het[(het["gene"] == g) & (het["judge"] == "Heterospecific")]["accession"])
    n_mid_bad = 0
    for acc, seq, lin in midori2_records(Path(midori_fa)):
        if acc in bad:                                 # MitoFish found this GenBank record mislabelled
            n_mid_bad += 1
            continue
        n_mid += 1
        add("MIDORI2", acc, seq, lin)
    if mitofish_fa and tables:
        for acc, seq, lin in mitofish_regions(Path(mitofish_fa), Path(tables), GENES[gene]["mitofish"]
                                              if gene in GENES else gene, stats):
            n_mf += 1
            add("MitoFish", acc, seq, lin)
    recs, n_shared = merge_shared(recs)
    fa = out_dir / f"fishmito_{gene}.fasta"
    with open(fa, "w") as fo, open(out_dir / "ids.tsv", "w") as ft:
        ft.write("\t".join(["id", "source", "accession"] + RANKS + ["species_all"]) + "\n")
        for i, r in enumerate(recs):
            fo.write(f">fm{i}\n{r['seq']}\n")
            ft.write("\t".join([f"fm{i}", r["source"], r["accession"]] + [r.get(k, "") for k in RANKS]
                               + [r.get("species_all", "")]) + "\n")
    if makeblastdb:
        subprocess.run([makeblastdb, "-in", str(fa), "-dbtype", "nucl", "-out", str(out_dir / f"fishmito_{gene}"),
                        "-parse_seqids"], check=True, capture_output=True)
    info = {"gene": gene, "versions": versions or {}, "n_midori2": n_mid, "n_mitofish": n_mf,
            "n_heterospecific_removed": stats.get("n_heterospecific_removed", 0),
            "n_midori2_heterospecific_removed": n_mid_bad, "n_ambiguous_removed": n_amb[0],
            "n_after_dedup": len(recs), "n_shared_seq_merged": n_shared,
            "built": dt.date.today().isoformat(),
            "references": ["Leray M, Knowlton N, Machida RJ (2022) MIDORI2. Environmental DNA 4:894-907",
                           "Zhu T, Sato Y, Sado T, Miya M, Iwasaki W (2023) MitoFish, MitoAnnotator, and MiFish "
                           "Pipeline: updates in 10 years. Mol Biol Evol 40:msad035"]}
    (out_dir / "fishmito.json").write_text(json.dumps(info, indent=1))
    return info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("latest")
    b = sub.add_parser("build")
    b.add_argument("--gene", required=True, choices=list(GENES))
    b.add_argument("--midori", required=True)
    b.add_argument("--mitofish", default="")
    b.add_argument("--tables", default="")
    b.add_argument("--out", required=True)
    b.add_argument("--makeblastdb", default="")
    b.add_argument("--midori_release", default="")
    b.add_argument("--mitofish_month", default="")
    a = ap.parse_args(argv)
    if a.cmd == "latest":
        print(json.dumps({"midori2": latest_midori2(), "mitofish": latest_mitofish()}))
        return 0
    info = build(Path(a.out), a.gene, Path(a.midori), Path(a.mitofish) if a.mitofish else None,
                 Path(a.tables) if a.tables else None, a.makeblastdb or None,
                 {"MIDORI2": a.midori_release, "MitoFish": a.mitofish_month})
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
