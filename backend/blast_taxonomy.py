#!/usr/bin/env python3
"""blast_taxonomy.py — names ASVs by BLASTn, the way fish eDNA reports are made.

  1. Local BLASTn against the fish mitochondrial reference (fish_mito_ref.py:
     MIDORI2 + MitoFish). Hits within 0.5 % identity of the best are ties; their lowest
     common ancestor is cut to the deepest rank whose identity threshold is met
     (default species 98.5, genus 95, family 90, order 85, class 80). A species is never
     given below its threshold.
  2. The limit (fish / vertebrates / all) mirrors NCBI's txid7898 / txid7742. An ASV whose
     best overall hit lies outside the limit (human DNA past the blocking primer) keeps
     its real name and is marked target = no — never forced onto a fish.
  3. Optional check against NCBI nt (remote BLASTn, entrez txid limit, as the lab's
     reports do). Agreement per ASV: agree / differ (cut to the common rank, flagged) /
     ncbi_only / local_better / ncbi_no_hit / not_checked. A failed check never fails
     the run.

    python3 blast_taxonomy.py --asv_fasta F --db fishmito_Cytb.fasta --limit fish --ncbi on --out DIR
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
DEFAULT_ID = {"Species": 98.5, "Genus": 95.0, "Family": 90.0, "Order": 85.0, "Class": 80.0}
FISH_CLASSES = {"Actinopteri", "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Dipnomorpha",
                "Hyperoartia", "Myxini", "Actinopterygii", "Dipnoi", "Petromyzonti",
                "Elasmobranchii", "Holocephali"}
NON_VERT_CHORDATES = {"Ascidiacea", "Thaliacea", "Appendicularia", "Leptocardii"}
# NCBI limits: Actinopterygii (ray-finned fishes) + Chondrichthyes (sharks, rays); Vertebrata.
ENTREZ = {"fish": "txid7898[ORGN] OR txid7777[ORGN]", "vertebrates": "txid7742[ORGN]", "all": ""}
ENTREZ_LABEL = {"fish": "ray-finned and cartilaginous fishes (NCBI txid 7898 and txid 7777)",
                "vertebrates": "vertebrates (NCBI txid 7742)", "all": "no taxonomic limit"}
TIE = 0.5
MAX_TARGETS = 5000      # ties must not be cut off: a common CytB haplotype can match hundreds of records
OUT_COLS = ["id", "sequence"] + RANKS + ["best_hit", "pident", "qcov", "source", "accession", "n_ties",
                                         "rank_reached", "target", "ncbi_best", "ncbi_pident",
                                         "ncbi_accession", "agreement", "check"]


def in_limit(lin: dict, limit: str) -> bool:
    if limit == "all":
        return True
    if limit == "fish":
        return lin.get("Class", "") in FISH_CLASSES
    return lin.get("Phylum", "") == "Chordata" and lin.get("Class", "") not in NON_VERT_CHORDATES


def _lca(lins: list[dict]) -> dict:
    """Lowest common ancestor. A rank some records leave empty (NCBI has no order for
    many fish families) is judged on the records that have it; disagreement stops
    the walk. Below the first agreed-but-empty rank the deeper ranks can still agree."""
    out = {r: "" for r in RANKS}
    for r in RANKS:
        vals = {l.get(r, "") for l in lins} - {""}
        if len(vals) > 1:
            break
        if vals:
            out[r] = vals.pop()
    return out


def _lca_ties(lins: list[dict], top: dict) -> dict:
    """LCA of tied hits that survives reclassification: the deepest rank every tie names
    alike is kept, with the ranks above it taken from the agreed value or, where the
    sources file the same taxon under different higher names (an order renamed), from the
    best hit. Falls back to the plain walk when the ties disagree at Class (homonyms)."""
    for i in range(len(RANKS) - 1, -1, -1):
        r = RANKS[i]
        vals = {l.get(r, "") for l in lins}
        if len(vals) == 1 and "" not in vals:
            cls = {l.get("Class", "") for l in lins} - {""}
            if len(cls) > 1:
                break
            out = {k: "" for k in RANKS}
            for k in RANKS[:i + 1]:
                v = {l.get(k, "") for l in lins} - {""}
                out[k] = v.pop() if len(v) == 1 else top.get(k, "")
            return out
    return _lca(lins)


def _species_of(lin: dict) -> list[str]:
    """Species a reference record stands for (a sequence shared by several species lists them all)."""
    sa = [x for x in (lin.get("species_all") or "").split(";") if x]
    return sa or ([lin["Species"]] if lin.get("Species") else [])


def _cut(lin: dict, pident: float, th: dict) -> dict:
    out = dict(lin)
    for r in RANKS[::-1]:
        if r in th and pident < th[r]:
            out[r] = ""
    return out


def _deepest(lin: dict) -> str:
    d = "none"
    for r in RANKS:
        if lin.get(r):
            d = r
    return d


def assign(hits: list[dict], lineages: dict, thresholds: dict, limit: str) -> dict:
    """One ASV's name from its BLAST hits (see module doc)."""
    empty = {**{r: "" for r in RANKS}, "best_hit": "", "pident": "", "qcov": "", "source": "",
             "accession": "", "n_ties": 0, "rank_reached": "none", "target": ""}
    hits = [h for h in hits if h["sseqid"] in lineages]
    if not hits:
        return empty
    key = lambda h: (float(h["bitscore"]), float(h["pident"]))
    best = max(hits, key=key)
    inl = [h for h in hits if in_limit(lineages[h["sseqid"]], limit)]
    best_in = max(inl, key=key) if inl else None
    off = best_in is None or float(best["pident"]) > float(best_in["pident"]) + TIE
    pool = hits if off else inl
    top = max(pool, key=key)
    pid = float(top["pident"])
    ties = [h for h in pool if float(h["pident"]) >= pid - TIE]
    src = lineages[top["sseqid"]]
    lin = _cut(_lca_ties([lineages[h["sseqid"]] for h in ties], src), pid, thresholds)
    tied_sp = sorted({x for h in ties for x in _species_of(lineages[h["sseqid"]])})
    return {**lin, "best_hit": src.get("Species", "") or " | ".join(_species_of(src)), "pident": round(pid, 2),
            "qcov": round(float(top.get("qcovhsp", 0) or 0), 1), "source": src.get("source", ""),
            "accession": src.get("accession", ""), "n_ties": len(tied_sp),
            "rank_reached": _deepest(lin), "target": "no" if off else "yes"}


_NOT_SPECIES = ("uncultured", "unidentified", "environmental", "hybrid", "synthetic", " x ", " sp.", " sp ",
                " cf.", " aff.", " nr.")


def valid_species(name: str) -> bool:
    """A binomial NCBI name that can be reported as a species: not 'sp.', uncultured or a hybrid."""
    n = f" {name or ''} "
    w = (name or "").split()
    return len(w) >= 2 and w[0][:1].isupper() and w[1][:1].islower() and not any(t in n for t in _NOT_SPECIES) \
        and not n.rstrip().endswith(" sp")


def genus_lineages(lineages: dict) -> dict:
    """{genus: {Phylum..Genus}} from the reference, to place an NCBI genus in the tree."""
    out = {}
    for l in lineages.values():
        g = l.get("Genus", "")
        if g and g not in out and l.get("Family"):
            out[g] = {r: l.get(r, "") for r in RANKS[:-1]}
    return out


def combine(local: dict, ncbi: list[dict] | None, thresholds: dict, genus_lin: dict | None = None) -> dict:
    """Merge the NCBI check (best hits of one ASV, sorted best first) into a local call.

    agree         the local species is NCBI's only (named) species at the top
    ncbi_tie      NCBI has two or more named species tied at the top: species cut, flagged
    local_better  NCBI's best is clearly (> 0.5 %) worse than the local hit: local stands
    differ        NCBI's best is another taxon as good or better: cut to the rank both share, flagged
    ncbi_only     no local name; NCBI's name used when its identity passes the thresholds
    outside_limit the local hit is non-target (e.g. human) and NCBI, searched inside the limit, has
                  no clearly closer fish: the local name stands (a fish record that matches human
                  DNA as well is a contaminated record)
    ncbi_no_hit / not_checked
    """
    gl = genus_lin or {}
    c = dict(local)
    c.update({"ncbi_best": "", "ncbi_pident": "", "ncbi_accession": "", "agreement": "not_checked",
              "check": ""})
    if ncbi is None:
        return c
    named = [h for h in ncbi if valid_species(h.get("species", ""))]
    if not ncbi:
        c["agreement"] = "ncbi_no_hit"
        return c
    b0 = ncbi[0]
    c.update({"ncbi_best": b0.get("species", ""), "ncbi_pident": round(float(b0["pident"]), 2),
              "ncbi_accession": b0.get("accession", "")})
    if local.get("rank_reached") == "none":
        c["agreement"] = "ncbi_only"
        if not named:
            return c
        b = named[0]
        bp = float(b["pident"])
        tops = [h for h in named if float(h["pident"]) >= bp - TIE]
        sp = {h["species"] for h in tops}
        gen = {h["species"].split()[0] for h in tops}
        g = b["species"].split()[0]
        up = gl.get(g, {})
        if bp >= thresholds.get("Species", 98.5) and len(sp) == 1:
            # NCBI was searched inside the limit (fish / vertebrates), so a name from it is a target
            c.update({**up, "Genus": g, "Species": b["species"], "rank_reached": "Species", "check": "yes",
                      "target": "yes",
                      "best_hit": b["species"], "pident": round(bp, 2), "source": "NCBI nt",
                      "accession": b.get("accession", "")})
        elif bp >= thresholds.get("Genus", 95.0) and len(gen) == 1:
            c.update({**up, "Genus": g, "rank_reached": "Genus", "check": "yes", "target": "yes",
                      "best_hit": b["species"],
                      "pident": round(bp, 2), "source": "NCBI nt", "accession": b.get("accession", "")})
        return c
    if not named:
        c["agreement"] = "ncbi_no_hit"
        return c
    b = named[0]
    bp = float(b["pident"])
    top_sp = {h["species"] for h in named if float(h["pident"]) >= bp - TIE}
    lsp = local.get("Species", "")
    lp = float(local.get("pident") or 0)
    lsps = {lsp} if lsp else set(x.strip() for x in (local.get("best_hit") or "").split("|") if x.strip())
    if lsps & top_sp:
        if len(top_sp) == 1 or not lsp:
            c["agreement"] = "agree"
            return c
        c["agreement"], c["check"] = "ncbi_tie", "yes"
        gens = {x.split()[0] for x in top_sp}
        c["Species"] = ""
        if len(gens) > 1:
            c["Genus"] = ""
        c["rank_reached"] = _deepest(c)
        return c
    if bp < lp - TIE:
        c["agreement"] = "local_better"
        return c
    # NCBI was searched INSIDE the limit (fish), so it cannot see the non-target organism the local
    # reference found (human DNA past the blocking primer). A fish record matching a human sequence
    # just as well is a contaminated GenBank record, not a fish: the local name stands. Only a fish
    # clearly closer than the local non-target hit overrules it.
    if str(local.get("target") or "") == "no" and bp <= lp + TIE:
        c["agreement"] = "outside_limit"
        return c
    c["agreement"], c["check"] = "differ", "yes"
    g = b["species"].split()[0]
    other = gl.get(g) or ({r: local.get(r, "") for r in RANKS[:-1]} if g == local.get("Genus") else None)
    keep = {r: "" for r in RANKS}
    if other:
        for r in RANKS[:-1]:
            if local.get(r) and local.get(r) == other.get(r):
                keep[r] = local[r]
            elif local.get(r) or other.get(r):
                break
    else:                                            # NCBI's genus is not in the reference: keep down to order only
        for r in RANKS[:3]:
            keep[r] = local.get(r, "")
    for r in RANKS:
        c[r] = keep[r]
    c["rank_reached"] = _deepest(c)
    return c


def load_ids(db_fasta: Path) -> dict:
    lin = {}
    with open(Path(db_fasta).parent / "ids.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin[r["id"]] = r
    return lin


def find_blastn(explicit: str = "") -> str | None:
    if explicit:
        return explicit
    bases = []
    for b in (os.environ.get("CONDA_PREFIX", ""), str(Path.home() / "miniconda3"), str(Path.home() / "miniforge3"),
              str(Path.home() / "mambaforge"), str(Path.home() / "anaconda3")):
        if b:
            p = Path(b)
            bases.append(p.parent.parent if p.parent.name == "envs" else p)
    for b in bases:
        for env in ("ngamp-amp", "ngamp-phylo", "ngamp-blast"):
            x = b / "envs" / env / "bin" / "blastn"
            if x.exists():
                return str(x)
    return shutil.which("blastn")


FMT = "6 qseqid sseqid pident length qlen bitscore qcovhsp evalue"


def _parse(txt: str, cols: list[str]) -> dict:
    out = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < len(cols):
            continue
        d = dict(zip(cols, f))
        if "sseqid" in d and "|" in d["sseqid"]:          # "lcl|fm12" with -parse_seqids databases
            d["sseqid"] = d["sseqid"].rstrip("|").split("|")[-1]
        out.setdefault(d["qseqid"], []).append(d)
    return out


def run_local(fasta: Path, db_fasta: Path, blastn: str, threads: int = 4) -> dict:
    db = str(Path(db_fasta).with_suffix(""))
    p = subprocess.run([blastn, "-task", "blastn", "-query", str(fasta), "-db", db, "-outfmt", FMT,
                        "-max_target_seqs", str(MAX_TARGETS), "-evalue", "1e-10", "-qcov_hsp_perc", "80",
                        "-num_threads", str(threads)], capture_output=True, text=True, check=True)
    return _parse(p.stdout, FMT.split()[1:])


NCBI_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _http(url: str, data: bytes | None = None) -> bytes:
    req = urllib.request.Request(url, data=data, headers={"User-Agent": "NextGen-Amplicon blast_taxonomy"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def run_ncbi(fasta: Path, entrez: str, timeout: int = 3600, http=_http, sleep=time.sleep,
             batch: int = 100, status: dict | None = None) -> dict | None:
    """NCBI nt check through the BLAST URL API (HTTPS, works behind web proxies; NCBI's
    rules: one status poll per minute per search). {qseqid: [hits best first]}; a query
    of a batch that failed is left out (= not checked). None when no batch finished.
    NGAMP_NCBI_EMAIL, when set, is sent as EMAIL as NCBI asks of scripted use."""
    import re
    import urllib.parse
    url = os.environ.get("NGAMP_NCBI_URL", NCBI_URL)
    email = os.environ.get("NGAMP_NCBI_EMAIL", "")
    st = status if status is not None else {}
    st.update({"batches": 0, "failed_batches": 0, "errors": []})
    recs = _fasta(Path(fasta))
    out = {}
    t_end = time.time() + timeout
    nb = (len(recs) + batch - 1) // batch
    for bi, i in enumerate(range(0, len(recs), batch), 1):
        st["batches"] += 1
        try:
            if time.time() > t_end:
                raise TimeoutError("NCBI time limit reached")
            q = "".join(f">{n}\n{s}\n" for n, s in recs[i:i + batch])
            form = {"CMD": "Put", "PROGRAM": "blastn", "MEGABLAST": "on", "DATABASE": "nt", "QUERY": q,
                    "HITLIST_SIZE": "20", "EXPECT": "1e-10", "TOOL": "NextGen-Amplicon"}
            if email:
                form["EMAIL"] = email
            if entrez:
                form["ENTREZ_QUERY"] = entrez
            r = http(url, urllib.parse.urlencode(form).encode()).decode(errors="replace")
            rid = re.search(r"RID = (\S+)", r).group(1)
            rtoe = int((re.search(r"RTOE = (\d+)", r) or [0, "30"])[1])
            sleep(max(10, min(rtoe, 60)))
            while True:
                stt = http(f"{url}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}").decode(errors="replace")
                m = re.search(r"Status=(\w+)", stt)
                if m and m.group(1) == "READY":
                    break
                if m and m.group(1) in ("FAILED", "UNKNOWN"):
                    raise RuntimeError(f"NCBI search {m.group(1)}")
                if time.time() > t_end:
                    raise TimeoutError("NCBI time limit reached")
                sleep(60)
            js = json.loads(http(f"{url}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}"))
            got = {}
            for rep in js.get("BlastOutput2", []):
                s = rep["report"]["results"]["search"]
                title = (s.get("query_title") or "").split()
                if not title:
                    continue
                qn, ql = title[0], float(s.get("query_len") or 1)
                rows = []
                for h in s.get("hits", []):
                    d, hs = h["description"][0], h["hsps"][0]
                    cov = 100 * (hs["query_to"] - hs["query_from"] + 1) / ql
                    if cov < 80:
                        continue
                    name = d.get("sciname", "")
                    rows.append({"species": name, "genus": name.split(" ")[0] if name else "",
                                 "pident": round(100 * hs["identity"] / hs["align_len"], 2),
                                 "accession": d.get("accession", ""), "bits": hs.get("bit_score", 0)})
                rows.sort(key=lambda r: (-r["bits"], -r["pident"]))
                got[qn] = rows
            out.update(got)
            print(f"[blast] NCBI batch {bi}/{nb} done ({len(got)} ASVs)", file=sys.stderr, flush=True)
        except Exception as ex:
            st["failed_batches"] += 1
            st["errors"].append(f"batch {bi}: {ex}")
            print(f"[blast] NCBI batch {bi}/{nb} failed: {ex}", file=sys.stderr, flush=True)
    return out if st["failed_batches"] < st["batches"] or not recs else None


def _ncbi_hits(ncbi: dict | None, q: str) -> list | None:
    """Hits of one ASV; None (= not checked) when NCBI failed or its batch did."""
    if ncbi is None or q not in ncbi:
        return None
    return ncbi[q]


def _fasta(p: Path) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    for line in Path(p).read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        out.append((name, "".join(seq)))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asv_fasta", required=True)
    ap.add_argument("--db", required=True, help="fishmito_<gene>.fasta (ids.tsv beside it)")
    ap.add_argument("--limit", default="fish", choices=["fish", "vertebrates", "all"])
    ap.add_argument("--ncbi", default="on", choices=["on", "off"])
    ap.add_argument("--id", default="98.5,95,90,85,80", help="species,genus,family,order,class %% identity")
    ap.add_argument("--blastn", default="")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    vals = [float(x) for x in a.id.split(",")]
    th = dict(zip(["Species", "Genus", "Family", "Order", "Class"], vals))
    info = {"limit": a.limit, "thresholds": th, "tie_window": TIE, "db": a.db, "ncbi": "off"}
    try:
        info["reference"] = json.loads((Path(a.db).parent / "fishmito.json").read_text())
    except Exception:
        info["reference"] = {}
    asvs = _fasta(Path(a.asv_fasta))
    blastn = find_blastn(a.blastn)
    rows = []
    if not blastn:
        info["error"] = "blastn not found (env ngamp-amp: download_databases.sh --fishmito)"
        (out / "blast_info.json").write_text(json.dumps(info, indent=1))
        print(f"[blast] {info['error']}", file=sys.stderr)
        return 3
    try:
        v = subprocess.run([blastn, "-version"], capture_output=True, text=True).stdout.splitlines()
        info["blastn"] = v[0] if v else ""
    except Exception:
        pass
    lin = load_ids(Path(a.db))
    t0 = time.time()
    local = run_local(Path(a.asv_fasta), Path(a.db), blastn, a.threads)
    info["local_seconds"] = round(time.time() - t0, 1)
    with open(out / "blast_local.tsv", "w") as fh:
        for q, hs in local.items():
            for h in hs:
                fh.write("\t".join(str(h.get(k, "")) for k in ["qseqid", "sseqid", "pident", "length", "qlen",
                                                               "bitscore", "qcovhsp", "evalue"]) + "\n")
    ncbi = None
    if a.ncbi == "on":
        t0 = time.time()
        nst = {}
        ncbi = run_ncbi(Path(a.asv_fasta), ENTREZ[a.limit], status=nst)
        info["ncbi"] = "failed" if ncbi is None else ("partial" if nst.get("failed_batches") else "done")
        info["ncbi_batches"] = {k: nst.get(k) for k in ("batches", "failed_batches", "errors")}
        info["ncbi_entrez"] = ENTREZ[a.limit]
        info["ncbi_limit_label"] = ENTREZ_LABEL[a.limit]
        info["ncbi_seconds"] = round(time.time() - t0, 1)
        if ncbi is not None:
            with open(out / "blast_ncbi.tsv", "w") as fh:
                for q, hs in ncbi.items():
                    for h in hs:
                        fh.write(f"{q}\t{h['accession']}\t{h['species']}\t{h['pident']}\n")
    gl = genus_lineages(lin)
    for q, seq in asvs:
        r = assign(local.get(q, []), lin, th, a.limit)
        r = combine(r, _ncbi_hits(ncbi, q), th, gl)
        rows.append({"id": q, "sequence": seq, **r})
    with open(out / "taxonomy_blast.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    agg = {}
    for r in rows:
        agg[r["agreement"]] = agg.get(r["agreement"], 0) + 1
    info.update({"n_asv": len(rows), "n_species": sum(1 for r in rows if r["Species"]),
                 "n_unassigned": sum(1 for r in rows if r["rank_reached"] == "none"),
                 "n_non_target": sum(1 for r in rows if r["target"] == "no"), "agreement": agg,
                 "n_check": sum(1 for r in rows if r["check"] == "yes")})
    (out / "blast_info.json").write_text(json.dumps(info, indent=1))
    print(f"[blast] {len(rows)} ASVs: {info['n_species']} to species, {info['n_unassigned']} unassigned, "
          f"{info['n_non_target']} non-target; NCBI check {info['ncbi']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
