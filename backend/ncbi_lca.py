#!/usr/bin/env python3
"""ncbi_lca.py — names ASVs by BLASTn against NCBI GenBank nt and the lowest common ancestor (LCA) of
the best hits, for markers no curated reference covers (coral ITS2 eDNA, other eukaryotic barcodes).
Follows DiBattista et al. 2025 (Coral Reefs 44:163-178): BLASTn vs nt, the top hits (default 10), of
these only hits with identity >= 90 % over 100 % of the ASV; the ASV is named down to the deepest rank
on which every such hit agrees. No hit passing -> Unassigned (the best hit is still listed).

  1. BLASTn (task blastn, not megablast: divergent ITS2 at 90-95 % identity is found) through the NCBI
     BLAST URL API, in batches (HTTPS; needs internet). A failed batch leaves its ASVs 'not_checked'.
  2. Lineage of every hit's taxid from NCBI taxonomy (E-utilities): superkingdom ... genus; species is
     the record's name when it is a proper binomial (not 'sp.', 'cf.', uncultured, environmental).
  3. LCA: a rank a record skips (NCBI has no order for many families) is judged on the records that
     have it; a record named no deeper (uncultured eukaryote) or a disagreement stops the walk.
     Species only when every passing hit names the same species.

    python3 ncbi_lca.py --asv_fasta asv.fasta --out DIR [--identity 90 --coverage 100 --hits 10]

Writes DIR/taxonomy_ncbi_lca.csv (one row per ASV, in the order of the FASTA) and DIR/ncbi_lca.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import blast_taxonomy as BT                                     # noqa: E402
import ncbi_unassigned as NU                                    # noqa: E402

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
COLS = ["id", "sequence"] + RANKS + ["rank_reached", "status", "n_hits", "n_passing", "n_uninformative_skipped",
                                     "best_hit", "best_pident", "best_qcov", "best_accession", "hits"]


def informative(hit: dict, lins: dict) -> bool:
    """A record that places the sequence below the domain: 'uncultured eukaryote' and other environmental
    records filed under Eukaryota alone say nothing about what the ASV is."""
    l = lineage_of(hit, lins)
    return any(l.get(r) for r in RANKS)


def lineage_of(hit: dict, lins: dict) -> dict:
    """Kingdom..Species of one NCBI record."""
    l = lins.get(str(hit.get("taxid", ""))) or {}
    sp = hit.get("species", "") if BT.valid_species(hit.get("species", "")) else ""
    out = {"Domain": l.get("superkingdom", ""), "Kingdom": l.get("kingdom", ""), "Phylum": l.get("phylum", ""),
           "Class": l.get("class", ""), "Order": l.get("order", ""), "Family": l.get("family", ""),
           "Genus": l.get("genus", "") or (sp.split()[0] if sp else ""), "Species": sp}
    if sp and out["Genus"] and sp.split()[0] != out["Genus"]:
        out["Species"] = ""                     # a record filed under another genus than its name says
    return out


def lca(lins: list[dict]) -> dict:
    """Deepest ranks every record agrees on (empty values judged on the records that have them);
    species only when every record names the same one. Kingdom is NCBI's kingdom (Metazoa, Fungi ...)
    or, for records filed under none (uncultured eukaryote), the domain (Eukaryota)."""
    out = {r: "" for r in ["Domain"] + RANKS}
    if not lins:
        return {r: "" for r in RANKS}
    rk = ["Domain"] + RANKS
    for i, r in enumerate(rk):
        vals = [l.get(r, "") for l in lins]
        # a record named no deeper than this (uncultured eukaryote, a genus-only record at Species)
        # is unresolved here: it stops the walk. One that skips this rank only (NCBI has no order
        # for many families) is judged on the ranks it has.
        if any(not any(l.get(x, "") for x in rk[i:]) for l in lins):
            break
        if r == "Species":
            if len(set(vals)) == 1 and vals[0]:
                out[r] = vals[0]
            break
        have = set(vals) - {""}
        if len(have) > 1:
            break
        if have:
            out[r] = have.pop()
    out["Kingdom"] = out["Kingdom"] or out["Domain"]
    out.pop("Domain")
    return out


def assign(hits: list[dict] | None, lins: dict, identity: float, coverage: float, top: int,
           skip_uninformative: bool = False) -> dict:
    """One ASV's row from its NCBI hits (best first): the top `top` hits that reach the identity over the
    coverage, their LCA. skip_uninformative: among the passing hits, records that name nothing below the
    domain (uncultured / environmental eukaryotes) are passed over when records that do name something
    pass too; an ASV whose passing hits are all such records is named Eukaryota, not left Unassigned."""
    row = {r: "" for r in RANKS}
    row.update({"rank_reached": "none", "n_hits": 0, "n_passing": 0, "n_uninformative_skipped": 0,
                "best_hit": "", "best_pident": "", "best_qcov": "", "best_accession": "", "hits": ""})
    if hits is None:
        row["status"] = "not_checked"
        return row
    if not hits:
        row["n_hits"] = 0
        row["status"] = "no_hit"
        return row
    # coverage as BLAST reports it (qcovs, a whole percent): 99.6 % of a 280 bp ASV is 100 %, as NCBI's own
    # 'Query cover' and BLAST+ -outfmt qcovs show it — the method's 100 % means a full-length hit
    ok_ = lambda h: (float(h.get("pident") or 0) >= identity
                     and int(float(h.get("qcov") or 0) + 0.5) >= coverage)
    if skip_uninformative:
        passing = [h for h in hits if ok_(h)]
        inf = [h for h in passing if informative(h, lins)]
        if inf:
            row["n_uninformative_skipped"] = sum(1 for h in passing[:passing.index(inf[min(len(inf), top) - 1]) + 1]
                                                if not informative(h, lins))
            use = inf[:top]
        else:
            use = passing[:top]
        shown = (use or hits)[:top]
        row["n_hits"] = len(hits)
    else:
        hits = hits[:top]
        use = [h for h in hits if ok_(h)]
        shown = hits
        row["n_hits"] = len(hits)
    b = (use or hits)[0]
    row.update({"best_hit": b.get("species", ""), "best_pident": b.get("pident", ""), "best_qcov": b.get("qcov", ""),
                "best_accession": b.get("accession", "")})
    row["n_passing"] = len(use)
    row["hits"] = "; ".join(f"{h.get('species', '')} {float(h.get('pident') or 0):.1f}%/{float(h.get('qcov') or 0):.0f}%"
                            for h in shown)
    if not use:
        row["status"] = "below_threshold"
        return row
    row.update(lca([lineage_of(h, lins) for h in use]))
    deep = [r for r in RANKS if row[r]]
    row["rank_reached"] = deep[-1] if deep else "none"
    row["status"] = "named" if deep else "no_agreement"
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--asv_fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--identity", type=float, default=90.0)
    ap.add_argument("--coverage", type=float, default=100.0)
    ap.add_argument("--hits", type=int, default=10)
    ap.add_argument("--task", default="blastn", choices=["blastn", "megablast"])
    ap.add_argument("--batch", type=int, default=25)
    ap.add_argument("--timeout", type=int, default=4 * 3600)
    ap.add_argument("--parallel", type=int, default=4, help="NCBI searches at once (started 10 s apart)")
    ap.add_argument("--from_hits", default="", help="ncbi_hits.json.gz of an earlier search: LCA again, no BLAST")
    ap.add_argument("--skip_uninformative", default="on", choices=["on", "off"],
                    help="pass over records named only 'Eukaryota' (uncultured eukaryote ...): up to 5x the hits are fetched")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    asvs = BT._fasta(Path(a.asv_fasta))
    info = {"date": date.today().isoformat(), "database": "NCBI GenBank nt (BLAST URL API)",
            "program": a.task, "identity": a.identity, "coverage": a.coverage, "hits": a.hits,
            "asvs": len(asvs), "method": "LCA of the top hits passing identity and coverage",
            "skip_uninformative": a.skip_uninformative == "on"}
    fetch = min(max(a.hits * 5, a.hits), 50) if a.skip_uninformative == "on" else a.hits
    info["fetched_per_asv"] = fetch
    print(f"[ncbi_lca] {len(asvs)} ASVs: {a.task} vs NCBI nt, top {a.hits} hits, >= {a.identity}% identity over "
          f">= {a.coverage}% of the ASV, LCA", flush=True)
    st = {}
    import gzip
    if a.from_hits:                                # the LCA again on hits kept from an earlier search
        with gzip.open(a.from_hits, "rt") as fh:
            hits = json.load(fh)
        info["from_hits"] = Path(a.from_hits).name
    else:
        hits = BT.run_ncbi(Path(a.asv_fasta), "", timeout=a.timeout, batch=a.batch, status=st, min_cov=0,
                           hitlist=fetch, megablast=(a.task == "megablast"), parallel=max(1, min(a.parallel, 5)),
                           search_timeout=20 * 60, retries=2)
        if hits is not None:                       # every hit kept: the names can be checked or redone offline
            with gzip.open(out / "ncbi_hits.json.gz", "wt") as fh:
                json.dump(hits, fh)
    info.update({"batches": st.get("batches", 0), "failed_batches": st.get("failed_batches", 0),
                 "errors": st.get("errors", [])[:10]})
    if hits is None:
        info["error"] = "NCBI could not be reached or every batch failed"
        (out / "ncbi_lca.json").write_text(json.dumps(info, indent=1))
        print("[ncbi_lca] NCBI BLAST failed: " + "; ".join(st.get("errors", [])[:3]), file=sys.stderr, flush=True)
        return 2
    taxids = {str(h.get("taxid", "")) for hs in hits.values() for h in hs[:fetch]}
    try:
        lins = NU.lineages(taxids)
    except Exception as ex:                       # names need lineages: no lineage, no LCA
        info["error"] = f"NCBI taxonomy (E-utilities) failed: {ex}"
        (out / "ncbi_lca.json").write_text(json.dumps(info, indent=1))
        print(f"[ncbi_lca] {info['error']}", file=sys.stderr, flush=True)
        return 3
    rows = []
    for name, seq in asvs:
        r = assign(BT._ncbi_hits(hits, name), lins, a.identity, a.coverage, a.hits,
                   skip_uninformative=(a.skip_uninformative == "on"))
        r.update({"id": name, "sequence": seq})
        rows.append(r)
    with open(out / "taxonomy_ncbi_lca.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS)
        w.writeheader()
        w.writerows(rows)
    cnt = {}
    for r in rows:
        cnt[r["status"]] = cnt.get(r["status"], 0) + 1
    info["status_counts"] = cnt
    info["rank_counts"] = {rk: sum(1 for r in rows if r["rank_reached"] == rk) for rk in RANKS}
    (out / "ncbi_lca.json").write_text(json.dumps(info, indent=1))
    print(f"[ncbi_lca] done: {cnt}; to species {info['rank_counts']['Species']}, genus "
          f"{info['rank_counts']['Genus']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
