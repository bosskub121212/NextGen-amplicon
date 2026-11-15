#!/usr/bin/env python3
"""NCBI nt check of the OTUs the reference left unassigned (VSEARCH OTU runs, e.g. fungal ITS
against UNITE). Primers made for one group also amplify other DNA where they happen to bind;
UNITE (fungi only) has no name for it, so it shows up as 'Unassigned'. Each such OTU is
BLASTed (megablast) against NCBI nt through the BLAST URL API, its closest record and the
record's lineage (E-utilities taxonomy) are kept, and it is named only as deep as the hit
supports. Results describe what is there, not how much: the reads of off-target DNA depend
on primer binding, not on abundance.

Writes, in the run folder:
  ncbi_unassigned.tsv   one row per checked OTU
  ncbi_unassigned.json  when / how many / failed batches / errors
On a finished run:  python3 ncbi_unassigned.py --run_dir results/<job>
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import blast_taxonomy as BT                                     # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# deepest rank kept for an identity over >= MIN_COV % of the OTU (below: the group only)
RANK_FLOORS = (("genus", 95.0), ("family", 90.0), ("phylum", 80.0))
MIN_COV = 80.0
TIE = 0.995                     # hits within 0.5 % of the best bit score are ties
COLS = ["otu", "reads", "length", "group", "named_as", "rank", "ncbi_hit", "pident", "qcov",
        "accession", "taxid", "phylum", "family", "genus", "ties", "title"]
GROUP_ORDER = ["Bacteria", "Archaea", "Fungi", "Plants", "Animals", "Other eukaryotes", "Viruses",
               "Unclassified record", "No NCBI match"]


def _id(x: str) -> str:
    return str(x).split(";")[0].split()[0]


def unassigned_otus(run_dir: Path) -> list[tuple[str, str, int]]:
    """(otu id, sequence, reads over all samples) of every OTU with no phylum, most reads first."""
    run_dir = Path(run_dir)
    tax = {}
    with open(run_dir / "taxonomy.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            tax[_id(r.get("ASV", ""))] = (r.get("Phylum") or "").strip()
    reads = {}
    with open(run_dir / "asv_table.csv", newline="") as fh:
        rd = csv.reader(fh)
        next(rd, None)
        for r in rd:
            if r:
                reads[_id(r[0])] = sum(int(round(float(v or 0))) for v in r[1:])
    seqs = {_id(n): s for n, s in BT._fasta(run_dir / "otus.fasta")}
    out = [(k, seqs.get(k, ""), reads.get(k, 0)) for k, ph in tax.items() if not ph]
    return sorted([x for x in out if x[1] and x[2] > 0], key=lambda x: -x[2])


def lineages(taxids, http=BT._http, sleep=time.sleep) -> dict:
    """{taxid: {"superkingdom", "kingdom", "phylum", "family", "genus", "name"}} from NCBI taxonomy."""
    import xml.etree.ElementTree as ET            # NCBI's own XML, fetched over HTTPS
    ids = sorted({str(t) for t in taxids if str(t).strip() and str(t) != "0"})
    out = {}
    for i in range(0, len(ids), 200):
        q = urllib.parse.urlencode({"db": "taxonomy", "id": ",".join(ids[i:i + 200]), "retmode": "xml",
                                    "tool": "NextGen-Amplicon"})
        root = ET.fromstring(http(f"{EUTILS}?{q}"))
        for tx in root.findall("Taxon"):
            tid = tx.findtext("TaxId", "")
            lin = {"name": tx.findtext("ScientificName", ""), "rank_self": tx.findtext("Rank", "")}
            for a in tx.findall("LineageEx/Taxon"):
                rk = a.findtext("Rank", "")
                if rk == "domain":
                    rk = "superkingdom"
                if rk in ("superkingdom", "kingdom", "phylum", "family", "genus"):
                    lin[rk] = a.findtext("ScientificName", "")
            if lin["rank_self"] in ("genus", "family", "phylum"):
                lin[lin["rank_self"]] = lin["name"]
            elif lin["rank_self"] in ("superkingdom", "domain"):
                lin["superkingdom"] = lin["name"]
            out[tid] = lin
        if i + 200 < len(ids):
            sleep(0.5)
    return out


def group_of(lin: dict | None) -> str:
    if not lin:
        return "Unclassified record"
    sk, k = lin.get("superkingdom", ""), lin.get("kingdom", "")
    if sk in ("Bacteria", "Archaea"):
        return sk
    if sk == "Viruses" or "virus" in (lin.get("name") or "").lower():
        return "Viruses"
    if sk == "Eukaryota":
        return {"Fungi": "Fungi", "Viridiplantae": "Plants", "Metazoa": "Animals"}.get(k, "Other eukaryotes")
    return "Unclassified record"


def name_hit(hits: list[dict], lins: dict) -> dict:
    """The OTU's row: closest record, group, and a name only as deep as identity, coverage and
    agreement of the tied hits allow."""
    if not hits:
        return {"group": "No NCBI match", "named_as": "No NCBI match", "rank": ""}
    best = hits[0]
    ties = [h for h in hits if h.get("bits", 0) >= TIE * best.get("bits", 0)]
    ls = [lins.get(str(h.get("taxid", ""))) or {} for h in ties]
    groups = {group_of(x) for x in ls}
    grp = groups.pop() if len(groups) == 1 else "Unclassified record"
    row = {"ncbi_hit": best.get("species", ""), "pident": best.get("pident", ""), "qcov": best.get("qcov", ""),
           "accession": best.get("accession", ""), "taxid": best.get("taxid", ""), "title": best.get("title", ""),
           "ties": len(ties), "group": grp, "named_as": grp, "rank": "group"}
    l0 = ls[0] if ls else {}
    for rk in ("phylum", "family", "genus"):
        row[rk] = l0.get(rk, "")
    pid, cov = float(best.get("pident") or 0), float(best.get("qcov") or 0)
    if cov >= MIN_COV:
        for rk, floor in RANK_FLOORS:
            names = {x.get(rk, "") for x in ls}
            if pid >= floor and len(names) == 1 and "" not in names:
                row["named_as"], row["rank"] = names.pop(), rk
                break
    return row


def run(run_dir, max_otus: int = 200, http=BT._http, sleep=time.sleep, timeout: int = 3600, log=print) -> dict:
    run_dir = Path(run_dir)
    otus = unassigned_otus(run_dir)
    info = {"date": date.today().isoformat(), "database": "NCBI nt (BLASTn megablast, URL API)",
            "unassigned_otus": len(otus), "checked": 0, "rank_floors": dict(RANK_FLOORS), "min_cov": MIN_COV}
    if len(otus) > max_otus:
        info["left_out"] = len(otus) - max_otus
        otus = otus[:max_otus]
    rows = []
    if otus:
        q = run_dir / "ncbi_unassigned_query.fasta"
        q.write_text("".join(f">{k}\n{s}\n" for k, s, _ in otus))
        st = {}
        log(f"  NCBI check of {len(otus)} unassigned OTUs (BLASTn megablast vs nt; needs internet)")
        hits = BT.run_ncbi(q, "", timeout=timeout, http=http, sleep=sleep, batch=50, status=st, min_cov=0)
        info.update({"batches": st.get("batches", 0), "failed_batches": st.get("failed_batches", 0),
                     "errors": st.get("errors", [])})
        if hits is None:
            info["error"] = "; ".join(st.get("errors", [])[:3]) or "NCBI did not answer"
            (run_dir / "ncbi_unassigned.json").write_text(json.dumps(info, indent=1))
            log(f"  [WARN] NCBI check failed: {info['error']}")
            return info
        try:
            lins = lineages({h.get("taxid") for v in hits.values() for h in v[:20]}, http=http, sleep=sleep)
        except Exception as ex:                    # names still shown, groups unknown
            lins, info["lineage_error"] = {}, str(ex)
        for k, s, n in otus:
            if k not in hits:
                continue                           # its batch failed: not checked
            r = name_hit(hits[k], lins)
            r.update({"otu": k, "reads": n, "length": len(s)})
            rows.append(r)
        q.unlink(missing_ok=True)
    info["checked"] = len(rows)
    with open(run_dir / "ncbi_unassigned.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (run_dir / "ncbi_unassigned.json").write_text(json.dumps(info, indent=1))
    log(f"  NCBI check: {len(rows)} OTUs named" + (f", {info['failed_batches']} batch(es) failed"
                                                 if info.get("failed_batches") else ""))
    return info


def load(run_dir) -> tuple[dict, dict | None]:
    """({otu id: row}, info) of a run; ({}, None) when the check was not run."""
    run_dir = Path(run_dir)
    jp = run_dir / "ncbi_unassigned.json"
    if not jp.exists():
        return {}, None
    info = json.loads(jp.read_text())
    rows = {}
    tp = run_dir / "ncbi_unassigned.tsv"
    if tp.exists():
        with open(tp, newline="") as fh:
            rows = {_id(r["otu"]): r for r in csv.DictReader(fh, delimiter="\t")}
    return rows, info


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--run_dir", required=True, help="finished VSEARCH OTU run (taxonomy.csv, asv_table.csv, otus.fasta)")
    p.add_argument("--max_otus", type=int, default=200, help="most-read unassigned OTUs checked")
    a = p.parse_args(argv)
    info = run(a.run_dir, a.max_otus)
    return 1 if info.get("error") else 0


if __name__ == "__main__":
    sys.exit(main())
