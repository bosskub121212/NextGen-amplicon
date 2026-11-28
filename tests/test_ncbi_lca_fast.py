#!/usr/bin/env python3
"""ITS named by BLASTn against NCBI nt + LCA (ncbi_lca.py; coral ITS2 eDNA, DiBattista et al. 2025):
the LCA rules, the script with NCBI replaced by canned answers, the pipeline / backend / form wiring,
and the report and Download results of such a run."""
from __future__ import annotations
import csv, json, re, sys, tempfile, traceback, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(ROOT / "tests"))
FAILS = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


CORAL = {"superkingdom": "Eukaryota", "kingdom": "Metazoa", "phylum": "Cnidaria", "class": "Anthozoa",
         "order": "Scleractinia", "family": "Poritidae", "genus": "Porites"}


def hit(sp, pid, cov, taxid):
    return {"species": sp, "pident": pid, "qcov": cov, "taxid": taxid, "accession": "X" + str(taxid), "bits": pid}


def _canned(L):
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        fa = t / "asv.fasta"
        fa.write_text(">q1\nACGT\n>q2\nGGGG\n>q3\nTTTT\n")
        L.BT.run_ncbi = lambda *a_, **k: {"q1": [hit("Porites lobata", 100, 100, 1)],
                                          "q2": [hit("Cylista ornata", 88, 100, 1)], "q3": []}
        L.NU.lineages = lambda ids: {"1": dict(CORAL)}
        rc = L.main(["--asv_fasta", str(fa), "--out", str(t / "o")])
        rows = list(csv.DictReader(open(t / "o" / "taxonomy_ncbi_lca.csv")))
        info = json.loads((t / "o" / "ncbi_lca.json").read_text())
        check("ncbi_lca.py: a row per ASV, in FASTA order, statuses, info", rc == 0
              and [r["status"] for r in rows] == ["named", "below_threshold", "no_hit"]
              and rows[0]["Species"] == "Porites lobata" and info["identity"] == 90 and info["hits"] == 10
              and info["program"] == "blastn", (rc, rows, info))
        L.BT.run_ncbi = lambda *a_, **k: None
        check("ncbi_lca.py: NCBI unreachable -> non-zero exit, said in ncbi_lca.json",
              L.main(["--asv_fasta", str(fa), "--out", str(t / "o2")]) != 0
              and "error" in json.loads((t / "o2" / "ncbi_lca.json").read_text()))



def _parallel_check():
    """run_ncbi with several searches at once: every batch done, every query answered."""
    import blast_taxonomy as BT, threading, urllib.parse
    rid_q, lock = {}, threading.Lock()

    def http(url, data=None):
        if data:
            f = urllib.parse.parse_qs(data.decode())
            qs = [l[1:] for l in f["QUERY"][0].splitlines() if l.startswith(">")]
            with lock:
                rid = f"R{len(rid_q) + 1}"
                rid_q[rid] = qs
            return f"RID = {rid}\nRTOE = 1\n".encode()
        rid = url.split("RID=")[-1]
        if "SearchInfo" in url:
            return b"Status=READY"
        rep = [{"report": {"results": {"search": {"query_title": q, "query_len": 4, "hits": [
            {"description": [{"sciname": "Porites lobata", "taxid": 1, "accession": "A"}],
             "hsps": [{"query_from": 1, "query_to": 4, "identity": 4, "align_len": 4, "bit_score": 8}]}]}}}}
            for q in rid_q[rid]]
        return json.dumps({"BlastOutput2": rep}).encode()
    with tempfile.TemporaryDirectory() as t:
        fa = Path(t) / "a.fasta"
        fa.write_text("".join(f">q{i}\nACGT\n" for i in range(1, 12)))
        st = {}
        got = BT.run_ncbi(fa, "", http=http, sleep=lambda x: None, batch=3, status=st, min_cov=0,
                          hitlist=10, megablast=False, parallel=3)
        rid_q.clear()
        _http0 = http

        def http_stuck(url, data=None):                # the first search never finishes
            if not data and "SearchInfo" in url and url.endswith("RID=R1"):
                return b"Status=WAITING"
            return _http0(url, data)
        st2 = {}
        got2 = BT.run_ncbi(fa, "", http=http_stuck, sleep=lambda x: None, batch=3, status=st2, min_cov=0,
                           hitlist=10, megablast=False, parallel=2, search_timeout=1e-6, retries=1)
        check("a search NCBI leaves waiting is given up and sent again", got2 is not None and len(got2) == 11
              and st2["failed_batches"] == 0, (st2, sorted(got2 or {})))
        check("NCBI searches 3 at a time: 4 batches, all 11 ASVs answered, no failure",
              got is not None and len(got) == 11 and st["batches"] == 4 and st["failed_batches"] == 0
              and got["q7"][0]["qcov"] == 100.0, (st, sorted(got or {})))


def main():
    print("ITS named at NCBI nt (BLASTn + LCA)")
    try:
        import ncbi_lca as L
        lins = {"1": dict(CORAL), "2": dict(CORAL), "3": {"superkingdom": "Eukaryota"},
                "4": dict(CORAL, family="Agariciidae", genus="Pavona")}
        a = L.assign([hit("Porites lobata", 100, 100, 1), hit("Porites lutea", 99.5, 100, 2)], lins, 90, 100, 10)
        check("LCA: two Porites species -> genus Porites", a["Genus"] == "Porites" and a["Species"] == ""
              and a["rank_reached"] == "Genus" and a["Order"] == "Scleractinia", a)
        a = L.assign([hit("Porites lobata", 100, 100, 1)] * 3, lins, 90, 100, 10)
        check("LCA: every hit the same species -> species", a["Species"] == "Porites lobata", a)
        a = L.assign([hit("Porites lobata", 100, 100, 1), hit("uncultured eukaryote", 100, 100, 3)], lins, 90, 100, 10)
        check("LCA: an uncultured eukaryote hit as good -> Eukaryota only", a["Kingdom"] == "Eukaryota"
              and not a["Phylum"], a)
        a = L.assign([hit("Porites lobata", 100, 100, 1), hit("Pavona decussata", 99, 100, 4)], lins, 90, 100, 10)
        check("LCA: two families -> order Scleractinia", a["Order"] == "Scleractinia" and not a["Family"], a)
        a = L.assign([hit("Cylista ornata", 88.3, 99, 1)], lins, 90, 100, 10)
        check("below 90 % identity -> Unassigned, best hit kept", a["status"] == "below_threshold"
              and a["best_hit"] == "Cylista ornata" and not a["Kingdom"], a)
        a = L.assign([hit("Porites lobata", 100, 98, 1)], lins, 90, 100, 10)
        b = L.assign([hit("Porites lobata", 100, 99.6, 1)], lins, 90, 100, 10)
        check("coverage under 100 % -> not used; 99.6 % is 100 % as BLAST reports it (qcovs)",
              a["status"] == "below_threshold" and b["status"] == "named", (a, b))
        a = L.assign([hit("Porites lobata", 100, 100, 1)] + [hit("Pavona decussata", 100, 100, 4)] * 2, lins, 90, 100, 1)
        check("only the top N hits are looked at", a["Species"] == "Porites lobata" and a["n_hits"] == 1, a)
        check("no hit / NCBI failed", L.assign([], lins, 90, 100, 10)["status"] == "no_hit"
              and L.assign(None, lins, 90, 100, 10)["status"] == "not_checked")
        uh = [hit("uncultured eukaryote", 100, 100, 3)] * 4 + [hit("Porites lobata", 98, 100, 1),
                                                                hit("Porites lutea", 97, 100, 2)]
        a = L.assign(uh, lins, 90, 100, 10, skip_uninformative=True)
        b = L.assign(uh, lins, 90, 100, 10)
        c = L.assign(uh[:4], lins, 90, 100, 10, skip_uninformative=True)
        check("uncultured-eukaryote records passed over when asked: named from the records that place it",
              a["Genus"] == "Porites" and a["n_uninformative_skipped"] == 4 and b["Kingdom"] == "Eukaryota"
              and not b["Phylum"] and c["Kingdom"] == "Eukaryota", (a, b, c))
        a = L.assign([hit("Porites sp.", 100, 100, 1), hit("Porites sp.", 100, 100, 2)], lins, 90, 100, 10)
        check("'sp.' records name the genus, never a species", a["Genus"] == "Porites" and not a["Species"], a)

        _orig = (L.BT.run_ncbi, L.NU.lineages)          # put back after: other tests use the real ones
        try:
            _canned(L)
        finally:
            L.BT.run_ncbi, L.NU.lineages = _orig

        _parallel_check()
        r = (ROOT / "backend" / "r_scripts" / "its_pipeline.R").read_text()
        check("its_pipeline.R: --tax_source ncbi_lca calls ncbi_lca.py; UNITE / FUNGuild skipped; failure stops",
              '"--tax_source"' in r and "ncbi_lca.py" in r and "if (use_ncbi) {" in r
              and "NCBI BLAST naming failed" in r and 'tax_method <- "ncbi_lca"' in r)
        check("its_pipeline.R: minLen, pool, chimera method and merge overlap from the form",
              "minOverlap=opt$merge_overlap, maxMismatch=0" in r and "method=chim_method" in r
              and "merge_overlap = opt$merge_overlap" in r)
        m = (ROOT / "backend" / "main.py").read_text()
        check("main.py: ITS Illumina passes minLen / pool / chimera / overlap and the NCBI LCA settings",
              '"--tax_source", "ncbi_lca"' in m and '"--merge_overlap", str(int(params.merge_overlap))' in m
              and '"--pool", params.pool or "FALSE"' in m and 'its_tax:            str   = "unite"' in m
              and '"its_tax must be unite or ncbi_lca"' in m)
        f = (ROOT / "frontend-beta" / "index.html").read_text()
        check("form: coral ITS2 preset (SCLER primers, NCBI + LCA, minLen 100, overlap 20, pool TRUE), "
              "no database needed for NCBI names",
              'f:"GARTCTTTGAACGCAAATGGC"' in f and 'its_tax:"ncbi_lca", its_min_len:100' in f
              and '["its_tax", "Names from (Illumina)"' in f and "var dbOk  = ncbiNames ? true" in f)

        # report + Download results of an ITS run named at NCBI
        import report_builder as RB, packages as PK
        from test_its_pacbio_pkg_fast import make_its, SEQ
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            d = make_its(t / "its")
            sm = json.loads((d / "summary.json").read_text())
            sm.update({"tax_method": "ncbi_lca", "lca": {"identity": 90, "coverage": 100, "hits": 10},
                       "minLen": 100, "pool": "TRUE", "merge_overlap": 20, "chimera_method": "consensus",
                       "primer_f": "GARTCTTTGAACGCAAATGGC"})
            (d / "summary.json").write_text(json.dumps(sm))
            (d / "ncbi_lca").mkdir()
            with open(d / "ncbi_lca" / "taxonomy_ncbi_lca.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, L.COLS)
                w.writeheader()
                for i, q in enumerate(SEQ):
                    w.writerow({"id": f"q{i + 1}", "sequence": q, "Genus": "Porites" if i == 0 else "",
                                "rank_reached": "Genus" if i == 0 else "none",
                                "status": "named" if i == 0 else "below_threshold", "n_hits": 10,
                                "n_passing": 10 if i == 0 else 0, "best_hit": "Porites lobata",
                                "best_pident": 100 if i == 0 else 85, "best_qcov": 100})
            (d / "ncbi_lca" / "ncbi_lca.json").write_text(json.dumps({"date": "2026-11-27"}))
            html = RB.build_profile_html(RB.Run(d, RB.KITOME, "x"))
            body = html.split("References")[0]
            refs = html.split("References")[-1]
            check("report: NCBI LCA method (90 % / 100 % / 10 hits, DiBattista 2025), the names section, "
                  "minLen 100, overlap 20, pool TRUE",
                  "lowest common ancestor" in body and "Names from NCBI GenBank (BLASTn + LCA)" in body
                  and "minimum length 100 bp" in body and "minimum overlap 20 bp" in body
                  and "pool = TRUE" in body and "SCLER5.8SForw" in body, re.findall(r"minimum [^,.]+", body)[:4])
            check("report: no UNITE / fungal wording, references BLAST and DiBattista, not UNITE",
                  "UNITE species hypothesis" not in body and "between fungi" not in body
                  and "Altschul" in refs and "DiBattista" in refs and "Abarenkov" not in refs
                  and "Wang, Q." not in refs)
            m = PK.build(d, None, "amplicon", None, False, t / "pk", job_label="U2BI-188", formats=("html",),
                         tools={}, combined=True)
            z = zipfile.ZipFile(m["files"][0]["path"])
            n = z.namelist()
            rep = z.read(next(x for x in n if x.endswith("_report.html"))).decode()
            check("Download results: 5.Taxonomy_Identification/BLAST/taxonomy_ncbi_lca.csv with the count "
                  "table's IDs, and the report names it",
                  any(x.endswith("5.Taxonomy_Identification/BLAST/taxonomy_ncbi_lca.csv") for x in n)
                  and "5.Taxonomy_Identification/BLAST/taxonomy_ncbi_lca.csv" in rep
                  and "Report_U2BI-188_Illumina_report.html" in " ".join(n), [x for x in n if "BLAST" in x])
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
