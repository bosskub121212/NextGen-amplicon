#!/usr/bin/env python3
"""Run report (Report PDF / Word) of a 12S / CytB run named by BLAST: the BLAST names and
method are described, fish vs non-target reads are reported, CytB is called CytB, no
bacterial wording, and the denoise cap is explained instead of blamed on read quality."""
from __future__ import annotations
import sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend")); sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)

def main():
    print("Fish (12S / CytB) run report with BLAST names")
    try:
        import report_builder as rb
        from test_amp_fast import make_amp_run
        with tempfile.TemporaryDirectory() as td:
            d = make_amp_run(Path(td) / "run")
            (d / "pipeline.log").write_text("  [denoise cap] S2-12S: 179261 -> 100000 reads (capped to bound denoising runtime)\n")
            run = rb.Run(d, set())
            html = rb.build_profile_html(run)
            check("BLAST rows loaded by sequence", len(run.blast) == 4, len(run.blast))
            check("CytB report title", "CytB Amplicon Sequencing Report" in html)
            check("BLAST method described", "nucleotide BLAST" in html and "98.5%" in html
                  and "taxonomy_table_dada2.csv" in html)
            check("no naive-Bayes paragraph as the naming method", "classified with the DADA2 implementation" not in html)
            check("fish / non-target section", "Fish and non-target reads (BLAST)" in html and "Fish %" in html)
            check("fish taxa table names Cyprinus carpio", "Cyprinus carpio" in html)
            check("human DNA reported as non-target", "Human DNA" in html)
            check("no Bacteria / Archaea sentence", "classified as Bacteria" not in html)
            check("no V3-V4 / 16S species note", "V3&ndash;V4 amplicon" not in html and "full-length 16S is" not in html)
            check("no rRNA operon note", "rRNA operons" not in html and "environmental DNA reflects" in html)
            check("denoise cap explained", "randomly subsampled" in html)
            check("BLAST+ and MitoFish referenced, not the naive Bayesian classifier",
                  "BLAST+: architecture" in html and "MitoFish, MitoAnnotator" in html
                  and "Naive Bayesian classifier for rapid" not in html)
            check("config rows say BLAST", "Names by" in html and "NCBI nt check" in html)
            check("chimera loss named, not blamed on read quality",
                  "lost most at chimera removal" in html and "comes from PCR" in html)
            check("composition explains phylum-level BLAST names", "places it only approximately" in html)
            import csv as _csv
            bp = d / "blast" / "taxonomy_blast.csv"
            rows = list(_csv.DictReader(open(bp)))
            for r in rows:
                if r["id"] == "q2":      # a fish cut to phylum by conflicting hits, 99% identity
                    r.update(Class="", Order="", Family="", Genus="", Species="", pident="99.0",
                             rank_reached="Phylum", check="yes")
                if r["id"] == "q4":      # a fish far from every reference
                    r.update(Phylum="Chordata", pident="70.0", rank_reached="Phylum", target="yes", check="yes")
            with open(bp, "w", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
            h3 = rb.build_profile_html(rb.Run(d, set()))
            check("phylum-only fish: distant vs conflicting hits told apart",
                  "Closest to a fish, below 80" in h3 and "name cut to phylum (conflicting" in h3)
            check("flagged list without repeats", "flagged for checking (" in h3 and "Fish-like" in h3)
            import blast_taxonomy as bt
            th = {"Species": 98.5, "Genus": 95, "Family": 90, "Order": 85, "Class": 80}
            c = bt.combine({"rank_reached": "none", "target": ""},
                           [{"species": "Cyprinus carpio", "pident": "99.5", "accession": "X"}], th, {})
            check("NCBI-only name is a target (searched inside the fish limit)",
                  c.get("Species") == "Cyprinus carpio" and c.get("target") == "yes", str(c))
            hum = {"Phylum": "Chordata", "Class": "Mammalia", "Genus": "Homo", "Species": "Homo sapiens",
                   "best_hit": "Homo sapiens", "pident": "98.8", "rank_reached": "Species", "target": "no"}
            c = bt.combine(hum, [{"species": "Parablennius sanguinolentus", "pident": "98.8", "accession": "A"}], th, {})
            check("human ASV vs a fish-limited NCBI record of the same identity: human kept, not flagged",
                  c["agreement"] == "outside_limit" and c["Species"] == "Homo sapiens" and not c["check"], str(c))
            c = bt.combine(dict(hum, pident="85.0"),
                           [{"species": "Parablennius sanguinolentus", "pident": "99.0", "accession": "A"}], th, {})
            check("a fish clearly closer than the non-target hit still overrules (flagged)",
                  c["agreement"] == "differ" and c["check"] == "yes", str(c))
            rows = list(_csv.DictReader(open(bp)))
            for r in rows:                # an older run: NCBI-only name left with an empty target
                if r["id"] == "q4":
                    r.update(Phylum="Chordata", Class="Actinopteri", Genus="Cyprinus", Species="Cyprinus carpio",
                             target="", agreement="ncbi_only", pident="99.5", rank_reached="Species")
            with open(bp, "w", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
            run4 = rb.Run(d, set())
            sec = rb._blast_section(run4, run4.samples, lambda h, r, c: str(r))
            check("older run: NCBI-only fish counted as fish", "No match" not in "".join(sec)[:0] and
                  any("Cyprinus carpio" in x for x in sec), "".join(sec)[:300])
            # same run without BLAST: DADA2 wording stays, no BLAST section
            (d / "blast" / "BLAST_SKIPPED.txt").write_text("blastn not found")
            run2 = rb.Run(d, set())
            h2 = rb.build_profile_html(run2)
            check("skipped BLAST: no BLAST section, DADA2 method", not run2.blast and "non-target reads (BLAST)" not in h2
                  and "classified with the DADA2 implementation" in h2)
            check("skipped BLAST: mitochondrial species note", "mitochondrial marker" in h2)
    except Exception:
        check("fish report builds", False, traceback.format_exc(limit=4))
    return 1 if FAILS else 0

if __name__ == "__main__":
    sys.exit(main())
