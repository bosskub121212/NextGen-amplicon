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
