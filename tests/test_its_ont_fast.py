#!/usr/bin/env python3
"""Nanopore ITS: routed to the VSEARCH OTU pipeline with UNITE (its_pipeline.R is DADA2
paired-end only), MinKNOW/Guppy chunk files (fastq_runid_<id>_<n>) merged into one
sample, UNITE headers parsed, fungal identity floors."""
from __future__ import annotations
import gzip, json, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)

def main():
    print("Nanopore ITS (VSEARCH + UNITE)")
    try:
        import qiime2_vsearch_pipeline as qv
        h = ("Abrothallus_subhalei|MT153946|SH1227328.10FU|reps|k__Fungi;p__Ascomycota;c__Dothideomycetes;"
             "o__Abrothallales;f__Abrothallaceae;g__Abrothallus;s__Abrothallus_subhalei")
        t = qv.parse_lineage(h)
        check("UNITE header parsed", t.get("superkingdom") == "Fungi" and t.get("genus") == "Abrothallus"
              and t.get("species") == "Abrothallus subhalei", str(t))
        t2 = qv.parse_lineage("Fungi_sp|X1|SH9|reps|k__Fungi;p__Ascomycota;c__unidentified;o__unidentified;"
                              "f__unidentified;g__unidentified;s__Fungi_sp")
        check("UNITE 'unidentified' / '_sp' left blank", t2 == {"superkingdom": "Fungi", "phylum": "Ascomycota"}, str(t2))
        check("SILVA trainset still parsed", qv.parse_lineage("Bacteria;Firmicutes;Bacilli;Lactobacillales;"
              "Lactobacillaceae;Lactobacillus;").get("genus") == "Lactobacillus")
        check("ITS identity floors (Vu et al. 2019)", qv._ITS_IDENTITY_FLOOR["species"] == 98.41
              and qv._ITS_IDENTITY_FLOOR["phylum"] == 75.0)
        with tempfile.TemporaryDirectory() as td:
            d = Path(td) / "in"; d.mkdir()
            names = [f"fastq_runid_cae8cbea-7b23_{i}.fastq.gz" for i in range(3)]
            for i, n in enumerate(names):
                with gzip.open(d / n, "wt") as fh:
                    fh.write(f"@r{i}\nACGT\n+\nIIII\n")
            (d / "sample_manifest.json").write_text(json.dumps(
                [{"sample": "PGU-ITS", "file1": n, "file2": ""} for n in names]))
            s = qv.find_fastq(d, Path(td) / "out" / "merged_reads")
            with gzip.open(s["PGU-ITS"], "rt") as fh:
                n_reads = sum(1 for l in fh) // 4
            check("chunk files of one sample merged", list(s) == ["PGU-ITS"] and n_reads == 3, f"{s} {n_reads}")
        import ont_wgs_pipeline as w
        fs = [f"fastq_runid_cae8cbea-7b23-4538-a575-f68dfe0df9cc_{i}.fastq.gz" for i in range(4)] \
            + ["fastq_runid_cae8cbea-7b23-4538-a575-f68dfe0df9cc_4_0.fastq.gz"]
        g = w.group_fastqs(fs, {f: "PGU-ITS" for f in fs})
        check("fastq_runid chunks grouped into the folder's sample", [(x["sample"], len(x["files"])) for x in g]
              == [("PGU-ITS", 5)], str([(x["sample"], len(x["files"])) for x in g]))
        g2 = w.group_fastqs(["L21_1.fastq.gz", "L21_2.fastq.gz"], {})
        check("replicates _1/_2 are not merged", len(g2) == 2)
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        i = m.find('if marker in ("ITS1", "ITS2", "ITS") and params.sequencerType in ("ont", "qiime2_vsearch")')
        j = m.find('elif marker in ("ITS1", "ITS2", "ITS"):')
        check("main: Nanopore ITS -> VSEARCH pipeline with UNITE before the DADA2 ITS branch",
              0 < i < j and "qiime2_vsearch_pipeline.py" in m[i:j] and "UNITE_ITS1" in m[i:j]
              and '"--marker",         region' in m[i:j])
        h_ = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
        check("form: ITS offers MinKNOW folder merging and the read-length window",
              'it.key === "ITS1" || it.key === "ITS2"' in h_ and "Minimum read length (Nanopore / VSEARCH)" in h_
              and "fastq_runid_" in h_)
    except Exception:
        check("its ont checks", False, traceback.format_exc(limit=4))
    return 1 if FAILS else 0

if __name__ == "__main__":
    sys.exit(main())
