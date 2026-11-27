#!/usr/bin/env python3
"""Download results for ITS Illumina (DADA2 + UNITE) and PacBio full-length 16S: the 8-folder result
reads their tables (old ITS runs with raw UNITE names, old PacBio runs with one row per sample)."""
from __future__ import annotations
import csv, json, sys, tempfile, traceback, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend")); sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
SEQ = ["ACGT" * 40, "TTGG" * 40, "CCAA" * 40, "GATC" * 40]
S = ["F1", "F2", "F3"]
COUNTS = [[500, 30, 200], [300, 120, 0], [0, 80, 60], [20, 5, 40]]


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def make_its(d: Path) -> Path:
    """An ITS2 Illumina run made before the pipeline cleaned UNITE names (prefixes, 'unidentified')."""
    d.mkdir(parents=True)
    with open(d / "asv_table.csv", "w") as fh:
        fh.write('"sequence",' + ",".join(f'"{s}"' for s in S) + "\n")
        for q, c in zip(SEQ, COUNTS):
            fh.write(f'"{q}",' + ",".join(map(str, c)) + "\n")
    tax = [["k__Fungi", "p__Ascomycota", "c__Eurotiomycetes", "o__Eurotiales", "f__Aspergillaceae", "g__Aspergillus", "s__Aspergillus_niger"],
           ["k__Fungi", "p__Basidiomycota", "c__Agaricomycetes", "o__Polyporales", "f__Polyporaceae", "g__unidentified", "s__unidentified"],
           ["k__Fungi", "p__Ascomycota", "c__Sordariomycetes", "o__Hypocreales", "f__Nectriaceae", "g__Fusarium", "s__Fusarium_sp"],
           ["k__Fungi", "NA", "NA", "NA", "NA", "NA", "NA"]]
    with open(d / "taxonomy_table.csv", "w") as fh:
        fh.write('"","Kingdom","Phylum","Class","Order","Family","Genus","Species","sequence"\n')
        for q, t in zip(SEQ, tax):
            fh.write(f'"{q}",' + ",".join(f'"{x}"' if x != "NA" else "NA" for x in t) + f',"{q}"\n')
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write('"","input","filtered","denoisedF","denoisedR","merged","nonchim"\n')
        for i, s in enumerate(S):
            fh.write(f'"{s}",1000,950,940,930,900,{sum(c[i] for c in COUNTS)}\n')
    (d / "summary.json").write_text(json.dumps({"pipeline": "its", "marker": "ITS2"}))
    (d / "run_params.json").write_text(json.dumps({"marker": "ITS2", "sequencerType": "illumina"}))
    return d


def make_pacbio(d: Path, old: bool) -> Path:
    """A PacBio run: old = before v2.17 (one row per SAMPLE, a column per ASV id, taxonomy keyed by id)."""
    d.mkdir(parents=True)
    ids = [f"ASV{i}" for i in range(1, len(SEQ) + 1)]
    keys = ids if old else SEQ
    with open(d / "asv_table.csv", "w") as fh:
        if old:
            fh.write("sample_id," + ",".join(ids) + "\n")
            for i, s in enumerate(S):
                fh.write(s + "," + ",".join(str(c[i]) for c in COUNTS) + "\n")
        else:
            fh.write(",".join(f'"{s}"' for s in S) + ',"sequence"\n')
            for q, c in zip(SEQ, COUNTS):
                fh.write(",".join(map(str, c)) + f',"{q}"\n')
    tax = [["Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Klebsiella", "pneumoniae"],
           ["Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus", "NA"],
           ["Bacteria", "Actinomycetota", "Actinomycetes", "Micrococcales", "Microbacteriaceae", "NA", "NA"],
           ["Bacteria", "NA", "NA", "NA", "NA", "NA", "NA"]]
    with open(d / "taxonomy_table.csv", "w") as fh:
        fh.write('"","Kingdom","Phylum","Class","Order","Family","Genus","Species"\n')
        for k, t in zip(keys, tax):
            fh.write(f'"{k}",' + ",".join(f'"{x}"' if x != "NA" else "NA" for x in t) + "\n")
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write('"","input","filtered","denoised","nonchim"\n')
        for i, s in enumerate(S):
            fh.write(f'"{s}",1000,900,880,{sum(c[i] for c in COUNTS)}\n')
    (d / "summary.json").write_text(json.dumps({"pipeline": "pacbio", "marker": "PACBIO"}))
    (d / "run_params.json").write_text(json.dumps({"marker": "PACBIO"}))
    return d


def main():
    print("Download results: ITS Illumina and PacBio (8-folder result)")
    try:
        import amp_data as AD
        import packages as PK
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            its = AD.load_run(make_its(td / "its"))
            g = {AD.label(a, "Genus") for a in its.asvs}
            sp = {AD.label(a, "Species") for a in its.asvs}
            check("ITS: engine its", its.engine == "its", its.engine)
            check("ITS (old run): UNITE prefixes dropped, 'unidentified' blank, species as binomial",
                  "Aspergillus" in g and "Polyporaceae unidentified" in g and "Aspergillus niger" in sp
                  and "Fusarium unidentified" in sp and not any("__" in x for x in g | sp), (g, sp))
            m = PK.build(its.dir, None, "amplicon", S, False, td / "pk_its", job_label="ITSjob",
                         formats=("html",), tools={}, combined=True)
            z = zipfile.ZipFile(m["files"][0]["path"])
            n = set(z.namelist())
            r = "ITSjob_results/"
            ct = z.read(r + "4.ASV/ASV_count_table_named.csv").decode()
            check("ITS: one combined result with 3.DADA2 / 4.ASV, taxa, alpha, beta, report",
                  all(r + x in n for x in ("4.ASV/ASV_count_table_named.csv",
                                           "5.Taxonomy_Identification/csv/Genus_relative_abundance_per_sample_wide.csv",
                                           "6.Diversity/Alpha/alpha_diversity.csv",
                                           "6.Diversity/Beta/bray_curtis_distance_matrix.csv",
                                           "Report_ITSjob_Illumina_report.html")), sorted(n)[:30])
            check("ITS: no UNITE prefix in the count table", "__" not in ct and "Aspergillus niger" in ct, ct[:300])
            for old in (True, False):
                tag = "old" if old else "new"
                pb = AD.load_run(make_pacbio(td / f"pb_{tag}", old))
                check(f"PacBio ({tag} table): engine pacbio, 3 samples, 4 ASVs named",
                      pb.engine == "pacbio" and pb.samples == S and len(pb.asvs) == 4
                      and {AD.label(a, "Genus") for a in pb.asvs} >= {"Klebsiella", "Bacillus"}
                      and sum(a["counts"]["F1"] for a in pb.asvs) == 820,
                      (pb.engine, pb.samples, [(a.get("id"), a["counts"], a["tax"].get("Genus")) for a in pb.asvs]))
                m = PK.build(pb.dir, None, "amplicon", S, False, td / f"pk_pb_{tag}", job_label="PB",
                             formats=("html",), tools={}, combined=True)
                z = zipfile.ZipFile(m["files"][0]["path"])
                n = set(z.namelist())
                w = list(csv.reader(z.read("PB_results/5.Taxonomy_Identification/csv/Genus_relative_abundance_per_sample_wide.csv")
                                    .decode().splitlines()))
                sums = [round(sum(float(x[i]) for x in w[1:]), 2) for i in range(1, 4)]
                check(f"PacBio ({tag} table): combined result, taxa sum to 100, report",
                      "PB_results/4.ASV/ASV_count_table_named.csv" in n and "PB_results/Report_PB_PacBio_report.html" in n
                      and sums == [100.0, 100.0, 100.0], (sums, sorted(n)[:20]))
                fa = [x for x in n if x.endswith(".fasta")]
                check(f"PacBio ({tag} table): ASV sequences written only when the run has them",
                      bool(fa) != old, fa)
            mm = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            check("Download results accepts ITS Illumina and PacBio; the customer report does not",
                  '_PKG_ONLY_MARKERS = ("ITS", "ITS1", "ITS2", "PACBIO")' in mm
                  and "_custom_check(job_id, packages=True)" in mm
                  and '_CUSTOM_MARKERS = ("12S", "16S", "ONT-16S", "ONT16S", "ONT", "18S", "18S-NEMA", "COX1", "COI")' in mm)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
