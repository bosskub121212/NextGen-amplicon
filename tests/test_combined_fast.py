#!/usr/bin/env python3
"""Download results -> One combined result: the chosen samples of a run in ONE result laid out like the
BI team's deliverable (QC, ASV / taxa table, taxa per rank with a column per sample, alpha, beta with
PCoA and PERMANOVA, Excel, one report)."""
from __future__ import annotations
import csv, json, sys, tempfile, traceback, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend")); sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def main():
    print("Download results: one combined result")
    try:
        import packages as PK
        import test_report_many_fast as TM
        import test_amp_fast as TA
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            d = TM.make_run(td / "emu", 6)
            S = (d / "asv_table.csv").read_text().splitlines()[0].split(",")[1:]
            (d / "metadata.csv").write_text("SampleID,Type\n" + "".join(f"{x},{'AB'[i % 2]}\n" for i, x in enumerate(S)))
            m = PK.build(d, None, "amplicon", S, False, td / "pk", job_label="U2BI-184_fastq", formats=("html",),
                         tools={}, combined=True)
            check("one zip for all samples", len(m["files"]) == 1 and m["files"][0]["name"] == "U2BI-184_fastq_results.zip"
                  and m["combined"] is True and m["files"][0]["samples"] == S, m["files"])
            z = zipfile.ZipFile(m["files"][0]["path"])
            names = set(z.namelist())
            r = "U2BI-184_fastq_results/"
            want = ["1.QC/read_tracking_summary.csv", "1.QC/Boxplot_read_counts_pipeline_steps.svg",
                    "4.Taxa_table/Taxa_table_count_table_named.csv", "4.Taxa_table/Taxa_table_counts_per_sample.csv",
                    "5.Taxa/csv/Genus_relative_abundance_per_sample_wide.csv",
                    "5.Taxa/csv/Genus_Top30_relative_abundance_per_sample_long.csv",
                    "5.Taxa/csv/Species_all_taxa_total_relative_abundance.csv",
                    "5.Taxa/png/Genus_Top30_relative_abundance_per_sample.svg", "5.Taxa/png/Genus_heatmap_top30.svg",
                    "6.Diversity/Alpha/alpha_diversity.csv", "6.Diversity/Alpha/Rarefaction_curve_data.csv",
                    "6.Diversity/Alpha/alpha_diversity_by_group.svg",
                    "6.Diversity/Beta/bray_curtis_distance_matrix.csv", "6.Diversity/Beta/PCoA_Bray_scores.csv",
                    "6.Diversity/Beta/PERMANOVA_Type.csv", "6.Diversity/Beta/PCoA_Bray.svg",
                    "7.Tables/U2BI-184_fastq_results.xlsx", "Report_U2BI-184_fastq_results.html", "README.txt"]
            try:                                   # the Excel file needs openpyxl (setup_amplicon_report.sh)
                import openpyxl  # noqa: F401
            except ImportError:
                want.remove("7.Tables/U2BI-184_fastq_results.xlsx")
                print("  (openpyxl not installed here: no Excel file expected)")
            miss = [w for w in want if r + w not in names]
            check("BI-style folders: QC, taxa table, taxa per rank, alpha, beta (PCoA, PERMANOVA), Excel, report",
                  not miss, miss)
            wide = list(csv.reader(z.read(r + "5.Taxa/csv/Genus_relative_abundance_per_sample_wide.csv").decode().splitlines()))
            check("taxa tables: a column per sample", wide[0] == ["Genus"] + S, wide[0])
            al = list(csv.reader(z.read(r + "6.Diversity/Alpha/alpha_diversity.csv").decode().splitlines()))
            check("alpha table carries the group", al[0][:2] == ["Sample", "Group"] and al[1][1] == "A", al[:2])
            rep = z.read(r + "Report_U2BI-184_fastq_results.html").decode()
            check("one report holding every sample", all(x in rep for x in S) and "Analysis by group" in rep)
            d2 = TA.make_amp_run(td / "d2")
            m2 = PK.build(d2, None, "amplicon", None, False, td / "pk2", job_label="Fish", formats=("html",), tools={},
                          combined=True)
            n2 = set(zipfile.ZipFile(m2["files"][0]["path"]).namelist())
            check("DADA2 run: ASV folder, cutadapt summary, ASV sequences",
                  "Fish_results/4.ASV/ASV_count_table_named.csv" in n2 and "Fish_results/2.Cutadapt/cutadapt_summary.csv" in n2
                  and any(x.startswith("Fish_results/4.ASV/") and x.endswith(".fasta") for x in n2), sorted(n2)[:20])
            m3 = PK.build(d, None, "amplicon", S[:2], True, td / "pk3", job_label="x", formats=("html",), tools={})
            check("per-sample packages still there (one zip per sample)", len(m3["files"]) == 2)
            fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            mm = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            check("Download results offers the combined result first", 'value="combined"' in fe
                  and "combined: !!opt.combined" in fe and "combined: bool = False" in mm
                  and 'combined=bool(req.combined and kind == "amplicon")' in mm)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")
    return FAILS


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
