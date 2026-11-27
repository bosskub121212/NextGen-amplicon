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
            with open(d / "asv_table.csv", "a") as fh:      # Emu's own row of reads it could not place
                fh.write("unassigned," + ",".join(str(50 + i) for i in range(len(S))) + "\n")
            (d / "metadata.csv").write_text("SampleID,Type\n" + "".join(f"{x},{'AB'[i % 2]}\n" for i, x in enumerate(S)))
            m = PK.build(d, None, "amplicon", S, False, td / "pk", job_label="U2BI-184_fastq", formats=("html",),
                         tools={}, combined=True)
            check("one zip for all samples", len(m["files"]) == 1 and m["files"][0]["name"] == "U2BI-184_fastq_results.zip"
                  and m["combined"] is True and m["files"][0]["samples"] == S, m["files"])
            z = zipfile.ZipFile(m["files"][0]["path"])
            names = set(z.namelist())
            r = "U2BI-184_fastq_results/"
            try:                                   # PNG needs cairosvg (setup_amplicon_report.sh); else the SVG stays
                import cairosvg  # noqa: F401
                PIC = ".png"
            except Exception:
                PIC = ".svg"
                print("  (cairosvg not installed here: Box_plot / Krona picture expected as SVG)")
            want = ["2.Filter/read_tracking_summary.csv", "2.Filter/Boxplot_read_counts_pipeline_steps.svg",
                    "2.Filter/" + S[0] + "/Box_plot" + PIC, "2.Filter/" + S[0] + "/read_tracking.tsv",
                    "4.Taxa_table/Taxa_table_count_table_named.csv", "4.Taxa_table/Taxa_table_counts_per_sample.csv",
                    "5.Taxonomy_Identification/csv/Genus_relative_abundance_per_sample_wide.csv",
                    "5.Taxonomy_Identification/csv/Genus_Top30_relative_abundance_per_sample_long.csv",
                    "5.Taxonomy_Identification/csv/Species_all_taxa_total_relative_abundance.csv",
                    "5.Taxonomy_Identification/Taxa_barplot/Genus_Top30_relative_abundance_per_sample.svg",
                    "5.Taxonomy_Identification/Heatmap/Genus_heatmap_top30.svg",
                    "5.Taxonomy_Identification/Krona/" + S[0] + "_krona_static" + PIC,
                    "6.Diversity/Alpha/alpha_diversity.csv", "6.Diversity/Alpha/Rarefaction_curve_data.csv",
                    "6.Diversity/Alpha/alpha_diversity_by_group.svg",
                    "6.Diversity/Beta/bray_curtis_distance_matrix.csv", "6.Diversity/Beta/PCoA_Bray_scores.csv",
                    "6.Diversity/Beta/PERMANOVA_Type.csv", "6.Diversity/Beta/PCoA_Bray.svg",
                    "Report_U2BI-184_fastq_ONT_report.html"]
            miss = [w for w in want if r + w not in names]
            check("BI team's folders: 2.Filter (tracking, per sample), taxa table, 5.Taxonomy_Identification, alpha, "
                  "beta (PCoA, PERMANOVA), report", not miss, miss)
            top = sorted({x[len(r):].split("/")[0] for x in names if x.startswith(r)})
            check("exactly the BI team's top level: 1.QC ... 6.Diversity + report (no README, Excel, 8.Reads, metadata)",
                  [x for x in top if x != "1.QC"] == ["2.Filter", "4.Taxa_table", "5.Taxonomy_Identification",
                                                     "6.Diversity", "Report_U2BI-184_fastq_ONT_report.html"], top)
            wide = list(csv.reader(z.read(r + "5.Taxonomy_Identification/csv/Genus_relative_abundance_per_sample_wide.csv").decode().splitlines()))
            check("taxa tables: a column per sample", wide[0] == ["Genus"] + S, wide[0])
            sums = [sum(float(x[i]) for x in wide[1:]) for i in range(1, len(S) + 1)]
            check("relative abundance: every sample sums to 100 % (Emu's unassigned row included)",
                  all(abs(v - 100) < 0.05 for v in sums), sums)
            ct = list(csv.reader(z.read(r + "4.Taxa_table/Taxa_table_count_table_named.csv").decode().splitlines()))
            un = [x for x in ct if "unassigned" in x]
            check("Emu count table: tax_id column, the unassigned row named", ct[0][1] == "Emu_tax_id"
                  and un and un[0][2] == "Unassigned", (ct[0], un))
            al = list(csv.reader(z.read(r + "6.Diversity/Alpha/alpha_diversity.csv").decode().splitlines()))
            check("alpha table carries the group", al[0][:2] == ["Sample", "Group"] and al[1][1] == "A", al[:2])
            rep = z.read(r + "Report_U2BI-184_fastq_ONT_report.html").decode()
            check("one report holding every sample", all(x in rep for x in S) and "Analysis by group" in rep)
            import re as _re, html as _h
            fol = [_h.unescape(x) for x in _re.findall(r'<div class="folder">Folder: (.*?)</div>', rep)]
            paths = [p_.strip() for f_ in fol for p_ in f_.split("·")]
            check("report 'Folder:' lines name this result's files (no r_tables/), and they exist",
                  fol and "r_tables/" not in rep.split("References")[0] and all(r + p_ in names for p_ in paths),
                  [p_ for p_ in paths if r + p_ not in names] or fol[:3])
            svg = z.read(r + "5.Taxonomy_Identification/Taxa_barplot/Order_Top30_relative_abundance_per_sample.svg").decode()
            body = svg[svg.index("<svg"):].split(">", 1)[1]          # the drawing, without the size attributes
            check("Taxa_barplot chart is the report's figure (same drawing)", body[:4000] in rep, body[:200])
            check("no R copies of the report's figures (other colours)",
                  not any("pipeline_figures/" in x or "/04_taxonomy" in x or "/04b_taxonomy" in x for x in names))
            d2 = TA.make_amp_run(td / "d2")
            m2 = PK.build(d2, None, "amplicon", None, False, td / "pk2", job_label="Fish", formats=("html",), tools={},
                          combined=True)
            n2 = set(zipfile.ZipFile(m2["files"][0]["path"]).namelist())
            check("DADA2 run: ASV folder, cutadapt summary, ASV sequences",
                  "Fish_results/4.ASV/ASV_count_table_named.csv" in n2 and "Fish_results/2.Cutadapt/cutadapt_summary.csv" in n2
                  and any(x.startswith("Fish_results/2.Cutadapt/") and x.endswith(".cutadapt.json") for x in n2)
                  and any(x.startswith("Fish_results/4.ASV/") and x.endswith(".fasta") for x in n2), sorted(n2)[:20])
            m3 = PK.build(d, None, "amplicon", S[:2], True, td / "pk3", job_label="x", formats=("html",), tools={})
            check("per-sample packages still there (one zip per sample)", len(m3["files"]) == 2)
            z3 = zipfile.ZipFile(m3["files"][0]["path"])
            n3 = z3.namelist()
            check("one-sample result: no README, no 2.Cutadapt and no Beta for Emu",
                  not any(x.endswith("README.txt") or "/2.Cutadapt/" in x or "/Beta/" in x for x in n3), sorted(n3)[:20])
            check("per-sample build: the expected one-sample beta note is not listed under every sample",
                  not any("one sample: no beta" in w for w in m3.get("warnings", [])), m3.get("warnings"))
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
