#!/usr/bin/env python3
"""Run report of a many-sample run (coral, 89 samples): the sample-wise figures go on landscape
pages (<section class="wide">), at most report_builder.CHUNK samples to a page, with readable
labels; the Word file gets landscape sections; small runs are unchanged."""
from __future__ import annotations
import json, random, sys, tempfile, traceback, re
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def make_run(d: Path, n: int, n_taxa: int = 40) -> Path:
    rnd = random.Random(3)
    d.mkdir(parents=True, exist_ok=True)
    S = [f"D1K{1 + i // 13}_{'CSGW'[i % 4]}_H_R{1 + i % 3}x{i}" for i in range(n)]
    taxa = [(str(100 + i), ["Bacteria", f"Phylum{i % 7}", f"Class{i % 9}", f"Order{i % 11}", f"Family{i % 13}",
                            f"Genus{i}", f"Genus{i} sp{i}"]) for i in range(n_taxa)]
    with open(d / "asv_table.csv", "w") as fh:
        fh.write("," + ",".join(S) + "\n")
        for tid, _ in taxa:
            fh.write(tid + "," + ",".join(str(rnd.randint(0, 900)) for _ in S) + "\n")
    with open(d / "taxonomy.csv", "w") as fh:
        fh.write("ASV,Kingdom,Phylum,Class,Order,Family,Genus,Species\n")
        for tid, t in taxa:
            fh.write(tid + "," + ",".join(t) + "\n")
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write("sample,input,qc_filtered,primer_trimmed,final\n")
        for s in S:
            fh.write(f"{s},30000,29000,28000,28000\n")
    (d / "summary.json").write_text(json.dumps({"marker": "ONT-16S", "platform": "ONT", "pipeline": "emu",
        "region": "V1-V9", "db_path": "/db/emu_db_mar2026", "emu_min_abundance": 0.0001,
        "emu_counts": "estimated_reads", "emu_chimera_filter": False}))
    (d / "run_params.json").write_text(json.dumps({"marker": "ONT-16S", "primer_f": "AGRGTTYGATYMTGGCTCAG",
        "primer_r": "RGYTACCTTGTTACGACTT", "ont_region": "V1-V9"}))
    return d


def main():
    print("Run report with many samples (landscape pages)")
    try:
        import report_builder as RB
        check("chunks: 89 -> 45 + 44, 45 -> one part", [len(c) for c in RB._chunks(range(89))] == [45, 44]
              and [len(c) for c in RB._chunks(range(45))] == [45])
        with tempfile.TemporaryDirectory() as td:
            h = RB.build_profile_html(RB.Run(make_run(Path(td) / "many", 60), set()))
            check("stacked bars and heatmap split in landscape parts",
                  h.count('<section class="wide">') >= 6 * 2 + 2 and "Top 30 Genus (Relative Abundance) &mdash; part 2 of 2" in h
                  and "Heatmap visualization of dominant genera &mdash; part 2 of 2" in h, h.count('class="wide"'))
            check("no archaeal reads -> no Archaeal community section", "Archaeal community" not in h)
            m_ = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            check("report request can leave the archaea section out",
                  "archaea:   bool = True" in m_ and "archaea=req.archaea" in m_)
            check("landscape page rule in the CSS", "@page wide { size: A4 landscape;" in h and "section.wide { page: wide; }" in h)
            check("sample labels vertical in the wide charts", 'transform="rotate(-90' in h)
            check("alpha on two landscape pages, PCoA on one", "Alpha diversity &mdash; Observed and Chao1" in h
                  and "Alpha diversity &mdash; Shannon and Simpson" in h and "Beta diversity &mdash; PCoA" in h)
            check("60 x 60 Bray-Curtis matrix not printed, said so", "too large to print legibly" in h
                  and "Bray&ndash;Curtis dissimilarity between samples; darker" not in h)
            steps = re.search(r"<svg[^>]*>(?:(?!</svg>).)*Pipeline step(?:(?!</svg>).)*</svg>", h, re.S)
            check("read-step chart: no 60-name legend", steps and "60 samples, one grey" in steps.group(0))
            # groups set after the run (Customise charts writes metadata.csv; any grouping column)
            d = make_run(Path(td) / "grp", 60)
            S = (d / "asv_table.csv").read_text().splitlines()[0].split(",")[1:]
            (d / "metadata.csv").write_text("SampleID,Type\n" + "".join(f"{x},{x.split('_')[1]}\n" for x in S))
            (d / "run_params.json").write_text(json.dumps({**json.loads((d / "run_params.json").read_text()),
                "metadata": [{"sampleId": x, "group": ""} for x in S]}))
            run = RB.Run(d, set())
            an = RB.Analysis(run) if hasattr(RB, "Analysis") else None
            gh = RB.build_profile_html(run)
            check("groups read from metadata.csv column 'Type' (over run_params)",
                  "grouped by the column <b>Type</b>" in gh and "Analysis by group" in gh, "")
            check("taxa, alpha (Kruskal-Wallis) and beta (PERMANOVA) by group",
                  "Taxa by group" in gh and "The ten most abundant genera of each group" in gh
                  and "Kruskal&ndash;Wallis" in gh and "PERMANOVA on Bray&ndash;Curtis dissimilarity (999 permutations)" in gh
                  and "mean relative abundance per group" in gh)
            sub = RB.samples_in_groups(d, "Type", ["C", "W"])
            rs = RB.Run(d, set(), "grp \u2014 C, W", samples=sub, group_col="Type", groups=["W", "C"])
            sh = RB.build_profile_html(rs)
            check("report of two groups: only their samples, groups in the order picked",
                  len(rs.samples) == len(sub) and len(rs.all_samples) == 60 and 0 < len(sub) < 60
                  and f"This report covers {len(sub)} of the run&rsquo;s 60 samples" in sh
                  and sh.index("<b>W</b>") < sh.index("<b>C</b>") and "Analysis by group" in sh
                  and all(x.split("_")[1] in ("C", "W") for x in rs.samples), len(sub))
            try:
                RB.Run(d, set(), samples=["nope"])
                check("no chosen sample in the run -> error", False)
            except ValueError:
                check("no chosen sample in the run -> error", True)
            m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            h_ = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            check("report request carries the grouping and the ticked groups",
                  "group_col: Optional[str]" in m and "rb.samples_in_groups(" in m
                  and 'id="repGcol"' in h_ and "group_col: oneRun ? R.gcol : null" in h_)
            check("wide bars ordered by group, with group brackets", "ordered by group" in gh
                  and 'font-weight="700" text-anchor="middle">C<' in gh)
        with tempfile.TemporaryDirectory() as td:          # 18S V4 against PR2
            d = Path(td) / "pr2"
            d.mkdir()
            S = ["KK-NE-1", "KK-NE-2", "KK-S-1"]
            L = [("Eukaryota", "TSAR", "Stramenopiles", "Gyrista", "Bacillariophyceae", "Bacillariophyceae_X",
                  "Naviculaceae", "Navicula", "Navicula_sp."),
                 ("Eukaryota", "Obazoa", "Opisthokonta", "Metazoa", "Arthropoda", "Crustacea", "Calanidae",
                  "Calanus", "Calanus_finmarchicus")]
            seqs = ["ACGT" * 95, "TTGCA" * 76]
            (d / "asv_table.csv").write_text("sequence," + ",".join(S) + "\n"
                                             + "".join(q + ",100,200,300\n" for q in seqs))
            (d / "taxonomy_table.csv").write_text("ASV,Kingdom,Supergroup,Phylum,Subdivision,Class,Order,Family,Genus,Species\n"
                                                  + "".join(q + "," + ",".join(l) + "\n" for q, l in zip(seqs, L)))
            (d / "read_tracking.csv").write_text("sample,input,filtered,nonchim\n" + "".join(f"{x},900,800,700\n" for x in S))
            pj = {"marker": "18S-nema", "dbPath": "/db/NemaBase/pr2_version_5.1.1_SSU_dada2.fasta.gz",
                  "primer_f": "CCAGCASCYGCGGTAATTCC", "primer_r": "ACTTTCGTTCTTGATYRA"}
            (d / "run_params.json").write_text(json.dumps(pj))
            (d / "summary.json").write_text(json.dumps({"marker": "18S-nema"}))
            ph = RB.build_profile_html(RB.Run(d, set()))
            check("18S PR2 run report: eukaryote title, PR2 ranks, '_X' blanked, species as binomials",
                  "18S rRNA V4 Eukaryote Amplicon Sequencing Report" in ph and "Division (reported here as Phylum)" in ph
                  and "Bacillariophyceae_X" not in ph and "Calanus finmarchicus" in ph and "16S" not in
                  re.sub(r"<[^>]+>", " ", ph).split("References")[0])
            import amp_data as AD, amp_package as AP, custom_report as CR
            r = AD.load_run(d)
            t, ref = CR._tax_text(r)
            refs = CR._refs(r, t + " " + ref)
            check("18S PR2 package: label, PR2 text and reference; no UNOISE / UNITE / Vu",
                  AP.marker_label(r)[0] == "18S rRNA (eukaryotes)" and "PR2 5.1.1" in ref
                  and any("(PR2)" in x for x in refs) and not any(k in x for x in refs for k in ("UNOISE2", "UNITE", "Vu D"))
                  and AD.label(next(a for a in r.asvs if a["tax"].get("Genus") == "Navicula"), "Species") == "Navicula unidentified",
                  (AP.marker_label(r), ref, refs))
        k = RB._kruskal([[0.92, 1.1, 3.4, 2.2, 7.7, 5.1, 6.4, 9.0], [3.5, 12.1, 9.9, 4.4, 8.8, 10.2, 7.3, 11.0, 6.1, 5.5],
                         [1, 1, 2, 2, 3]])
        check("Kruskal-Wallis with ties: H and p as scipy.stats.kruskal",
              k and abs(k[0] - 11.0456973) < 1e-5 and k[1] == 2 and abs(k[2] - 0.0039944529) < 1e-8, k)
        check("chi-square tail (incomplete gamma) as scipy chi2.sf(4.6, 3)", abs(RB._gammaincc(1.5, 2.3) - 0.2035421) < 1e-6)
        with tempfile.TemporaryDirectory() as td:
            small = RB.build_profile_html(RB.Run(make_run(Path(td) / "few", 6), set()))
            check("a 6-sample run keeps its portrait figures", 'class="wide"' not in small
                  and "Bray&ndash;Curtis dissimilarity between samples; darker" in small)
            try:
                import docx, docx_export as DX
                p = Path(td) / "r.docx"
                DX.html_to_docx(h, p)
                secs = docx.Document(str(p)).sections
                o = ["L" if s_.page_width > s_.page_height else "P" for s_ in secs]
                check("Word: landscape sections for the wide parts, portrait around them",
                      o[0] == "P" and "L" in o and o[-1] == "P" and all(a != b for a, b in zip(o, o[1:])), o)
            except ImportError:
                print("  skip  Word check (python-docx / cairosvg not installed)")
    except Exception:
        check("many-sample report", False, traceback.format_exc(limit=5))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
