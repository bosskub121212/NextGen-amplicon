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
            check("landscape page rule in the CSS", "@page wide { size: A4 landscape;" in h and "section.wide { page: wide; }" in h)
            check("sample labels vertical in the wide charts", 'transform="rotate(-90' in h)
            check("alpha on two landscape pages, PCoA on one", "Alpha diversity &mdash; Observed and Chao1" in h
                  and "Alpha diversity &mdash; Shannon and Simpson" in h and "Beta diversity &mdash; PCoA" in h)
            check("60 x 60 Bray-Curtis matrix not printed, said so", "too large to print legibly" in h
                  and "Bray&ndash;Curtis dissimilarity between samples; darker" not in h)
            steps = re.search(r"<svg[^>]*>(?:(?!</svg>).)*Pipeline step(?:(?!</svg>).)*</svg>", h, re.S)
            check("read-step chart: no 60-name legend", steps and "60 samples, one grey" in steps.group(0))
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
