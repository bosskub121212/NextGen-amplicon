#!/usr/bin/env python3
"""Download results: a PNG beside every PDF figure; NanoPlot in the combined result too; the run report
without the 'Reagent-associated genera' table (BI team, v2.42.0)."""
from __future__ import annotations
import gzip, os, shutil, stat, sys, tempfile, traceback, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(ROOT / "tests"))
FAILS = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def tiny_pdf(path: Path, pages: int = 1):
    """A valid PDF of `pages` pages, each a filled square (no plotting library needed)."""
    objs = ["<< /Type /Catalog /Pages 2 0 R >>",
            "<< /Type /Pages /Kids [" + " ".join(f"{3 + 2 * i} 0 R" for i in range(pages)) + f"] /Count {pages} >>"]
    for i in range(pages):
        c = b"0.2 0.4 0.6 rg 20 20 160 160 re f"
        objs.append(f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 200 200] /Contents {4 + 2 * i} 0 R >>")
        objs.append(f"<< /Length {len(c)} >>\nstream\n{c.decode()}\nendstream")
    out, offs = b"%PDF-1.4\n", []
    for k, o in enumerate(objs, 1):
        offs.append(len(out))
        out += f"{k} 0 obj\n{o}\nendobj\n".encode()
    x = len(out)
    out += f"xref\n0 {len(objs) + 1}\n0000000000 65535 f \n".encode()
    out += "".join(f"{o:010d} 00000 n \n" for o in offs).encode()
    out += f"trailer\n<< /Size {len(objs) + 1} /Root 1 0 R >>\nstartxref\n{x}\n%%EOF\n".encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(out)


def main():
    print("PNG of the PDF figures, NanoPlot in every result, no reagent table")
    try:
        import amp_combined as AC
        import report_builder as RB
        tool = AC._pdf_tool()
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            # --- PNG beside each PDF ----------------------------------------------------------------
            r = t / "res"
            tiny_pdf(r / "6.Diversity/Beta/figures/12_upgma_bray.pdf")
            tiny_pdf(r / "1.QC/QualityProfile.pdf", pages=3)
            tiny_pdf(r / "Report_x.pdf")
            tiny_pdf(r / "2.Filter/S/filtered_reads/x.pdf")
            notes = []
            n = AC.pngs_for_pdfs(r, notes)
            if tool:
                check(f"PNG beside each PDF figure ({tool}); a page each for a multi-page PDF",
                      n == 2 and (r / "6.Diversity/Beta/figures/12_upgma_bray.png").read_bytes()[:4] == b"\x89PNG"
                      and all((r / f"1.QC/QualityProfile_p{i}.png").exists() for i in (1, 2, 3)), (n, notes))
                check("the report and the reads left alone", not (r / "Report_x.png").exists()
                      and not (r / "2.Filter/S/filtered_reads/x.png").exists())
                check("a PDF that has its PNG is not done again", AC.pngs_for_pdfs(r, []) == 0)
            else:
                check("no PDF renderer: said in the build notes", n == 0 and notes and "pymupdf" in notes[0], notes)
            sh = (ROOT / "setup_amplicon_report.sh").read_text(encoding="utf-8")
            check("setup_amplicon_report.sh installs PyMuPDF", "pip install --quiet pymupdf" in sh)

            # --- combined result: NanoPlot of each sample, PNG of the pipeline's PDFs ----------------
            import packages as PK
            from test_amp_fast import make_emu_run
            d = make_emu_run(t / "emu")
            S = (d / "asv_table.csv").read_text().splitlines()[0].split(",")[1:]
            for s in S:
                (d / "merged_reads").mkdir(exist_ok=True)
                with gzip.open(d / "merged_reads" / f"{s}.fastq.gz", "wt") as fh:
                    fh.write("@r1\nACGTACGT\n+\nIIIIIIII\n")
            tiny_pdf(d / "r_plots" / "01_read_counts.pdf")
            nano = t / "bin" / "NanoPlot"
            nano.parent.mkdir()
            nano.write_text("#!/usr/bin/env python3\nimport sys, os\na = sys.argv[1:]\n"
                            "if '--version' in a:\n    print('NanoPlot 1.42.0'); sys.exit(0)\n"
                            "o = a[a.index('-o') + 1]; os.makedirs(o, exist_ok=True)\n"
                            "open(os.path.join(o, 'NanoStats.txt'), 'w').write('Mean read length: 8\\n')\n"
                            "open(os.path.join(o, 'NanoPlot-report.html'), 'w').write('<html></html>')\n")
            nano.chmod(nano.stat().st_mode | stat.S_IEXEC)
            m = PK.build(d, None, "amplicon", S, False, t / "pk", job_label="J", formats=("html",),
                         tools={"NanoPlot": str(nano)}, combined=True)
            z = zipfile.ZipFile(m["files"][0]["path"])
            names = set(z.namelist())
            check("combined result: NanoPlot of each sample in 1.QC/<sample>/NanoPlot",
                  all(f"J_results/1.QC/{s}/NanoPlot/NanoStats.txt" in names for s in S), sorted(x for x in names if "1.QC" in x))
            if tool:
                check("combined result: the pipeline's PDF figure has its PNG",
                      "J_results/1.QC/01_read_counts.pdf" in names and "J_results/1.QC/01_read_counts.png" in names)
            AC_doc = AC.__doc__
            check("layout doc: NanoPlot under 1.QC, PNG copies said", "NanoPlot" in AC_doc and "PNG beside it" in AC_doc)
            m2 = PK.build(d, None, "amplicon", S[:1], False, t / "pk2", job_label="J", formats=("html",),
                          tools={"NanoPlot": None}, combined=True)
            check("NanoPlot missing: said in the build's notes, not silently left out",
                  any("NanoPlot is not installed" in w for w in m2.get("warnings", [])), m2.get("warnings"))

            # --- the run report: no reagent table, no Salter --------------------------------------
            from test_report_many_fast import make_run
            d3 = make_run(t / "bac", 4)
            html = RB.build_profile_html(RB.Run(d3, RB.KITOME, "x"))
            check("run report: no 'Reagent-associated genera' table, note or Salter reference",
                  "Reagent-associated genera" not in html and "reagent-contaminant list" not in html
                  and "Salter" not in html and "Most abundant" in html)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
