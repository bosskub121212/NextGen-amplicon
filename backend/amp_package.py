#!/usr/bin/env python3
"""Per-sample package of an amplicon run (ONT 16S / Emu, 16S Illumina / DADA2, 12S·CytB) in the
same standard as the ONT WGS package: numbered folders + Report_<S>.pdf / .docx in the package
root, the report laid out like the WGS report (no cover; "<marker> amplicon sequencing report /
Sample: S"; numbered sections, each ending "Results <folder>").

  <S>/
    Report_<S>.pdf / .docx
    1.QcReads/     read tracking table + box plot, FastQC / fastp (Illumina) or NanoPlot (ONT),
                   DADA2 quality profiles, filtered_reads/ (optional)
    2.Taxonomy/    outfiles/ (abundance tables per rank, ASV FASTA, BIOM), Bar charts/, Heatmap/,
                   Krona/, Emu/ (Emu's own per-sample table), BLAST/ (12S·CytB: local + NCBI check)
    3.Diversity/   alpha_diversity.tsv, Rarefaction_curves.png
    4.Tables/      <S>.xlsx
  A folder whose step did not run holds NOT_AVAILABLE.txt.

Built from a one-sample view of the run (amp_data.RunData.only): nothing of another sample —
names, taxa, counts, SV numbers or figures — reaches the package.

    build_sample(run_dir, upload_dir, sample, out_dir, formats, name, tools, include_reads) -> manifest
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

import amp_data as AD
import custom_report as CR
import qc_tools as Q

FOLDERS = ["1.QcReads", "2.Taxonomy", "3.Diversity", "4.Tables"]
e = CR.e

CSS = """
@page { size: A4; margin: 22mm 20mm 18mm 20mm;
  @top-right { content: "Page | " counter(page); font-size: 9pt; color: #555; font-family: Carlito, Calibri, "DejaVu Sans", sans-serif; } }
body { font-family: Carlito, Calibri, "DejaVu Sans", sans-serif; font-size: 11pt; line-height: 1.45; color: #111; }
h1 { font-size: 20pt; text-align: center; font-weight: normal; margin: 0 0 8mm; }
.sample { font-size: 12pt; font-weight: bold; margin-bottom: 2mm; } .marker { margin-bottom: 5mm; }
h2 { font-size: 12.5pt; margin: 0 0 3mm; page-break-before: always; } h2.first { page-break-before: avoid; }
h3 { font-size: 11.5pt; margin: 4mm 0 2mm; }
h4 { font-size: 11pt; margin: 3mm 0 1.5mm; text-decoration: underline; font-weight: normal; }
p { margin: 0 0 3mm; text-align: justify; }
p.res { font-weight: bold; font-size: 12pt; margin-top: 4mm; page-break-before: avoid; }
p.cap { text-align: center; font-size: 9pt; color: #444; font-style: italic; }
p.note { color: #b00000; text-align: center; font-style: italic; }
p.nofig { text-align: center; color: #777; border: 1px dashed #bbb; padding: 12pt; }
div.fig { text-align: center; margin: 2mm 0; page-break-inside: avoid; page-break-before: avoid; }
div.keep { page-break-inside: avoid; } p.cap { page-break-before: avoid; }
table { border-collapse: collapse; width: 100%; margin: 2mm 0; font-size: 9pt; table-layout: fixed; }
th, td { border: 1px solid #c8c8c8; padding: 1.5mm 2mm; text-align: left; vertical-align: top; overflow-wrap: anywhere; }
th { background: #f2f2f2; } table.num { font-size: 7.6pt; } table.num td + td + td { text-align: right; }
ul.files { font-size: 9.5pt; list-style: none; padding-left: 0; } .indent { margin-left: 6mm; }
.ul { text-decoration: underline; font-weight: bold; }
"""


def results(folder: str) -> str:
    return f'<p class="res">Results {e(folder)}</p>'


def marker_label(run: AD.RunData) -> tuple[str, str]:
    """(title word, method line) from what the run was."""
    p, sm = run.params, run.summary
    mk = str(p.get("marker") or sm.get("marker") or "").upper()
    if run.kind == "emu":
        return "16S rRNA", "Full-length 16S rRNA gene, Oxford Nanopore reads, Emu"
    if mk == "12S":
        g = "CytB" if str(p.get("mito_gene", "")).upper() == "CYTB" else "12S rRNA"
        return g, f"{g} gene (mitochondrial), Illumina paired reads, DADA2"
    return ("16S rRNA" if mk == "16S" else (mk or "Amplicon")), \
        f"{'16S rRNA gene' if mk == '16S' else (mk or 'amplicon')}, Illumina paired reads, DADA2"


def _put(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _move(src: Path, dst: Path):
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))


def build_sample(run_dir: Path, upload_dir: Path | None, sample: str, out_dir: Path, formats=("pdf", "docx"),
                 name: str | None = None, tools: dict | None = None, include_reads: bool = True) -> dict:
    tools = dict(tools or {})
    full = AD.load_run(Path(run_dir))
    run = full.only(sample)                                # refuses an unknown sample
    s = sample
    name = name or CR.safe(s)
    out_dir = Path(out_dir)
    root = out_dir / name
    shutil.rmtree(root, ignore_errors=True)
    F = {f: root / f for f in FOLDERS}
    for d in F.values():
        d.mkdir(parents=True)
    work = out_dir / f".{name}_work"
    shutil.rmtree(work, ignore_errors=True)
    img = work / "img"
    for sub in ("Taxonomic profiles", "Heatmap", "Krona plot", "Rarefaction curves"):
        (img / sub).mkdir(parents=True, exist_ok=True)
    emu = run.kind == "emu"
    unit = ("taxon" if emu else "ASV") if run.counts_are_reads else "ppm"
    notes = []

    # 1. reads: QC tools on the raw reads, filtered reads
    raw = CR._find_reads(upload_dir, s, run.peer_names) if not emu else ()
    filt = CR._filtered(run, s)
    reads_rel = None
    if include_reads and filt:
        for q in filt:
            _put(q, F["1.QcReads"] / "filtered_reads" / q.name)
        reads_rel = "1.QcReads/filtered_reads"
    qc = {}
    if emu:
        fs = CR._ont_raw(run, upload_dir, s)
        if fs:
            qc["stats"] = Q.read_stats(fs)
            qc["nanoplot"] = Q.nanoplot(fs, F["1.QcReads"] / "NanoPlot", CR._tool(tools, "NanoPlot"))
        else:
            qc["missing"] = True
    else:
        r1, r2 = raw if raw else (None, None)
        if r1:
            qc["fastqc"] = Q.fastqc([x for x in (r1, r2) if x], F["1.QcReads"] / "FastQC", CR._tool(tools, "fastqc"))
            qc["fastp"] = Q.fastp_check(r1, r2, F["1.QcReads"] / "fastp", CR._tool(tools, "fastp"))
        else:
            qc["missing"] = True
    t = run.tracking[s]
    with open(F["1.QcReads"] / "read_tracking.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        cols = [c for c in run.track_cols if c != "perc"]
        w.writerow(["Sample"] + cols + ["perc"])
        w.writerow([s] + ["" if t.get(c) is None else t.get(c) for c in cols] + [t.get("perc", "")])

    # 2. taxonomy tables, Emu's own table, BLAST / NCBI check
    CR._write_outfiles(run, s, F["2.Taxonomy"] / "outfiles")
    ed = Path(run_dir) / "emu_output" / s
    if emu and ed.is_dir():
        for f in sorted(ed.iterdir()):
            if f.is_file() and f.suffix.lower() in (".tsv", ".csv", ".txt") and not f.name.endswith(".log"):
                _put(f, F["2.Taxonomy"] / "Emu" / f.name)
    asvs = run.sample_asvs(s)
    if not emu and any(a.get("best_hit") or a.get("ncbi_best") for a in asvs):
        (F["2.Taxonomy"] / "BLAST").mkdir(exist_ok=True)
        with open(F["2.Taxonomy"] / "BLAST" / "blast_ncbi_check.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["ASV", "Name given", "Local best hit", "Local % identity", "NCBI best hit",
                        "NCBI % identity", "Target group", "Check", "Reads"])
            for a in asvs:
                w.writerow([a["sv"], AD.label(a, "Species"), a.get("best_hit", ""), a.get("pident", ""),
                            a.get("ncbi_best", ""), a.get("ncbi_pident", ""), a.get("target", ""),
                            a.get("check", ""), a["counts"][s]])

    # 4. Excel
    xlsx = F["4.Tables"] / f"{name}.xlsx"
    try:
        import amplicon_excel as AX
        AX.write_sample(xlsx, run, s)
    except ImportError:
        pass

    # figures: R (phyloseq / ggplot2 / microbiomeutilities / vegan / dada2) and Krona
    CR._prepare(run, s, work / "1", raw, filt)
    with open(work / "jobs.tsv", "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow([str(work / "1"), str(img), s, unit])
    rs = CR._tool(tools, "Rscript")
    log = out_dir / f"{name}.plots.log"                    # beside the package, never inside it
    if rs:
        try:
            pr = subprocess.run([rs, str(CR.RSCRIPT), "--jobs", str(work / "jobs.tsv")], capture_output=True,
                                text=True, timeout=3 * 3600)
            log.write_text((pr.stdout or "") + (pr.stderr or ""))
        except subprocess.TimeoutExpired:
            log.write_text("custom_report_plots.R stopped after 3 hours\n")
        except OSError as ex:
            log.write_text(f"custom_report_plots.R could not run: {ex}\n")
    else:
        notes.append("R (Rscript) not found: no figures")
    _move(img / "Box_plot.png", F["1.QcReads"] / "Box_plot.png")
    for f in ("plotQualityProfile.pdf", "plotQualityProfile_filtered.pdf"):
        _move(img / f, F["1.QcReads"] / f)
    for _, _, low in CR.RANK_FILES:
        _move(img / "Taxonomic profiles" / f"plot_bar_{low}.png", F["2.Taxonomy"] / "Bar charts" / f"plot_bar_{low}.png")
    _move(img / "Heatmap" / "Heat_map_label_ASV_number_top50ASVs.png",
          F["2.Taxonomy"] / "Heatmap" / "Heat_map_top50.png")
    _move(img / "Krona plot" / "krona_static.png", F["2.Taxonomy"] / "Krona" / "krona_static.png")
    _move(img / "Rarefaction curves" / "Rarefaction_curves.png", F["3.Diversity"] / "Rarefaction_curves.png")
    kt = CR._tool(tools, "ktImportText")
    kh = F["2.Taxonomy"] / "Krona" / f"{name}-krona.html"
    if kt:
        try:
            subprocess.run([kt, "-o", str(kh), f"{work / '1' / 'krona.tsv'},{s}"], capture_output=True, timeout=600)
        except subprocess.TimeoutExpired:
            pass
    if not kh.exists() and (F["2.Taxonomy"] / "Krona" / "krona_static.png").exists():
        kh.write_text("<!doctype html><meta charset=utf-8><title>" + e(s) + "</title><body style='text-align:center;"
                      "font-family:sans-serif'><h3>" + e(s) + " — taxonomy (KronaTools not installed: still image)</h3>"
                      + CR._fig(F["2.Taxonomy"] / "Krona" / "krona_static.png", 20) + "</body>")
    plots, alpha = {}, None
    try:
        plots = json.loads((work / "1" / "plots_info.json").read_text())
    except Exception:
        pass
    try:
        with open(work / "1" / "alpha.tsv") as fh:
            alpha = next(csv.DictReader(fh, delimiter="\t"), None)
    except OSError:
        pass
    a = alpha or run.alpha.get(s) or {}
    if a:
        with open(F["3.Diversity"] / "alpha_diversity.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            keys = ["Observed", "Chao1", "se.chao1", "Shannon", "Simpson"]
            w.writerow(["Sample"] + keys)
            w.writerow([s] + [("" if (k in ("Chao1", "se.chao1") and not run.counts_are_reads) else a.get(k, ""))
                              for k in keys])

    # report
    doc = build_html(run, s, root, qc, plots, a, reads_rel, xlsx)
    man = {"dir": str(root), "sample": s, "folders": FOLDERS, "notes": notes}
    stem = f"Report_{name}"
    if "html" in formats:
        (root / f"{stem}.html").write_text(doc, encoding="utf-8")
        man["html"] = str(root / f"{stem}.html")
    if "pdf" in formats:
        try:
            from weasyprint import HTML as _WH
            _WH(string=doc, base_url=str(root)).write_pdf(str(root / f"{stem}.pdf"))
            man["pdf"] = str(root / f"{stem}.pdf")
        except Exception as ex:
            man["pdf_error"] = str(ex)
    if "docx" in formats:
        try:
            import docx_export
            docx_export.html_to_docx(doc, root / f"{stem}.docx", title=f"Amplicon sequencing report — {s}")
            man["docx"] = str(root / f"{stem}.docx")
        except Exception as ex:
            man["docx_error"] = str(ex)
    shutil.rmtree(work, ignore_errors=True)
    for f in root.rglob("*.log"):           # tool logs (NanoPlot's) carry server paths: not for the customer
        f.unlink()
    for f in FOLDERS:
        d = root / f
        if not any(d.rglob("*")):
            (d / "NOT_AVAILABLE.txt").write_text("This step was not run for this sample, or its tool was not "
                                                  "installed — see the report.\n")
    return man


def nanostats(p: Path) -> dict:
    """Mean / median read quality from NanoPlot's NanoStats.txt ({} when absent)."""
    out = {}
    try:
        for line in Path(p).read_text(errors="replace").splitlines():
            k, _, v = line.partition(":")
            k = k.strip().lower()
            if k in ("mean read quality", "median read quality"):
                try:
                    out["mean" if k.startswith("mean") else "median"] = float(v.strip().replace(",", ""))
                except ValueError:
                    pass
    except OSError:
        pass
    return out


def build_html(run: AD.RunData, s: str, root: Path, qc: dict, plots: dict, a: dict, reads_rel, xlsx: Path) -> str:
    emu = run.kind == "emu"
    t = run.tracking[s]
    word, line = marker_label(run)
    rel = lambda x: f"1.QcReads/{x}"
    xl = (xlsx, f"4.Tables/{xlsx.name}")
    o = [f"<h1>{e(word)} amplicon sequencing report</h1>", f'<p class="sample">Sample: {e(s)}</p>',
         f'<p class="marker">{e(line)}</p>']

    # 1. QC
    o.append('<h2 class="first">1. Quality control and clean data</h2>')
    o += (CR._emu_steps if emu else CR._dada2_steps)(run, 1, s, "1.QcReads", root / "1.QcReads", qc, plots, "",
                                                     rel=rel, reads_rel=reads_rel, xlsx=xl)
    if t.get("raw_data") is None:
        o.append('<p class="note">raw_data is unknown for this sample: ' + (
            "the pipeline's read tracking has no row for it.</p>" if emu else
            "cutadapt's own read counts were not found, so the count before primer removal cannot be given.</p>"))
    if s in getattr(run, "dropped", []):
        o.append('<p class="note">No reads of this sample were left after filtering: the pipeline dropped it.</p>')
    ns = nanostats(root / "1.QcReads" / "NanoPlot" / "NanoStats.txt")
    if ns.get("mean") is not None:
        o.append("<p>NanoPlot's own summary (1.QcReads/NanoPlot/NanoStats.txt) gives a mean read quality of "
                 f"Q{ns['mean']:.1f}" + (f" and a median of Q{ns['median']:.1f}" if ns.get("median") is not None else "")
                 + ": NanoPlot averages the reads' error rates, which weighs low-quality reads more than the mean of "
                 "the reads' quality values given above. Both describe the same reads.</p>")
    o.append(results("1.QcReads"))

    # 2. Taxonomy
    tax_t, _ = CR._tax_text(run)
    o += ["<h2>2. Taxonomic profiles</h2>", "<h3>Taxonomic assignment</h3>", f"<p>{tax_t}</p>"]
    if t.get("Nonchim") == 0:
        o.append('<p class="note">' + ("No read of this sample was assigned by Emu" if emu else
                                       "No reads of this sample passed filtering, merging and chimera removal")
                 + "; the profiles below are empty. Check the library or sequencing depth.</p>")
    if not emu and (run.dir / "blast" / "taxonomy_blast.csv").exists():
        nt = run.non_target_share(s)
        nts = sorted({AD.label(x, "Species") for x in run.sample_asvs(s) if x.get("target") == "no"})
        o.append("<p><b>Non-target reads:</b> " + (f"{nt:.1f}% of this sample's reads came from outside the limit "
                                                    f"({e(', '.join(nts))}). They are kept in the tables and marked."
                                                    if nt > 0 else "none.") + "</p>")
        chk = [x for x in run.sample_asvs(s) if x.get("check") == "yes"]
        if chk:
            o.append("<p><b>ASVs to check</b> (local reference and NCBI disagree; named only to the shared rank):</p>"
                     + CR._table(["ASV", "Name given", "Local best hit (% id)", "NCBI best hit (% id)", "Reads"],
                                 [[x["sv"], e(AD.label(x, "Species")), f"{e(x.get('best_hit', ''))} ({x.get('pident', '')})",
                                   f"{e(x.get('ncbi_best', ''))} ({x.get('ncbi_pident', '')})", f"{x['counts'][s]:,}"]
                                  for x in chk]))
    T = root / "2.Taxonomy"
    o += ["<h3>Bar charts of taxonomic profiles</h3>",
          f"<p>Relative abundance of the taxa detected in this sample at each rank (phyloseq "
          f"{e(plots.get('phyloseq', ''))}, ggplot2 {e(plots.get('ggplot2', ''))}). "
          + ("Taxa" if emu else "ASVs") + " without a name at a rank are shown as '&lt;higher rank&gt; unidentified'. "
          "At a rank with more than 30 taxa the chart colours the 30 most abundant and sums the rest into one grey "
          "'Others' segment (its title gives the count); the tables in 2.Taxonomy/outfiles and the Excel file list "
          "every taxon.</p>"]
    for rank, _, low in CR.RANK_FILES:              # label and its chart never split over two pages
        o.append(f'<div class="keep"><p><u>At {rank} level</u></p>'
                 + CR._fig(T / "Bar charts" / f"plot_bar_{low}.png", 9) + "</div>")
    unit = "taxa" if emu else "ASVs"
    hm = plots.get("heatmap_tool", "pheatmap")
    o += [f'<div class="keep"><h3>Top 50 {unit} heatmap ({e(hm.split(" ")[0])})</h3>'
          + CR._fig(T / "Heatmap" / "Heat_map_top50.png", 13) + "</div>",
          '<div class="keep"><h3>Krona plot</h3>'
          "<p>Krona shows the taxonomy as nested rings that can be zoomed in a web browser (the HTML file in "
          "2.Taxonomy/Krona); the figure is a still of it.</p>" + CR._fig(T / "Krona" / "krona_static.png", 12)
          + "</div>"]
    files = ([] if emu else ["outfiles/expr.asv.fasta"]) + ["outfiles/expr.abundance.all.txt",
                                                            "outfiles/expr.relative_abundance.all.txt",
                                                            "outfiles/expr.biom"]
    o.append("<p>Tables per rank (counts and relative abundance), "
             + ("" if emu else "the ASV sequences (FASTA), ") + "and a BIOM file are in 2.Taxonomy/outfiles"
             + (", Emu's own per-sample table in 2.Taxonomy/Emu" if (T / "Emu").exists() else "")
             + (", the BLAST / NCBI check of every ASV in 2.Taxonomy/BLAST" if (T / "BLAST").exists() else "")
             + ".</p>")
    if emu:
        o.append("<p>Emu assigns each read to a reference species, so there are no ASV sequences; the taxon IDs are "
                 "those of the Emu database.</p>")
    o.append(results("2.Taxonomy"))

    # 3. Diversity
    reads_ok = run.counts_are_reads
    o += ["<h2>3. Diversity</h2>", "<h3>Rarefaction</h3>",
          f"<p>How the number of {unit} found grows with the number of reads sampled (vegan "
          f"{e(plots.get('vegan', ''))}); a curve that levels off means the sequencing depth captured the sample's "
          "diversity." + (" Emu's counts are estimated reads and taxa below the minimum abundance were removed, so "
                          "the curve and Chao1 are approximate." if emu and reads_ok else "") + "</p>"]
    if reads_ok:
        o.append(CR._fig(root / "3.Diversity" / "Rarefaction_curves.png", 13))
    else:
        o.append('<p class="note">Not made: the counts of this run are parts per million, not reads.</p>')
    o += ["<h3>Alpha diversity</h3>",
          "<p>Richness (Observed, Chao1) and diversity weighing evenness (Shannon, Simpson) within this sample.</p>"]
    if a:
        o.append(CR._table(["Sample", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson"],
                           [[e(s)] + [("—" if (k in ("Chao1", "se.chao1") and not reads_ok) else
                                       (f"{float(a[k]):.2f}" if str(a.get(k, "")).strip() not in ("", "NA") else "—"))
                                      for k in ("Observed", "Chao1", "se.chao1", "Shannon", "Simpson")]], "num"))
    else:
        o.append("<p>Alpha diversity was not computed for this sample.</p>")
    o.append(results("3.Diversity"))

    # 4. Tables
    o += ["<h2>4. Tables</h2>",
          ("<p>The Excel file 4.Tables/" + e(xlsx.name) + " holds the read tracking"
           + (" and the taxonomic profiles at every rank" if emu else
              ", the taxonomic profiles at every rank and the sequence of each ASV") + ".</p>") if xlsx.exists() else
          '<p class="note">The Excel file was not made on this machine (Python package openpyxl is missing; install '
          "it with: bash setup_amplicon_report.sh).</p>", results("4.Tables")]

    _, ref_text = CR._tax_text(run)
    o.append("<h2>References</h2><ol>" + "".join(f"<li>{e(r)}</li>" for r in CR._refs(run, ref_text)) + "</ol>")
    return ("<!doctype html><html><head><meta charset=\"utf-8\"><title>Amplicon sequencing report</title>"
            f"<style>{CSS}</style></head><body>" + "\n".join(o) + "</body></html>")
