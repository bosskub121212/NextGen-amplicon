"""amp_combined.py — the result of a run for download (Download results): several samples side by side
in one result ("One combined result") or one sample per result — laid out exactly as the BI team's
deliverable to the customer (U2BI-184, v2.43.0):

<name>/
  Report_<job>_<ONT|PacBio|Illumina>_report.pdf / .docx   the run report for these samples
  1.QC/          each sample's raw-read QC (<sample>/NanoPlot, or <sample>/FastQC + fastp); the
                 pipeline's read count / quality figures
  2.Filter/ | 2.Cutadapt/   read tracking (every sample, every step: csv + figures); per sample:
                 Box_plot.png, read_tracking.tsv, filtered_reads/ (and cutadapt's counts on Illumina)
  3.DADA2/ | 3.Emu/ | 3.VSEARCH/   what the denoising / classification step wrote (error models;
                 Emu's table per sample; OTU centroids)
  4.ASV/ | 4.OTU/ | 4.Taxa_table/  the feature table: counts raw and with lineage, taxonomy, sequences,
                 features and reads per sample, OTU distribution
  5.Taxonomy_Identification/   csv/ (every rank: relative abundance and counts per sample, Top 30
                 wide and long, totals); Taxa_barplot/ (bars per rank, Top 30); Heatmap/ (genus);
                 Krona/ (a chart per sample + all, and a picture of each); species summary;
                 NCBI/ or BLAST/ when the run has them
  6.Diversity/   Alpha/ (indices, figures, rarefaction curves + data, Shannon rarefaction,
                 rank abundance); Beta/ (Bray-Curtis matrix + heatmap, PCoA and PERMANOVA from 3
                 samples, prevalence; figures/ the R beta figures)
  7.Functional_prediction/   PICRUSt2 (16S runs that have it)
Every figure the pipeline wrote as PDF has a PNG beside it, same name (a page each: _p1, _p2 ...).
"""
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import amp_data as AD
import custom_report as CR

RANK_ORDER = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _svg_size(svg: str):
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    return (float(m.group(1)), float(m.group(2))) if m else (800.0, 500.0)


def save_figure(svg: str, stem: Path) -> list[Path]:
    """<stem>.png (cairosvg, twice the drawn size) and <stem>.svg (scales without loss); a PDF in
    place of the PNG when cairosvg is not installed."""
    if not svg:
        return []
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    w, h = _svg_size(svg)
    fixed = re.sub(r'width="100%"( style="max-width:[\d.]+px")?', f'width="{w:g}" height="{h:g}"', svg, count=1)
    if 'xmlns="http://www.w3.org/2000/svg"' not in fixed:
        fixed = fixed.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
    out = [stem.with_suffix(".svg")]
    out[0].write_text(fixed, encoding="utf-8")
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=fixed.encode("utf-8"), write_to=str(stem.with_suffix(".png")),
                         output_width=int(w * 2), background_color="white")
        out.append(stem.with_suffix(".png"))
    except Exception:
        try:
            from weasyprint import HTML
            HTML(string=f"<html><head><style>@page {{ size: {w:g}px {h:g}px; margin: 0 }} body {{ margin: 0 }}"
                        f"</style></head><body>{fixed}</body></html>").write_pdf(str(stem.with_suffix(".pdf")))
            out.append(stem.with_suffix(".pdf"))
        except Exception:
            pass
    return out


def _put(src: Path, dst: Path):
    src, dst = Path(src), Path(dst)
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return True
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    return True


def _csv(path: Path, head, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(rows)


def _picrust(run_dir: Path, dest: Path, sel: list, whole: bool):
    """PICRUSt2's predicted functions (KO, EC, MetaCyc pathways) for these samples — the report shows them.
    A result of every sample takes the files as they are; a part of the run gets each table cut to its
    samples (a column per sample, or a row per sample for the NSTI tables)."""
    import gzip
    src = run_dir / "PICRUSt2"
    if not src.is_dir():
        return
    keep = set(sel)
    for sub in ("KO_metagenome_out", "EC_metagenome_out", "pathways_out"):
        for f in sorted((src / sub).glob("*.tsv.gz")) if (src / sub).is_dir() else []:
            d = dest / sub / f.name
            if whole:
                _put(f, d)
                continue
            try:
                with gzip.open(f, "rt", newline="") as fh:
                    rows = list(csv.reader(fh, delimiter="\t"))
                if not rows:
                    continue
                head = rows[0]
                cols = [i for i, h in enumerate(head) if i == 0 or h in keep or h in ("description", "taxon")]
                if len(cols) > 1 and any(head[i] in keep for i in cols):
                    out = [[r[i] if i < len(r) else "" for i in cols] for r in rows]
                elif head and head[0].lower() == "sample":
                    out = [head] + [r for r in rows[1:] if r and r[0] in keep]
                else:
                    continue                        # no sample in it: not this part's to carry
                d.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(d, "wt", newline="") as fh:
                    csv.writer(fh, delimiter="\t").writerows(out)
            except (OSError, csv.Error, EOFError):
                continue


def _report_files(root: Path, D: dict, feat: str, rparts, rsfx) -> dict:
    """{file of the run folder the report cites: [the files of this result that hold it]} — only files
    that are there."""
    rel = lambda p: str(Path(p).relative_to(root)).replace(os.sep, "/")
    have = lambda *ps: [rel(p) for p in ps if Path(p).exists()]
    k = len(rparts)
    m = {"read_tracking.csv": have(D["filt"] / "read_tracking_summary.csv"),
         "qc/read_qc.json": have(D["filt"] / "read_qc.json"),
         "asv_table.csv": have(D["feat"] / f"{feat}_count_table_named.csv"),
         "asv_lengths.csv": [],
         "r_tables/rarefaction.csv": have(D["alpha"] / "Rarefaction_curve_data.csv",
                                          D["alpha"] / "Rarefaction_curves.png"),
         "r_tables/taxonomy_heatmap_genus.csv": have(D["taxa"] / "Heatmap" / f"Genus_heatmap_top30{rsfx(k)}.png",
                                                     D["taxa"] / "csv" / "Genus_relative_abundance_per_sample_wide.csv"),
         "r_tables/alpha_diversity.csv": have(D["alpha"] / "alpha_diversity.csv", D["alpha"] / "alpha_diversity.png"),
         "alpha_diversity.csv": [],
         "r_tables/pca_scores.csv": have(D["beta"] / "PCoA_Bray_scores.csv", D["beta"] / "PCoA_Bray.png"),
         "beta_stats.txt": have(*sorted(D["beta"].glob("PERMANOVA_*.csv"))) if D["beta"].is_dir() else [],
         "r_tables/beta_heatmap.csv": have(D["beta"] / "bray_curtis_distance_matrix.csv",
                                           D["beta"] / "Bray_Curtis_heatmap.png"),
         "asvs.fasta": have(D["feat"] / f"{feat}_sequences.fasta"),
         "otus.fasta": have(D["s3"] / "otus.fasta"),
         "blast/taxonomy_blast.csv": have(D["taxa"] / "BLAST" / "taxonomy_blast.csv"),
         "taxonomy_table_dada2.csv": have(D["feat"] / "taxonomy_table_dada2.csv"),
         "metadata.csv": [],
         "PICRUSt2": [rel(D["func"] / "PICRUSt2")] if (D["func"] / "PICRUSt2").is_dir() else []}
    for rank in RANK_ORDER:
        m[f"r_tables/taxonomy_{rank.lower()}.csv"] = have(
            D["taxa"] / "csv" / f"{rank}_relative_abundance_per_sample_wide.csv",
            D["taxa"] / "Taxa_barplot" / f"{rank}_Top30_relative_abundance_per_sample{rsfx(k)}.png")
    return m


def _taxa_tables(run: AD.RunData, d: Path, group: dict):
    """Per rank: every taxon's relative abundance (%) and reads per sample, the Top 30 wide and long,
    and each taxon's total share — the BI team's Top_taxa_all_levels/csv."""
    S = list(run.samples)
    # share of every read of the sample, Emu's unassigned included (it is a row of the tables), as the report
    # and its charts count it: a rank sums to 100 (dividing by the assigned reads gave 100.4% with that row)
    tot = {s: sum(a["counts"].get(s, 0) for a in run.asvs) or 1 for s in S}
    allt = sum(tot.values()) or 1
    for rank in RANK_ORDER:
        per = {s: run.rank_table(rank, s) for s in S}
        names = sorted({k for t in per.values() for k in t},     # ties by name: the same order every build
                       key=lambda k: (-round(sum(per[s].get(k, 0) / tot[s] for s in S), 12), k))
        pct = lambda k, s: 100 * per[s].get(k, 0) / tot[s]
        _csv(d / f"{rank}_relative_abundance_per_sample_wide.csv", [rank] + S,
             [[k] + [f"{pct(k, s):.4f}" for s in S] for k in names])
        _csv(d / f"{rank}_read_counts_per_sample_wide.csv", [rank] + S,
             [[k] + [str(per[s].get(k, 0)) for s in S] for k in names])
        top = [k for k in names if not k.startswith("Unassigned")][:30]
        _csv(d / f"{rank}_Top30_relative_abundance_per_sample_wide.csv", [rank] + S,
             [[k] + [f"{pct(k, s):.4f}" for s in S] for k in top]
             + [["Other"] + [f"{100 - sum(pct(k, s) for k in top):.4f}" for s in S]])
        _csv(d / f"{rank}_Top30_relative_abundance_per_sample_long.csv",
             ["Sample"] + (["Group"] if group else []) + [rank, "Relative_abundance_pct", "Reads"],
             [[s] + ([group.get(s, "")] if group else []) + [k, f"{pct(k, s):.4f}", per[s].get(k, 0)]
              for s in S for k in top])
        _csv(d / f"{rank}_all_taxa_total_relative_abundance.csv", [rank, "Reads", "Relative_abundance_pct",
                                                                    "Samples_present"],
             [[k, sum(per[s].get(k, 0) for s in S), f"{100 * sum(per[s].get(k, 0) for s in S) / allt:.4f}",
               sum(1 for s in S if per[s].get(k, 0) > 0)] for k in names])


PNG_DPI = 200


def _pymupdf():
    try:
        import pymupdf as fitz           # PyMuPDF >= 1.24 (the name "fitz" is being retired)
    except ImportError:
        import fitz
    return fitz


def _pdf_tool():
    """What turns a PDF page into a PNG here: PyMuPDF, else poppler's pdftoppm, else Ghostscript."""
    try:
        _pymupdf()
        return "fitz"
    except Exception:
        pass
    for t in ("pdftoppm", "gs"):
        if shutil.which(t):
            return t
    return None


def pdf_to_png(pdf: Path, tool: str, dpi: int = PNG_DPI, max_pages: int = 50) -> list[Path]:
    """<stem>.png beside a one-page PDF, <stem>_p1.png, _p2.png ... for more pages; white background."""
    import tempfile
    pdf = Path(pdf)
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        if tool == "fitz":
            fitz = _pymupdf()
            with fitz.open(str(pdf)) as doc:
                for i, page in enumerate(doc):
                    if i >= max_pages:
                        break
                    page.get_pixmap(dpi=dpi, alpha=False).save(str(t / f"p-{i + 1:03d}.png"))
        elif tool == "pdftoppm":
            subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-l", str(max_pages), str(pdf), str(t / "p")],
                           capture_output=True, timeout=600, check=True)
        elif tool == "gs":
            subprocess.run(["gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", f"-r{dpi}",
                            "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-dLastPage={max_pages}",
                            f"-sOutputFile={t / 'p-%03d.png'}", str(pdf)],
                           capture_output=True, timeout=600, check=True)
        else:
            return []
        pages = sorted(t.glob("p-*.png"), key=lambda x: int(re.sub(r"\D", "", x.stem) or 0))
        if not pages:
            raise RuntimeError("no page rendered")
        out = [pdf.with_suffix(".png")] if len(pages) == 1 else \
              [pdf.with_name(f"{pdf.stem}_p{i + 1}.png") for i in range(len(pages))]
        for src, dst in zip(pages, out):
            shutil.move(str(src), str(dst))
    return out


def pngs_for_pdfs(root: Path, notes: list) -> int:
    """A PNG beside every PDF figure of the result (the R pipeline's figures, NanoPlot's, a figure saved
    as PDF when cairosvg is missing), for the BI team to paste into their reports. The report itself
    (the PDF at the top) and the reads are left alone; a PDF that already has its PNG is skipped."""
    pdfs = [p for p in sorted(Path(root).rglob("*.pdf"))
            if p.parent != Path(root) and "filtered_reads" not in p.parts
            and not p.with_suffix(".png").exists() and not p.with_name(p.stem + "_p1.png").exists()]
    if not pdfs:
        return 0
    tool = _pdf_tool()
    if not tool:
        notes.append(f"PNG copies of the {len(pdfs)} PDF figures not made: no PDF renderer on the server "
                     "(pip install pymupdf, or sudo apt-get install poppler-utils)")
        return 0
    n, bad = 0, []
    for p in pdfs:
        try:
            pdf_to_png(p, tool)
            n += 1
        except Exception as ex:
            bad.append(f"{p.relative_to(root)} ({type(ex).__name__})")
    if bad:
        notes.append("PNG not made for " + ", ".join(bad[:10]) + (" ..." if len(bad) > 10 else ""))
    return n


SUN_PAL = ["#E08FBE", "#D4A373", "#44B9C6", "#A8A6E3", "#9FB06B", "#C69ADF", "#6FB3D9", "#E39C8E",
           "#5FC4A8", "#C4A96A", "#8EA9DE", "#E3B95F", "#7FC77A", "#D98CA4", "#62B6A8", "#B9A0D8"]


def _arc(cx, cy, r0, r1, a0, a1):
    """Ring sector from angle a0 to a1 (radians, clockwise from 12 o'clock)."""
    import math
    if a1 - a0 >= 2 * math.pi - 1e-6:                 # a full ring: two halves (one arc cannot close)
        m = a0 + math.pi
        return _arc(cx, cy, r0, r1, a0, m) + " " + _arc(cx, cy, r0, r1, m, a1)
    pt = lambda r, a: (cx + r * math.sin(a), cy - r * math.cos(a))
    big = 1 if a1 - a0 > math.pi else 0
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = pt(r1, a0), pt(r1, a1), pt(r0, a1), pt(r0, a0)
    return (f"M{x1:.2f},{y1:.2f} A{r1:.2f},{r1:.2f} 0 {big} 1 {x2:.2f},{y2:.2f} L{x3:.2f},{y3:.2f} "
            f"A{r0:.2f},{r0:.2f} 0 {big} 0 {x4:.2f},{y4:.2f} Z")


def svg_sunburst(run: AD.RunData, s: str) -> str:
    """The sample's Krona as a still picture: rings Phylum (inside) ... Species (outside), each taxon's
    angle its share of the sample's reads; a lineage stops at its first unnamed rank."""
    import math
    from html import escape as esc
    tree = {"v": 0, "c": {}}
    for a in run.sample_asvs(s):
        n = a["counts"].get(s, 0)
        if n <= 0:
            continue
        path = ["Unassigned"] if a.get("not_taxon") else []
        if not path:
            for r in RANK_ORDER:
                v = str(a["tax"].get(r) or "").strip()
                if not AD._v(v):
                    break
                path.append(v)
        path = path or ["Unassigned"]
        tree["v"] += n
        node = tree
        for name in path:
            node = node["c"].setdefault(name, {"v": 0, "c": {}})
            node["v"] += n
    if not tree["v"]:
        return ""
    W, H, cx, cy, hole, ring = 1000, 1040, 500, 530, 58, 70
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
           f'font-family="DejaVu Sans, Arial, sans-serif"><rect width="{W}" height="{H}" fill="#fff"/>',
           f'<text x="{W/2}" y="38" font-size="26" text-anchor="middle" fill="#222">{esc(s)}</text>']
    labels = []

    def draw(node, depth, a0, a1, colour_at):
        kids = sorted(node["c"].items(), key=lambda kv: -kv[1]["v"])
        a = a0
        for i, (name, ch) in enumerate(kids):
            span = (a1 - a0) * ch["v"] / node["v"]
            col = "#BDBDBD" if name == "Unassigned" else SUN_PAL[(colour_at + i) % len(SUN_PAL)]
            r0, r1 = hole + depth * ring, hole + (depth + 1) * ring
            out.append(f'<path d="{_arc(cx, cy, r0, r1, a, a + span)}" fill="{col}" stroke="#fff" '
                       f'stroke-width="2"><title>{esc(name)}: {ch["v"]:,} reads '
                       f'({100 * ch["v"] / tree["v"]:.2f}%)</title></path>')
            rm = (r0 + r1) / 2
            m = a + span / 2
            arc_len = span * rm
            # as Krona does: the name along the ring where the sector is wide, across it (outward) where it
            # is narrow; never wider than its own sector, so no two labels overlap; cut with "…" to fit
            along = arc_len >= ring
            chord = 2 * math.sqrt(max((r1 - 4) ** 2 - rm ** 2, 0))  # a straight line along a curved ring
            room, high = (min(arc_len - 10, chord), ring) if along else (ring - 10, arc_len)
            fit = min(int(room / 6.4), 30)
            if high >= 13 and fit >= 5:
                txt = name if len(name) <= fit else name[:fit - 1] + "…"
                deg = math.degrees(m)
                rot = (deg if math.cos(m) >= 0 else deg + 180) if along else \
                      (deg - 90 if math.sin(m) >= 0 else deg + 90)
                x, y = cx + rm * math.sin(m), cy - rm * math.cos(m)
                if depth == 0 and span >= 2 * math.pi - 1e-6:          # one phylum: the name in the middle
                    rot, y = 0, cy - rm
                labels.append(f'<text x="{x:.1f}" y="{y:.1f}" dy="4" transform="rotate({rot:.1f} {x:.1f} {y:.1f})" '
                              f'font-size="11" text-anchor="middle" fill="#111">{esc(txt)}</text>')
            if ch["c"] and depth + 1 < len(RANK_ORDER):
                draw(ch, depth + 1, a, a + span, colour_at + i * 3 + 1)
            a += span

    draw(tree, 0, 0.0, 2 * math.pi, 0)
    out += labels
    out.append(f'<text x="{W/2}" y="{H - 14}" font-size="12" text-anchor="middle" fill="#666">'
               f'{tree["v"]:,} reads · rings from the centre: {", ".join(RANK_ORDER)}</text></svg>')
    return "".join(out)


def _png_only(svg: str, stem: Path, notes: list):
    """<stem>.png alone (the BI layout has no SVG beside these); the SVG stays when no PNG can be made."""
    got = save_figure(svg, stem) if svg else []
    if Path(stem).with_suffix(".png").exists():
        Path(stem).with_suffix(".svg").unlink(missing_ok=True)
    elif svg:
        notes.append(f"{Path(stem).name}.png not made (cairosvg missing): SVG instead")
    return got


class _OneTrack:
    """report_builder.Run seen with one sample's read tracking (Box_plot.png of that sample)."""
    def __init__(self, run, row):
        self._r, self.tracking = run, [row]

    def __getattr__(self, k):
        return getattr(self._r, k)


def _strip_num(name: str) -> str:
    """'03b_shannon_rarefaction.pdf' -> 'shannon_rarefaction.pdf' (the BI team's names)."""
    return re.sub(r"^\d+[a-z]?_", "", name)


def build_combined(run_dir: Path, upload_dir: Path | None, samples, out_dir: Path, name: str,
                   include_reads: bool = True, ncbi_names: bool = False, tools: dict | None = None,
                   formats=("pdf", "docx"), job_label: str = "", progress=None, read_qc: bool = True) -> dict:
    """read_qc: also NanoPlot (nanopore) or FastQC + fastp (Illumina) of each sample's raw reads, in
    1.QC/<sample>/ (1.QC/ for one sample). Every result has it since v2.42.0 (the combined one left it
    out before: the BI team found no NanoPlot in 1.QC)."""
    import report_builder as RB
    import qc_tools as Q
    tools = dict(tools or {})
    progress = progress or (lambda i, n, lab: None)
    run_dir = Path(run_dir)
    full = AD.load_run(run_dir)
    if ncbi_names:
        full.use_ncbi_names()
    sel = [s for s in full.samples if s in set(samples)] if samples else list(full.samples)
    subset = len(sel) < len(full.samples)
    run = full.subset(sel) if subset else full
    root = Path(out_dir) / name
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    notes, N = [], 6
    engine = getattr(run, "engine", run.kind)
    feat = {"emu": "Taxa_table", "vsearch": "OTU"}.get(engine, "ASV")
    step3 = {"emu": "Emu", "vsearch": "VSEARCH"}.get(engine, "DADA2")
    rrun = RB.Run(run_dir, RB.KITOME, job_label or name, samples=sel if subset else None)
    long_read = run.kind == "emu" or engine == "pacbio" or (engine == "vsearch" and not getattr(rrun, "paired", False))
    tech = "ONT" if run.kind == "emu" or (engine == "vsearch" and long_read) else "PacBio" if engine == "pacbio" \
        else "Illumina"
    D = {"qc": root / "1.QC", "filt": root / ("2.Filter" if long_read else "2.Cutadapt"), "s3": root / f"3.{step3}",
         "feat": root / f"4.{feat}", "taxa": root / "5.Taxonomy_Identification",
         "alpha": root / "6.Diversity" / "Alpha", "beta": root / "6.Diversity" / "Beta",
         "func": root / "7.Functional_prediction"}
    sdir = lambda base, s: base / CR.safe(s)             # a folder per sample, as the BI team lays them out
    an = RB.Analysis(rrun)
    S = an.samples
    group = {s: an.groups[s] for s in S if s in an.groups}
    many = len(S) > RB.MANY
    if group:                                   # bars of one group next to each other
        S_ord = [x for _g, v in an.group_list() for x in v if x in S]
        S_ord += [x for x in S if x not in S_ord]
    else:
        S_ord = list(S)
    parts = RB._chunks(S_ord) if many else [S_ord]
    sfx = lambda k: f"_part{k}" if len(parts) > 1 else ""

    def fig(stem, f):
        try:
            save_figure(f(), stem)
        except Exception as ex:
            notes.append(f"figure {Path(stem).name} not made ({ex})")

    # 1. QC
    progress(0, N, "QC and read tracking")
    cols = [c for c in run.track_cols if c != "perc"]
    _csv(D["filt"] / "read_tracking_summary.csv", ["Sample"] + (["Group"] if group else []) + cols + ["perc_retained"],
         [[s] + ([group.get(s, "")] if group else []) + ["" if run.tracking.get(s, {}).get(c) is None
                                                         else run.tracking[s][c] for c in cols]
          + [run.tracking.get(s, {}).get("perc", "")] for s in sel])
    # Files the pipeline wrote for the whole run hold every sample's data: they go in only when this
    # result holds every sample (a one-sample result may be for one customer among several).
    whole = not subset
    if whole:
        _put(run_dir / "read_tracking.csv", D["filt"] / "read_tracking.csv")
    fig(D["filt"] / "Boxplot_read_counts_pipeline_steps", lambda: RB.svg_steps(rrun, an))
    fig(D["filt"] / "Reads_surviving_each_step", lambda: RB.svg_readflow(rrun))
    rows_by = {str(next(iter(t.values()))): t for t in (rrun.tracking or [])}
    tcols = [c for c in run.track_cols if c != "perc"]
    for s in sel:
        t = run.tracking.get(s, {})
        sd = sdir(D["filt"], s)
        sd.mkdir(parents=True, exist_ok=True)
        with open(sd / "read_tracking.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["Sample"] + tcols + ["perc"])
            w.writerow([s] + ["" if t.get(c) is None else t.get(c) for c in tcols] + [t.get("perc", "")])
        if s in rows_by:
            try:
                _png_only(RB.svg_steps(_OneTrack(rrun, rows_by[s]), an), sd / "Box_plot", notes)
            except Exception as ex:
                notes.append(f"{s}: Box_plot.png not made ({ex})")
        _put(run_dir / "cutadapt" / f"{s}.cutadapt.json", sd / f"{s}.cutadapt.json")   # Illumina runs
        if include_reads:
            for q in CR._filtered(run, s):
                _put(q, sd / "filtered_reads" / q.name)
    if whole:
        for f in ("qc_readcount_boxplot.pdf", "read_tracking_plot.pdf"):
            _put(run_dir / f, D["qc"] / f)
        for p in list(run_dir.glob("*QualityProfile*.pdf")) + list(run_dir.glob("Dada2/*QualityProfile*.pdf")) \
                + list(run_dir.glob("r_plots/01_*.pdf")):
            _put(p, D["qc"] / p.name)
        _put(run_dir / "qc" / "read_qc.json", D["filt"] / "read_qc.json")
    def _nano(fs, sd, s):
        tool = CR._tool(tools, "NanoPlot")
        if not tool:
            notes.append(f"{s}: NanoPlot is not installed on the server (setup_amplicon_report.sh): no NanoPlot")
        elif not Q.nanoplot(fs, sd / "NanoPlot", tool, threads=min(4, os.cpu_count() or 2)):
            notes.append(f"{s}: NanoPlot did not finish: no NanoPlot")

    if read_qc:
        for k, s in enumerate(sel):
            if len(sel) > 1:
                progress(0, N, f"Read QC of the raw reads: {s} ({k + 1}/{len(sel)})")
            sd = sdir(D["qc"], s)
            try:
                if run.kind == "emu" or engine == "vsearch":
                    fs = CR._ont_raw(run, upload_dir, s)
                    if fs:
                        _nano(fs, sd, s)
                    else:
                        notes.append(f"{s}: raw reads not found in the upload folder: no NanoPlot")
                elif engine == "pacbio":           # long single reads: NanoPlot, not FastQC
                    r1 = CR._find_reads(upload_dir, s, run.peer_names)[0]
                    if r1:
                        _nano([r1], sd, s)
                    else:
                        notes.append(f"{s}: raw reads not found in the upload folder: no NanoPlot")
                else:
                    raw = CR._find_reads(upload_dir, s, run.peer_names)
                    r1, r2 = raw if raw else (None, None)
                    if r1:
                        Q.fastqc([x for x in (r1, r2) if x], sd / "FastQC", CR._tool(tools, "fastqc"))
                        Q.fastp_check(r1, r2, sd / "fastp", CR._tool(tools, "fastp"))
                    else:
                        notes.append(f"{s}: raw reads not found in the upload folder: no FastQC / fastp")
            except Exception as ex:
                notes.append(f"{s}: read QC not made ({ex})")

    # 2. cutadapt
    if (run_dir / "cutadapt").is_dir() and any((run_dir / "cutadapt" / f"{s}.cutadapt.json").exists() for s in sel):
        rows = []
        for s in sel:
            try:
                rc = json.loads((run_dir / "cutadapt" / f"{s}.cutadapt.json").read_text()).get("read_counts", {})
                rows.append([s, rc.get("input"), rc.get("read1_with_adapter"), rc.get("read2_with_adapter"),
                             rc.get("output")])
            except Exception:
                pass
        _csv(D["filt"] / "cutadapt_summary.csv", ["Sample", "Read_pairs_in", "R1_with_primer", "R2_with_primer",
                                                 "Read_pairs_out"], rows)

    # 3. denoising / classification step
    for p in run_dir.glob("Dada2/*.pdf"):
        _put(p, D["s3"] / p.name)
    if engine == "emu":
        for s in sel:
            for p in (run_dir / "emu_output" / s).glob("*"):
                if p.is_file() and p.suffix in (".tsv", ".csv", ".txt") and not p.name.endswith(".log"):
                    _put(p, sdir(D["s3"], s) / p.name)
    if engine == "vsearch" and (run_dir / "otus.fasta").exists():
        ids = {str(a.get("id", "")).split(";")[0] for a in run.asvs}
        keep, cur = [], False
        for line in (run_dir / "otus.fasta").read_text(errors="replace").splitlines():
            if line.startswith(">"):
                cur = line[1:].split(";")[0].split()[0] in ids
            if cur:
                keep.append(line)
        if keep:
            D["s3"].mkdir(parents=True, exist_ok=True)
            (D["s3"] / "otus.fasta").write_text("\n".join(keep) + "\n")

    # 4. feature table
    progress(1, N, "Feature and taxa tables")
    tot = {s: sum(a["counts"].get(s, 0) for a in run.asvs) for s in sel}
    emu = run.kind == "emu"     # Emu: its tax_id beside the ID, so the unassigned rows are named, not blank
    # Emu's tax_id / VSEARCH's OTU name beside the ID: the pipeline's own name of each row
    idcol = "Emu_tax_id" if engine == "emu" else "OTU_id" if engine == "vsearch" else None
    _csv(D["feat"] / f"{feat}_count_table_named.csv", ["ID"] + ([idcol] if idcol else []) + AD.RANKS + sel
         + ["Total"] + (["Sequence"] if not emu else []),
         [[a["sv"]] + ([a.get("id", "")] if idcol else [])
          + [("Unassigned" if a.get("not_taxon") and r == "Phylum" else a["tax"].get(r, "")) for r in AD.RANKS]
          + [a["counts"].get(s, 0) for s in sel] + [a["total"]]
          + ([a.get("sequence", "")] if not emu else []) for a in run.asvs])
    _csv(D["feat"] / f"{feat}_counts_per_sample.csv", ["Sample"] + (["Group"] if group else [])
         + ["Reads", f"{'Taxa' if engine == 'emu' else feat + 's'}"],
         [[s] + ([group.get(s, "")] if group else []) + [tot[s], sum(1 for a in run.asvs if a["counts"].get(s, 0) > 0)]
          for s in sel])
    if whole:
        for f in ("asv_table.csv", "taxonomy_table.csv", "taxonomy.csv", "asvs.fasta", "asv_lengths.csv",
                  "asv_length_distribution.pdf", "observed_asvs.pdf", "phylo_tree.nwk"):
            _put(run_dir / f, D["feat"] / f)
    # the sequences under the IDs of the count table (SV001 ...): the pipeline's asvs.fasta, when the result
    # carries it, numbers them its own way. Old COX1 / PacBio tables carry ASV ids, no sequences: no FASTA then
    if run.kind != "emu" and any(a.get("sequence") for a in run.asvs):
        with open(D["feat"] / f"{feat}_sequences.fasta", "w") as fh:
            for a in run.asvs:
                if a.get("sequence"):
                    fh.write(f">{a['sv']} {AD.label(a, 'Species')}\n{a['sequence']}\n")
    for p in (run_dir.glob("r_plots/18_*.pdf") if whole else []):
        _put(p, D["feat"] / p.name)
    if (run_dir / "taxonomy_table_dada2.csv").exists():    # BLAST runs: the DADA2 naming kept beside, these ASVs
        with open(run_dir / "taxonomy_table_dada2.csv", newline="", encoding="utf-8-sig") as fh:
            rd = list(csv.reader(fh))
        sq = {a.get("sequence") for a in run.asvs if a.get("sequence")}
        rows = [r for r in rd[1:] if r and r[0] in sq]
        if rd and rows:
            D["feat"].mkdir(parents=True, exist_ok=True)
            with open(D["feat"] / "taxonomy_table_dada2.csv", "w", newline="") as fh:
                csv.writer(fh).writerows([rd[0]] + rows)
    # 5. taxa
    _taxa_tables(run, D["taxa"] / "csv", group)
    # the report's own figures, drawn by the same calls as in report_builder.build_profile_html, so a PNG
    # here is the figure of the report (the BI team found the bars here differed from the report's)
    gby = an.groups if len({an.groups.get(x) for x in S} - {None}) >= 2 else {}
    rparts = (RB._group_chunks(S, gby) if gby else RB._chunks(S)) if many else [S]
    rsfx = lambda k: f"_part{k}" if len(rparts) > 1 else ""
    for rank in RANK_ORDER:
        if rank not in rrun.rank_cols:
            continue
        for k, part in enumerate(rparts, 1):
            fig(D["taxa"] / "Taxa_barplot" / f"{rank}_Top30_relative_abundance_per_sample{rsfx(k)}",
                (lambda rank=rank, part=part: RB.svg_stacked(an, rank, 30, samples=part, wide=True,
                                                             groups=gby or None)) if many else
                (lambda rank=rank: RB.svg_stacked(an, rank, 30)))
    for k, part in enumerate(rparts, 1):
        fig(D["taxa"] / "Heatmap" / f"Genus_heatmap_top30{rsfx(k)}",
            (lambda part=part: RB.svg_heatmap(an, "Genus", 30, samples=part, wide=True)) if many else
            (lambda: RB.svg_heatmap(an, "Genus", 30)))
    # R figures that draw what the report already shows (bars per rank, group bars, heatmaps) are left
    # out: same data, other colours and order — they read as a second, different result
    if whole:
        # the BI team's places: prevalence with the beta figures, the species summary in 5., no number prefix
        for p in sorted(run_dir.glob("r_plots/22*.pdf")):
            nm = _strip_num(p.name)
            if nm.startswith("prevalence"):
                _put(p, D["beta"] / p.name)
            else:
                _put(p, D["taxa"] / nm)
        if not any(run_dir.glob("r_plots/22_*.pdf")):
            _put(run_dir / "prevalence_abundance.pdf", D["beta"] / "prevalence_abundance.pdf")
    ids = {str(a.get("id", "")).split(";")[0] for a in run.asvs}
    seqs = {a.get("sequence") for a in run.asvs if a.get("sequence")}
    if (run_dir / "ncbi_unassigned.tsv").exists():           # the OTUs of these samples only
        with open(run_dir / "ncbi_unassigned.tsv", newline="") as fh:
            rd = list(csv.DictReader(fh, delimiter="\t"))
        rows = [r for r in rd if str(r.get("otu", "")).split(";")[0] in ids]
        if rows:
            (D["taxa"] / "NCBI").mkdir(parents=True, exist_ok=True)
            with open(D["taxa"] / "NCBI" / "unassigned_otus_ncbi.tsv", "w", newline="") as fh:
                w = csv.DictWriter(fh, list(rd[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
    if (run_dir / "blast" / "taxonomy_blast.csv").exists():  # the ASVs of these samples only
        with open(run_dir / "blast" / "taxonomy_blast.csv", newline="") as fh:
            rd = list(csv.DictReader(fh))
        rows = [r for r in rd if r.get("sequence") in seqs]
        if rows:
            (D["taxa"] / "BLAST").mkdir(parents=True, exist_ok=True)
            with open(D["taxa"] / "BLAST" / "taxonomy_blast.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, list(rd[0].keys()))
                w.writeheader()
                w.writerows(rows)
        _put(run_dir / "blast" / "blast_info.json", D["taxa"] / "BLAST" / "blast_info.json")
    kd = D["taxa"] / "Krona"
    kt = CR._tool(tools, "ktImportText")
    if kt:
        kd.mkdir(parents=True, exist_ok=True)
        args = []
        for s in sel:
            p = kd / f".{CR.safe(s)}.krona.tsv"
            with open(p, "w") as fh:
                for a in run.sample_asvs(s):
                    fh.write("\t".join([str(a["counts"][s])] + [AD.label(a, r) for r in AD.RANKS]) + "\n")
            args.append((s, f"{p},{s}"))
        jobs = [(kd / f"{CR.safe(s)}-krona.html", [arg]) for s, arg in args]
        if len(args) > 1:
            jobs.append((kd / "all-krona.html", [arg for _s, arg in args]))
        for out_html, a_ in jobs:
            try:
                subprocess.run([kt, "-o", str(out_html)] + a_, capture_output=True, timeout=900)
            except (subprocess.TimeoutExpired, OSError) as ex:
                notes.append(f"Krona {out_html.name} not made ({ex})")
        for p in kd.glob(".*.krona.tsv"):
            p.unlink()
    else:
        notes.append("KronaTools (ktImportText) not found: no Krona chart")
    for s in sel:                                  # a picture of each Krona, to paste into a report
        try:
            _png_only(svg_sunburst(run, s), kd / f"{CR.safe(s)}_krona_static", notes)
        except Exception as ex:
            notes.append(f"{s}: Krona picture not made ({ex})")

    # 6. diversity
    progress(2, N, "Diversity")
    al = an.alpha()
    keys = ["Observed", "Chao1", "ACE", "Shannon", "Simpson"]
    _csv(D["alpha"] / "alpha_diversity.csv", ["Sample"] + (["Group"] if group else []) + ["Reads"] + keys,
         [[s] + ([group.get(s, "")] if group else []) + [round(an.comm_tot.get(s, 0))]
          + [(round(al.get(s, {}).get(k), 4) if isinstance(al.get(s, {}).get(k), float) else al.get(s, {}).get(k, "")) for k in keys] for s in S])
    fig(D["alpha"] / "alpha_diversity", lambda: RB.svg_alpha(an, al))
    if group and len(set(group.values())) >= 2:
        try:
            gl = [(g, [x for x in v if x in S]) for g, v in an.group_list()]
            fig(D["alpha"] / "alpha_diversity_by_group", lambda: RB.svg_alpha_groups(gl, al))
        except Exception as ex:
            notes.append(f"alpha by group not made ({ex})")
    fig(D["alpha"] / "Rarefaction_curves", lambda: RB.svg_rare(an))
    try:
        rf = an.rarefaction()
        _csv(D["alpha"] / "Rarefaction_curve_data.csv", ["Sample", "Depth", "Expected_features"],
             [[s, d_, f"{v:.3f}"] for s, pts in rf.items() for d_, v in pts])
    except Exception as ex:
        notes.append(f"rarefaction data not written ({ex})")
    # R's alpha and rarefaction figures repeat the report's; its Shannon rarefaction and rank-abundance do not
    for p in (list(run_dir.glob("r_plots/03b_*.pdf")) + list(run_dir.glob("r_plots/03c_*.pdf")) if whole else []):
        _put(p, D["alpha"] / _strip_num(p.name))
    try:
        if len(S) < 2:
            raise RuntimeError("one sample: no beta diversity (it compares samples)")
        Dm = an.bray()
        _csv(D["beta"] / "bray_curtis_distance_matrix.csv", [""] + S,
             [[s] + [f"{x:.4f}" for x in row] for s, row in zip(S, Dm)])
        if len(S) >= 3:
            coords, pct = RB._pcoa(Dm)
            _csv(D["beta"] / "PCoA_Bray_scores.csv", ["Sample"] + (["Group"] if group else [])
                 + [f"PCoA1 ({pct[0]:.1f}%)", f"PCoA2 ({pct[1]:.1f}%)"],
                 [[s] + ([group.get(s, "")] if group else []) + [f"{c[0]:.5f}", f"{c[1]:.5f}"]
                  for s, c in zip(S, coords)])
            fig(D["beta"] / "PCoA_Bray", lambda: RB.svg_pcoa(an, wide=many))
        else:
            notes.append("PCoA needs 3 samples or more: not drawn")
        if len(S) >= 3:
            fig(D["beta"] / "UPGMA_bray", lambda: RB.svg_upgma(an))
        if len(S) <= 40:                                # the report's Bray-Curtis heatmap (up to 40 samples)
            fig(D["beta"] / "Bray_Curtis_heatmap", lambda: RB.svg_bc_heat(an))
        pm = an.permanova_groups()
        if pm:
            _csv(D["beta"] / f"PERMANOVA_{an.group_col or 'group'}.csv",
                 ["Grouping", "Groups", "Samples", "R2", "F", "p", "Permutations"],
                 [[an.group_col or "group", pm.get("k"), pm.get("n"), pm.get("R2"), pm.get("F"), pm.get("p"),
                   pm.get("perms")]])
    except Exception as ex:
        notes.append(f"beta diversity tables not made ({ex})")
    if whole:
        for f in ("beta_stats.txt", "beta_nmds_jaccard.pdf", "pca_scores.csv", "nmds_bray.csv"):
            _put(run_dir / f, D["beta"] / "figures" / f)
        _put(run_dir / "UPGMA", D["beta"] / "figures" / "UPGMA")
        for pat in ("05b*", "05c*", "12*", "17*", "20*"):
            for p in run_dir.glob(f"r_plots/{pat}.pdf"):
                _put(p, D["beta"] / "figures" / p.name)

    # 7. PICRUSt2 (the BI layout has no Excel and no 8.Reads: the filtered reads are in 2.Filter/<sample>/)
    progress(3, N, "Functional prediction")
    _picrust(run_dir, D["func"] / "PICRUSt2", sel, whole)

    # report
    progress(4, N, "Writing the report")
    doc = RB.build_profile_html(rrun, title=None, app="NextGen-Amplicon", company="",
                                files=_report_files(root, D, feat, rparts, rsfx),
                                ids={k: a["sv"] for a in run.asvs
                                     for k in {a.get("sequence"), a.get("id"), str(a.get("id") or "").split(";")[0]}
                                     if k})
    stem = root / f"Report_{CR.safe(job_label or name)}_{tech}_report"
    root.mkdir(parents=True, exist_ok=True)
    if "html" in formats:
        stem.with_suffix(".html").write_text(doc, encoding="utf-8")
    if "pdf" in formats:
        try:
            from weasyprint import HTML
            HTML(string=doc).write_pdf(str(stem.with_suffix(".pdf")))
        except Exception as ex:
            notes.append(f"PDF not made ({ex})")
    if "docx" in formats:
        try:
            import docx_export as DX
            DX.html_to_docx(doc, stem.with_suffix(".docx"))
        except Exception as ex:
            notes.append(f"Word file not made ({ex})")
    progress(5, N, "PNG copies of the PDF figures")
    pngs_for_pdfs(root, notes)
    for f in root.rglob("*.log"):           # tool logs (NanoPlot's) carry server paths: not for the customer
        f.unlink()
    progress(N, N, "Done")
    return {"dir": str(root), "samples": sel, "notes": notes, "report": stem.name}
