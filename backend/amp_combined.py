"""amp_combined.py — one result for several samples of a run (Download results -> "One combined
result"), laid out like the BI team's own deliverable: the samples side by side in one report and
in every table and figure, instead of a package per sample.

<name>/
  Report_<name>.pdf / .docx   the run report for these samples (the same as Build report)
  1.QC/          read tracking (every sample, every step), reads-per-step figure, the pipeline's read
                 count / quality figures
  2.Cutadapt/    cutadapt's per-sample counts (Illumina runs)
  3.DADA2/ | 3.Emu/ | 3.VSEARCH/   what the denoising / classification step wrote (error models;
                 Emu's per-sample tables; OTU centroids)
  4.ASV/ | 4.OTU/ | 4.Taxa_table/  the feature table: counts raw and with lineage, taxonomy, sequences,
                 features and reads per sample (+ figure), length distribution
  5.Taxa/        csv/ every rank: relative abundance and counts per sample (all taxa, Top 30 wide and
                 long, totals); png/ the report's figures: bars per rank (Top 30), genus heatmap; Krona/ one chart
                 with a ring set per sample; additional_figures/ (R figures the report does not
                 have); NCBI/ or BLAST/ when the run has them
  6.Diversity/   Alpha/ (indices per sample with the group, figure, rarefaction curves + data);
                 Beta/ (Bray-Curtis matrix, PCoA scores + figure, PERMANOVA when there are groups;
                 additional_figures/ with the R beta figures)
  7.Tables/      <name>.xlsx: read tracking, profiles at every rank, alpha (every sample)
  metadata.csv   the samples' groups (when the job has them)
  8.Reads/       quality-filtered reads (when chosen)
  README.txt
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
    m = {"read_tracking.csv": have(D["qc"] / "read_tracking_summary.csv"),
         "asv_table.csv": have(D["feat"] / f"{feat}_count_table_named.csv"),
         "asv_lengths.csv": [],
         "r_tables/rarefaction.csv": have(D["alpha"] / "Rarefaction_curve_data.csv",
                                          D["alpha"] / "Rarefaction_curves.png"),
         "r_tables/taxonomy_heatmap_genus.csv": have(D["taxa"] / "png" / f"Genus_heatmap_top30{rsfx(k)}.png",
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
         "metadata.csv": have(root / "metadata.csv"),
         "PICRUSt2": [rel(D["tab"] / "PICRUSt2")] if (D["tab"] / "PICRUSt2").is_dir() else []}
    for rank in RANK_ORDER:
        m[f"r_tables/taxonomy_{rank.lower()}.csv"] = have(
            D["taxa"] / "csv" / f"{rank}_relative_abundance_per_sample_wide.csv",
            D["taxa"] / "png" / f"{rank}_Top30_relative_abundance_per_sample{rsfx(k)}.png")
    return m


def _readme_body(root: Path, name: str) -> str:
    """The folder list of the docstring, keeping only the folders this result has (a Nanopore run has no
    2.Cutadapt, one sample no Beta, a run without its reads no 8.Reads) and only the 3. / 4. name it uses."""
    out, keep = [], True
    for ln in __doc__.split("\n\n", 1)[1].splitlines():
        m = re.match(r"^  (\d\.[^/ ]+/(?: \| \d\.[^/ ]+/)*)(\s+)(.*)$", ln)
        if m:
            have = [f for f in re.findall(r"(\d\.[^/ ]+)/", m.group(1)) if (root / f).is_dir()]
            keep = bool(have)
            if keep:
                head = have[0] + "/"
                ln = "  " + head + " " * max(1, 15 - len(head)) + m.group(3)
        elif ln.startswith("  metadata.csv"):
            keep = (root / "metadata.csv").exists()
        elif not ln.startswith("    "):
            keep = True
        if keep:
            out.append(ln.replace("<name>", name))
    body = "\n".join(out) + "\n"
    if not (root / "6.Diversity" / "Beta").is_dir():
        body = re.sub(r";\n\s+Beta/ \(.*?figures\)", "", body, flags=re.S)
    return body


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
        names = sorted({k for t in per.values() for k in t},
                       key=lambda k: -sum(per[s].get(k, 0) / tot[s] for s in S))
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


def build_combined(run_dir: Path, upload_dir: Path | None, samples, out_dir: Path, name: str,
                   include_reads: bool = True, ncbi_names: bool = False, tools: dict | None = None,
                   formats=("pdf", "docx"), job_label: str = "", progress=None, read_qc: bool = False) -> dict:
    """read_qc: also NanoPlot (nanopore) or FastQC + fastp (Illumina) of each sample's raw reads, in
    1.QC/<sample>/ — the one-sample results do; the combined one leaves it out (slow on many samples)."""
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
    D = {"qc": root / "1.QC", "cut": root / "2.Cutadapt", "s3": root / f"3.{step3}", "feat": root / f"4.{feat}",
         "taxa": root / "5.Taxa", "alpha": root / "6.Diversity" / "Alpha", "beta": root / "6.Diversity" / "Beta",
         "tab": root / "7.Tables", "reads": root / "8.Reads"}

    rrun = RB.Run(run_dir, RB.KITOME, job_label or name, samples=sel if subset else None)
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
    _csv(D["qc"] / "read_tracking_summary.csv", ["Sample"] + (["Group"] if group else []) + cols + ["perc_retained"],
         [[s] + ([group.get(s, "")] if group else []) + ["" if run.tracking.get(s, {}).get(c) is None
                                                         else run.tracking[s][c] for c in cols]
          + [run.tracking.get(s, {}).get("perc", "")] for s in sel])
    # Files the pipeline wrote for the whole run hold every sample's data: they go in only when this
    # result holds every sample (a one-sample result may be for one customer among several).
    whole = not subset
    if whole:
        _put(run_dir / "read_tracking.csv", D["qc"] / "read_tracking.csv")
    fig(D["qc"] / "Boxplot_read_counts_pipeline_steps", lambda: RB.svg_steps(rrun, an))
    fig(D["qc"] / "Reads_surviving_each_step", lambda: RB.svg_readflow(rrun))
    if whole:
        for f in ("qc_readcount_boxplot.pdf", "read_tracking_plot.pdf"):
            _put(run_dir / f, D["qc"] / f)
        for p in list(run_dir.glob("*QualityProfile*.pdf")) + list(run_dir.glob("Dada2/*QualityProfile*.pdf")) \
                + list(run_dir.glob("r_plots/01_*.pdf")):
            _put(p, D["qc"] / p.name)
        _put(run_dir / "qc" / "read_qc.json", D["qc"] / "read_qc.json")
    if read_qc:
        for s in sel:
            sd = D["qc"] / (CR.safe(s) if len(sel) > 1 else "")
            try:
                if run.kind == "emu" or engine == "vsearch":
                    fs = CR._ont_raw(run, upload_dir, s)
                    if fs:
                        Q.nanoplot(fs, sd / "NanoPlot", CR._tool(tools, "NanoPlot"))
                    else:
                        notes.append(f"{s}: raw reads not found in the upload folder: no NanoPlot")
                elif engine == "pacbio":           # long single reads: NanoPlot, not FastQC
                    r1 = CR._find_reads(upload_dir, s, run.peer_names)[0]
                    if r1:
                        Q.nanoplot([r1], sd / "NanoPlot", CR._tool(tools, "NanoPlot"))
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
    for s in sel:
        _put(run_dir / "cutadapt" / f"{s}.cutadapt.json", D["cut"] / f"{s}.cutadapt.json")
    if (run_dir / "cutadapt").is_dir() and D["cut"].exists():
        rows = []
        for s in sel:
            try:
                rc = json.loads((run_dir / "cutadapt" / f"{s}.cutadapt.json").read_text()).get("read_counts", {})
                rows.append([s, rc.get("input"), rc.get("read1_with_adapter"), rc.get("read2_with_adapter"),
                             rc.get("output")])
            except Exception:
                pass
        _csv(D["cut"] / "cutadapt_summary.csv", ["Sample", "Read_pairs_in", "R1_with_primer", "R2_with_primer",
                                                 "Read_pairs_out"], rows)

    # 3. denoising / classification step
    for p in run_dir.glob("Dada2/*.pdf"):
        _put(p, D["s3"] / p.name)
    if engine == "emu":
        for s in sel:
            for p in (run_dir / "emu_output" / s).glob("*"):
                if p.is_file() and p.suffix in (".tsv", ".csv", ".txt") and not p.name.endswith(".log"):
                    _put(p, D["s3"] / (s if len(sel) > 1 else "") / p.name)
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
    if (run_dir / "metadata.csv").exists():                 # the groups of these samples (the report cites it)
        try:
            with open(run_dir / "metadata.csv", newline="", encoding="utf-8-sig") as fh:
                rd = list(csv.reader(fh))
            rows = [r for r in rd[1:] if r and r[0] in set(sel)]
            if rd and rows:
                with open(root / "metadata.csv", "w", newline="") as fh:
                    csv.writer(fh).writerows([rd[0]] + rows)
        except (OSError, csv.Error):
            pass

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
            fig(D["taxa"] / "png" / f"{rank}_Top30_relative_abundance_per_sample{rsfx(k)}",
                (lambda rank=rank, part=part: RB.svg_stacked(an, rank, 30, samples=part, wide=True,
                                                             groups=gby or None)) if many else
                (lambda rank=rank: RB.svg_stacked(an, rank, 30)))
    for k, part in enumerate(rparts, 1):
        fig(D["taxa"] / "png" / f"Genus_heatmap_top30{rsfx(k)}",
            (lambda part=part: RB.svg_heatmap(an, "Genus", 30, samples=part, wide=True)) if many else
            (lambda: RB.svg_heatmap(an, "Genus", 30)))
    # R figures that draw what the report already shows (bars per rank, group bars, heatmaps) are left
    # out: same data, other colours and order — they read as a second, different result
    if whole:
        for p in list(run_dir.glob("r_plots/22*.pdf")) + [run_dir / "prevalence_abundance.pdf"]:
            _put(p, D["taxa"] / "additional_figures" / p.name)
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
    kt = CR._tool(tools, "ktImportText")
    if kt:
        kd = D["taxa"] / "Krona"
        kd.mkdir(parents=True, exist_ok=True)
        args = []
        for s in sel:
            p = kd / f".{CR.safe(s)}.krona.tsv"
            with open(p, "w") as fh:
                for a in run.sample_asvs(s):
                    fh.write("\t".join([str(a["counts"][s])] + [AD.label(a, r) for r in AD.RANKS]) + "\n")
            args.append(f"{p},{s}")
        try:
            subprocess.run([kt, "-o", str(kd / f"{name}-krona.html")] + args, capture_output=True, timeout=900)
        except (subprocess.TimeoutExpired, OSError) as ex:
            notes.append(f"Krona not made ({ex})")
        for p in kd.glob(".*.krona.tsv"):
            p.unlink()
    else:
        notes.append("KronaTools (ktImportText) not found: no Krona chart")

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
        _put(p, D["alpha"] / "additional_figures" / p.name)
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
            _put(run_dir / f, D["beta"] / "additional_figures" / f)
        _put(run_dir / "UPGMA", D["beta"] / "additional_figures" / "UPGMA")
        for pat in ("05b*", "05c*", "12*", "17*", "20*"):
            for p in run_dir.glob(f"r_plots/{pat}.pdf"):
                _put(p, D["beta"] / "additional_figures" / p.name)

    # 7. Excel
    progress(3, N, "Excel")
    try:
        import amplicon_excel as AX
        D["tab"].mkdir(parents=True, exist_ok=True)
        AX.write_all(D["tab"] / f"{name}.xlsx", run)
    except ImportError:
        notes.append("openpyxl missing: no Excel file")

    _picrust(run_dir, D["tab"] / "PICRUSt2", sel, whole)

    # 8. reads
    if include_reads:
        for s in sel:
            for q in CR._filtered(run, s):
                _put(q, D["reads"] / q.name)

    # report
    progress(4, N, "Writing the report")
    doc = RB.build_profile_html(rrun, title=None, app="NextGen-Amplicon", company="",
                                files=_report_files(root, D, feat, rparts, rsfx),
                                ids={k: a["sv"] for a in run.asvs
                                     for k in {a.get("sequence"), a.get("id"), str(a.get("id") or "").split(";")[0]}
                                     if k})
    stem = root / f"Report_{name}"
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
    for f in root.rglob("*.log"):           # tool logs (NanoPlot's) carry server paths: not for the customer
        f.unlink()
    body = _readme_body(root, name)
    (root / "README.txt").write_text(
        f"{name}: {len(sel)} sample{'s' if len(sel) != 1 else ''} in one result ({', '.join(sel)}).\n"
        + (f"Groups ({an.group_col}): " + ", ".join(f"{g} ({len(v)})" for g, v in an.group_list()) + "\n" if group else "")
        + "\n" + body
        + ("\nNot made in this result:\n" + "".join(f"  - {n}\n" for n in notes) if notes else "")
        + ("\nNCBI names: the OTUs the reference left unassigned carry the name of their NCBI check in the "
           "tables, Excel and Krona; the report lists them in its NCBI section.\n" if ncbi_names else ""),
        encoding="utf-8")
    progress(N, N, "Done")
    return {"dir": str(root), "samples": sel, "notes": notes}
