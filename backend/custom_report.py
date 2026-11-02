#!/usr/bin/env python3
"""custom_report.py — the per-sample customer report for 12S / CytB eDNA runs.

Same headings, order and figure placement as the lab's "Custom Metagenomic Sequencing
Report" (one block per sample), with every sentence written from what this run really
did: its parameters, tool versions, reference and checks. The cover's customer fields are
left blank for the BI team; there is no company name or page art.

    build(run_dir, upload_dir, out_dir) -> {"html", "pdf", "docx", "zip", "samples"}

out_dir gets SampleN_<name>/Images_and_Outfiles/... (every file a caption names),
<job>_all_samples.xlsx, the report (HTML, PDF, Word) and a zip of it all.
Figures: r_scripts/custom_report_plots.R (phyloseq, ggplot2, microbiomeutilities or
pheatmap, vegan, dada2). Krona: KronaTools ktImportText. FastQC / fastp: qc_tools.py.
"""
from __future__ import annotations

import base64
import csv
import html
import json
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path

import amp_data as AD
import qc_tools as Q

HERE = Path(__file__).resolve().parent
RSCRIPT = HERE / "r_scripts" / "custom_report_plots.R"
RANK_FILES = [("Phylum", "phy", "phylum"), ("Class", "cls", "class"), ("Order", "ord", "order"),
              ("Family", "fam", "family"), ("Genus", "gen", "genus"), ("Species", None, "species")]
PRIMER_NAMES = {("TTCCTAGCCATACAYTAYAC", "GGTGGCKCCTCAGAAGGACATTTGKCCYCA"): ("Fish-CytB-L14912", "Fish-CytB-H15149"),
                ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCCAGTTTG"): ("MiFish-U-F", "MiFish-U-R"),
                ("TTAGATACCCCACTATGC", "TAGAACAGGCTCCTCTAG"): ("Vert-12SV5-F", "Vert-12SV5-R")}
REFS = [
    "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics.",
    "Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34:i884–i890.",
    "Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17:10–12.",
    "Callahan BJ, McMurdie PJ, Rosen MJ, et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature Methods 13:581–583.",
    "McMurdie PJ, Holmes S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE 8:e61217.",
    "Camacho C, Coulouris G, Avagyan V, et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics 10:421.",
    "Leray M, Knowlton N, Machida RJ. (2022). MIDORI2: a collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. Environmental DNA 4:894–907.",
    "Zhu T, Sato Y, Sado T, Miya M, Iwasaki W. (2023). MitoFish, MitoAnnotator, and MiFish Pipeline: updates in 10 years. Molecular Biology and Evolution 40:msad035.",
    "Ondov BD, Bergman NH, Phillippy AM. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12:385.",
    "Oksanen J, et al. vegan: Community Ecology Package. R package.",
    "Shetty SA, Lahti L. microbiomeutilities: utilities for microbiome analytics. R package.",
    "Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer.",
]
e = html.escape


def safe(name: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", name).strip("_") or "sample"


def _png_size(p: Path) -> tuple[int, int]:
    b = p.read_bytes()[:32]
    if b[:8] == b"\x89PNG\r\n\x1a\n":
        return struct.unpack(">II", b[16:24])
    return (800, 600)


def _fig(p: Path, width_cm: float = 15.0) -> str:
    """A PNG inside an SVG: weasyprint (PDF) and docx_export (Word) both carry it."""
    if not p.exists():
        return f'<p class="nofig">Figure not available: {e(p.name)}</p>'
    w, h = _png_size(p)
    w_cm = min(width_cm, 15.0)
    h_cm = w_cm * h / max(w, 1)
    if h_cm > 19:
        w_cm, h_cm = w_cm * 19 / h_cm, 19
    data = base64.b64encode(p.read_bytes()).decode()
    return (f'<div class="fig"><svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{w_cm:.2f}cm" height="{h_cm:.2f}cm" viewBox="0 0 {w} {h}"><image width="{w}" height="{h}" '
            f'xlink:href="data:image/png;base64,{data}" href="data:image/png;base64,{data}"/></svg></div>')


def _cap(rel: str) -> str:
    return f'<p class="cap">&ldquo;{e(rel)}&rdquo;</p>'


def _table(head: list[str], rows: list[list], cls: str = "") -> str:
    th = "".join(f"<th>{e(str(h))}</th>" for h in head)
    tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>'


def _find_reads(upload_dir: Path, sample: str) -> tuple[Path | None, Path | None]:
    if not upload_dir or not Path(upload_dir).is_dir():
        return None, None
    try:
        import read_qc
        found = read_qc.find_samples(Path(upload_dir))
    except Exception:
        found = {}
    for name, files in found.items():
        if name in (sample, safe(sample)) or safe(name) == safe(sample):
            fs = sorted(files)
            r1 = next((f for f in fs if re.search(r"(_R?1[._]|_1\.f)", Path(f).name)), fs[0] if fs else None)
            r2 = next((f for f in fs if f != r1 and re.search(r"(_R?2[._]|_2\.f)", Path(f).name)), None)
            return (Path(r1) if r1 else None), (Path(r2) if r2 else None)
    return None, None


def _tool(tools: dict, name: str):
    if name in tools:
        return tools[name]
    return shutil.which(name) if name == "Rscript" else Q.find_tool(name)


def _write_outfiles(run: AD.RunData, s: str, d: Path):
    d.mkdir(parents=True, exist_ok=True)
    asvs = run.sample_asvs(s)
    tot = sum(a["counts"][s] for a in asvs) or 1
    with open(d / "expr.asv.fasta", "w") as fh:
        for a in asvs:
            fh.write(f">{a['sv']} {AD.label(a, 'Species')}\n{a['sequence']}\n")
    for rank, short, _ in RANK_FILES[:5]:
        t = run.rank_table(rank, s)
        for kind, f in (("abundance", 1), ("relative_abundance", 1 / tot)):
            with open(d / f"expr.{kind}.ab{short}.txt", "w") as fh:
                fh.write(f"{rank}\t{s}\n" + "".join(f"{k}\t{v * f:.6g}\n" for k, v in t.items()))
    for kind, f in (("abundance", 1), ("relative_abundance", 1 / tot)):
        with open(d / f"expr.{kind}.all.txt", "w") as fh:
            fh.write("\t".join(["ASV"] + AD.RANKS + [s]) + "\n")
            for a in asvs:
                fh.write("\t".join([a["sv"]] + [a["tax"].get(r, "") for r in AD.RANKS]
                                   + [f"{a['counts'][s] * f:.6g}"]) + "\n")
    biom = {"id": s, "format": "Biological Observation Matrix 1.0.0", "format_url": "http://biom-format.org",
            "type": "OTU table", "generated_by": "custom_report.py", "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "matrix_type": "sparse", "matrix_element_type": "int", "shape": [len(asvs), 1],
            "rows": [{"id": a["sv"], "metadata": {"taxonomy": [a["tax"].get(r, "") for r in AD.RANKS]}} for a in asvs],
            "columns": [{"id": s, "metadata": None}],
            "data": [[i, 0, a["counts"][s]] for i, a in enumerate(asvs)]}
    (d / "expr.biom").write_text(json.dumps(biom))


def _prepare(run: AD.RunData, s: str, wd: Path, raw: tuple, filt: tuple):
    wd.mkdir(parents=True, exist_ok=True)
    t = run.tracking[s]
    with open(wd / "track.tsv", "w") as fh:
        fh.write("step\tcount\n")
        for k, lab in (("trimmed", "trimmed"), ("filtered", "filtered"), ("denoisedR", "denoisedR"),
                       ("denoisedF", "denoisedF"), ("merged", "merged"), ("Nonchim", "nonchim")):
            fh.write(f"{lab}\t{t[k]}\n")
    with open(wd / "asv.tsv", "w") as fh:
        fh.write("\t".join(["SV", "label"] + AD.RANKS + ["count"]) + "\n")
        for a in run.sample_asvs(s):
            fh.write("\t".join([a["sv"], f"{a['sv']}:{AD.label(a, 'Species')}"]
                               + [AD.label(a, r) for r in AD.RANKS] + [str(a["counts"][s])]) + "\n")
    with open(wd / "krona.tsv", "w") as fh:
        for a in run.sample_asvs(s):
            fh.write("\t".join([str(a["counts"][s])] + [AD.label(a, r) for r in AD.RANKS]) + "\n")
    (wd / "reads.json").write_text(json.dumps({"raw": [str(x) for x in raw if x], "filtered": [str(x) for x in filt if x]}))


CSS = """
@page { size: A4; margin: 2cm 2cm 2.2cm 2cm;
  @bottom-center { content: "Page " counter(page) " of " counter(pages); font-size: 9pt; color: #555; } }
body { font-family: "DejaVu Sans", Arial, sans-serif; font-size: 10.5pt; line-height: 1.45; color: #111; }
h1 { font-size: 15pt; text-align: center; margin: 0 0 12pt; }
h2 { font-size: 13pt; margin: 18pt 0 6pt; page-break-before: always; }
h3 { font-size: 11pt; margin: 12pt 0 4pt; }
h4 { font-size: 10.5pt; margin: 10pt 0 3pt; text-decoration: underline; font-weight: normal; }
p { margin: 3pt 0 6pt; text-align: justify; }
p.res { margin-left: 0; } .ul { text-decoration: underline; font-weight: bold; }
p.cap { text-align: center; font-style: italic; font-size: 9.5pt; }
p.note { color: #b00000; text-align: center; font-style: italic; }
p.nofig { text-align: center; color: #777; border: 1px dashed #bbb; padding: 18pt; }
div.fig { text-align: center; margin: 6pt 0; }
table { border-collapse: collapse; width: 100%; margin: 6pt 0; font-size: 9pt; table-layout: fixed; }
table.num { font-size: 7.6pt; } td, th { overflow-wrap: anywhere; }
th, td { border: 1px solid #999; padding: 3pt 5pt; text-align: left; vertical-align: top; }
th { background: #eee; }
table.cover td { height: 18pt; } table.cover td:first-child { width: 45%; font-weight: bold; }
table.num td + td + td { text-align: right; }
ul.files { font-size: 9.5pt; list-style: none; padding-left: 0; }
.wf td { border: 1px solid #555; } .wf td:first-child { width: 28%; font-weight: bold; text-align: center; }
.indent { margin-left: 18pt; }
"""


def _cover(run: AD.RunData, n: int, ref_text: str) -> list[str]:
    blank = lambda k: f"<tr><td>{e(k)}</td><td></td></tr>"
    o = ["<h1>Custom Metagenomic Sequencing Report Summary</h1>",
         '<table class="cover"><tbody>' + "".join(blank(k) for k in
            ("Client Name", "Order Number", "Customer Email", "Number of Samples to Process",
             "Sequencing Service Provided")) + "</tbody></table>",
         '<table class="cover"><tbody>' + "".join(blank(k) for k in
            ("Download Link:", "Zip file name:", "Password zip file:", "Expiration of link download:")) + "</tbody></table>",
         "<h3>Downloading Your Files</h3>",
         "<p>All files of this service have been uploaded to our repository and can be downloaded for up to 1 month "
         "after the service is completed. Please download every file before the link expires: open the download "
         "link above, enter the user name and password provided, and click Download. If you cannot reach the files, "
         "please contact us. All fields are case sensitive.</p>",
         '<table class="cover"><tbody>' + blank("Date and Time of Release of Report:") + blank("Report by:")
         + blank("Date:") + "</tbody></table>",
         '<table class="cover"><tbody>'
         f"<tr><td>Metagenomic Data Analysis Workflow Used</td><td>See Workflow Analysis</td></tr>"
         f"<tr><td>Sources of Metagenomic Reference Sequencing</td><td>Custom primers</td></tr>"
         f"<tr><td>Reference or Database Sequencing Used</td><td>{e(ref_text)}</td></tr>"
         f"<tr><td>Number of samples in this run</td><td>{n}</td></tr></tbody></table>"]
    o.append('<h2>Data basis analysis:</h2><ul>' + "".join(f"<li>{e(x)}</li>" for x in (
        "Clean data", "Data quality control", "Venn diagram", "Taxonomic profiles (Bar chart & Box plot)",
        "Krona plot", "TOP 50 ASVs Heatmap of relative abundance", "Rarefaction curves", "Alpha diversity")) + "</ul>")
    wf = [("Raw Data", []),
          ("Clean Data", ["Quality control checks on raw sequence data", "Adapter trimming", "Primer removal"]),
          ("DADA2 Pipeline", ["Filter the forward and reverse reads", "Learn the Error Rates", "Merge paired reads",
                              "Remove chimeras", "Assign taxonomy"]),
          ("R / Phyloseq", ["Plot bar chart & box plot of taxonomic profiles", "Krona plot", "Plot Top 50 ASVs Heatmap",
                            "Plot Rarefaction curves", "Alpha diversity"])]
    rows = []
    for i, (k, v) in enumerate(wf):
        rows.append(f"<tr><td>{e(k)}</td><td>{'<br>'.join('- ' + e(x) for x in v)}</td></tr>")
        if i < len(wf) - 1:
            rows.append('<tr><td style="border:none">&#8595;</td><td style="border:none"></td></tr>')
    o.append('<h3>Workflow Analysis</h3><table class="wf"><tbody>' + "".join(rows) + "</tbody></table>")
    return o


LOCAL_LIMIT = {"fish": "fishes (ray-finned, cartilaginous, jawless, lungfish and coelacanth classes)",
               "vertebrates": "vertebrates (Chordata except tunicates and lancelets)", "all": "no taxonomic group"}
NCBI_LIMIT = {"fish": "ray-finned and cartilaginous fishes (NCBI txid 7898 and txid 7777)",
              "vertebrates": "vertebrates (NCBI txid 7742)", "all": "no taxonomic limit"}


def _n(k: int, one: str, many: str) -> str:
    return f"{k} {one if k == 1 else many}"


def _ncbi_sentence(bi: dict, lim: str) -> str:
    ag = bi.get("agreement", {})
    lab = bi.get("ncbi_limit_label") or NCBI_LIMIT.get(lim, "no taxonomic limit")
    parts = [f"{ag.get('agree', 0)} agreed at species",
             f"{ag.get('ncbi_tie', 0)} tied in NCBI between two or more species (cut to the shared rank and flagged)",
             f"{ag.get('differ', 0)} differed (cut to the shared rank and flagged for checking)",
             f"{ag.get('local_better', 0)} matched better in the local reference (local name kept)",
             f"{ag.get('ncbi_only', 0)} were named from NCBI only",
             f"{ag.get('ncbi_no_hit', 0)} had no NCBI hit",
             f"{ag.get('not_checked', 0)} were not checked"]
    t = (f" Every ASV was also checked against the NCBI nucleotide collection (nt) by remote BLASTn limited to "
         f"{lab}: " + "; ".join(parts) + ".")
    if bi.get("ncbi") == "partial":
        nb = bi.get("ncbi_batches", {})
        t += (f" {nb.get('failed_batches', 0)} of {nb.get('batches', 0)} NCBI batches got no answer; their ASVs are "
              "counted as not checked and rest on the local reference.")
    return t


def _dada2_text(p: dict) -> tuple[str, str]:
    db = p.get("dbPath", "") or ""
    return (f"Classify sequences with the DADA2 naive Bayesian classifier (assignTaxonomy, minBoot="
            f"{p.get('minBoot', 50)}) against {e(Path(db).name or 'the reference database')}.", Path(db).name or "—")


def _tax_text(run: AD.RunData) -> tuple[str, str]:
    p, bi = run.params, run.blast_info
    gene = "CytB" if str(p.get("mito_gene", "")).upper() == "CYTB" else "12S rRNA"
    if bi and (run.dir / "blast" / "taxonomy_blast.csv").exists():
        rf = bi.get("reference", {})
        v = rf.get("versions", {})
        has_mf = bool(v.get("MitoFish")) and int(rf.get("n_mitofish", 1) or 0) > 0
        lim = bi.get("limit", "fish")
        th = bi.get("thresholds", {})
        src_ref = f"MIDORI2 {v.get('MIDORI2', '')}" + (f" + MitoFish {v.get('MitoFish', '')}" if has_mf else "")
        ref = f"{src_ref} ({gene}, BLASTn)"
        t = (f"Classify sequences by Nucleotide BLAST ({e(bi.get('blastn', 'BLASTn'))}) against a fish mitochondrial "
             f"{gene} reference built from MIDORI2 {e(v.get('MIDORI2', ''))} (GenBank, quality-controlled; Leray et al. "
             "2022)" + (f" and MitoFish {e(v.get('MitoFish', ''))} (Zhu et al. 2023); records MitoFish judges "
                        "heterospecific (mislabelled) were removed from both sources" if has_mf else "")
             + "; records with more than 1% ambiguous bases were removed, and a sequence filed under several species "
             "is kept once and named only to the rank those species share. "
             f"Names were limited to {LOCAL_LIMIT.get(lim, 'no taxonomic group')}. "
             "Each ASV is named by the lowest rank shared by all hits within 0.5% identity of the best hit, and only "
             f"down to the deepest rank whose identity threshold is met (species ≥ {th.get('Species', 98.5)}%, genus ≥ "
             f"{th.get('Genus', 95)}%, family ≥ {th.get('Family', 90)}%, order ≥ {th.get('Order', 85)}%, class ≥ "
             f"{th.get('Class', 80)}%). ASVs whose best hit lies outside the limit keep their real name and are reported "
             "as non-target.")
        if bi.get("ncbi") in ("done", "partial"):
            t += _ncbi_sentence(bi, lim)
            ref += f"; NCBI nt check ({bi.get('ncbi_limit_label') or NCBI_LIMIT.get(lim, '')})"
        elif bi.get("ncbi") == "failed":
            t += " The NCBI nt check could not be completed (no answer from NCBI); the names rest on the local reference."
        return t, ref
    t, ref = _dada2_text(p)
    if str(p.get("tax_method", "dada2")) in ("blast", "both"):
        why = ""
        sk = run.dir / "blast" / "BLAST_SKIPPED.txt"
        if sk.exists():
            why = sk.read_text(errors="replace").strip().splitlines()[0][:300] if sk.read_text().strip() else ""
        t += (" BLAST naming was requested but not used" + (f" ({e(why)})" if why else "")
              + "; the names in this report come from the DADA2 classifier above.")
    return t, ref


MISSING_RAW = "The raw reads of this sample were not found in the job's upload folder; {tool} was not run on them."


def _num(x, fmt="{:,}") -> str:
    return "unknown" if x is None else fmt.format(x)


def _xlsx_cap(img: Path, rel, job_label: str, i: int) -> str:
    """Caption for the sample's Excel file, only when it was really written."""
    name = f"{job_label}_sample{i}.xlsx"
    if (Path(img) / name).exists():
        return _cap(rel(name))
    return ('<p class="note">The Excel file was not made on this machine (Python package openpyxl is missing; '
            "install it with: bash setup_amplicon_report.sh).</p>")


def _sample_block(run, i, s, folder, img, qc, plots, alpha_r, job_label) -> list[str]:
    p = run.params
    t = run.tracking[s]
    rel = lambda x: f"{folder}/Images_and_Outfiles/{x}"
    f_, r_ = p.get("primer_f", ""), p.get("primer_r", "")
    fn, rn = PRIMER_NAMES.get((f_, r_), ("Forward primer", "Reverse primer"))
    o = ['<h2>Report Summary:</h2>', f"<h3>{i}. Sample{i}_{e(s)}</h3>",
         "<h3>Clean data</h3>", "<h4>1. Quality control checks on raw sequence data by FastQC</h4>",
         "<p>FastQC" + (f" {e(qc['fastqc']['version'])}" if qc.get("fastqc") else "") + " reports per-base quality, "
         "sequence content, duplication and adapter content for each raw FASTQ file.</p>",
         f'<p><span class="ul">Result:</span> {e(MISSING_RAW.format(tool="FastQC") if qc.get("missing") else Q.sentence_fastqc(qc.get("fastqc")))}</p>',
         "<h4>2. Adapter trimming by fastp</h4>",
         "<p>fastp" + (f" {e(qc['fastp'].get('version', ''))}" if qc.get("fastp") else "") + " detects adapters (for read "
         "pairs, from the overlap of the two reads). It was run on the raw reads as a check only: its trimmed reads are "
         "not used further, because primer removal and read-through trimming for the analysis are done by cutadapt "
         "(next step).</p>",
         f'<p><span class="ul">Result:</span> {e(MISSING_RAW.format(tool="fastp") if qc.get("missing") else Q.sentence_fastp(qc.get("fastp")))}</p>',
         "<h4>3. Remove primer sequences by cutadapt" + (f" v.{e(t['cutadapt_version'])}" if t.get("cutadapt_version") else "")
         + "</h4>",
         "<p>cutadapt searched each read pair for the forward primer at the start of read 1 and the reverse primer at "
         f"the start of read 2 (error rate {p.get('cutadaptErrorRate', 0.1)}, minimum overlap {p.get('cutadaptOverlap', 3)} "
         "bases), removed them and their reverse complements at the read ends, and dropped reads shorter than 50 bases."
         "</p>"]
    if f_ and r_:
        hits = t.get("primer_hits", 0)
        res = (f'The primers {e(fn)}="{e(f_)}" and {e(rn)}="{e(r_)}" were removed'
               + (f" ({hits:,} of {t['raw_data']:,} read pairs, {100 * hits / max(t['raw_data'], 1):.1f}%, carried the "
                  "forward primer)" if hits and t.get("raw_data") else "") + ".")
    else:
        res = "No primer sequences were given; primer removal was skipped."
    o.append(f'<p><span class="ul">Result:</span> {res}</p>')
    trF, trR = int(p.get("truncLen_F") or 0), int(p.get("truncLen_R") or 0)
    trunc = (f"trimming the forward reads at {trF} and the reverse reads at {trR} based on the forward / reverse "
             "read quality profiles" if trF > 0 else
             "no truncation (truncLen = 0): reads were kept at their full length after primer removal")
    dver = plots.get("dada2", "")
    paired = any((x.get("denoisedR") or 0) > 0 for x in run.tracking.values()) or bool(p.get("truncLen_R"))
    pool = str(p.get("pool", "FALSE")).upper()
    pool_t = {"TRUE": "all samples were pooled for denoising (pool=TRUE)",
              "PSEUDO": "samples were denoised with pseudo-pooling (pool=\"pseudo\")"}.get(
        pool, "each sample was then denoised on its own (pool=FALSE)")
    o += ["<h3>DADA2 pipeline</h3>", f"<p>DADA2{(' ' + e(dver)) if dver else ''} in R.</p>",
          "<h4>1. Pre-processing of reads for amplicon analysis</h4>",
          f'<p class="indent">1) Filter the forward and reverse reads -- using filterAndTrim with maxN=0 (DADA2 requires no '
          f"Ns), truncQ=2, rm.phix={'TRUE' if paired else 'FALSE'}, maxEE={float(p.get('maxEE_F', 2)):g},{float(p.get('maxEE_R', 2)):g} and "
          f"{trunc}.</p>",
          '<p class="indent">Output:</p><ul class="files">' + "".join(
              f"<li>- &ldquo;{e(x)}&rdquo;</li>" for x in (rel("plotQualityProfile.pdf"),
                                                            rel("plotQualityProfile_filtered.pdf"),
                                                            f"{folder}/filtered_reads")) + "</ul>",
          '<p class="indent">2) Learn the Error Rates -- learnErrors estimates the error model of this sequencing run '
          "from the data, alternating between estimating error rates and inferring the sample composition until they "
          f"converge. The error model was learned from the reads of this run's samples together; {pool_t}.</p>",
          '<p class="indent">3) Merge paired reads -- each denoised forward read is aligned with the reverse complement of '
          "its denoised reverse read; pairs are merged when they overlap by at least 12 bases with no mismatch in the "
          "overlap.</p>",
          f'<p class="indent">4) Chimera Removal -- removeBimeraDenovo (method "{e(str(p.get("chimeraMethod", "consensus")))}") '
          "removes sequence variants that are bimeras of two more abundant parents.</p>",
          '<p><span class="ul">Result:</span> <b>Data quality control</b></p>',
          _fig(img / "Box_plot.png", 14), _cap(rel("Box_plot.png")),
          _table(["SampleID", "Sample_name", "raw_data", "trimmed", "filtered", "denoisedF", "denoisedR", "merged",
                  "Nonchim", "perc"],
                 [[e(s), e(s)] + [_num(t[k]) for k in ("raw_data", "trimmed", "filtered", "denoisedF", "denoisedR",
                                                       "merged", "Nonchim")] + [_num(t["perc"], "{:.2f}")]], "num"),
          _xlsx_cap(img, rel, job_label, i)]
    if t.get("raw_data") is None:
        o.append('<p class="note">raw_data is unknown for this sample: cutadapt\'s own read counts were not found, so '
                 "the count before primer removal cannot be given (trimmed is the count after it).</p>")
    if s in getattr(run, "dropped", []):
        o.append('<p class="note">No reads of this sample were left after filtering: the pipeline dropped it from the '
                 "ASV table.</p>")
    if t["Nonchim"] == 0:
        o.append('<p class="note">No reads of this sample passed filtering, merging and chimera removal; the profiles '
                 "below are empty. Check the library or sequencing depth.</p>")
    tax_t, _ = _tax_text(run)
    nt = run.non_target_share(s)
    nts = sorted({AD.label(a, "Species") for a in run.sample_asvs(s) if a.get("target") == "no"})
    o += ["<h3>2. Taxonomic assignment</h3>", f"<p>{tax_t}</p>",
          ("<p><b>Non-target reads:</b> " + (f"{nt:.1f}% of this sample's reads came from outside the limit "
                                              f"({e(', '.join(nts))}) — for example human DNA that the blocking primer did "
                                              "not stop. They are kept in the tables and marked." if nt > 0 else
                                              "none; no read was assigned outside the target group.") + "</p>")]
    chk = [a for a in run.sample_asvs(s) if a.get("check") == "yes"]
    if chk:
        o.append("<p><b>ASVs to check</b> (local reference and NCBI disagree; named only to the shared rank):</p>"
                 + _table(["ASV", "Name given", "Local best hit (% id)", "NCBI best hit (% id)", "Reads"],
                          [[a["sv"], e(AD.label(a, "Species")), f"{e(a.get('best_hit', ''))} ({a.get('pident', '')})",
                            f"{e(a.get('ncbi_best', ''))} ({a.get('ncbi_pident', '')})", f"{a['counts'][s]:,}"]
                           for a in chk]))
    hm = plots.get("heatmap_tool", "pheatmap")
    o += ["<h3>Explore microbiome profiles using phyloseq package in R</h3>",
          f"<p>{e(plots.get('phyloseq') or 'phyloseq')} holds the ASV counts, their taxonomy and the sample data as one object and hands them to the R "
          "packages that draw and test them (ggplot2 for the figures, vegan for rarefaction and diversity).</p>",
          "<h4>1. Venn diagram</h4>", "<p>Result:</p>",
          '<p class="note">It is not possible to draw a Venn diagram in a single sample.</p>',
          "<h4>2. Bar chart of taxonomic profiles by ggplot2 package</h4>",
          "<p>Relative abundance of the taxa detected in this sample at each taxonomic rank. ASVs without a name at a "
          "rank are shown as '&lt;higher rank&gt; unidentified'.</p>"]
    for rank, _, low in RANK_FILES:
        o += [f"<p><u>At {rank} level</u></p>", _fig(img / "Taxonomic profiles" / f"plot_bar_{low}.png", 9),
              _cap(rel(f"Taxonomic profiles/plot_bar_{low}.png"))]
    files = ["expr.asv.fasta"] + [f"expr.abundance.ab{x}.txt" for _, x, _ in RANK_FILES[:5]] + ["expr.abundance.all.txt"] \
        + [f"expr.relative_abundance.ab{x}.txt" for _, x, _ in RANK_FILES[:5]] + ["expr.relative_abundance.all.txt",
                                                                                 "expr.biom"]
    o.append(f"<p><b><u>All File</u></b> {e(folder)}</p><ul class=\"files\">"
             + "".join(f"<li>- Images_and_Outfiles/Taxonomic profiles/outfiles/{f}</li>" for f in files) + "</ul>")
    o += [f"<h4>3. Top 50 ASVs Heatmap by {e(hm.split(' ')[0])} package</h4>", "<p>Result:</p>",
          _fig(img / "Heatmap" / "Heat_map_label_ASV_number_top50ASVs.png", 13),
          _cap(rel("Heatmap/Heat_map_label_ASV_number_top50ASVs.png")),
          "<h4>4. Krona plot</h4>",
          "<p>Krona shows the taxonomy as nested rings that can be zoomed in a web browser; the interactive chart is the "
          "HTML file below, the figure is a still of it.</p>", "<p>Result:</p>",
          _fig(img / "Krona plot" / "krona_static.png", 12), _cap(rel("Krona plot/SampleID-krona.html")),
          "<h4>5. Statistical analysis</h4>",
          '<p class="indent">1) Rarefaction analysis shows how the number of ASVs found grows with the number of reads '
          "sampled; a curve that levels off means the sequencing depth captured the sample's diversity.</p>",
          "<p>Result:</p>", _fig(img / "Rarefaction curves" / "Rarefaction_curves.png", 13),
          _cap(rel("Rarefaction curves/Rarefaction_curves.png")),
          '<p class="indent">2) Alpha diversity measures the diversity within one sample: richness (Observed, Chao1), '
          "and diversity weighing evenness (Shannon, Simpson).</p>",
          "<p>Result: A comparison between groups (Wilcoxon test) needs two or more samples per group and is not "
          "made within a single-sample section.</p>", "<p><b>Estimate richness</b></p>"]
    a = alpha_r.get(s) or run.alpha[s]
    o.append(_table(["Sample_ID", "Sample_name", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson"],
                    [[e(s), e(s)] + [f"{float(a[k]):.2f}" for k in ("Observed", "Chao1", "se.chao1", "Shannon",
                                                                    "Simpson")]], "num"))
    o += ["<p>Please see more results in the Excel file, which includes the taxonomic profiles and the sequence of "
          "each taxon.</p>", _xlsx_cap(img, rel, job_label, i)]
    return o


def build(run_dir: Path, upload_dir: Path | None, out_dir: Path, tools: dict | None = None,
          formats=("html", "pdf", "docx"), job_label: str = "report") -> dict:
    tools = dict(tools or {})
    run = AD.load_run(Path(run_dir))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_work"
    work.mkdir(exist_ok=True)
    job_label = safe(job_label)
    folders, jobs = {}, []
    for i, s in enumerate(run.samples, 1):
        folder = f"Sample{i}_{safe(s)}"
        folders[s] = folder
        img = out / folder / "Images_and_Outfiles"
        for sub in ("Taxonomic profiles/outfiles", "Heatmap", "Krona plot", "Rarefaction curves"):
            (img / sub).mkdir(parents=True, exist_ok=True)
        raw = _find_reads(upload_dir, s)
        filt = tuple(p for p in (run.dir / "filtered" / f"{s}_F_filt.fastq.gz", run.dir / "filtered" / f"{s}_R_filt.fastq.gz")
                     if p.exists())
        fr = out / folder / "filtered_reads"
        fr.mkdir(exist_ok=True)
        for p in filt:
            try:
                os.link(p, fr / p.name)
            except OSError:
                shutil.copy(p, fr / p.name)
        _write_outfiles(run, s, img / "Taxonomic profiles" / "outfiles")
        _prepare(run, s, work / str(i), raw, filt)
        jobs.append([str(work / str(i)), str(img), s])
        try:
            import amplicon_excel as AX
            AX.write_sample(img / f"{job_label}_sample{i}.xlsx", run, s)
        except ImportError:
            pass
    try:
        import amplicon_excel as AX
        AX.write_all(out / f"{job_label}_all_samples.xlsx", run)
    except ImportError:
        pass
    # figures (R) and Krona
    with open(work / "jobs.tsv", "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(jobs)
    rs = _tool(tools, "Rscript")
    if rs:
        try:
            subprocess.run([rs, str(RSCRIPT), "--jobs", str(work / "jobs.tsv")], capture_output=True, text=True,
                           timeout=3 * 3600)
        except subprocess.TimeoutExpired:
            pass                                       # missing figures are shown as "Figure not available"
    kt = _tool(tools, "ktImportText")
    for wd, img, s in jobs:
        html_out = Path(img) / "Krona plot" / "SampleID-krona.html"
        if kt:
            try:
                subprocess.run([kt, "-o", str(html_out), f"{Path(wd) / 'krona.tsv'},{s}"], capture_output=True,
                               timeout=600)
            except subprocess.TimeoutExpired:
                pass
        if not html_out.exists():
            html_out.write_text("<!doctype html><meta charset=utf-8><title>" + e(s) + "</title><body style='text-align:center;"
                                "font-family:sans-serif'><h3>" + e(s) + " — taxonomy (KronaTools not installed: still "
                                "image)</h3>" + _fig(Path(img) / "Krona plot" / "krona_static.png", 20) + "</body>")
    # FastQC / fastp on the raw reads
    qc = {}
    for i, s in enumerate(run.samples, 1):
        r1, r2 = _find_reads(upload_dir, s)
        q = {} if r1 else {"missing": True}
        if r1:
            q["fastqc"] = Q.fastqc([x for x in (r1, r2) if x], work / str(i) / "fastqc", _tool(tools, "fastqc"))
            q["fastp"] = Q.fastp_check(r1, r2, work / str(i) / "fastp", _tool(tools, "fastp"))
        qc[s] = q
    plots, alpha_r = {}, {}
    for wd, img, s in jobs:
        try:
            plots = {**plots, **json.loads((Path(wd) / "plots_info.json").read_text())}
        except Exception:
            pass
        ap = Path(wd) / "alpha.tsv"
        if ap.exists():
            with open(ap) as fh:
                r = next(csv.DictReader(fh, delimiter="\t"), None)
                if r:
                    alpha_r[s] = r
    _, ref_text = _tax_text(run)
    body = _cover(run, len(run.samples), ref_text)
    for i, s in enumerate(run.samples, 1):
        body += _sample_block(run, i, s, folders[s], out / folders[s] / "Images_and_Outfiles", qc[s], plots, alpha_r,
                              job_label)
    refs = [r for r in REFS if not (("MitoFish" in r and "MitoFish" not in ref_text)
                                     or ("BLAST+" in r and "BLASTn" not in ref_text))]
    body += ['<h2>References</h2><ol class="files">' + "".join(f"<li>[{k}] {e(r)}</li>" for k, r in enumerate(refs, 1))
             + "</ol>"]
    doc = ("<!doctype html><html><head><meta charset=\"utf-8\"><title>Custom Metagenomic Sequencing Report</title>"
           f"<style>{CSS}</style></head><body>" + "\n".join(body) + "</body></html>")
    res = {"samples": run.samples, "folders": folders}
    hp = out / f"{job_label}_report.html"
    hp.write_text(doc, encoding="utf-8")
    res["html"] = str(hp)
    if "pdf" in formats:
        try:
            from weasyprint import HTML as _WH
            _WH(string=doc, base_url=str(out)).write_pdf(str(out / f"{job_label}_report.pdf"))
            res["pdf"] = str(out / f"{job_label}_report.pdf")
        except Exception as ex:
            res["pdf_error"] = str(ex)
    if "docx" in formats:
        try:
            import docx_export
            docx_export.html_to_docx(doc, out / f"{job_label}_report.docx", title="Custom Metagenomic Sequencing Report")
            res["docx"] = str(out / f"{job_label}_report.docx")
        except Exception as ex:
            res["docx_error"] = str(ex)
    shutil.rmtree(work, ignore_errors=True)
    z = shutil.make_archive(str(out.parent / f"{out.name}"), "zip", root_dir=str(out))
    res["zip"] = z
    (out / "report_info.json").write_text(json.dumps({k: v for k, v in res.items() if k != "folders"}, indent=1))
    return res
