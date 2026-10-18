#!/usr/bin/env python3
"""shotgun_report.py — the HTML report for a shotgun metagenome run (beta).

Called by report_builder.build_profile_html() when Run.kind == "shotgun". It
reuses report_builder's analysis and figures (composition, heatmap, alpha and
beta diversity, Controls, MaAsLin 3) but has its own Method and QC text: no
primers, amplicons or ASVs, and abundances are reads per species.
"""
from __future__ import annotations

import report_builder as RB
from report_builder import esc, fmt, Analysis, _tbl, _fig

REFS = [
    ("Wood DE, Lu J, Langmead B. (2019).", "Improved metagenomic analysis with Kraken 2.",
     "Genome Biology 20:257."),
    ("Lu J, Breitwieser FP, Thielen P, Salzberg SL. (2017).",
     "Bracken: estimating species abundance in metagenomics data.", "PeerJ Computer Science 3:e104."),
    ("Chen S, Zhou Y, Chen Y, Gu J. (2018).", "fastp: an ultra-fast all-in-one FASTQ preprocessor.",
     "Bioinformatics 34:i884–i890."),
    ("Constantinides B, Hunt M, Crook DW. (2023).", "Hostile: accurate decontamination of microbial "
     "host sequences.", "Bioinformatics 39:btad728."),
    ("De Coster W, Rademakers R. (2023).", "NanoPack2: population-scale evaluation of long-read "
     "sequencing data.", "Bioinformatics 39:btad311."),
]


def build(run, title=None, app="NextGen-Amplicon", company="", subtitle=None) -> str:
    an = Analysis(run)
    S = an.samples
    sm, p = run.summary, run.params
    platform = sm.get("platform") or ("ont" if p.get("sequencerType") == "ont" else "illumina")
    plat_txt = "Oxford Nanopore" if platform == "ont" else "Illumina paired-end"
    rep_title = "Shotgun Metagenome Sequencing Report"
    project = title if title and title.strip() and title.strip() != rep_title else None
    subtitle = subtitle or (f"{run.label} · {len(S)} sample{'s' if len(S) != 1 else ''} · {plat_txt}, "
                            "Kraken2 + Bracken · beta")
    fig, tab = [0], [0]

    def F(svg, cap, folder=None):
        if not svg:
            return ""
        fig[0] += 1
        return _fig(svg, cap, folder, fig[0])

    def T(head, rows, cap):
        tab[0] += 1
        return f'<p class="tcap"><b>Table {tab[0]}.</b> {cap}</p>' + _tbl(head, rows)

    tr = run.tracking or []
    num = lambda r, k: float(r.get(k) or 0)
    inp = sum(num(r, "input") for r in tr)
    fin = sum(num(r, "final") for r in tr) or run.total
    clas = sum(num(r, "classified") for r in tr)
    after = sum(num(r, "after_qc") for r in tr)

    out = [f'<div class="eyebrow">{esc(run.label)} &middot; shotgun metagenome &middot; '
           f'{len(S)} sample{"s" if len(S) != 1 else ""}</div>', f"<h1>{esc(rep_title)}</h1>"]
    if project:
        out.append(f'<p class="project">{esc(project)}</p>')
    out.append(f'<p class="sub">{esc(subtitle)}</p>')
    kpis = [(fmt(inp) if inp else "—", "read pairs sequenced" if platform == "illumina" else "reads sequenced"),
            (f"{100 * clas / after:.1f}%" if after else "—", "of QC-passed reads classified"),
            (fmt(fin), "reads assigned to species"),
            (str(run.n_asvs), "species detected")]
    out.append('<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>' for v, l in kpis) + "</div>")
    out.append('<div class="note"><p style="margin:0"><b>Beta pipeline.</b> Shotgun profiling in this app is '
               "new. Results are produced with established tools and settings, but the pipeline has been "
               "validated on simulated and public mock-community data only, not yet on this laboratory's "
               "own samples.</p></div>")
    if an.groups or an.desc:
        out.append(_tbl(["Sample", "Group", "Description"],
                        [[esc(s), esc(an.groups.get(s, "—")), f'<span class="dim">{esc(an.desc.get(s, ""))}</span>']
                         for s in S]))

    # ═══ METHOD ═══
    out.append("<h2>Method</h2>")
    tools = sm.get("tools") or {}
    v = lambda k: f" {tools[k]}" if tools.get(k) else ""
    host = sm.get("host") or p.get("sg_host") or "none"
    host_txt = {"none": "Host reads were not removed before classification",
                "human": "Human reads were removed with Hostile" + v("hostile") + " against the T2T-CHM13 "
                         "human genome with HLA alleles",
                "mouse": "Mouse reads were removed with Hostile" + v("hostile") + " against GRCm39"}.get(
        host, f"Host reads were removed with Hostile{v('hostile')} against {esc(host)}")
    qc_txt = (f"Reads were trimmed and filtered with fastp{v('fastp')} (adapter detection for pairs, quality "
              "trimming, minimum length 50 bp)" if platform == "illumina" else
              f"Reads shorter than {p.get('sg_min_len', 500)} bp or with mean quality below Q"
              f"{p.get('sg_min_q', 10)} were removed with "
              + ("chopper" + v("chopper") if tools.get("chopper") else "SeqKit" + v("seqkit")))
    blen = sm.get("bracken_read_length") or {}
    out.append(f"<p>{qc_txt}. {host_txt}. Every remaining read was classified with Kraken 2{v('kraken2')} "
               f"(k-mer matching, confidence threshold {sm.get('confidence', p.get('sg_confidence', 0.1))}) "
               f"against the index <b>{esc(sm.get('kraken_db_name') or '')}</b>. Species abundances were "
               f"re-estimated with Bracken{v('bracken')} at the "
               f"{'species' if sm.get('bracken_level', 'S') == 'S' else 'genus'} level (minimum "
               f"{sm.get('bracken_threshold', 10)} reads"
               + (f"; k-mer distribution for {', '.join(sorted({str(x) for x in blen.values()}))}-bp reads"
                  if blen else "") + "), which redistributes reads that Kraken 2 could only place at "
               "genus or family level. Reads still classified as host (human / mouse) were counted and "
               "left out of the community tables. Diversity, composition and the statistics below use "
               "the Bracken species table.</p>")
    if platform == "ont" and blen and max(int(x) for x in blen.values()) < 1000:
        out.append("<p class=\"dim\">Bracken's k-mer distributions are built for a fixed read length; for "
                   "Nanopore reads the longest available one was used, which is shorter than the reads. "
                   "Species estimates remain usable; genus-level figures are the more robust ones.</p>")
    prov = RB._provenance_rows(run)
    if prov:
        out.append(T(["Software / reference", "Version", "Checksum (MD5)"], prov,
                     "Software and reference versions recorded when the run finished (provenance.json)."))

    # ═══ RESULT ═══
    out.append('<h2 class="pbreak">Result</h2>')
    out.append("<h3>Read processing</h3>")
    rows = []
    for r in tr:
        i0 = num(r, "input")
        rows.append([esc(r.get("sample", "")), fmt(i0), fmt(num(r, "after_qc")),
                     fmt(num(r, "host_removed")) if r.get("host_removed") not in (None, "") else "—",
                     f"{100 * num(r, 'classified') / max(1, num(r, 'after_qc') - num(r, 'host_removed')):.1f}%",
                     fmt(num(r, "host_classified")), fmt(num(r, "final"))])
    out.append(T(["Sample", "Input", "After QC", "Host removed", "Classified", "Host (Kraken 2)",
                  "Species reads"], rows,
                 "Reads at each step. <i>Classified</i> is the share of host-free reads Kraken 2 placed "
                 "anywhere in the taxonomy; <i>species reads</i> are Bracken's estimate of reads per "
                 "species, the basis of every figure below. A low classified share usually means organisms "
                 "the index does not contain (environmental samples), not a failed run."))
    low = [r.get("sample") for r in tr if num(r, "final") < 10000]
    if low:
        out.append(f"<p><b>Low depth:</b> {esc(', '.join(low))} — fewer than 10,000 reads reached species "
                   "level; rare taxa in these samples are likely to be missed.</p>")
    dom = sm.get("domains") or {}
    if dom:
        tot = sum(dom.values()) or 1
        out.append(T(["Domain", "Reads", "Share"],
                     [[esc(k), fmt(v), f"{100 * v / tot:.1f}%"] for k, v in dom.items()],
                     "Species reads by domain across all samples (viruses as the index names them)."))

    out.append('<h3 class="pbreak">Taxonomic composition</h3>')
    for rank in ("Phylum", "Genus", "Species"):
        if rank in run.rank_cols:
            out.append(F(RB.svg_stacked(an, rank, 30),
                         f"Top 30 {rank.lower()} by relative abundance (Bracken species reads summed to "
                         f"{rank.lower()})."))
    if "Genus" in run.rank_cols:
        out.append(F(RB.svg_heatmap(an, "Genus", 30), "Top 30 genera — relative abundance heatmap."))
    def _dom(rank, s):
        return ", ".join(f"{n} ({v:.1f}%)" if isinstance(v, (int, float)) else str(n)
                         for n, v in (d if isinstance(d, tuple) else (d, None) for d in an.dominant(rank, s)))
    dom_rows = [[esc(s), esc(_dom("Phylum", s)), esc(_dom("Species", s))] for s in S]
    out.append(T(["Sample", "Top phyla", "Top species"], dom_rows, "Dominant taxa per sample."))

    out.append('<h3 class="pbreak">Diversity</h3>')
    out.append(F(RB.svg_alpha(an, an.alpha()), "Alpha diversity per sample (species level): observed "
                 "species, Chao1, Shannon and Simpson."))
    out.append(F(RB.svg_rare(an), "Rarefaction of species against sequencing depth."))
    if len(S) >= 3:
        pc = RB.svg_pcoa(an)
        if pc:
            out.append(F(pc, f"PCoA of Bray–Curtis dissimilarity (axis 1 {an.pcoa_pct[0]:.1f}%, axis 2 "
                             f"{an.pcoa_pct[1]:.1f}%)."))
        out.append(F(RB.svg_bc_heat(an), "Bray–Curtis dissimilarity between samples."))

    out += RB._controls_section(run, an, T)
    out += RB._maaslin_section(run, T)

    out.append('<h2 class="pbreak">Notes on reading this report</h2>')
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in [
        "Shotgun read counts reflect DNA, not cells: larger genomes contribute more reads per cell, so "
        "relative abundance is biased towards organisms with big genomes.",
        "Kraken 2 can only name organisms whose genomes, or close relatives, are in the index. Species "
        "near the detection limit (a few hundred reads) should be confirmed before being reported.",
        "Full tables: asv_table.csv (species × sample, estimated reads; first column = NCBI taxid), "
        "taxonomy.csv (lineage per taxid), read_tracking.csv, and per-sample Kraken 2 / Bracken reports in "
        "shotgun/&lt;sample&gt;/.",
    ]) + "</ul>")
    out.append("<h2>References</h2>")
    refs = REFS[:2] + ([REFS[2]] if platform == "illumina" else [REFS[4]]) + \
        ([REFS[3]] if host != "none" else [])
    if (run.path / "r_tables" / "maaslin3_results.csv").exists():
        refs.append(RB.REFERENCES[25])
    out.append('<ul class="refs">' + "".join(f"<li>{esc(a)} <i>{esc(t)}</i> {esc(j)}</li>"
                                             for a, t, j in refs) + "</ul>")
    css = RB._css(app, company, f"{run.label} · shotgun") + RB._CSS_FULL.replace("@@TITLE@@", esc(rep_title))
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{esc(rep_title)}</title>'
            f"<style>{css}</style></head><body>" + "".join(out) + "</body></html>")
