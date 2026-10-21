#!/usr/bin/env python3
"""shotgun_report.py — the HTML report for a shotgun metagenome run (beta).

Called by report_builder.build_profile_html() when Run.kind == "shotgun". It
reuses report_builder's analysis and figures (composition, heatmap, alpha and
beta diversity, Controls, MaAsLin 3) but has its own Method and QC text: no
primers, amplicons or ASVs, and abundances are reads per species.
"""
from __future__ import annotations

import csv
from pathlib import Path

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
REFS_AMR = [
    ("Li H. (2018).", "Minimap2: pairwise alignment for nucleotide sequences.", "Bioinformatics 34:3094–3100."),
    ("Feldgarden M, Brover V, Gonzalez-Escalona N, et al. (2021).", "AMRFinderPlus and the Reference Gene "
     "Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, "
     "and virulence.", "Scientific Reports 11:12728."),
]
REFS_MAG = {
    "megahit": ("Li D, Liu CM, Luo R, Sadakane K, Lam TW. (2015).", "MEGAHIT: an ultra-fast single-node "
                "solution for large and complex metagenomics assembly via succinct de Bruijn graph.",
                "Bioinformatics 31:1674–1676."),
    "flye": ("Kolmogorov M, Bickhart DM, Behsaz B, et al. (2020).", "metaFlye: scalable long-read "
             "metagenome assembly using repeat graphs.", "Nature Methods 17:1103–1110."),
    "metabat2": ("Kang DD, Li F, Kirton E, et al. (2019).", "MetaBAT 2: an adaptive binning algorithm for "
                 "robust and efficient genome reconstruction from metagenome assemblies.", "PeerJ 7:e7359."),
    "checkm2": ("Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. (2023).", "CheckM2: a rapid, scalable and "
                "accurate tool for assessing microbial genome quality using machine learning.",
                "Nature Methods 20:1203–1212."),
    "mimag": ("Bowers RM, Kyrpides NC, Stepanauskas R, et al. (2017).", "Minimum information about a single "
              "amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea.",
              "Nature Biotechnology 35:725–731."),
    "gtdb": ("Parks DH, Chuvochina M, Rinke C, et al. (2022).", "GTDB: an ongoing census of bacterial and "
             "archaeal diversity through a phylogenetically consistent, rank normalized and complete "
             "genome-based taxonomy.", "Nucleic Acids Research 50:D785–D794."),
}


def _csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(v, d=1, dash="—"):
    try:
        return f"{float(v):.{d}f}"
    except (TypeError, ValueError):
        return dash


def _short(name: str) -> str:
    """Escherichia coli -> E. coli (genus alone stays whole)."""
    p = (name or "").split()
    return f"{p[0][0]}. {' '.join(p[1:])}" if len(p) >= 2 and p[0][:1].isupper() else name or ""


def _shade(v: float, vmax: float) -> str:
    """Cell background for a gene x sample RPKM matrix (log-ish ramp)."""
    if v <= 0 or vmax <= 0:
        return ""
    import math
    t = min(1.0, math.log10(1 + v) / math.log10(1 + vmax))
    a = 0.10 + 0.55 * t
    return f' style="background:rgba(196,84,52,{a:.2f})"'


def _resistome(run, sm, S, T, F) -> tuple[list[str], list[str]]:
    """(method sentences, result HTML) for the read-based resistome."""
    info = sm.get("amr") or {}
    rows = _csv(run.path / "amr_genes.csv")
    if not info:
        return [], []
    tools = sm.get("tools") or {}
    mm = f" {tools['minimap2']}" if tools.get("minimap2") else ""
    meth = [f"Antimicrobial resistance (AMR) genes were detected directly in the reads: reads were aligned "
            f"with minimap2{mm} to the nucleotide sequences of the NCBI AMRFinderPlus Reference Gene "
            f"Catalog (database {RB.esc(info.get('db_version', ''))}); alignments of at least "
            f"{_f(info.get('min_identity'), 0)}% identity were kept, reads were pooled per gene family "
            f"(alleles such as blaTEM-1 and blaTEM-116 count as blaTEM) and a family was called present "
            f"when its reads covered at least {_f(info.get('min_breadth'), 0)}% of the gene. Abundance is "
            "RPKM, reads per kilobase of gene per million reads. The organism carrying each gene is the "
            "Kraken 2 call shared by most of the reads that carry it."]
    out = ['<h3 class="pbreak">Antimicrobial resistance genes (resistome)</h3>']
    amr = [r for r in rows if r.get("type") == "AMR"]
    other = [r for r in rows if r.get("type") != "AMR"]
    if not amr:
        out.append("<p>No acquired AMR gene reached the identity and coverage thresholds in any sample.</p>")
    else:
        per = {smp: [r for r in amr if r["sample"] == smp] for smp in S}
        out.append(T(["Sample", "AMR gene families", "Drug classes"],
                     [[RB.esc(smp), str(len(v)),
                       RB.esc(", ".join(sorted({(r.get("class") or "—").title() for r in v})) or "—")]
                      for smp, v in per.items()],
                     "AMR gene families found in the reads of each sample."))
        if len(S) > 1:
            genes, label = {}, {}
            for r in amr:
                k = r.get("family_id") or r["gene"]
                genes.setdefault(k, r.get("class", ""))
                label.setdefault(k, r["gene"])
            order = sorted(genes, key=lambda g: (genes[g], label[g]))[:60]
            val = {(r["sample"], r.get("family_id") or r["gene"]): float(r.get("rpkm") or 0) for r in amr}
            vmax = max(val.values()) if val else 0
            body = []
            for g in order:
                cells = []
                for smp in S:
                    v = val.get((smp, g), 0)
                    cells.append(f'<td{_shade(v, vmax)}>{_f(v, 1) if v else ""}</td>')
                body.append(f"<tr><td><i>{RB.esc(label[g])}</i></td><td>{RB.esc(genes[g].title())}</td>"
                            + "".join(cells) + "</tr>")
            head = "".join(f"<th>{RB.esc(smp)}</th>" for smp in S)
            out.append(f'<p class="tcap"><b>Table {T.n()}.</b> AMR genes by sample, RPKM (darker = more). '
                       "Empty cells: not found." + (" First 60 gene families by drug class; all in "
                                                    "amr_rpkm.csv." if len(genes) > 60 else "") + "</p>"
                       f"<table><thead><tr><th>Gene</th><th>Class</th>{head}</tr></thead>"
                       f"<tbody>{''.join(body)}</tbody></table>")
        det = sorted(amr, key=lambda r: (S.index(r["sample"]) if r["sample"] in S else 0,
                                         -float(r.get("rpkm") or 0)))
        lim = 80
        out.append(T(["Sample", "Gene", "Class", "Reads", "Gene covered", "Identity", "RPKM", "Likely host"],
                     [[RB.esc(r["sample"]), f"<i>{RB.esc(r['gene'])}</i>"
                       + (f' <span class="dim">({RB.esc(r["allele"])})</span>'
                          if r.get("allele") and r["allele"] != r["gene"] else ""),
                       RB.esc((r.get("subclass") or r.get("class") or "").title()),
                       RB.fmt(r.get("reads")), f"{_f(r.get('breadth'), 0)}%", f"{_f(r.get('identity'), 1)}%",
                       _f(r.get("rpkm"), 1),
                       (f"<i>{RB.esc(_short(r['host_taxon']))}</i> {_f(r.get('host_share'), 0)}%"
                        if r.get("host_taxon") else '<span class="dim">unresolved</span>')]
                      for r in det[:lim]],
                     "AMR genes in the reads (allele in brackets when one matched best). <i>Gene "
                     "covered</i>: share of the gene's length covered by reads; <i>likely host</i>: the "
                     "organism Kraken 2 assigns to most of the reads carrying the gene, with the share of "
                     "those reads." + (f" First {lim} rows; all in amr_genes.csv."
                                                           if len(det) > lim else "")))
    if other:
        lim = 40
        o = sorted(other, key=lambda r: (S.index(r["sample"]) if r["sample"] in S else 0,
                                         r.get("type", ""), -float(r.get("rpkm") or 0)))
        out.append(T(["Sample", "Gene", "Type", "Product", "Reads", "RPKM", "Likely host"],
                     [[RB.esc(r["sample"]), f"<i>{RB.esc(r['gene'])}</i>", RB.esc(r.get("type", "").title()),
                       f'<span class="dim">{RB.esc(r.get("product", ""))}</span>', RB.fmt(r.get("reads")),
                       _f(r.get("rpkm"), 1),
                       f"<i>{RB.esc(_short(r['host_taxon']))}</i>" if r.get("host_taxon") else
                       '<span class="dim">unresolved</span>'] for r in o[:lim]],
                     "Stress-response and virulence genes from the same catalogue (AMRFinderPlus "
                     "&ldquo;plus&rdquo; genes)." + (f" First {lim} rows; all in amr_genes.csv."
                                                    if len(o) > lim else "")))
    out.append('<div class="note acc"><p style="margin:0"><b>Reading the resistome.</b> Genes are found from '
               "reads, without assembly: a gene is reported when its reads cover enough of it, so genes of "
               "rare organisms can be missed, and resistance caused by point mutations (for example "
               "<i>gyrA</i>, <i>parC</i>, 23S rRNA) is not detected. Several genes are intrinsic to their "
               "species and expected whenever it is present (for example <i>blaEC</i> in <i>E. coli</i>, "
               "<i>mexAB</i>/<i>mexEF</i> efflux pumps in <i>P. aeruginosa</i>, <i>fosX</i> in "
               "<i>Listeria</i>, <i>lsa(A)</i> in <i>E. faecalis</i>). The likely host is reliable for "
               "Nanopore reads, which carry kilobases of flanking DNA; for short reads and for mobile genes "
               "shared by many species it is tentative. Gene presence is not a phenotype: confirm "
               "resistance by culture and susceptibility testing.</p></div>")
    return meth, out


def _mags(run, sm, S, T) -> tuple[list[str], list[str], list[str]]:
    """(method sentences, result HTML, reference keys) for MAGs."""
    info = sm.get("mags") or {}
    if not info:
        return [], [], []
    if info.get("error"):
        return [], ['<h3 class="pbreak">Metagenome-assembled genomes</h3>',
                    f"<p>The MAG step did not finish ({RB.esc(info['error'])}); the profile above is "
                    "unaffected. See logs/mags_traceback.txt.</p>"], []
    tv = info.get("tools") or {}
    v = lambda k: f" {tv[k]}" if tv.get(k) else ""
    ill = info.get("assembler", "").startswith("MEGAHIT")
    refs = ["megahit" if ill else "flye"]
    asm = (f"MEGAHIT{v('megahit')}" if ill else f"Flye{v('flye')} in metagenome mode (--meta)")
    meth = [f"Metagenome-assembled genomes (MAGs) were built per sample: the QC-passed reads were assembled "
            f"with {asm}, mapped back to the contigs with minimap2 to measure contig depth, and contigs of at "
            f"least {info.get('min_contig', 1500):,} bp were binned with MetaBAT 2"
            + (f" ({tv['metabat2']})" if tv.get("metabat2") else "") + "."]
    if info.get("binner"):
        refs.append("metabat2")
    if info.get("checkm2"):
        meth.append(f"Bin completeness and contamination were estimated with CheckM2{v('checkm2')}; bins "
                    "are graded after MIMAG: high quality &gt;90% complete and &lt;5% contaminated, medium "
                    "&ge;50% and &lt;10%, otherwise low.")
        refs += ["checkm2", "mimag"]
    meth.append("Each bin was named by Kraken 2 on its contigs (the taxon holding most of the bin's length)"
                + (", and matched to GTDB species representatives with sourmash" if info.get("gtdb") else "")
                + "; AMR genes in the contigs were called with AMRFinderPlus and tied to their bin.")
    if info.get("gtdb"):
        refs.append("gtdb")
    out = ['<h3 class="pbreak">Metagenome-assembled genomes</h3>']
    rows = []
    for smp in S:
        a = (info.get("samples") or {}).get(smp)
        if not a:
            continue
        x = a.get("assembly") or {}
        rows.append([RB.esc(smp), RB.fmt(x.get("contigs")), _f((x.get("bp") or 0) / 1e6, 1) if x else "—",
                     RB.fmt(x.get("n50")), RB.fmt(x.get("longest")), str(a.get("bins", 0)),
                     RB.esc("; ".join(a.get("notes") or [])) or ""])
    if rows:
        out.append(T(["Sample", "Contigs", "Assembly (Mb)", "N50 (bp)", "Longest (bp)", "Bins", "Notes"],
                     rows, "Assembly per sample."))
    bins = _csv(run.path / "mags_summary.csv")
    if bins:
        out.append(T(["Bin", "Quality", "Complete", "Contam.", "Size (Mb)", "Contigs", "N50", "Depth",
                      "Taxon (Kraken 2)"] + (["GTDB (ANI)"] if info.get("gtdb") else []) + ["AMR"],
                     [[RB.esc(b["bin"]), RB.esc(b.get("quality", "")),
                       f"{_f(b.get('completeness'))}%" if b.get("completeness") not in ("", None) else "—",
                       f"{_f(b.get('contamination'))}%" if b.get("contamination") not in ("", None) else "—",
                       _f(float(b.get("size_bp") or 0) / 1e6, 2), RB.fmt(b.get("contigs")), RB.fmt(b.get("n50")),
                       f"{_f(b.get('mean_depth'))}×", f"<i>{RB.esc(b.get('kraken2_taxon') or '—')}</i>"]
                      + ([f"<i>{RB.esc(b.get('gtdb_species') or '—')}</i>"
                          + (f" ({_f(b.get('gtdb_ani'))}%)" if b.get("gtdb_ani") not in ("", None, "0.0") else "")]
                         if info.get("gtdb") else [])
                      + [str(len([g for g in (b.get("amr_genes") or "").split(";") if g])) or "0"]
                      for b in bins[:80]],
                     "Genome bins, best first. <i>Depth</i>: mean read depth of the bin's contigs. Genomes in "
                     "mags/&lt;sample&gt;/, all columns in mags_summary.csv, AMR genes per contig in "
                     "mags_amr.csv. <i>AMR</i>: number of AMR genes in the bin." + (" First 80 bins."
                                                                                   if len(bins) > 80 else "")))
        withamr = [b for b in bins if b.get("amr_genes")]
        if withamr:
            out.append(T(["Bin", "Taxon", "AMR genes"],
                         [[RB.esc(b["bin"]), f"<i>{RB.esc(b.get('kraken2_taxon') or '—')}</i>",
                           "<i>" + RB.esc(", ".join(b["amr_genes"].split(";"))) + "</i>"] for b in withamr[:80]],
                         "AMR genes found by AMRFinderPlus in the contigs of each bin."))
    elif rows:
        out.append("<p>No bins were formed. Binning needs depth: roughly 10× or more of an organism's genome; "
                   "low-depth or very diverse samples often assemble into short contigs only.</p>")
    out.append('<div class="note acc"><p style="margin:0"><b>Reading MAGs.</b> A bin is a draft genome '
               "pieced together from a mixture; high-quality MAGs are reliable for species identity and gene "
               "content, low-quality ones are fragments. Plasmids and other mobile elements often stay "
               "unbinned, so AMR genes on them appear in mags_amr.csv without a bin."
               + ("" if ill else " Nanopore assemblies here are not polished: expect small indels in genes.")
               + "</p></div>")
    return meth, out, refs


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

    def _tn():
        tab[0] += 1
        return tab[0]
    T.n = _tn

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
               "left out of the community tables"
               + (f"; species below {sm.get('min_abundance_pct')}% of the reads in every sample "
                  f"({sm.get('species_below_min_abundance', 0)} here) were left out as well — mostly reads "
                  "Bracken spreads onto close relatives of abundant species"
                  if sm.get("min_abundance_pct") else "")
               + ". Diversity, composition and the statistics below use the Bracken species table.</p>")
    if platform == "ont" and blen and max(int(x) for x in blen.values()) < 1000:
        out.append("<p class=\"dim\">Bracken's k-mer distributions are built for a fixed read length; for "
                   "Nanopore reads the longest available one was used, which is shorter than the reads. "
                   "Species estimates remain usable; genus-level figures are the more robust ones.</p>")
    # method text now; the result tables later, so that table numbers follow page order
    def _noT(*a, **k):
        return ""
    _noT.n = lambda: 0
    amr_meth, _ = _resistome(run, sm, S, _noT, F)
    mag_meth, _, mag_refs = _mags(run, sm, S, _noT)
    for m in (amr_meth, mag_meth):
        if m:
            out.append("<p>" + " ".join(m) + "</p>")
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

    out += _resistome(run, sm, S, T, F)[1]
    out += _mags(run, sm, S, T)[1]
    out += RB._controls_section(run, an, T)
    out += RB._maaslin_section(run, T)

    # Close relatives around a dominant species: shared k-mers, not extra organisms
    spill = []
    pct = an.rank_pct("Species") if "Species" in run.rank_cols else {}
    for smp in S:
        ranked = sorted(((t, v[smp]) for t, v in pct.items() if t != "Unclassified"), key=lambda x: -x[1])
        if ranked and ranked[0][1] >= 50:
            top, share = ranked[0]
            g = top.split()[0]
            rel = [t for t, v in ranked[1:] if t.split()[0] == g and 0 < v < 1]
            if len(rel) >= 5:
                spill.append(f"{esc(smp)}: <i>{esc(top)}</i> {share:.1f}% with {len(rel)} other <i>{esc(g)}</i> "
                             "species below 1%")
    out.append('<h2 class="pbreak">Notes on reading this report</h2>')
    extra = ([f"<b>Close relatives of a dominant species.</b> {'; '.join(spill)}. Reads of one organism "
              "that match its relatives' genomes as well are spread over them by Kraken 2 and Bracken; such "
              "small same-genus entries are expected even from a pure culture and are not evidence of "
              "additional organisms. Raise the minimum abundance (pipeline step 3) to hide them."] if spill else [])
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in extra + [
        "Shotgun read counts reflect DNA, not cells: larger genomes contribute more reads per cell, so "
        "relative abundance is biased towards organisms with big genomes.",
        "Kraken 2 can only name organisms whose genomes, or close relatives, are in the index. Species "
        "near the detection limit (a few hundred reads) should be confirmed before being reported.",
        "Full tables: asv_table.csv (species × sample, estimated reads; first column = NCBI taxid), "
        "taxonomy.csv (lineage per taxid), read_tracking.csv, and per-sample Kraken 2 / Bracken reports in "
        "shotgun/&lt;sample&gt;/."
        + (" Resistome: amr_genes.csv (one row per gene and sample) and amr_rpkm.csv (gene × sample)."
           if amr_meth else "")
        + (" MAGs: mags_summary.csv, mags_amr.csv and the genomes in mags/." if mag_meth else ""),
    ]) + "</ul>")
    out.append("<h2>References</h2>")
    refs = REFS[:2] + ([REFS[2]] if platform == "illumina" else [REFS[4]]) + \
        ([REFS[3]] if host != "none" else [])
    if amr_meth:
        refs += REFS_AMR
    refs += [REFS_MAG[k] for k in mag_refs]
    if (run.path / "r_tables" / "maaslin3_results.csv").exists():
        refs.append(RB.REFERENCES[25])
    out.append('<ul class="refs">' + "".join(f"<li>{esc(a)} <i>{esc(t)}</i> {esc(j)}</li>"
                                             for a, t, j in refs) + "</ul>")
    css = RB._css(app, company, f"{run.label} · shotgun") + RB._CSS_FULL.replace("@@TITLE@@", esc(rep_title))
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{esc(rep_title)}</title>'
            f"<style>{css}</style></head><body>" + "".join(out) + "</body></html>")
