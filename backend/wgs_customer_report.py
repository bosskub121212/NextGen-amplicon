"""wgs_customer_report.py — the lab-style ONT-WGS report of one sample, with its package.

    build_sample(result_dir, sample, out_dir, formats=("html", "pdf", "docx")) -> manifest

Builds the eight-folder package (wgs_package.build_package) and writes into its root
Report_<sample>.html / .pdf / .docx, laid out like the lab's "Whole genome sequencing
report": no cover, "Sample: X", sections 1–8 in folder order, each ending with
"Results <folder>". Figures are also saved as PNG in their folders. Every sentence names the
tool that really ran here (Flye + Medaka, CheckM2, sourmash/GTDB, ANIb/ANIm/TETRA, Bakta,
AMRFinderPlus, ...) with this run's numbers; a step that did not run says so.
"""
from __future__ import annotations

import csv
import html
import json
import math
from pathlib import Path

import wgs_package as WP

esc = html.escape
ACC, INK, MUT = "#3c9fae", "#1f2937", "#6b7280"

REFS = {
    "seqkit": "Shen W, Le S, Li Y, Hu F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE 11:e0163962.",
    "flye": "Kolmogorov M, Yuan J, Lin Y, Pevzner PA. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology 37:540–546.",
    "medaka": "Oxford Nanopore Technologies. medaka: sequence correction provided by ONT Research (software).",
    "checkm2": "Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nature Methods 20:1203–1212.",
    "sourmash": "Brown CT, Irber L. (2016). sourmash: a library for MinHash sketching of DNA. Journal of Open Source Software 1:27.",
    "gtdb": "Parks DH, Chuvochina M, Rinke C, et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Research 50:D785–D794.",
    "skani": "Shaw J, Yu YW. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods 20:1661–1665.",
    "jspecies": "Richter M, Rosselló-Móra R. (2009). Shifting the genomic gold standard for the prokaryotic species definition. PNAS 106:19126–19131.",
    "fastme": "Lefort V, Desper R, Gascuel O. (2015). FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program. Molecular Biology and Evolution 32:2798–2800.",
    "mafft": "Katoh K, Standley DM. (2013). MAFFT multiple sequence alignment software version 7. Molecular Biology and Evolution 30:772–780.",
    "iqtree": "Minh BQ, Schmidt HA, Chernomor O, et al. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution 37:1530–1534.",
    "mlst": "Jolley KA, Bray JE, Maiden MCJ. (2018). Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Research 3:124.",
    "bakta": "Schwengers O, Jelonek L, Dieckmann MA, et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics 7:000685.",
    "amrfinder": "Feldgarden M, Brover V, Gonzalez-Escalona N, et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific Reports 11:12728.",
    "vfdb": "Chen L, Yang J, Yu J, et al. (2005). VFDB: a reference database for bacterial virulence factors. Nucleic Acids Research 33:D325–D328.",
    "plasmidfinder": "Carattoli A, Zankari E, García-Fernández A, et al. (2014). In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrobial Agents and Chemotherapy 58:3895–3903.",
    "mobsuite": "Robertson J, Nash JHE. (2018). MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. Microbial Genomics 4:e000206.",
    "minimap2": "Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094–3100.",
    "samtools": "Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience 10:giab008.",
    "kraken2": "Wood DE, Lu J, Langmead B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology 20:257.",
    "bracken": "Lu J, Breitwieser FP, Thielen P, Salzberg SL. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104.",
    "centrifuge": "Kim D, Song L, Breitwieser FP, Salzberg SL. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome Research 26:1721–1729.",
    "centrifuger": "Song L, Langmead B. (2024). Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. Genome Biology 25:106.",
    "gottcha": "Freitas TAK, Li PE, Scholz MB, Chain PSG. (2015). Accurate read-based metagenome characterization using a hierarchical suite of unique signatures. Nucleic Acids Research 43:e69.",
    "tygs": "Meier-Kolthoff JP, Göker M. (2019). TYGS is an automated high-throughput platform for state-of-the-art genome-based taxonomy. Nature Communications 10:2182.",
    "gbdp": "Meier-Kolthoff JP, Auch AF, Klenk HP, Göker M. (2013). Genome sequence-based species delimitation with confidence intervals and improved distance functions. BMC Bioinformatics 14:60.",
}


# ── small helpers ─────────────────────────────────────────────────────────────
def _f(x, d=2):
    try:
        return f"{float(x):,.{d}f}"
    except (TypeError, ValueError):
        return "—"


def _i(x):
    try:
        return f"{int(float(x)):,}"
    except (TypeError, ValueError):
        return "—"


def ital(sp: str) -> str:
    return f"<i>{esc(sp)}</i>" if sp else "—"


def _tbl(head, rows, cls="t") -> str:
    th = "".join(f"<th>{esc(str(h))}</th>" for h in head)
    tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>'


def _read_tsv(p: Path) -> list[dict]:
    try:
        with open(p, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return []


def results(folder: str) -> str:
    return f'<p class="res">Results {esc(folder)}</p>'


# ── figures (SVG; also saved as PNG in the package) ───────────────────────────
def svg_hist(bins, vals, log=False, title="", w=640, h=200) -> str:
    if not bins or not any(vals):
        return ""
    edges = list(bins) + [bins[-1] * 2 if bins[-1] else 1]
    pl, pr, pt, pb = 54, 14, 26, 40
    W, H = w - pl - pr, h - pt - pb
    lo = math.log10(max(edges[1], 1) / 2) if log else 0
    hi = math.log10(edges[-1]) if log else len(vals)

    def X(i, e):
        if not log:
            return pl + W * i / len(vals)
        return pl + W * (math.log10(max(e, 1) if e else 10 ** lo) - lo) / (hi - lo)
    mx = max(vals) or 1
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
         f'font-family="DejaVu Sans" font-size="9"><rect width="{w}" height="{h}" fill="#fff"/>',
         f'<text x="{w / 2}" y="15" text-anchor="middle" font-size="11" fill="{INK}">{esc(title)}</text>']
    for i, v in enumerate(vals):
        x0, x1 = X(i, edges[i]), X(i + 1, edges[i + 1])
        bh = H * v / mx
        o.append(f'<rect x="{x0 + 0.5:.1f}" y="{pt + H - bh:.1f}" width="{max(0.8, x1 - x0 - 1):.1f}" '
                 f'height="{bh:.1f}" fill="{ACC}"/>')
        lab = "0" if edges[i] == 0 else (f"{edges[i] // 1000}k" if edges[i] >= 1000 else str(edges[i]))
        o.append(f'<text x="{x0:.1f}" y="{pt + H + 12}" text-anchor="middle" fill="{MUT}" font-size="8">{lab}</text>')
    for k in range(5):
        y = pt + H - H * k / 4
        o.append(f'<line x1="{pl}" y1="{y:.1f}" x2="{pl + W}" y2="{y:.1f}" stroke="#e5e7eb"/>'
                 f'<text x="{pl - 4}" y="{y + 3:.1f}" text-anchor="end" fill="{MUT}" font-size="8">'
                 f'{_human(mx * k / 4)}</text>')
    o.append(f'<text x="{pl + W / 2}" y="{h - 6}" text-anchor="middle" fill="{INK}">Read length (bp)'
             f'{" — log scale" if log else ""}</text>'
             f'<text x="12" y="{pt + H / 2}" transform="rotate(-90 12 {pt + H / 2})" text-anchor="middle" '
             f'fill="{INK}">Number of bases</text></svg>')
    return "".join(o)


def _human(v):
    for u, d in (("G", 1e9), ("M", 1e6), ("k", 1e3)):
        if v >= d:
            return f"{v / d:.0f}{u}"
    return f"{v:.0f}"


def svg_bars(labels, vals, title, ylab, logy=False, w=330, h=240) -> str:
    if not vals:
        return ""
    pl, pr, pt, pb = 50, 10, 26, 60
    W, H = w - pl - pr, h - pt - pb
    f = (lambda v: math.log10(v + 1)) if logy else (lambda v: v)
    mx = max(f(v) for v in vals) or 1
    bw = W / len(vals)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
         f'font-family="DejaVu Sans" font-size="9"><rect width="{w}" height="{h}" fill="#fff"/>',
         f'<text x="{w / 2}" y="15" text-anchor="middle" font-size="11" font-weight="bold" fill="{INK}">{esc(title)}</text>']
    for i, (lab, v) in enumerate(zip(labels, vals)):
        bh = H * f(v) / mx
        x = pl + i * bw
        o.append(f'<rect x="{x + bw * 0.15:.1f}" y="{pt + H - bh:.1f}" width="{bw * 0.7:.1f}" height="{bh:.1f}" '
                 f'fill="#bdbdbd" stroke="#555" stroke-width=".6"/>'
                 f'<text x="{x + bw / 2:.1f}" y="{pt + H - bh - 3:.1f}" text-anchor="middle" font-size="8">{v:,}</text>'
                 f'<text x="{x + bw / 2:.1f}" y="{pt + H + 10}" text-anchor="end" font-size="8.5" '
                 f'transform="rotate(-45 {x + bw / 2:.1f} {pt + H + 10})">{esc(lab)}</text>')
    o.append(f'<line x1="{pl}" y1="{pt + H}" x2="{pl + W}" y2="{pt + H}" stroke="{INK}"/>'
             f'<text x="12" y="{pt + H / 2}" transform="rotate(-90 12 {pt + H / 2})" text-anchor="middle">'
             f'{esc(ylab)}</text></svg>')
    return "".join(o)


def svg_protein_hist(lengths, w=330, h=240) -> str:
    if not lengths:
        return ""
    step = 100
    nb = min(40, max(lengths) // step + 1)
    counts = [0] * nb
    for L in lengths:
        counts[min(nb - 1, L // step)] += 1
    return svg_bars([f"{i * step}" if i % 5 == 0 else "" for i in range(nb)], counts,
                    "Distribution of protein size", "Count", w=w, h=h).replace(
        'text-anchor="middle" font-size="8">', 'text-anchor="middle" font-size="0">')


def _faa_lengths(p: Path) -> list[int]:
    out, cur = [], 0
    try:
        with open(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        out.append(cur)
                    cur = 0
                else:
                    cur += len(line.strip())
    except OSError:
        return []
    if cur:
        out.append(cur)
    return out


def save_png(svg: str, dest: Path) -> bool:
    if not svg:
        return False
    try:
        import cairosvg
    except ImportError:
        return False
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        cairosvg.svg2png(bytestring=svg.encode(), write_to=str(dest), scale=2)
        return True
    except Exception:
        return False


def _fig(svg: str, cap: str = "") -> str:
    if not svg:
        return ""
    return f'<div class="fig">{svg}' + (f'<p class="cap">{cap}</p>' if cap else "") + "</div>"


# ── sections ──────────────────────────────────────────────────────────────────
def s_qc(S, res, rp, root, used) -> str:
    raw, flt = res.get("reads_raw") or {}, res.get("reads_filtered") or {}
    o = ['<h2 class="first">1.QC</h2>']
    if not raw:
        return "".join(o) + "<p>Read statistics were not recorded for this sample.</p>" + results("1.QcReads")
    used.add("seqkit")
    pc = lambda a, b: f" ({100 * a / b:.2f} %)" if a and b else ""
    o += ["<h3>Before Trimming</h3>",
          f"<p>Reads #: {raw.get('reads', '—')}<br>Total bases: {raw.get('bases', '—')}<br>"
          f"Reads Length: {raw.get('mean_len', '—')}</p>"]
    if flt:
        o += ["<h3>After Trimming</h3>",
              f"<p>Reads #: {flt.get('reads', '—')}{pc(flt.get('reads'), raw.get('reads'))}<br>"
              f"Total bases: {flt.get('bases', '—')}{pc(flt.get('bases'), raw.get('bases'))}<br>"
              f"Mean Reads Length: {flt.get('mean_len', '—')}</p>",
              f"<p class='m'>Reads shorter than {rp.get('min_read_len', 1000)} bp or with a mean quality below "
              f"Q{rp.get('min_read_q', 10)} were removed with seqkit.</p>"]
    for log, name, title in ((False, "read_length_histogram", "Weighted histogram of read lengths"),
                             (True, "read_length_histogram_log", "Weighted histogram of read lengths after log transformation")):
        sv = svg_hist(raw.get("len_bins") or [], raw.get("len_hist_bp") or [], log=log, title=title)
        save_png(sv, root / "1.QcReads" / f"{name}.png")
        o.append(f'<div class="band">{esc(title)}</div>' + _fig(sv))
    return "".join(o) + results("1.QcReads")


def s_assembly(S, res, root, used) -> str:
    o = ['<h2>2.Genome Assembly</h2>']
    fa = root / "2.Assembly" / f"{S}_assembly.fasta"
    if not fa.exists():
        return "".join(o) + "<p>The assembly was not produced for this sample.</p>" + results("2.Assembly")
    st = dict(WP.assembly_stats(fa))
    a = res.get("assembly") or {}
    used.update({"flye"} | ({"medaka"} if a.get("polish") else set()))
    pol = a.get("polish") or ""
    pol_t = (f" and polished with {esc(pol[0].upper() + pol[1:])}" if pol and not pol.lower().startswith("no")
             else "")
    o.append(f"<p>{esc(S)} was assembled using Flye{pol_t}. There were {st['# contigs']} contigs, an estimated "
             f"genome length of {st['Total length']:,} bp, and an average G+C content of {st['GC (%)']:.2f}%. "
             f"The N50 length, which is defined as the shortest sequence length at 50% of the genome is "
             f"{st['N50']:,} bp. The L50 count which is defined as the smallest number of contigs whose length sum "
             f"produces N50 is {st['L50']} (Table 1).</p>")
    o.append(_tbl(["Statistics without reference", f"{esc(S)}_assembly"],
                  [[esc(k), f"{v:,}" if isinstance(v, int) else f"{v:.2f}"] for k, v in WP.assembly_stats(fa)],
                  "t narrow"))
    o.append('<p class="cap">Table 1. Assembly statistics (computed by the app in the way QUAST reports them).</p>')
    return "".join(o) + results("2.Assembly")


def s_completeness(S, res, root, used) -> str:
    o = ['<h2>3. Estimates of genome completeness and contamination</h2>']
    rows = _read_tsv(root / "3.Completeness" / "checkm2_quality_report.tsv")
    if not rows:
        return "".join(o) + "<p>CheckM2 was not run for this sample.</p>" + results("3.Completeness")
    used.add("checkm2")
    head = list(rows[0].keys())
    body = [[esc(f"{S}_assembly" if k == "Name" else r.get(k, "")) for k in head] for r in rows]
    o.append("<p>Genome completeness and contamination were estimated with CheckM2.</p>")
    o.append(_tbl([h.replace("_", " ") for h in head], body, "t tiny"))
    return "".join(o) + results("3.Completeness")


CONCLUSION = {"known": "belongs to known species", "borderline": "borderline (95–96 % ANI)",
              "novel": "potential new species"}


def s_taxonomy(S, res, sd, root, used, tygs) -> str:
    o = ['<h2>4. Genome-based taxonomy</h2>', "<h3>Table 4.1 Genome-based identification of strain</h3>"]
    ts = res.get("typestrain") or {}
    idn = res.get("identification") or {}
    if tygs and tygs.get("identification"):
        used.update({"tygs", "gbdp"})
        o.append(_tbl(["Strain", "Conclusion", "Identification result", "Remark"],
                      [[esc(r.get("strain", "")), esc(r.get("conclusion", "")), ital(r.get("species", "")),
                        esc(r.get("remark", ""))] for r in tygs["identification"]]))
        if tygs.get("dddh"):
            o.append("<p>The following table contains the pairwise dDDH values between the genome and the selected "
                     "type-strain genomes, with their confidence intervals (C.I.) for the three GBDP formulas:</p>")
            o.append(_tbl(["Query", "Subject", "d0", "C.I. d0", "d4", "C.I. d4", "d6", "C.I. d6", "Diff. G+C %"],
                          [[esc(r.get(k, "")) for k in ("query", "subject", "d0", "ci_d0", "d4", "ci_d4", "d6",
                                                        "ci_d6", "gc_diff")] for r in tygs["dddh"]], "t tiny"))
            o.append("<ul class='m'><li>formula d0 (GGDC formula 1): length of all HSPs divided by total genome length</li>"
                     "<li>formula d4 (GGDC formula 2): sum of all identities found in HSPs divided by overall HSP length</li>"
                     "<li>formula d6 (GGDC formula 3): sum of all identities found in HSPs divided by total genome length</li>"
                     "</ul><p class='m'>Note: formula d4 is independent of genome length and is thus robust against the "
                     "use of incomplete draft genomes. Values from the Type (Strain) Genome Server (TYGS, DSMZ).</p>")
        o.append(results("4.Taxonomy/TYGS"))
    else:
        v = ts.get("verdict") or {}
        concl = CONCLUSION.get(v.get("level"), idn.get("call", "") or "undetermined")
        sp = v.get("species") or idn.get("species", "")
        meth = str(idn.get("method", "") or "")
        rem = (f"ANI {_f(idn.get('ani'))} % ({meth or 'ANI'}) to {idn.get('reference', '')} "
               f"({idn.get('accession', '')})" if idn.get("ani") else meth)
        jsp = any(t.get(k) is not None for t in (ts.get("type_strains") or []) for k in ("anib", "anim", "tetra"))
        if "skani" in meth.lower() or (ts.get("type_strains") and not jsp):
            used.add("skani")
        if "sourmash" in meth.lower():
            used.update({"sourmash", "gtdb"})
        if jsp:
            used.add("jspecies")
        o.append(_tbl(["Strain", "Conclusion", "Identification result", "Remark"],
                      [[esc(f"{S}_assembly"), esc(concl), ital(sp), esc(rem)]]))
        basis = ("average nucleotide identity (ANIb, ANIm, TETRA as in JSpecies) against type-strain genomes "
                 "(NCBI assembly from type material)" if jsp else
                 f"{meth}" if meth else "average nucleotide identity")
        o.append(f"<p class='m'>Identification by {esc(basis)}; ANI ≥ 95–96 % = same species. "
                 "dDDH (TYGS) not imported for this sample — upload the assembly to TYGS and import the result to "
                 "add the dDDH comparison and the TYGS trees.</p>")
    # trees
    o.append("<h3>Phylogenetic tree</h3>")
    trees = (tygs or {}).get("trees") or {}
    from wgs_report import svg_tree
    shown = False
    for key, lab in (("16S", "16S rRNA gene tree"), ("genome", "Genome tree")):
        nwk = trees.get(key)
        if nwk:
            sv = svg_tree(nwk, S, w=640)
            save_png(sv, root / "4.Taxonomy" / "TYGS" / f"{key}_tree.png")
            o.append(_fig(sv, f"{lab} from TYGS (FastME from GBDP distances)."))
            shown = True
    if not shown:
        for key, mkey, lab in (("rrna_tree", "rrna_tree_method", "16S rRNA gene tree"),
                               ("genome_tree", "genome_tree_method", "Genome tree")):
            if ts.get(key):
                used.update({"fastme"} if key == "genome_tree" else {"mafft", "iqtree"})
                sv = svg_tree(ts[key], S, w=640)
                save_png(sv, root / "4.Taxonomy" / "Phylogeny" / f"{key.replace('rrna', '16S')}.png")
                o.append(_fig(sv, f"{lab} of the isolate and its closest type strains — {esc(ts.get(mkey, ''))}."))
                shown = True
    if not shown:
        o.append("<p>Phylogenetic trees were not built for this sample.</p>")
    o.append(results("4.Taxonomy/TYGS" if trees else "4.Taxonomy/Phylogeny"))
    # read-based
    o.append("<h3>Read-based Taxonomy Classification</h3>")
    rt = sd / "readtax"
    if (rt / "read_taxonomy.tsv").exists():
        import readtax
        info = readtax.load_info(rt)
        used.update(info.get("tools_ref", []))
        sv = readtax.svg_heatmap(rt / "read_taxonomy.tsv")
        save_png(sv, root / "4.Taxonomy" / "Read-based Taxonomy" / "heatmap_read_taxonomy.png")
        o.append(f"<p>{esc(readtax.method_sentence(info))}</p>" + _fig(sv))
    else:
        o.append("<p>Read-based taxonomy was not run for this sample.</p>")
    o.append(results("4.Taxonomy/Read-based Taxonomy"))
    # contig-based
    o.append("<h3>Contig-based Community Profiling</h3>")
    hits = res.get("sourmash") or []
    if hits:
        used.update({"sourmash", "gtdb"})
        ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
        agg = {r: {} for r in ranks}
        for h in hits:
            parts = (h.get("gtdb_lineage") or "").split(";")
            for r, p in zip(ranks, parts):
                name = p.split("__", 1)[-1] or "unclassified"
                agg[r][name] = agg[r].get(name, 0) + float(h.get("f_query") or 0)
        rows = []
        for r in ranks:
            top = sorted(agg[r].items(), key=lambda kv: -kv[1])[:5]
            rows.append([r] + [f"{ital(n) if r in ('genus', 'species') else esc(n)} ({v:.1f} %)" for n, v in top]
                        + ["NA"] * (5 - len(top)))
        o.append("<p>Composition of the assembly by sourmash gather against GTDB representative genomes: share of "
                 "the assembly's k-mers explained by each genome, summed per rank.</p>")
        o.append(_tbl(["Level", "Top1", "Top2", "Top3", "Top4", "Top5"], rows, "t tiny"))
        sp = sorted(agg["species"].items(), key=lambda kv: -kv[1])[:10]
        sv = svg_bars([n for n, _ in sp], [round(v, 1) for _, v in sp], f"Rank: Species ({len(agg['species'])})",
                      "% of assembly k-mers", w=520, h=200)
        save_png(sv, root / "4.Taxonomy" / "Contig-based Taxonomy" / "species_composition.png")
        o.append(_fig(sv))
    else:
        o.append("<p>Assembly composition (sourmash / GTDB) was not computed for this sample.</p>")
    o.append(results("4.Taxonomy/Contig-based Taxonomy"))
    # MLST
    o.append("<h3>MLST</h3>")
    ml = res.get("mlst") or {}
    if ml.get("scheme") and ml.get("scheme") != "-":
        used.add("mlst")
        o.append(_tbl(["Scheme", "Sequence type (ST)", "Alleles"],
                      [[esc(ml.get("scheme", "")), esc(str(ml.get("st", ""))), esc(ml.get("alleles", ""))]]))
        o.append("<p class='m'>mlst against PubMLST schemes; ST '-' = new or incomplete allele combination.</p>")
    else:
        o.append("<p>No MLST scheme matched this genome.</p>")
    return "".join(o) + results("4.Taxonomy/MLST")


def s_ani(S, res, root, used) -> str:
    o = ['<h2>5. ANI analysis</h2>']
    ts = (res.get("typestrain") or {}).get("type_strains") or []
    if not ts:
        return "".join(o) + "<p>The type-strain ANI comparison was not run for this sample.</p>" + results("5.ANI")
    cols = [("anib", "ANIb [%]", 2), ("anib_cov_q", "Aligned [%]", 1), ("anim", "ANIm [%]", 2),
            ("tetra", "TETRA", 4), ("skani_ani", "skani ANI [%]", 2)]
    cols = [c for c in cols if any(t.get(c[0]) is not None for t in ts)]      # only what was computed
    have = {c[0] for c in cols}
    if have & {"anib", "anim", "tetra"}:
        used.add("jspecies")
    if "skani_ani" in have:
        used.add("skani")
    rows = [[f"{ital(t.get('organism', ''))} {esc(str(t.get('strain') or ''))} [T]", esc(t.get("accession", ""))]
            + [_f(t.get(k), nd) for k, _, nd in cols] for t in ts]
    o.append(_tbl(["Genome", "Accession"] + [lab for _, lab, _ in cols], rows))
    notes = []
    if "anib" in have:
        notes.append("ANIb: mean identity of 1020-nt fragments (BLAST+), both directions; aligned = share of the "
                     "isolate genome in accepted fragments")
    if "anim" in have:
        notes.append("ANIm: nucmer")
    if "tetra" in have:
        notes.append("TETRA: tetranucleotide signature correlation (species ≥ 0.999)")
    if "skani_ani" in have:
        notes.append("skani: sparse-chaining ANI")
    o.append("<p class='m'>" + esc(". ".join(notes)) + ". Species threshold: ANI ≥ 95–96 %.</p>")
    return "".join(o) + results("5.ANI")


def s_annotation(S, res, sd, root, used) -> str:
    o = ['<h2>6. Genome annotation</h2>']
    bj = sd / "bakta" / f"{S}.json"
    summ = (res.get("bakta") or {}).get("summary") or {}
    if not summ and not bj.exists():
        return "".join(o) + "<p>Genome annotation was not run for this sample.</p>" + results("6.Annotation")
    used.add("bakta")
    ver, dbv, tt = "", "", ""
    try:
        d = json.loads(bj.read_text())
        ver = (d.get("version") or {}).get("bakta", "")
        db = (d.get("version") or {}).get("db") or {}
        dbv = f"database v{db.get('version', '')} ({db.get('type', '')})" if db else ""
        tt = str((d.get("genome") or {}).get("translation_table", "") or "")
    except Exception:
        pass
    if not tt:
        rows = _read_tsv(root / "3.Completeness" / "checkm2_quality_report.tsv")
        tt = rows[0].get("Translation_Table_Used", "") if rows else ""
    o.append(f"<p>{esc(S)} genome was annotated using Bakta{(' v' + esc(ver)) if ver else ''}"
             f"{(', ' + esc(dbv)) if dbv else ''}. This genome is in the superkingdom Bacteria"
             f"{(' and was annotated using genetic code ' + esc(tt)) if tt else ''}.</p>")
    feats = [(k, summ.get(key)) for k, key in (("CDS", "CDSs"), ("rRNA", "rRNAs"), ("tRNA", "tRNAs"),
                                               ("tmRNA", "tmRNAs"), ("ncRNA", "ncRNAs"))]
    feats = [(k, int(float(v))) for k, v in feats if v not in (None, "")]
    sv1 = svg_bars([k for k, _ in feats], [v for _, v in feats], "Feature count", "Count (log10)", logy=True)
    sv2 = svg_protein_hist(_faa_lengths(sd / "bakta" / f"{S}.faa"))
    for sv, n in ((sv1, "feature_count.png"), (sv2, "protein_size_distribution.png")):
        save_png(sv, root / "6.Annotation" / n)
    o.append(f'<div class="fig2">{sv1}{sv2}</div>')
    keep = ["CDSs", "hypotheticals", "pseudogenes", "rRNAs", "tRNAs", "tmRNAs", "ncRNAs", "CRISPR arrays", "oriCs"]
    o.append(_tbl(["Feature", "Count"], [[esc(k), esc(str(summ[k]))] for k in keep if k in summ], "t narrow"))
    return "".join(o) + results("6.Annotation")


def s_amr(S, res, root, used) -> str:
    o = ['<h2>7.Antimicrobial Resistance Genes</h2>']
    am = res.get("amrfinder")
    if am is None and not (root / "7.AMR" / "amrfinderplus.tsv").exists():
        o.append("<p>AMRFinderPlus was not run for this sample.</p>")
    else:
        used.add("amrfinder")
        am = am or []
        o.append("<p>A summary of the AMR, stress-response and virulence genes found in this genome by AMRFinderPlus "
                 "(--plus) and the corresponding resistance class is provided in the table.</p>")
        if am:
            o.append(_tbl(["Gene", "Product", "Type", "Class", "Subclass", "Method", "Contig", "Position"],
                          [[esc(g.get("gene", "")), esc(g.get("name", "")), esc(g.get("type", "")),
                            esc(g.get("class", "")), esc(g.get("subclass", "")), esc(g.get("method", "")),
                            esc(g.get("contig", "")), f"{esc(str(g.get('start', '')))}–{esc(str(g.get('stop', '')))}"]
                           for g in am], "t tiny"))
        else:
            o.append("<p>No AMR genes were found.</p>")
    vf = res.get("vfdb")
    if vf is not None:
        used.add("vfdb")
        o.append("<h3>Virulence factors (VFDB)</h3>")
        o.append(_tbl(["Gene", "Product", "Coverage %", "Identity %", "Contig"],
                      [[esc(g.get("gene", "")), esc(g.get("product", "")), esc(g.get("coverage", "")),
                        esc(g.get("identity", "")), esc(g.get("contig", ""))] for g in vf], "t tiny")
                 if vf else "<p>No virulence factor genes were found.</p>")
    pf, mob = res.get("plasmidfinder"), res.get("mobsuite") or {}
    if pf is not None or mob:
        used.update(({"plasmidfinder"} if pf is not None else set()) | ({"mobsuite"} if mob else set()))
        o.append("<h3>Plasmids</h3>")
        pls = mob.get("plasmids") or []
        o.append(_tbl(["Plasmid", "Contigs", "Size (bp)", "Replicons", "Mobility"],
                      [[esc(str(p.get("id", p.get("cluster", "")))), esc(", ".join(p.get("contigs", []))
                                                                         if isinstance(p.get("contigs"), list)
                                                                         else str(p.get("contigs", ""))),
                        _i(p.get("size")), esc(str(p.get("replicons", ""))), esc(str(p.get("mobility", "")))]
                       for p in pls], "t tiny") if pls else "<p>No plasmid was reconstructed (MOB-suite).</p>")
    return "".join(o) + results("7.AMR")


def s_mapping(S, res, root, used) -> str:
    o = ['<h2>8. Reads Mapped to Reference</h2>']
    m = res.get("mapping") or {}
    if not m:
        return "".join(o) + ("<p>Mapping to a reference genome was not performed for this sample (no closest "
                             "reference genome, or the step was skipped).</p>") + results("8.readsMappingToRef")
    used.update({"minimap2", "samtools"})
    cols = [("Reference", f"{ital(m.get('reference', ''))}"), ("Accession", esc(m.get("accession", ""))),
            ("Ref Len", _i(m.get("ref_length"))), ("Mapped Reads %", _f(m.get("reads_mapped_pct"))),
            ("Avg Fold", f"{_f(m.get('mean_depth'))}X" if m.get("mean_depth") is not None else None),
            ("Base Coverage ≥1x %", _f(m.get("breadth_1x"))), ("≥10x %", _f(m.get("breadth_10x"))),
            ("Absent regions", _i(m.get("absent_regions")) if m.get("absent_regions") is not None else None),
            ("SNVs", _i(m.get("snps")) if m.get("snps") is not None else None),
            ("INDELs", _i(m.get("indels")) if m.get("indels") is not None else None)]
    cols = [(h, v) for h, v in cols if v not in (None, "—")]
    o.append(_tbl([h for h, _ in cols], [[v for _, v in cols]], "t tiny"))
    from wgs_report import svg_coverage
    sv = svg_coverage(m)
    if sv:
        save_png(sv, root / "8.readsMappingToRef" / "coverage_plot.png")
        o.append(_fig(sv, f"Read depth along the reference (minimap2 map-ont, {_i(m.get('bin_size'))}-bp windows)."))
    return "".join(o) + results("8.readsMappingToRef")


CSS = """
@page { size: A4; margin: 22mm 20mm 18mm 20mm;
  @top-right { content: "Page | " counter(page); font-size: 9pt; color: #555; font-family: Carlito, Calibri, "DejaVu Sans", sans-serif; } }
body { font-family: Carlito, Calibri, "DejaVu Sans", sans-serif; font-size: 11pt; line-height: 1.45; color: #111; }
h1 { font-size: 20pt; text-align: center; font-weight: normal; margin: 0 0 8mm; }
.sample { font-size: 12pt; font-weight: bold; margin-bottom: 5mm; }
h2 { font-size: 12.5pt; margin: 0 0 3mm; page-break-before: always; } h2.first { page-break-before: avoid; }
h3 { font-size: 11.5pt; margin: 4mm 0 2mm; }
p { margin: 0 0 3mm; text-align: justify; } p.m { font-size: 9.5pt; color: #333; }
p.res { font-weight: bold; font-size: 12pt; margin-top: 4mm; page-break-before: avoid; }
.band { background: #3cc8c8; color: #fff; font-size: 8pt; padding: 2.5mm 3mm; margin-top: 3mm; }
.fig { text-align: center; margin: 2mm 0; } .fig2 { display: flex; gap: 3mm; justify-content: center; }
p.cap { text-align: center; font-size: 9pt; color: #444; }
table.t { border-collapse: collapse; width: 100%; font-size: 9.5pt; margin: 2mm 0; }
table.t th, table.t td { border: 1px solid #c8c8c8; padding: 1.5mm 2mm; text-align: left; vertical-align: top; }
table.t th { background: #f2f2f2; } table.tiny { font-size: 7.6pt; } table.narrow { width: auto; min-width: 90mm; }
"""


def build_html(result_dir: Path, sample: str, root: Path, tygs: dict | None = None) -> str:
    sd = Path(result_dir) / "samples" / sample
    try:
        res = json.loads((sd / "result.json").read_text())
    except Exception:
        res = {}
    try:
        rp = json.loads((Path(result_dir) / "run_params.json").read_text())
    except Exception:
        rp = {}
    used: set = set()
    body = ['<h1>Whole genome sequencing report</h1>', f'<div class="sample">Sample: {esc(sample)}</div>',
            s_qc(sample, res, rp, root, used), s_assembly(sample, res, root, used),
            s_completeness(sample, res, root, used), s_taxonomy(sample, res, sd, root, used, tygs),
            s_ani(sample, res, root, used), s_annotation(sample, res, sd, root, used), s_amr(sample, res, root, used),
            s_mapping(sample, res, root, used)]
    order = list(REFS)
    refs = [REFS[k] for k in order if k in used]
    body.append('<h2>Reference</h2><ol class="m">' + "".join(f"<li>{esc(r)}</li>" for r in refs) + "</ol>")
    return ("<!doctype html><html><head><meta charset='utf-8'><title>Whole genome sequencing report</title>"
            f"<style>{CSS}</style></head><body>" + "\n".join(body) + "</body></html>")


def build_sample(result_dir: Path, sample: str, out_dir: Path, formats=("html", "pdf", "docx"),
                 name: str | None = None) -> dict:
    """Package + report of one sample. Returns the package manifest with report paths / errors."""
    man = WP.build_package(result_dir, sample, out_dir, name=name)
    root = Path(man["dir"])
    sd = Path(result_dir) / "samples" / sample
    tygs = None
    try:
        import tygs as TY
        tygs = TY.load(sd / "tygs")
    except Exception:
        tygs = None
    doc = build_html(result_dir, sample, root, tygs)
    stem = f"Report_{name or WP.safe(sample)}"
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
            docx_export.html_to_docx(doc, root / f"{stem}.docx", title="Whole genome sequencing report")
            man["docx"] = str(root / f"{stem}.docx")
        except Exception as ex:
            man["docx_error"] = str(ex)
    if "html" not in formats:
        (root / f"{stem}.html").unlink()
        man.pop("html", None)
    man["files"] = sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())
    return man
