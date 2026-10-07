#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wgs_report.py — HTML (and, through the app's /report endpoint, PDF) report for an
ONT-WGS isolate run. Reads only what ont_wgs_pipeline.py leaves in the result
folder (summary.json + wgs_results.json), so it can be rebuilt at any time:

    python3 wgs_report.py RESULTS_DIR [-o report.html]

Styling comes from report_builder.py so amplicon and genome reports look like
the same product.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from report_builder import (ACC, FAINT, INK, MUT, RULE, STOP, WARN, _css, _tbl,  # noqa: E402
                            esc, fmt)

GOOD = "#2E6B46"


# ── formatting ────────────────────────────────────────────────────────────────
def bp(n) -> str:
    try:
        n = float(n)
    except (TypeError, ValueError):
        return "—"
    if n >= 1e9:
        return f"{n / 1e9:.2f} Gb"
    if n >= 1e6:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.1f} kb"
    return f"{int(n)} bp"


def pct(v, d=1) -> str:
    try:
        return f"{float(v):.{d}f}%"
    except (TypeError, ValueError):
        return "—"


def num(v, d=1) -> str:
    try:
        return f"{float(v):.{d}f}"
    except (TypeError, ValueError):
        return "—"


def ital(species: str) -> str:
    """Binomials in italics, GTDB placeholder suffixes (_A) kept upright."""
    if not species:
        return "—"
    return f"<i>{esc(species)}</i>"


def call_badge(call: str) -> str:
    col = GOOD if call == "species" else (WARN if call.startswith(("genus", "16S"))
                                          else STOP)
    return f'<span style="color:{col};font-weight:700">{esc(call)}</span>'


# ── figures ───────────────────────────────────────────────────────────────────
def svg_length_hist(rp: dict, w=330, h=150) -> str:
    bins, vals = rp.get("len_bins", []), rp.get("len_hist_bp", [])
    if not bins or not any(vals):
        return ""
    labels = []
    for b in bins:
        labels.append("0" if b == 0 else (f"{b // 1000}k" if b >= 1000 else str(b)))
    mx = max(vals) or 1
    pl, pb, pt = 34, 26, 10
    bw = (w - pl - 6) / len(vals)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
           f'viewBox="0 0 {w} {h}" font-family="DejaVu Sans" font-size="7">']
    for i, v in enumerate(vals):
        bh = (h - pb - pt) * v / mx
        x = pl + i * bw
        out.append(f'<rect x="{x + 1:.1f}" y="{h - pb - bh:.1f}" width="{bw - 2:.1f}" '
                   f'height="{bh:.1f}" fill="{ACC}"/>')
        out.append(f'<text x="{x + bw / 2:.1f}" y="{h - pb + 9}" text-anchor="middle" '
                   f'fill="{MUT}">{labels[i]}</text>')
    out.append(f'<line x1="{pl}" y1="{h - pb}" x2="{w - 4}" y2="{h - pb}" stroke="{INK}" '
               f'stroke-width=".6"/>')
    out.append(f'<text x="{pl}" y="{pt - 2}" fill="{MUT}">bases per read-length bin '
               f'(max {bp(mx)})</text>')
    out.append(f'<text x="{(w + pl) / 2}" y="{h - 4}" text-anchor="middle" fill="{MUT}">'
               f'read length (bp, bin start)</text>')
    out.append("</svg>")
    return "".join(out)


def svg_gc_hist(rp: dict, asm_gc=None, w=250, h=150) -> str:
    vals = rp.get("gc_hist", [])
    if not vals or not any(vals):
        return ""
    lo, hi = 15, 80
    vals = vals[lo:hi + 1]
    mx = max(vals) or 1
    pl, pb, pt = 10, 26, 10
    bw = (w - pl - 6) / len(vals)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
           f'viewBox="0 0 {w} {h}" font-family="DejaVu Sans" font-size="7">']
    for i, v in enumerate(vals):
        bh = (h - pb - pt) * v / mx
        out.append(f'<rect x="{pl + i * bw:.1f}" y="{h - pb - bh:.1f}" width="{bw:.2f}" '
                   f'height="{bh:.1f}" fill="{GOOD}"/>')
    for g in range(20, 81, 10):
        x = pl + (g - lo + .5) * bw
        out.append(f'<text x="{x:.1f}" y="{h - pb + 9}" text-anchor="middle" '
                   f'fill="{MUT}">{g}</text>')
    if asm_gc:
        x = pl + (float(asm_gc) - lo + .5) * bw
        out.append(f'<line x1="{x:.1f}" y1="{pt}" x2="{x:.1f}" y2="{h - pb}" '
                   f'stroke="{STOP}" stroke-dasharray="2,2" stroke-width=".8"/>')
    out.append(f'<line x1="{pl}" y1="{h - pb}" x2="{w - 4}" y2="{h - pb}" stroke="{INK}" '
               f'stroke-width=".6"/>')
    out.append(f'<text x="{pl}" y="{pt - 2}" fill="{MUT}">bases by read GC</text>')
    out.append(f'<text x="{(w + pl) / 2}" y="{h - 4}" text-anchor="middle" fill="{MUT}">'
               f'read GC (%) — dashed: assembly GC</text>')
    out.append("</svg>")
    return "".join(out)


def _wedge(cx, cy, r0, r1, a0, a1) -> str:
    def pt(r, a):
        return f"{cx + r * math.cos(a):.2f},{cy + r * math.sin(a):.2f}"
    return f"M{pt(r0, a0)} L{pt(r1, a0)} L{pt(r1, a1)} L{pt(r0, a1)} Z"


def svg_genome_map(r: dict, size=330) -> str:
    """Circular map of the chromosome: GC content and GC skew per window, rRNA
    operons, AMR and virulence genes. Plasmids are listed beside it."""
    contigs = r.get("contigs", [])
    if not contigs:
        return ""
    chrom = max(contigs, key=lambda c: c["length"])
    gw = r.get("gc_windows", {}).get(chrom["contig"])
    if not gw or not gw.get("gc"):
        return ""
    L = chrom["length"]
    cx = cy = size / 2
    R = size / 2 - 30
    tau = 2 * math.pi
    ang = lambda pos: tau * pos / L - math.pi / 2   # noqa: E731
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
           f'viewBox="0 0 {size} {size}" font-family="DejaVu Sans" font-size="7">']
    # backbone + Mb ticks
    out.append(f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" stroke="{INK}" '
               f'stroke-width="1.4"/>')
    step = 500_000 if L > 2_000_000 else 250_000
    for pos in range(0, L, step):
        a = ang(pos)
        major = pos % 1_000_000 == 0
        r1 = R + (6 if major else 3)
        out.append(f'<line x1="{cx + R * math.cos(a):.1f}" y1="{cy + R * math.sin(a):.1f}" '
                   f'x2="{cx + r1 * math.cos(a):.1f}" y2="{cy + r1 * math.sin(a):.1f}" '
                   f'stroke="{INK}" stroke-width=".6"/>')
        if major:
            out.append(f'<text x="{cx + (R + 14) * math.cos(a):.1f}" '
                       f'y="{cy + (R + 14) * math.sin(a) + 2.5:.1f}" text-anchor="middle" '
                       f'fill="{MUT}">{pos // 1_000_000} Mb</text>')
    # features on the backbone
    feats = []
    for f in r.get("rrna_features", []):
        if f["contig"] == chrom["contig"] and f["name"].startswith("16S"):
            feats.append((f["start"], ACC, 9))
    for x in r.get("amrfinder", []):
        if x.get("contig") == chrom["contig"] and x.get("type") == "AMR":
            feats.append((int(float(x.get("start") or 0)), STOP, 11))
    for x in r.get("vfdb", []):
        if x.get("contig") == chrom["contig"]:
            feats.append((int(float(x.get("start") or 0)), WARN, 7))
    for pos, col, ln in feats:
        a = ang(pos)
        out.append(f'<line x1="{cx + (R - ln / 2) * math.cos(a):.1f}" '
                   f'y1="{cy + (R - ln / 2) * math.sin(a):.1f}" '
                   f'x2="{cx + (R + ln / 2) * math.cos(a):.1f}" '
                   f'y2="{cy + (R + ln / 2) * math.sin(a):.1f}" stroke="{col}" '
                   f'stroke-width="1.3"/>')
    # GC content (deviation from mean) ring
    gc = gw["gc"]
    n = len(gc)
    mean = sum(gc) / n
    dev = max(1e-6, max(abs(g - mean) for g in gc))
    r_gc = R - 26
    for i, g in enumerate(gc):
        a0, a1 = tau * i / n - math.pi / 2, tau * (i + 1) / n - math.pi / 2
        h = 16 * (g - mean) / dev
        col = GOOD if h >= 0 else "#8FA8B2"
        out.append(f'<path d="{_wedge(cx, cy, r_gc, r_gc + h, a0, a1)}" fill="{col}"/>')
    # GC skew ring
    sk = gw["skew"]
    smax = max(1e-6, max(abs(s) for s in sk))
    r_sk = R - 62
    for i, s in enumerate(sk):
        a0, a1 = tau * i / n - math.pi / 2, tau * (i + 1) / n - math.pi / 2
        h = 14 * s / smax
        col = "#6B4E8C" if h >= 0 else "#C9A227"
        out.append(f'<path d="{_wedge(cx, cy, r_sk, r_sk + h, a0, a1)}" fill="{col}"/>')
    out.append(f'<text x="{cx}" y="{cy - 4}" text-anchor="middle" font-size="9" '
               f'font-weight="bold" fill="{INK}">{esc(chrom["contig"])}</text>')
    out.append(f'<text x="{cx}" y="{cy + 8}" text-anchor="middle" font-size="8" '
               f'fill="{MUT}">{bp(L)} · {"circular" if chrom.get("circular") else "linear"}'
               f'</text>')
    out.append("</svg>")
    return "".join(out)


def map_legend() -> str:
    items = [(ACC, "16S rRNA gene"), (STOP, "AMR gene (AMRFinderPlus)"),
             (WARN, "virulence factor (VFDB)"), (GOOD, "GC content above mean"),
             ("#8FA8B2", "GC content below mean"), ("#6B4E8C", "GC skew +"),
             ("#C9A227", "GC skew −")]
    return "".join(f'<div style="font-size:7.6pt;margin:.6mm 0"><span style="display:'
                   f'inline-block;width:9px;height:9px;background:{c};margin-right:4px">'
                   f'</span>{esc(t)}</div>' for c, t in items)


def tbl(head, rows, caption=None, left=(), nowrap=()) -> str:
    """_tbl with chosen columns left-aligned (text) and/or kept on one line."""
    def st(i):
        css = []
        if i in left:
            css.append("text-align:left")
        if i in nowrap:
            css.append("white-space:nowrap")
        return f' style="{";".join(css)}"' if css else ""
    h = "".join(f"<th{st(i)}>{c}</th>" for i, c in enumerate(head))
    b = "".join("<tr>" + "".join(f"<td{st(i)}>{c}</td>" for i, c in enumerate(r)) + "</tr>"
                for r in rows)
    cap = f"<caption>{caption}</caption>" if caption else ""
    return f"<table><thead><tr>{h}</tr></thead><tbody>{b}</tbody>{cap}</table>"


def vf_product(p: str) -> str:
    """VFDB products read '(gene) description [gene (VFxxxx) - category (VFCxxxx)]
    [organism]' — keep the description and the category, drop the organism."""
    import re
    m = re.match(r"^\(\S+\)\s*(.*?)\s*\[[^\]]*?-\s*([^\[\]]+?)\s*\(VFC\d+\)\]", p or "")
    if m:
        return f"{m.group(1)} — {m.group(2)}"
    return re.sub(r"\s*\[[^\]]*\]\s*$", "", p or "")


# ── sections ──────────────────────────────────────────────────────────────────
def kpis(items) -> str:
    return '<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>'
        for v, l in items) + "</div>"


def section_sample(r: dict) -> str:
    S = r["sample"]
    if r.get("status") != "ok":
        return (f'<h2 class="pbreak">{esc(S)}</h2><div class="note"><b>Sample failed.</b> '
                f'{esc(r.get("error", ""))}</div>')
    idn, a = r.get("identification", {}), r.get("assembly", {})
    ck, ml = r.get("checkm2", {}), r.get("mlst", {})
    rq, fq = r.get("reads_raw", {}), r.get("reads_filtered", {})
    ri = r.get("run_info", {})
    out = [f'<h2 class="pbreak">{esc(S)} — {ital(idn.get("species"))}</h2>']
    src = []
    if r.get("barcode"):
        src.append(f"barcode <span class='mono'>{esc(r['barcode'])}</span>")
    if r.get("folder"):
        src.append(f"folder <span class='mono'>{esc(r['folder'])}</span>")
    src.append(f"{r.get('n_files', 0)} FASTQ file(s) merged")
    out.append(f'<p class="sub">{" · ".join(src)}</p>')
    ani = idn.get("ani")
    out.append(kpis([
        (bp(a.get("total_length")), "genome size"),
        (f"{a.get('contigs', '—')} / {a.get('circular', '—')}", "contigs / circular"),
        (f"{num(a.get('mean_depth'), 0)}x", "mean depth"),
        (pct(ani, 2) if ani is not None else "—", "ANI to closest reference"),
        (f"{num(ck.get('completeness'))} / {num(ck.get('contamination'))}"
         if ck else "—", "CheckM2 complete / contam. (%)"),
        (f"ST{esc(ml.get('st'))}" if ml.get("st") not in (None, "", "-") else "—",
         f"MLST {ml.get('scheme', '')}"),
    ]))
    if r.get("warnings"):
        out.append('<div class="note"><b>Check before reporting:</b><ul>' + "".join(
            f"<li>{esc(w)}</li>" for w in r["warnings"]) + "</ul></div>")

    # sequencing data
    out.append("<h3>Sequencing data</h3>")
    rows = []
    for lab, d in (("raw (merged)", rq), ("after filtering", fq)):
        rows.append([lab, fmt(d.get("reads")), bp(d.get("bases")), fmt(d.get("mean_len")),
                     fmt(d.get("n50")), num(d.get("mean_q")), pct(d.get("pct_q20")),
                     pct(d.get("gc"))])
    out.append(_tbl(["Reads", "Count", "Bases", "Mean len", "N50", "Mean Q", "≥Q20",
                     "GC"], rows,
                    f"Basecall model {esc(ri.get('basecall_model') or 'not recorded')}; "
                    f"flow cell {esc(ri.get('flow_cell_id') or '—')}. Filter: reads "
                    f"≥ min length and ≥ min mean Q (see Methods)."))
    peaks = rq.get("gc_peaks", [])
    out.append('<figure><div style="display:flex;gap:6mm;align-items:flex-end">'
               f'{svg_length_hist(rq)}{svg_gc_hist(rq, a.get("gc"))}</div>'
               '<figcaption>Left: sequenced bases by read length. Right: bases by read GC '
               '— a pure culture gives one peak'
               + (f"; this sample shows {len(peaks)} peaks ("
                  + ", ".join(f"{p['gc']}%" for p in peaks) + ")" if len(peaks) > 1 else
                  "; one peak here") + ".</figcaption></figure>")

    # assembly
    out.append("<h3>Assembly</h3>")
    out.append(f"<p>Flye (<span class='mono'>--{esc(r.get('flye_mode', ''))}</span>) "
               f"followed by {esc(a.get('polish', ''))}. N50 {bp(a.get('n50'))}, "
               f"GC {pct(a.get('gc'), 2)}. "
               + ("The chromosome closed into a single circular contig." if
                  a.get("chromosome_closed") else
                  "<b>The chromosome did not close</b> — the genome is in pieces; gene "
                  "content is still reliable, synteny and copy number less so.")
               + "</p>")
    crow = []
    for c in r.get("contigs", [])[:40]:
        crow.append([f"<span class='mono'>{esc(c['contig'])}</span>", esc(c.get("type", "")),
                     fmt(c["length"]), pct(c.get("gc"), 2), f"{num(c.get('depth'), 0)}x",
                     "yes" if c.get("circular") else "no", esc(c.get("replicons", "")) or
                     "—"])
    out.append(tbl(["Contig", "Type", "Length", "GC", "Depth", "Circular", "Replicons"],
                   crow, left=(1, 6), nowrap=(0,), caption=("Type: chromosome = largest replicon ≥ 1 Mb; plasmid = MOB-suite "
                           "call or a circular contig carrying a PlasmidFinder replicon.")
                    + (f" {len(r['contigs']) - 40} more contigs in assembly_contigs.csv."
                       if len(r.get("contigs", [])) > 40 else "")))
    gm = svg_genome_map(r)
    if gm:
        out.append(f'<figure><div style="display:flex;gap:6mm;align-items:center">{gm}'
                   f'<div>{map_legend()}</div></div><figcaption>Chromosome map. Outer '
                   f'ring: position (Mb) with 16S, AMR and virulence genes; middle: GC '
                   f'content per window relative to the genome mean; inner: GC skew '
                   f'(its sign change marks the replication origin and terminus).'
                   f'</figcaption></figure>')

    # identification
    out.append("<h3>Species identification</h3>")
    if idn.get("method"):
        line = (f"<p><b>{ital(idn.get('species'))}</b> — {call_badge(idn.get('call', ''))}"
                f" by {esc(idn['method'])}")
        if ani is not None:
            line += f", ANI {pct(ani, 2)}"
        if idn.get("af") is not None:
            line += f" over {pct(float(idn['af']), 1)} of the assembly"
        if idn.get("reference"):
            line += (f" to <span class='mono'>{esc(idn.get('accession', ''))}</span> "
                     f"({esc(idn['reference'])})")
        out.append(line + ".</p>")
        if idn.get("lineage"):
            out.append(f"<p class='dim' style='font-size:8pt'>GTDB lineage: "
                       f"{esc(idn['lineage'])}</p>")
    sm = r.get("sourmash", [])
    if sm:
        rows = [[str(h["rank"]), ital(h.get("gtdb_species") or h.get("reference_name")),
                 f"<span class='mono'>{esc(h['accession'])}</span>", pct(h.get("ani"), 2),
                 pct(h.get("f_query"), 1)] for h in sm[:6]]
        out.append(tbl(["#", "GTDB species", "Reference", "ANI (est.)", "% of assembly"],
                       rows, left=(1, 2), caption= "sourmash gather against GTDB species representatives "
                        "(k=31). A second species holding a real share of the assembly "
                        "points to a mixed culture."))
    out.append('<div class="note acc">ANI ≥ 95 % to a species representative is the '
               'accepted species boundary for prokaryotes; 16S rRNA identity alone cannot '
               'separate many close species (e.g. <i>Cronobacter sakazakii</i> / <i>C. '
               'malonaticus</i>, the <i>Bacillus cereus</i> group, <i>Enterobacter cloacae'
               '</i> complex), so the genome call takes precedence.</div>')

    # 16S
    rr = r.get("rrna", {})
    out.append("<h3>16S rRNA genes</h3>")
    if rr.get("copies"):
        out.append(f"<p>{rr['copies']} full-length 16S copies in the assembly "
                   f"({rr.get('distinct_sequences', '?')} distinct sequence(s)); "
                   f"{rr.get('operons_23S', 0)} 23S and {rr.get('operons_5S', 0)} 5S genes"
                   + (f"; {rr['partial']} partial 16S fragment(s) not used" if rr.get("partial")
                      else "") + ". Sequences delivered as "
                   f"<span class='mono'>{esc(S)}_16S.fasta</span>.</p>")
        hits = rr.get("hits", [])
        if rr.get("consensus_species"):
            if rr.get("resolved"):
                out.append(f"<p>All copies agree on {ital(rr['consensus_species'])}.</p>")
            else:
                out.append(f"<p><b>16S alone resolves only to genus:</b> "
                           f"{esc(rr['consensus_species'])}. The ‘best match’ column differs "
                           f"by ≤ 0.2 % between these species and must not be read as an "
                           f"identification — the genome result above decides.</p>")
        if hits:
            rows = [[f"<span class='mono'>{esc(h['copy'])}</span>", ital(h.get("best_species")),
                     pct(h.get("best_identity"), 2), str(h.get("n_tied", "")),
                     esc(h.get("tied_species", ""))[:160]] for h in hits]
            out.append(tbl(["Copy", "Best 16S match", "Identity", "Tied", "Species within "
                            "0.2 % of best"], rows, left=(1, 4), nowrap=(0,), caption=
                            f"vsearch global alignment against {esc(rr.get('db', ''))}."))
    else:
        out.append("<p class='dim'>No full-length 16S gene recovered.</p>")

    # MLST
    out.append("<h3>MLST</h3>")
    if ml:
        out.append(_tbl(["Scheme", "ST", "Alleles"], [[esc(ml.get("scheme")),
                    esc(ml.get("st")), f"<span class='mono'>{esc(ml.get('alleles'))}</span>"]],
                        "'-' = no scheme or novel allele combination; '~' = novel allele; "
                        "'?' = partial match."))
    else:
        out.append("<p class='dim'>Not run.</p>")

    # AMR
    amr = r.get("amrfinder", [])
    out.append("<h3>Antimicrobial resistance</h3>")
    core = [x for x in amr if x.get("type") == "AMR"]
    if amr or r.get("amr_organism") is not None:
        org = r.get("amr_organism")
        out.append(f"<p>AMRFinderPlus (<span class='mono'>--plus</span>"
                   + (f", organism <i>{esc(org)}</i>: point mutations screened" if org else
                      ", no organism-specific point-mutation set for this species")
                   + f"): {len(core)} AMR determinant(s).</p>")
        if core:
            rows = [[f"<b>{esc(x['gene'])}</b>", esc(x.get("class", "")),
                     esc(x.get("subclass", "")), esc(x.get("method", "")),
                     pct(x.get("identity")), pct(x.get("coverage")),
                     f"<span class='mono'>{esc(x.get('contig', ''))}</span>"] for x in core]
            out.append(tbl(["Gene", "Class", "Subclass", "Method", "Identity", "Coverage",
                            "Contig"], rows, left=(1, 2, 3), nowrap=(6,), caption=
                            "Genotype only: presence of a gene predicts, but does not prove, "
                            "phenotypic resistance — confirm by AST where it matters."))
        other = [x for x in amr if x.get("type") in ("STRESS", "VIRULENCE")]
        if other:
            rows = [[esc(x.get("type", "")).title(), f"<b>{esc(x['gene'])}</b>",
                     esc(x.get("name", ""))[:70], esc(x.get("subclass", "")),
                     pct(x.get("identity"))] for x in other[:40]]
            out.append(tbl(["Type", "Gene", "Name", "Subclass", "Identity"], rows,
                           left=(1, 2, 3), caption=
                            "Stress-response (metal, biocide, heat, acid) and virulence genes "
                            "from the AMRFinderPlus 'plus' set."))
    else:
        out.append("<p class='dim'>Not run.</p>")

    # VFDB
    vf = r.get("vfdb", [])
    out.append("<h3>Virulence factors (VFDB)</h3>")
    if vf:
        rows = [[f"<b>{esc(x['gene'])}</b>", esc(vf_product(x.get("product", "")))[:110],
                 pct(x.get("identity")), pct(x.get("coverage")),
                 f"<span class='mono'>{esc(x.get('contig', ''))}</span>"] for x in vf[:60]]
        out.append(tbl(["Gene", "Product — VFDB category", "Identity", "Coverage", "Contig"],
                       rows, f"abricate + VFDB, ≥ 80 % identity and coverage. {len(vf)} hit(s)"
                       + (" — first 60 shown; all in virulence_vfdb.csv." if len(vf) > 60
                          else ".") + " Many are genes shared with related non-pathogens "
                       "(flagella, chemotaxis, LPS); presence alone is not proof of "
                       "virulence.", left=(1,), nowrap=(0, 4)))
    else:
        out.append("<p class='dim'>No VFDB hits at ≥ 80 % identity / coverage "
                   "(or screen not run).</p>")

    # plasmids
    out.append("<h3>Plasmids</h3>")
    mob = r.get("mobsuite", {}).get("plasmids", [])
    pf = r.get("plasmidfinder", [])
    if mob:
        rows = [[esc(p.get("plasmid", "")), fmt(p.get("size")), esc(p.get("rep_types", "")),
                 esc(p.get("relaxase", "")), esc(p.get("mobility", "")),
                 esc(p.get("host_range", ""))] for p in mob]
        out.append(_tbl(["Plasmid", "Size", "Replicon", "Relaxase", "Mobility",
                         "Predicted host range"], rows, "MOB-suite mob_recon / mob_typer."))
    if pf:
        rows = [[f"<b>{esc(x['gene'])}</b>", f"<span class='mono'>{esc(x['contig'])}</span>",
                 pct(x.get("identity")), pct(x.get("coverage"))] for x in pf]
        out.append(_tbl(["Replicon", "Contig", "Identity", "Coverage"], rows,
                        "PlasmidFinder replicons (Enterobacterales and Gram-positive "
                        "schemes; other taxa may carry plasmids it does not recognise)."))
    if not mob and not pf:
        out.append(f"<p>{a.get('plasmids', 0)} plasmid(s) called; no replicon matched "
                   "PlasmidFinder.</p>")

    # annotation
    bk = r.get("bakta", {}).get("summary", {})
    out.append("<h3>Genome annotation</h3>")
    if bk:
        keys = [k for k in ("CDSs", "hypotheticals", "pseudogenes", "tRNAs", "tmRNAs",
                            "rRNAs", "ncRNAs", "ncRNA regions", "CRISPR arrays", "sORFs",
                            "gaps", "oriCs", "oriVs", "oriTs") if k in bk]
        out.append(_tbl(["Feature", "Count"], [[esc(k), fmt(bk[k])] for k in keys],
                        f"Bakta. GenBank/GFF3/protein files: samples/{esc(S)}/bakta/."))
    else:
        out.append("<p class='dim'>Not run (Bakta or its database not installed).</p>")
    return "\n".join(out)


def build_html(out_dir: Path, title: str | None = None, company: str = "",
               app: str = "NextGen-Amplicon", subtitle: str | None = None) -> str:
    out_dir = Path(out_dir)
    summ = json.loads((out_dir / "summary.json").read_text())
    results = json.loads((out_dir / "wgs_results.json").read_text())
    ok = [r for r in results if r.get("status") == "ok"]
    job = summ.get("job_name", out_dir.name)
    if summ.get("app_version") and "v" not in app:
        app = f"{app} v{summ['app_version']}"
    title = title or "Whole-genome sequencing report — bacterial isolates"
    ri = next((r.get("run_info", {}) for r in ok if r.get("run_info")), {})
    body = [f'<div class="eyebrow">{esc(job)} &middot; ONT-WGS &middot; '
            f'{esc(summ.get("timestamp", "")[:10])}</div>', f"<h1>{esc(title)}</h1>",
            f'<p class="sub">{esc(subtitle) if subtitle else ""}'
            f'{len(results)} sample(s) · Oxford Nanopore '
            f'{esc(ri.get("chemistry", ""))} {esc(ri.get("accuracy", "").upper())} · '
            f'flow cell {esc(ri.get("flow_cell_id") or "—")}</p>']
    closed = sum(1 for r in ok if r["assembly"].get("chromosome_closed"))
    sp_ok = sum(1 for r in ok if r["identification"].get("call") == "species")
    body.append(kpis([(f"{len(ok)}/{len(results)}", "samples assembled"),
                      (str(closed), "closed chromosomes"),
                      (str(sp_ok), "species-level IDs (ANI ≥ 95 %)"),
                      (f"{sum(len([x for x in r.get('amrfinder', []) if x.get('type') == 'AMR']) for r in ok)}",
                       "AMR determinants (all samples)")]))
    rows = []
    for r in results:
        if r.get("status") != "ok":
            rows.append([f"<b>{esc(r['sample'])}</b>", esc(r.get("barcode", "")),
                         f"<span style='color:{STOP}'>failed</span>"] + ["—"] * 6)
            continue
        idn, a, ck = r["identification"], r["assembly"], r.get("checkm2", {})
        nam = len([x for x in r.get("amrfinder", []) if x.get("type") == "AMR"])
        rows.append([f"<b>{esc(r['sample'])}</b>", esc(r.get("barcode", "")),
                     ital(idn.get("species")) + f"<br><span class='dim' style='font-size:7pt'>"
                     f"{esc(idn.get('call', ''))}</span>",
                     pct(idn.get("ani"), 2) if idn.get("ani") is not None else "—",
                     bp(a.get("total_length")),
                     f"{a.get('contigs')} ({a.get('circular')} circ.)",
                     f"{num(a.get('mean_depth'), 0)}x",
                     (f"{num(ck.get('completeness'))}/{num(ck.get('contamination'))}"
                      if ck else "—"),
                     f"{esc(r.get('mlst', {}).get('st') or '—')} · {nam} · "
                     f"{a.get('plasmids', 0)}"])
    body.append("<h2>Summary</h2>")
    body.append(_tbl(["Sample", "Barcode", "Species", "ANI", "Genome", "Contigs", "Depth",
                      "Compl./Cont. %", "ST · AMR · plasmids"], rows,
                     "Species from whole-genome ANI against GTDB species representatives "
                     "where available, otherwise from 16S rRNA (marked)."))
    if summ.get("warnings"):
        gen = [w for w in summ["warnings"] if not w.startswith("[")]
        if gen:
            body.append('<div class="note"><b>Run notes:</b><ul>' + "".join(
                f"<li>{esc(w)}</li>" for w in gen) + "</ul></div>")
    for r in results:
        body.append(section_sample(r))

    # methods
    v = summ.get("tool_versions", {})
    st = summ.get("settings", {})
    dbs = summ.get("databases", {})

    def tv(t):
        return f" {esc(v[t])}" if t in v else " (not installed)"
    body.append('<h2 class="pbreak">Methods</h2>')
    body.append(
        "<p>MinKNOW chunk files were merged per sample (file list in "
        "<span class='mono'>sample_map.csv</span>). Reads were profiled and filtered with "
        f"seqkit{tv('seqkit')} (length ≥ {esc(st.get('min_read_len'))} bp, mean "
        f"Q ≥ {esc(st.get('min_read_q'))}). Genomes were assembled de novo with "
        f"Flye{tv('flye')}"
        + (f" using the longest {esc(st.get('asm_coverage'))}x of reads for the initial "
           f"disjointigs" if st.get("asm_coverage") else "")
        + f", and polished with Medaka{tv('medaka')} "
        "(bacterial methylation-aware model where compatible with the basecaller). "
        f"Completeness and contamination: CheckM2{tv('checkm2')}. Species: sourmash"
        f"{tv('sourmash')} gather against GTDB species representatives, with ANI to the "
        f"closest reference genome confirmed by skani{tv('skani')} when the reference "
        f"could be retrieved from NCBI. 16S rRNA genes: barrnap{tv('barrnap')}, "
        f"classified with vsearch{tv('vsearch')}. MLST: mlst{tv('mlst')} (PubMLST "
        f"schemes). AMR, stress and virulence genes: AMRFinderPlus{tv('amrfinder')} "
        f"with --plus. Virulence factors: abricate{tv('abricate')} with VFDB; plasmid "
        f"replicons with PlasmidFinder. Plasmid reconstruction: MOB-suite"
        f"{tv('mob_recon')}. Annotation: Bakta{tv('bakta')}.</p>")
    dbrows = [[esc(k), f"<span class='mono'>{esc(Path(p).name if p else '—')}</span>"]
              for k, p in dbs.items()]
    if dbrows:
        body.append(_tbl(["Database", "Version / file"], dbrows))
    body.append(
        "<h3>Interpretation thresholds</h3><ul>"
        "<li>Species: ANI ≥ 95 % with ≥ ~60 % of the genome aligned; 90–95 % "
        "suggests a related but different species.</li>"
        "<li>Isolate purity: CheckM2 contamination ≤ 5 %, one read-GC peak, and no "
        "second species with ≥ 5 % of the assembly.</li>"
        "<li>Assembly: a closed (circular) chromosome at ≥ 50x is complete; "
        "fragmented assemblies remain valid for gene content.</li>"
        "<li>AMR genotypes predict phenotype imperfectly and need AST to confirm.</li>"
        "</ul>")
    body.append(
        "<h3>Files delivered</h3><ul>"
        "<li><span class='mono'>wgs_summary.csv</span> — one row per sample</li>"
        "<li><span class='mono'>samples/&lt;sample&gt;/&lt;sample&gt;_assembly.fasta</span>"
        " — polished genome; <span class='mono'>_16S.fasta</span> — 16S copies; "
        "<span class='mono'>bakta/</span> — annotation (GenBank, GFF3, proteins)</li>"
        "<li><span class='mono'>assembly_contigs.csv, species_id_gtdb.csv, "
        "rrna_16S_hits.csv, mlst.csv, amr_stress_virulence_amrfinder.csv, "
        "virulence_vfdb.csv, plasmids_mobsuite.csv, read_qc.csv, sample_map.csv</span></li>"
        "</ul>")
    css = _css(app, company, f"{job} · ONT-WGS")
    return (f'<!DOCTYPE html><html><head><meta charset="utf-8"><title>{esc(title)}</title>'
            f'<style>{css}\n@media screen {{ body {{ max-width: 190mm; margin: 10mm auto; '
            f'padding: 0 6mm; }} }}</style></head><body>' + "\n".join(body)
            + "</body></html>")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("-o", "--out", default="")
    a = ap.parse_args(argv)
    html = build_html(Path(a.results))
    outp = Path(a.out) if a.out else Path(a.results) / "wgs_report.html"
    outp.write_text(html, encoding="utf-8")
    print(outp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
