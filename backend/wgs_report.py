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
    col = GOOD if call == "species" else (
        WARN if call.startswith(("genus", "16S", "species (")) else
        ACC if call.startswith("potential novel") else STOP)
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


# ── type strains, trees, mapping ─────────────────────────────────────────────
def _split_label(name: str):
    """'Genus species [subsp. x] strainT ACC' → (binomial, strain, is_type, acc)."""
    import re
    w = name.split()
    acc = w.pop() if w and re.match(r"^(GC[AF]_\d+|NR_\d+)", w[-1]) else ""
    n = 4 if len(w) >= 4 and w[2] in ("subsp.", "pv.", "bv.", "serovar") else 2
    sp, rest = " ".join(w[:n]), " ".join(w[n:])
    is_t = rest.endswith("T")
    if is_t:
        rest = rest[:-1]
    return sp, rest, is_t, acc


def svg_tree(newick: str, isolate: str, w: int = 640, title: str = "") -> str:
    """Rectangular phylogram; the isolate in accent colour, type strains with ᵀ."""
    if not newick:
        return ""
    import wgs_taxonomy as WT
    try:
        root = WT.parse_newick(newick)
    except Exception:
        return ""
    lay = WT.tree_layout(root)
    rows, mx = lay["rows"], lay["max_x"] or 1.0
    rh, top, left = 13, 16, 8
    lab_w = 330
    tw = w - lab_w - left - 10
    h = top + rows * rh + 26

    def X(x):
        return left + tw * x / mx

    def Y(y):
        return top + y * rh + rh / 2
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
         f'viewBox="0 0 {w} {h}" font-family="DejaVu Sans" font-size="7.4">']
    if title:
        o.append(f'<text x="{left}" y="10" fill="{MUT}" font-size="7">{esc(title)}</text>')
    for e in lay["edges"]:
        if "y0" in e:
            o.append(f'<line x1="{X(e["x0"]):.1f}" y1="{Y(e["y0"]):.1f}" x2="{X(e["x0"]):.1f}" '
                     f'y2="{Y(e["y1"]):.1f}" stroke="{INK}" stroke-width=".7"/>')
            sup = e.get("support")
            if sup is not None and sup >= 50 and e["x0"] > 0:
                o.append(f'<text x="{X(e["x0"]) - 2:.1f}" y="{Y(e["yn"]) - 2:.1f}" '
                         f'text-anchor="end" fill="{MUT}" font-size="6">{sup:.0f}</text>')
        else:
            o.append(f'<line x1="{X(e["x0"]):.1f}" y1="{Y(e["y"]):.1f}" x2="{X(e["x1"]):.1f}" '
                     f'y2="{Y(e["y"]):.1f}" stroke="{INK}" stroke-width=".7"/>')
    for lf in lay["leaves"]:
        x, y = X(lf["x"]) + 3, Y(lf["y"]) + 2.6
        nm = lf["name"]
        if nm.startswith(isolate):
            o.append(f'<circle cx="{X(lf["x"]):.1f}" cy="{Y(lf["y"]):.1f}" r="2.4" fill="{ACC}"/>'
                     f'<text x="{x + 2:.1f}" y="{y:.1f}" fill="{ACC}" font-weight="bold">'
                     f'{esc(nm)}</text>')
            continue
        sp, strain, is_t, acc = _split_label(nm)
        t = (f'<tspan font-style="italic">{esc(sp)}</tspan>'
             + (f' {esc(strain)}' if strain else "")
             + ('<tspan baseline-shift="super" font-size="5.5">T</tspan>' if is_t else "")
             + (f' <tspan fill="{MUT}" font-size="6.4">{esc(acc)}</tspan>' if acc else ""))
        o.append(f'<text x="{x:.1f}" y="{y:.1f}" fill="{INK}">{t}</text>')
    # scale bar
    nice = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    sb = next((v for v in nice if tw * v / mx >= 30), nice[-1])
    yb = h - 10
    o.append(f'<line x1="{left}" y1="{yb}" x2="{left + tw * sb / mx:.1f}" y2="{yb}" '
             f'stroke="{INK}" stroke-width=".8"/><text x="{left + tw * sb / mx + 4:.1f}" '
             f'y="{yb + 2.5}" fill="{MUT}" font-size="6.5">{sb:g}</text>')
    o.append("</svg>")
    return "".join(o)


def svg_coverage(m: dict, w: int = 640, h: int = 150) -> str:
    bins = m.get("coverage_bins") or []
    seqs = m.get("ref_seqs") or []
    if not bins:
        return ""
    if not seqs:
        lens = {}
        for c, st_, _ in bins:
            lens[c] = max(lens.get(c, 0), st_ + m.get("bin_size", 1000))
        seqs = list(lens.items())
    off, tot = {}, 0
    for c, L in seqs:
        off[c] = tot
        tot += L
    tot = tot or 1
    mean = m.get("mean_depth") or (sum(d for _, _, d in bins) / len(bins))
    cap = max(1.0, 3 * (mean or 1))
    pl, pr, pt, pb = 34, 8, 12, 30
    W, Hh = w - pl - pr, h - pt - pb

    def X(c, pos):
        return pl + W * (off.get(c, 0) + pos) / tot
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
         f'viewBox="0 0 {w} {h}" font-family="DejaVu Sans" font-size="7">']
    pts = []
    for c, st_, d in bins:
        y = pt + Hh * (1 - min(d, cap) / cap)
        pts.append(f"{X(c, st_):.1f},{y:.1f}")
    o.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{ACC}" '
             f'stroke-width=".8"/>')
    ym = pt + Hh * (1 - min(mean, cap) / cap)
    o.append(f'<line x1="{pl}" y1="{ym:.1f}" x2="{w - pr}" y2="{ym:.1f}" stroke="{MUT}" '
             f'stroke-dasharray="2,2" stroke-width=".5"/>')
    for c, s0, e0 in m.get("absent") or []:
        x0, x1 = X(c, s0), X(c, e0)
        o.append(f'<rect x="{x0:.1f}" y="{pt + Hh + 2}" width="{max(.8, x1 - x0):.1f}" '
                 f'height="5" fill="{STOP}"/>')
    for c, L in seqs:
        if off[c] > 0 and L / tot > 0.002:
            o.append(f'<line x1="{X(c, 0):.1f}" y1="{pt}" x2="{X(c, 0):.1f}" y2="{pt + Hh}" '
                     f'stroke="{RULE}" stroke-width=".4"/>')
    step = 0.5 if tot <= 3e6 else (1 if tot <= 8e6 else 2)
    v = 0.0
    while v * 1e6 <= tot:
        x = pl + W * v * 1e6 / tot
        o.append(f'<text x="{x:.1f}" y="{h - pb + 17}" text-anchor="middle" '
                 f'fill="{MUT}">{v:g}</text>')
        v += step
    o.append(f'<line x1="{pl}" y1="{pt + Hh}" x2="{w - pr}" y2="{pt + Hh}" stroke="{INK}" '
             f'stroke-width=".6"/>')
    o.append(f'<text x="{pl - 3}" y="{pt + 4}" text-anchor="end" fill="{MUT}">{cap:.0f}x</text>'
             f'<text x="{pl - 3}" y="{pt + Hh}" text-anchor="end" fill="{MUT}">0</text>')
    o.append(f'<text x="{pl}" y="{pt - 3}" fill="{MUT}">read depth along the reference '
             f'(dashed: mean {mean:.0f}x; red: ≥ 1 kb with ≤ 2x — absent from the isolate)'
             f'</text><text x="{pl + W / 2:.0f}" y="{h - 3}" text-anchor="middle" '
             f'fill="{MUT}">reference position (Mb)</text>')
    o.append("</svg>")
    return "".join(o)


def svg_dotplot(m: dict, size: int = 300) -> str:
    blocks = m.get("synteny_blocks") or []
    rs, qs = m.get("ref_seqs") or [], m.get("query_seqs") or []
    if not blocks or not rs or not qs:
        return ""
    ro, qo, rt, qt = {}, {}, 0, 0
    for c, L in rs:
        ro[c] = rt
        rt += L
    for c, L in qs:
        qo[c] = qt
        qt += L
    pl, pb, pt, pr = 30, 26, 8, 8
    W, Hh = size - pl - pr, size - pt - pb

    def X(c, p):
        return pl + W * (ro.get(c, 0) + p) / (rt or 1)

    def Y(c, p):
        return pt + Hh - Hh * (qo.get(c, 0) + p) / (qt or 1)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
         f'viewBox="0 0 {size} {size}" font-family="DejaVu Sans" font-size="7">',
         f'<rect x="{pl}" y="{pt}" width="{W}" height="{Hh}" fill="none" stroke="{INK}" '
         f'stroke-width=".5"/>']
    for c, L in qs:
        if qo[c] > 0 and L / (qt or 1) > 0.003:
            o.append(f'<line x1="{pl}" y1="{Y(c, 0):.1f}" x2="{pl + W}" y2="{Y(c, 0):.1f}" '
                     f'stroke="{RULE}" stroke-width=".4"/>')
    for c, L in rs:
        if ro[c] > 0 and L / (rt or 1) > 0.003:
            o.append(f'<line x1="{X(c, 0):.1f}" y1="{pt}" x2="{X(c, 0):.1f}" y2="{pt + Hh}" '
                     f'stroke="{RULE}" stroke-width=".4"/>')
    for tn, ts, te, qn, qs_, qe, strand in blocks:
        if strand == "+":
            x0, y0, x1, y1, col = X(tn, ts), Y(qn, qs_), X(tn, te), Y(qn, qe), ACC
        else:
            x0, y0, x1, y1, col = X(tn, ts), Y(qn, qe), X(tn, te), Y(qn, qs_), STOP
        o.append(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
                 f'stroke="{col}" stroke-width="1"/>')
    o.append(f'<text x="{pl + W / 2:.0f}" y="{size - 6}" text-anchor="middle" fill="{MUT}">'
             f'reference ({rt / 1e6:.2f} Mb)</text>'
             f'<text x="9" y="{pt + Hh / 2:.0f}" text-anchor="middle" fill="{MUT}" '
             f'transform="rotate(-90 9 {pt + Hh / 2:.0f})">isolate ({qt / 1e6:.2f} Mb)</text>')
    o.append("</svg>")
    return "".join(o)


def _cell(v, ok, warn=None, d=2, suffix=""):
    if v is None or v == "":
        return "—"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return esc(str(v))
    col = GOOD if ok(f) else (WARN if warn and warn(f) else INK)
    wt = "700" if ok(f) else "400"
    return f'<span style="color:{col};font-weight:{wt}">{f:.{d}f}{suffix}</span>'


def section_typestrain(r: dict) -> str:
    t = r.get("typestrain") or {}
    if not t or (not t.get("type_strains") and not t.get("rrna_type")):
        return ""
    S = r["sample"]
    out = ["<h3>Comparison with type strains (JSpecies / TYGS-style)</h3>"]
    v = t.get("verdict") or {}
    if v:
        col = {"known": GOOD, "borderline": WARN, "novel": ACC}.get(v.get("level"), MUT)
        out.append(f'<div class="note" style="border-left:3px solid {col}"><b style="color:'
                   f'{col}">{esc(v.get("call", "").capitalize())}</b> — {esc(v.get("text", ""))}'
                   f'</div>')
    ts = t.get("type_strains") or []
    if ts:
        rows = []
        for i, p in enumerate(ts, 1):
            sp, strain = p.get("organism", ""), (p.get("strain") or "")
            rows.append([
                str(i), f"<i>{esc(sp)}</i> {esc(strain)}<sup>T</sup>",
                f"<span class='mono'>{esc(p.get('accession', ''))}</span>",
                _cell(p.get("anib"), lambda x: x >= 95, lambda x: x >= 90),
                _cell(p.get("anib_cov_q"), lambda x: False, d=1),
                _cell(p.get("anim"), lambda x: x >= 95, lambda x: x >= 90),
                _cell(p.get("tetra"), lambda x: x >= 0.999, lambda x: x >= 0.989, d=4),
                _cell(p.get("skani_ani"), lambda x: x >= 95, lambda x: x >= 90),
                _cell(p.get("rrna_identity"), lambda x: x >= 98.65, d=2)])
        out.append(tbl(["#", "Type strain", "Genome", "ANIb %", "aligned %", "ANIm %",
                        "TETRA", "skani %", "16S %"], rows, left=(1, 2), nowrap=(2,),
                       caption="ANIb: mean identity of 1020-nt fragments (blastn, ≥ 70 % "
                       "coverage, ≥ 30 % identity), average of both directions; aligned = "
                       "share of the isolate genome in accepted fragments. ANIm: nucmer "
                       "--mum + delta-filter -1. TETRA: tetranucleotide z-score correlation. "
                       "Green = above the species threshold (ANI ≥ 95 %, TETRA ≥ 0.999, "
                       "16S ≥ 98.65 %); amber = transition zone. Type-strain genomes: NCBI "
                       "'assembly from type material'."))
    gt, rt = t.get("genome_tree"), t.get("rrna_tree")
    if gt:
        out.append(f'<figure>{svg_tree(gt, S)}<figcaption>Genome tree of the isolate and '
                   f'its closest type strains — {esc(t.get("genome_tree_method", ""))}. '
                   f'Comparable to the TYGS genome tree, but built from ANI distances '
                   f'rather than GBDP; no branch support.</figcaption></figure>')
    if rt:
        out.append(f'<figure>{svg_tree(rt, S)}<figcaption>16S rRNA gene tree — '
                   f'{esc(t.get("rrna_tree_method", ""))}; numbers are bootstrap support '
                   f'(≥ 50 shown). Type-strain sequences: NCBI RefSeq 16S.</figcaption>'
                   f'</figure>')
    r16 = t.get("rrna_type") or []
    if r16:
        rows = [[f"<i>{esc(x['species'])}</i> {esc(x.get('strain', ''))}<sup>T</sup>",
                 f"<span class='mono'>{esc(x['accession'])}</span>",
                 _cell(x.get("identity"), lambda y: y >= 98.65, d=2),
                 str(x.get("aln_len", ""))] for x in r16[:10]]
        out.append(tbl(["Closest type strains by 16S", "Record", "Identity %", "Aligned bp"],
                       rows, left=(0, 1), caption="Best of all 16S copies against NCBI "
                       "RefSeq 16S (type material). ≥ 98.65 % is compatible with the same "
                       "species but does not prove it."))
    out.append('<p class="dim" style="font-size:7.6pt">Files: <span class="mono">taxonomy/'
               'type_strain_ani.csv, type_strain_16S.csv, ani_matrix_skani.tsv, '
               f'{esc(S)}_genome_tree.nwk, {esc(S)}_16S_tree.nwk, {esc(S)}_16S_alignment.fasta'
               '</span>. dDDH (TYGS/GGDC) is not computed here — submit the assembly to TYGS '
               'when a formal dDDH value is required.</p>')
    return "".join(out)


def section_mapping(r: dict) -> str:
    m = r.get("mapping") or {}
    if not m:
        return ""
    out = ["<h3>Mapping to the closest reference genome</h3>"]
    ref = (f"<i>{esc(m.get('reference', ''))}</i>"
           + (f" {esc(m.get('strain', ''))}" if m.get("strain") and m.get("strain") not in
              m.get("reference", "") else "")
           + ("<sup>T</sup>" if m.get("type_strain") else "")
           + f" (<span class='mono'>{esc(m.get('accession', ''))}</span>, "
             f"{bp(m.get('ref_length'))}, {m.get('ref_contigs', '?')} sequence(s))")
    out.append(f"<p>Reference: {ref}.</p>")
    rows = [
        ["Reads mapped", pct(m.get("reads_mapped_pct"), 1),
         "Assembly aligned to reference", pct(m.get("asm_query_aligned_pct"), 1)],
        ["Mean depth on reference", f"{num(m.get('mean_depth'), 0)}x",
         "Reference covered by assembly", pct(m.get("asm_ref_aligned_pct"), 1)],
        ["Reference covered ≥ 1x / ≥ 10x",
         f"{pct(m.get('breadth_1x'), 1)} / {pct(m.get('breadth_10x'), 1)}",
         "Identity in aligned blocks", pct(m.get("alignment_identity"), 2)],
        ["Reference regions absent (≥ 1 kb)",
         f"{m.get('absent_regions', '—')} ({bp(m.get('absent_bp'))})",
         "SNPs / indels vs reference",
         f"{fmt(m.get('snps'))} / {fmt(m.get('indels'))} "
         f"({num(m.get('snps_per_100kb'), 0)} SNPs per 100 kb)"],
        ["", "", "Isolate-specific regions (≥ 5 kb)",
         f"{m.get('isolate_unique_regions', '—')} ({bp(m.get('isolate_unique_bp'))})"],
    ]
    out.append(tbl(["Reads → reference", "", "Assembly → reference", ""], rows,
                   left=(0, 2)))
    cv = svg_coverage(m)
    if cv:
        out.append(f"<figure>{cv}<figcaption>Read depth along the reference "
                   f"(minimap2 map-ont, {fmt(m.get('bin_size'))}-bp windows). Gaps mark "
                   f"reference genes the isolate lacks; peaks mark repeats or plasmid "
                   f"copy number.</figcaption></figure>")
    dp = svg_dotplot(m)
    if dp:
        out.append(f'<figure style="display:flex;gap:5mm;align-items:center">{dp}'
                   f'<figcaption style="max-width:70mm">Synteny dot plot of the polished '
                   f'assembly against the reference (minimap2 asm20, blocks ≥ 2 kb). '
                   f'<span style="color:{ACC}">Blue</span>: same orientation; '
                   f'<span style="color:{STOP}">red</span>: inverted. A straight diagonal '
                   f'means conserved gene order (a diagonal split in two on a circular '
                   f'chromosome only reflects a different start position); other breaks are '
                   f'rearrangements, insertions or replicons the other genome lacks.'
                   f'</figcaption></figure>')
    out.append('<p class="dim" style="font-size:7.6pt">SNP counts against a reference '
               'of the same species at 97–99 % ANI measure species-level divergence '
               '(thousands of SNPs are normal) — they are not an outbreak SNP distance, '
               'which needs isolates compared with each other. Files: <span class="mono">'
               'mapping/reference_coverage_bins.tsv, reference_absent_regions.csv, '
               'isolate_unique_regions.csv, assembly_vs_reference.paf</span>.</p>')
    return "".join(out)


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
                           "call or a circular contig carrying a PlasmidFinder replicon; "
                           "plasmid (putative) = circular replicon < 1 Mb with no MOB-suite / "
                           "PlasmidFinder support (common for Bacillales megaplasmids, which "
                           "these Enterobacterales-centred databases miss; may also be a "
                           "circular prophage).")
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
        if idn.get("gtdb_species") and idn.get("gtdb_species") != idn.get("species"):
            out.append(f"<p class='dim' style='font-size:8pt'>GTDB name: "
                       f"{ital(idn['gtdb_species'])}</p>")
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
        if rr.get("off_genus"):
            out.append(f"<p class='dim'>Also matched at the same identity: "
                       f"{esc(', '.join(rr['off_genus'][:4]))} — a minority of reference names "
                       f"outside the majority genus, most likely mislabelled database entries; "
                       f"ignored for the genus call.</p>")
        if hits:
            rows = [[f"<span class='mono'>{esc(h['copy'])}</span>", ital(h.get("best_species")),
                     pct(h.get("best_identity"), 2), str(h.get("n_tied", "")),
                     esc(h.get("tied_species", ""))[:160]] for h in hits]
            out.append(tbl(["Copy", "Best 16S match", "Identity", "Tied", "Species within "
                            "0.2 % of best"], rows, left=(1, 4), nowrap=(0,), caption=
                            f"vsearch global alignment against {esc(rr.get('db', ''))}."))
    else:
        out.append("<p class='dim'>No full-length 16S gene recovered.</p>")

    out.append(section_typestrain(r))
    out.append(section_mapping(r))

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
    put = [c for c in r.get("contigs", []) if c.get("type") == "plasmid (putative)"
           or (c.get("type") in ("contig", "plasmid?") and c.get("circular")
               and c.get("length", 0) < 1_000_000)]
    if put:
        out.append(f"<p>{len(put)} putative plasmid(s) — circular replicons below 1 Mb "
                   "without a MOB-suite or PlasmidFinder match: " + ", ".join(
                       f"<span class='mono'>{esc(c['contig'])}</span> ({fmt(c['length'])} bp)"
                       for c in put) + ".</p>")
    if not mob and not pf and not put:
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
    sp_ok = sum(1 for r in ok if str(r["identification"].get("call", "")).startswith("species"))
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
                     f"{a.get('plasmids', 0)}"
                     + (f" (+{a['putative_plasmids']}?)" if a.get("putative_plasmids") else "")])
    body.append("<h2>Summary</h2>")
    body.append(_tbl(["Sample", "Barcode", "Species", "ANI", "Genome", "Contigs", "Depth",
                      "Compl./Cont. %", "ST · AMR · plasmids"], rows,
                     "Species from ANIb against the closest type-strain genome (NCBI type "
                     "material) where available, else ANI to the GTDB species "
                     "representative, else 16S rRNA (marked)."))
    if summ.get("warnings"):
        gen = list(dict.fromkeys(w for w in summ["warnings"] if not w.startswith("[")))
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
    if any(r.get("typestrain") for r in ok):
        iq = v.get("iqtree3") or v.get("iqtree2") or v.get("iqtree")
        body.append(
            "<p><b>Comparison with type strains</b> (after the DSMZ TYGS workflow and "
            "JSpeciesWS, rebuilt from open tools): all 16S copies were searched with vsearch "
            "against NCBI RefSeq 16S (type material). Type-strain genomes ('assembly from type "
            "material', one per species, most complete assembly first) of the candidate "
            "genera were downloaded from NCBI Datasets and screened with skani; for the "
            f"{esc(st.get('ts_max', 10))} closest, ANIb (Goris et al. 2007: 1020-nt "
            f"fragments, blastn{tv('blastn')} -task blastn -xdrop_gap_final 150 -dust no "
            "-evalue 1e-15, fragments with ≥ 70 % coverage and ≥ 30 % identity, mean "
            "pident; as implemented in pyANI), ANIm (nucmer"
            f"{tv('nucmer')} --mum, delta-filter -1) and TETRA (Teeling et al. 2004) were "
            "computed in both directions. Genome tree: FastME"
            f"{tv('fastme')} (BioNJ + SPR) on skani distances (1 − ANI/100); 16S tree: "
            f"MAFFT{tv('mafft')} (--auto), columns with ≥ 50 % bases, IQ-TREE"
            f"{' ' + esc(iq) if iq else ' (not installed)'} with ModelFinder and 1000 "
            "ultrafast bootstraps; both midpoint-rooted. Species verdict: ANIb ≥ 96 % (or "
            "≥ 95 % with TETRA ≥ 0.99) to a type strain = known species; 95–96 % = "
            "transition zone; < 95 % = potential novel species. dDDH (GBDP) is not computed "
            "because the DSMZ software is not distributed.</p>")
    if any(r.get("mapping") for r in ok):
        body.append(
            f"<p><b>Reference mapping</b>: filtered reads were mapped with minimap2"
            f"{tv('minimap2')} (-ax map-ont) to the closest genome found above (the most "
            f"complete assembly within 0.3 % ANI of the best), sorted and summarised with "
            f"samtools{tv('samtools')} (depth, flagstat; ≥ 1 kb at ≤ 2x = absent region). "
            "The polished assembly was aligned to the same genome with minimap2 -cx asm20 "
            "--cs; SNPs and indels were counted from primary alignments and isolate-specific "
            "regions are assembly stretches ≥ 5 kb without an alignment.</p>")
    dbrows = [[esc(k), f"<span class='mono'>{esc(Path(p).name if p else ('bundled default' if k == 'amrfinder_db' else 'not configured'))}</span>"]
              for k, p in dbs.items()]
    if dbrows:
        body.append(_tbl(["Database", "Version / file"], dbrows))
    body.append(
        "<h3>Interpretation thresholds</h3><ul>"
        "<li>Species: ANI (ANIb / ANIm / skani) ≥ 95–96 % to the type strain with "
        "≥ ~60 % of the genome aligned; TETRA ≥ 0.989–0.999 supports it (Richter & "
        "Rosselló-Móra 2009); 16S ≥ 98.65 % is compatible with, but not proof of, the "
        "same species (Kim et al. 2014). dDDH ≥ 70 % (TYGS) is the formal criterion when "
        "ANI falls in 95–96 %.</li>"
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
        "virulence_vfdb.csv, plasmids_mobsuite.csv, read_qc.csv, sample_map.csv, "
        "type_strain_comparison.csv, reference_mapping.csv</span></li>"
        "<li><span class='mono'>samples/&lt;sample&gt;/taxonomy/</span> — ANI tables, "
        "trees (Newick), 16S alignment; <span class='mono'>mapping/</span> — coverage, "
        "absent and isolate-specific regions, assembly-vs-reference PAF</li>"
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
