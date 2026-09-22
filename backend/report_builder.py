#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_builder.py — turn one or more NextGen-Amplicon result folders into a PDF report.

    python3 report_builder.py RESULTS_DIR [RESULTS_DIR ...] -o report.pdf

One folder produces a profile report; two or more produce a comparison, with every
figure and table gaining a column per run. Nothing about the runs needs to be
declared on the command line — the folder is read and the rest follows from it.

Why this exists as its own program rather than more R inside the pipeline: a report
is written after the fact, often over runs finished days apart, and it should be
possible to rebuild it without re-running anything. It reads only the CSV and JSON a
finished run leaves behind.

Four things the raw tables get wrong for a reader, which this corrects:

  1. SILVA 144 marks an unresolved genus by appending "--other"
     (Chryseobacterium--other). Printed verbatim that reads as a software fault, and
     it splits one genus across two rows. Collapsed here, with the count of affected
     reads kept and reported rather than hidden.

  2. Haemotropic mycoplasmas arrive under two superseded genus names,
     Eperythrozoon and Haemobartonella — both transferred into Mycoplasma in 2001.
     Left alone, one organism appears twice and neither row reaches its true
     abundance. Merged into a single line.

  3. Host mitochondrial DNA is co-amplified by "universal" 16S primers and can be a
     quarter of the reads. It sits inside the taxonomy table as an ordinary family,
     so a naive abundance table silently counts it as a bacterium.

  4. A well-known set of genera are routinely recovered from extraction kits and PCR
     reagents. They are flagged, never removed, and every figure says plainly that
     without a sequenced blank control the split is an estimate.

Requires: weasyprint (pip install weasyprint). Without it the HTML is still written
and the path reported, so the report can be printed from a browser instead.
"""

from __future__ import annotations

import argparse
import csv
import json
import html
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path

__version__ = "1.0.0"

# ── palette ────────────────────────────────────────────────────────────────────
INK, ACC, WARN, STOP = "#16242B", "#1E6E8C", "#96660A", "#A83232"
MUT, FAINT, RULE, KIT_C = "#4E6672", "#7D949E", "#D6E0E4", "#C9A227"
SERIES = ["#1E6E8C", "#8FA8B2", "#2E6B46", "#96660A", "#6B4E8C"]

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Genus labels meaning "we could not resolve this", in every spelling seen so far.
UNRESOLVED_RE = re.compile(r"--+(other|unclassified|uncultured|unknown)$", re.I)
BLANK_RE = re.compile(
    r"^(na|nan|null|unknown|unassigned|unclassified|uncultured|undetermined)$", re.I
)

# Both were folded into Mycoplasma in 2001; reference databases still carry them.
HAEMOPLASMA = {"Eperythrozoon", "Haemobartonella"}
HAEMOPLASMA_LABEL = "Mycoplasma (haemotropic)"

# Genera repeatedly recovered from DNA extraction kits and PCR reagents.
# Salter et al. 2014, BMC Biol 12:87, plus later additions. Override with --kit-list.
KITOME = {
    "Sphingomonas", "Alteriyabuuchia", "Ralstonia", "Methylobacterium", "Cutibacterium",
    "Propionibacterium", "Bradyrhizobium", "Acinetobacter", "Pseudomonas", "Burkholderia",
    "Delftia", "Herbaspirillum", "Novosphingobium", "Sphingosinithalassobacter",
    "Stenotrophomonas", "Chryseobacterium", "Brevundimonas", "Paracoccus", "Afipia",
    "Massilia", "Curvibacter", "Roseateles", "Pelomonas", "Diaphorobacter",
    "Pseudacidovorax", "Variovorax", "Comamonas", "Janthinobacterium", "Rhodococcus",
}


def esc(s) -> str:
    return html.escape(str(s))


def fmt(v) -> str:
    """Thousands-separated integer, tolerant of the strings CSV columns hand back."""
    try:
        return f"{int(round(float(v))):,}"
    except (TypeError, ValueError):
        return str(v) if v not in (None, "") else "\u2014"


# ══════════════════════════════════════════════════════════════════════════════
#  Reading a result folder
# ══════════════════════════════════════════════════════════════════════════════
class Run:
    """One result folder, normalized.

    Everything downstream reads this, so the differences between pipeline versions
    and reference databases are absorbed in one place.
    """

    def __init__(self, path: Path, kitome: set[str], label: str | None = None):
        self.path = Path(path)
        if not self.path.is_dir():
            raise FileNotFoundError(f"not a directory: {self.path}")
        self.kitome = kitome
        self.label = label or self.path.name

        self.params = self._json("run_params.json")
        self.summary = self._json("summary.json")
        self.qc = self._json("qc_report.json")

        self._load_tables()
        self._classify()

    # ── file helpers ──────────────────────────────────────────────────────────
    def _json(self, name) -> dict:
        p = self.path / name
        if not p.exists():
            return {}
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _rows(self, name):
        p = self.path / name
        if not p.exists():
            return []
        with p.open(newline="", encoding="utf-8-sig") as fh:
            return list(csv.DictReader(fh))

    # ── tables ────────────────────────────────────────────────────────────────
    def _load_tables(self):
        tax_rows = self._rows("taxonomy_table.csv")
        if not tax_rows:
            raise ValueError(f"{self.path.name}: taxonomy_table.csv missing or empty")
        key = list(tax_rows[0].keys())[0]  # unnamed first column = the ASV sequence
        self.tax = {r[key]: r for r in tax_rows}
        # SILVA 144 inserts a Kingdom rank below Domain, so the column set differs
        # between references. Take whatever ranks are actually present.
        self.rank_cols = [r for r in RANKS if r in tax_rows[0]]

        asv_rows = self._rows("asv_table.csv")
        if not asv_rows:
            raise ValueError(f"{self.path.name}: asv_table.csv missing or empty")
        cols = list(asv_rows[0].keys())
        seq_col = "sequence" if "sequence" in cols else cols[-1]
        self.samples = [c for c in cols if c != seq_col and c != ""]
        # Per-sample counts, and the total each ASV carries across all samples.
        self.per_sample = {s: 0.0 for s in self.samples}
        self.abund: dict[str, float] = {}
        for r in asv_rows:
            seq = r[seq_col]
            tot = 0.0
            for s in self.samples:
                try:
                    v = float(r[s] or 0)
                except (TypeError, ValueError):
                    v = 0.0
                self.per_sample[s] += v
                tot += v
            self.abund[seq] = tot
        self.total = sum(self.abund.values())

        self.tracking = self._rows("read_tracking.csv")
        self.alpha = self._rows("alpha_diversity.csv")
        self.rarefaction = self._rows("rarefaction.csv")
        self.orientation = self._json("orientation_report.json")

    # ── derived views ─────────────────────────────────────────────────────────
    def _clean_genus(self, g):
        """One label per genus: no blanks, no '--other', no superseded haemoplasma names."""
        if not g:
            return None
        g = str(g).strip()
        if not g or BLANK_RE.match(g):
            return None
        if g in HAEMOPLASMA:
            return HAEMOPLASMA_LABEL
        return UNRESOLVED_RE.sub("", g)

    def _classify(self):
        # Host organelle reads hide inside the taxonomy as an ordinary family.
        self.host_seqs = [
            s for s, r in self.tax.items()
            if any(v and ("Mitochondria" in v or "Chloroplast" in v) for v in r.values())
        ]
        top = "Domain" if "Domain" in next(iter(self.tax.values()), {}) else "Kingdom"
        self.euk_seqs = [
            s for s, r in self.tax.items()
            if s not in self.host_seqs and (r.get(top) or "") == "Eukaryota"
        ]
        host_set = set(self.host_seqs) | set(self.euk_seqs)
        self.bact_seqs = [s for s in self.tax if s not in host_set]

        self.host = sum(self.abund.get(s, 0) for s in self.host_seqs)
        self.euk = sum(self.abund.get(s, 0) for s in self.euk_seqs)
        self.bacteria = sum(self.abund.get(s, 0) for s in self.bact_seqs)

        self.genera = Counter()
        self.unresolved = Counter()   # reads that arrived with a '--other' suffix
        for s in self.bact_seqs:
            raw = self.tax[s].get("Genus")
            g = self._clean_genus(raw)
            self.genera[g or "(unassigned)"] += self.abund.get(s, 0)
            if raw and UNRESOLVED_RE.search(str(raw)):
                self.unresolved[g] += self.abund.get(s, 0)

        self.kit = sum(v for g, v in self.genera.items() if g in self.kitome)
        self.no_genus = self.genera.get("(unassigned)", 0.0)
        self.sample_derived = self.bacteria - self.kit - self.no_genus

        # The literal string "NA" survives the CSV round-trip and would otherwise be
        # printed as a phylum called NA.
        self.phyla = Counter()
        for s in self.bact_seqs:
            v = (self.tax[s].get("Phylum") or "").strip()
            self.phyla["(unassigned)" if not v or BLANK_RE.match(v) else v] += self.abund.get(s, 0)

        self.species = OrderedDict()
        for s, r in self.tax.items():
            v = r.get("Species")
            if v and not BLANK_RE.match(str(v).strip()):
                g = self._clean_genus(r.get("Genus")) or "?"
                self.species[f"{g} {v}"] = (self.abund.get(s, 0), len(s))
        self.species = OrderedDict(
            sorted(self.species.items(), key=lambda kv: -kv[1][0])
        )

    # ── reporting helpers ─────────────────────────────────────────────────────
    def rank_pct(self, rank) -> float | None:
        if rank not in self.rank_cols or not self.total:
            return None
        hit = [s for s, r in self.tax.items()
               if r.get(rank) and not BLANK_RE.match(str(r[rank]).strip())]
        return round(100 * sum(self.abund.get(s, 0) for s in hit) / self.total, 1)

    def top_genera(self, n=10):
        out = []
        for g, v in self.genera.most_common(n + 1):
            if g == "(unassigned)":
                continue
            out.append({
                "name": g, "reads": v,
                "pct": round(100 * v / self.bacteria, 1) if self.bacteria else 0.0,
                "kit": g in self.kitome,
                "unresolved": self.unresolved.get(g, 0),
            })
            if len(out) == n:
                break
        return out

    def top_asvs(self, n=30):
        """The n most abundant ASVs with their lineage, deepest rank first.

        ASV identifiers are positional (ASV1, ASV2, ...) in abundance order, which is
        what the pipeline's own outputs use; the sequence itself stays in asvs.fasta.
        """
        order = sorted(self.abund.items(), key=lambda kv: -kv[1])
        out = []
        for i, (seq, reads) in enumerate(order[:n], start=1):
            r = self.tax.get(seq, {})
            lineage, deepest = [], "\u2014"
            for rank in self.rank_cols:
                v = r.get(rank)
                if v and not BLANK_RE.match(str(v).strip()):
                    v = UNRESOLVED_RE.sub(" sp.", str(v))
                    lineage.append(v)
                    deepest = f"{rank[0]}: {v}"
            host = seq in self.host_seqs
            out.append({
                "id": f"ASV{i}", "len": len(seq), "reads": reads,
                "pct": 100 * reads / self.total if self.total else 0,
                "phylum": (r.get("Phylum") or "\u2014"),
                "genus": (self._clean_genus(r.get("Genus")) or "\u2014"),
                "deepest": deepest, "lineage": " > ".join(lineage) or "unassigned",
                "host": host,
            })
        return out

    def lengths(self) -> Counter:
        c = Counter()
        for s, v in self.abund.items():
            c[len(s)] += v
        return c

    @property
    def db_name(self) -> str:
        p = self.params.get("dbPath") or ""
        if p:
            return Path(p).name
        return self.params.get("taxDatabase") or "(pipeline default)"

    @property
    def trunc(self) -> str:
        f, r = self.params.get("truncLen_F"), self.params.get("truncLen_R")
        return f"{f} / {r}" if f and r else "—"

    @property
    def ceiling(self):
        f, r = self.params.get("truncLen_F"), self.params.get("truncLen_R")
        return f + r - 12 if f and r else None

    @property
    def n_asvs(self) -> int:
        return self.summary.get("n_asvs") or len(self.tax)


# ══════════════════════════════════════════════════════════════════════════════
#  Figures — inline SVG, so the PDF carries no external assets
# ══════════════════════════════════════════════════════════════════════════════
def svg_ranks(runs) -> str:
    """Share of reads assigned at each rank, one bar per run."""
    present = [r for r in RANKS if any(x.rank_pct(r) is not None for x in runs)]
    W, lab, plot = 520, 74, 380
    barh, gap = 7, 2.5
    rowh = max(22, len(runs) * (barh + gap) + 11)
    H = len(present) * rowh + 26
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for gx in range(0, 101, 25):
        x = lab + plot * gx / 100
        p.append(f'<line x1="{x:.1f}" y1="10" x2="{x:.1f}" y2="{len(present)*rowh+10:.1f}" '
                 f'stroke="{RULE}" stroke-width="1"/>')
        p.append(f'<text x="{x:.1f}" y="{len(present)*rowh+22}" font-size="7.5" fill="{FAINT}" '
                 f'text-anchor="middle" font-family="DejaVu Sans Mono">{gx}%</text>')
    for i, rank in enumerate(present):
        y = 10 + i * rowh
        p.append(f'<text x="{lab-8}" y="{y+11}" font-size="9" fill="{INK}" text-anchor="end">{rank}</text>')
        for j, run in enumerate(runs):
            v = run.rank_pct(rank)
            if v is None:
                continue
            bw = plot * v / 100
            yy = y + 3 + j * (barh + gap)
            col = SERIES[j % len(SERIES)]
            p.append(f'<rect x="{lab}" y="{yy:.1f}" width="{max(bw,0.6):.1f}" height="{barh}" fill="{col}"/>')
            p.append(f'<text x="{lab+bw+4:.1f}" y="{yy+barh-0.5:.1f}" font-size="7.5" fill="{col}" '
                     f'font-family="DejaVu Sans Mono">{v}%</text>')
    p.append("</svg>")
    return "".join(p)


def svg_composition(runs) -> str:
    """Host / no-genus / reagent / sample-derived, as a share of each run's reads."""
    segs = [("host", "Host / non-bacterial", "#B4C5CC"),
            ("nogen", "Bacteria, no genus", FAINT),
            ("kit", "Reagent-associated", KIT_C),
            ("samp", "Sample-derived", ACC)]
    W, x0, bw, barh, pitch = 520, 6, 452, 28, 56
    H = len(runs) * pitch + 22
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for i, run in enumerate(runs):
        y = 18 + i * pitch
        tot = run.total or 1
        vals = {"host": run.host + run.euk, "nogen": run.no_genus,
                "kit": run.kit, "samp": max(run.sample_derived, 0)}
        p.append(f'<text x="{x0}" y="{y-5}" font-size="8.5" fill="{INK}" '
                 f'font-family="DejaVu Sans Mono">{esc(run.label)}</text>')
        x = x0
        for k, _l, c in segs:
            w = bw * vals[k] / tot
            p.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{barh}" fill="{c}"/>')
            if w > 32:
                p.append(f'<text x="{x+w/2:.1f}" y="{y+18}" font-size="8" fill="#fff" '
                         f'text-anchor="middle" font-family="DejaVu Sans Mono">{fmt(vals[k])}</text>')
            x += w
        p.append(f'<text x="{x0+bw+6}" y="{y+18}" font-size="7.5" fill="{FAINT}" '
                 f'font-family="DejaVu Sans Mono">{fmt(tot)}</text>')
    lx = x0
    for _k, l, c in segs:
        p.append(f'<rect x="{lx}" y="{H-13}" width="9" height="9" fill="{c}"/>')
        p.append(f'<text x="{lx+12}" y="{H-5}" font-size="7.5" fill="{MUT}">{l}</text>')
        lx += 12 + len(l) * 4.2 + 14
    p.append("</svg>")
    return "".join(p)


def svg_lengths(run) -> str:
    """ASV length histogram against the merge ceiling.

    The ceiling is truncLen_F + truncLen_R - minOverlap: the longest insert mergePairs()
    can join. When the longest surviving ASV sits on that number, longer amplicons were
    discarded with no error anywhere, so drawing the two together is the whole check.
    """
    h = run.lengths()
    if not h:
        return ""
    binned = Counter()
    for L, v in h.items():
        binned[(L // 10) * 10] += v
    ceil = run.ceiling
    lo = min(binned) - 5
    hi = max(max(binned) + 15, (ceil + 10) if ceil else 0)
    W, H, x0, y0, pw, ph = 520, 156, 34, 20, 400, 100
    mx = max(binned.values())
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for gy in (0, 0.5, 1.0):
        y = y0 + ph - ph * gy
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}" stroke="{RULE}"/>')
        p.append(f'<text x="{x0-5}" y="{y+3:.1f}" font-size="7" fill="{FAINT}" text-anchor="end" '
                 f'font-family="DejaVu Sans Mono">{fmt(mx*gy)}</text>')
    span = max(hi - lo, 1)
    for L, v in sorted(binned.items()):
        x = x0 + pw * (L - lo) / span
        w = max(pw * 10 / span - 1, 1)
        bh = ph * v / mx
        p.append(f'<rect x="{x:.1f}" y="{y0+ph-bh:.1f}" width="{w:.1f}" height="{bh:.1f}" '
                 f'fill="{ACC}" opacity="0.85"/>')
    if ceil:
        cx = x0 + pw * (ceil - lo) / span
        p.append(f'<line x1="{cx:.1f}" y1="{y0-1}" x2="{cx:.1f}" y2="{y0+ph}" stroke="{STOP}" '
                 f'stroke-width="1.2" stroke-dasharray="3 2"/>')
        # above the plot, never inside it — at plot height this label lands on the bars
        p.append(f'<text x="{cx+3:.1f}" y="{y0-3}" font-size="7.5" fill="{STOP}" '
                 f'text-anchor="start" font-family="DejaVu Sans Mono">merge ceiling {ceil} bp</text>')
    step = 50 if span > 150 else 20
    t = int((lo + step) // step * step)
    while t < hi:
        x = x0 + pw * (t - lo) / span
        p.append(f'<text x="{x:.1f}" y="{y0+ph+13}" font-size="7" fill="{FAINT}" '
                 f'text-anchor="middle" font-family="DejaVu Sans Mono">{t}</text>')
        t += step
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-4}" font-size="8" fill="{MUT}" '
             f'text-anchor="middle">ASV length (bp)</text>')
    p.append("</svg>")
    return "".join(p)


def svg_genera(runs, n=10) -> str:
    """Top genera per run: name and value on one line, bar beneath, one shared axis.

    An earlier layout interleaved bar, value and the next genus name inside one narrow
    row; at print size the value read as belonging to the label below it.
    """
    cols = [r.top_genera(n) for r in runs]
    mx = max((g["pct"] for c in cols for g in c), default=1) or 1
    ncol = len(runs)
    gap = 24
    colw = (520 - gap * (ncol - 1)) / ncol
    rowh = 26
    H = n * rowh + 34
    p = [f'<svg viewBox="0 0 520 {H}" width="100%" style="max-width:520px">']
    for ci, (run, top) in enumerate(zip(runs, cols)):
        ox = ci * (colw + gap)
        col = SERIES[ci % len(SERIES)]
        p.append(f'<text x="{ox:.1f}" y="8" font-size="8" fill="{INK}" '
                 f'font-family="DejaVu Sans Mono">{esc(run.label)}</text>')
        p.append(f'<line x1="{ox:.1f}" y1="12" x2="{ox+colw:.1f}" y2="12" stroke="{RULE}"/>')
        for i, g in enumerate(top):
            y = 22 + i * rowh
            nm = g["name"]
            if len(nm) > int(colw / 4.6):
                nm = nm[: int(colw / 4.6) - 1] + "…"
            p.append(f'<text x="{ox:.1f}" y="{y+6}" font-size="8" fill="{INK}">'
                     f'{esc(nm)}{" *" if g["kit"] else ""}</text>')
            p.append(f'<text x="{ox+colw:.1f}" y="{y+6}" font-size="7.5" fill="{MUT}" '
                     f'text-anchor="end" font-family="DejaVu Sans Mono">{g["pct"]}%</text>')
            bw = colw * g["pct"] / mx
            p.append(f'<rect x="{ox:.1f}" y="{y+9}" width="{bw:.1f}" height="5.5" '
                     f'fill="{KIT_C if g["kit"] else col}"/>')
    p.append(f'<text x="0" y="{H-4}" font-size="7" fill="{FAINT}">'
             f'* reagent-associated genus &#183; bars share one axis across all columns</text>')
    p.append("</svg>")
    return "".join(p)


# ══════════════════════════════════════════════════════════════════════════════
#  Document
# ══════════════════════════════════════════════════════════════════════════════
def _css(app: str, company: str, sample: str) -> str:
    foot = f"{app}  \\00B7  {company}  \\00B7  "
    return f"""
@page {{
  size: A4; margin: 17mm 16mm 20mm 16mm;
  @bottom-right {{ content: "{foot}" counter(page) " / " counter(pages);
    font-family: "DejaVu Sans Mono", monospace; font-size: 6.6pt; color: {FAINT}; }}
  @bottom-left {{ content: "{sample}";
    font-family: "DejaVu Sans Mono", monospace; font-size: 6.6pt; color: {FAINT}; }}
}}
* {{ box-sizing: border-box; }}
body {{ font-family: Carlito, "DejaVu Sans", sans-serif; font-size: 9.4pt;
        line-height: 1.55; color: {INK}; margin: 0; }}
h1 {{ font-family: "Bitstream Charter", "DejaVu Serif", serif; font-size: 22pt;
      font-weight: 700; margin: 0 0 2mm; line-height: 1.12; letter-spacing: -.01em; }}
h2 {{ font-family: "Bitstream Charter", "DejaVu Serif", serif; font-size: 13pt;
      font-weight: 700; margin: 9mm 0 1.5mm; padding-bottom: 1mm;
      border-bottom: 1px solid {RULE}; break-after: avoid; }}
h3 {{ font-family: "Bitstream Charter", "DejaVu Serif", serif; font-size: 10.5pt;
      font-weight: 700; margin: 5mm 0 1mm; break-after: avoid; }}
p {{ margin: 0 0 2.6mm; }}
.eyebrow {{ font-family: "DejaVu Sans Mono", monospace; font-size: 7pt;
            letter-spacing: .14em; text-transform: uppercase; color: {ACC}; margin-bottom: 2mm; }}
.sub {{ color: {MUT}; font-size: 10pt; margin-bottom: 5mm; }}
.dim {{ color: {FAINT}; }} .hl {{ color: {ACC}; font-weight: 700; }}
.mono {{ font-family: "DejaVu Sans Mono", monospace; }}
i {{ font-family: "Bitstream Charter", "DejaVu Serif", serif; font-style: italic; }}
.kpis {{ display: flex; gap: 3mm; margin: 4mm 0 5mm; }}
.kpi {{ flex: 1; border: 1px solid {RULE}; border-top: 2px solid {ACC}; padding: 2.5mm 3mm; }}
.kpi .v {{ font-family: "DejaVu Sans Mono", monospace; font-size: 13pt; font-weight: 700; line-height: 1.1; }}
.kpi .l {{ font-size: 7.2pt; color: {MUT}; margin-top: .8mm; line-height: 1.3; }}
table {{ width: 100%; border-collapse: collapse; font-size: 8.6pt; margin: 2mm 0; }}
thead th {{ font-family: "DejaVu Sans Mono", monospace; font-size: 6.8pt; letter-spacing: .09em;
            text-transform: uppercase; color: {FAINT}; text-align: right; font-weight: 400;
            padding: 0 2mm 1.2mm 0; border-bottom: 1px solid {INK}; }}
thead th:first-child {{ text-align: left; }}
tbody td {{ padding: 1.3mm 2mm 1.3mm 0; border-bottom: .5px solid {RULE}; text-align: right;
            font-variant-numeric: tabular-nums; }}
tbody td:first-child {{ text-align: left; }}
tbody tr {{ break-inside: avoid; }}
caption {{ caption-side: bottom; text-align: left; font-size: 7.6pt; color: {FAINT}; padding-top: 1.2mm; }}
figure {{ margin: 3mm 0 2mm; break-inside: avoid; }}
figcaption {{ font-size: 7.8pt; color: {MUT}; margin-top: 1.5mm; }}
.note {{ border-left: 2.5px solid {WARN}; background: #FAF6EC; padding: 2.5mm 3mm;
         margin: 3mm 0; font-size: 8.8pt; break-inside: avoid; }}
.note.acc {{ border-left-color: {ACC}; background: #EDF4F7; }}
.note b {{ color: {INK}; }}
ul {{ margin: 1mm 0 3mm; padding-left: 5mm; }} li {{ margin-bottom: 1.4mm; }}
.pbreak {{ break-before: page; }}
"""


def _tbl(head, rows, caption=None) -> str:
    h = "".join(f"<th>{c}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    cap = f"<caption>{caption}</caption>" if caption else ""
    return f"<table><thead><tr>{h}</tr></thead><tbody>{b}</tbody>{cap}</table>"


def build_profile_html(run, title=None, app="NextGen-Amplicon", company="",
                       subtitle=None, top_asv=30) -> str:
    """Result report for a single run: what was sequenced, whether it worked, and what
    was in it. This is the routine deliverable; the comparison layout below is for the
    less common case of holding two analyses against each other."""
    sample_label = ", ".join(run.samples[:3]) + ("\u2026" if len(run.samples) > 3 else "")
    marker = run.summary.get("marker") or run.params.get("marker") or "16S"
    title = title or f"{sample_label} \u2014 {marker} Amplicon Analysis"
    subtitle = subtitle or ("Sequence variants, community composition and quality control "
                            "for this sequencing run.")
    L = run.lengths()
    inp = 0.0
    for t in run.tracking:
        try:
            inp += float(t.get("input") or 0)
        except (TypeError, ValueError):
            pass
    out = [f'<div class="eyebrow">{esc(run.label)} &middot; {esc(marker)} &middot; '
           f'{len(run.samples)} sample{"s" if len(run.samples)!=1 else ""}</div>',
           f"<h1>{esc(title)}</h1>", f'<p class="sub">{esc(subtitle)}</p>']
    kpis = [(fmt(inp) if inp else "\u2014", "read pairs sequenced"),
            (fmt(run.total), "reads analysed"),
            (str(run.n_asvs), "sequence variants"),
            (f"{run.rank_pct('Genus') or 0}%", "reads with a genus")]
    out.append('<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>'
        for v, l in kpis) + "</div>")

    # ── 1. configuration ──────────────────────────────────────────────────────
    out.append("<h2>1. Run configuration</h2>")
    cfg = [("Sample(s)", sample_label), ("Marker", marker),
           ("Reference database", run.db_name),
           ("Forward primer", run.params.get("primer_f") or "\u2014"),
           ("Reverse primer", run.params.get("primer_r") or "\u2014"),
           ("truncLen F / R", run.trunc),
           ("maxEE F / R", f'{run.params.get("maxEE_F")} / {run.params.get("maxEE_R")}'),
           ("Minimum bootstrap", run.params.get("minBoot")),
           ("Chimera method", run.params.get("chimeraMethod"))]
    out.append(_tbl(["Parameter", "Value"],
                    [[k, esc(v if v is not None else "\u2014")] for k, v in cfg]))

    # ── 2. quality control ────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">2. Quality control</h2>')
    out.append("<h3>Read tracking</h3>")
    if run.tracking:
        cols = [c for c in ("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
                if c in run.tracking[0]]
        rows = []
        for t in run.tracking:
            base = float(t.get("input") or 0) or 1
            rows.append([esc(next(iter(t.values())))] + [fmt(t.get(c) or 0) for c in cols]
                        + [f'{100*float(t.get("nonchim") or 0)/base:.1f}%'])
        out.append(_tbl(["Sample"] + [c.capitalize() for c in cols] + ["Retained"], rows))
    ff = svg_readflow(run)
    if ff:
        out.append("<figure>" + ff + "<figcaption>Reads surviving each stage, as a share of "
                   "input. The right-hand figure is what that stage cost. Merging and chimera "
                   "removal are the two stages that normally dominate the loss."
                   "</figcaption></figure>")

    o = (run.orientation or {}).get("checked") or {}
    if o:
        out.append("<h3>Read orientation</h3>")
        out.append(_tbl(["Check", "Pairs", "Share"], [
            ["Already in the expected orientation", fmt(o.get("already_correct", 0)),
             f'{o.get("pct_already_correct", 0)}%'],
            ["Reversed (R1 and R2 swapped)", fmt(o.get("flipped", 0)),
             f'{o.get("pct_flipped", 0)}%'],
            ["Ambiguous", fmt(o.get("ambiguous", 0)), f'{o.get("pct_ambiguous", 0)}%'],
        ]))
        if run.orientation.get("applied"):
            out.append('<div class="note acc"><p style="margin:0"><b>Orientation was repaired '
                       "before analysis.</b> Read pairs in the reverse orientation cannot be "
                       "merged as they stand, so without this step they would have been lost "
                       "silently at the merge stage.</p></div>")

    if L:
        out.append("<h3>Amplicon length</h3>")
        out.append("<figure>" + svg_lengths(run) + "<figcaption>"
                   + (f"Distribution of ASV lengths. The dashed line is the merge ceiling "
                      f"(truncLen_F + truncLen_R &minus; 12 = {run.ceiling} bp), the longest "
                      f"fragment that can be assembled from a read pair. "
                      + ceiling_verdict(run)[1]
                      if run.ceiling else "Distribution of ASV lengths.")
                   + "</figcaption></figure>")

        if ceiling_verdict(run)[0]:
            out.append('<div class="note"><p style="margin:0"><b>The truncation settings are '
                       "cutting the amplicon short.</b> Every figure in this report is computed "
                       "over the fragments that survived merging, so a clipped window biases the "
                       "whole profile toward whatever happens to be short \u2014 host organelle "
                       "amplicons in particular. Treat the composition below as provisional until "
                       "the run is repeated with a wider window.</p></div>")

    qc = (run.qc or {}).get("offtarget") or {}
    if qc:
        out.append("<h3>Off-target reads</h3>")
        out.append(_tbl(["Category", "Reads", "Share"], [
            ["Host organelle (mitochondria / chloroplast)",
             fmt(qc.get("organelle_reads", 0)), f'{qc.get("organelle_pct", 0)}%'],
            ["Unassigned at domain level",
             fmt(qc.get("unassigned_reads", 0)), f'{qc.get("unassigned_pct", 0)}%'],
            ["Total off-target", fmt(qc.get("reads", 0)), f'{qc.get("pct", 0)}%'],
        ], "&ldquo;Universal&rdquo; 16S primers co-amplify host mitochondrial DNA, which is "
           "normal for host-associated and ectoparasite samples."))

    # ── 3. composition ────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">3. Community composition</h2>')
    out.append(f"<p>Shares below are of the {fmt(run.bacteria)} bacterial reads, with host "
               "organelle reads set aside.</p>")
    ph = [(n, v) for n, v in run.phyla.most_common(8) if n != "(unassigned)"]
    if ph:
        out.append("<h3>Phylum</h3>")
        out.append("<figure>" + svg_bars(ph, run.bacteria) + "</figure>")
    tg = run.top_genera(15)
    if tg:
        out.append("<h3>Genus</h3>")
        out.append("<figure>" + svg_bars([(g["name"], g["reads"]) for g in tg], run.bacteria,
                                         kit={g["name"] for g in tg if g["kit"]})
                   + "<figcaption>* reagent-associated genus. Unresolved sub-genus labels are "
                     "counted under the parent genus.</figcaption></figure>")
        out.append(_tbl(["Genus", "Reads", "% bacteria", "Note"],
                        [[f'<i>{esc(g["name"])}</i>', fmt(g["reads"]), f'{g["pct"]}%',
                          f'<span class="dim">'
                          + "; ".join(([f'reagent-associated'] if g["kit"] else [])
                                      + ([f'{fmt(g["unresolved"])} unresolved below genus']
                                         if g["unresolved"] else []))
                          + "</span>"] for g in tg]))
    out.append('<div class="note"><p style="margin:0"><b>Reagent-associated genera are flagged, '
               "not removed.</b> They are assigned from the published list of genera commonly "
               "recovered from DNA extraction kits and PCR reagents. Unless a blank extraction "
               "control was sequenced alongside this sample, the split is an estimate rather "
               "than a measurement.</p></div>")

    # ── 4. top ASVs ───────────────────────────────────────────────────────────
    out.append(f'<h2 class="pbreak">4. Most abundant sequence variants</h2>')
    out.append(f"<p>The {top_asv} most abundant of {run.n_asvs} ASVs. Identifiers are positional "
               "in abundance order; the sequences themselves are in "
               "<span class=\"mono\">asvs.fasta</span>, and the full table is in "
               "<span class=\"mono\">asv_table.csv</span>.</p>")
    rows = []
    for a in run.top_asvs(top_asv):
        rows.append([f'<span class="mono">{a["id"]}</span>',
                     f'{a["len"]}',
                     fmt(a["reads"]), f'{a["pct"]:.2f}%',
                     (f'<span class="dim">host organelle</span>' if a["host"]
                      else f'<i>{esc(a["genus"])}</i>' if a["genus"] != "\u2014"
                      else f'<span class="dim">{esc(a["deepest"])}</span>'),
                     f'<span style="font-size:7.4pt" class="dim">{esc(a["lineage"][:78])}</span>'])
    out.append(_tbl(["ASV", "bp", "Reads", "Share", "Genus", "Lineage"], rows))

    # ── 5. diversity ──────────────────────────────────────────────────────────
    if run.alpha or run.rarefaction:
        out.append('<h2 class="pbreak">5. Diversity</h2>')
    if run.alpha:
        cols = list(run.alpha[0].keys())
        out.append(_tbl(cols, [[esc(r.get(c, "")) for c in cols] for r in run.alpha],
                        "Observed = ASVs seen. Chao1 estimates how many there would be at "
                        "unlimited depth. Shannon and Simpson weight richness by evenness; "
                        "Simpson approaches 1 when no single organism dominates."))
    rr = svg_rarefaction(run)
    if rr:
        out.append("<figure>" + rr + "<figcaption>Rarefaction. A curve that has flattened means "
                   "the sequencing depth was enough to see what is present; one still rising at "
                   "the right-hand edge means deeper sequencing would recover more."
                   "</figcaption></figure>")

    # ── 6. species ────────────────────────────────────────────────────────────
    out.append("<h2>6. Species-level assignment</h2>")
    if run.species:
        out.append(_tbl(["Species", "Reads", "Share", "ASV length"],
                        [[f"<i>{esc(k)}</i>", fmt(v[0]),
                          f"{100*v[0]/run.total:.2f}%", f"{v[1]} bp"]
                         for k, v in run.species.items()],
                        f'{fmt(sum(v[0] for v in run.species.values()))} of {fmt(run.total)} '
                        f'reads \u2014 {run.rank_pct("Species")}% of the analysed data.'))
    else:
        out.append("<p>No read reached a species-level assignment.</p>")
    out.append("<p>Species are assigned by exact 100% match of the whole fragment against a "
               "species-annotated reference, not by a confidence score. A low rate is the "
               "expected outcome for a short amplicon: the V3&ndash;V4 region covers roughly "
               "440 bp of a 1,540 bp gene, and uncultured environmental bacteria frequently "
               "have no species-annotated entry to match against. Full-length 16S "
               "(27F/1492R on PacBio or ONT) is the appropriate marker where species-level "
               "identification is required.</p>")

    # ── 7. notes ──────────────────────────────────────────────────────────────
    out.append("<h2>7. Notes on reading this report</h2>")
    notes = []
    if run.total:
        notes.append(f"{100*(run.host+run.euk)/run.total:.0f}% of analysed reads are host or "
                     "other non-bacterial material co-amplified by the primers, and are "
                     "excluded from the composition figures.")
    if run.bacteria:
        notes.append(f"About {100*run.kit/run.bacteria:.0f}% of the bacterial fraction matches "
                     "the reagent contaminant profile. Relative abundances are more honestly "
                     "computed over the remaining sample-derived reads.")
    if run.unresolved:
        notes.append(f"{fmt(sum(run.unresolved.values()))} reads carry a genus the reference "
                     "could not resolve further; they are counted under the parent genus here "
                     "and appear as separate rows in the raw taxonomy table.")
    if run.genera.get(HAEMOPLASMA_LABEL):
        notes.append("Haemotropic mycoplasmas are reported as one group: reference databases "
                     "still split them between the superseded genus names Eperythrozoon and "
                     "Haemobartonella, both transferred into Mycoplasma in 2001.")
    notes.append("A blank extraction control sequenced alongside the sample would turn the "
                 "reagent estimate into a measurement.")
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")

    css = _css(app, company, f"{sample_label} \u00b7 {marker}")
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            + "".join(out) + "</body></html>")


def build_html(runs, title=None, app="NextGen-Amplicon", company="", subtitle=None) -> str:
    multi = len(runs) > 1
    ref = runs[0]
    sample_label = ", ".join(ref.samples[:3]) + ("…" if len(ref.samples) > 3 else "")
    marker = ref.summary.get("marker") or ref.params.get("marker") or "16S"
    title = title or (f"{sample_label} — Reference Comparison" if multi
                      else f"{sample_label} — Community Profile")
    subtitle = subtitle or (
        "Two or more analyses of the same sequencing run, compared side by side. "
        "Where the runs share a pipeline configuration, every difference below is "
        "attributable to what changed between them."
        if multi else
        "Community profile from a single amplicon sequencing run."
    )
    names = [r.label for r in runs]

    # ── header + KPIs ─────────────────────────────────────────────────────────
    inp = 0
    for t in ref.tracking:
        try:
            inp += float(t.get("input") or 0)
        except (TypeError, ValueError):
            pass
    L = ref.lengths()
    kpis = [(fmt(inp) if inp else "—", "read pairs sequenced"),
            (fmt(ref.total), "reads analysed"),
            (str(ref.n_asvs), "sequence variants"),
            (f"{max(L)} bp" if L else "—",
             f"longest ASV (ceiling {ref.ceiling})" if ref.ceiling else "longest ASV")]
    out = [f'<div class="eyebrow">{esc(sample_label)} &middot; {esc(marker)} &middot; '
           f'{len(ref.samples)} sample{"s" if len(ref.samples)!=1 else ""}</div>',
           f"<h1>{esc(title)}</h1>", f'<p class="sub">{esc(subtitle)}</p>',
           '<div class="kpis">' + "".join(
               f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>'
               for v, l in kpis) + "</div>"]

    # ── 1. method ─────────────────────────────────────────────────────────────
    out.append("<h2>1. Run configuration</h2>")
    keys = [("Reference database", "db_name"), ("truncLen F / R", "trunc")]
    prm = [("maxEE F / R", lambda r: f'{r.params.get("maxEE_F")} / {r.params.get("maxEE_R")}'),
           ("Minimum bootstrap", lambda r: r.params.get("minBoot")),
           ("Chimera method", lambda r: r.params.get("chimeraMethod")),
           ("Sequencer", lambda r: r.params.get("sequencerType"))]
    rows = [[k] + [esc(getattr(r, a)) for r in runs] for k, a in keys]
    for k, fn in prm:
        vals = [fn(r) for r in runs]
        if any(v not in (None, "", "None") for v in vals):
            rows.append([k] + [esc(v if v is not None else "—") for v in vals])
    out.append(_tbl(["Parameter"] + [esc(n) for n in names], rows,
                    "Rows that differ between columns are the variables under comparison."
                    if multi else None))

    if ref.tracking:
        cols = [c for c in ("input", "filtered", "merged", "nonchim") if c in ref.tracking[0]]
        trows = []
        for t in ref.tracking:
            name = next(iter(t.values()))
            base = float(t.get("input") or 0) or 1
            trows.append([esc(name)] + [fmt(t.get(c) or 0) for c in cols]
                         + [f'{100*float(t.get("nonchim") or 0)/base:.1f}%'])
        out.append("<h3>Read tracking</h3>")
        out.append(_tbl(["Sample"] + [c.capitalize() for c in cols] + ["Retained"], trows,
                        "Identical across runs that share a pipeline configuration — "
                        "the reference database is applied after this point." if multi else None))

    if L:
        out.append("<figure>" + svg_lengths(ref) + "<figcaption>"
                   + ("ASV length distribution. The dashed line marks the merge ceiling "
                      f"(truncLen_F + truncLen_R &minus; 12 = {ref.ceiling} bp), the longest insert "
                      f"that can be joined. The longest ASV here is {max(L)} bp"
                      + (", so the window is not clipping the amplicon."
                         if ref.ceiling and max(L) < ref.ceiling - 2
                         else " — on the ceiling, so longer amplicons were discarded.")
                      if ref.ceiling else "ASV length distribution.")
                   + "</figcaption></figure>")

    # ── 2. depth ──────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">2. Taxonomic depth</h2>')
    out.append(f"<p>Share of the {fmt(ref.total)} analysed reads carrying an assignment at "
               f"each rank.</p>")
    out.append("<figure>" + svg_ranks(runs) + "<figcaption>"
               + " &middot; ".join(f'<span style="color:{SERIES[i%len(SERIES)]}">&#9632;</span> {esc(n)}'
                                   for i, n in enumerate(names)) + "</figcaption></figure>")
    drows = []
    for rank in RANKS:
        vals = [r.rank_pct(rank) for r in runs]
        if all(v is None for v in vals):
            continue
        row = [rank] + [f"{v}%" if v is not None else "—" for v in vals]
        if multi and vals[0] is not None and vals[-1] is not None:
            row.append(f"{vals[-1]-vals[0]:+.1f}")
        drows.append(row)
    out.append(_tbl(["Rank"] + [esc(n) for n in names] + (["Difference"] if multi else []), drows))

    # ── 3. composition ────────────────────────────────────────────────────────
    out.append("<h2>3. Composition of the analysed reads</h2>")
    out.append("<figure>" + svg_composition(runs) + "<figcaption>"
               "Host organelle reads sit inside the taxonomy table as an ordinary family; "
               "counting them as bacteria would overstate every abundance below."
               "</figcaption></figure>")
    out.append('<div class="note"><p style="margin:0"><b>The reagent-associated fraction is an '
               "estimate, not a measurement.</b> It is assigned from the published list of genera "
               "commonly recovered from DNA extraction kits and PCR reagents. Unless a blank "
               "extraction control was sequenced alongside this sample, the split cannot be "
               "verified.</p></div>")

    # ── 4. genera ─────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">4. Genus-level profile</h2>')
    out.append("<p>Most abundant genera as a share of each run's bacterial reads, host "
               "material excluded.</p>")
    out.append("<figure>" + svg_genera(runs) + "<figcaption>Bars share one axis. Gold marks "
               "reagent-associated genera.</figcaption></figure>")
    seen, order = set(), []
    for g in runs[-1].top_genera(16):
        order.append(g["name"]); seen.add(g["name"])
    grows = []
    for nm in order:
        row = [f"<i>{esc(nm)}</i>" if " " not in nm or nm.startswith("Candidatus") else esc(nm)]
        for r in runs:
            v = r.genera.get(nm, 0)
            row.append(f"{fmt(v)} &middot; {100*v/r.bacteria:.1f}%" if v and r.bacteria
                       else '<span class="dim">not named</span>')
        notes = []
        if nm in runs[-1].kitome:
            notes.append("reagent")
        if runs[-1].unresolved.get(nm):
            notes.append(f'{fmt(runs[-1].unresolved[nm])} unresolved below genus')
        row.append(f'<span class="dim">{"; ".join(notes)}</span>')
        grows.append(row)
    out.append(_tbl(["Genus"] + [esc(n) for n in names] + ["Note"], grows,
                    "&ldquo;Not named&rdquo; means that reference assigns those reads to a "
                    "different genus." if multi else None))

    if runs[-1].unresolved:
        tot_u = sum(runs[-1].unresolved.values())
        out.append(f'<div class="note"><p style="margin:0"><b>{fmt(tot_u)} reads carry a genus the '
                   "reference could not resolve further.</b> SILVA 144 marks these by appending "
                   "<span class=\"mono\">--other</span> to the genus name. They are counted under "
                   "the parent genus in this report; in the raw "
                   "<span class=\"mono\">taxonomy_table.csv</span> they appear as separate rows."
                   "</p></div>")

    haemo = [r.genera.get(HAEMOPLASMA_LABEL, 0) for r in runs]
    if any(haemo):
        r = runs[-1]
        out.append(f'<div class="note acc"><p style="margin:0"><b>Haemotropic <i>Mycoplasma</i> '
                   f"spp. &mdash; {fmt(haemo[-1])} reads, "
                   f"{100*haemo[-1]/r.bacteria:.1f}% of bacterial reads.</b> Reference databases "
                   "split these across the genus labels <i>Eperythrozoon</i> and "
                   "<i>Haemobartonella</i>; members of both were transferred into <i>Mycoplasma</i> "
                   "in 2001, so they are one group under superseded names and are merged here."
                   "</p></div>")

    # ── 5. phyla ──────────────────────────────────────────────────────────────
    out.append("<h2>5. Phylum-level summary</h2>")
    for r in runs:
        rows = [[esc(n), fmt(v), f"{100*v/r.bacteria:.1f}%"] for n, v in r.phyla.most_common(8)]
        out.append(f"<h3>{esc(r.label)}</h3>" + _tbl(["Phylum", "Reads", "%"], rows))
    if multi:
        out.append('<p class="dim" style="font-size:8pt">Reference releases differ in higher-rank '
                   "naming (SILVA 144 follows LPSN, which splits Bacillota into Bacillota and "
                   "Clostridiota), so these tables are not row-comparable below their totals.</p>")

    # ── 6. species ────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">6. Species-level assignment</h2>')
    out.append("<p>Species are assigned by exact 100% match of the full fragment against a "
               "species-annotated reference, not by a confidence score.</p>")
    allsp = OrderedDict()
    for r in runs:
        for k, v in r.species.items():
            allsp.setdefault(k, {})[r.label] = v
    if allsp:
        rows = []
        for k, per in allsp.items():
            first = next(iter(per.values()))
            row = [f"<i>{esc(k)}</i>"]
            for r in runs:
                row.append(fmt(per[r.label][0]) if r.label in per else '<span class="dim">—</span>')
            row.append(f"{first[1]} bp")
            rows.append(row)
        out.append(_tbl(["Species"] + [esc(n) for n in names] + ["ASV length"], rows,
                        f"{fmt(sum(v[0] for v in runs[-1].species.values()))} of "
                        f"{fmt(runs[-1].total)} reads — "
                        f"{runs[-1].rank_pct('Species')}% of the analysed data."))
    else:
        out.append("<p>No read reached a species-level assignment in any run.</p>")
    out.append("<p>A low species rate is the expected outcome for a short amplicon. The "
               "V3&ndash;V4 region covers roughly 440 bp of a 1,540 bp gene, and uncultured "
               "environmental bacteria frequently have no species-annotated entry in any 16S "
               "reference to match against. The figure reflects the marker and the sample rather "
               "than the analysis; a different database or different settings do not change it. "
               "Full-length 16S (27F/1492R on PacBio or ONT) is the appropriate marker where "
               "species-level identification is required.</p>")

    # ── 7. summary ────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">7. Summary</h2>')
    srows = [["Reads analysed"] + [fmt(r.total) for r in runs],
             ["Sequence variants"] + [str(r.n_asvs) for r in runs]]
    for rank in ("Family", "Genus", "Species"):
        vals = [r.rank_pct(rank) for r in runs]
        if any(v is not None for v in vals):
            srows.append([f"Reads placed to {rank.lower()}"] +
                         [f"{v}%" if v is not None else "—" for v in vals])
    srows += [["Host / non-bacterial reads"] + [fmt(r.host + r.euk) for r in runs],
              ["Reagent-associated (estimate)"] + [fmt(r.kit) for r in runs],
              ["Sample-derived bacteria"] + [fmt(max(r.sample_derived, 0)) for r in runs]]
    out.append(_tbl([""] + [esc(n) for n in names], srows))

    bullets = []
    if multi:
        bullets.append("The runs compared here start from the same sequencing data; differences "
                       "in the tables above follow from what was changed between them.")
    r = runs[-1]
    if r.total:
        bullets.append(f"{100*(r.host+r.euk)/r.total:.0f}% of analysed reads are host or other "
                       "non-bacterial material co-amplified by the primers.")
    if r.bacteria:
        bullets.append(f"About {100*r.kit/r.bacteria:.0f}% of the bacterial fraction matches the "
                       "reagent contaminant profile. Relative abundances are more honestly "
                       "computed over the sample-derived fraction than over the total.")
    bullets.append("A blank extraction control sequenced alongside the sample would turn the "
                   "reagent estimate into a measurement; a host-blocking primer would recover "
                   "most of the reads currently spent on host DNA.")
    out.append("<ul>" + "".join(f"<li>{b}</li>" for b in bullets) + "</ul>")

    css = _css(app, company, f"{sample_label} · {marker}")
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            + "".join(out) + "</body></html>")


# ══════════════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════════════
def _app_version(start: Path) -> str:
    """Stamp the report with the app version, looked up next to this script."""
    for d in [start] + list(start.parents)[:3]:
        p = d / "version.json"
        if p.exists():
            try:
                return "NextGen-Amplicon v" + json.loads(p.read_text())["version"]
            except Exception:
                pass
    return "NextGen-Amplicon"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Build a PDF report from NextGen-Amplicon result folders.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  report_builder.py results/RunA -o profile.pdf\n"
               "  report_builder.py results/RunA results/RunB -o comparison.pdf\n"
               "  report_builder.py results/* --label-from-folder -o all.pdf\n")
    ap.add_argument("folders", nargs="+", type=Path, help="one or more result folders")
    ap.add_argument("-o", "--out", type=Path, default=Path("report.pdf"), help="output PDF")
    ap.add_argument("--html", type=Path, help="also write the HTML here")
    ap.add_argument("--title", help="report title")
    ap.add_argument("--subtitle", help="one-line description under the title")
    ap.add_argument("--labels", help="comma-separated names for the runs, in order")
    ap.add_argument("--company", default="", help="company name for the page footer")
    ap.add_argument("--app", help="program name for the page footer")
    ap.add_argument("--top-asv", type=int, default=30, metavar="N",
                    help="how many of the most abundant ASVs to table in a single-run report "
                         "(default 30; 50 and 100 are the other usual choices)")
    ap.add_argument("--compare", action="store_true",
                    help="force the comparison layout even for a single folder")
    ap.add_argument("--kit-list", type=Path,
                    help="text file of reagent-contaminant genera, one per line, replacing the "
                         "built-in list")
    ap.add_argument("--version", action="version", version=f"report_builder {__version__}")
    a = ap.parse_args(argv)

    kitome = set(KITOME)
    if a.kit_list:
        kitome = {l.strip() for l in a.kit_list.read_text(encoding="utf-8").splitlines()
                  if l.strip() and not l.startswith("#")}

    labels = [s.strip() for s in a.labels.split(",")] if a.labels else []
    runs = []
    for i, f in enumerate(a.folders):
        try:
            runs.append(Run(f, kitome, labels[i] if i < len(labels) else None))
        except Exception as e:
            print(f"  [skip] {f}: {e}", file=sys.stderr)
    if not runs:
        print("No readable result folder given.", file=sys.stderr)
        return 1

    app = a.app or _app_version(Path(__file__).resolve().parent)
    # One folder is the ordinary case and gets the result report; the comparison layout
    # only makes sense with something to compare against.
    if len(runs) == 1 and not a.compare:
        doc = build_profile_html(runs[0], title=a.title, app=app, company=a.company,
                                 subtitle=a.subtitle, top_asv=a.top_asv)
    else:
        doc = build_html(runs, title=a.title, app=app, company=a.company, subtitle=a.subtitle)

    html_path = a.html or a.out.with_suffix(".html")
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(doc, encoding="utf-8")

    try:
        from weasyprint import HTML as _H
    except ImportError:
        print(f"weasyprint is not installed, so no PDF was produced.\n"
              f"  HTML written to {html_path} — open it in a browser and print to PDF,\n"
              f"  or install it with:  pip install weasyprint", file=sys.stderr)
        return 2
    a.out.parent.mkdir(parents=True, exist_ok=True)
    _H(string=doc, base_url=str(html_path.parent)).write_pdf(str(a.out))
    kind = "profile" if (len(runs) == 1 and not a.compare) else f"comparison of {len(runs)}"
    print(f"{a.out}  ({kind})")
    return 0




# ══════════════════════════════════════════════════════════════════════════════
#  Figures used by the single-run profile
# ══════════════════════════════════════════════════════════════════════════════
def ceiling_verdict(run):
    """Is truncLen cutting the amplicon short? Returns (clipped: bool, sentence).

    A fixed tolerance is the wrong test. A distribution that decays to a single read
    two bases below the ceiling is not clipped; one that is still carrying real
    abundance when it hits the ceiling is, even if its maximum is a base or two under.
    So look at the mass piled against the edge, not only at the maximum.
    """
    h = run.lengths()
    if not h or not run.ceiling:
        return False, ""
    mx, ceil = max(h), run.ceiling
    edge = sum(v for L, v in h.items() if L >= ceil - 2)
    share = 100 * edge / run.total if run.total else 0
    if mx >= ceil - 1 or share > 2:
        return True, (f"The longest ASV is {mx} bp against a {ceil} bp ceiling, with "
                      f"{fmt(edge)} reads ({share:.1f}%) stacked against that edge \u2014 "
                      "longer fragments could not be assembled and were discarded. Raising "
                      "truncLen, if the read lengths allow it, would recover them.")
    return False, (f"The longest ASV recovered is {mx} bp, {ceil-mx} bp clear of the "
                   f"{ceil} bp ceiling, so the truncation settings are not cutting the "
                   "amplicon short.")


def svg_readflow(run) -> str:
    """Where the reads went, stage by stage.

    A funnel rather than a table because the shape is the point: a step that loses far
    more than its neighbours is visible instantly and does not have to be worked out
    from four numbers.
    """
    stages = [("Input", "input"), ("Filtered", "filtered"),
              ("Merged", "merged"), ("Non-chimeric", "nonchim")]
    tot = Counter()
    for t in run.tracking:
        for _lbl, k in stages:
            try:
                tot[k] += float(t.get(k) or 0)
            except (TypeError, ValueError):
                pass
    base = tot["input"] or 1
    present = [(l, k) for l, k in stages if tot[k] > 0]
    if len(present) < 2:
        return ""
    W, x0, bw, rowh = 520, 92, 330, 30
    H = len(present) * rowh + 14
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for i, (lbl, k) in enumerate(present):
        y = 6 + i * rowh
        frac = tot[k] / base
        w = max(bw * frac, 1.5)
        p.append(f'<text x="{x0-8}" y="{y+15}" font-size="9" fill="{INK}" text-anchor="end">{lbl}</text>')
        p.append(f'<rect x="{x0}" y="{y}" width="{bw}" height="21" fill="{RULE}" opacity="0.5"/>')
        p.append(f'<rect x="{x0}" y="{y}" width="{w:.1f}" height="21" fill="{ACC}"/>')
        p.append(f'<text x="{x0+6}" y="{y+15}" font-size="9" fill="#fff" '
                 f'font-family="DejaVu Sans Mono">{fmt(tot[k])}</text>')
        p.append(f'<text x="{x0+bw+8}" y="{y+15}" font-size="8.5" fill="{MUT}" '
                 f'font-family="DejaVu Sans Mono">{100*frac:.1f}%</text>')
        if i:
            drop = tot[present[i-1][1]] - tot[k]
            if drop > 0:
                p.append(f'<text x="{W-4}" y="{y+15}" font-size="8" fill="{FAINT}" '
                         f'text-anchor="end" font-family="DejaVu Sans Mono">'
                         f'-{fmt(drop)}</text>')
    p.append("</svg>")
    return "".join(p)


def svg_bars(pairs, total, color=ACC, kit=None, unit="%") -> str:
    """Generic labelled horizontal bar list: [(name, value), ...]."""
    if not pairs:
        return ""
    kit = kit or set()
    W, colw, rowh = 520, 520, 24
    mx = max(v for _n, v in pairs) or 1
    H = len(pairs) * rowh + 8
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for i, (nm, v) in enumerate(pairs):
        y = 4 + i * rowh
        isk = nm in kit
        p.append(f'<text x="0" y="{y+7}" font-size="8.2" fill="{INK}">'
                 f'{esc(nm)}{" *" if isk else ""}</text>')
        p.append(f'<text x="{colw}" y="{y+7}" font-size="8" fill="{MUT}" text-anchor="end" '
                 f'font-family="DejaVu Sans Mono">{fmt(v)} &#183; {100*v/total:.1f}{unit}</text>')
        p.append(f'<rect x="0" y="{y+10}" width="{colw*v/mx:.1f}" height="5.5" '
                 f'fill="{KIT_C if isk else color}"/>')
    p.append("</svg>")
    return "".join(p)


def svg_rarefaction(run) -> str:
    """Observed richness against sequencing depth.

    A curve that has flattened means the depth was sufficient to see what is there; one
    still climbing at the right-hand edge means deeper sequencing would find more.
    """
    rows = run.rarefaction
    if not rows or len(rows) < 3:
        return ""
    cols = list(rows[0].keys())
    depth_col = cols[0]
    series = cols[1:]
    if not series:
        return ""
    W, H, x0, y0, pw, ph = 520, 170, 40, 12, 400, 118
    try:
        xs = [float(r[depth_col]) for r in rows]
    except (TypeError, ValueError):
        return ""
    mx_x = max(xs) or 1
    mx_y = 0.0
    for sname in series:
        for r in rows:
            try:
                mx_y = max(mx_y, float(r[sname] or 0))
            except (TypeError, ValueError):
                pass
    mx_y = mx_y or 1
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for gy in (0, 0.5, 1.0):
        y = y0 + ph - ph * gy
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}" stroke="{RULE}"/>')
        p.append(f'<text x="{x0-5}" y="{y+3:.1f}" font-size="7" fill="{FAINT}" text-anchor="end" '
                 f'font-family="DejaVu Sans Mono">{fmt(mx_y*gy)}</text>')
    for si, sname in enumerate(series):
        pts = []
        for r in rows:
            try:
                x = x0 + pw * float(r[depth_col]) / mx_x
                y = y0 + ph - ph * float(r[sname] or 0) / mx_y
                pts.append(f"{x:.1f},{y:.1f}")
            except (TypeError, ValueError):
                pass
        if pts:
            col = SERIES[si % len(SERIES)]
            p.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="1.6"/>')
            lx, ly = pts[-1].split(",")
            p.append(f'<text x="{float(lx)+4:.1f}" y="{float(ly)+3:.1f}" font-size="7.5" '
                     f'fill="{col}">{esc(sname)}</text>')
    for f in (0, 0.5, 1.0):
        x = x0 + pw * f
        p.append(f'<text x="{x:.1f}" y="{y0+ph+13}" font-size="7" fill="{FAINT}" '
                 f'text-anchor="middle" font-family="DejaVu Sans Mono">{fmt(mx_x*f)}</text>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-4}" font-size="8" fill="{MUT}" text-anchor="middle">'
             f'Sequencing depth (reads)</text>')
    p.append(f'<text x="8" y="{y0+8}" font-size="8" fill="{MUT}">ASVs observed</text>')
    p.append("</svg>")
    return "".join(p)


if __name__ == "__main__":
    sys.exit(main())
