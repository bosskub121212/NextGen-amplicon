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

  2. Haemotropic mycoplasmas arrive under two genus names, Eperythrozoon and
     Haemobartonella. Left alone, one group appears twice and neither row reaches
     its true abundance. Merged into a single line. (Those two names are the
     CURRENT ones: the transfer of these organisms into Mycoplasma proposed in
     2001 was never validly published, and Judicial Opinion 92 (2014) denied the
     Request to place the new combinations on the Approved Lists. Earlier versions
     of this file described them as superseded, which is the wrong way round.)

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

__version__ = "2.0.0"

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

# Two current genus names for one group of haemotropic mycoplasmas. The 2001
# transfer into Mycoplasma was never validly published and Opinion 92 (2014)
# denied it, so these are not superseded names — but they do split one biological
# group across two rows, which is why they are merged under a neutral label.
HAEMOPLASMA = {"Eperythrozoon", "Haemobartonella"}
HAEMOPLASMA_LABEL = "Haemoplasma group"

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
        # ONT 16S (Emu): one row per TAXON, not per sequence variant; abundances
        # are Emu's relative abundances; read stages are input/QC/primers/chimera.
        self.emu = (self.summary.get("pipeline") == "emu"
                    or (self.path / "emu_output").is_dir())

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
        if not tax_rows and self.emu:
            tax_rows = self._rows("taxonomy.csv")          # emu_pipeline's name for it
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
        if self.emu:
            seq_col = cols[0]                              # tax_id, first column
        else:
            seq_col = "sequence" if "sequence" in cols else cols[-1]
        self.samples = [c for c in cols if c != seq_col and c != ""]
        self.tracking = self._rows("read_tracking.csv")
        # Emu reports relative abundance (the pipeline stores it as parts per
        # million). Scaled to each sample's classified reads it becomes Emu's
        # estimated read count, so every "reads" figure below means reads.
        scale = {s: 1.0 for s in self.samples}
        # From v2.15 Emu writes real estimated read counts (--keep-counts) and
        # says so in summary.json; older runs hold relative abundance x 1e6 and
        # are rescaled to each sample's classified reads.
        self.emu_counts = str(self.summary.get("emu_counts") or ("ppm" if self.emu else ""))
        if self.emu and self.emu_counts != "estimated_reads":
            final = {}
            for t in self.tracking:
                try:
                    final[str(t.get("sample", ""))] = float(t.get("final") or 0)
                except (TypeError, ValueError):
                    pass
            tot = {s: 0.0 for s in self.samples}
            for r in asv_rows:
                for s_ in self.samples:
                    try:
                        tot[s_] += float(r[s_] or 0)
                    except (TypeError, ValueError):
                        pass
            for s_ in self.samples:
                if final.get(s_) and tot[s_]:
                    scale[s_] = final[s_] / tot[s_]
        # Per-sample counts, and the total each ASV carries across all samples.
        self.per_sample = {s: 0.0 for s in self.samples}
        self.abund: dict[str, float] = {}
        self.mat: dict[str, dict[str, float]] = {}      # feature -> sample -> reads
        for r in asv_rows:
            seq = r[seq_col]
            tot = 0.0
            row = self.mat.setdefault(seq, {s: 0.0 for s in self.samples})
            for s in self.samples:
                try:
                    v = float(r[s] or 0) * scale[s]
                except (TypeError, ValueError):
                    v = 0.0
                self.per_sample[s] += v
                row[s] += v
                tot += v
            self.abund[seq] = self.abund.get(seq, 0.0) + tot
        self.total = sum(self.abund.values())

        self.alpha = self._rows("alpha_diversity.csv")
        self.rarefaction = self._rows("rarefaction.csv")
        self.orientation = self._json("orientation_report.json")

    # ── derived views ─────────────────────────────────────────────────────────
    def _clean_genus(self, g):
        """One label per genus: no blanks, no '--other', haemoplasma names merged."""
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
        # Eukaryotic reads are off-target only for a bacterial marker; for 18S,
        # 12S, ITS or COI they ARE the community and must not be set aside.
        marker = str(self.summary.get("marker") or self.params.get("marker") or "16S").upper()
        self.euk_target = any(k in marker for k in ("18S", "12S", "ITS", "COX", "COI", "NEMA", "EUK"))
        self.euk_seqs = [] if self.euk_target else [
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
                v = str(v).strip()
                if self.emu:
                    # Emu writes the binomial ("Brevundimonas aurantiaca"), and
                    # "<Genus> sp." for genus-only taxa — not a species call
                    if v.endswith(" sp.") or v == "sp.":
                        continue
                    name = v if v.split(" ")[0] == g or " " in v else f"{g} {v}"
                    self.species[name] = (self.species.get(name, (0, None))[0]
                                          + self.abund.get(s, 0), None)
                    continue
                self.species[f"{g} {v}"] = (self.abund.get(s, 0), len(s))
        self.species = OrderedDict(
            sorted(self.species.items(), key=lambda kv: -kv[1][0])
        )

    # ── reporting helpers ─────────────────────────────────────────────────────
    def rank_pct(self, rank) -> float | None:
        if rank not in self.rank_cols or not self.total:
            return None
        hit = [s for s, r in self.tax.items()
               if r.get(rank) and not BLANK_RE.match(str(r[rank]).strip())
               # Emu writes "<Genus> sp." where it has no species
               and not (self.emu and rank == "Species" and str(r[rank]).strip().endswith("sp."))]
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
                "id": f"T{i}" if self.emu else f"ASV{i}",
                "len": None if self.emu else len(seq), "reads": reads,
                "pct": 100 * reads / self.total if self.total else 0,
                "phylum": (r.get("Phylum") or "\u2014"),
                "genus": (self._clean_genus(r.get("Genus")) or "\u2014"),
                "deepest": deepest, "lineage": " > ".join(lineage) or "unassigned",
                "host": host, "seq": seq,
            })
        return out

    def lengths(self) -> Counter:
        c = Counter()
        if self.emu:            # rows are taxa: there is no sequence length
            return c
        for s, v in self.abund.items():
            c[len(s)] += v
        return c

    @property
    def db_name(self) -> str:
        p = self.params.get("dbPath") or ""
        if p:
            pp = Path(p)
            # an Emu database is a folder; its file is always species_taxid.fasta
            return pp.parent.name if self.emu and pp.suffix else pp.name
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
        return self.summary.get("n_asvs") or self.summary.get("n_taxa") or len(self.tax)

    # wording and read stages that differ between DADA2 and Emu runs
    @property
    def unit(self) -> str:
        return "taxa" if self.emu else "sequence variants"

    @property
    def stages(self):
        if self.emu:
            last = "Classified" if _skipped(self).get("chimera") else "Non-chimeric"
            return [("Input", "input"), ("QC (Q / length)", "qc_filtered"),
                    ("Primers trimmed", "primer_trimmed"), (last, "final")]
        return [("Input", "input"), ("Filtered", "filtered"),
                ("Merged", "merged"), ("Non-chimeric", "nonchim")]

    @property
    def retained_key(self) -> str:
        return "final" if self.emu else "nonchim"

    def emu_config(self):
        p, sm = self.params, self.summary
        return [("Classifier", "Emu — minimap2 alignment + expectation–maximisation "
                 "over full-length reads"),
                ("Reference database", self.db_name),
                ("Region", sm.get("region") or p.get("ont_region") or "—"),
                ("Forward primer", p.get("primer_f") or "—"),
                ("Reverse primer", p.get("primer_r") or "—"),
                ("Minimum abundance", sm.get("emu_min_abundance") or p.get("ont_min_abundance")
                 or "—")]


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
    # Reports go to customers as they are: no program or company name in the
    # footer — only the job label (left) and the page number (right).
    foot = ""
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
    if ref.emu:
        kpis = [(fmt(inp) if inp else "—", "reads sequenced"),
                (fmt(ref.total), "reads classified"),
                (" / ".join(str(r.n_asvs) for r in runs), "taxa detected"),
                (" / ".join(f"{r.rank_pct('Genus') or 0}%" for r in runs), "reads with a genus")]
    else:
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
    if ref.emu:                      # Emu runs: no truncLen / maxEE / bootstrap
        keys = [("Reference database", "db_name")]
        prm = [(k, (lambda r, k=k: dict(r.emu_config()).get(k)))
               for k, _v in ref.emu_config() if k != "Reference database"]
    rows = [[k] + [esc(getattr(r, a)) for r in runs] for k, a in keys]
    for k, fn in prm:
        vals = [fn(r) for r in runs]
        if any(v not in (None, "", "None") for v in vals):
            rows.append([k] + [esc(v if v is not None else "—") for v in vals])
    out.append(_tbl(["Parameter"] + [esc(n) for n in names], rows,
                    "Rows that differ between columns are the variables under comparison."
                    if multi else None))

    if ref.tracking:
        cols = [k for _l, k in ref.stages if k in ref.tracking[0]]
        heads = [l for l, k in ref.stages if k in ref.tracking[0]]
        trows = []
        for t in ref.tracking:
            name = next(iter(t.values()))
            base = float(t.get("input") or 0) or 1
            trows.append([esc(name)] + [fmt(t.get(c) or 0) for c in cols]
                         + [f'{100*float(t.get(ref.retained_key) or 0)/base:.1f}%'])
        out.append("<h3>Read tracking</h3>")
        out.append(_tbl(["Sample"] + heads + ["Retained"], trows,
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
        out.append(f'<div class="note acc"><p style="margin:0"><b>Haemotropic mycoplasmas '
                   f"&mdash; {fmt(haemo[-1])} reads, "
                   f"{100*haemo[-1]/r.bacteria:.1f}% of bacterial reads.</b> Reference databases "
                   "label these <i>Eperythrozoon</i> and <i>Haemobartonella</i>. Those are the "
                   "names LPSN currently treats as correct: the transfer of these organisms into "
                   "<i>Mycoplasma</i> proposed in 2001 was never validly published, and Judicial "
                   "Opinion 92 (2014) denied the Request to place the new combinations on the "
                   "Approved Lists. They are one group under two current names, and are merged "
                   "here for that reason. Haemoplasmas are blood-borne pathogens of dogs; their "
                   "route of transmission is unresolved, so recovery from an ectoparasite sample "
                   "is not by itself evidence of arthropod vectoring."
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
    if ref.emu:
        out.append("<p>Species from Emu over full-length reads. A database built to genus "
                   "level (SILVA) has no species to give — an empty column there is the "
                   "reference, not the data.</p>")
    else:
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
            if not ref.emu:
                row.append(f"{first[1]} bp")
            rows.append(row)
        out.append(_tbl(["Species"] + [esc(n) for n in names] + ([] if ref.emu else ["ASV length"]),
                        rows,
                        f"{fmt(sum(v[0] for v in runs[-1].species.values()))} of "
                        f"{fmt(runs[-1].total)} reads — "
                        f"{runs[-1].rank_pct('Species')}% of the analysed data."))
    else:
        out.append("<p>No read reached a species-level assignment in any run.</p>")
    if not ref.emu:
      out.append("<p>A low species rate is the expected outcome for a short amplicon. The "
               "V3&ndash;V4 region covers roughly 440 bp of a 1,540 bp gene, and uncultured "
               "environmental bacteria frequently have no species-annotated entry in any 16S "
               "reference to match against. The figure reflects the marker and the sample rather "
               "than the analysis; a different database or different settings do not change it. "
               "Full-length 16S (27F/1492R on PacBio or ONT) is the appropriate marker where "
               "species-level identification is required.</p>")

    # ── 7. summary ────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">7. Summary</h2>')
    srows = [["Reads classified" if ref.emu else "Reads analysed"] + [fmt(r.total) for r in runs],
             ["Taxa" if ref.emu else "Sequence variants"] + [str(r.n_asvs) for r in runs]]
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
#  Full single-run report (2.13): per-sample analysis, computed in pure Python
#  from asv_table + taxonomy so that DADA2 and Emu runs get the same figures
#  whether or not viz_pipeline.R left r_tables/ behind.
# ══════════════════════════════════════════════════════════════════════════════
import math

PAL30 = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
         "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
         "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#393B79", "#637939", "#8C6D31",
         "#843C39", "#7B4173", "#5254A3", "#8CA252", "#BD9E39", "#AD494A", "#A55194",
         "#6B6ECF", "#E7BA52"]
OTHER_C, UNCL_C = "#98A4AA", "#CBD3D7"
MONO = 'font-family="DejaVu Sans Mono"'
PLURAL = {"Phylum": "phyla", "Class": "classes", "Order": "orders", "Family": "families",
          "Genus": "genera", "Species": "species"}


def _rc(seq: str) -> str:
    comp = str.maketrans("ACGTRYKMBVDHNSWacgtrykmbvdhnsw", "TGCAYRMKVBHDNSWtgcayrmkvbhdnsw")
    return (seq or "").translate(comp)[::-1]


def _pal(i: int) -> str:
    return PAL30[i % len(PAL30)]


def _short(s, n):
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def _col_heads(S, n=12):
    """Sample names short enough for a column header and still distinct: a prefix every
    name shares (up to a separator, e.g. "Bacteria_") is dropped, and if cutting to n
    characters still makes two headers equal, the end of the name is kept instead."""
    names = [str(s) for s in S]
    if len(names) > 1:
        pre = names[0]
        for x in names[1:]:
            while not x.startswith(pre):
                pre = pre[:-1]
        cut = max(pre.rfind(c) for c in "_-. ") + 1
        if cut > 0 and all(len(x) > cut for x in names):
            names = [x[cut:] for x in names]
    out = [_short(x, n) for x in names]
    if len(set(out)) < len(out):
        out = [x if len(x) <= n else "…" + x[-(n - 1):] for x in names]
    return out


# ── numerical helpers ─────────────────────────────────────────────────────────
def _alpha_of(counts) -> dict:
    """Observed, Chao1 (bias-corrected), ACE, Shannon (ln) and Simpson (1 - D),
    the same estimators phyloseq::estimate_richness / vegan::estimateR use."""
    c = [int(round(x)) for x in counts if round(x) > 0]
    S, N = len(c), sum(c)
    if not N:
        return {"Observed": 0, "Chao1": 0, "ACE": 0, "Shannon": 0, "Simpson": 0}
    f1, f2 = c.count(1), c.count(2)
    chao1 = S + f1 * (f1 - 1) / (2 * (f2 + 1))
    rare = [x for x in c if x <= 10]
    s_rare, n_rare = len(rare), sum(rare)
    ace = float(S)
    if n_rare and n_rare != f1:
        cace = 1 - f1 / n_rare
        g2 = max(s_rare / cace * sum(i * (i - 1) * rare.count(i) for i in range(1, 11))
                 / (n_rare * (n_rare - 1)) - 1, 0) if n_rare > 1 else 0
        ace = (S - s_rare) + s_rare / cace + f1 / cace * g2
    ps = [x / N for x in c]
    return {"Observed": S, "Chao1": chao1, "ACE": ace,
            "Shannon": -sum(p * math.log(p) for p in ps),
            "Simpson": 1 - sum(p * p for p in ps)}


def _hurlbert(counts, n) -> float:
    """Expected number of features in a random subsample of n reads (vegan::rarefy)."""
    N = sum(counts)
    if n >= N:
        return float(sum(1 for x in counts if x > 0))
    lg = math.lgamma
    lc = lg(N + 1) - lg(n + 1) - lg(N - n + 1)
    s = 0.0
    for x in counts:
        if x <= 0:
            continue
        if N - x < n:
            s += 1.0
            continue
        s += 1.0 - math.exp(lg(N - x + 1) - lg(n + 1) - lg(N - x - n + 1) - lc)
    return s


def _bray(a, b) -> float:
    """Bray-Curtis between two composition vectors (proportions)."""
    num = sum(abs(x - y) for x, y in zip(a, b))
    den = sum(a) + sum(b)
    return num / den if den else 0.0


def _jacobi(A, sweeps=80):
    """Eigen-decomposition of a small symmetric matrix. Returns (values, vectors as columns)."""
    n = len(A)
    A = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(A[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        if off < 1e-20:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(A[p][q]) < 1e-14:
                    continue
                th = (A[q][q] - A[p][p]) / (2 * A[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p], A[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k], A[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p], V[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [A[i][i] for i in range(n)], V


def _pcoa(D):
    """Classical multidimensional scaling. Returns (coords [[ax1, ax2], ...], [pct1, pct2])."""
    n = len(D)
    D2 = [[D[i][j] ** 2 for j in range(n)] for i in range(n)]
    rm = [sum(r) / n for r in D2]
    gm = sum(rm) / n
    B = [[-0.5 * (D2[i][j] - rm[i] - rm[j] + gm) for j in range(n)] for i in range(n)]
    vals, vecs = _jacobi(B)
    order = sorted(range(n), key=lambda k: -vals[k])
    pos = sum(v for v in vals if v > 1e-12) or 1.0
    coords = [[0.0, 0.0] for _ in range(n)]
    pct = [0.0, 0.0]
    for ax in range(2):
        if ax >= len(order):
            break
        k = order[ax]
        lam = vals[k]
        if lam <= 1e-12:
            continue
        pct[ax] = 100 * lam / pos
        # sign convention: largest absolute loading positive, so reruns look the same
        col = [vecs[i][k] for i in range(n)]
        if col[max(range(n), key=lambda i: abs(col[i]))] < 0:
            col = [-x for x in col]
        for i in range(n):
            coords[i][ax] = col[i] * math.sqrt(lam)
    return coords, pct


def _upgma(D):
    """Average-linkage clustering. Nodes: ("L", i) or ("N", left, right, height)."""
    n = len(D)
    if n == 0:
        return None
    nodes = {i: ("L", i) for i in range(n)}
    size = {i: 1 for i in range(n)}
    dist = {(i, j): D[i][j] for i in range(n) for j in range(i + 1, n)}
    nxt = n
    while len(nodes) > 1:
        (a, b), d = min(dist.items(), key=lambda kv: (kv[1], kv[0]))
        nodes[nxt] = ("N", nodes.pop(a), nodes.pop(b), d)
        sa, sb = size.pop(a), size.pop(b)
        size[nxt] = sa + sb
        for k in list(nodes):
            if k == nxt:
                continue
            dak = dist.pop((min(a, k), max(a, k)))
            dbk = dist.pop((min(b, k), max(b, k)))
            dist[(k, nxt)] = (sa * dak + sb * dbk) / (sa + sb)
        dist.pop((min(a, b), max(a, b)), None)
        nxt += 1
    return next(iter(nodes.values()))


def _leaves(node):
    if node is None:
        return []
    if node[0] == "L":
        return [node[1]]
    return _leaves(node[1]) + _leaves(node[2])


def _height(node):
    return 0.0 if node is None or node[0] == "L" else node[3]


def _dendro(node, leafpos):
    """Segments of a dendrogram in (position, height) space."""
    if node[0] == "L":
        return leafpos[node[1]], 0.0, []
    pl, hl, sl = _dendro(node[1], leafpos)
    pr, hr, sr = _dendro(node[2], leafpos)
    h = node[3]
    return (pl + pr) / 2, h, sl + sr + [((pl, hl), (pl, h)), ((pr, hr), (pr, h)),
                                        ((pl, h), (pr, h))]


def _euclid(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


# ── per-sample analysis of one run ────────────────────────────────────────────
class Analysis:
    """Everything the full report draws, computed once per run."""

    def __init__(self, run):
        self.run = run
        self.samples = list(run.samples)
        self.comm = [q for q in run.bact_seqs if q in run.mat]   # host organelle / eukaryote reads set aside
        self.mat = run.mat
        self.comm_tot = {s: sum(self.mat[q][s] for q in self.comm) for s in self.samples}
        self._rank = {}
        self.groups = self._groups()

    # sample → group, from the run's metadata
    def _groups(self):
        g, desc = {}, {}
        for m in (self.run.params.get("metadata") or []):
            sid = str(m.get("sampleId") or m.get("sample") or "").strip()
            if sid:
                if m.get("group"):
                    g[sid] = str(m["group"]).strip()
                if m.get("description"):
                    desc[sid] = str(m["description"]).strip()
        if not g:
            for r in self.run._rows("metadata.csv"):
                vals = list(r.values())
                sid = str(r.get("sampleId") or r.get("sample") or (vals[0] if vals else "")).strip()
                grp = r.get("group") or r.get("Group")
                if sid and grp:
                    g[sid] = str(grp).strip()
                if sid and r.get("description"):
                    desc[sid] = str(r["description"]).strip()
        self.desc = desc
        return {s: g[s] for s in self.samples if s in g}

    # ── taxonomy ──────────────────────────────────────────────────────────────
    def label(self, seq, rank) -> str:
        run, r = self.run, self.run.tax.get(seq, {})
        if rank == "Species":
            sp = str(r.get("Species") or "").strip()
            if not sp or BLANK_RE.match(sp) or sp.endswith(" sp.") or sp == "sp.":
                return "Unclassified"
            g = run._clean_genus(r.get("Genus"))
            if run.emu or " " in sp or not g:
                return sp
            return f"{g} {sp}"
        v = str(r.get(rank) or "").strip()
        if not v or BLANK_RE.match(v):
            return "Unclassified"
        if rank == "Genus":
            return run._clean_genus(v) or "Unclassified"
        return UNRESOLVED_RE.sub("", v) or "Unclassified"

    def rank_counts(self, rank):
        """{taxon: {sample: reads}} over the community reads."""
        if rank not in self._rank:
            m = {}
            for q in self.comm:
                t = self.label(q, rank)
                row = m.setdefault(t, {s: 0.0 for s in self.samples})
                for s in self.samples:
                    row[s] += self.mat[q][s]
            self._rank[rank] = m
        return self._rank[rank]

    def rank_pct(self, rank):
        m = self.rank_counts(rank)
        return {t: {s: (100 * v[s] / self.comm_tot[s] if self.comm_tot[s] else 0.0)
                    for s in self.samples} for t, v in m.items()}

    def top(self, rank, n=30):
        """(top taxa by mean relative abundance, pct table incl. Other / Unclassified)."""
        pct = self.rank_pct(rank)
        ns = len(self.samples) or 1
        named = sorted((t for t in pct if t != "Unclassified"),
                       key=lambda t: -sum(pct[t].values()) / ns)
        keep = named[:n]
        out = {t: pct[t] for t in keep}
        rest = named[n:]
        if rest:
            out["Other"] = {s: sum(pct[t][s] for t in rest) for s in self.samples}
        if "Unclassified" in pct:
            out["Unclassified"] = pct["Unclassified"]
        return keep, out

    def dominant(self, rank, s, k=3):
        pct = self.rank_pct(rank)
        items = sorted(((t, v[s]) for t, v in pct.items() if t != "Unclassified"),
                       key=lambda kv: -kv[1])
        return [(t, v) for t, v in items[:k] if v > 0]

    # ── features per sample ───────────────────────────────────────────────────
    def n_features(self, s):
        return sum(1 for q in self.run.abund if self.mat[q][s] >= 0.5)

    def alpha(self):
        if not hasattr(self, "_alpha"):
            self._alpha = {s: _alpha_of([self.mat[q][s] for q in self.run.abund])
                           for s in self.samples}
        return self._alpha

    def rarefaction(self, points=30):
        """{sample: [(depth, expected features), ...]} — Hurlbert expectation."""
        if hasattr(self, "_rare"):
            return self._rare
        cnt = {s: [int(round(self.mat[q][s])) for q in self.run.abund] for s in self.samples}
        mx = max((sum(c) for c in cnt.values()), default=0)
        out = {}
        for s, c in cnt.items():
            N = sum(c)
            if N < 2:
                continue
            depths = sorted({max(1, int(N * (i / points) ** 1.6)) for i in range(1, points + 1)} | {N})
            out[s] = [(d, _hurlbert(c, d)) for d in depths]
        self._rare, self.rare_max = out, mx
        return out

    def plateau(self, s):
        """Share of richness gained over the last 10% of the sample's depth."""
        c = [int(round(self.mat[q][s])) for q in self.run.abund]
        N = sum(c)
        if N < 20:
            return None
        full = _hurlbert(c, N)
        at90 = _hurlbert(c, int(N * 0.9))
        return 100 * (full - at90) / full if full else None

    # ── between samples ───────────────────────────────────────────────────────
    def bray(self):
        if hasattr(self, "_bray"):
            return self._bray
        feats = list(self.run.abund)
        tot = {s: sum(self.mat[q][s] for q in feats) for s in self.samples}
        vec = {s: [self.mat[q][s] / tot[s] if tot[s] else 0.0 for q in feats] for s in self.samples}
        n = len(self.samples)
        D = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                D[i][j] = D[j][i] = _bray(vec[self.samples[i]], vec[self.samples[j]])
        self._bray = D
        return D

    def archaea(self):
        """[(sample, archaeal reads, % of community, dominant archaeal taxon)]"""
        run = self.run
        arch = [q for q in self.comm
                if any(str(run.tax[q].get(k) or "").strip() == "Archaea"
                       for k in ("Domain", "Kingdom", "Superkingdom", "superkingdom"))]
        rows = []
        for s in self.samples:
            reads = sum(self.mat[q][s] for q in arch)
            best = Counter()
            for q in arch:
                if self.mat[q][s] > 0:
                    g = self.label(q, "Genus")
                    best[g if g != "Unclassified" else self.label(q, "Phylum")] += self.mat[q][s]
            top = best.most_common(1)
            rows.append((s, reads, 100 * reads / self.comm_tot[s] if self.comm_tot[s] else 0.0,
                         top[0][0] if top else "—"))
        return rows

    def permanova(self):
        """(R2, p) from beta_stats.txt when the test had residual degrees of freedom."""
        p = self.run.path / "beta_stats.txt"
        if not p.exists():
            return None
        txt = p.read_text(encoding="utf-8", errors="replace")
        if "No residual component" in txt:
            return None
        m = re.search(r"\n\s*(?:Model|grp_factor|group\w*)\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
                      r"\s+([\d.eE+-]+)\s+([\d.eE+-]+)", txt)
        if not m:
            return None
        try:
            return float(m.group(3)), float(m.group(5)), float(m.group(4))
        except ValueError:
            return None

# ── figures for the full report ───────────────────────────────────────────────
def _bcf(v):
    """Two decimals, three where rounding would print a misleading 1.00 or 0.00."""
    return f"{v:.3f}" if (0.99 <= v < 1 or 0 < v <= 0.01) else f"{v:.2f}"


def _svg_open(W, H):
    return [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px" '
            f'xmlns="http://www.w3.org/2000/svg">']


def _xlabels(p, labels, xs, y, maxlen=22, size=8.4):
    rot = len(labels) > 3 or max((len(str(l)) for l in labels), default=0) > 10
    for l, x in zip(labels, xs):
        t = esc(_short(l, maxlen))
        if rot:
            p.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{INK}" '
                     f'text-anchor="end" transform="rotate(-40 {x:.1f} {y:.1f})">{t}</text>')
        else:
            p.append(f'<text x="{x:.1f}" y="{y+4:.1f}" font-size="{size}" fill="{INK}" '
                     f'text-anchor="middle">{t}</text>')


def svg_stacked(an, rank, n=30, samples=None, table=None, italic=None) -> str:
    """Top-n taxa per sample as 100% stacked bars, with a legend down the right."""
    samples = samples or an.samples
    if table is None:
        keep, table = an.top(rank, n)
    else:
        keep = [t for t in table if t not in ("Other", "Unclassified")]
    order = keep + [t for t in ("Other", "Unclassified") if t in table]
    col = {t: _pal(i) for i, t in enumerate(keep)}
    col["Other"], col["Unclassified"] = OTHER_C, UNCL_C
    italic = rank in ("Genus", "Species") if italic is None else italic
    longest = max((len(str(s)) for s in samples), default=4)
    W, x0, y0, ph = 520, 40, 22, 430
    leg_x = 322
    pw = leg_x - x0 - 16
    bot = 20 + min(longest, 22) * 4.0
    H = int(y0 + ph + bot + 10)
    p = _svg_open(W, H)
    p.append(f'<text x="{x0}" y="11" font-size="9.4" fill="{MUT}">Relative abundance — '
             f'{esc(rank)} level</text>')
    for gy in range(0, 101, 25):
        y = y0 + ph - ph * gy / 100
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
        p.append(f'<text x="{x0-4}" y="{y+2.5:.1f}" font-size="8" fill="{FAINT}" text-anchor="end" {MONO}>{gy}</text>')
    p.append(f'<text x="9" y="{y0+ph/2:.1f}" font-size="8.4" fill="{MUT}" text-anchor="middle" '
             f'transform="rotate(-90 9 {y0+ph/2:.1f})">Relative abundance (%)</text>')
    ns = max(len(samples), 1)
    slot = pw / ns
    bw = min(slot * 0.66, 84)
    xs = []
    for i, s in enumerate(samples):
        xc = x0 + slot * (i + 0.5)
        xs.append(xc)
        tot = sum(table[t].get(s, 0) for t in order) or 100
        acc = 0.0
        for t in order:
            v = table[t].get(s, 0) * 100 / tot
            if v <= 0:
                continue
            h = ph * v / 100
            y = y0 + ph - ph * (acc + v) / 100
            p.append(f'<rect x="{xc-bw/2:.1f}" y="{y:.2f}" width="{bw:.1f}" height="{h:.2f}" '
                     f'fill="{col[t]}" stroke="#fff" stroke-width=".25"><title>{esc(t)} {v:.2f}%</title></rect>')
            acc += v
    _xlabels(p, samples, xs, y0 + ph + 10)
    # legend
    rows = len(order)
    lh = min(14.0, (H - 30) / max(rows, 1))
    fs = min(8.6, lh * 0.70)
    p.append(f'<text x="{leg_x}" y="11" font-size="8.8" fill="{INK}" font-weight="700">{esc(rank)}</text>')
    for i, t in enumerate(order):
        y = 18 + i * lh
        p.append(f'<rect x="{leg_x}" y="{y:.1f}" width="{lh*0.75:.1f}" height="{lh*0.75:.1f}" fill="{col[t]}"/>')
        sty = ' font-style="italic"' if italic and t not in ("Other", "Unclassified") else ""
        p.append(f'<text x="{leg_x+lh:.1f}" y="{y+lh*0.68:.1f}" font-size="{fs:.1f}" fill="{INK}"{sty}>'
                 f'{esc(_short(t, 34))}</text>')
    p.append("</svg>")
    return "".join(p)


def _heat_col(v, vmax):
    """White → deep blue on a log scale of percent."""
    if v <= 0:
        return "#F4F7F8"
    f = math.log10(v + 0.01) - math.log10(0.01)
    f = max(0.0, min(1.0, f / (math.log10(vmax + 0.01) - math.log10(0.01))))
    a = (0xF4, 0xF7, 0xF8)
    b = (0x0B, 0x2E, 0x59)
    mid = (0x5B, 0x9B, 0xC8)
    if f < 0.5:
        c = [a[i] + (mid[i] - a[i]) * f * 2 for i in range(3)]
    else:
        c = [mid[i] + (b[i] - mid[i]) * (f - 0.5) * 2 for i in range(3)]
    return "#%02X%02X%02X" % tuple(int(x) for x in c)


def svg_heatmap(an, rank="Genus", n=30) -> str:
    """Top-n taxa × samples, clustered: rows by profile (UPGMA, Euclidean on log %),
    columns by Bray-Curtis when there are three or more samples."""
    keep, table = an.top(rank, n)
    if not keep:
        return ""
    samples = list(an.samples)
    # columns
    col_tree = None
    if len(samples) >= 3:
        vec = [[table[t][s] for t in keep] for s in samples]
        D = [[_bray(a, b) for b in vec] for a in vec]
        col_tree = _upgma(D)
        samples = [samples[i] for i in _leaves(col_tree)]
    # rows
    lv = [[math.log10(table[t][s] + 0.01) for s in an.samples] for t in keep]
    row_tree = _upgma([[_euclid(a, b) for b in lv] for a in lv]) if len(keep) > 1 else None
    rows = [keep[i] for i in _leaves(row_tree)] if row_tree else keep
    vmax = max((table[t][s] for t in keep for s in samples), default=1) or 1
    W = 520
    dend_w = 70
    lab_w = 160
    top_h = 52 if col_tree else 8
    cw = min(92, (W - dend_w - lab_w - 10) / max(len(samples), 1))
    ch = 15.0
    x0 = dend_w + 4
    y0 = 20 + top_h
    longest = max((len(str(s)) for s in samples), default=4)
    H = int(y0 + ch * len(rows) + 22 + min(longest, 22) * 4.0 + 38)
    p = _svg_open(W, H)
    p.append(f'<text x="0" y="11" font-size="9.4" fill="{MUT}">Top {len(keep)} {esc(rank.lower())} '
             f'— relative abundance (%)</text>')
    for j, s in enumerate(samples):
        for i, t in enumerate(rows):
            v = table[t][s]
            x, y = x0 + j * cw, y0 + i * ch
            c = _heat_col(v, vmax)
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="{c}" '
                     f'stroke="#fff" stroke-width=".5"/>')
            if cw >= 24:
                dark = c != "#F4F7F8" and int(c[1:3], 16) < 0x80
                txt = f"{v:.1f}" if v >= 0.05 else ("0" if v == 0 else "&lt;.1")
                p.append(f'<text x="{x+cw/2:.1f}" y="{y+ch*0.72:.1f}" font-size="7.2" '
                         f'fill="{"#fff" if dark else MUT}" text-anchor="middle" {MONO}>{txt}</text>')
    # row labels
    lx = x0 + cw * len(samples) + 4
    for i, t in enumerate(rows):
        sty = ' font-style="italic"' if rank in ("Genus", "Species") else ""
        p.append(f'<text x="{lx:.1f}" y="{y0+i*ch+ch*0.75:.1f}" font-size="8" fill="{INK}"{sty}>'
                 f'{esc(_short(t, 34))}</text>')
    # row dendrogram (root at the left)
    if row_tree:
        hmax = _height(row_tree) or 1
        pos = {k: y0 + (rows.index(keep[k]) + 0.5) * ch for k in range(len(keep))}
        _pp, _hh, segs = _dendro(row_tree, pos)
        for (p1, h1), (p2, h2) in segs:
            xa = x0 - 3 - (dend_w - 6) * h1 / hmax
            xb = x0 - 3 - (dend_w - 6) * h2 / hmax
            p.append(f'<line x1="{xa:.1f}" y1="{p1:.1f}" x2="{xb:.1f}" y2="{p2:.1f}" stroke="{MUT}" stroke-width=".7"/>')
    # column dendrogram (root at the top)
    if col_tree:
        hmax = _height(col_tree) or 1
        pos = {k: x0 + (samples.index(an.samples[k]) + 0.5) * cw for k in range(len(an.samples))}
        _pp, _hh, segs = _dendro(col_tree, pos)
        for (p1, h1), (p2, h2) in segs:
            ya = y0 - 3 - (top_h - 6) * h1 / hmax
            yb = y0 - 3 - (top_h - 6) * h2 / hmax
            p.append(f'<line x1="{p1:.1f}" y1="{ya:.1f}" x2="{p2:.1f}" y2="{yb:.1f}" stroke="{MUT}" stroke-width=".7"/>')
    yb = y0 + ch * len(rows)
    _xlabels(p, samples, [x0 + (j + 0.5) * cw for j in range(len(samples))], yb + 9)
    # colour key
    ky = H - 14
    kx = 0
    p.append(f'<text x="{kx}" y="{ky-4}" font-size="7.6" fill="{MUT}">% (log scale)</text>')
    stops = [0, 0.01, 0.1, 1, 10, vmax]
    stops = sorted({round(s, 2) for s in stops if s <= vmax})
    for k, sv in enumerate(stops):
        p.append(f'<rect x="{kx + k*34}" y="{ky}" width="34" height="7" fill="{_heat_col(sv, vmax)}" stroke="{RULE}" stroke-width=".4"/>')
        p.append(f'<text x="{kx + k*34 + 17}" y="{ky+14}" font-size="7.2" fill="{FAINT}" text-anchor="middle" {MONO}>{sv:g}</text>')
    p.append("</svg>")
    return "".join(p)


def svg_steps(run, an) -> str:
    """Reads per pipeline step, one line per sample (box per step from 4 samples up)."""
    if not run.tracking:
        return ""
    if run.emu:
        steps = [(l, k) for l, k in run.stages if k in run.tracking[0]]
    else:
        steps = [(c.capitalize() if c != "nonchim" else "Non-chimeric", c)
                 for c in ("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
                 if c in run.tracking[0]]
    if len(steps) < 2:
        return ""
    data = []
    for t in run.tracking:
        name = str(next(iter(t.values())))
        vals = []
        for _l, k in steps:
            try:
                vals.append(float(t.get(k) or 0))
            except (TypeError, ValueError):
                vals.append(0.0)
        data.append((name, vals))
    mx = max((v for _n, vs in data for v in vs), default=1) or 1
    nice = 10 ** int(math.log10(mx))
    mx = math.ceil(mx / nice) * nice
    W, H, x0, y0, pw, ph = 520, 290, 52, 14, 330, 224
    p = _svg_open(W, H)
    for k in range(5):
        v = mx * k / 4
        y = y0 + ph - ph * k / 4
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
        p.append(f'<text x="{x0-4}" y="{y+2.5:.1f}" font-size="8" fill="{FAINT}" text-anchor="end" {MONO}>{fmt(v)}</text>')
    xs = [x0 + pw * (i + 0.5) / len(steps) for i in range(len(steps))]
    if len(data) >= 4:
        for i, x in enumerate(xs):
            col = sorted(vs[i] for _n, vs in data)
            q = lambda f: col[min(len(col) - 1, int(round(f * (len(col) - 1))))]
            yq1, yq3 = y0 + ph - ph * q(.25) / mx, y0 + ph - ph * q(.75) / mx
            ymed = y0 + ph - ph * q(.5) / mx
            ylo, yhi = y0 + ph - ph * col[0] / mx, y0 + ph - ph * col[-1] / mx
            p.append(f'<line x1="{x:.1f}" y1="{ylo:.1f}" x2="{x:.1f}" y2="{yhi:.1f}" stroke="{FAINT}" stroke-width=".8"/>')
            p.append(f'<rect x="{x-11:.1f}" y="{yq3:.1f}" width="22" height="{max(yq1-yq3,0.8):.1f}" fill="#EAF1F4" stroke="{MUT}" stroke-width=".7"/>')
            p.append(f'<line x1="{x-11:.1f}" y1="{ymed:.1f}" x2="{x+11:.1f}" y2="{ymed:.1f}" stroke="{INK}" stroke-width="1.1"/>')
    for si, (name, vs) in enumerate(data):
        c = SERIES[si] if len(data) <= len(SERIES) else _pal(si)
        pts = " ".join(f"{x:.1f},{y0+ph-ph*v/mx:.1f}" for x, v in zip(xs, vs))
        p.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1" opacity=".75"/>')
        for x, v in zip(xs, vs):
            p.append(f'<circle cx="{x:.1f}" cy="{y0+ph-ph*v/mx:.1f}" r="2.1" fill="{c}"/>')
    for (l, _k), x in zip(steps, xs):
        p.append(f'<text x="{x:.1f}" y="{y0+ph+12}" font-size="8.2" fill="{INK}" text-anchor="middle">{esc(l)}</text>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-24}" font-size="8.4" fill="{MUT}" text-anchor="middle">Pipeline step</text>')
    p.append(f'<text x="9" y="{y0+ph/2:.1f}" font-size="8.4" fill="{MUT}" text-anchor="middle" '
             f'transform="rotate(-90 9 {y0+ph/2:.1f})">Read count</text>')
    lx = x0 + pw + 14
    for si, (name, _vs) in enumerate(data[:18]):
        c = SERIES[si] if len(data) <= len(SERIES) else _pal(si)
        y = y0 + 4 + si * 12.5
        p.append(f'<circle cx="{lx+3}" cy="{y:.1f}" r="2.6" fill="{c}"/>')
        p.append(f'<text x="{lx+10}" y="{y+2.6:.1f}" font-size="8" fill="{INK}">{esc(_short(name, 24))}</text>')
    p.append("</svg>")
    return "".join(p)


def svg_alpha(an, alpha) -> str:
    """Four panels — Observed, Chao1, Shannon, Simpson — one bar per sample."""
    idx = ["Observed", "Chao1", "Shannon", "Simpson"]
    samples = an.samples
    longest = max((len(str(s)) for s in samples), default=4)
    pw_, ph_ = 205, 130
    bot = 16 + min(longest, 20) * 3.8
    cellh = ph_ + bot + 26
    W, H = 520, int(cellh * 2 + 14)
    p = _svg_open(W, H)
    for k, name in enumerate(idx):
        ox = 40 + (k % 2) * 262
        oy = 26 + (k // 2) * cellh
        vals = [alpha[s][name] for s in samples]
        mx = max(vals, default=1) or 1
        if name == "Simpson":
            mx = 1.0
        p.append(f'<text x="{ox-36}" y="{oy-13}" font-size="9.4" fill="{INK}" font-weight="700">{name}</text>')
        for g in (0, .5, 1):
            y = oy + ph_ - ph_ * g
            v = mx * g
            p.append(f'<line x1="{ox}" y1="{y:.1f}" x2="{ox+pw_}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
            lab = f"{v:.2f}" if name in ("Shannon", "Simpson") else fmt(v)
            p.append(f'<text x="{ox-4}" y="{y+2.5:.1f}" font-size="7.6" fill="{FAINT}" text-anchor="end" {MONO}>{lab}</text>')
        slot = pw_ / max(len(samples), 1)
        bw = min(slot * 0.55, 52)
        xs = []
        for i, (s, v) in enumerate(zip(samples, vals)):
            xc = ox + slot * (i + 0.5)
            xs.append(xc)
            h = ph_ * v / mx
            c = SERIES[i] if len(samples) <= len(SERIES) else _pal(i)
            p.append(f'<rect x="{xc-bw/2:.1f}" y="{oy+ph_-h:.1f}" width="{bw:.1f}" height="{h:.1f}" fill="{c}"/>')
            if len(samples) <= 12:
                lab = f"{v:.2f}" if name in ("Shannon", "Simpson") else fmt(v)
                p.append(f'<text x="{xc:.1f}" y="{oy+ph_-h-2:.1f}" font-size="7.2" fill="{MUT}" text-anchor="middle" {MONO}>{lab}</text>')
        _xlabels(p, samples, xs, oy + ph_ + 8, maxlen=20, size=6.8)
    p.append("</svg>")
    return "".join(p)


def svg_rare(an) -> str:
    """Rarefaction: expected features against subsampled depth, one curve per sample."""
    curves = an.rarefaction()
    if not curves:
        return ""
    samples = [s for s in an.samples if s in curves]
    mx_x = max(c[-1][0] for c in curves.values()) or 1
    mx_y = max(max(v for _d, v in c) for c in curves.values()) or 1
    W, H, x0, y0, pw, ph = 520, 300, 50, 12, 320, 244
    p = _svg_open(W, H)
    for k in range(5):
        y = y0 + ph - ph * k / 4
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+pw}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
        p.append(f'<text x="{x0-4}" y="{y+2.5:.1f}" font-size="8" fill="{FAINT}" text-anchor="end" {MONO}>{fmt(mx_y*k/4)}</text>')
        x = x0 + pw * k / 4
        p.append(f'<text x="{x:.1f}" y="{y0+ph+12}" font-size="8" fill="{FAINT}" text-anchor="middle" {MONO}>{fmt(mx_x*k/4)}</text>')
    for i, s in enumerate(samples):
        c = SERIES[i] if len(samples) <= len(SERIES) else _pal(i)
        pts = " ".join(f"{x0+pw*d/mx_x:.1f},{y0+ph-ph*v/mx_y:.1f}" for d, v in curves[s])
        p.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1.5"/>')
        d, v = curves[s][-1]
        p.append(f'<circle cx="{x0+pw*d/mx_x:.1f}" cy="{y0+ph-ph*v/mx_y:.1f}" r="2" fill="{c}"/>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-14}" font-size="8.4" fill="{MUT}" text-anchor="middle">Sequencing depth (reads)</text>')
    p.append(f'<text x="9" y="{y0+ph/2:.1f}" font-size="8.4" fill="{MUT}" text-anchor="middle" '
             f'transform="rotate(-90 9 {y0+ph/2:.1f})">Observed {"taxa" if an.run.emu else "ASVs"}</text>')
    lx = x0 + pw + 14
    for i, s in enumerate(samples[:20]):
        c = SERIES[i] if len(samples) <= len(SERIES) else _pal(i)
        y = y0 + 6 + i * 12.5
        p.append(f'<line x1="{lx}" y1="{y:.1f}" x2="{lx+10}" y2="{y:.1f}" stroke="{c}" stroke-width="2"/>')
        p.append(f'<text x="{lx+14}" y="{y+2.6:.1f}" font-size="8" fill="{INK}">{esc(_short(s, 24))}</text>')
    p.append("</svg>")
    return "".join(p)


def _group_colors(an):
    groups = [an.groups.get(s) for s in an.samples]
    uniq = []
    for g in groups:
        if g and g not in uniq:
            uniq.append(g)
    use_groups = bool(uniq) and len(uniq) < len(an.samples)
    if use_groups:
        cmap = {g: _pal(i) for i, g in enumerate(uniq)}
        return [cmap.get(g, FAINT) for g in groups], cmap
    return [SERIES[i] if len(an.samples) <= len(SERIES) else _pal(i)
            for i in range(len(an.samples))], None


def svg_pcoa(an) -> str:
    D = an.bray()
    if len(D) < 3:
        return ""
    coords, pct = _pcoa(D)
    cols, cmap = _group_colors(an)
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    def rng(v):
        lo, hi = min(v), max(v)
        pad = (hi - lo) * 0.18 or 0.1
        return lo - pad, hi + pad
    (xl, xh), (yl, yh) = rng(xs), rng(ys)
    W, H, x0, y0, pw, ph = 520, 300, 46, 14, 330, 240
    p = _svg_open(W, H)
    p.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="none" stroke="{RULE}"/>')
    X = lambda v: x0 + pw * (v - xl) / (xh - xl)
    Y = lambda v: y0 + ph - ph * (v - yl) / (yh - yl)
    if xl < 0 < xh:
        p.append(f'<line x1="{X(0):.1f}" y1="{y0}" x2="{X(0):.1f}" y2="{y0+ph}" stroke="{RULE}" stroke-dasharray="3 2"/>')
    if yl < 0 < yh:
        p.append(f'<line x1="{x0}" y1="{Y(0):.1f}" x2="{x0+pw}" y2="{Y(0):.1f}" stroke="{RULE}" stroke-dasharray="3 2"/>')
    for k in range(5):
        vx = xl + (xh - xl) * k / 4
        vy = yl + (yh - yl) * k / 4
        p.append(f'<text x="{X(vx):.1f}" y="{y0+ph+11}" font-size="7.6" fill="{FAINT}" text-anchor="middle" {MONO}>{vx:.2f}</text>')
        p.append(f'<text x="{x0-4}" y="{Y(vy)+2.5:.1f}" font-size="7.6" fill="{FAINT}" text-anchor="end" {MONO}>{vy:.2f}</text>')
    placed = [(X(cx) - 5, Y(cy) + 4, 10) for cx, cy in coords]     # points are obstacles too
    for s, (cx, cy), c in zip(an.samples, coords, cols):
        p.append(f'<circle cx="{X(cx):.1f}" cy="{Y(cy):.1f}" r="4.2" fill="{c}" stroke="#fff" stroke-width=".8"/>')
        anchor = "end" if X(cx) > x0 + pw * 0.7 else "start"
        dx = -6 if anchor == "end" else 6
        lx_, ly_ = X(cx) + dx, Y(cy) - 5
        w_ = len(_short(s, 24)) * 3.7
        x1_ = lx_ - w_ if anchor == "end" else lx_
        for _ in range(12):      # step a label down until it clears the ones already placed
            if not any(abs(ly_ - py) < 9 and x1_ < px + pw_ and px < x1_ + w_ for px, py, pw_ in placed):
                break
            ly_ += 10
        placed.append((x1_, ly_, w_))
        p.append(f'<text x="{lx_:.1f}" y="{ly_:.1f}" font-size="8" fill="{INK}" text-anchor="{anchor}">{esc(_short(s, 24))}</text>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-14}" font-size="8.6" fill="{MUT}" text-anchor="middle">PCoA 1 ({pct[0]:.1f}%)</text>')
    p.append(f'<text x="10" y="{y0+ph/2:.1f}" font-size="8.6" fill="{MUT}" text-anchor="middle" '
             f'transform="rotate(-90 10 {y0+ph/2:.1f})">PCoA 2 ({pct[1]:.1f}%)</text>')
    lx = x0 + pw + 14
    items = list(cmap.items()) if cmap else list(zip(an.samples, cols))
    p.append(f'<text x="{lx}" y="{y0+4}" font-size="8.2" fill="{INK}" font-weight="700">{"Group" if cmap else "Sample"}</text>')
    for i, (nm, c) in enumerate(items[:20]):
        y = y0 + 16 + i * 12.5
        p.append(f'<circle cx="{lx+3}" cy="{y:.1f}" r="3" fill="{c}"/>')
        p.append(f'<text x="{lx+10}" y="{y+2.6:.1f}" font-size="8" fill="{INK}">{esc(_short(nm, 24))}</text>')
    p.append("</svg>")
    an.pcoa_pct = pct
    return "".join(p)


def svg_upgma(an) -> str:
    D = an.bray()
    if len(D) < 2:
        return ""
    tree = _upgma(D)
    order = _leaves(tree)
    rowh = 22
    W = 520
    x0, pw = 20, 320
    y0 = 10
    H = int(y0 + rowh * len(order) + 34)
    hmax = _height(tree) or 1
    p = _svg_open(W, H)
    pos = {k: y0 + (order.index(k) + 0.5) * rowh for k in order}
    _pp, _hh, segs = _dendro(tree, pos)
    X = lambda h: x0 + pw - pw * h / hmax
    for (p1, h1), (p2, h2) in segs:
        p.append(f'<line x1="{X(h1):.1f}" y1="{p1:.1f}" x2="{X(h2):.1f}" y2="{p2:.1f}" stroke="{INK}" stroke-width="1"/>')
    cols, _cm = _group_colors(an)
    for k in order:
        p.append(f'<circle cx="{x0+pw:.1f}" cy="{pos[k]:.1f}" r="2.6" fill="{cols[k]}"/>')
        p.append(f'<text x="{x0+pw+7:.1f}" y="{pos[k]+2.8:.1f}" font-size="8.6" fill="{INK}">{esc(an.samples[k])}</text>')
    ya = y0 + rowh * len(order) + 6
    p.append(f'<line x1="{x0}" y1="{ya}" x2="{x0+pw}" y2="{ya}" stroke="{FAINT}" stroke-width=".7"/>')
    for k in range(5):
        h = hmax * k / 4
        p.append(f'<line x1="{X(h):.1f}" y1="{ya}" x2="{X(h):.1f}" y2="{ya+3}" stroke="{FAINT}"/>')
        p.append(f'<text x="{X(h):.1f}" y="{ya+11}" font-size="7.6" fill="{FAINT}" text-anchor="middle" {MONO}>{h:.2f}</text>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-3}" font-size="8.2" fill="{MUT}" text-anchor="middle">Bray–Curtis dissimilarity (average linkage)</text>')
    p.append("</svg>")
    return "".join(p)


def svg_bc_heat(an) -> str:
    D = an.bray()
    n = len(D)
    if n < 2:
        return ""
    order = _leaves(_upgma(D)) if n > 2 else list(range(n))
    longest = max((len(str(s)) for s in an.samples), default=4)
    lab = min(longest, 22) * 4.8 + 8
    cs = min(84, (520 - lab - 70) / n)
    x0, y0 = lab, 10
    W = 520
    H = int(y0 + cs * n + 16 + min(longest, 22) * 4.0 + 12)
    p = _svg_open(W, H)
    for i, a in enumerate(order):
        p.append(f'<text x="{x0-5:.1f}" y="{y0+(i+0.5)*cs+2.6:.1f}" font-size="8" fill="{INK}" text-anchor="end">{esc(_short(an.samples[a], 22))}</text>')
        for j, b in enumerate(order):
            v = D[a][b]
            c = _lin_col(v) if a != b else "#FFFFFF"
            p.append(f'<rect x="{x0+j*cs:.1f}" y="{y0+i*cs:.1f}" width="{cs:.1f}" height="{cs:.1f}" fill="{c}" stroke="#fff" stroke-width=".8"/>')
            if cs >= 18:
                dark = a != b and int(c[1:3], 16) < 0x80
                p.append(f'<text x="{x0+(j+0.5)*cs:.1f}" y="{y0+(i+0.5)*cs+2.6:.1f}" font-size="{min(9.5, cs/3.6):.1f}" '
                         f'fill="{"#fff" if dark else INK}" text-anchor="middle" {MONO}>{_bcf(v)}</text>')
    _xlabels(p, [an.samples[k] for k in order], [x0 + (j + 0.5) * cs for j in range(n)], y0 + cs * n + 9, size=7.2)
    # key
    kx = x0 + cs * n + 16
    for k in range(6):
        v = k / 5
        p.append(f'<rect x="{kx}" y="{y0+k*14}" width="10" height="14" fill="{_lin_col(v)}" stroke="{RULE}" stroke-width=".4"/>')
        p.append(f'<text x="{kx+14}" y="{y0+k*14+9}" font-size="7.6" fill="{FAINT}" {MONO}>{v:.1f}</text>')
    p.append("</svg>")
    return "".join(p)


def _lin_col(f):
    f = max(0.0, min(1.0, f))
    a, mid, b = (0xF4, 0xF7, 0xF8), (0x5B, 0x9B, 0xC8), (0x0B, 0x2E, 0x59)
    c = ([a[i] + (mid[i] - a[i]) * f * 2 for i in range(3)] if f < 0.5
         else [mid[i] + (b[i] - mid[i]) * (f - 0.5) * 2 for i in range(3)])
    return "#%02X%02X%02X" % tuple(int(x) for x in c)


REFERENCES = [
    ("Martin, M. (2011).", "Cutadapt removes adapter sequences from high-throughput sequencing reads.",
     "EMBnet.journal, 17(1), 10–12."),
    ("Callahan, B. J., et al. (2016).", "DADA2: High-resolution sample inference from Illumina amplicon data.",
     "Nature Methods, 13, 581–583."),
    ("Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007).",
     "Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy.",
     "Applied and Environmental Microbiology, 73(16), 5261–5267."),
    ("Quast, C., et al. (2013).", "The SILVA ribosomal RNA gene database project.",
     "Nucleic Acids Research, 41, D590–D596."),
    ("Curry, K. D., et al. (2022).", "Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data.",
     "Nature Methods, 19, 845–853."),
    ("De Coster, W., & Rademakers, R. (2023).", "NanoPack2: population-scale evaluation of long-read sequencing data.",
     "Bioinformatics, 39(5), btad311."),
    ("Rognes, T., et al. (2016).", "VSEARCH: a versatile open source tool for metagenomics.", "PeerJ, 4, e2584."),
    ("Li, H. (2018).", "Minimap2: pairwise alignment for nucleotide sequences.", "Bioinformatics, 34(18), 3094–3100."),
    ("McMurdie, P. J., & Holmes, S. (2013).", "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data.",
     "PLoS ONE, 8(4), e61217."),
    ("Oksanen, J., et al. (2022).", "vegan: Community Ecology Package.", "R package version 2.6."),
    ("Shannon, C. E. (1948).", "A mathematical theory of communication.", "Bell System Technical Journal, 27, 379–423."),
    ("Simpson, E. H. (1949).", "Measurement of diversity.", "Nature, 163, 688."),
    ("Chao, A. (1984).", "Nonparametric estimation of the number of classes in a population.",
     "Scandinavian Journal of Statistics, 11, 265–270."),
    ("Hurlbert, S. H. (1971).", "The nonconcept of species diversity: a critique and alternative parameters.",
     "Ecology, 52(4), 577–586."),
    ("Bray, J. R., & Curtis, J. T. (1957).", "An ordination of the upland forest communities of southern Wisconsin.",
     "Ecological Monographs, 27(4), 325–349."),
    ("Lozupone, C., & Knight, R. (2008).", "Species divergence and the measurement of microbial diversity.",
     "FEMS Microbiology Reviews, 32(4), 557–578."),
    ("Salter, S. J., et al. (2014).", "Reagent and laboratory contamination can critically impact sequence-based microbiome analyses.",
     "BMC Biology, 12, 87."),
]
REF_DADA2 = [0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16]
REF_EMU = [5, 0, 6, 7, 4, 8, 9, 10, 11, 12, 13, 14, 15, 16]


def _fig(svg, caption, folder=None, n=None):
    if not svg:
        return ""
    lab = f"<b>Figure {n}.</b> " if n else ""
    fol = f'<div class="folder">Folder: {folder}</div>' if folder else ""
    return f"<figure>{svg}<figcaption>{lab}{caption}</figcaption>{fol}</figure>"


def _pct(v, d=1):
    return f"{v:.{d}f}%"


def _range_txt(pairs, f=fmt):
    """'from X (A) to Y (B)' over [(sample, value)]."""
    lo = min(pairs, key=lambda kv: kv[1])
    hi = max(pairs, key=lambda kv: kv[1])
    return f"from {f(lo[1])} ({esc(lo[0])}) to {f(hi[1])} ({esc(hi[0])})"


def _skipped(run):
    """Which optional Emu steps the run log says were skipped: {"qc": bool, "chimera": bool}.
    None for a step when there is no log to tell (the report then does not guess)."""
    lp = run.path / "pipeline.log"
    if not lp.exists():
        return {"qc": None, "chimera": None, "chimera_off": False}
    t = lp.read_text(encoding="utf-8", errors="replace")
    off = bool(re.search(r"Chimera removal off", t))
    return {"qc": bool(re.search(r"chopper not found|QC filtering skipped", t)),
            "chimera": bool(off or re.search(r"vsearch not found|Chimera removal skipped|"
                                             r"not found \u2014 skipping chimera removal", t)),
            "chimera_off": off}


def _methods(run, an):
    p, out = run.params, []
    pf, pr = p.get("primer_f") or "", p.get("primer_r") or ""
    if run.emu:
        sk = _skipped(run)
        reg = run.summary.get("region") or p.get("ont_region") or "V1-V9"
        lo, hi = {"V1-V9": (1000, 1700), "V7-V8": (250, 450)}.get(reg, (50, 2000))
        out.append("<h3>Read quality filtering</h3>")
        if sk["qc"]:
            out.append("<p>Raw Oxford Nanopore reads were demultiplexed by MinKNOW; reads from the same "
                       "sample spread over several barcode files were combined before analysis. No "
                       "additional read-quality filter was applied in this run (chopper was not "
                       "available); reads were used as passed by MinKNOW, and length was controlled at "
                       "the primer-trimming step below.</p>")
        else:
            out.append("<p>Raw Oxford Nanopore reads were demultiplexed by MinKNOW; reads from the same "
                       "sample spread over several barcode files were combined before analysis. Reads were "
                       "filtered with chopper (NanoPack2) to a mean quality of Q10 or better and a length "
                       "of 100&ndash;2,000 bases.</p>")
            out.append('<p class="cmd">chopper -q 10 -l 100 --maxlength 2000</p>')
        out.append("<h3>Primer trimming</h3>")
        if pf or pr:
            out.append(f"<p>Primer sequences were located and removed with cutadapt, searching both "
                       f"strands. The forward primer was trimmed from the 5&prime; end and the reverse "
                       f"complement of the reverse primer from the 3&prime; end; reads in which no primer "
                       f"was found, or that fell outside {lo:,}&ndash;{hi:,} bp after trimming (the "
                       f"expected length of the {esc(reg)} amplicon), were discarded.</p>")
            out.append("<ul>" + (f"<li>Forward primer: <b class=\"mono\">{esc(pf)}</b></li>" if pf else "")
                       + (f"<li>Reverse primer: <b class=\"mono\">{esc(pr)}</b></li>" if pr else "") + "</ul>")
            out.append(f'<p class="cmd">cutadapt -e 0.2 -O 10 -g {esc(pf)} -a {esc(_rc(pr))} --revcomp '
                       f'-n 2 --discard-untrimmed -m {lo} -M {hi}</p>')
        else:
            out.append("<p>No primer sequences were given, so reads were analysed as delivered.</p>")
        out.append("<h3>Chimera removal</h3>")
        out.append("<p>No separate chimera filter was applied, following the Emu method (Curry et "
                   "al. 2022). Reference-based chimera detection (UCHIME) was evaluated and not used: "
                   "on full-length nanopore reads it flags genuine reads of taxa that are sparsely "
                   "represented in the reference as chimeras, removing them selectively and distorting "
                   "relative abundance. Emu&rsquo;s expectation&ndash;maximisation and its minimum "
                   "abundance threshold limit the influence of the few true chimeras.</p>"
                   if sk.get("chimera_off") else
                   "<p>Chimera removal with VSEARCH was not performed in this run (VSEARCH was not "
                   "available). Emu&rsquo;s alignment-based abundance estimation is relatively tolerant "
                   "of chimeras, but low-abundance taxa should be read with that in mind.</p>"
                   if sk["chimera"] else
                   "<p>Chimeric reads were removed with VSEARCH in reference mode "
                   "(<span class=\"mono\">--uchime_ref</span>) against the sequences of the "
                   "classification database.</p>")
        out.append("<h3>Taxonomic classification</h3>")
        out.append(f"<p>Each full-length read was classified with Emu, which aligns reads to the "
                   f"reference with minimap2 and estimates the relative abundance of every reference "
                   f"taxon by expectation&ndash;maximisation, so a read matching several close relatives "
                   f"is shared between them in proportion to the evidence. The reference was "
                   f"<b>{esc(run.db_name)}</b>; taxa below a relative abundance of "
                   f"{esc(run.summary.get('emu_min_abundance') or p.get('ont_min_abundance') or 0.0001)} "
                   f"were dropped. "
                   + ("Read counts in this report are Emu's estimated read counts per taxon."
                      if run.emu_counts == "estimated_reads" else
                      "Read counts in this report are derived from Emu's relative abundance "
                      "&times; each sample's classified reads.") + "</p>")
        out.append('<p class="cmd">emu abundance --type map-ont --db ' + esc(run.db_name)
                   + f' --min-abundance {esc(run.summary.get("emu_min_abundance") or 0.0001)}'
                   + (' --keep-counts' if run.emu_counts == "estimated_reads" else "") + '</p>')
    else:
        paired = any("denoisedR" in t for t in run.tracking[:1])
        e = p.get("cutadaptErrorRate", 0.1)
        o = p.get("cutadaptOverlap", 3)
        disc = " --discard-untrimmed" if p.get("discardUntrimmed") else ""
        out.append("<h3>Primer trimming</h3>")
        if pf or pr:
            out.append(f"<p>Paired-end reads were primer-trimmed with cutadapt before denoising. The "
                       f"forward primer was removed from the 5&prime; end of R1 and the reverse primer "
                       f"from the 5&prime; end of R2, with the reverse complements removed from the "
                       f"3&prime; ends where the read ran through the amplicon; an error rate of {e} and a "
                       f"minimum overlap of {o} bases were allowed, and reads shorter than 50 bp after "
                       f"trimming were discarded" + ("; reads with no primer were discarded." if disc else ".")
                       + "</p>" if paired else
                       f"<p>Reads were primer-trimmed with cutadapt (error rate {e}, minimum overlap {o}), "
                       f"and reads shorter than 50 bp after trimming were discarded.</p>")
            out.append("<ul>" + (f"<li>Forward primer: <b class=\"mono\">{esc(pf)}</b></li>" if pf else "")
                       + (f"<li>Reverse primer: <b class=\"mono\">{esc(pr)}</b></li>" if pr else "") + "</ul>")
            cmd = (f"cutadapt -g {pf} -a {_rc(pr)} -G {pr} -A {_rc(pf)} -e {e} -O {o} -m 50{disc} "
                   f"-o R1_cut.fastq.gz -p R2_cut.fastq.gz R1.fastq.gz R2.fastq.gz" if paired else
                   f"cutadapt -g {pf} -a {_rc(pr)} -e {e} -O {o} -m 50{disc} -o cut.fastq.gz in.fastq.gz")
            out.append(f'<p class="cmd">{esc(cmd)}</p>')
        else:
            out.append("<p>No primer sequences were given, so primer trimming was skipped.</p>")
        out.append("<h3>Denoising and ASV inference</h3>")
        out.append(f"<p>Trimmed reads were processed with DADA2. Quality filtering used "
                   f"<span class=\"mono\">filterAndTrim()</span> with truncation at {run.trunc} bases "
                   f"(forward / reverse), maxEE = {p.get('maxEE_F')} / {p.get('maxEE_R')}, truncQ = 2, "
                   f"maxN = 0 and PhiX removal. Error models were learned separately for forward and "
                   f"reverse reads ({fmt(p.get('nbases') or 1e8)} bases), and exact amplicon sequence "
                   f"variants (ASVs) were inferred with <span class=\"mono\">dada()</span> "
                   f"(pool = {esc(p.get('pool') or 'FALSE')}). Denoised pairs were merged with "
                   f"<span class=\"mono\">mergePairs()</span> (minimum overlap 12 bp, no mismatches), and "
                   f"chimeras were removed with <span class=\"mono\">removeBimeraDenovo()</span> "
                   f"(method = {esc(p.get('chimeraMethod') or 'consensus')}). The non-chimeric ASV table "
                   f"was used for everything that follows.</p>")
        out.append("<h3>Taxonomic assignment</h3>")
        out.append(f"<p>ASVs were classified with the DADA2 implementation of the naive Bayesian "
                   f"classifier (<span class=\"mono\">assignTaxonomy()</span>) against "
                   f"<b>{esc(run.db_name)}</b>, with a minimum bootstrap confidence of "
                   f"{esc(p.get('minBoot', 50))}. Species were added where the reference carries them, by "
                   f"exact matching (<span class=\"mono\">addSpecies()</span>). Organelle sequences "
                   f"(mitochondria, chloroplasts) and eukaryotic ASVs were identified from the assigned "
                   f"lineage and set aside from the community figures.</p>")
    out.append("<h3>Diversity analysis</h3>")
    unit = "taxa" if run.emu else "ASVs"
    out.append(f"<p>Alpha diversity (Observed {unit}, Chao1, ACE, Shannon and Simpson) was computed "
               f"from the per-sample counts. Rarefaction curves give the expected number of {unit} in a "
               f"random subsample of each depth (Hurlbert&rsquo;s formula, as in vegan&rsquo;s "
               f"<span class=\"mono\">rarefy()</span>). Beta diversity was measured as Bray&ndash;Curtis "
               f"dissimilarity on relative abundances and visualised by principal coordinates analysis "
               f"(PCoA), average-linkage (UPGMA) clustering and a dissimilarity heatmap. Taxonomic "
               f"composition is shown as the top 30 taxa at each rank per sample; the remaining taxa are "
               f"pooled as &ldquo;Other&rdquo; and reads without an assignment at that rank as "
               f"&ldquo;Unclassified&rdquo;."
               + (" Because there is one sample per condition, comparisons are descriptive; no "
                  "statistical test between groups is possible without replicates."
                  if len(an.samples) > 1 and an.groups and len(set(an.groups.values())) >= len(an.samples)
                  else " PERMANOVA (vegan adonis2) is reported where groups have replicates."
                  if an.groups and len(set(an.groups.values())) < len(an.samples) else "")
               + "</p>")
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  PICRUSt2 (predicted function) — read from <run>/PICRUSt2 when it exists
# ══════════════════════════════════════════════════════════════════════════════
REF_PICRUST = [
    ("Douglas, G. M., et al. (2020).", "PICRUSt2 for prediction of metagenome functions.",
     "Nature Biotechnology, 38, 685–688."),
    ("Barbera, P., et al. (2019).", "EPA-ng: massively parallel evolutionary placement of genetic sequences.",
     "Systematic Biology, 68(2), 365–369."),
    ("Czech, L., Barbera, P., & Stamatakis, A. (2020).", "Genesis and Gappa: processing, analyzing and "
     "visualizing phylogenetic (placement) data.", "Bioinformatics, 36(10), 3263–3265."),
    ("Louca, S., & Doebeli, M. (2018).", "Efficient comparative phylogenetics on large trees.",
     "Bioinformatics, 34(6), 1053–1055."),
    ("Ye, Y., & Doak, T. G. (2009).", "A parsimony approach to biological pathway reconstruction/inference "
     "for genomes and metagenomes.", "PLoS Computational Biology, 5(8), e1000465."),
    ("Caspi, R., et al. (2020).", "The MetaCyc database of metabolic pathways and enzymes — a 2019 update.",
     "Nucleic Acids Research, 48, D445–D453."),
    ("Kanehisa, M., & Goto, S. (2000).", "KEGG: Kyoto Encyclopedia of Genes and Genomes.",
     "Nucleic Acids Research, 28(1), 27–30."),
]


def _read_tsv_any(p: Path):
    """Rows of a PICRUSt2 table (.tsv or .tsv.gz) as lists of strings; [] when absent."""
    import gzip
    for f in (p, p.with_suffix("") if p.suffix == ".gz" else Path(str(p) + ".gz")):
        if f.exists():
            op = gzip.open if f.suffix == ".gz" else open
            with op(f, "rt", encoding="utf-8", errors="replace") as fh:
                return [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    return []


def load_picrust(run, samples):
    """PICRUSt2 results of a run, or None. Pathway abundances are turned into % of each
    sample's total predicted pathway abundance."""
    d = run.path / "PICRUSt2"
    pw = _read_tsv_any(d / "pathways_out" / "path_abun_unstrat_descrip.tsv.gz")
    has_desc = bool(pw)
    if not pw:
        pw = _read_tsv_any(d / "pathways_out" / "path_abun_unstrat.tsv.gz")
    if len(pw) < 2:
        return None
    head = pw[0]
    c0 = 2 if has_desc else 1
    cols = head[c0:]
    order = [s for s in samples if s in cols] + [s for s in cols if s not in samples]
    idx = {s: cols.index(s) + c0 for s in order}
    raw = {}
    for r in pw[1:]:
        try:
            raw[r[0]] = ((r[1] if has_desc else ""), {s: float(r[idx[s]] or 0) for s in order})
        except (ValueError, IndexError):
            continue
    tot = {s: sum(v[s] for _d, v in raw.values()) or 1.0 for s in order}
    pct = {k: (dsc, {s: 100 * v[s] / tot[s] for s in order}) for k, (dsc, v) in raw.items()}
    nsti = {}
    for sub in ("KO_metagenome_out", "EC_metagenome_out"):
        rows = _read_tsv_any(d / sub / "weighted_nsti.tsv.gz")
        for r in rows[1:]:
            try:
                nsti[r[0]] = float(r[1])
            except (ValueError, IndexError):
                pass
        if nsti:
            break

    def n_feat(sub):
        rows = _read_tsv_any(d / sub / "pred_metagenome_unstrat.tsv.gz")
        return max(len(rows) - 1, 0) if rows else None

    summ = {}
    try:
        summ = json.loads((d / "picrust2_summary.json").read_text())
    except Exception:
        pass
    return {"samples": order, "pct": pct, "n_path": len(raw), "nsti": nsti,
            "n_ko": n_feat("KO_metagenome_out"), "n_ec": n_feat("EC_metagenome_out"),
            "version": summ.get("picrust2_version") or "", "summary": summ}


def svg_pathway_heat(pc, n=30) -> str:
    """Top-n MetaCyc pathways (by mean %) × samples, ordered by mean abundance."""
    S = pc["samples"]
    keep = sorted(pc["pct"], key=lambda k: -sum(pc["pct"][k][1].values()))[:n]
    if not keep:
        return ""
    # Every pathway sits within a narrow band (≈0.5–2%), so one colour scale for the whole
    # map is near-uniform; colour is scaled within each row (lowest → highest sample).
    rng_ = {}
    for k in keep:
        vs = [pc["pct"][k][1][s] for s in S]
        rng_[k] = (min(vs), max(vs))
    W = 560
    lab_w = 262
    x0, y0 = 46, 20                 # left margin: rotated sample names run left of their column
    cw = min(60, (W - lab_w - x0 - 6) / max(len(S), 1))
    ch = 14.0
    longest = max((len(str(s)) for s in S), default=4)
    H = int(y0 + ch * len(keep) + 22 + min(longest, 22) * 4.0 + 38)
    p = _svg_open(W, H)
    p.append(f'<text x="0" y="11" font-size="9.4" fill="{MUT}">Top {len(keep)} MetaCyc pathways '
             f'— % of predicted pathway abundance (colour scaled per row)</text>')
    for j, s in enumerate(S):
        for i, k in enumerate(keep):
            v = pc["pct"][k][1][s]
            x, y = x0 + j * cw, y0 + i * ch
            lo_, hi_ = rng_[k]
            f_ = (v - lo_) / (hi_ - lo_) if hi_ > lo_ else 0.5
            c = _heat_col(0.01 * 10 ** (3 * f_) - 0.01 + 0.001, 9.99)   # f 0→1 maps to light→dark
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="{c}" '
                     f'stroke="#fff" stroke-width=".5"/>')
            if cw >= 24:
                rgb = [int(c[i:i+2], 16) for i in (1, 3, 5)]
                dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 165
                txt = f"{v:.2f}" if v >= 0.005 else ("0" if v == 0 else "&lt;.01")
                p.append(f'<text x="{x+cw/2:.1f}" y="{y+ch*0.72:.1f}" font-size="7" '
                         f'fill="{"#fff" if dark else INK}" text-anchor="middle" {MONO}>{txt}</text>')
    lx = x0 + cw * len(S) + 5
    for i, k in enumerate(keep):
        dsc = pc["pct"][k][0] or k
        p.append(f'<text x="{lx:.1f}" y="{y0+i*ch+ch*0.74:.1f}" font-size="7.6" fill="{INK}">'
                 f'{esc(_short(dsc, 62))}</text>')
    yb = y0 + ch * len(keep)
    _xlabels(p, S, [x0 + (j + 0.5) * cw for j in range(len(S))], yb + 9)
    ky = H - 14
    p.append(f'<text x="0" y="{ky-4}" font-size="7.6" fill="{MUT}">Within each pathway</text>')
    for k, (f_, lab) in enumerate(((0, "lowest"), (0.25, ""), (0.5, ""), (0.75, ""), (1, "highest"))):
        c = _heat_col(0.01 * 10 ** (3 * f_) - 0.01 + 0.001, 9.99)
        p.append(f'<rect x="{k*34}" y="{ky}" width="34" height="7" fill="{c}" '
                 f'stroke="{RULE}" stroke-width=".4"/>')
        if lab:
            p.append(f'<text x="{k*34 + 17}" y="{ky+14}" font-size="7.2" fill="{FAINT}" '
                     f'text-anchor="middle">{lab}</text>')
    p.append("</svg>")
    return "".join(p)


def _nsti_word(v):
    return ("close to reference genomes" if v <= 0.06 else
            "acceptable" if v <= 0.15 else "distant from reference genomes — read with caution")


def _picrust_methods(pc):
    ver = f" v{esc(pc['version'])}" if pc.get("version") else ""
    feats = " and ".join(x for x in (("KEGG orthologs (KO)" if pc.get("n_ko") is not None else ""),
                                     ("Enzyme Commission numbers (EC)" if pc.get("n_ec") is not None else ""))
                         if x) or "gene families"
    out = ["<h3>Functional prediction (PICRUSt2)</h3>",
           f"<p>The functional potential of each community was predicted from the ASVs with "
           f"PICRUSt2{ver} (Douglas et al. 2020). ASV sequences were aligned to the PICRUSt2 reference "
           f"(HMMER) and placed into its reference phylogeny of bacterial genomes with EPA-ng and gappa. "
           f"For every ASV, the 16S rRNA gene copy number and the copy numbers of {feats} were inferred "
           f"from its placement by hidden-state prediction (castor), and the nearest sequenced taxon "
           f"index (NSTI) was recorded as its distance to the closest reference genome. ASV abundances "
           f"were divided by their predicted 16S copy number and multiplied by the predicted gene copy "
           f"numbers to give a predicted metagenome per sample; ASVs with NSTI above 2 were excluded. "
           f"MetaCyc pathway abundances were inferred from the predicted EC profile with MinPath "
           f"(pathway_pipeline), and descriptions were added from KEGG and MetaCyc.</p>",
           '<p class="cmd">place_seqs.py -s study_seqs.fna -o placed_seqs.tre<br>'
           'hsp.py -i 16S -t placed_seqs.tre -o marker_predicted_and_nsti.tsv.gz -n<br>'
           'hsp.py -i EC|KO -t placed_seqs.tre -o EC|KO_predicted.tsv.gz<br>'
           'metagenome_pipeline.py -i asv_table.tsv -m marker_predicted_and_nsti.tsv.gz '
           '-f EC|KO_predicted.tsv.gz -o EC|KO_metagenome_out --max_nsti 2.0<br>'
           'pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out<br>'
           'add_descriptions.py -m KO|EC|METACYC</p>']
    return out


def _picrust_result(pc, F, T, top=25):
    S = pc["samples"]
    out = ['<h3 class="pbreak">Predicted functional profile (PICRUSt2)</h3>']
    nz = pc["nsti"]
    parts = []
    if pc.get("n_ko") is not None:
        parts.append(f"{fmt(pc['n_ko'])} KEGG orthologs")
    if pc.get("n_ec") is not None:
        parts.append(f"{fmt(pc['n_ec'])} EC numbers")
    parts.append(f"{fmt(pc['n_path'])} MetaCyc pathways")
    txt = f"<p>PICRUSt2 predicted {', '.join(parts[:-1]) + ' and ' if len(parts) > 1 else ''}{parts[-1]} across the {len(S)} sample{'s' if len(S) != 1 else ''}. "
    if nz:
        vals = [(s, nz[s]) for s in S if s in nz] or list(nz.items())
        hi = max(v for _s, v in vals)
        txt += (f"Weighted NSTI ranged {_range_txt(vals, lambda v: f'{v:.3f}')}; "
                + ("every sample is below 0.15, so the ASVs are close enough to sequenced genomes for the "
                   "prediction to be informative. " if hi <= 0.15 else
                   "samples above 0.15 contain many ASVs distant from any sequenced genome, and their "
                   "predictions are less reliable. "))
    out.append(txt + "</p>")
    if nz:
        out.append(T(["Sample", "Weighted NSTI", "Assessment"],
                     [[esc(s), f"{nz[s]:.3f}", _nsti_word(nz[s])] for s in S if s in nz],
                     "Weighted nearest sequenced taxon index per sample: the mean branch distance of the "
                     "sample&rsquo;s reads to the closest reference genome. Lower is better; values above "
                     "0.15 call for caution."))
    keep = sorted(pc["pct"], key=lambda k: -sum(pc["pct"][k][1].values()))[:top]
    per = len(S) <= 6
    rows = []
    for k in keep:
        dsc, v = pc["pct"][k]
        mean = sum(v.values()) / len(S)
        rows.append([f'<div style="text-align:left">{esc(dsc or k)}</div>',
                     f'<span class="mono" style="font-size:6.8pt">{esc(k)}</span>', f"{mean:.2f}"]
                    + ([f"{v[s]:.2f}" for s in S] if per else []))
    out.append(T(["Pathway", "MetaCyc ID", "Mean %"] + ([f"{esc(h)} %" for h in _col_heads(S)] if per else []),
                 rows, f"The {len(keep)} most abundant predicted MetaCyc pathways, as a percentage of each "
                 f"sample&rsquo;s total predicted pathway abundance"
                 + (", with the value for each sample." if per else ".")))
    out.append('<div class="note"><p style="margin:0"><b>These are predictions, not measurements.</b> '
               "PICRUSt2 infers gene content from the closest sequenced genomes of each ASV; it cannot "
               "see strain-level gene gain or loss, plasmids or genes absent from the reference genomes. "
               "Highly conserved housekeeping pathways dominate every profile, so differences between "
               "samples are more informative than the ranking itself. Shotgun metagenomics is needed to "
               "confirm a functional difference.</p></div>")
    out.append(F(svg_pathway_heat(pc, 30), "Predicted abundance of the 30 most abundant MetaCyc pathways "
                 "per sample. Numbers are the percentage of the sample&rsquo;s total predicted pathway "
                 "abundance; colour is scaled within each pathway from the sample with the lowest value "
                 "(light) to the highest (dark), so it shows which sample is enriched, not how large the "
                 "difference is.", "PICRUSt2/pathways_out"))
    return out


def build_profile_html(run, title=None, app="NextGen-Amplicon", company="",
                       subtitle=None, top_asv=30) -> str:
    """Full result report for one run, laid out as Method / Result / References with a
    figure for every analysis: read tracking, composition at every rank, a clustered
    genus heatmap, archaea, alpha diversity and rarefaction, and beta diversity."""
    an = Analysis(run)
    S = an.samples
    unit = "taxa" if run.emu else "ASVs"
    marker = run.summary.get("marker") or run.params.get("marker") or "16S"
    region = run.summary.get("region") or ("V1–V9" if run.emu else "V3–V4"
                                          if (run.params.get("primer_f") or "").startswith("CCTACGG")
                                          else "")
    rep_title = (f"16S {region} Amplicon Sequencing Report".replace("  ", " ")
                 if "16S" in marker.upper() else f"{marker} Amplicon Sequencing Report")
    sample_label = ", ".join(S[:3]) + ("…" if len(S) > 3 else "")
    project = title if title and title.strip() and title.strip() != rep_title else None
    title = rep_title          # the report type heads the page; a user title becomes the project line
    subtitle = subtitle or (f"{run.label} · {len(S)} sample{'s' if len(S) != 1 else ''} · "
                            + ("Oxford Nanopore full-length 16S, classified with Emu" if run.emu
                               else "Illumina paired-end, DADA2 amplicon sequence variants"))
    fig = [0]
    tab = [0]

    def F(svg, cap, folder=None):
        if not svg:
            return ""
        fig[0] += 1
        return _fig(svg, cap, folder, fig[0])

    def T(head, rows, cap):
        tab[0] += 1
        return (f'<p class="tcap"><b>Table {tab[0]}.</b> {cap}</p>' + _tbl(head, rows))

    rt = run.path / "r_tables"
    def folder(*names):
        have = [n for n in names if (run.path / n).exists()]
        return " &middot; ".join(esc(n) for n in have) if have else None

    # totals
    inp = sum(float(t.get("input") or 0) for t in run.tracking) if run.tracking else 0
    fin = sum(float(t.get(run.retained_key) or 0) for t in run.tracking) if run.tracking else run.total
    kreads = (run.qc or {}).get("kingdom_reads") or {}
    arch_rows = an.archaea()
    n_arch = sum(r[1] for r in arch_rows)

    out = [f'<div class="eyebrow">{esc(run.label)} &middot; {esc(marker)} &middot; '
           f'{len(S)} sample{"s" if len(S) != 1 else ""}</div>',
           f"<h1>{esc(title)}</h1>"]
    if project:
        out.append(f'<p class="project">{esc(project)}</p>')
    out.append(f'<p class="sub">{esc(subtitle)}</p>')
    kpis = [(fmt(inp) if inp else "—", "reads sequenced" if run.emu else "read pairs sequenced"),
            (fmt(fin), "reads classified" if run.emu else "non-chimeric reads"),
            (str(run.n_asvs), f"{unit} detected"),
            (f"{run.rank_pct('Genus') or 0}%", "reads with a genus")]
    out.append('<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>'
        for v, l in kpis) + "</div>")
    srows = [[esc(s), esc(an.groups.get(s, "—")), f'<span class="dim">{esc(an.desc.get(s, ""))}</span>']
             for s in S]
    out.append(_tbl(["Sample", "Group", "Description"], srows) if (an.groups or an.desc)
               else _tbl(["Sample", "Reads analysed"], [[esc(s), fmt(an.comm_tot[s])] for s in S]))

    # ═══ METHOD ═══
    out.append("<h2>Method</h2>")
    out += _methods(run, an)
    pic = load_picrust(run, S)
    if pic:
        out += _picrust_methods(pic)
    cfg = ([("Marker", marker)] + run.emu_config()) if run.emu else [
        ("Marker", marker), ("Reference database", run.db_name),
        ("Forward primer", run.params.get("primer_f") or "—"),
        ("Reverse primer", run.params.get("primer_r") or "—"),
        ("truncLen F / R", run.trunc),
        ("maxEE F / R", f'{run.params.get("maxEE_F")} / {run.params.get("maxEE_R")}'),
        ("Pooling", run.params.get("pool") or "—"),
        ("Minimum bootstrap", run.params.get("minBoot")),
        ("Chimera method", run.params.get("chimeraMethod"))]
    out.append(T(["Parameter", "Value"], [[k, esc(v if v is not None else "—")] for k, v in cfg],
                 "Run configuration."))

    # ═══ RESULT ═══
    out.append('<h2 class="pbreak">Result</h2>')
    out.append("<h3>Data quality control</h3>")
    if run.tracking:
        per = []
        for t in run.tracking:
            nm = str(next(iter(t.values())))
            i0 = float(t.get("input") or 0)
            f0 = float(t.get(run.retained_key) or 0)
            per.append((nm, i0, f0, 100 * f0 / i0 if i0 else 0))
        steps = run.stages if run.emu else [("filtering", "filtered"), ("denoising", "denoisedF"),
                                            ("merging", "merged"), ("chimera removal", "nonchim")]
        losses = []
        prev = "input"
        for lbl, k in steps:
            if k not in run.tracking[0]:
                continue
            a = sum(float(t.get(prev) or 0) for t in run.tracking)
            b = sum(float(t.get(k) or 0) for t in run.tracking)
            if a:
                losses.append((lbl, 100 * (a - b) / a))
            prev = k
        worst = max(losses, key=lambda kv: kv[1]) if losses else None
        out.append("<p><b>Overall pipeline outcome.</b> "
                   + (f"All {len(S)} samples were processed through "
                      + ", ".join(x for x, k in (("quality filtering", "qc"), ("primer trimming", None),
                                                 ("chimera removal", "chimera"))
                                  if not (k and _skipped(run)[k]))
                      + ", and the remaining reads were classified with Emu. "
                      if run.emu else
                      f"All {len(S)} paired-end samples were processed through primer trimming, quality "
                      f"filtering, denoising, read-pair merging and chimera removal in DADA2. ")
                   + f"From {fmt(inp)} input {'reads' if run.emu else 'read pairs'}, {fmt(fin)} "
                   f"({100*fin/inp if inp else 0:.1f}%) were retained "
                   + (f"and assigned to {run.n_asvs} taxa. " if run.emu else
                      f"and resolved into {run.n_asvs} ASVs. ")
                   + (f"{fmt(kreads.get('Bacteria', 0))} reads were classified as Bacteria and "
                      f"{fmt(kreads.get('Archaea', 0))} as Archaea. " if kreads else "")
                   + (f"{fmt(run.host + run.euk)} reads ({100*(run.host+run.euk)/run.total:.1f}%) were host "
                      f"organelle or eukaryotic and are excluded from the community figures. "
                      if run.total and (run.host + run.euk) > 0 else
                      "No mitochondrial, chloroplast or eukaryotic reads were detected. ")
                   + "</p>")
        qi = (f"<p><b>QC interpretation.</b> Input depth ranged "
              f"{_range_txt([(n, i) for n, i, _f, _p in per])} per sample. The final retained reads "
              f"ranged {_range_txt([(n, f) for n, _i, f, _p in per])}, a retention of "
              f"{min(x[3] for x in per):.1f}&ndash;{max(x[3] for x in per):.1f}%. ")
        if worst:
            qi += (f"The largest loss across the run was at the {esc(worst[0].lower())} step "
                   f"({worst[1]:.1f}% of the reads entering that step). ")
        low = [x for x in per if x[3] < 50]
        if low:
            qi += (f"{', '.join(esc(x[0]) for x in low)} retained less than half of "
                   f"{'its' if len(low) == 1 else 'their'} input, which deserves a look at the raw read "
                   f"quality before the results are interpreted. ")
        if run.emu:
            def _same(a, b):
                return all(float(t.get(a) or 0) == float(t.get(b) or 0) for t in run.tracking)
            sk = _skipped(run)
            if sk["qc"]:
                qi += ("No read-quality filter was applied (chopper was not available), so the QC "
                       "column equals the input. ")
            elif "qc_filtered" in run.tracking[0] and _same("input", "qc_filtered"):
                qi += ("The quality filter removed no reads: either every read already met Q10 and "
                       "100&ndash;2,000 bp (MinKNOW&rsquo;s own pass filter is Q10), or chopper was not "
                       "found and the step was skipped &mdash; the run log states which. ")
            if sk.get("chimera_off"):
                qi += "No chimera filter was applied (see Method), so the last two columns are equal. "
            elif sk["chimera"]:
                qi += "Chimera removal was not performed (VSEARCH was not available). "
            elif "primer_trimmed" in run.tracking[0] and _same("primer_trimmed", "final"):
                qi += ("Chimera removal removed no reads: either VSEARCH found no chimeras or it was "
                       "not found and the step was skipped &mdash; the run log states which. ")
        fins = [x[2] for x in per if x[2] > 0]
        if fins and max(fins) / min(fins) > 3:
            qi += ("Depth differs more than three-fold between samples, so richness values are best "
                   "compared through the rarefaction curves rather than as raw counts. ")
        elif fins:
            qi += "Depth is reasonably even across samples. "
        out.append(qi + "</p>")

        if run.emu:
            cols = [k for _l, k in run.stages if k in run.tracking[0]]
            heads = [l for l, k in run.stages if k in run.tracking[0]]
        else:
            cols = [c for c in ("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
                    if c in run.tracking[0]]
            heads = cols
        rows = [[esc(next(iter(t.values())))] + [fmt(t.get(c) or 0) for c in cols]
                + [f'{100*float(t.get(run.retained_key) or 0)/(float(t.get("input") or 0) or 1):.1f}%']
                for t in run.tracking]
        out.append(T(["Sample"] + heads + ["Retained"], rows,
                     "Read counts per sample at each processing step"
                     + ((" (raw, QC-filtered, primer-trimmed, classified)." if _skipped(run).get("chimera")
                         else " (raw, QC-filtered, primer-trimmed, non-chimeric).") if run.emu else
                        " (raw, filtered, denoisedF, denoisedR, merged, non-chimeric).")))
        out.append(F(svg_steps(run, an), "Read count distribution per pipeline step. Each line is one "
                     "sample" + ("; boxes show the spread across samples." if len(S) >= 4 else "."),
                     folder("read_tracking.csv")))
        out.append(F(svg_readflow(run), "Reads surviving each stage across the whole run, as a share "
                     "of input; the right-hand figure is what that stage cost.", None))

    o = (run.orientation or {}).get("checked") or {}
    if o:
        out.append("<h3>Read orientation</h3>")
        out.append(T(["Check", "Pairs", "Share"], [
            ["Already in the expected orientation", fmt(o.get("already_correct", 0)), f'{o.get("pct_already_correct", 0)}%'],
            ["Reversed (R1 and R2 swapped)", fmt(o.get("flipped", 0)), f'{o.get("pct_flipped", 0)}%'],
            ["Ambiguous", fmt(o.get("ambiguous", 0)), f'{o.get("pct_ambiguous", 0)}%']],
            "Read-pair orientation check." + (" Reversed pairs were re-oriented before analysis."
                                              if run.orientation.get("applied") else "")))

    if run.lengths():
        out.append("<h3>ASVs and amplicon length</h3>")
        L = run.lengths()
        lens = sorted(L)
        mode = max(L, key=lambda k: L[k])
        cv = ceiling_verdict(run)
        out.append(f"<p>The final ASV table contained {run.n_asvs} ASVs carrying {fmt(run.total)} reads. "
                   f"ASV lengths ranged from {lens[0]} to {lens[-1]} bp, with most reads at {mode} bp. "
                   + (cv[1] if cv[1] else "") + "</p>")
        out.append(F(svg_lengths(run), "Distribution of ASV lengths, weighted by reads."
                     + (f" The dashed line is the merge ceiling (truncLen_F + truncLen_R &minus; 12 = "
                        f"{run.ceiling} bp)." if run.ceiling else ""), folder("asv_lengths.csv", "asv_table.csv")))
        if cv[0]:
            out.append('<div class="note"><p style="margin:0"><b>The truncation settings are cutting the '
                       "amplicon short.</b> Longer fragments could not be merged, which biases the profile "
                       "toward short amplicons. Treat the composition as provisional until the run is "
                       "repeated with a wider window.</p></div>")

    out.append(T(["Sample", "Reads", f"n_{'taxa' if run.emu else 'ASV'}"],
                 [[esc(s), fmt(sum(an.mat[q][s] for q in run.abund)), fmt(an.n_features(s))] for s in S],
                 f"Reads and number of {unit} per sample after processing."))
    out.append(F(svg_rare(an), f"Rarefaction curves: expected number of {unit} against sequencing depth "
                 "for each sample.", folder("r_tables/rarefaction.csv")))

    qc = (run.qc or {}).get("offtarget") or {}
    if qc and qc.get("reads"):
        out.append("<h3>Off-target reads</h3>")
        out.append(T(["Category", "Reads", "Share"], [
            ["Host organelle (mitochondria / chloroplast)", fmt(qc.get("organelle_reads", 0)), f'{qc.get("organelle_pct", 0)}%'],
            ["Unassigned at domain level", fmt(qc.get("unassigned_reads", 0)), f'{qc.get("unassigned_pct", 0)}%'],
            ["Total off-target", fmt(qc.get("reads", 0)), f'{qc.get("pct", 0)}%']],
            "Off-target reads. &ldquo;Universal&rdquo; 16S primers co-amplify host mitochondrial DNA."))

    # ── composition ──────────────────────────────────────────────────────────
    out.append('<h3 class="pbreak">Taxonomic composition of samples</h3>')
    ph_all = an.rank_pct("Phylum")
    mean_ph = sorted(((t, sum(v.values()) / len(S)) for t, v in ph_all.items() if t != "Unclassified"),
                     key=lambda kv: -kv[1])[:4]
    txt = ("<p>At the phylum level the community was dominated by "
           + ", ".join(f"{esc(t)} ({v:.1f}% on average)" for t, v in mean_ph) + ". ")
    for s in (S if len(S) <= 8 else []):          # beyond 8 samples the table below says it
        dg = an.dominant("Genus", s, 3)
        dp = an.dominant("Phylum", s, 1)
        if dg:
            txt += (f"{esc(s)} was {'dominated by' if dg[0][1] >= 30 else 'led by'} <i>{esc(dg[0][0])}</i> "
                    f"({dg[0][1]:.1f}%)"
                    + (", followed by " + " and ".join(f"<i>{esc(t)}</i> ({v:.1f}%)" for t, v in dg[1:])
                       if len(dg) > 1 else "")
                    + (f", within {esc(dp[0][0])}" if dp else "") + ". ")
    if len(S) > 1:
        tops = {an.dominant("Genus", s, 1)[0][0] for s in S if an.dominant("Genus", s, 1)}
        txt += ("The same genus was the most abundant in every sample. " if len(tops) == 1 else
                f"The most abundant genus differed between samples ({len(tops)} different leading "
                f"genera), indicating distinct community structures. ")
    out.append(txt + "</p>")
    gu = an.rank_pct("Genus").get("Unclassified", {})
    if gu:
        out.append(f"<p>Reads without a genus assignment made up "
                   f"{min(gu.values()):.1f}&ndash;{max(gu.values()):.1f}% of each sample; they are shown "
                   f"as &ldquo;Unclassified&rdquo; and are part of every percentage below.</p>")
    dom_rows = []
    for s in S:
        dom_rows.append([esc(s),
                         "; ".join(f"{esc(t)} {v:.1f}%" for t, v in an.dominant("Phylum", s, 3)),
                         "; ".join(f"<i>{esc(t)}</i> {v:.1f}%" for t, v in an.dominant("Genus", s, 3))])
    out.append(T(["Sample", "Top phyla", "Top genera"], dom_rows, "Dominant taxa per sample."))

    for rank in RANKS:
        if rank not in run.rank_cols:
            continue
        keep, _t = an.top(rank, 30)
        if not keep:
            continue
        out.append(f'<h3 class="pbreak">Top 30 {rank} (Relative Abundance)</h3>')
        out.append(F(svg_stacked(an, rank, 30),
                     f"Relative abundance of the {len(keep)} most abundant "
                     f"{PLURAL.get(rank, rank.lower())} per sample"
                     + ". Percentages are of each sample&rsquo;s community reads"
                     + (" with host organelle reads set aside" if run.host + run.euk > 0 else "") + ".",
                     folder(f"r_tables/taxonomy_{rank.lower()}.csv")))

    out.append('<h3 class="pbreak">Heatmap visualization of dominant genera</h3>')
    out.append("<p>The genus-level heatmap summarises the 30 dominant genera across samples. Genera "
               "(rows) are clustered by the similarity of their abundance profiles"
               + (" and samples (columns) by Bray&ndash;Curtis dissimilarity" if len(S) >= 3 else "")
               + "; colour is relative abundance on a log scale, and each cell gives the percentage. "
               "Blocks of dark cells confined to one sample mark genera characteristic of it.</p>")
    out.append(F(svg_heatmap(an, "Genus", 30), "Top 30 genera &mdash; relative abundance heatmap with "
                 "clustering dendrogram.", folder("r_tables/taxonomy_heatmap_genus.csv")))

    out.append("<h3>Archaeal community</h3>")
    if n_arch > 0:
        out.append(f"<p>Archaeal reads were {'rare' if n_arch / max(fin, 1) < 0.01 else 'present'} in this "
                   f"dataset: {fmt(n_arch)} reads in total, "
                   f"{100*n_arch/max(sum(an.comm_tot.values()),1):.3f}% of the community reads. "
                   + ("The primers are designed for Bacteria and amplify most Archaea poorly, so this is a "
                      "minimum rather than an estimate of archaeal load." if not run.emu else "") + "</p>")
    else:
        out.append("<p>No archaeal reads were detected in any sample.</p>")
    out.append(T(["Sample", "Archaeal reads", "% of sample", "Dominant taxon"],
                 [[esc(s), fmt(r), f"{p_:.3f}%", f"<i>{esc(t)}</i>" if t != "—" else t]
                  for s, r, p_, t in arch_rows], "Archaeal reads per sample."))

    # reagent-associated genera, top features, species
    tg = run.top_genera(15)
    if tg:
        out.append('<h3 class="pbreak">Reagent-associated genera</h3>')
        out.append(T(["Genus", "Reads", "% bacteria", "Note"],
                     [[f'<i>{esc(g["name"])}</i>', fmt(g["reads"]), f'{g["pct"]}%',
                       '<span class="dim">' + "; ".join((["reagent-associated"] if g["kit"] else [])
                       + ([f'{fmt(g["unresolved"])} unresolved below genus'] if g["unresolved"] else []))
                       + "</span>"] for g in tg],
                     "The 15 most abundant genera over the whole run, with genera commonly recovered "
                     "from extraction kits and PCR reagents flagged."))
        out.append('<div class="note"><p style="margin:0"><b>Reagent-associated genera are flagged, not '
                   "removed.</b> They match the published list of genera recovered from DNA extraction "
                   "kits and PCR reagents (Salter et al. 2014). Many are also genuine environmental "
                   "organisms; without a sequenced blank extraction control the two cannot be told "
                   "apart.</p></div>")

    out.append(f"<h3>Most abundant {'taxa' if run.emu else 'sequence variants'}</h3>")
    per = len(S) <= 6                      # a share column per sample while it still fits
    def _lin(q):
        r = run.tax.get(q, {})
        parts = []
        for rk in ("Phylum", "Class", "Order", "Family"):
            v = str(r.get(rk) or "").strip()
            if v and not BLANK_RE.match(v):
                parts.append(UNRESOLVED_RE.sub("", v))
        return " &rsaquo; ".join(esc(x) for x in parts) or "unassigned"
    rows = []
    for a in run.top_asvs(top_asv):
        q = a["seq"]
        name = an.label(q, "Species") if run.emu else "Unclassified"
        if name == "Unclassified":
            name = a["genus"] if a["genus"] != "\u2014" else ""
        nm = ('<span class="dim">host organelle</span>' if a["host"]
              else f'<i>{esc(name)}</i>' if name else f'<span class="dim">{esc(a["deepest"])}</span>')
        rows.append([f'<span class="mono">{a["id"]}</span>'] + ([] if run.emu else [f'{a["len"]}'])
                    + [fmt(a["reads"]), f'{a["pct"]:.2f}%']
                    + ([f'{100*an.mat[q][s_]/an.comm_tot[s_]:.2f}' if an.comm_tot[s_] else "0"
                        for s_ in S] if per else [])
                    + [f'<div style="text-align:left">{nm}</div>',
                       f'<div style="text-align:left;font-size:7pt" class="dim">{_lin(q)}</div>'])
    out.append(T((["Taxon"] if run.emu else ["ASV", "bp"]) + ["Reads", "All"]
                 + ([f"{esc(h)} %" for h in _col_heads(S)] if per else [])
                 + (["Species / genus"] if run.emu else ["Genus"]) + ["Phylum &rsaquo; Class &rsaquo; Order &rsaquo; Family"],
                 rows,
                 f"The {top_asv} most abundant of {run.n_asvs} {unit}: reads and share over the whole run"
                 + (", and the share of each sample" if per else "")
                 + (". Identifiers are positional in abundance order; sequences are in asvs.fasta."
                    if not run.emu else ".")))

    out.append("<h3>Species-level assignment</h3>")
    spc = an.rank_counts("Species")
    named = sorted((t for t in spc if t != "Unclassified"), key=lambda t: -sum(spc[t].values()))
    if named:
        lim = 50
        tot_all = sum(an.comm_tot.values()) or 1
        rows = []
        for t in named[:lim]:
            v = sum(spc[t].values())
            rows.append([f"<i>{esc(t)}</i>", fmt(v), f"{100*v/tot_all:.2f}%"]
                        + ([f'{100*spc[t][s_]/an.comm_tot[s_]:.2f}' if an.comm_tot[s_] else "0"
                            for s_ in S] if per else []))
        cov = 100 * sum(sum(spc[t].values()) for t in named) / tot_all
        out.append(T(["Species", "Reads", "All"] + ([f"{esc(h)} %" for h in _col_heads(S)] if per else []),
                     rows,
                     f"Species-level assignments: {len(named)} reference taxa below genus carrying "
                     f"{cov:.1f}% of the community reads"
                     + (f" ({sum(1 for t in named if ' sp.' not in t)} named species; the others are "
                        f"unnamed reference strains such as &ldquo;sp. &hellip;&rdquo;)"
                        if any(' sp.' in t for t in named) else "")
                     + (f"; the {lim} most abundant are listed" if len(named) > lim else "")
                     + (", with the share of each sample." if per else ".")))
    else:
        out.append("<p>No read reached a species-level assignment with this reference.</p>")
    out.append("<p>" + ("Emu resolves species from the whole ~1,500 bp gene, but species whose 16S genes "
                        "are near-identical (for example within <i>Stenotrophomonas</i> or the "
                        "<i>Bacillus cereus</i> group) cannot be separated by any 16S method; genus is "
                        "the robust level to report." if run.emu else
                        "Species are assigned only by exact match of the whole fragment. A low rate is "
                        "expected for a ~440 bp V3&ndash;V4 amplicon; full-length 16S is the appropriate "
                        "marker where species identification is required.") + "</p>")

    # ── alpha ────────────────────────────────────────────────────────────────
    alpha = an.alpha()
    out.append('<h3 class="pbreak">Alpha Diversity</h3>')
    obs = [(s, alpha[s]["Observed"]) for s in S]
    sha = [(s, alpha[s]["Shannon"]) for s in S]
    sim = [(s, alpha[s]["Simpson"]) for s in S]
    f2 = lambda v: f"{v:.2f}"
    f3 = lambda v: f"{v:.3f}"
    at = (f"<p>Observed {unit} ranged {_range_txt(obs)} across samples. Shannon diversity ranged "
          f"{_range_txt(sha, f2)} and Simpson diversity {_range_txt(sim, f3)}. ")
    if len(S) > 1:
        hi_s = max(sha, key=lambda kv: kv[1])[0]
        lo_s = min(sha, key=lambda kv: kv[1])[0]
        at += (f"{esc(hi_s)} therefore combined the richest and most even community, while {esc(lo_s)} "
               f"was the most dominated by a few {unit}. " if hi_s == max(obs, key=lambda kv: kv[1])[0]
               else f"{esc(hi_s)} had the most even community and {esc(lo_s)} the least even. ")
    plats = [(s, an.plateau(s)) for s in S]
    plats = [(s, v) for s, v in plats if v is not None]
    if plats:
        flat = [s for s, v in plats if v < 1]
        rising = [s for s, v in plats if v >= 5]
        at += ("Rarefaction curves " + ("reached a plateau in every sample (less than 1% of richness "
               "gained over the last 10% of reads), so sequencing depth was sufficient to capture the "
               "community. " if len(flat) == len(plats) else
               f"were still rising for {', '.join(esc(s) for s in rising)}, so deeper sequencing would "
               f"recover further rare {unit} there. " if rising else
               "were close to a plateau in every sample, so depth captured nearly all of the community. "))
    if run.emu:
        at += ("Emu reports abundances rather than counts of distinct reads and drops taxa below its "
               "minimum abundance, so singletons are absent and Chao1 and ACE stay close to the "
               "observed count. ")
    out.append(at + "</p>")
    out.append(T(["Sample", "Reads", f"Observed {unit}", "Chao1", "ACE", "Shannon", "Simpson"],
                 [[esc(s), fmt(sum(an.mat[q][s] for q in run.abund)), fmt(alpha[s]["Observed"]),
                   f'{alpha[s]["Chao1"]:.1f}', f'{alpha[s]["ACE"]:.1f}', f'{alpha[s]["Shannon"]:.3f}',
                   f'{alpha[s]["Simpson"]:.3f}'] for s in S],
                 "Alpha diversity indices per sample. Shannon uses the natural logarithm; Simpson is "
                 "1 &minus; &Sigma;p&sup2;, approaching 1 when no single organism dominates."))
    out.append(F(svg_alpha(an, alpha), "Alpha diversity indices per sample.",
                 folder("r_tables/alpha_diversity.csv", "alpha_diversity.csv")))

    # ── beta ─────────────────────────────────────────────────────────────────
    if len(S) >= 2:
        D = an.bray()
        pairs = [(S[i], S[j], D[i][j]) for i in range(len(S)) for j in range(i + 1, len(S))]
        lo = min(pairs, key=lambda x: x[2])
        hi = max(pairs, key=lambda x: x[2])
        out.append('<h3 class="pbreak">Beta diversity</h3>')
        bt = (f"<p>Beta diversity compares community composition between samples. Bray&ndash;Curtis "
              f"dissimilarity (0 = identical, 1 = no shared {unit}) "
              + (f"ranged from {lo[2]:.2f} to {hi[2]:.2f}. " if len(pairs) > 1 else
                 f"between {esc(lo[0])} and {esc(lo[1])} was {_bcf(lo[2])}. "))
        if len(pairs) > 1:
            bt += (f"The most similar pair was {esc(lo[0])} and {esc(lo[1])} ({_bcf(lo[2])}); the most "
                   f"different was {esc(hi[0])} and {esc(hi[1])} ({_bcf(hi[2])}). ")
        mean_d = sum(x[2] for x in pairs) / len(pairs)
        bt += ("Overall the samples have highly distinct compositions. " if mean_d > 0.7 else
               "The samples share a substantial part of their composition. " if mean_d < 0.4 else
               "The samples differ moderately in composition. ")
        pm = an.permanova()
        if pm:
            bt += (f"PERMANOVA on group gave R&sup2; = {pm[0]:.3f}, p = {pm[1]:.3g}"
                   + (" &mdash; a significant group effect. " if pm[1] < 0.05 else
                      " &mdash; no significant group effect. "))
        elif an.groups and len(set(an.groups.values())) >= len(S):
            bt += ("With one sample per group no permutation test is possible, so the ordination and "
                   "clustering below are exploratory. ")
        else:
            bt += ("No replicated groups were defined for this run, so the comparison is descriptive. ")
        out.append(bt + "</p>")
        if len(S) <= 8:
            feats = list(run.abund)
            rows = []
            for i in range(len(S)):
                for j in range(i + 1, len(S)):
                    a_, b_ = S[i], S[j]
                    A = {q for q in feats if an.mat[q][a_] >= 0.5}
                    B = {q for q in feats if an.mat[q][b_] >= 0.5}
                    sh = A & B
                    ra = sum(an.mat[q][a_] for q in sh) / (sum(an.mat[q][a_] for q in A) or 1)
                    rb_ = sum(an.mat[q][b_] for q in sh) / (sum(an.mat[q][b_] for q in B) or 1)
                    jac = 1 - len(sh) / (len(A | B) or 1)
                    rows.append([f"{esc(a_)} &ndash; {esc(b_)}", fmt(len(sh)), fmt(len(A - B)), fmt(len(B - A)),
                                 f"{100*ra:.1f}% / {100*rb_:.1f}%", f"{jac:.3f}", _bcf(D[i][j])])
            out.append(T(["Pair", f"Shared {unit}", "Only first", "Only second",
                          "Reads in shared", "Jaccard", "Bray&ndash;Curtis"], rows,
                         f"{unit.capitalize()} shared between each pair of samples. &ldquo;Reads in shared&rdquo; is "
                         f"the share of each sample&rsquo;s reads that fall in the shared {unit}; Jaccard counts "
                         f"presence/absence only, Bray&ndash;Curtis weights by abundance."))
        # A taxon that dominates one sample and turns up at a trace level in another is
        # the signature of barcode / index cross-talk, not of a shared organism.
        gp = an.rank_pct("Genus")
        xt = []
        for g, v in gp.items():
            if g == "Unclassified":
                continue
            for a_ in S:
                for b_ in S:
                    if a_ != b_ and v[a_] >= 10 and 0 < v[b_] <= 0.5:
                        xt.append((g, a_, v[a_], b_, v[b_]))
        if xt:
            xt.sort(key=lambda x: -x[2])
            out.append("<p><b>Trace levels of another sample&rsquo;s dominant taxa.</b> "
                       + "; ".join(f"<i>{esc(g)}</i> is {va:.1f}% of {esc(a_)} but "
                                   f"{(f'{vb:.2f}' if vb >= 0.01 else '&lt;0.01')}% of {esc(b_)}"
                                   for g, a_, va, b_, vb in xt[:5])
                       + ". Traces like these can be genuine, but they are also what barcode / index "
                       "cross-talk between samples multiplexed in one run looks like, so they should not "
                       "be taken as members of the second sample&rsquo;s community without independent "
                       "evidence.</p>")
        pc = svg_pcoa(an)
        if pc:
            out.append(F(pc, f"Principal coordinates analysis (PCoA) of Bray&ndash;Curtis dissimilarity. "
                         f"Axis 1 explains {an.pcoa_pct[0]:.1f}% and axis 2 {an.pcoa_pct[1]:.1f}% of the "
                         f"variation.", folder("r_tables/pca_scores.csv", "beta_stats.txt")))
        if len(S) >= 3:
          out.append(F(svg_upgma(an), "UPGMA clustering tree based on Bray&ndash;Curtis dissimilarity.",
                     folder("beta_stats.txt")))
        out.append(F(svg_bc_heat(an), "Bray&ndash;Curtis dissimilarity between samples; darker cells and "
                     "higher values are more dissimilar.", folder("r_tables/beta_heatmap.csv", "beta_stats.txt")))

    # ── group view ───────────────────────────────────────────────────────────
    grp = {}
    for s in S:
        if s in an.groups:
            grp.setdefault(an.groups[s], []).append(s)
    if len(grp) >= 2 and any(len(v) >= 2 for v in grp.values()):
        out.append('<h3 class="pbreak">Composition by group</h3>')
        out.append("<p>Mean relative abundance of the samples in each group "
                   + "; ".join(f"{esc(g)}: {', '.join(esc(x) for x in v)}" for g, v in grp.items())
                   + ".</p>")
        for rank in ("Phylum", "Genus"):
            pct = an.rank_pct(rank)
            gt = {t: {g: sum(pct[t][s] for s in v) / len(v) for g, v in grp.items()} for t in pct}
            named = sorted((t for t in gt if t != "Unclassified"), key=lambda t: -sum(gt[t].values()))
            tbl_ = {t: gt[t] for t in named[:30]}
            if named[30:]:
                tbl_["Other"] = {g: sum(gt[t][g] for t in named[30:]) for g in grp}
            if "Unclassified" in gt:
                tbl_["Unclassified"] = gt["Unclassified"]
            out.append(F(svg_stacked(an, rank, 30, samples=list(grp), table=tbl_),
                         f"Top 30 {rank.lower()} &mdash; mean relative abundance per group."))

    if pic:
        out += _picrust_result(pic, F, T)

    # ── notes / references ───────────────────────────────────────────────────
    out.append('<h2 class="pbreak">Notes on reading this report</h2>')
    notes = []
    if run.total and run.host + run.euk > 0:
        notes.append(f"{100*(run.host+run.euk)/run.total:.1f}% of analysed reads are host organelle or "
                     "other non-bacterial material and are excluded from the composition figures.")
    if run.bacteria:
        notes.append(f"About {100*run.kit/run.bacteria:.0f}% of the bacterial reads belong to genera on "
                     "the reagent-contaminant list. A blank extraction control sequenced alongside the "
                     "samples would turn that estimate into a measurement.")
    if run.unresolved:
        notes.append(f"{fmt(sum(run.unresolved.values()))} reads carry a genus label the reference marks "
                     "as unresolved (&ldquo;--other&rdquo;); they are counted under the parent genus here.")
    if run.genera.get(HAEMOPLASMA_LABEL):
        notes.append("Haemotropic mycoplasmas are reported as one group (Eperythrozoon and "
                     "Haemobartonella combined).")
    notes.append("Relative abundance from 16S amplicons reflects gene copies, not cells: taxa with many "
                 "rRNA operons are over-represented.")
    notes.append(f"Full tables are in the result folder: {'taxonomy.csv' if run.emu else 'taxonomy_table.csv'} "
                 f"and asv_table.csv (per-{'taxon' if run.emu else 'ASV'} counts), r_tables/ (per-rank "
                 f"tables, diversity) and read_tracking.csv.")
    if pic:
        notes.append("Predicted functions (PICRUSt2) are in PICRUSt2/: KO_metagenome_out and "
                     "EC_metagenome_out (gene families per sample), pathways_out (MetaCyc pathways); the "
                     "*_descrip.tsv.gz files carry descriptions.")
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")
    out.append("<h2>References</h2>")
    out.append('<ul class="refs">' + "".join(
        f"<li>{esc(a)} <i>{esc(t)}</i> {esc(j)}</li>"
        for a, t, j in ([REFERENCES[i] for i in (REF_EMU if run.emu else REF_DADA2)]
                        + (REF_PICRUST if pic else []))) + "</ul>")

    css = _css(app, company, f"{run.label} · {marker}") + _CSS_FULL.replace("@@TITLE@@", esc(title))
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            + "".join(out) + "</body></html>")


_CSS_FULL = """
@page { @top-right { content: "@@TITLE@@"; font-family: Carlito, "DejaVu Sans", sans-serif;
        font-size: 7.4pt; font-weight: 700; color: #4E6672; border-bottom: 2.5px solid #9CC57A;
        padding-bottom: 1mm; } }
@page :first { @top-right { content: none; border: none; } }
h2 { font-size: 15pt; margin-top: 7mm; }
h3 { font-size: 11pt; margin-top: 5mm; }
figure { margin: 3mm 0 4mm; }
figcaption { font-size: 8pt; color: #16242B; }
.folder { font-weight: 700; font-size: 8.4pt; margin-top: 1.5mm; color: #16242B; }
.tcap { font-size: 8.4pt; margin: 4mm 0 1mm; break-after: avoid; }
.cmd { font-family: "DejaVu Sans Mono", monospace; font-size: 7.6pt; background: #F1F5F6;
       border-left: 2.5px solid #1E6E8C; padding: 2mm 3mm; word-break: break-all; }
.refs li { font-size: 8.4pt; margin-bottom: 1.6mm; }
.project { font-family: "Bitstream Charter", "DejaVu Serif", serif; font-size: 13pt; font-weight: 700;
           color: #1E6E8C; margin: 0 0 1mm; }
"""


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
    stages = run.stages
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
