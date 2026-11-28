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

  4. (Removed in v2.44.0 at the BI team's request: genera on a published list of extraction-kit /
     PCR-reagent contaminants were flagged as 'reagent-associated'. A list is no measurement: a
     genus on it is often a real member of the sample. Contamination is judged from sequenced
     blanks only — the Controls section.)

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
# Runs with more samples than MANY get their sample-wise figures on landscape pages, at most
# CHUNK samples to a page, so that every bar and label stays readable (an 89-sample coral run
# was printed as one A4 strip of hair-thin bars).
MANY, CHUNK = 24, 45


def _skey(x) -> str:
    """A sample name folded for matching (case, punctuation, R's make.names dots)."""
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def samples_in_groups(path, group_col: str, groups) -> list[str]:
    """Samples of metadata.csv whose value in group_col is one of groups (all grouped samples when
    groups is empty), in the file's order."""
    p = Path(path) / "metadata.csv"
    if not p.exists() or not group_col:
        return []
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or group_col not in rows[0]:
        return []
    idk = list(rows[0].keys())[0]
    want = {str(g).strip() for g in (groups or [])}
    return [str(r[idk]).strip() for r in rows
            if (r.get(group_col) or "").strip() and (not want or (r.get(group_col) or "").strip() in want)]


def _gammaincc(a, x):
    """Regularised upper incomplete gamma Q(a, x) (Numerical Recipes gser / gcf)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        ap, s_, d = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1
            d *= x / ap
            s_ += d
            if abs(d) < abs(s_) * 1e-12:
                break
        return max(0.0, 1.0 - s_ * math.exp(-x + a * math.log(x) - gln))
    b, c, d = x + 1 - a, 1e300, 1.0 / (x + 1 - a)
    h = d
    for i in range(1, 500):
        an_ = -i * (i - a)
        b += 2
        d = an_ * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an_ / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return min(1.0, math.exp(-x + a * math.log(x) - gln) * h)


def _kruskal(groups_vals):
    """Kruskal–Wallis H (ties corrected), df and p (chi-square) over [[values of a group], ...];
    None when it cannot be computed."""
    gv = [list(v) for v in groups_vals if len(v)]
    if len(gv) < 2:
        return None
    allv = sorted((x, gi) for gi, vs in enumerate(gv) for x in vs)
    N = len(allv)
    if N - len(gv) < 1:
        return None
    ranks, i, ties = [0.0] * N, 0, 0.0
    while i < N:
        j = i
        while j + 1 < N and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for q in range(i, j + 1):
            ranks[q] = r
        t = j - i + 1
        ties += t ** 3 - t
        i = j + 1
    R = [0.0] * len(gv)
    for (x, gi), r in zip(allv, ranks):
        R[gi] += r
    H = 12.0 / (N * (N + 1)) * sum(R[g] ** 2 / len(gv[g]) for g in range(len(gv))) - 3 * (N + 1)
    C = 1 - ties / (N ** 3 - N) if N > 1 else 1
    if C <= 0:
        return None
    H /= C
    df = len(gv) - 1
    return H, df, _gammaincc(df / 2.0, H / 2.0)


def _group_chunks(samples, groups, k=CHUNK):
    """Samples ordered by group (first-seen order; ungrouped last) and cut into pages of at most
    k, a page break falling between groups where the groups allow it."""
    order, seen = [], []
    for s_ in samples:
        g = groups.get(s_)
        if g is not None and g not in seen:
            seen.append(g)
    blocks = [[s_ for s_ in samples if groups.get(s_) == g] for g in seen]
    rest = [s_ for s_ in samples if s_ not in groups]
    if rest:
        blocks.append(rest)
    parts, cur = [], []
    for b in blocks:
        if len(b) > k:                       # a group larger than a page: its own near-equal parts
            if cur:
                parts.append(cur)
                cur = []
            parts.extend(_chunks(b, k))
            continue
        if len(cur) + len(b) > k:
            parts.append(cur)
            cur = []
        cur = cur + b
    if cur:
        parts.append(cur)
    return parts


def _chunks(items, k=CHUNK):
    """items in near-equal consecutive parts of at most k (89 -> 45 + 44, not 45 + 44 + 0)."""
    items = list(items)
    parts = max(1, math.ceil(len(items) / k))
    size = math.ceil(len(items) / parts) if items else 1
    return [items[i:i + size] for i in range(0, len(items), size)]

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


_DNA_RE = re.compile(r"^[ACGTUNRYSWKMBDHV]{20,}$", re.I)


def _is_dna(s) -> bool:
    return bool(_DNA_RE.match(str(s)))


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

    def __init__(self, path: Path, kitome: set[str], label: str | None = None,
                 samples=None, group_col: str | None = None, groups=None):
        """samples: report only these samples (a report for some groups of the run); everything
        — totals, tables, diversity — is then computed from them alone. group_col: the grouping
        (metadata.csv column) to analyse, over the one chosen in Customise charts; groups: the
        groups of it that were picked (their order is kept)."""
        self.path = Path(path)
        if not self.path.is_dir():
            raise FileNotFoundError(f"not a directory: {self.path}")
        self.only = list(samples) if samples else None
        self.group_pick = group_col or ""
        self.groups_pick = list(groups or [])
        self.kitome = kitome
        self.label = label or self.path.name

        self.params = self._json("run_params.json")
        self.summary = self._json("summary.json")
        self.qc = self._json("qc_report.json")
        self.prov = self._json("provenance.json")
        # Which pipeline made this folder. ONT 16S (Emu): one row per TAXON, not
        # per sequence variant; VSEARCH: one row per OTU; the DADA2 family (16S,
        # ITS, COX1, PacBio) one row per ASV.
        self.kind = self._detect_kind()
        self.emu = self.kind == "emu"

        self._load_tables()
        self._classify()
        self.blast, self.blast_info = self._load_blast()
        self.denoise_cap = self._denoise_cap()
        try:                                         # NCBI check of unassigned OTUs (VSEARCH runs)
            import ncbi_unassigned as _NU
            self.ncbi_un, self.ncbi_un_info = _NU.load(self.path)
        except Exception:
            self.ncbi_un, self.ncbi_un_info = {}, None

    def _detect_kind(self) -> str:
        pl = str(self.summary.get("pipeline") or "").lower()
        mk = str(self.params.get("marker") or self.summary.get("marker") or "").upper()
        if pl == "emu" or (self.path / "emu_output").is_dir():
            return "emu"
        if pl == "shotgun":                     # shotgun_pipeline.py (beta): rows are species
            return "shotgun"
        if pl == "qiime2_vsearch" or self.params.get("sequencerType") == "qiime2_vsearch":
            return "vsearch"
        if pl == "its" or mk in ("ITS", "ITS1", "ITS2"):
            return "its"
        if pl == "cox1" or mk == "COX1":
            return "cox1"
        if pl.startswith("pacbio") or mk == "PACBIO":
            return "pacbio"
        return "dada2"

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
        if not tax_rows and self.kind in ("emu", "vsearch", "shotgun"):
            tax_rows = self._rows("taxonomy.csv")          # emu / vsearch / shotgun name for it
        if not tax_rows:
            raise ValueError(f"{self.path.name}: taxonomy_table.csv missing or empty")
        key = list(tax_rows[0].keys())[0]  # first column = the ASV sequence (or feature id)
        tax_rows = self._clean_tax(tax_rows, key)
        self.tax = {r[key]: r for r in tax_rows}
        # SILVA 144 inserts a Kingdom rank below Domain, so the column set differs
        # between references. Take whatever ranks are actually present.
        self.rank_cols = [r for r in RANKS if r in tax_rows[0]]

        asv_rows = self._rows("asv_table.csv")
        if not asv_rows:
            raise ValueError(f"{self.path.name}: asv_table.csv missing or empty")
        cols = list(asv_rows[0].keys())
        if cols and cols[0].strip().lower() in ("sample_id", "sampleid", "sample"):
            # COX1 and PacBio before v2.17: one row per SAMPLE, a column per ASV id.
            # Turn it the right way round so everything below reads one layout.
            sc = cols[0]
            samples = [str(r[sc]) for r in asv_rows]
            asv_rows = [dict({"__id": fid}, **{samples[i]: asv_rows[i].get(fid) for i in range(len(samples))})
                        for fid in cols[1:]]
            cols = ["__id"] + samples
            seq_col = "__id"
        elif self.kind in ("emu", "vsearch", "shotgun"):
            seq_col = cols[0]                              # tax_id / OTU id, first column
        else:
            seq_col = "sequence" if "sequence" in cols else cols[-1]
        self.samples = [c for c in cols if c != seq_col and c != ""]
        self.tracking = self._rows("read_tracking.csv")
        self.all_samples = list(self.samples)
        if self.only:
            want = {_skey(x) for x in self.only}
            self.samples = [x for x in self.samples if _skey(x) in want]
            if not self.samples:
                raise ValueError(f"{self.path.name}: none of the chosen samples is in this run")
            keep = {_skey(x) for x in self.samples}
            gone = {_skey(x) for x in self.all_samples} - keep
            self.tracking = [t for t in self.tracking if _skey(next(iter(t.values()), "")) not in gone]
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
        if self.only:                    # features none of the chosen samples holds are not in this report
            for q in [q for q, v in self.abund.items() if v <= 0]:
                self.abund.pop(q, None)
                self.mat.pop(q, None)
        self.total = sum(self.abund.values())

        self.alpha = self._rows("alpha_diversity.csv")
        self.rarefaction = self._rows("rarefaction.csv")
        if self.only:
            gone = {_skey(x) for x in self.all_samples} - {_skey(x) for x in self.samples}
            drop = lambda r: any(_skey(v) in gone for v in list(r.values())[:2])
            self.alpha = [r for r in self.alpha if not drop(r)]
            self.rarefaction = [r for r in self.rarefaction if not drop(r)]
        self.orientation = self._json("orientation_report.json")

    def _load_blast(self):
        """BLAST naming of a 12S / CytB run (blast_taxonomy.py): {sequence: row} and blast_info.json.
        Only when the run asked for it and it ran — a skipped BLAST leaves the DADA2 names, which
        the report then describes as such."""
        if str(self.params.get("tax_method") or "dada2") not in ("blast", "both"):
            return {}, {}
        b = self.path / "blast"
        if (b / "BLAST_SKIPPED.txt").exists() or not (b / "taxonomy_blast.csv").exists():
            return {}, {}
        rows = {}
        for r in self._rows("blast/taxonomy_blast.csv"):
            q = str(r.get("sequence") or "").strip()
            if q:
                rows[q] = r
        return rows, self._json("blast/blast_info.json")

    def _denoise_cap(self) -> dict:
        """Samples the DADA2 pipeline subsampled before denoising ('[denoise cap] S: a -> b reads'):
        {sample: (before, after)}. Their drop at the denoising step is this cap, not read loss."""
        p = self.path / "pipeline.log"
        out = {}
        if not p.exists():
            return out
        try:
            for m in re.finditer(r"\[denoise cap\] (.+?): (\d+) -> (\d+) reads",
                                 p.read_text(encoding="utf-8", errors="replace")):
                out[m.group(1).strip()] = (int(m.group(2)), int(m.group(3)))
        except OSError:
            pass
        return out

    @property
    def mito_gene(self) -> str:
        """The mitochondrial gene of a 12S-pipeline run: 'CytB' or '12S' ('' for other markers)."""
        mk = str(self.summary.get("marker") or self.params.get("marker") or "").upper()
        if "12S" not in mk:
            return ""
        g = str(self.params.get("mito_gene") or "").upper()
        db = str(self.params.get("dbPath") or "").lower()
        return "CytB" if g == "CYTB" or "cytb" in db else "12S"

    def _clean_tax(self, rows, key):
        """Names as the reference meant them, whatever the run's age.

        UNITE prefixes ("g__Aspergillus") and MIDORI2 taxids ("Musca_7370") are
        stripped; MIDORI2 placeholders ("order_Sciaenidae") are blanked. COX1 runs
        before v2.17 named MIDORI2's six levels Kingdom..Genus although they are
        Phylum..Species; such a table is shifted back here, recognised by binomial
        names sitting in its "Genus" column.
        """
        if not rows:
            return rows
        ranks = [c for c in rows[0] if c != key]
        midori = self.kind == "cox1" or "midori" in str(self.params.get("dbPath") or "").lower()
        # PR2 (18S eukaryotes): 'Supergroup' / 'Subdivision' columns or the file name. Its placeholders
        # for a rank PR2 leaves undefined ('Chlorophyta_X', '_XX') are blanks; species are 'Genus_epithet'.
        self.is_pr2 = ("Supergroup" in ranks or "Subdivision" in ranks
                       or "pr2" in str(self.params.get("dbPath") or self.summary.get("db_path") or "").lower())
        out = []
        for r in rows:
            r = dict(r)
            for c in ranks:
                v = r.get(c)
                if v is None:
                    continue
                v = re.sub(r"^[a-z]__", "", str(v).strip())
                if midori:
                    v = re.sub(r"_[0-9]+$", "", v)
                    if re.match(r"^(kingdom|phylum|class|order|family|genus|species)_", v):
                        v = ""
                if re.fullmatch(r"(?i)unidentified|.*_incertae_sedis", v):
                    v = ""
                if self.is_pr2:
                    if re.search(r"_X+$", v) or re.search(r"(?i)_sp\.?$", v):
                        v = ""
                    elif c == "Species":
                        v = v.replace("_", " ")
                r[c] = v
            out.append(r)
        if "Species" in ranks and "Genus" in ranks:
            # VSEARCH naming before v2.33.3 cut MIDORI2 species ("Musca domestica") at the space and
            # wrote the genus again as the species; a one-word species equal to its genus is no species
            for r in out:
                sp, g = str(r.get("Species") or "").strip(), str(r.get("Genus") or "").strip()
                if sp and g and sp == g:
                    r["Species"] = ""
                    self.species_is_genus = getattr(self, "species_is_genus", 0) + 1
        if self.kind == "cox1" and "Species" not in rows[0] and "Genus" in rows[0]:
            gv = [str(r.get("Genus") or "") for r in out if r.get("Genus")]
            if gv and sum(" " in g for g in gv) > 0.5 * len(gv):
                std = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
                new = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
                out = [{**{k: v for k, v in r.items() if k not in std},
                        **{n: r.get(o, "") for o, n in zip(std, new)}} for r in out]
                self.rank_shift_fixed = True
        return out

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
        # shotgun: fungi and protists are part of the community too
        self.euk_target = any(k in marker for k in ("18S", "12S", "ITS", "COX", "COI", "NEMA", "EUK",
                                                    "SHOTGUN")) or self.kind == "shotgun"
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
                v = v.replace("_", " ")
                name = v if (" " in v or g == "?") else f"{g} {v}"
                self.species[name] = (self.species.get(name, (0, None))[0] + self.abund.get(s, 0),
                                      len(s) if _is_dna(s) else None)
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
            pid = getattr(self, "pkg_ids", None) or {}
            out.append({
                "id": (pid.get(seq) or pid.get(str(seq).split(";")[0])
                       or (f"T{i}" if self.emu or self.kind == "shotgun" else f"OTU{i}" if self.kind == "vsearch"
                           else f"ASV{i}")),
                "len": len(seq) if _is_dna(seq) else None, "reads": reads,
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
            if _is_dna(s):      # feature ids (old COX1/PacBio, VSEARCH OTUs) have none
                c[len(s)] += v
        return c

    @property
    def db_name(self) -> str:
        # what the run actually used (provenance.json), else what was asked for
        for d in (self.prov.get("databases") or []):
            if d.get("role") == "taxonomy" and d.get("file"):
                f = d["file"].split("/")[0] if self.emu else d["file"]
                lab = d.get("label") or ""
                return f"{lab} ({f})" if lab and lab != f and f not in lab else (lab or f)
        p = self.params.get("dbPath") or ""
        if p:
            pp = Path(p)
            # an Emu database is a folder; its file is always species_taxid.fasta
            return pp.parent.name if self.emu and pp.suffix else pp.name
        return (self.summary.get("taxonomy_db") or self.params.get("taxDatabase")
                or "(pipeline default)")

    def config(self, marker):
        """Run configuration rows for the DADA2-family and VSEARCH pipelines."""
        p, sm = self.params, self.summary
        prim = [("Forward primer", p.get("primer_f") or sm.get("primer_f") or "—"),
                ("Reverse primer", p.get("primer_r") or sm.get("primer_r") or "—")]
        lc = sm.get("lca") or {}
        mb = (("Classifier", "VSEARCH best hit, rank kept above an identity cut-off")
              if sm.get("tax_method") == "vsearch" else
              ("Naming", f"BLASTn vs NCBI nt, LCA of the top {lc.get('hits', 10)} hits at >= {lc.get('identity', 90)}% "
                         f"identity over {lc.get('coverage', 100)}% of the ASV")
              if sm.get("tax_method") == "ncbi_lca" else
              ("Minimum bootstrap", sm.get("minBoot") or p.get("minBoot")))
        if self.kind == "its":
            return ([("Marker", marker), ("Reference database", self.db_name)] + prim
                    + [("Length truncation", "none (ITS length varies between taxa)"),
                       ("maxEE F / R", f'{p.get("maxEE_F")} / {p.get("maxEE_R")}')]
                    + ([("Minimum length · pool · pair overlap",
                         f'{sm.get("minLen")} bp · {sm.get("pool")} · {sm.get("merge_overlap")} bp, no mismatch')]
                       if sm.get("minLen") else []) + [mb])
        if self.kind == "cox1":
            return ([("Marker", marker), ("Reference database", self.db_name)] + prim
                    + [("truncLen F / R", f'{p.get("truncLen_cox1_f", "—")} / {p.get("truncLen_cox1_r", "—")}'),
                       ("maxEE F / R", f'{p.get("maxEE_F")} / {p.get("maxEE_R")}'),
                       ("ASV length kept", sm.get("length_range")
                        or f'{p.get("cox1_min_len", "—")}-{p.get("cox1_max_len", "—")} bp'),
                       ("Genetic code (stop-codon filter)", sm.get("codon_table") or p.get("codon_table")),
                       ("LULU curation", {"ran": "yes", "off": "no",
                                          "missing": "skipped (LULU not installed)",
                                          "no_blast": "skipped (BLAST+ not found)",
                                          "error": "skipped (LULU failed on this run)",
                                          "empty_matchlist": "yes (no ASV pair at >= 84%: nothing merged)"}
                        .get(self.lulu_state, self.lulu_state)),
                       mb])
        if self.kind == "pacbio":
            return ([("Marker", "16S full length (PacBio HiFi)"), ("Reference database", self.db_name)]
                    + prim
                    + [("Read length kept", sm.get("length_range")
                        or f'{p.get("pb_min_len", "—")}-{p.get("pb_max_len", "—")} bp'),
                       ("maxEE", p.get("pb_maxEE")), ("Pooling", sm.get("pool") or p.get("pool") or "—"), mb])
        if self.kind == "vsearch":
            return ([("Marker", marker), ("Reference database", self.db_name)] + prim
                    + [("Reads", "paired-end, merged with VSEARCH" if self.paired else "single-end"),
                       ("Read length kept", "50 bp and longer after merging; expected errors ≤ 1"
                        if self.paired else f'{p.get("ontMinLen", "—")}-{p.get("ontMaxLen", "—")} bp'),
                       ("UNOISE alpha / minsize", f'{p.get("otuUnoiseAlpha", 2.0)} / {p.get("otuUnoiseMinsize", 8)}')])
        rows = [("Marker", marker), ("Reference database", self.db_name)] + prim + [
            ("truncLen F / R", self.trunc),
            ("maxEE F / R", f'{p.get("maxEE_F")} / {p.get("maxEE_R")}'),
            ("Pooling", p.get("pool") or "—"), mb,
            ("Chimera method", p.get("chimeraMethod"))]
        if self.blast:
            bi = self.blast_info or {}
            th = bi.get("thresholds") or {}
            rows += [("Names by", "BLASTn against MIDORI2 + MitoFish (DADA2 classifier kept beside)"),
                     ("BLAST limited to", {"fish": "fish (NCBI txid 7898 + 7777)", "vertebrates":
                                           "vertebrates (NCBI txid 7742)", "all": "no limit"}.get(
                                               str(bi.get("limit") or p.get("blast_limit") or "fish"), "fish")),
                     ("% identity (species / genus / family / order / class)",
                      " / ".join(str(th.get(k, "—")) for k in ("Species", "Genus", "Family", "Order", "Class"))
                      if th else str(p.get("blast_id") or "98.5,95,90,85,80").replace(",", " / ")),
                     ("NCBI nt check", "on" if str(bi.get("ncbi") or "off") != "off" else "off")]
        return rows

    def _trunc_pair(self):
        """The truncLen the pipeline actually used: COX1 has its own pair (truncLen_cox1_f/r), and the
        ITS pipelines do not truncate — truncLen_F/R ride along in their run_params unused."""
        p = self.params
        if self.kind == "cox1":
            return p.get("truncLen_cox1_f"), p.get("truncLen_cox1_r")
        if self.kind == "its":
            return None, None
        return p.get("truncLen_F"), p.get("truncLen_R")

    @property
    def trunc(self) -> str:
        f, r = self._trunc_pair()
        return f"{f} / {r}" if f and r else "—"

    @property
    def ceiling(self):
        f, r = self._trunc_pair()
        try:
            return int(f) + int(r) - 12 if f and r else None
        except (TypeError, ValueError):
            return None

    @property
    def app_version(self) -> tuple:
        v = str((self.prov or {}).get("app_version") or "")
        try:
            return tuple(int(x) for x in v.split(".")[:3])
        except ValueError:
            return ()

    @property
    def cutadapt_two_rounds(self) -> bool:
        """cutadapt -n 2 (both primers of a read-through pair removed) from v2.32.1 on."""
        return bool(self.app_version) and self.app_version >= (2, 32, 1)

    @property
    def lulu_state(self) -> str:
        """'ran', 'missing' (asked for, package not installed: skipped) or 'off' — COX1 runs."""
        sm = self.summary
        asked = sm.get("lulu_requested", sm.get("lulu", self.params.get("run_lulu", True)))
        if not asked:
            return "off"
        st = sm.get("lulu_status")             # v2.45.0: what happened, not only what was asked for
        if st:
            return {"ran": "ran", "missing": "missing", "off": "off"}.get(st, st)
        if "lulu_requested" in sm and not sm.get("lulu"):
            return "missing"
        try:
            if "LULU package not installed" in (self.path / "pipeline.log").read_text(errors="replace"):
                return "missing"
        except OSError:
            pass
        return "ran"

    @property
    def n_asvs(self) -> int:
        return self.summary.get("n_asvs") or self.summary.get("n_taxa") or len(self.tax)

    # wording and read stages that differ between DADA2 and Emu runs
    @property
    def unit(self) -> str:
        return {"emu": "taxa", "vsearch": "OTUs", "shotgun": "species"}.get(self.kind, "sequence variants")

    @property
    def feat(self) -> str:
        """Plural short name of a table row: taxa / OTUs / ASVs."""
        return {"emu": "taxa", "vsearch": "OTUs", "shotgun": "species"}.get(self.kind, "ASVs")

    @property
    def feat1(self) -> str:
        return {"emu": "taxon", "vsearch": "OTU", "shotgun": "species"}.get(self.kind, "ASV")

    @property
    def paired(self) -> bool:
        t0 = self.tracking[0] if self.tracking else {}
        return any(k in t0 for k in ("denoisedR", "denoised_R", "merged"))

    @property
    def reads_word(self) -> str:
        return "read pairs" if self.paired else "reads"

    @property
    def stages(self):
        if self.emu:
            last = "Classified" if _skipped(self).get("chimera") else "Non-chimeric"
            return [("Input", "input"), ("QC (Q / length)", "qc_filtered"),
                    ("Primers trimmed", "primer_trimmed"), (last, "final")]
        allst = {
            # Illumina pairs: trimmed, merged, then expected-error filtered; nanopore reads:
            # chopper's quality filter first, then primer trimming
            "vsearch": ([("Input", "input"), ("Primers trimmed", "primer_trimmed"),
                         ("Merged", "merged"), ("QC filtered", "qc_filtered"),
                         ("Mapped to OTUs", "final")] if self.paired else
                        [("Input", "input"), ("QC filtered", "qc_filtered"),
                         ("Primers trimmed", "primer_trimmed"), ("Mapped to OTUs", "final")]),
            # the length filter runs before chimera removal in cox1_pipeline.R, although
            # read_tracking.csv lists its column after nonchim
            "cox1": [("Input", "input"), ("Filtered", "filtered"), ("Merged", "merged"),
                     ("Length filtered", "length_filtered"), ("Non-chimeric", "nonchim"),
                     ("No stop codon", "codon_filtered")]
                    + ([("LULU curated", "lulu_curated")] if self.lulu_state in ("ran", "empty_matchlist") else []),
            "pacbio": [("Primer-oriented", "input"), ("Filtered", "filtered"),
                       ("Denoised", "denoised"), ("Non-chimeric", "nonchim")],
        }.get(self.kind, [("Input", "input"), ("Filtered", "filtered"),
                          ("Merged", "merged"), ("Non-chimeric", "nonchim")])
        t0 = self.tracking[0] if self.tracking else {}
        have = [(l, k) for l, k in allst if k in t0]
        return have or allst

    @property
    def retained_key(self) -> str:
        if self.emu:
            return "final"
        st = self.stages
        return st[-1][1] if st else "nonchim"

    @property
    def bacterial(self) -> bool:
        """A bacterial / archaeal marker: reagent kitome, archaea and organelle notes apply."""
        mk = str(self.summary.get("marker") or self.params.get("marker") or "16S").upper()
        return self.kind in ("emu", "pacbio") or (
            self.kind in ("dada2", "vsearch") and not any(k in mk for k in ("18S", "12S", "NEMA", "EUK", "ITS", "COX")))

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
    """Host / no-genus / bacteria named to genus, as a share of each run's reads."""
    segs = [("host", "Host / non-bacterial", "#B4C5CC"),
            ("nogen", "Bacteria, no genus", FAINT),
            ("samp", "Bacteria, genus named", ACC)]
    W, x0, bw, barh, pitch = 520, 6, 452, 28, 56
    H = len(runs) * pitch + 22
    p = [f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px">']
    for i, run in enumerate(runs):
        y = 18 + i * pitch
        tot = run.total or 1
        vals = {"host": run.host + run.euk, "nogen": run.no_genus,
                "samp": max(run.bacteria - run.no_genus, 0)}
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
                     f'{esc(nm)}</text>')
            p.append(f'<text x="{ox+colw:.1f}" y="{y+6}" font-size="7.5" fill="{MUT}" '
                     f'text-anchor="end" font-family="DejaVu Sans Mono">{g["pct"]}%</text>')
            bw = colw * g["pct"] / mx
            p.append(f'<rect x="{ox:.1f}" y="{y+9}" width="{bw:.1f}" height="5.5" '
                     f'fill="{col}"/>')
    p.append(f'<text x="0" y="{H-4}" font-size="7" fill="{FAINT}">'
             f'bars share one axis across all columns</text>')
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

    # ── 4. genera ─────────────────────────────────────────────────────────────
    out.append('<h2 class="pbreak">4. Genus-level profile</h2>')
    out.append("<p>Most abundant genera as a share of each run's bacterial reads, host "
               "material excluded.</p>")
    out.append("<figure>" + svg_genera(runs) + "<figcaption>Bars share one axis.</figcaption></figure>")
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
              ["Bacteria named to genus"] + [fmt(max(r.bacteria - r.no_genus, 0)) for r in runs]]
    out.append(_tbl([""] + [esc(n) for n in names], srows))

    bullets = []
    if multi:
        bullets.append("The runs compared here start from the same sequencing data; differences "
                       "in the tables above follow from what was changed between them.")
    r = runs[-1]
    if r.total:
        bullets.append(f"{100*(r.host+r.euk)/r.total:.0f}% of analysed reads are host or other "
                       "non-bacterial material co-amplified by the primers.")
    if r.total and (r.host + r.euk) / r.total >= 0.05:
        bullets.append("A host-blocking primer would recover most of the reads currently spent on host DNA.")
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
        """sample -> group. metadata.csv first: it is the file Customise charts writes, so groups
        set after the run win over the ones given with it (run_params.json). Its 'group' / 'Group'
        column, else its first other column holding 2+ values (a grouping named, say, 'Type')."""
        g, desc = {}, {}
        self.group_col = ""
        rows = self.run._rows("metadata.csv")
        if rows:
            keys = list(rows[0].keys())
            idk = keys[0] if keys else ""
            skip = {idk.lower(), "sampleid", "sample", "description", "sample_id", "#sampleid"}
            try:                     # the grouping chosen in Customise charts, when it is a column here
                pick = str(json.loads((self.run.path / "edit_charts" / "settings.json").read_text())
                           .get("groupCol") or "")
            except Exception:
                pick = ""
            col = pick if pick in keys and pick != idk else ""
            forced = getattr(self.run, "group_pick", "")
            if forced and forced in keys and forced != idk:
                col = forced
            col = col or next((k for k in keys if k in ("group", "Group")), "")
            if not col:
                col = next((k for k in keys if k.lower() not in skip
                            and len({(r.get(k) or "").strip() for r in rows} - {""}) >= 2), "")
            for r in rows:
                sid = str(r.get("sampleId") or r.get("sample") or r.get(idk) or "").strip()
                if sid and col and (r.get(col) or "").strip():
                    g[sid] = str(r[col]).strip()
                if sid and r.get("description"):
                    desc[sid] = str(r["description"]).strip()
            self.group_col = col if g else ""
        if not g:
            for m in (self.run.params.get("metadata") or []):
                sid = str(m.get("sampleId") or m.get("sample") or "").strip()
                if sid:
                    if m.get("group"):
                        g[sid] = str(m["group"]).strip()
                    if m.get("description"):
                        desc[sid] = str(m["description"]).strip()
            self.group_col = "group" if g else ""
        self.desc = desc
        return {s: g[s] for s in self.samples if s in g}

    def group_list(self):
        """[(group, [samples])]: the groups picked for the report in their order, else in
        first-seen order of the run's samples."""
        out = OrderedDict((g, []) for g in getattr(self.run, "groups_pick", []) or [])
        for s_ in self.samples:
            if s_ in self.groups:
                out.setdefault(self.groups[s_], []).append(s_)
        return [(g, v) for g, v in out.items() if v]

    def permanova_groups(self, perms=999, seed=1):
        """PERMANOVA (Anderson 2001) of Bray–Curtis on the current groups, computed here so that it
        follows metadata.csv as it is now (the run's beta_stats.txt used the groups of the run).
        {"R2", "F", "p", "n", "k", "perms"} or None without two groups and a replicate."""
        if getattr(self, "_perm", "x") != "x":
            return self._perm
        import random as _r
        S = self.samples
        idx = [i for i, s_ in enumerate(S) if s_ in self.groups]
        labs = [self.groups[S[i]] for i in idx]
        n, levels = len(idx), sorted(set(labs))
        k = len(levels)
        self._perm = None
        if k < 2 or n - k < 1 or max(labs.count(l) for l in levels) < 2:
            return None
        D = self.bray()
        code = {l: i for i, l in enumerate(levels)}
        lab = [code[l] for l in labs]
        pairs = [(i, j, D[idx[i]][idx[j]] ** 2) for i in range(n) for j in range(i + 1, n)]
        sst = sum(v for _i, _j, v in pairs) / n
        if sst <= 0:
            return None

        def ssw(L):
            acc = [0.0] * k
            cnt = [0] * k
            for c in L:
                cnt[c] += 1
            for i, j, v in pairs:
                if L[i] == L[j]:
                    acc[L[i]] += v
            return sum(acc[c] / cnt[c] for c in range(k) if cnt[c])

        def fstat(w):
            b = sst - w
            return (b / (k - 1)) / (w / (n - k)) if w > 0 else float("inf")
        w0 = ssw(lab)
        f0 = fstat(w0)
        rnd, L, hits = _r.Random(seed), list(lab), 1
        for _ in range(perms):
            rnd.shuffle(L)
            if fstat(ssw(L)) >= f0 - 1e-12:
                hits += 1
        self._perm = {"R2": (sst - w0) / sst, "F": f0, "p": hits / (perms + 1), "n": n, "k": k, "perms": perms}
        return self._perm

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


def _xlabels(p, labels, xs, y, maxlen=22, size=8.4, vertical=False):
    rot = len(labels) > 3 or max((len(str(l)) for l in labels), default=0) > 10
    for l, x in zip(labels, xs):
        t = esc(_short(l, maxlen))
        if vertical:                      # reads bottom-to-top, hanging below the axis
            xx = x + size * 0.34
            p.append(f'<text x="{xx:.1f}" y="{y:.1f}" font-size="{size}" fill="{INK}" '
                     f'text-anchor="end" transform="rotate(-90 {xx:.1f} {y:.1f})">{t}</text>')
        elif rot:
            p.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{INK}" '
                     f'text-anchor="end" transform="rotate(-40 {x:.1f} {y:.1f})">{t}</text>')
        else:
            p.append(f'<text x="{x:.1f}" y="{y+4:.1f}" font-size="{size}" fill="{INK}" '
                     f'text-anchor="middle">{t}</text>')


def svg_stacked(an, rank, n=30, samples=None, table=None, italic=None, wide=False, groups=None) -> str:
    """Top-n taxa per sample as 100% stacked bars, with a legend down the right. wide: for a
    landscape page (one part of a many-sample run; colours stay those of the whole run)."""
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
    if wide:
        W, x0, y0, ph = 1000, 46, (40 if groups else 24), 440
        leg_x = W - 225
        pw = leg_x - x0 - 18
        bot = 14 + min(longest, 22) * 4.8
    else:
        W, x0, y0, ph = 520, 40, (40 if groups else 22), 430
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
    bw = min(slot * (0.74 if wide else 0.66), 84)
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
    _xlabels(p, samples, xs, y0 + ph + 8, size=8, vertical=wide)
    if groups:                                  # a labelled bracket over each group's bars
        i = 0
        while i < len(samples):
            g = groups.get(samples[i], "")
            j = i
            while j + 1 < len(samples) and groups.get(samples[j + 1], "") == g:
                j += 1
            if g:
                xa, xb = xs[i] - slot * 0.42, xs[j] + slot * 0.42
                p.append(f'<line x1="{xa:.1f}" y1="{y0-6}" x2="{xb:.1f}" y2="{y0-6}" stroke="{INK}" stroke-width="1"/>')
                p.append(f'<text x="{(xa+xb)/2:.1f}" y="{y0-10}" font-size="8.6" fill="{INK}" font-weight="700" '
                         f'text-anchor="middle">{esc(_short(g, max(4, int((xb-xa)/5))))}</text>')
            i = j + 1
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


def svg_heatmap(an, rank="Genus", n=30, samples=None, wide=False) -> str:
    """Top-n taxa × samples, clustered: rows by profile (UPGMA, Euclidean on log %),
    columns by Bray-Curtis when there are three or more samples. samples: one part of a
    many-sample run, in the order given (rows and colour scale stay those of the whole run)."""
    keep, table = an.top(rank, n)
    if not keep:
        return ""
    part = samples is not None
    samples = list(samples) if part else list(an.samples)
    # columns
    col_tree = None
    if len(samples) >= 3 and not part:
        vec = [[table[t][s] for t in keep] for s in samples]
        D = [[_bray(a, b) for b in vec] for a in vec]
        col_tree = _upgma(D)
        samples = [samples[i] for i in _leaves(col_tree)]
    # rows
    lv = [[math.log10(table[t][s] + 0.01) for s in an.samples] for t in keep]
    row_tree = _upgma([[_euclid(a, b) for b in lv] for a in lv]) if len(keep) > 1 else None
    rows = [keep[i] for i in _leaves(row_tree)] if row_tree else keep
    vmax = max((table[t][s] for t in keep for s in an.samples), default=1) or 1
    W = 1000 if wide else 520
    dend_w = 70
    lab_w = 175 if wide else 160
    top_h = 52 if col_tree else 8
    cw = min(92, (W - dend_w - lab_w - 10) / max(len(samples), 1))
    ch = 12.6 if wide else 15.0
    x0 = dend_w + 4
    y0 = 20 + top_h
    longest = max((len(str(s)) for s in samples), default=4)
    H = int(y0 + ch * len(rows) + 22 + min(longest, 22) * (4.8 if wide else 4.0) + 38)
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
    _xlabels(p, samples, [x0 + (j + 0.5) * cw for j in range(len(samples))], yb + 9,
             size=8 if wide else 8.4, vertical=wide)
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
    if run.kind not in ("dada2", "its"):
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
    many = len(data) > MANY
    for si, (name, vs) in enumerate(data):
        c = FAINT if many else (SERIES[si] if len(data) <= len(SERIES) else _pal(si))
        pts = " ".join(f"{x:.1f},{y0+ph-ph*v/mx:.1f}" for x, v in zip(xs, vs))
        p.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="{.6 if many else 1}" '
                 f'opacity="{.45 if many else .75}"><title>{esc(name)}</title></polyline>')
        for x, v in zip(xs, vs):
            p.append(f'<circle cx="{x:.1f}" cy="{y0+ph-ph*v/mx:.1f}" r="{1.3 if many else 2.1}" fill="{c}"/>')
    for (l, _k), x in zip(steps, xs):
        p.append(f'<text x="{x:.1f}" y="{y0+ph+12}" font-size="8.2" fill="{INK}" text-anchor="middle">{esc(l)}</text>')
    p.append(f'<text x="{x0+pw/2:.1f}" y="{H-24}" font-size="8.4" fill="{MUT}" text-anchor="middle">Pipeline step</text>')
    p.append(f'<text x="9" y="{y0+ph/2:.1f}" font-size="8.4" fill="{MUT}" text-anchor="middle" '
             f'transform="rotate(-90 9 {y0+ph/2:.1f})">Read count</text>')
    lx = x0 + pw + 14
    if many:                                    # 89 names do not fit a legend: the table lists them
        p.append(f'<text x="{lx}" y="{y0+8}" font-size="8" fill="{MUT}">{len(data)} samples, one grey</text>')
        p.append(f'<text x="{lx}" y="{y0+20}" font-size="8" fill="{MUT}">line each; boxes: spread</text>')
        p.append(f'<text x="{lx}" y="{y0+32}" font-size="8" fill="{MUT}">per step (see the table)</text>')
    for si, (name, _vs) in enumerate([] if many else data[:18]):
        c = SERIES[si] if len(data) <= len(SERIES) else _pal(si)
        y = y0 + 4 + si * 12.5
        p.append(f'<circle cx="{lx+3}" cy="{y:.1f}" r="2.6" fill="{c}"/>')
        p.append(f'<text x="{lx+10}" y="{y+2.6:.1f}" font-size="8" fill="{INK}">{esc(_short(name, 24))}</text>')
    p.append("</svg>")
    return "".join(p)


def svg_alpha(an, alpha, idx=None, wide=False, samples=None) -> str:
    """Four panels — Observed, Chao1, Shannon, Simpson — one bar per sample. wide: the given
    indices one above the other across a landscape page (many samples)."""
    samples = samples or an.samples
    if wide:
        return _svg_alpha_wide(samples, alpha, idx or ["Observed", "Chao1"])
    idx = ["Observed", "Chao1", "Shannon", "Simpson"]
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


def _svg_alpha_wide(samples, alpha, idx) -> str:
    longest = max((len(str(s)) for s in samples), default=4)
    W, ox, pw_, ph_ = 1000, 58, 922, 168
    bot = 12 + min(longest, 22) * 4.2
    cellh = ph_ + bot + 30
    H = int(cellh * len(idx) + 8)
    p = _svg_open(W, H)
    slot = pw_ / max(len(samples), 1)
    bw = min(slot * 0.62, 40)
    for k, name in enumerate(idx):
        oy = 24 + k * cellh
        vals = [alpha[s][name] for s in samples]
        mx = 1.0 if name == "Simpson" else (max(vals, default=1) or 1)
        p.append(f'<text x="0" y="{oy-11}" font-size="10" fill="{INK}" font-weight="700">{name}</text>')
        for g in (0, .25, .5, .75, 1):
            y = oy + ph_ - ph_ * g
            v = mx * g
            p.append(f'<line x1="{ox}" y1="{y:.1f}" x2="{ox+pw_}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
            lab = f"{v:.2f}" if name in ("Shannon", "Simpson") else fmt(v)
            p.append(f'<text x="{ox-4}" y="{y+2.5:.1f}" font-size="7.8" fill="{FAINT}" text-anchor="end" {MONO}>{lab}</text>')
        xs = []
        for i, (s, v) in enumerate(zip(samples, vals)):
            xc = ox + slot * (i + 0.5)
            xs.append(xc)
            h = ph_ * v / mx
            p.append(f'<rect x="{xc-bw/2:.1f}" y="{oy+ph_-h:.1f}" width="{bw:.1f}" height="{h:.1f}" fill="{ACC}">'
                     f'<title>{esc(s)} {v:.3g}</title></rect>')
        _xlabels(p, samples, xs, oy + ph_ + 6, size=min(7.2, slot * 0.8), vertical=True)
    p.append("</svg>")
    return "".join(p)


def svg_alpha_groups(glist, alpha, idx=("Observed", "Chao1", "Shannon", "Simpson")) -> str:
    """One panel per index: a box (quartiles, median, whiskers to the extremes) and the samples as
    dots, per group."""
    names = [g for g, _ in glist]
    longest = max((len(str(g)) for g in names), default=4)
    rot = len(names) > 6 or longest > 9
    pw_, ph_ = 205, 130
    bot = (16 + min(longest, 18) * 3.8) if rot else 20
    cellh = ph_ + bot + 26
    W, H = 520, int(cellh * ((len(idx) + 1) // 2) + 14)
    p = _svg_open(W, H)
    for k, name in enumerate(idx):
        ox = 40 + (k % 2) * 262
        oy = 26 + (k // 2) * cellh
        vals = [[alpha[s_][name] for s_ in mem] for _g, mem in glist]
        mx = 1.0 if name == "Simpson" else (max((v for vs in vals for v in vs), default=1) or 1)
        p.append(f'<text x="{ox-36}" y="{oy-13}" font-size="9.4" fill="{INK}" font-weight="700">{name}</text>')
        for gg in (0, .5, 1):
            y = oy + ph_ - ph_ * gg
            lab = f"{mx*gg:.2f}" if name in ("Shannon", "Simpson") else fmt(mx * gg)
            p.append(f'<line x1="{ox}" y1="{y:.1f}" x2="{ox+pw_}" y2="{y:.1f}" stroke="{RULE}" stroke-width=".6"/>')
            p.append(f'<text x="{ox-4}" y="{y+2.5:.1f}" font-size="7.6" fill="{FAINT}" text-anchor="end" {MONO}>{lab}</text>')
        slot = pw_ / max(len(glist), 1)
        bw = min(slot * 0.5, 40)
        xs = []
        Y = lambda v: oy + ph_ - ph_ * v / mx
        for gi, vs in enumerate(vals):
            xc = ox + slot * (gi + 0.5)
            xs.append(xc)
            c = _pal(gi)
            if not vs:
                continue
            q = sorted(vs)
            qq = lambda f: q[0] if len(q) == 1 else (q[int(math.floor(f * (len(q) - 1)))]
                                                     + (q[min(len(q) - 1, int(math.floor(f * (len(q) - 1))) + 1)]
                                                        - q[int(math.floor(f * (len(q) - 1)))])
                                                     * (f * (len(q) - 1) - math.floor(f * (len(q) - 1))))
            p.append(f'<line x1="{xc:.1f}" y1="{Y(q[0]):.1f}" x2="{xc:.1f}" y2="{Y(q[-1]):.1f}" stroke="{MUT}" stroke-width=".8"/>')
            p.append(f'<rect x="{xc-bw/2:.1f}" y="{Y(qq(.75)):.1f}" width="{bw:.1f}" '
                     f'height="{max(Y(qq(.25))-Y(qq(.75)), .8):.1f}" fill="{c}" fill-opacity=".25" stroke="{c}" stroke-width="1"/>')
            p.append(f'<line x1="{xc-bw/2:.1f}" y1="{Y(qq(.5)):.1f}" x2="{xc+bw/2:.1f}" y2="{Y(qq(.5)):.1f}" stroke="{INK}" stroke-width="1.3"/>')
            for m_, v in enumerate(vs):        # dots spread across the box, deterministically
                jx = (((m_ * 37) % 11) / 10 - 0.5) * bw * 0.8
                p.append(f'<circle cx="{xc+jx:.1f}" cy="{Y(v):.1f}" r="1.7" fill="{c}" fill-opacity=".85"/>')
        _xlabels(p, names, xs, oy + ph_ + 8, maxlen=18, size=7.4) if rot else [
            p.append(f'<text x="{x:.1f}" y="{oy+ph_+12:.1f}" font-size="7.8" fill="{INK}" text-anchor="middle">'
                     f'{esc(_short(g, 12))}</text>') for g, x in zip(names, xs)]
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
             f'transform="rotate(-90 9 {y0+ph/2:.1f})">Observed {an.run.feat}</text>')
    lx = x0 + pw + 14
    if len(samples) > MANY:
        p.append(f'<text x="{lx}" y="{y0+8}" font-size="8" fill="{MUT}">{len(samples)} samples, one curve</text>')
        p.append(f'<text x="{lx}" y="{y0+20}" font-size="8" fill="{MUT}">each; depths and counts</text>')
        p.append(f'<text x="{lx}" y="{y0+32}" font-size="8" fill="{MUT}">per sample are in the tables</text>')
    for i, s in enumerate([] if len(samples) > MANY else samples[:20]):
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


def svg_pcoa(an, wide=False) -> str:
    """wide: a landscape page for many samples — labels only where they do not collide,
    no sample legend."""
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
    if wide:
        W, H, x0, y0, pw, ph = 1000, 540, 56, 14, (760 if cmap else 920), 480
    else:
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
    if wide:
        for s, (cx, cy), c in zip(an.samples, coords, cols):
            p.append(f'<circle cx="{X(cx):.1f}" cy="{Y(cy):.1f}" r="4" fill="{c}" stroke="#fff" stroke-width=".8">'
                     f'<title>{esc(s)}</title></circle>')
        for s, (cx, cy) in zip(an.samples, coords):
            t = _short(s, 24)
            w_ = len(t) * 3.9
            for dx, dy, anc in ((6, -5, "start"), (6, 9, "start"), (-6, -5, "end"), (-6, 9, "end")):
                lx_, ly_ = X(cx) + dx, Y(cy) + dy
                x1_ = lx_ - w_ if anc == "end" else lx_
                if x1_ < x0 or x1_ + w_ > x0 + pw or ly_ < y0 + 8 or ly_ > y0 + ph - 2:
                    continue
                if not any(abs(ly_ - py) < 8.5 and x1_ < px + q_ and px < x1_ + w_ for px, py, q_ in placed):
                    placed.append((x1_, ly_, w_))
                    p.append(f'<text x="{lx_:.1f}" y="{ly_:.1f}" font-size="7.2" fill="{INK}" '
                             f'text-anchor="{anc}">{esc(t)}</text>')
                    break
    for s, (cx, cy), c in ([] if wide else zip(an.samples, coords, cols)):
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
    items = list(cmap.items()) if cmap else ([] if wide else list(zip(an.samples, cols)))
    if items:
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
    rowh = min(22.0, max(6.0, 660.0 / len(order)))      # many samples: the tree on one page
    fs = min(8.6, rowh * 0.86)
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
        p.append(f'<text x="{x0+pw+7:.1f}" y="{pos[k]+fs*0.33:.1f}" font-size="{fs:.1f}" fill="{INK}">{esc(an.samples[k])}</text>')
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
    # 17
    ("Abarenkov, K., et al. (2024).", "The UNITE database for molecular identification and taxonomic "
     "communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered.",
     "Nucleic Acids Research, 52(D1), D791–D797."),
    ("Leray, M., Knowlton, N., & Machida, R. J. (2022).", "MIDORI2: A collection of quality controlled, "
     "preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic "
     "mitochondrial sequences.", "Environmental DNA, 4(4), 894–907."),
    ("Frøslev, T. G., et al. (2017).", "Algorithm for post-clustering curation of DNA amplicon data yields "
     "reliable biodiversity estimates.", "Nature Communications, 8, 1188."),
    ("Callahan, B. J., et al. (2019).", "High-throughput amplicon sequencing of the full-length 16S rRNA gene "
     "with single-nucleotide resolution.", "Nucleic Acids Research, 47(18), e103."),
    ("Edgar, R. C. (2016).", "UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing.",
     "bioRxiv, 081257."),
    # 22
    ("Guillou, L., et al. (2013).", "The Protist Ribosomal Reference database (PR2): a catalog of unicellular "
     "eukaryote Small Sub-Unit rRNA sequences with curated taxonomy.", "Nucleic Acids Research, 41, D597–D604."),
    ("Chuvochina, M., et al. (2026).", "SILVA in 2026: a global core biodata resource for rRNA within the "
     "DSMZ digital diversity.", "Nucleic Acids Research, 54(D1), D334–D341."),
    ("Parks, D. H., et al. (2022).", "GTDB: an ongoing census of bacterial and archaeal diversity through a "
     "phylogenetically consistent, rank normalized and complete genome-based taxonomy.",
     "Nucleic Acids Research, 50(D1), D785–D794."),
    # 25
    ("Nickols, W. A., et al. (2024).", "MaAsLin 3: Refining and extending generalized multivariable linear "
     "models for meta-omic association discovery.", "bioRxiv, 2024.12.13.628459."),
    ("Camacho, C., et al. (2009).", "BLAST+: architecture and applications.", "BMC Bioinformatics, 10, 421."),
    ("Zhu, T., Sato, Y., Sado, T., Miya, M., & Iwasaki, W. (2023).", "MitoFish, MitoAnnotator, and MiFish "
     "Pipeline: updates in 10 years.", "Molecular Biology and Evolution, 40(3), msad035."),
    ("Vu, D., et al. (2019).", "Large-scale generation and analysis of filamentous fungal DNA barcodes boosts "
     "coverage for kingdom fungi and reveals thresholds for fungal species and higher taxon delimitation.",
     "Studies in Mycology, 92, 135–154."),
    # 29 — ITS named by BLASTn against NCBI nt and the LCA of the top hits (ncbi_lca.py)
    # cited as the request form gives it (journal, volume, pages, DOI); no title is made up here
    ("DiBattista et al. (2025).", "", "Coral Reefs, 44, 163–178. https://doi.org/10.1007/s00338-024-02594-6"),
    ("Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990).", "Basic local alignment "
     "search tool.", "Journal of Molecular Biology, 215(3), 403–410."),
]


def ncbi_lca(run) -> bool:
    """ITS named by BLASTn against NCBI nt + LCA (its_pipeline.R --tax_source ncbi_lca) instead of UNITE."""
    return str((getattr(run, "summary", None) or {}).get("tax_method") or "") == "ncbi_lca"
REF_DADA2 = [0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16]
REF_EMU = [5, 0, 6, 7, 4, 8, 9, 10, 11, 12, 13, 14, 15]


def _ref_indices(run) -> list[int]:
    """The reference list for this run: its tools, its database, its statistics."""
    if run.emu:
        return REF_EMU
    if ncbi_lca(run):                    # no classifier, no reference file: BLAST at NCBI + LCA
        return [0, 1, 30, 29, 8, 9, 10, 11, 12, 13, 14, 15]
    db = " ".join([run.db_name] + [d.get("label", "") for d in (run.prov.get("databases") or [])]).lower()
    dbref = ([23] if "silva 14" in db or "v144" in db else [3] if "silva" in db else []) \
        + ([17] if "unite" in db or "sh_general" in db or run.kind == "its" else []) \
        + ([18] if "midori" in db or run.kind == "cox1" else []) \
        + ([22] if "pr2" in db else []) + ([24] if "gtdb" in db else [])
    if run.kind == "vsearch":
        tools = [0, 6, 21] + ([5] if not run.paired else [])
        if str(run.params.get("marker") or run.summary.get("marker") or "").upper().startswith("ITS"):
            dbref = dbref + [28]
    else:
        tools = [0, 1] + ([20] if run.kind == "pacbio" else []) + [2] + ([19] if run.kind == "cox1" else [])
    if getattr(run, "blast", None):
        tools = [t for t in tools if t != 2] + [26]       # BLASTn names, not the naive Bayesian classifier
        dbref = [d for d in dbref if d != 18] + [18, 27]   # MIDORI2 + MitoFish
    out = tools + dbref + [8, 9, 10, 11, 12, 13, 14, 15]
    if (run.path / "r_tables" / "maaslin3_results.csv").exists():
        out.append(25)
    return list(dict.fromkeys(out))      # 16 (Salter 2014) went with the reagent table (v2.42.0)


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


def _provenance_rows(run):
    """[[name, version, md5]] from provenance.json; [] for runs made before v2.17."""
    pv = run.prov or {}
    if not pv:
        return []
    rows = []
    if pv.get("app_version"):
        rows.append(["NextGen-Amplicon", esc(pv["app_version"]), ""])
    names = {"R": "R", "dada2": "DADA2", "DECIPHER": "DECIPHER", "phyloseq": "phyloseq",
             "vegan": "vegan", "Biostrings": "Biostrings", "ShortRead": "ShortRead", "lulu": "LULU",
             "cutadapt": "cutadapt", "vsearch": "VSEARCH", "emu": "Emu", "minimap2": "minimap2",
             "chopper": "chopper"}
    for k, v in (pv.get("tools") or {}).items():
        rows.append([esc(names.get(k, k)), esc(v), ""])
    for d in pv.get("databases") or []:
        md5 = d.get("md5") or ""
        ok = d.get("md5_ok")
        mark = " &#10003;" if ok else (" &#10007; differs from the published file" if ok is False else "")
        role = {"taxonomy": "Reference", "species": "Species reference"}.get(d.get("role"), "Reference")
        rows.append([f'{role}: {esc(d.get("label") or d.get("file") or "")}',
                     f'<span class="dim" style="font-size:7pt">{esc(d.get("file") or "")}</span>',
                     f'<span class="mono" style="font-size:7pt">{esc(md5)}</span>{mark}'])
    return rows


def _tophit_text(run, marker):
    cut = run.summary.get("tax_cutoffs") or {}
    order = ["Species", "Genus", "Family", "Order", "Class", "Phylum"]
    cuts = ", ".join(f"{r.lower()} {cut[r]:g}%" for r in order if isinstance(cut.get(r), (int, float)))
    src = ("Ratnasingham &amp; Hebert 2013; common COI metabarcoding practice" if marker == "COI"
           else "Vu et al. 2019")
    return (f"<p>Each ASV took the lineage of its best match in <b>{esc(run.db_name)}</b> "
            f"(VSEARCH <span class=\"mono\">--usearch_global</span>, both strands), and a rank was kept "
            f"only where the identity reached its threshold" + (f" ({cuts}; {src})" if cuts else "")
            + ". This reference has too many lineages for the DADA2 classifier to hold in memory "
              "(it keeps one k-mer profile per lineage), so best-hit naming was used instead. ")


def _methods_other(run, an):
    """Method text for ITS, COX1, PacBio and VSEARCH runs."""
    p, sm, out = run.params, run.summary, []
    pf, pr = p.get("primer_f") or sm.get("primer_f") or "", p.get("primer_r") or sm.get("primer_r") or ""
    mb = sm.get("minBoot") or p.get("minBoot") or 50
    prim_ul = ("<ul>" + (f"<li>Forward primer: <b class=\"mono\">{esc(pf)}</b></li>" if pf else "")
               + (f"<li>Reverse primer: <b class=\"mono\">{esc(pr)}</b></li>" if pr else "") + "</ul>"
               if pf or pr else "")
    e, o = p.get("cutadaptErrorRate", 0.1), p.get("cutadaptOverlap", 3)
    if run.kind == "its":
        out.append("<h3>Primer trimming</h3>")
        nl = ncbi_lca(run)
        out.append(f"<p>Primers were removed with cutadapt from both reads, including their reverse "
                   f"complements where a short ITS amplicon reads through into the opposite primer "
                   f"(error rate {e}, minimum overlap {o}). This matters more for ITS than for 16S because "
                   f"ITS length varies widely between {'taxa' if nl else 'fungi'}.</p>" + prim_ul)
        out.append("<h3>Denoising and ASV inference</h3>")
        pool = str(sm.get("pool") or "FALSE")
        out.append(f"<p>Reads were quality-filtered with DADA2 <span class=\"mono\">filterAndTrim()</span> "
                   f"without length truncation (maxEE = {p.get('maxEE_F')} / {p.get('maxEE_R')}, no ambiguous "
                   f"bases, minimum length {esc(sm.get('minLen') or 50)} bp), error models were learned per read "
                   f"direction, ASVs were inferred with <span class=\"mono\">dada()</span>"
                   + {"TRUE": " with all samples pooled (pool = TRUE)", "pseudo": " with pseudo-pooling"}.get(pool, "")
                   + f", read pairs were merged (minimum overlap {esc(sm.get('merge_overlap') or 12)} bp, no "
                   f"mismatch) and chimeras removed with <span class=\"mono\">removeBimeraDenovo()</span> "
                   f"({esc(sm.get('chimera_method') or 'consensus')}), following the DADA2 ITS workflow.</p>")
        out.append("<h3>Taxonomic assignment</h3>")
        lc = sm.get("lca") or {}
        if nl:
            out.append(f"<p>Every ASV was compared with NCBI GenBank nt by BLASTn (blastn, NCBI BLAST service) "
                       f"and named by the lowest common ancestor (LCA) of its top {esc(lc.get('hits', 10))} hits "
                       f"that reach at least {esc(lc.get('identity', 90))}% identity over "
                       f"{esc(lc.get('coverage', 100))}% of the ASV, after DiBattista et al. (2025): the ASV is "
                       f"named down to the deepest rank all such hits share (NCBI taxonomy), and to species only "
                       f"when they all name the same one. "
                       + ("Among the hits that pass, records NCBI files under nothing below the domain "
                          "(&ldquo;uncultured eukaryote&rdquo; and other environmental sequences) were passed over "
                          "when records that do place the sequence pass too (up to "
                          f"{esc((run._json('ncbi_lca/ncbi_lca.json') or {}).get('fetched_per_asv', 50))} records "
                          "fetched per ASV for this); an ASV whose passing hits are all such records is named "
                          "Eukaryota. "
                          if (run._json("ncbi_lca/ncbi_lca.json") or {}).get("skip_uninformative") else "")
                       + f"An ASV with no hit passing is Unassigned. No curated "
                       f"reference was used"
                       + (" (SILVA, UNITE and Greengenes hold no coral ITS2; primers SCLER5.8SForw / "
                          "SCLER28SRev)" if str(pf).upper().startswith("GARTCTTTGAACGCAAATGGC") else "")
                       + ". The search needs NCBI's service: the names are those of GenBank on the day of the "
                       f"run ({esc((run._json('ncbi_lca/ncbi_lca.json') or {}).get('date', ''))}).</p>")
        else:
            out.append((_tophit_text(run, "ITS") if sm.get("tax_method") == "vsearch" else
                    f"<p>ASVs were classified with the DADA2 naive Bayesian classifier against "
                    f"<b>{esc(run.db_name)}</b> (minimum bootstrap {esc(mb)}, both orientations tried). ")
                   + "UNITE placeholders (&ldquo;unidentified&rdquo;, <i>incertae sedis</i>) are reported as "
                     "unassigned at that rank.</p>")
    elif run.kind == "cox1":
        out.append("<h3>Primer trimming</h3>")
        out.append(f"<p>Primers were removed with cutadapt (error rate {e}, minimum overlap {o}); reads "
                   f"without a primer were discarded.</p>" + prim_ul)
        out.append(f"<p>Primers were removed in two rounds (<span class=\"mono\">-n 2</span>), so a read that "
                   f"runs through the amplicon loses both primers.</p>" if run.cutadapt_two_rounds else "")
        out.append("<h3>Denoising and ASV inference</h3>")
        out.append(f"<p>Reads were truncated at {p.get('truncLen_cox1_f', '—')} / "
                   f"{p.get('truncLen_cox1_r', '—')} bases and filtered with DADA2 (maxEE = "
                   f"{p.get('maxEE_F')} / {p.get('maxEE_R')}); ASVs were inferred, read pairs merged "
                   f"(minimum overlap {sm.get('merge_overlap') or p.get('merge_overlap') or 20} bp) and "
                   f"chimeras removed.</p>")
        out.append("<h3>COI-specific curation</h3>")
        out.append(f"<p>ASVs outside {esc(sm.get('length_range') or 'the expected length range')} were "
                   f"removed, then every ASV was translated with genetic code "
                   f"{esc(sm.get('codon_table') or p.get('codon_table') or 5)} and those with a stop codon "
                   f"in every reading frame were discarded as likely nuclear pseudogenes (NUMTs"
                   + (f"; {fmt(sm['numt_removed'])} removed" if sm.get("numt_removed") is not None else "")
                   + "). "
                   + ("Remaining ASVs were curated with LULU, which merges erroneous daughter sequences "
                      "into their parent by co-occurrence"
                      + (f" ({fmt(sm['lulu_removed'])} merged)" if sm.get("lulu_removed") is not None else "")
                      + "." if run.lulu_state == "ran" else
                      "LULU curation was asked for but the LULU package is not installed on this machine, so "
                      "it was skipped: ASVs that LULU would merge into a parent are still counted separately."
                      if run.lulu_state == "missing" else
                      "LULU curation was asked for but BLAST+ (which LULU needs to compare the ASVs) was not "
                      "found, so it was skipped: ASVs that LULU would merge into a parent are still counted "
                      "separately." if run.lulu_state == "no_blast" else
                      "LULU curation was asked for but failed on this run, so it was skipped: ASVs that LULU "
                      "would merge into a parent are still counted separately." if run.lulu_state == "error" else
                      "LULU curation ran but no two ASVs matched at 84% identity or more, so nothing was merged."
                      if run.lulu_state == "empty_matchlist" else "LULU curation was not used.")
                   + "</p>")
        out.append("<h3>Taxonomic assignment</h3>")
        out.append((_tophit_text(run, "COI") if sm.get("tax_method") == "vsearch" else
                    f"<p>ASVs were classified with the DADA2 naive Bayesian classifier against "
                    f"<b>{esc(run.db_name)}</b> (minimum bootstrap {esc(mb)}). ")
                   + "MIDORI2 ranks run from phylum to species; NCBI taxon identifiers were stripped "
                     "from the names.</p>")
    elif run.kind == "pacbio":
        out.append("<h3>Primer trimming and orientation</h3>")
        out.append("<p>HiFi reads were searched for the primer pair in both orientations with cutadapt "
                   "(<span class=\"mono\">--revcomp</span>), so every read leaves this step in the "
                   "forward orientation; reads without both primers were discarded.</p>" + prim_ul)
        out.append("<h3>Denoising and ASV inference</h3>")
        out.append(f"<p>Reads of {esc(sm.get('length_range') or 'the expected length')} were filtered "
                   f"(maxEE = {p.get('pb_maxEE', 3)}), and ASVs were inferred with DADA2 using the PacBio "
                   f"error model (<span class=\"mono\">PacBioErrfun</span>, band size 32, pooling "
                   f"{esc(sm.get('pool') or p.get('pool') or 'pseudo')}), following Callahan et al. "
                   f"(2019). Chimeras were removed de novo.</p>")
        out.append("<h3>Taxonomic assignment</h3>")
        out.append(f"<p>ASVs were classified with the DADA2 naive Bayesian classifier against "
                   f"<b>{esc(run.db_name)}</b> (minimum bootstrap {esc(mb)}); species were added by exact "
                   f"match of the full-length sequence where the reference carries them.</p>")
    elif run.kind == "vsearch":
        out.append("<h3>Read preparation</h3>")
        out.append("<p>" + ("Read pairs were primer-trimmed with cutadapt and merged with VSEARCH "
                            "(<span class=\"mono\">--fastq_mergepairs</span>); merged reads with more than "
                            "one expected error were discarded (<span class=\"mono\">--fastq_filter "
                            "--fastq_maxee 1</span>). " if run.paired else
                            "Reads were quality-filtered with chopper and primer-trimmed with cutadapt, "
                            "searching both strands" + ("; both primers were removed from each read "
                            "(<span class=\"mono\">-n 2</span>)" if sm.get("cutadapt_times") == 2 else "")
                            + ". ")
                   + ("Merged reads shorter than 50 bp were discarded." if run.paired else
                      f"Reads outside {p.get('ontMinLen', '—')}&ndash;{p.get('ontMaxLen', '—')} bp were "
                      f"discarded.") + "</p>" + prim_ul)
        out.append("<h3>OTU inference</h3>")
        out.append(f"<p>Reads from all samples were pooled and dereplicated, chimeras were removed "
                   f"against the reference (<span class=\"mono\">--uchime_ref</span>), and the remaining "
                   f"unique sequences were denoised with UNOISE3 (<span class=\"mono\">--cluster_unoise</span>, "
                   f"alpha {p.get('otuUnoiseAlpha', 2.0)}, minimum size {p.get('otuUnoiseMinsize', 8)}). "
                   + (f"OTUs within {100 * float(sm['otu_merge_id']):g}% identity of a more abundant one were "
                      f"merged into it (<span class=\"mono\">--cluster_size</span>, the most abundant sequence "
                      f"kept as centroid" + (f"; {sm['otus_unoise']} OTUs before merging" if sm.get("otus_unoise") else "")
                      + "). " if sm.get("otu_merge_id") else "")
                   + f"Every read was then mapped back to the OTU centroids "
                   f"(<span class=\"mono\">--usearch_global</span>) to count it.</p>")
        out.append("<h3>Taxonomic assignment</h3>")
        out.append(f"<p>Each OTU took the taxonomy of its best hit in <b>{esc(run.db_name)}</b> "
                   f"(<span class=\"mono\">--usearch_global</span>); ranks were kept only down to the "
                   f"level the identity supports"
                   + (" &mdash; for fungal ITS species 98.41%, genus 94.3%, family 88.5%, order 81.2%, "
                      "class 80.9% and phylum 75% identity (Vu et al. 2019)"
                      if str(run.params.get("marker") or run.summary.get("marker") or "").upper().startswith("ITS")
                      else "") + ".</p>")
    return out


def _methods(run, an):
    p, out = run.params, []
    pf, pr = p.get("primer_f") or "", p.get("primer_r") or ""
    if run.kind in ("its", "cox1", "pacbio", "vsearch"):
        out += _methods_other(run, an)
    elif run.emu:
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
        _dbp = str(run.params.get("ont_db_path") or run.params.get("dbPath") or run.summary.get("db_path") or "")
        _dbp = (Path(_dbp).parent.name if Path(_dbp).suffix else Path(_dbp).name) or run.db_name
        _mf = run.summary.get("emu_minimap2_f") if run.summary.get("emu_align_speed") == "fast" else None
        if _mf:
            out.append(f"<p>The alignment used Emu's minimap2 settings with one change: minimizers found in "
                       f"more than {esc(_mf)} reference sequences were ignored "
                       f"(<span class=\"mono\">-f {esc(_mf)}</span>). Full-length 16S reads share most of "
                       f"their minimizers with tens of thousands of references, and those shared minimizers "
                       f"take most of the alignment time without telling the taxa apart; the alignment ran "
                       f"about ten times faster. Emu's expectation&ndash;maximisation then ran on that "
                       f"alignment unchanged.</p>")
            out.append('<p class="cmd">minimap2 -ax map-ont -N 50 -p 0.9 -f ' + esc(_mf)
                       + ' species_taxid.fasta reads.fastq.gz -o reads.sam</p>')
        out.append('<p class="cmd">emu abundance --type map-ont --db ' + esc(_dbp)
                   + f' --min-abundance {esc(run.summary.get("emu_min_abundance") or 0.0001)}'
                   + (' --keep-counts' if run.emu_counts == "estimated_reads" else "")
                   + (' reads.sam' if _mf else "") + '</p>')
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
            n2 = " -n 2" if run.cutadapt_two_rounds else ""
            cmd = (f"cutadapt -g {pf} -a {_rc(pr)} -G {pr} -A {_rc(pf)} -e {e} -O {o}{n2} -m 50{disc} "
                   f"-o R1_cut.fastq.gz -p R2_cut.fastq.gz R1.fastq.gz R2.fastq.gz" if paired else
                   f"cutadapt -g {pf} -a {_rc(pr)} -e {e} -O {o}{n2} -m 50{disc} -o cut.fastq.gz in.fastq.gz")
            if n2:
                out.append("<p>cutadapt ran two rounds (<span class=\"mono\">-n 2</span>): a read that runs "
                           "through an amplicon shorter than itself holds both primers and the adapter behind "
                           "them, and both rounds are needed to remove them.</p>")
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
        if run.blast:
            out += _blast_methods(run)
            return out + _diversity_methods(run, an)
        pr2 = getattr(run, "is_pr2", False)
        out.append(f"<p>ASVs were classified with the DADA2 implementation of the naive Bayesian "
                   f"classifier (<span class=\"mono\">assignTaxonomy()</span>) against "
                   f"<b>{esc(run.db_name)}</b>, with a minimum bootstrap confidence of "
                   f"{esc(p.get('minBoot', 50))}. "
                   + ("PR2 (Guillou et al. 2013) is the reference for eukaryotic 18S; its ranks are Domain, "
                      "Supergroup, Division (reported here as Phylum), Subdivision, Class, Order, Family, Genus and "
                      "Species, and its species come from the same classifier and bootstrap threshold. Names PR2 "
                      "uses for a rank it leaves undefined (ending in &ldquo;_X&rdquo;) are reported as "
                      "unclassified at that rank. ASVs whose lineage names them mitochondria or chloroplasts were "
                      "set aside from the community figures; ASVs PR2 places outside the Eukaryota stay in the "
                      "tables under their own Domain.</p>" if pr2 else
                      "Species were added where the reference carries them, by exact matching "
                      "(<span class=\"mono\">addSpecies()</span>). Organelle sequences (mitochondria, chloroplasts)"
                      + (" were identified from the assigned lineage and set aside from the community figures.</p>"
                         if run.euk_target else
                         " and eukaryotic ASVs were identified from the assigned lineage and set aside from the "
                         "community figures.</p>")))
    return out + _diversity_methods(run, an)


def _blast_methods(run) -> list:
    """How a 12S / CytB run was named: BLASTn, ties, identity thresholds, limit, NCBI check."""
    bi = run.blast_info or {}
    th = bi.get("thresholds") or {"Species": 98.5, "Genus": 95, "Family": 90, "Order": 85, "Class": 80}
    lim = str(bi.get("limit") or run.params.get("blast_limit") or "fish")
    lim_txt = {"fish": "ray-finned and cartilaginous fishes (NCBI txid 7898 and 7777)",
               "vertebrates": "vertebrates (NCBI txid 7742)", "all": "no taxonomic limit"}.get(lim, lim)
    ver = (bi.get("reference") or {}).get("versions") or {}
    ref = " + ".join(f"{k} {v}" for k, v in ver.items()) or "MIDORI2 + MitoFish"
    ncbi = str(bi.get("ncbi") or "off") != "off"
    out = [f"<p>Each ASV was named by nucleotide BLAST ({esc(bi.get('blastn') or 'BLASTn')}) against a "
           f"fish mitochondrial reference built from {esc(ref)}. Hits within 0.5% identity of the best "
           f"hit count as ties, and the name is their lowest common ancestor, cut to the deepest rank "
           f"whose identity threshold the best hit meets: species {th.get('Species')}%, genus "
           f"{th.get('Genus')}%, family {th.get('Family')}%, order {th.get('Order')}% and class "
           f"{th.get('Class')}%. A species is never named below its threshold; an ASV that meets none of "
           f"them keeps only its phylum. Target taxa are {esc(lim_txt)}; an ASV whose best hit lies outside "
           f"them (other animals, algae, human DNA) keeps its own name and is counted as non-target, never "
           f"forced onto a fish.</p>"]
    if ncbi:
        out.append("<p>Every ASV was also searched against NCBI nt (remote BLASTn, same taxonomic limit). "
                   "Where NCBI and the local reference disagree, the name is cut to the rank they share and "
                   "the ASV is flagged for checking.</p>")
    _td, _tb = _cite(run, "taxonomy_table_dada2.csv"), _cite(run, "blast/taxonomy_blast.csv")
    out.append("<p>The DADA2 naive Bayesian classification against "
               f"<b>{esc(run.db_name)}</b> is kept beside it" + (f" ({esc(_td)})" if _td else "")
               + "; the BLAST names are used in every table and figure of this report"
               + (f" ({esc(_tb)})" if _tb else "") + ".</p>")
    return out


def _cite(run, name: str):
    """A file of the run folder as the reader has it: its name in a Download results package when the report
    is built for one (run.pkg_files, from build_profile_html(files=...)) — None when the package has none —
    else the name in the run folder."""
    m = getattr(run, "pkg_files", None)
    if m is None:
        return name if (Path(run.path) / name).exists() else None
    v = m.get(name)
    return v[0] if v else None


def _blast_name(r: dict):
    """(deepest named rank, name) of a BLAST row; ('', '') when nothing was named."""
    for k in ("Species", "Genus", "Family", "Order", "Class", "Phylum"):
        v = str(r.get(k) or "").strip()
        if v and v.upper() != "NA":
            return k, v
    return "", ""


def _group_section(run, an, alpha, S, F, T) -> list:
    """Analysis by group (metadata.csv): which taxa each group holds, its alpha diversity
    (Kruskal–Wallis) and how the groups differ in composition (PERMANOVA, distances)."""
    gl = [(g, [x for x in v if x in S]) for g, v in an.group_list()]
    gl = [(g, v) for g, v in gl if v]
    if len(gl) < 2 or not any(len(v) >= 2 for _g, v in gl):
        return []
    col = getattr(an, "group_col", "") or "group"
    ung = [x for x in S if x not in an.groups]
    out = ['<h2 class="pbreak">Analysis by group</h2>',
           f"<p>The samples were grouped by the column <b>{esc(col)}</b> of the job&rsquo;s "
           f"{esc(_cite(run, 'metadata.csv') or 'metadata')} "
           f"({len(gl)} groups"
           + (f"; {len(ung)} sample{'s' if len(ung) != 1 else ''} without a value "
              f"{'are' if len(ung) != 1 else 'is'} left out of the comparisons" if ung else "")
           + "). Group means are the means of the samples&rsquo; relative abundances, so every sample "
           "counts equally whatever its depth.</p>"]
    tot = {x: an.comm_tot.get(x, 0) for x in S}
    out.append(T(["Group", "Samples", "Mean reads", "Samples in the group"],
                 [[esc(g), str(len(v)), fmt(sum(tot[x] for x in v) / len(v)), _left(esc(", ".join(v)))] for g, v in gl],
                 f"Groups of the column &ldquo;{esc(col)}&rdquo;."))

    # taxa by group
    out.append("<h3>Taxa by group</h3>")
    gen = an.rank_pct("Genus")
    names = [t for t in gen if t not in ("Unclassified", "Other")]
    mean = {t: {g: sum(gen[t].get(x, 0) for x in v) / len(v) for g, v in gl} for t in names}
    pres = {t: {g: sum(1 for x in v if gen[t].get(x, 0) > 0) for g, v in gl} for t in names}
    in_all = [t for t in names if all(pres[t][g] > 0 for g, _v in gl)]
    only = {g: sorted((t for t in names if pres[t][g] * 2 >= len(v) and pres[t][g] > 0
                       and all(pres[t][h] == 0 for h, _w in gl if h != g)), key=lambda t: -mean[t][g])
            for g, v in gl}
    det = [t for t in names if any(pres[t][g] for g, _v in gl)]
    out.append(f"<p>{len(det)} genera were detected; {len(in_all)} of them in every group. "
               + ("Genera found only in one group (in at least half of its samples and in no sample of another "
                  "group) are listed below; " if any(only.values()) else
                  "No genus was found consistently in one group and absent from all the others; ")
               + "the top genera of each group follow, with the number of its samples holding each.</p>")
    out.append(T(["Group", "Top genera: mean % (samples holding it)"],
                 [[esc(g), _left("; ".join(f"<i>{esc(t)}</i> {mean[t][g]:.1f}% ({pres[t][g]}/{len(v)})"
                                           for t in sorted(names, key=lambda t: -mean[t][g])[:10] if mean[t][g] > 0))]
                  for g, v in gl], "The ten most abundant genera of each group."))
    if any(only.values()):
        out.append(T(["Group", "Genera", "Found only in this group: mean % (samples holding it)"],
                     [[esc(g), str(len(only[g])),
                       _left("; ".join(f"<i>{esc(t)}</i> {mean[t][g]:.2f}% ({pres[t][g]}/{len(v)})" for t in only[g][:12])
                        + (" &hellip;" if len(only[g]) > 12 else "") or "&mdash;")] for g, v in gl],
                     "Genera present in at least half of one group&rsquo;s samples and in no sample of any other group."))
    for rank in ("Phylum", "Genus", "Species"):
        if rank not in run.rank_cols:
            continue
        pct = an.rank_pct(rank)
        gt = {t: {g: sum(pct[t].get(x, 0) for x in v) / len(v) for g, v in gl} for t in pct}
        named = sorted((t for t in gt if t != "Unclassified"), key=lambda t: -sum(gt[t].values()))
        tbl_ = {t: gt[t] for t in named[:30]}
        if named[30:]:
            tbl_["Other"] = {g: sum(gt[t][g] for t in named[30:]) for g, _v in gl}
        if "Unclassified" in gt:
            tbl_["Unclassified"] = gt["Unclassified"]
        out.append(F(svg_stacked(an, rank, 30, samples=[g for g, _v in gl], table=tbl_),
                     f"Top 30 {PLURAL.get(rank, rank.lower())} &mdash; mean relative abundance per group."))

    # alpha by group
    idx = ["Observed", "Chao1", "Shannon", "Simpson"]
    out.append('<h3 class="pbreak">Alpha diversity by group</h3>')
    kw = {i_: _kruskal([[alpha[x][i_] for x in v] for _g, v in gl]) for i_ in idx}
    sig = [i_ for i_ in idx if kw[i_] and kw[i_][2] < 0.05]
    out.append("<p>Alpha diversity of the samples in each group, compared with the Kruskal&ndash;Wallis test "
               "(rank-based, no assumption of normality). "
               + (f"The groups differ in {', '.join(sig)} (p &lt; 0.05); the test says that at least one group "
                  "differs, not which one &mdash; the boxes show where the differences lie. " if sig else
                  "No index differs significantly between the groups (all p &ge; 0.05). ")
               + "Richness (Observed, Chao1) grows with sequencing depth, so groups of very different depth "
               "are compared with care."
               + (" Chao1 equals Observed in every sample here: the counts hold no singletons (Emu removes taxa "
                  "below its minimum abundance), so Chao1 has nothing to extrapolate from and its test repeats "
                  "Observed&rsquo;s." if all(abs(alpha[x]["Chao1"] - alpha[x]["Observed"]) < 1e-9
                                             for _g, v in gl for x in v) else "")
               + (" With two groups the Kruskal&ndash;Wallis test is equivalent to the Mann&ndash;Whitney U test; "
                  "it depends only on the samples&rsquo; ranks, so two indices that order the samples alike "
                  "(Shannon and Simpson often do) give the same result." if len(gl) == 2 else "")
               + "</p>")
    out.append(F(svg_alpha_groups(gl, alpha, idx), "Alpha diversity per group: box = quartiles, line = median, "
                 "dots = samples."))

    def ms(vs, d):
        if not vs:
            return "&mdash;"
        m = sum(vs) / len(vs)
        sd = (sum((v - m) ** 2 for v in vs) / (len(vs) - 1)) ** 0.5 if len(vs) > 1 else 0.0
        return f"{m:.{d}f} &plusmn; {sd:.{d}f}"
    out.append(T(["Group", "n"] + idx,
                 [[esc(g), str(len(v))] + [ms([alpha[x][i_] for x in v], 3 if i_ in ("Shannon", "Simpson") else 1)
                                           for i_ in idx] for g, v in gl]
                 + [["Kruskal&ndash;Wallis", ""] + [(f"H = {kw[i_][0]:.2f}, p = {kw[i_][2]:.3g}" if kw[i_] else "&mdash;")
                                                    for i_ in idx]],
                 "Mean &plusmn; SD per group; last row: Kruskal&ndash;Wallis test across the groups."))

    # beta by group
    out.append("<h3>Beta diversity by group</h3>")
    D = an.bray()
    ix = {x: i for i, x in enumerate(an.samples)}
    pg = an.permanova_groups()
    rows = []
    for g, v in gl:
        wi = [D[ix[a]][ix[b]] for i, a in enumerate(v) for b in v[i + 1:]]
        bt = [D[ix[a]][ix[b]] for a in v for h, w in gl if h != g for b in w]
        rows.append([esc(g), (f"{sum(wi)/len(wi):.3f}" if wi else "&mdash;"),
                     (f"{sum(bt)/len(bt):.3f}" if bt else "&mdash;")])
    out.append('<div style="break-inside: avoid">')
    out.append("<p>" + (f"PERMANOVA on Bray&ndash;Curtis dissimilarity ({pg['perms']} permutations) gave "
                        f"R&sup2; = {pg['R2']:.3f}, pseudo-F = {pg['F']:.2f}, p = {pg['p']:.3g}: the grouping explains "
                        f"{100*pg['R2']:.1f}% of the variation in composition"
                        + (", a significant difference between the groups. " if pg["p"] < 0.05 else
                           ", which is not a significant difference. ")
                        if pg else "PERMANOVA needs two groups and replicated samples. ")
               + "PERMANOVA is also sensitive to groups that differ in spread rather than in centre, so the table "
               "gives each group&rsquo;s mean dissimilarity within itself and to the other groups; the PCoA under "
               "Beta diversity is coloured by group.</p>")
    out.append(T(["Group", "Mean Bray&ndash;Curtis within the group", "Mean to the other groups"], rows,
                 "Lower within than between: the group&rsquo;s samples resemble each other more than they "
                 "resemble other groups.") + "</div>")
    return out


def _left(html_):
    """A table cell of running text: left-aligned (the report's cells are right-aligned numbers)."""
    return f'<div style="text-align:left">{html_}</div>'


def _ncbi_un_section(run, S, T) -> list:
    """OTUs the reference could not name, checked at NCBI nt (ncbi_unassigned.py): what they are by
    group per sample, and why their reads are no profile of those organisms."""
    info, rows = run.ncbi_un_info or {}, run.ncbi_un or {}
    out = ["<h3>Sequences the reference could not name (NCBI check)</h3>"]
    if info.get("error") and not info.get("checked"):
        return out + [f"<p>The OTUs the reference could not name were to be compared with NCBI nt, but NCBI could "
                      f"not be reached ({esc(str(info['error'])[:160])}); they stay 'Unassigned'.</p>"]
    if not rows:
        return out + ["<p>Every OTU had a hit in the reference; nothing was sent to NCBI.</p>"]
    groups = ["Bacteria", "Archaea", "Fungi", "Plants", "Animals", "Other eukaryotes", "Viruses",
              "Unclassified record", "No NCBI match"]
    tot = {s_: sum(r.get(s_, 0.0) for r in run.mat.values()) or 1.0 for s_ in S}
    per = {}
    named = {}
    for q, row in run.mat.items():
        r = rows.get(str(q).split(";")[0])
        if not r:
            continue
        g = r.get("group") or "Unclassified record"
        for s_ in S:
            per.setdefault(g, {}).setdefault(s_, 0.0)
            per[g][s_] += row.get(s_, 0.0)
        if r.get("rank") in ("genus", "family", "phylum"):
            named.setdefault(g, Counter())[r.get("named_as", "")] += sum(row.get(s_, 0.0) for s_ in S)
    fl = info.get("rank_floors") or {}
    out.append(f"<p>{info.get('checked', len(rows))} OTUs had no hit in the reference and were compared with NCBI nt "
               f"(BLASTn megablast, {esc(info.get('date', ''))}). Each is named to genus at &ge;{fl.get('genus', 95):g}% "
               f"identity, family &ge;{fl.get('family', 90):g}%, phylum &ge;{fl.get('phylum', 80):g}% over at least "
               f"{info.get('min_cov', 80):g}% of the OTU when the equally good records agree, otherwise to its group"
               + (f"; {info['left_out']} OTUs with fewer reads were not sent" if info.get("left_out") else "")
               + ". The primers are made for the target group; other DNA is amplified only where they happen to "
               "bind, so these reads show what is present, not how much of it, and organisms the primers do not "
               "bind are missed. A bacterial community profile needs 16S rRNA sequencing.</p>")
    g_in = [g for g in groups if g in per]
    out.append(T(["Group"] + [esc(s_) for s_ in S] + ["Named (most reads first)"],
                 [[esc(g)] + [f"{100 * per[g][s_] / tot[s_]:.2f}%" for s_ in S]
                  + [esc(", ".join(n for n, _ in named.get(g, Counter()).most_common(6)))] for g in g_in],
                 "Reads of the OTUs the reference could not name, by the group of their closest NCBI record "
                 "(% of each sample's counted reads)"
                 + (f"; the list per OTU is {esc(_cite(run, 'ncbi_unassigned.tsv'))}." if _cite(run, "ncbi_unassigned.tsv")
                    else ".")))
    return out


def _ncbi_lca_section(run, S, T) -> list:
    """ITS named at NCBI (ncbi_lca.py): how far the names reach, in reads, and the most abundant ASVs
    with the hit behind each name — what a reader needs to judge a BLAST + LCA name."""
    import csv as _csv
    p = Path(run.path) / "ncbi_lca" / "taxonomy_ncbi_lca.csv"
    if not ncbi_lca(run) or not p.exists():
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        rows = {r["sequence"]: r for r in _csv.DictReader(fh)}
    reads = {q: sum(run.mat.get(q, {}).get(s_, 0.0) for s_ in S) for q in rows}
    used = {q: n for q, n in reads.items() if n > 0}
    tot = sum(used.values()) or 1
    st, rk = Counter(), Counter()
    for q, n in used.items():
        st[rows[q].get("status", "")] += n
        rk[rows[q].get("rank_reached", "none")] += n
    lc = (run.summary or {}).get("lca") or {}
    pc = lambda n: f"{100 * n / tot:.1f}%"
    deep = lambda *ks: sum(rk[k] for k in ks)
    out = ["<h3>Names from NCBI GenBank (BLASTn + LCA)</h3>",
           f"<p>{len([q for q in used if rows[q].get('status') == 'named'])} of {len(used)} ASVs were named "
           f"({pc(st['named'])} of the reads): to species {pc(deep('Species'))}, to genus or deeper "
           f"{pc(deep('Species', 'Genus'))}, to family or deeper {pc(deep('Species', 'Genus', 'Family'))}. "
           f"Unassigned: {pc(st['below_threshold'])} of the reads had NCBI hits, none at "
           f"&ge; {esc(lc.get('identity', 90))}% identity over {esc(lc.get('coverage', 100))}% of the ASV"
           + (f"; {pc(st['no_hit'])} had no NCBI hit" if st["no_hit"] else "")
           + (f"; {pc(st['not_checked'])} could not be checked (NCBI batch failed)" if st["not_checked"] else "")
           + ".</p>"]
    trs = []
    for a in run.top_asvs(15):
        r = rows.get(a["seq"])
        if not r:
            continue
        rank = r.get("rank_reached") or "none"
        name = (r.get(rank) if rank in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
                else "") or "Unassigned"
        nm = f"<i>{esc(name)}</i>" if rank in ("Genus", "Species") else esc(name)
        rl = "domain" if rank == "Kingdom" and name in ("Eukaryota", "Bacteria", "Archaea") else rank.lower()
        trs.append([f'<span class="mono">{esc(a["id"])}</span>', f'{a["pct"]:.2f}%',
                    nm + (f' <span class="dim">({esc(rl)})</span>' if rank not in ("none", "") else ""),
                    f"<i>{esc(r.get('best_hit') or '—')}</i>",
                    (f"{float(r['best_pident']):.1f}% / {float(r['best_qcov']):.0f}%"
                     if r.get("best_pident") not in (None, "") else "—"),
                    f"{r.get('n_passing') or 0} of {r.get('n_hits') or 0}"])
    cite = _cite(run, "ncbi_lca/taxonomy_ncbi_lca.csv")
    out.append(T(["ASV", "Reads", "Named as", "Best NCBI hit", "Identity / coverage", "Hits used / fetched"], trs,
                 "The most abundant ASVs: the name the LCA of the passing hits supports, and the best hit "
                 "behind it" + (f"; every ASV with all its hits is in {esc(cite)}" if cite else "") + "."))
    return out


def _blast_section(run, S, T) -> list:
    """Target (fish) and non-target reads per sample, the fish named, and what was not fish —
    the first question of a fish eDNA run, answered before the community figures."""
    bi = run.blast_info or {}
    th = bi.get("thresholds") or {}
    cls_th = th.get("Class", 80)
    lim = str(bi.get("limit") or run.params.get("blast_limit") or "fish")
    tword = {"fish": "fish", "vertebrates": "vertebrate", "all": "target"}.get(lim, "target")
    tot = {s_: 0.0 for s_ in S}
    part = {k: {s_: 0.0 for s_ in S} for k in ("yes", "far", "no", "none")}
    fish, other = {}, {}
    flagged = []
    agree = Counter()
    for q, row in run.mat.items():
        r = run.blast.get(q) or {}
        t = str(r.get("target") or "").strip()
        t = t if t in ("yes", "no") else "none"
        if t == "none" and r.get("agreement") == "ncbi_only" and (r.get("Genus") or r.get("Species")):
            t = "yes"          # named by NCBI, searched inside the limit (runs before v2.25.2 left it blank)
        reads = sum(row.get(s_, 0.0) for s_ in S)
        if r.get("agreement"):
            agree[r["agreement"]] += 1
        rank, name = _blast_name(r)
        try:
            pid = float(r.get("pident") or "nan")
        except ValueError:
            pid = float("nan")
        far = not (pid == pid and pid >= float(cls_th))
        k_ = "far" if (t == "yes" and far) else t
        for s_ in S:
            tot[s_] += row.get(s_, 0.0)
            part[k_][s_] += row.get(s_, 0.0)
        if t == "yes":
            if rank in ("", "Phylum"):
                rank, name = "", (f"Closest to a fish, below {cls_th}% identity (unconfirmed)" if far else
                                  "Fish, name cut to phylum (conflicting hits or NCBI disagreement)")
            e = fish.setdefault(name, {"rank": rank, "reads": 0.0, "pid": [], "in": set()})
        elif t == "no" and not name:
            bh = str(r.get("best_hit") or r.get("ncbi_best") or "").strip()
            name, rank = (f"Not named (closest: {bh[:70]})" if bh else "Not named (no close hit)"), ""
            e = other.setdefault(name, {"rank": rank, "reads": 0.0, "pid": [], "in": set()})
        elif t == "no" and name:
            if rank == "Phylum":
                name = (f"{name} (best hit below {cls_th}% identity)" if far else
                        f"{name} (local and NCBI best hits disagree; flagged)"
                        if r.get("agreement") == "differ" else f"{name} (cut to phylum: conflicting hits)")
            e = other.setdefault(name, {"rank": rank, "reads": 0.0, "pid": [], "in": set()})
        else:
            e = None
        if e is not None:
            e["reads"] += reads
            if pid == pid:
                e["pid"].append(pid)
            e["in"].update(s_ for s_ in S if row.get(s_, 0.0) > 0)
        if str(r.get("check") or "").lower() == "yes":
            flagged.append((name or "unnamed", reads))
    pct = lambda k, s_: (100 * part[k][s_] / tot[s_]) if tot[s_] else 0.0
    all_tot = sum(tot.values()) or 1.0
    out = [f'<h3 class="pbreak">{tword.capitalize()} and non-target reads (BLAST)</h3>']
    fs = [(s_, pct("yes", s_)) for s_ in S]
    hi, lo = max(fs, key=lambda kv: kv[1]), min(fs, key=lambda kv: kv[1])
    far_tot = sum(part["far"].values())
    out.append(f"<p>{tword.capitalize()} DNA identified at {cls_th}% identity or more made up "
               f"{100*sum(part['yes'].values())/all_tot:.1f}% of the reads over the whole run, from "
               f"{lo[1]:.1f}% ({esc(lo[0])}) to {hi[1]:.1f}% ({esc(hi[0])}) per sample."
               + (f" A further {100*far_tot/all_tot:.1f}% had a {tword} as best hit but below {cls_th}% "
                  f"identity: too distant to confirm as {tword}, so it is reported apart." if far_tot else "")
               + f" The rest is DNA of other organisms the primers also amplified "
               f"({100*sum(part['no'].values())/all_tot:.1f}%) and sequences without a match in the reference "
               f"({100*sum(part['none'].values())/all_tot:.1f}%). Composition figures further on include all "
               f"of it; the {tword} share is the part this marker targets.</p>")
    out.append(T(["Sample", "Reads", f"{tword.capitalize()} %"]
                 + ([f"{tword.capitalize()}-like, &lt;{cls_th}% %"] if far_tot else [])
                 + ["Non-target %", "No match %"],
                 [[esc(s_), fmt(tot[s_]), f"{pct('yes', s_):.1f}"]
                  + ([f"{pct('far', s_):.1f}"] if far_tot else [])
                  + [f"{pct('no', s_):.1f}", f"{pct('none', s_):.1f}"] for s_ in S],
                 f"Reads per sample by BLAST result: {tword} (within the BLAST limit, identity at least "
                 f"{cls_th}%)" + (f", best hit a {tword} but below {cls_th}%" if far_tot else "")
                 + ", other organisms, and no match in the reference."))
    def _rows(d, n):
        rows = []
        for name, e in sorted(d.items(), key=lambda kv: -kv[1]["reads"])[:n]:
            pr = (f"{min(e['pid']):.1f}&ndash;{max(e['pid']):.1f}" if e["pid"] and min(e["pid"]) != max(e["pid"])
                  else f"{e['pid'][0]:.1f}" if e["pid"] else "&mdash;")
            nm = f"<i>{esc(name)}</i>" if e["rank"] in ("Species", "Genus") else esc(name)
            rows.append([nm, esc(e["rank"].lower()) or "&mdash;", fmt(e["reads"]),
                         f"{100*e['reads']/all_tot:.2f}%", pr, str(len(e["in"]))])
        return rows
    if fish:
        out.append(T(["Taxon", "Named to", "Reads", "All", "Best-hit identity %", "Samples"], _rows(fish, 40),
                     f"{tword.capitalize()} taxa, named to the deepest rank whose identity threshold was met. "
                     f"Sequences whose best hit is a {tword} below {cls_th}% identity get no name: they may be a "
                     f"{tword} missing from the reference, or another organism whose closest reference "
                     f"happens to be a {tword}."))
    hs_t = sum(e["reads"] for k, e in fish.items() if k.startswith("Homo"))
    if hs_t:
        out.append(f"<p>Human DNA (<i>Homo sapiens</i>) made up {fmt(hs_t)} reads ({100*hs_t/all_tot:.1f}%). "
                   f"With the BLAST limit set to {tword}s it counts as {tword} DNA in the figures above; it "
                   f"usually comes from handling and is not part of the {tword} community.</p>")
    if other:
        hs = sum(e["reads"] for k, e in other.items() if k.startswith("Homo"))
        out.append(T(["Taxon", "Named to", "Reads", "All", "Best-hit identity %", "Samples"], _rows(other, 15),
                     "The most abundant non-target taxa. A name at phylum level comes from a distant best hit "
                     "and is approximate."))
        if hs:
            out.append(f"<p>Human DNA (<i>Homo sapiens</i>) made up {fmt(hs)} reads ({100*hs/all_tot:.2f}%); "
                       f"it is non-target and does not affect the {tword} names.</p>")
    if agree:
        lab = {"agree": "agree", "differ": "differ (cut to the shared rank, flagged)", "local_better":
               "local reference closer", "ncbi_only": "named by NCBI only", "ncbi_tie": "NCBI ties",
               "ncbi_no_hit": "no NCBI hit within the limit", "not_checked": "not checked",
               "outside_limit": f"non-target kept (NCBI, searched for {tword} only, had no closer {tword})"}
        out.append("<p><b>NCBI check.</b> " + "; ".join(f"{lab.get(k, esc(k))}: {v} {run.feat}"
                                                        for k, v in agree.most_common()) + ". "
                   + (f"{len(flagged)} {run.feat} are flagged for checking ("
                      + ", ".join(f"{esc(n)}" + (f" &times;{c}" if c > 1 else "")
                                  for n, c in Counter(n for n, _r in flagged).most_common(6))
                      + ")" + (f"; see the check column of {esc(_cite(run, 'blast/taxonomy_blast.csv'))}."
                               if _cite(run, "blast/taxonomy_blast.csv") else ".")
                      if flagged else "No ASV is flagged for checking.") + "</p>")
    return out


def _diversity_methods(run, an) -> list:
    out = []
    out.append("<h3>Diversity analysis</h3>")
    unit = run.feat
    out.append(f"<p>Alpha diversity (Observed {unit}, Chao1, ACE, Shannon and Simpson) was computed "
               f"from the per-sample counts. Rarefaction curves give the expected number of {unit} in a "
               f"random subsample of each depth (Hurlbert&rsquo;s formula, as in vegan&rsquo;s "
               f"<span class=\"mono\">rarefy()</span>). "
               + ("Beta diversity was measured as Bray&ndash;Curtis dissimilarity on relative abundances "
                  "(heatmap; PCoA and UPGMA need three samples or more). "
                  if len(an.samples) == 2 else
                  "Beta diversity was measured as Bray&ndash;Curtis "
                  "dissimilarity on relative abundances and visualised by principal coordinates analysis "
                  "(PCoA), average-linkage (UPGMA) clustering and a dissimilarity heatmap. "
                  if len(an.samples) > 1 else
                  "Beta diversity compares samples and was not computed for this single sample. ")
               + f"Taxonomic "
               f"composition is shown as the top 30 taxa at each rank per sample; the remaining taxa are "
               f"pooled as &ldquo;Other&rdquo; and reads without an assignment at that rank as "
               f"&ldquo;Unclassified&rdquo;."
               + (" Because there is one sample per condition, comparisons are descriptive; no "
                  "statistical test between groups is possible without replicates."
                  if len(an.samples) > 1 and an.groups and len(set(an.groups.values())) >= len(an.samples)
                  else (" Where samples are grouped (metadata.csv, also set after the run under Customise charts), "
                        if getattr(run, "pkg_files", None) is None else
                        f" Where samples are grouped ({esc(_cite(run, 'metadata.csv') or 'the job metadata')}), ")
                       + "groups are compared by PERMANOVA on Bray&ndash;Curtis (999 permutations, "
                       "Anderson 2001) and their alpha diversity by the Kruskal&ndash;Wallis test."
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


def _picrust_result(pc, F, T, top=25, run=None):
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
                 "difference is.", (f"{_cite(run, 'PICRUSt2')}/pathways_out"
                                    if run is not None and _cite(run, "PICRUSt2") else
                                    None if run is not None else "PICRUSt2/pathways_out")))
    return out


def svg_qprofile(qc, trunc=None, sugg=None) -> str:
    """Median quality per position with the 25-75% band, R1 and R2 side by side."""
    pooled = qc.get("pooled") or {}
    panels = [(k, (pooled.get(k) or {}).get("per_position") or []) for k in ("R1", "R2")]
    panels = [(k, p) for k, p in panels if p]
    if not panels:
        return ""
    W, H, ml, mr, mt, mb, gap = 560, 210, 34, 8, 18, 30, 26
    pw = (W - ml - mr - gap * (len(panels) - 1)) / len(panels)
    ph = H - mt - mb
    p = _svg_open(W, H)
    for j, (name, prof) in enumerate(panels):
        x0 = ml + j * (pw + gap)
        L = len(prof)
        X = lambda i: x0 + pw * i / max(1, L - 1)
        Y = lambda q: mt + ph * (1 - min(41, max(0, q)) / 41)
        for qv in (20, 30):
            p.append(f'<line x1="{x0:.1f}" x2="{x0+pw:.1f}" y1="{Y(qv):.1f}" y2="{Y(qv):.1f}" '
                     f'stroke="{RULE}" stroke-dasharray="3 3"/>')
            if j == 0:
                p.append(f'<text x="{x0-4:.1f}" y="{Y(qv)+3:.1f}" font-size="7.5" fill="{FAINT}" '
                         f'text-anchor="end">Q{qv}</text>')
        band = [f"{X(i):.1f},{Y(r[3]):.1f}" for i, r in enumerate(prof)] + \
               [f"{X(i):.1f},{Y(r[1]):.1f}" for i, r in reversed(list(enumerate(prof)))]
        p.append(f'<polygon points="{" ".join(band)}" fill="{ACC}" opacity="0.18"/>')
        p.append('<polyline fill="none" stroke="' + ACC + '" stroke-width="1.4" points="'
                 + " ".join(f"{X(i):.1f},{Y(r[2]):.1f}" for i, r in enumerate(prof)) + '"/>')
        p.append(f'<rect x="{x0:.1f}" y="{mt}" width="{pw:.1f}" height="{ph}" fill="none" stroke="{RULE}"/>')
        p.append(f'<text x="{x0:.1f}" y="{mt-5}" font-size="8.6" fill="{INK}" font-weight="700">{name}</text>')
        for tick in range(0, L + 1, 50):
            p.append(f'<text x="{X(min(tick, L-1)):.1f}" y="{H-mb+12}" font-size="7.4" fill="{FAINT}" '
                     f'text-anchor="middle">{tick}</text>')
        for val, col, lab in ((trunc[j] if trunc else None, STOP, "used"),
                              (sugg[j] if sugg else None, "#2E6B46", "suggested")):
            if val:
                xv = X(min(L - 1, val))
                p.append(f'<line x1="{xv:.1f}" x2="{xv:.1f}" y1="{mt}" y2="{mt+ph}" stroke="{col}" '
                         f'stroke-width="1.1" stroke-dasharray="4 2"/>')
                p.append(f'<text x="{xv-3:.1f}" y="{mt+10 + (0 if lab == "used" else 10)}" font-size="7.4" '
                         f'fill="{col}" text-anchor="end">{lab} {val}</text>')
    p.append(f'<text x="{W/2:.0f}" y="{H-4}" font-size="8" fill="{MUT}" text-anchor="middle">'
             f'position in read (bp)</text>')
    p.append("</svg>")
    return "".join(p)


def svg_longread(qc) -> str:
    """Read length and mean read quality, pooled over samples (NanoPlot-style)."""
    lh = [0] * 40
    qh = [0] * 41
    for ent in (qc.get("samples") or {}).values():
        for f in ent.get("files", []):
            for i, v in enumerate(f.get("len_hist_100bp") or []):
                lh[i] += v
            for i, v in enumerate(f.get("meanq_hist") or []):
                qh[i] += v
    if not sum(lh):
        return ""
    W, H, mt, mb = 560, 170, 16, 28
    p = _svg_open(W, H)
    for j, (hist, lab, unit, step) in enumerate(((lh, "Read length", "kb", 100), (qh, "Mean read quality", "Q", 1))):
        x0, pw, ph = 30 + j * 280, 240, H - mt - mb
        last = max(i for i, v in enumerate(hist) if v) + 1
        hist = hist[:max(last, 5)]
        mx = max(hist) or 1
        bw = pw / len(hist)
        for i, v in enumerate(hist):
            h = ph * v / mx
            p.append(f'<rect x="{x0+i*bw:.1f}" y="{mt+ph-h:.1f}" width="{max(bw-0.6,0.6):.1f}" '
                     f'height="{h:.1f}" fill="{ACC}"/>')
        p.append(f'<text x="{x0}" y="{mt-4}" font-size="8.6" font-weight="700" fill="{INK}">{lab}</text>')
        for i in range(0, len(hist) + 1, max(1, len(hist) // 5)):
            t = f"{i*step/1000:g}" if unit == "kb" else f"{i}"
            p.append(f'<text x="{x0+i*bw:.1f}" y="{H-mb+11}" font-size="7.4" fill="{FAINT}" '
                     f'text-anchor="middle">{t}</text>')
        p.append(f'<text x="{x0+pw/2:.0f}" y="{H-4}" font-size="7.8" fill="{MUT}" '
                 f'text-anchor="middle">{"kb" if unit == "kb" else "Phred Q"}</text>')
    p.append("</svg>")
    return "".join(p)


def _rawqc_caption(qc: dict) -> str:
    """What Table 'Raw read quality' looked at: the first n reads of each file, or (MinKNOW chunks,
    read_qc.py since v2.43.3) n reads per sample spread over all of its files."""
    n = fmt(qc.get("reads_per_file_examined", 0))
    ch = [(k, v["files"][0]["files_in_sample"]) for k, v in (qc.get("samples") or {}).items()
          if v.get("files") and v["files"][0].get("files_in_sample")]
    if not ch:
        return f"Raw reads before any trimming: the first {n} reads of each file"
    many = len(ch) > 4
    which = (f"{len(ch)} samples of {min(c for _k, c in ch)}–{max(c for _k, c in ch)} files each" if many else
             "; ".join(f"{esc(k)}: {c} files" for k, c in ch))
    return (f"Raw reads before any trimming: {n} reads per sample, an even share from each of its MinKNOW "
            f"chunk files ({which})")


def _readqc_section(run, F, T):
    qc = run._json("qc/read_qc.json")
    if not qc or not qc.get("samples"):
        return []
    out = ["<h3>Raw read quality</h3>"]
    rows = []
    keep = {_skey(x) for x in run.samples}
    for name, ent in qc["samples"].items():
        if run.only and _skey(name) not in keep:
            continue
        fs = ent.get("files", [])
        if not fs:
            continue
        f1 = fs[0]
        rows.append([esc(name), fmt(f1.get("reads_examined", 0)),
                     " / ".join(str(int(f.get("median_len") or 0)) for f in fs),
                     " / ".join(f'{f.get("mean_q") or 0:.1f}' for f in fs),
                     " / ".join(f'{f.get("pct_q30") or 0:.1f}%' for f in fs)]
                    + ([fmt(f1.get("n50", 0))] if qc.get("platform") == "long" else []))
    paired = qc.get("platform") == "paired"
    out.append(T(["Sample", "Reads examined", "Median length" + (" R1 / R2" if paired else ""),
                  "Mean Q" + (" R1 / R2" if paired else ""), "Bases ≥ Q30"]
                 + (["N50"] if qc.get("platform") == "long" else []), rows,
                 _rawqc_caption(qc)
                 + f" (read_qc.py" + (f"; {esc(_cite(run, 'qc/read_qc.json'))}" if _cite(run, "qc/read_qc.json")
                                             else "") + ")."))
    s = qc.get("suggestion") or {}
    used = (run.params.get("truncLen_F"), run.params.get("truncLen_R")) if run.kind == "dada2" else None
    if paired:
        pf = len(run.params.get("primer_f") or "")
        pr = len(run.params.get("primer_r") or "")
        # the figure is in raw-read positions; truncLen counts after the primer
        tr = (used[0] + pf if used and used[0] else None, used[1] + pr if used and used[1] else None) if used else None
        sg = (s.get("truncLen_F", 0) + pf if s.get("truncLen_F") else None,
              s.get("truncLen_R", 0) + pr if s.get("truncLen_R") else None) if s else None
        out.append(F(svg_qprofile(qc, tr, sg),
                     "Quality along the read, all samples pooled: the line is the median, the band the 25th to 75th "
                     "percentile. Dashed lines mark where reads were cut (red) and where the data suggest cutting "
                     "(green); both are drawn in raw-read positions, i.e. including the primer.",
                     _cite(run, "qc/read_qc.json")))
        if s and run.kind == "dada2":
            txt = (f"<p><b>truncLen.</b> The data suggest {s.get('truncLen_F')} / {s.get('truncLen_R')} bp "
                   f"after primer removal. {esc(s.get('reason', ''))}")
            if used and all(used):
                diff = abs(int(used[0]) - int(s.get('truncLen_F') or 0)) + abs(int(used[1]) - int(s.get('truncLen_R') or 0))
                txt += (f" This run used {used[0]} / {used[1]} bp"
                        + (", which matches." if diff <= 10 else ".") )
            out.append(txt + "</p>")
    elif qc.get("platform") == "long":
        out.append(F(svg_longread(qc), "Read length and mean read quality of the raw reads, all samples pooled.",
                     _cite(run, "qc/read_qc.json")))
    if qc.get("low_depth"):
        out.append(f'<p class="dim">Fewer than 1,000 reads in: {", ".join(esc(x) for x in qc["low_depth"])}.</p>')
    if (qc.get("tools") or {}).get("multiqc") and getattr(run, "pkg_files", None) is None:   # run folder only
        out.append('<p class="dim">FastQC and MultiQC reports: qc/fastqc/ (multiqc_report.html).</p>')
    return out


def _controls_section(run, an, T):
    """Negative controls and mock communities, whatever the pipeline (controls.py)."""
    try:
        import controls
    except Exception:
        return []
    negs, mocks = controls.identify(run)
    if not negs and not mocks:
        return []
    out = ['<h3 class="pbreak">Controls</h3>']
    if negs:
        ns = controls.negative_summary(an, negs, mocks)
        out.append("<p><b>Negative controls.</b> "
                   + "; ".join(f"{esc(n['sample'])}: {fmt(n['reads'])} reads"
                               + (f" ({100*n['ratio_to_median']:.1f}% of the median sample)"
                                  if n.get("ratio_to_median") is not None else "")
                               for n in ns["negatives"])
                   + ". A blank with more than a few percent of a typical sample's reads points to "
                     "contamination or cross-talk worth checking before the results are used. "
                   + (f"Genera at least twice as abundant in the blanks as in the samples — the signature of a "
                      f"reagent or handling contaminant — are marked: "
                      f"{', '.join('<i>' + esc(g) + '</i>' for g in ns['flagged'])}." if ns["flagged"] else
                      "No genus is markedly more abundant in the blanks than in the samples.")
                   + "</p>")
        rows = [[f"<i>{esc(r['genus'])}</i>" + (' <b style="color:#A83232">&#9888;</b>' if r["flag"] else ""),
                 f"{r['pct_in_blanks']:.2f}%", f"{r['mean_pct_in_samples']:.2f}%",
                 f"{r['samples_present']} / {r['n_samples']}"] for r in ns["genera"]]
        out.append(T(["Genus", "% of blank reads", "Mean % in samples", "Samples with ≥ 0.1%"], rows,
                     "The most abundant genera in the negative controls (reads pooled over blanks) and how "
                     "much of them the samples carry. &#9888; = at least 1% of blank reads and at least twice "
                     "as abundant in the blanks as in the samples. Nothing has been removed"
                     + (" (Decontam was not run)." if run.params.get("decontam") in (None, "", "none") else ".")))
    if mocks:
        std = run.params.get("mockStandard") or "zymo_d6300"
        ms = controls.mock_summary(an, mocks, std)
        if ms.get("no_composition"):
            out.append(f"<p><b>Mock community.</b> {esc(', '.join(mocks))} ({esc(ms['standard'])}): no "
                       "theoretical composition is on file for this standard, so it is shown only in the "
                       "composition figures.</p>")
        for m in ms.get("samples", []):
            verdict = ("matches the standard closely" if m["bray_curtis"] <= 0.15 else
                       "is close to the standard" if m["bray_curtis"] <= 0.30 else
                       "departs clearly from the standard")
            out.append(f"<p><b>Mock community {esc(m['sample'])}</b> ({esc(ms['standard'])}) {verdict}: "
                       f"Bray&ndash;Curtis {m['bray_curtis']:.2f} to the theoretical composition, "
                       f"{m['other_pct']:.2f}% of reads in genera that are not in the standard"
                       + (f", not found: {', '.join('<i>' + esc(g) + '</i>' for g in m['missing'])}"
                          if m["missing"] else "")
                       + ". Ratios far from 1 for one genus usually reflect primer bias or 16S copy number "
                         "rather than a pipeline error; genera that are missing or reads in unexpected "
                         "genera are the warning signs.</p>")
            rows = [[f"<i>{esc(r['genus'])}</i>", f"{r['expected']:.1f}%",
                     (f"{r['observed']:.1f}%" if not r["family_only"] else
                      f'{r["observed"]:.1f}% <span class="dim">(family only)</span>'),
                     f"{r['ratio']:.2f}" if r["ratio"] is not None else "—"] for r in m["genera"]]
            out.append(T(["Genus", "Expected", "Observed", "Observed / expected"], rows,
                         f"Mock community {esc(m['sample'])}: theoretical and observed genus composition."))
    return out


def _maaslin_section(run, T):
    """MaAsLin 3 results (r_tables/maaslin3_results.csv), or why it did not run."""
    d = run.path / "maaslin3"
    res = run._rows("r_tables/maaslin3_results.csv")
    skipped = (d / "SKIPPED.txt").read_text().strip() if (d / "SKIPPED.txt").exists() else ""
    if not res and not skipped:
        return []
    out = ['<h3 class="pbreak">Differential abundance (MaAsLin 3)</h3>']
    if not res:
        out.append(f"<p>MaAsLin 3 was not run: {esc(skipped)}.</p>")
        return out
    try:
        design = json.loads((d / "design.json").read_text())
    except Exception:
        design = {}
    fixed = design.get("fixed") or sorted({r.get("metadata", "") for r in res})
    rank = design.get("rank") or (res[0].get("rank") if res else "Genus")
    thr = float(design.get("max_significance") or 0.1)

    def q(r):
        try:
            return float(r.get("qval_individual") or "nan")
        except ValueError:
            return float("nan")
    sig = sorted([r for r in res if q(r) == q(r) and q(r) < thr], key=q)
    out.append(
        f"<p>Each {esc(rank.lower())} was tested for an association with "
        + ", ".join(f"<b>{esc(v)}</b>" for v in fixed)
        + (f" (random effect: {', '.join(esc(x) for x in design['random'])})" if design.get("random") else "")
        + f", all variables in one model, in two ways: <i>abundance</i> (how much of it, where it is "
          f"present; log of total-sum-scaled counts) and <i>prevalence</i> (whether it is present at all). "
          f"{len(res)} tests on {design.get('n_features', '?')} features and {design.get('n_samples', '?')} "
          f"samples; p-values were corrected by Benjamini&ndash;Hochberg. "
        + (f"{len(sig)} association{'s' if len(sig) != 1 else ''} reached q &lt; {thr:g}." if sig else
           f"No association reached q &lt; {thr:g}.")
        + (f" Not tested (no level with two or more samples): {', '.join(esc(x) for x in design['not_tested'])}."
           if design.get("not_tested") else "") + "</p>")
    if sig:
        rows = []
        for r in sig[:40]:
            try:
                coef = float(r.get("coef") or 0)
            except ValueError:
                coef = 0.0
            rows.append([f"<i>{esc(r.get('feature', ''))}</i>", esc(r.get("metadata", "")),
                         esc(r.get("value", "")), esc(r.get("model", "")),
                         f'<span style="color:{"#2E6B46" if coef > 0 else "#A83232"}">{coef:+.2f}</span>',
                         f"{q(r):.3g}"])
        out.append(T(["Feature", "Variable", "Level", "Model", "Coefficient", "q"], rows,
                     f"Associations with q &lt; {thr:g}, strongest first"
                     + (f" (the first 40 of {len(sig)})" if len(sig) > 40 else "")
                     + ". A positive coefficient means more (abundance) or more often present (prevalence) "
                       "at the level shown than at the reference level, or with a higher value of a "
                       "numeric variable. All results: r_tables/maaslin3_results.csv; figures in maaslin3/."))
    out.append('<div class="note"><p style="margin:0"><b>Read with the design in mind.</b> MaAsLin 3 '
               "adjusts every association for the other variables listed. With few samples per group the "
               "test has little power, so the absence of significant features is not evidence that groups "
               "are the same.</p></div>")
    return out


def build_profile_html(run, title=None, app="NextGen-Amplicon", company="",
                       subtitle=None, top_asv=30, archaea=True, files=None, ids=None) -> str:
    """Full result report for one run, laid out as Method / Result / References with a
    figure for every analysis: read tracking, composition at every rank, a clustered
    genus heatmap, archaea, alpha diversity and rarefaction, and beta diversity.

    files: for a Download results package — {file of the run folder: [its files in the package]}, so
    the 'Folder:' lines under the figures name what the customer has (5.Taxonomy_Identification/csv/..., not r_tables/...).
    ids: {sequence / feature id of the run: its ID in that package's count table (SV001 ...)}, so the
    report's table of the most abundant features uses the IDs the customer can look up."""
    run.pkg_files = files                      # the package's file names, for the text that cites a file (_cite)
    run.pkg_ids = ids
    if run.kind == "shotgun":                  # its own layout: no amplicon steps or markers
        import shotgun_report
        return shotgun_report.build(run, title=title, app=app, company=company, subtitle=subtitle)
    an = Analysis(run)
    S = an.samples
    unit = run.feat
    marker = run.summary.get("marker") or run.params.get("marker") or "16S"
    if run.kind == "pacbio":
        marker = "16S"
    elif run.kind == "cox1":
        marker = "COI"
    elif run.mito_gene == "CytB":
        marker = "CytB"
    elif getattr(run, "is_pr2", False) and "18S" in str(marker).upper():
        marker = "18S"                 # the job's marker key is "18S-nema" whatever the reference
    region = run.summary.get("region") or ("V1–V9" if run.emu else
                                          (run.params.get("pb_region") or "V1–V9") if run.kind == "pacbio"
                                          else "V3–V4"
                                          if (run.params.get("primer_f") or "").startswith("CCTACGG")
                                          else "")
    rep_title = (f"16S {region} Amplicon Sequencing Report".replace("  ", " ")
                 if "16S" in marker.upper() else
                 ("18S rRNA V4 Eukaryote Amplicon Sequencing Report" if getattr(run, "is_pr2", False)
                  else "18S rRNA Nematode Amplicon Sequencing Report")
                 if "18S" in marker.upper() or "NEMA" in marker.upper() else f"{marker} Amplicon Sequencing Report")
    project = title if title and title.strip() and title.strip() != rep_title else None
    title = rep_title          # the report type heads the page; a user title becomes the project line
    tech = {
        "emu": "Oxford Nanopore full-length 16S, classified with Emu",
        "vsearch": ("Illumina paired-end" if run.paired else "single-end reads")
                   + ", VSEARCH operational taxonomic units",
        "its": "Illumina paired-end, DADA2 amplicon sequence variants, "
               + ("named by BLASTn against NCBI nt (LCA)" if ncbi_lca(run) else "UNITE reference"),
        "cox1": "Illumina paired-end, DADA2 amplicon sequence variants, MIDORI2 reference",
        "pacbio": "PacBio HiFi full-length 16S, DADA2 amplicon sequence variants",
    }.get(run.kind, "Illumina paired-end, DADA2 amplicon sequence variants")
    if run.blast:
        tech += ", named by BLASTn (MIDORI2 + MitoFish)"
    subtitle = subtitle or (f"{run.label} · {len(S)} sample{'s' if len(S) != 1 else ''} · " + tech)
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
        if files is not None:
            have = list(dict.fromkeys(x for n in names for x in (files.get(n) or [])))
        else:
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
    fin_word = {"emu": "reads classified", "vsearch": "reads mapped to OTUs",
                "cox1": "reads after curation"}.get(run.kind, "non-chimeric reads")
    kpis = [(fmt(inp) if inp else "—", f"{run.reads_word} sequenced"),
            (fmt(fin), fin_word),
            (str(run.n_asvs), f"{unit} detected"),
            (f"{run.rank_pct('Genus') or 0}%", "reads with a genus")]
    out.append('<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>'
        for v, l in kpis) + "</div>")
    if run.only:
        gl_ = an.group_list()
        out.append('<div class="note acc"><b>This report covers ' + f"{len(S)} of the run&rsquo;s "
                   f"{len(run.all_samples)} samples</b>"
                   + (f": the group{'s' if len(gl_) != 1 else ''} "
                      + ", ".join(f"<b>{esc(g)}</b> ({len(v)})" for g, v in gl_)
                      + f" of the grouping &ldquo;{esc(an.group_col)}&rdquo;" if gl_ else "")
                   + ". Every number, table and figure below is computed from these samples only."
                   + (" The groups are compared in the chapter <b>Analysis by group</b>, after Beta diversity."
                      if len(gl_) >= 2 else "") + "</div>")
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
    cfg = ([("Marker", marker)] + run.emu_config()) if run.emu else run.config(marker)
    out.append(T(["Parameter", "Value"], [[k, esc(v if v is not None else "—")] for k, v in cfg],
                 "Run configuration."))
    prov = _provenance_rows(run)
    if prov:
        out.append(T(["Software / reference", "Version", "Checksum (MD5)"], prov,
                     "Software versions and reference files this run used, recorded when it "
                     "finished (provenance.json in the result folder)."
                     + (" &#10003; marks a reference whose checksum matches the file its publisher "
                        "released." if any("&#10003;" in r[2] for r in prov) else "")))

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
        steps = ([(l.lower(), k) for l, k in run.stages[1:]] if run.kind not in ("dada2", "its") else
                 [("filtering", "filtered"), ("denoising", "denoisedF"),
                  ("merging", "merged"), ("chimera removal", "nonchim")])
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
        flow = {
            "vsearch": f"All {len(S)} samples were processed through primer trimming"
                       + (", read-pair merging" if run.paired else "")
                       + ", quality filtering, dereplication, chimera removal and UNOISE denoising "
                         "in VSEARCH, and every read was then mapped back to the OTUs. ",
            "cox1": f"All {len(S)} paired-end samples were processed through primer trimming, quality "
                    f"filtering, denoising, read-pair merging and chimera removal in DADA2, followed by "
                    f"a length filter, removal of sequences with stop codons (likely pseudogenes) and "
                    f"LULU curation. ",
            "pacbio": f"All {len(S)} samples were oriented and primer-trimmed with cutadapt, then "
                      f"quality-filtered, denoised with the PacBio error model and checked for chimeras "
                      f"in DADA2. ",
            "its": f"All {len(S)} paired-end samples were processed through primer trimming (both "
                   f"orientations, no length truncation), quality filtering, denoising, read-pair "
                   f"merging and chimera removal in DADA2. ",
        }.get(run.kind, f"All {len(S)} paired-end samples were processed through primer trimming, "
                        f"quality filtering, denoising, read-pair merging and chimera removal in DADA2. ")
        out.append("<p><b>Overall pipeline outcome.</b> "
                   + (f"All {len(S)} samples were processed through "
                      + ", ".join(x for x, k in (("quality filtering", "qc"), ("primer trimming", None),
                                                 ("chimera removal", "chimera"))
                                  if not (k and _skipped(run)[k]))
                      + ", and the remaining reads were classified with Emu. "
                      if run.emu else flow)
                   + f"From {fmt(inp)} input {run.reads_word}, {fmt(fin)} "
                   f"({100*fin/inp if inp else 0:.1f}%) were retained "
                   + (f"and assigned to {run.n_asvs} taxa. " if run.emu else
                      f"and resolved into {run.n_asvs} {run.feat}. ")
                   + (f"{fmt(kreads.get('Bacteria', 0))} reads were classified as Bacteria and "
                      f"{fmt(kreads.get('Archaea', 0))} as Archaea. " if kreads and run.bacterial else "")
                   + (f"{fmt(run.host + run.euk)} reads ({100*(run.host+run.euk)/run.total:.1f}%) were host "
                      f"organelle or eukaryotic and are excluded from the community figures. "
                      if run.total and (run.host + run.euk) > 0 else
                      "No mitochondrial, chloroplast or eukaryotic reads were detected. "
                      if run.bacterial else "")
                   + "</p>")
        qi = (f"<p><b>QC interpretation.</b> Input depth ranged "
              f"{_range_txt([(n, i) for n, i, _f, _p in per])} per sample. The final retained reads "
              f"ranged {_range_txt([(n, f) for n, _i, f, _p in per])}, a retention of "
              f"{min(x[3] for x in per):.1f}&ndash;{max(x[3] for x in per):.1f}%. ")
        cap = {k: v for k, v in (run.denoise_cap or {}).items() if k in {x[0] for x in per}}
        if worst and not (cap and worst[0].lower() == "denoising"):
            qi += (f"The largest loss across the run was at the {esc(worst[0].lower())} step "
                   f"({worst[1]:.1f}% of the reads entering that step). ")
        low = [x for x in per if x[3] < 50]
        if cap:
            n_cap = sorted({b for _a, b in cap.values()})
            qi += (f"Before denoising, {len(cap)} sample{'s' if len(cap) != 1 else ''} with more than "
                   f"{fmt(n_cap[0])} filtered reads {'were' if len(cap) != 1 else 'was'} randomly subsampled to "
                   f"{fmt(n_cap[0])} reads ({', '.join(esc(k) for k in cap)}), which bounds the denoising run "
                   f"time and leaves the composition unchanged. The drop at the denoising step and the "
                   f"retention below are mostly this subsampling, not reads lost to quality. ")
            # judge loss on the reads that actually entered denoising
            low = []
            for t in run.tracking:
                nm = str(next(iter(t.values())))
                ent = cap[nm][1] if nm in cap else float(t.get("filtered") or 0)
                f0 = float(t.get(run.retained_key) or 0)
                if ent and 100 * f0 / ent < 50:
                    low.append((nm, ent, f0, 100 * f0 / ent))
        if low:
            byname = {str(next(iter(t.values()))): t for t in run.tracking}
            def _worst_step(nm):
                t, prev, best = byname.get(nm) or {}, "input", None
                for lbl, k in steps:
                    if k not in t:
                        continue
                    a = cap[nm][1] if (k == "denoisedF" and nm in cap) else float(t.get(prev) or 0)
                    b = float(t.get(k) or 0)
                    if a and (best is None or (a - b) / a > best[1]):
                        best = (lbl, (a - b) / a)
                    prev = k
                return best
            ws = {x[0]: _worst_step(x[0]) for x in low}
            chim = [n for n, w in ws.items() if w and w[0] == "chimera removal"]
            rest = [n for n in ws if n not in chim]
            if chim:
                qi += (f"{', '.join(esc(n) for n in chim)} lost most at chimera removal ("
                       + ", ".join(f"{100*ws[n][1]:.0f}%" for n in chim)
                       + " of merged reads): a high chimera share comes from PCR (many cycles or little "
                         "template), not from sequencing quality. ")
            if rest:
                qi += (f"{', '.join(esc(n) for n in rest)} retained less than half of "
                       f"{'its' if len(rest) == 1 else 'their'} "
                       + ("reads that entered denoising" if cap else "input")
                       + ", which deserves a look at the raw read quality before the results are interpreted. ")
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

        if run.kind not in ("dada2", "its"):
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
                        " (raw, filtered, denoisedF, denoisedR, merged, non-chimeric)."
                        if run.kind in ("dada2", "its") else
                        " (" + ", ".join(l.lower() for l, _k in run.stages) + ").")))
        out.append(F(svg_steps(run, an), "Read count distribution per pipeline step. Each line is one "
                     "sample" + ("; boxes show the spread across samples." if len(S) >= 4 else "."),
                     folder("read_tracking.csv")))
        out.append(F(svg_readflow(run), "Reads surviving each stage across the whole run, as a share "
                     "of input; the right-hand figure is what that stage cost.", None))

    out += _readqc_section(run, F, T)

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

    out.append(T(["Sample", "Reads", f"n_{run.feat1}"],
                 [[esc(s), fmt(sum(an.mat[q][s] for q in run.abund)), fmt(an.n_features(s))] for s in S],
                 f"Reads and number of {unit} per sample after processing."))
    out.append(F(svg_rare(an), f"Rarefaction curves: expected number of {unit} against sequencing depth "
                 "for each sample.", folder("r_tables/rarefaction.csv")))

    qc = (run.qc or {}).get("offtarget") or {}
    if qc and qc.get("reads"):
        out.append("<h3>Off-target reads</h3>")
        _org = float(qc.get("organelle_reads", 0) or 0)
        _dom = max(0.0, float(qc.get("reads", 0) or 0) - _org)
        _tot = (float(qc.get("reads", 0) or 0) / (float(qc.get("pct", 0) or 0) / 100)) if qc.get("pct") else 0
        _dpct = f"{100 * _dom / _tot:.2f}%" if _tot else "&mdash;"
        _bact = run.bacterial
        out.append(T(["Category", "Reads", "Share"], [
            ["Host organelle (mitochondria / chloroplast)", fmt(_org), f'{qc.get("organelle_pct", 0)}%'],
            ["Outside the target domain (" + ("eukaryotes" if _bact else "bacteria and archaea") + ")",
             fmt(_dom), _dpct],
            ["Total off-target", fmt(qc.get("reads", 0)), f'{qc.get("pct", 0)}%'],
            ["Not assigned at domain level (not counted above)", fmt(qc.get("unassigned_reads", 0)),
             f'{qc.get("unassigned_pct", 0)}%']],
            "Off-target reads. &ldquo;Universal&rdquo; 16S primers co-amplify host mitochondrial DNA."
            if _bact else
            "Off-target reads: DNA outside the domain this marker targets that the primers also amplified, "
            "and organelle sequences."))

    if run.blast:
        out += _blast_section(run, S, T)
    out += _ncbi_lca_section(run, S, T)
    if getattr(run, "ncbi_un_info", None) is not None:
        out += _ncbi_un_section(run, S, T)

    # ── composition ──────────────────────────────────────────────────────────
    out.append('<h3 class="pbreak">Taxonomic composition of samples</h3>')
    ph_all = an.rank_pct("Phylum")
    mean_ph = sorted(((t, sum(v.values()) / len(S)) for t, v in ph_all.items() if t != "Unclassified"),
                     key=lambda kv: -kv[1])[:4]
    txt = ("<p>At the phylum level the community was dominated by "
           + ", ".join(f"{esc(t)} ({v:.1f}% on average)" for t, v in mean_ph) + ". ")
    if run.blast:
        txt += ("Names here are the BLAST names; a sequence named only to phylum mostly had no reference within "
                "the class threshold (its phylum is a distant best hit and places it only approximately) or had "
                "conflicting best hits. ")
    gphy: dict = {}
    for r_ in run.tax.values():
        g_, ph_ = run._clean_genus(r_.get("Genus")), (r_.get("Phylum") or "").strip()
        if g_ and ph_:
            gphy.setdefault(g_, Counter())[ph_] += 1
    for s in (S if len(S) <= 8 else []):          # beyond 8 samples the table below says it
        dg = an.dominant("Genus", s, 3)
        dp = [(gphy[dg[0][0]].most_common(1)[0][0], 0)] if dg and dg[0][0] in gphy else []
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

    many = len(S) > MANY
    gby = an.groups if len({an.groups.get(x) for x in S} - {None}) >= 2 else {}
    parts = (_group_chunks(S, gby) if gby else _chunks(S)) if many else [S]
    for rank in RANKS:
        if rank not in run.rank_cols:
            continue
        keep, _t = an.top(rank, 30)
        if not keep:
            continue
        cap = (f"Relative abundance of the {len(keep)} most abundant "
               f"{PLURAL.get(rank, rank.lower()) if len(keep) != 1 else rank.lower()} per sample"
               + ". Percentages are of each sample&rsquo;s community reads"
               + (" with host organelle reads set aside" if run.host + run.euk > 0 else "") + ".")
        if not many:
            out.append(f'<h3 class="pbreak">Top 30 {rank} (Relative Abundance)</h3>')
            out.append(F(svg_stacked(an, rank, 30), cap, folder(f"r_tables/taxonomy_{rank.lower()}.csv")))
            continue
        for i, part in enumerate(parts):          # landscape, at most CHUNK samples a page
            out.append('<section class="wide">'
                       f'<h3 class="pbreak">Top 30 {rank} (Relative Abundance)'
                       + (f" &mdash; part {i + 1} of {len(parts)}" if len(parts) > 1 else "") + "</h3>"
                       + F(svg_stacked(an, rank, 30, samples=part, wide=True, groups=gby or None),
                           cap + f" Samples {esc(part[0])} to {esc(part[-1])} ({len(part)} of {len(S)}"
                           + (", ordered by group" if gby else "") + "); "
                           "the taxa and colours are the same on every part.",
                           folder(f"r_tables/taxonomy_{rank.lower()}.csv") if i == len(parts) - 1 else None)
                       + "</section>")

    if not many:
        out.append('<h3 class="pbreak">Heatmap visualization of dominant genera</h3>')
        out.append("<p>The genus-level heatmap summarises the 30 dominant genera across samples. Genera "
                   "(rows) are clustered by the similarity of their abundance profiles"
                   + (" and samples (columns) by Bray&ndash;Curtis dissimilarity" if len(S) >= 3 else "")
                   + "; colour is relative abundance on a log scale, and each cell gives the percentage. "
                   "Blocks of dark cells confined to one sample mark genera characteristic of it.</p>")
        out.append(F(svg_heatmap(an, "Genus", 30), "Top 30 genera &mdash; relative abundance heatmap with "
                     "clustering dendrogram.", folder("r_tables/taxonomy_heatmap_genus.csv")))
    else:
        for i, part in enumerate(parts):
            out.append('<section class="wide">'
                       '<h3 class="pbreak">Heatmap visualization of dominant genera'
                       + (f" &mdash; part {i + 1} of {len(parts)}" if len(parts) > 1 else "") + "</h3>"
                       + ("<p>The 30 dominant genera of the run (rows, clustered by the similarity of their "
                          "abundance profiles) across the samples (columns, in the order of the stacked bars, "
                          f"{len(parts)} parts); colour is relative abundance on a log scale, the same scale on "
                          "every part. Dark cells confined to a few samples mark genera characteristic of them; "
                          "the samples&rsquo; similarity is shown under Beta diversity.</p>" if i == 0 else "")
                       + F(svg_heatmap(an, "Genus", 30, samples=part, wide=True),
                           f"Top 30 genera &mdash; relative abundance, samples {esc(part[0])} to {esc(part[-1])}.",
                           folder("r_tables/taxonomy_heatmap_genus.csv") if i == len(parts) - 1 else None)
                       + "</section>")

    # Archaea: only when some were found (a page of zeros says nothing), and never when the
    # report was asked without it (a project about Bacteria only)
    show_arch = run.bacterial and archaea and n_arch > 0
    if show_arch:
        out.append("<h3>Archaeal community</h3>")
    if not show_arch:
        pass                                   # 18S / ITS / COI, none found, or left out
    elif n_arch > 0:
        out.append(f"<p>Archaeal reads were {'rare' if n_arch / max(fin, 1) < 0.01 else 'present'} in this "
                   f"dataset: {fmt(n_arch)} reads in total, "
                   f"{100*n_arch/max(sum(an.comm_tot.values()),1):.3f}% of the community reads. "
                   + ("The primers are designed for Bacteria and amplify most Archaea poorly, so this is a "
                      "minimum rather than an estimate of archaeal load." if not run.emu else "") + "</p>")
    else:
        out.append("<p>No archaeal reads were detected in any sample.</p>")
    if show_arch:
        rows_a = [x for x in arch_rows if x[1] > 0]
        out.append(T(["Sample", "Archaeal reads", "% of sample", "Dominant taxon"],
                     [[esc(s), fmt(r), f"{p_:.3f}%", f"<i>{esc(t)}</i>" if t != "—" else t]
                      for s, r, p_, t in rows_a],
                     "Archaeal reads per sample" + (f" ({len(rows_a)} of {len(arch_rows)} samples hold any; the "
                                                    "others have none)" if len(rows_a) < len(arch_rows) else "")
                     + "."))

    # top features, species. (The 'Reagent-associated genera' table was dropped in v2.42.0 at the BI
    # team's request: a list of kit genera is no measurement; a sequenced blank (see Controls) is.)
    out.append(f"<h3>Most abundant {run.unit}</h3>")
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
        name = an.label(q, "Species") if run.kind in ("emu", "vsearch", "its", "cox1") else "Unclassified"
        if name == "Unclassified":
            name = a["genus"] if a["genus"] != "\u2014" else ""
        nm = ('<span class="dim">host organelle</span>' if a["host"]
              else f'<i>{esc(name)}</i>' if name else f'<span class="dim">{esc(a["deepest"])}</span>')
        rows.append([f'<span class="mono">{a["id"]}</span>'] + ([] if run.kind in ("emu", "vsearch") else [f'{a["len"] or "—"}'])
                    + [fmt(a["reads"]), f'{a["pct"]:.2f}%']
                    + ([f'{100*an.mat[q][s_]/an.comm_tot[s_]:.2f}' if an.comm_tot[s_] else "0"
                        for s_ in S] if per else [])
                    + [f'<div style="text-align:left">{nm}</div>',
                       f'<div style="text-align:left;font-size:7pt" class="dim">{_lin(q)}</div>'])
    out.append(T((["Taxon"] if run.emu else ["OTU"] if run.kind == "vsearch" else ["ASV", "bp"]) + ["Reads", "All"]
                 + ([f"{esc(h)} %" for h in _col_heads(S)] if per else [])
                 + (["Species / genus"] if run.kind in ("emu", "vsearch", "its", "cox1") else ["Genus"])
                 + ["Phylum &rsaquo; Class &rsaquo; Order &rsaquo; Family"],
                 rows,
                 f"The {top_asv} most abundant of {run.n_asvs} {unit}: reads and share over the whole run"
                 + (", and the share of each sample" if per else "")
                 + ((f". Identifiers are those of {esc(_cite(run, 'asv_table.csv'))}"
                     if getattr(run, "pkg_ids", None) and _cite(run, "asv_table.csv") else
                     ". Identifiers are positional in abundance order")
                    + ((f"; sequences are in {esc(_cite(run, 'asvs.fasta'))}." if _cite(run, "asvs.fasta") else ".")
                       if run.kind not in ("emu", "vsearch") else
                       (f"; centroid sequences are in {esc(_cite(run, 'otus.fasta'))}." if _cite(run, "otus.fasta") else ".")
                       if run.kind == "vsearch" else "."))))

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
                        "A species is named only when every NCBI hit that passes the identity and coverage "
                        "names it; GenBank names are only as good as the records deposited under them, and a "
                        "species missing from GenBank is named at a higher rank or not at all." if ncbi_lca(run) else
                        "Species come from the UNITE species hypothesis the ASV was classified into; ITS "
                        "separates most fungal species, but names in UNITE are only as good as the "
                        "sequences deposited under them." if run.kind == "its" else
                        "Species come from the MIDORI2 reference the ASV was classified into. COI separates "
                        "most animal species, but a name is only as reliable as the reference coverage of "
                        "that group; an unexpected species should be checked against BOLD." if run.kind == "cox1" else
                        "Species are assigned by exact match of the full-length sequence; the whole gene "
                        "separates many but not all species." if run.kind == "pacbio" else
                        "Species come from the closest reference sequence above the identity threshold; "
                        "treat them as the nearest named relative." if run.kind == "vsearch" else
                        (f"Species are named only when the best BLAST hit reaches "
                         f"{(run.blast_info.get('thresholds') or {}).get('Species', 98.5)}% identity and no other "
                         f"species lies within 0.5%; otherwise the name stops at the deepest rank whose threshold "
                         f"is met. A species missing here is usually absent from the reference, not from the "
                         f"sample.") if run.blast else
                        "Species are assigned only by exact match of the whole fragment; for a mitochondrial "
                        "marker the reference covers only part of the local fauna, so a missing species is "
                        "usually missing from the reference." if run.mito_gene else
                        "PR2 species are the classifier&rsquo;s best match at the bootstrap threshold. The ~380 bp "
                        "18S V4 fragment separates many genera but few species, so read species names as tentative "
                        "and the genus as the reliable rank." if getattr(run, "is_pr2", False) else
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
    if not many:
        out.append(F(svg_alpha(an, alpha), "Alpha diversity indices per sample.",
                     folder("r_tables/alpha_diversity.csv", "alpha_diversity.csv")))
    else:
        for k, idx_ in enumerate((["Observed", "Chao1"], ["Shannon", "Simpson"])):
            out.append('<section class="wide">'
                       + f'<h3 class="pbreak">Alpha diversity &mdash; {" and ".join(idx_)}</h3>'
                       + F(svg_alpha(an, alpha, idx_, wide=True, samples=[x for pt in parts for x in pt]),
                           f"{' and '.join(idx_)} per sample" + (", ordered by group" if gby else "") + ".",
                           folder("r_tables/alpha_diversity.csv", "alpha_diversity.csv") if k else None)
                       + "</section>")

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
        pg = an.permanova_groups() if an.groups else None
        pm = (pg["R2"], pg["p"]) if pg else (None if an.groups else an.permanova())
        if pm:
            bt += (f"PERMANOVA on group gave R&sup2; = {pm[0]:.3f}, p = {pm[1]:.3g}"
                   + (" &mdash; a significant group effect (see Analysis by group). " if pm[1] < 0.05 else
                      " &mdash; no significant group effect (see Analysis by group). "))
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
        pc = svg_pcoa(an, wide=many)
        if pc:
            fig_ = F(pc, f"Principal coordinates analysis (PCoA) of Bray&ndash;Curtis dissimilarity. "
                     f"Axis 1 explains {an.pcoa_pct[0]:.1f}% and axis 2 {an.pcoa_pct[1]:.1f}% of the "
                     f"variation." + (" Sample names are printed where they do not overlap." if many else ""),
                     folder("r_tables/pca_scores.csv", "beta_stats.txt"))
            out.append('<section class="wide"><h3 class="pbreak">Beta diversity &mdash; PCoA</h3>' + fig_
                       + "</section>" if many else fig_)
        if len(S) >= 3:
          tree_ = F(svg_upgma(an), "UPGMA clustering tree based on Bray&ndash;Curtis dissimilarity.",
                    folder("beta_stats.txt"))
          out.append(f'<div class="pbreak">{tree_}</div>' if many else tree_)
        if len(S) <= 40:
            out.append(F(svg_bc_heat(an), "Bray&ndash;Curtis dissimilarity between samples; darker cells and "
                         "higher values are more dissimilar.", folder("r_tables/beta_heatmap.csv", "beta_stats.txt")))
        else:
            have = [n_ for n_ in ("r_tables/beta_heatmap.csv", "beta_stats.txt") if (run.path / n_).exists()]
            out.append(f"<p>With {len(S)} samples the {len(S)} &times; {len(S)} Bray&ndash;Curtis matrix is too "
                       "large to print legibly; the tree above orders the samples by it"
                       + (f", and the matrix itself is in {', '.join(have)}" if have else "") + ".</p>")

    # ── group view ───────────────────────────────────────────────────────────
    out += _group_section(run, an, alpha, S, F, T)

    out += _controls_section(run, an, T)
    out += _maaslin_section(run, T)

    if pic:
        out += _picrust_result(pic, F, T, run=run)

    # ── notes / references ───────────────────────────────────────────────────
    out.append('<h2 class="pbreak">Notes on reading this report</h2>')
    notes = []
    if run.total and run.host + run.euk > 0:
        notes.append(f"{100*(run.host+run.euk)/run.total:.1f}% of analysed reads are host organelle or "
                     "other non-bacterial material and are excluded from the composition figures.")
    if run.unresolved:
        notes.append(f"{fmt(sum(run.unresolved.values()))} reads carry a genus label the reference marks "
                     "as unresolved (&ldquo;--other&rdquo;); they are counted under the parent genus here.")
    if run.genera.get(HAEMOPLASMA_LABEL):
        notes.append("Haemotropic mycoplasmas are reported as one group (Eperythrozoon and "
                     "Haemobartonella combined).")
    if getattr(run, "rank_shift_fixed", False):
        notes.append("This run was made before v2.17, when COX1 rank names were one level off (MIDORI2 "
                     "starts at phylum). The ranks are shown corrected here; the run's own "
                     "taxonomy_table.csv still carries the old column names.")
    notes.append({
        "its": ("Relative abundance from ITS amplicons in environmental DNA reflects rDNA copies and the DNA "
                "each organism shed, not colonies or cover: copy number differs between taxa, and ITS length "
                "differences favour shorter amplicons." if ncbi_lca(run) else
                "Relative abundance from ITS amplicons reflects rDNA copies, not cells: copy number varies "
                "more than tenfold between fungi, and ITS length differences favour shorter amplicons."),
        "cox1": "Relative abundance from COI amplicons reflects mitochondrial copies and primer fit, not "
                "biomass or individuals; COI primers are degenerate and amplify some groups better than "
                "others.",
        "vsearch": "Relative abundance from ITS amplicons reflects rDNA copies, not cells: copy number varies "
                   "more than tenfold between fungi, and ITS length differences favour shorter amplicons."
                   if str(run.params.get("marker") or run.summary.get("marker") or "").upper().startswith("ITS")
                   else "Relative abundance from rRNA amplicons reflects gene copies, not cells: taxa with "
                        "many rRNA operons are over-represented.",
    }.get(run.kind, "Relative abundance in environmental DNA reflects the DNA each organism shed and how "
                    "well the primers fit it, not the number or biomass of organisms."
                    if run.mito_gene else
                    "Relative abundance from rRNA amplicons reflects gene copies, not cells: taxa with "
                    "many rRNA operons are over-represented."))
    if files is not None:
        _trk = _cite(run, "read_tracking.csv")
        notes.append("Full tables are in this result's folders: "
                     + (f"{esc(_trk.split('/')[0])} (read tracking), " if _trk else "")
                     + (f"{esc(_cite(run, 'asv_table.csv'))} (the count table with every taxon's lineage), "
                        if _cite(run, "asv_table.csv") else "")
                     + "5.Taxonomy_Identification/csv (every rank, a column per sample) and 6.Diversity "
                     + ("(alpha and beta)." if _cite(run, "r_tables/beta_heatmap.csv") else "(alpha)."))
    else:
        _tx = "taxonomy.csv" if run.kind in ("emu", "vsearch") else "taxonomy_table.csv"
        _have = [f"{x} ({what})" for x, what in ((_tx, "the names"), ("asv_table.csv", f"per-{run.feat1} counts"),
                                                  ("r_tables", "per-rank tables, diversity"),
                                                  ("read_tracking.csv", "reads per step"))
                 if (Path(run.path) / x).exists()]
        if _have:
            notes.append("Full tables are in the result folder: " + ", ".join(_have)
                         + ("; provenance.json records the software and reference versions."
                            if run.prov and (Path(run.path) / "provenance.json").exists() else "."))
    _pc = _cite(run, "PICRUSt2") if pic else None
    if pic and _pc:
        notes.append(f"Predicted functions (PICRUSt2) are in {esc(_pc)}/: KO_metagenome_out and "
                     "EC_metagenome_out (gene families per sample), pathways_out (MetaCyc pathways); the "
                     "*_descrip.tsv.gz files carry descriptions.")
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")
    out.append("<h2>References</h2>")
    out.append('<ul class="refs">' + "".join(
        f"<li>{esc(a)} <i>{esc(t)}</i> {esc(j)}</li>"
        for a, t, j in ([REFERENCES[i] for i in _ref_indices(run)]
                        + (REF_PICRUST if pic else []))) + "</ul>")

    css = _css(app, company, f"{run.label} · {marker}") + _CSS_FULL.replace("@@TITLE@@", esc(title))
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            + "".join(out) + "</body></html>")


_CSS_FULL = """
@page wide { size: A4 landscape; margin: 14mm 14mm 16mm 14mm; }
section.wide { page: wide; }
section.wide h3 { margin-top: 0; }
section.wide figure { margin: 1mm 0 0; }
section.wide p { margin-bottom: 1.5mm; }
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
