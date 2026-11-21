"""Multi-marker report: one set of samples sequenced with several markers (e.g. 12S fish,
18S eukaryotes and COI animals of the same eDNA samples), each analysed as its own job.

It answers two questions, side by side for every marker:
  1. Overview per marker — what each marker targets, the reference it was named against,
     its samples, reads and how far its names reach.
  2. What was found at each sampling site — the taxa of each marker, site by site, with the
     replicates they were found in and the controls that also carried them.

Samples are matched between the jobs by their name without the marker part each job adds
(KK-NE-1-Riaz12S, KK-NE-1_TAR and KK-NE-1_COI-XT are all KK-NE-1). The site is the name
without its replicate number (KK-NE-1 -> KK-NE); a name ending in CT / NC / NTC / NEG / BLANK
is that site's control.
"""
from __future__ import annotations

import re
from collections import OrderedDict

import report_builder as rb

esc, fmt = rb.esc, rb.fmt
SEP = "-_. "
CONTROL_RE = re.compile(r"^(.*?)[-_. ]+(CT|NC|NTC|NEG|NEGATIVE|BLANK|EB|CTRL|CONTROL)[-_. ]?\d*$", re.I)
REP_RE = re.compile(r"^(.*?)[-_. ]+(?:R|REP)?(\d{1,3})$", re.I)
MAX_ROWS = 80          # taxa listed per marker in the site table
MAX_NAMES = 12         # names listed per marker in a site summary


# ── markers ───────────────────────────────────────────────────────────────────
def marker_of(run) -> str:
    """Short marker name of a run: 12S, CytB, 18S, COI, ITS, 16S."""
    mk = str(run.summary.get("marker") or run.params.get("marker") or "").upper()
    if run.kind == "cox1" or "COI" in mk or "COX" in mk:
        return "COI"
    if run.mito_gene:
        return run.mito_gene
    if "18S" in mk or "NEMA" in mk or "EUK" in mk:
        return "18S"
    if run.kind == "its" or mk.startswith("ITS"):
        return "ITS"
    if run.kind in ("emu", "pacbio") or "16S" in mk:
        return "16S"
    return mk or run.kind.upper()


def target_of(run, m: str) -> str:
    if m == "12S" or m == "CytB":
        lim = str((run.blast_info or {}).get("limit") or run.params.get("blast_limit") or "fish")
        return {"vertebrates": "vertebrates", "all": "vertebrates (no BLAST limit)"}.get(lim, "fish")
    if m == "18S":
        return "eukaryotes (protists, algae, fungi, animals)" if getattr(run, "is_pr2", False) else "nematodes"
    if m == "COI":
        return "animals (metazoans) and other eukaryotes"
    if m == "ITS":
        return "fungi"
    if m == "16S":
        return "bacteria and archaea"
    return "—"


def naming_of(run) -> str:
    if run.blast:
        return "BLASTn, MIDORI2 + MitoFish"
    if str(run.summary.get("tax_method") or "").lower() == "vsearch":
        return "VSEARCH best hit, identity cut-offs"
    return {"vsearch": "VSEARCH best hit", "emu": "Emu"}.get(run.kind, "DADA2 classifier")


def ref_of(run) -> str:
    """The reference's label, without the file name db_name adds in brackets."""
    n = str(run.db_name)
    m = re.match(r"^(.+?) \((.+)\)$", n)
    return m.group(1) if m else n


def _group_ranks(run, m: str):
    """Ranks tried, in order, for the broad organism group of a feature."""
    if m in ("12S", "CytB"):
        return ("Class", "Phylum")
    if getattr(run, "is_pr2", False):
        return ("Subdivision", "Phylum", "Supergroup")
    return ("Phylum", "Kingdom", "Domain")


def _val(v) -> str:
    v = str(v or "").strip()
    if not v or rb.BLANK_RE.match(v) or v.upper() == "NA":
        return ""
    return rb.UNRESOLVED_RE.sub("", v)


def _feature_name(run, q, m):
    """(group, rank, name) of one feature: the deepest name the run gives it."""
    if run.blast:
        r = run.blast.get(q) or {}
    else:
        r = run.tax.get(q) or {}
    group = next((_val(r.get(k)) for k in _group_ranks(run, m) if _val(r.get(k))), "") or "Unassigned"
    for rank in ("Species", "Genus", "Family", "Order", "Class", "Subdivision", "Phylum"):
        v = _val(r.get(rank))
        if not v:
            continue
        if rank == "Species":
            v = v.replace("_", " ")
            g = run._clean_genus(r.get("Genus"))
            if g and v == g:            # a species that is only its genus name (old VSEARCH naming)
                continue
            if " " not in v and g:
                v = f"{g} {v}"
        elif rank == "Genus":
            v = run._clean_genus(v) or v
        return group.replace("_", " "), rank, v.replace("_", " ")
    return "Unassigned", "", "Not named"


# ── samples, sites ────────────────────────────────────────────────────────────
def _common_suffix(names) -> str:
    names = [str(n) for n in names]
    if len(names) < 2:
        return ""
    s = names[0]
    for n in names[1:]:
        i = 0
        while i < min(len(s), len(n)) and s[-1 - i] == n[-1 - i]:
            i += 1
        s = s[len(s) - i:] if i else ""
        if not s:
            return ""
    # cut back to a separator, so a shared last replicate digit is not taken for the marker tag
    j = next((k for k, c in enumerate(s) if c in SEP), None)
    if j is not None and j > 0:
        s = s[j:]
    if not re.search(r"[A-Za-z]", s) or len(s) < 2:
        return ""
    if any(len(n) <= len(s) for n in names):
        return ""
    return s


def sample_keys(samples) -> dict:
    """{sample: key}: the sample name without the marker tag all of the job's names share."""
    suf = _common_suffix(samples)
    out = {}
    for s in samples:
        k = str(s)[: len(str(s)) - len(suf)] if suf else str(s)
        out[s] = k.rstrip(SEP) or str(s)
    return out


def site_of(key: str):
    """(site, replicate label, is_control) of a sample key."""
    m = CONTROL_RE.match(key)
    if m and m.group(1):
        return m.group(1).rstrip(SEP), m.group(2).upper(), True
    m = REP_RE.match(key)
    if m and m.group(1):
        return m.group(1).rstrip(SEP), m.group(2), False
    return key, "", False


# ── one marker ────────────────────────────────────────────────────────────────
class MarkerData:
    def __init__(self, run):
        self.run = run
        self.m = marker_of(run)
        self.target = target_of(run, self.m)
        self.keys = sample_keys(run.samples)                    # sample -> key
        self.by_fold = {rb._skey(k): s for s, k in self.keys.items()}
        comm = [q for q in run.bact_seqs if q in run.mat]
        self.comm_tot = {s: sum(run.mat[q].get(s, 0.0) for q in comm) for s in run.samples}
        self.taxa = OrderedDict()   # (group, rank, name) -> {sample: % of the sample's reads}
        self.reads = {}             # (group, rank, name) -> reads in all samples
        for q in comm:
            t = _feature_name(run, q, self.m)
            row = run.mat[q]
            d = self.taxa.setdefault(t, {})
            for s in run.samples:
                v = row.get(s, 0.0)
                if v > 0 and self.comm_tot.get(s):
                    d[s] = d.get(s, 0.0) + 100.0 * v / self.comm_tot[s]
            self.reads[t] = self.reads.get(t, 0.0) + sum(row.get(s, 0.0) for s in run.samples)

    def sample(self, key):
        return self.by_fold.get(rb._skey(key))

    def site_stats(self, t, reps):
        """(mean %, found in, of) of taxon t over the samples (keys) of one site."""
        ss = [self.sample(k) for k in reps]
        ss = [s for s in ss if s is not None]
        if not ss:
            return None
        vals = [self.taxa[t].get(s, 0.0) for s in ss]
        return sum(vals) / len(ss), sum(v > 0 for v in vals), len(ss)

    def overview(self):
        run = self.run
        inp = 0.0
        for t in run.tracking:
            try:
                inp += float(t.get("input") or 0)
            except (TypeError, ValueError):
                pass
        named = [t for t in self.taxa if t[1]]
        return {
            "inp": inp, "fin": sum(self.comm_tot.values()), "n_feat": len(run.mat),
            "n_taxa": len(named), "n_sp": sum(t[1] == "Species" for t in named),
            "n_gen": sum(t[1] in ("Species", "Genus") for t in named),
            "genus_pct": (100.0 * sum(v for t, v in self.reads.items() if t[1] in ("Species", "Genus"))
                          / max(1e-9, sum(self.reads.values()))),
        }

    def groups(self):
        """[(group, mean % of a sample's reads, taxa)] largest first."""
        S = [s for s in self.run.samples if self.comm_tot.get(s)]
        g = OrderedDict()
        for t, d in self.taxa.items():
            e = g.setdefault(t[0], [0.0, 0])
            e[0] += sum(d.get(s, 0.0) for s in S) / max(1, len(S))
            e[1] += 1 if t[1] else 0
        return sorted(((k, v[0], v[1]) for k, v in g.items()), key=lambda x: -x[1])


def pipeline_of(run) -> str:
    """The pipeline that made the job, as the figure headings name it."""
    m = marker_of(run)
    return {"cox1": "COX1 pipeline (DADA2 + codon filter)", "its": "ITS pipeline (DADA2)",
            "vsearch": "VSEARCH OTU pipeline", "emu": "ONT 16S pipeline (Emu)",
            "pacbio": "PacBio pipeline (DADA2)"}.get(run.kind, f"DADA2 {m} pipeline")


# ── figures side by side ──────────────────────────────────────────────────────
def _figure_pages(md) -> list:
    """One landscape page per topic, the markers' charts side by side, each headed with the
    marker, the job and the pipeline that made it."""
    ans = {}
    for d in md:
        try:
            ans[id(d)] = rb.Analysis(d.run)
        except Exception:
            ans[id(d)] = None

    def grp_rank(d):
        cols = list(next(iter(d.run.tax.values()), {}).keys()) + list(d.run.rank_cols)
        return next((r for r in _group_ranks(d.run, d.m) if r in cols), "Phylum")

    def safe(f):
        try:
            return f() or ""
        except Exception:
            return ""

    topics = [
        ("Reads through the pipeline steps",
         "Reads per sample after each step of each pipeline; the steps differ between pipelines (COX1 adds a "
         "length and a stop-codon filter).",
         lambda d, an: rb.svg_steps(d.run, an)),
        ("Composition &mdash; organism groups",
         "The 15 most abundant organism groups per sample (class for 12S, PR2 subdivision for 18S, phylum "
         "for COI), as % of each sample&rsquo;s reads for that marker.",
         lambda d, an: rb.svg_stacked(an, grp_rank(d), 15)),
        ("Composition &mdash; genera",
         "The 15 most abundant genera per sample, as % of each sample&rsquo;s reads for that marker.",
         lambda d, an: rb.svg_stacked(an, "Genus", 15)),
        ("Alpha diversity",
         "Observed ASVs, Chao1, Shannon and Simpson per sample. Each marker counts its own ASVs, so the "
         "values compare samples within a marker, not markers with each other.",
         lambda d, an: rb.svg_alpha(an, an.alpha())),
        ("Rarefaction",
         "ASVs found against sequencing depth; a curve that levels off means the depth was enough to see "
         "what that marker can see in the sample.",
         lambda d, an: rb.svg_rare(an)),
        ("Beta diversity &mdash; PCoA (Bray&ndash;Curtis)",
         "Samples close together hold a similar community for that marker.",
         lambda d, an: rb.svg_pcoa(an)),
    ]
    out = []
    for k, (title, note, f) in enumerate(topics):
        cells = []
        for d in md:
            an = ans[id(d)]
            svg = safe(lambda: f(d, an)) if an is not None else ""
            head = (f'<div class="sh"><b>{esc(d.lab)}</b> &middot; {esc(d.run.label)}</div>'
                    f'<div class="sp">{esc(pipeline_of(d.run))} &middot; {esc(ref_of(d.run))}</div>')
            cells.append(f"<td>{head}" + (svg or '<p class="dim">Not available for this job.</p>') + "</td>")
        first = ('<h2 style="margin-top:0">2. Figures side by side</h2><p>The same figure from each '
                 "marker&rsquo;s job, one topic to a page; each column is headed with the marker, the job and "
                 "the pipeline that made it.</p>") if k == 0 else ""
        out.append(f'<section class="wide side">{first}<h3>2{"ABCDEF"[k]}. {title}</h3><p class="dim">{note}</p>'
                   f'<table class="side"><tr>{"".join(cells)}</tr></table></section>')
    return out


# ── report ────────────────────────────────────────────────────────────────────
def _tbl(head, rows, left=(0,), cls="mmt"):
    """report_builder's table with the given columns (text) aligned left."""
    al = lambda i: ' style="text-align:left"' if i in left else ""
    h = "".join(f"<th{al(i)}>{c}</th>" for i, c in enumerate(head))
    b = "".join("<tr>" + "".join(f"<td{al(i)}>{c}</td>" for i, c in enumerate(r)) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'


def _n(n, one, many=None):
    return f"{n} {one if n == 1 else (many or one + 's')}"


def _name_html(rank, name):
    if rank in ("Species", "Genus"):
        return f"<i>{esc(name)}</i>"
    if not rank:
        return f'<span class="dim">{esc(name)}</span>'
    return f'{esc(name)} <span class="dim">({esc(rank.lower())})</span>'


def _cell(st):
    if st is None:
        return '<span class="dim">n/a</span>'
    mean, n, of = st
    if n == 0:
        return '<span class="dim">&ndash;</span>'
    p = f"{mean:.2f}" if mean < 10 else f"{mean:.1f}"
    return f"{p}%<br><span class=\"dim\">{n}/{of}</span>"


def build_multimarker_html(runs, title=None, app="NextGen-Amplicon", company="", subtitle=None) -> str:
    md = [MarkerData(r) for r in runs]
    names = [d.m for d in md]
    dup = {n for n in names if names.count(n) > 1}
    for d in md:                           # two runs of one marker: tell them apart by job
        d.lab = f"{d.m} ({d.run.label})" if d.m in dup else d.m

    # every sample key of every marker, folded so that case and punctuation do not split them
    keys = OrderedDict()
    for d in md:
        for s, k in d.keys.items():
            keys.setdefault(rb._skey(k), k)
    info = {f: site_of(k) for f, k in keys.items()}
    sites = OrderedDict()
    for f in sorted(keys, key=lambda f: (info[f][0].lower(), info[f][2], info[f][1].zfill(4))):
        sites.setdefault(info[f][0], {"reps": [], "ct": []})[
            "ct" if info[f][2] else "reps"].append(keys[f])
    controls = [k for v in sites.values() for k in v["ct"]]
    n_field = sum(len(v["reps"]) for v in sites.values())

    rep_title = "Multi-marker Amplicon Comparison Report"
    project = title if title and title.strip() and title.strip() != rep_title else None
    title = rep_title
    mk = " / ".join(d.lab for d in md)
    subtitle = subtitle or (f"{mk} · {len(sites)} site{'s' if len(sites) != 1 else ''} · "
                            f"{n_field} field sample{'s' if n_field != 1 else ''}"
                            + (f" + {len(controls)} control{'s' if len(controls) != 1 else ''}" if controls else ""))
    tab = [0]

    def T(cap):
        tab[0] += 1
        return f'<p class="tcap"><b>Table {tab[0]}.</b> {cap}</p>'

    out = [f'<div class="eyebrow">{esc(mk)} &middot; {len(keys)} samples</div>', f"<h1>{esc(title)}</h1>"]
    if project:
        out.append(f'<p class="project">{esc(project)}</p>')
    out.append(f'<p class="sub">{esc(subtitle)}</p>')
    ov = {id(d): d.overview() for d in md}
    allsp = {t[2] for d in md for t in d.taxa if t[1] == "Species"}
    kpis = [(str(len(md)), "markers"), (str(len(sites)), "sampling sites"),
            (str(sum(ov[id(d)]["n_taxa"] for d in md)), "taxa named (all markers)"),
            (str(len(allsp)), "named to species")]
    out.append('<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="v">{v}</div><div class="l">{esc(l)}</div></div>' for v, l in kpis) + "</div>")
    out.append(
        "<p>The same samples were sequenced with " + ", ".join(f"<b>{esc(d.lab)}</b>" for d in md)
        + ". Each marker was analysed as its own job; this report sets their results side by side. "
        "Each marker amplifies a different part of the community, so a taxon missing from one marker "
        "is often outside what that marker&rsquo;s primers and reference cover rather than absent from the "
        "sample.</p>")

    # ── 1. overview per marker ────────────────────────────────────────────────
    out.append("<h2>1. Overview per marker</h2>")
    rows = []
    for d in md:
        o = ov[id(d)]
        rows.append([f"<b>{esc(d.lab)}</b>",
                     esc(d.target), esc(ref_of(d.run)) + f"<br><span class=\"dim\">{esc(naming_of(d.run))}</span>",
                     str(len(d.run.samples)), fmt(o["inp"]) if o["inp"] else "—", fmt(o["fin"]),
                     fmt(o["n_feat"]), str(o["n_taxa"]), str(o["n_sp"]), f"{o['genus_pct']:.1f}%"])
    out.append(T("The markers. <i>Reads</i>: sequenced, and in the final table after quality filtering, "
                 "merging and chimera removal (organelle reads of 18S set aside). <i>Taxa</i>: distinct names at "
                 "their deepest rank; <i>species</i>: of these, named to species. <i>Genus %</i>: share of the "
                 "final reads named to genus or species."))
    out.append(_tbl(["Marker", "Target", "Reference / naming", "Samples", "Reads<br>sequenced", "Final<br>reads",
                     "ASVs", "Taxa", "Species", "Genus<br>%"], rows, left=(0, 1, 2), cls="mmt ov"))

    rows = []
    for d in md:
        gs = d.groups()
        for i, (g, p, n) in enumerate(gs[:6]):
            rows.append([f"<b>{esc(d.lab)}</b>" if i == 0 else "", esc(g), f"{p:.1f}%", str(n)])
        if len(gs) > 6:
            rest = sum(p for _g, p, _n in gs[6:])
            rows.append(["", f'<span class="dim">{len(gs) - 6} other groups</span>', f"{rest:.1f}%",
                         str(sum(n for _g, _p, n in gs[6:]))])
    out.append(T("What each marker sees: its main organism groups (fish and other vertebrates by class for "
                 "12S; PR2 subdivision for 18S; phylum otherwise), as the mean share of a sample&rsquo;s reads, "
                 "and the number of taxa named in each."))
    out.append(_tbl(["Marker", "Organism group", "Mean % of reads", "Taxa"], rows, left=(0, 1)))

    rows = []
    for f, k in keys.items():
        site, rep, ct = info[f]
        cells = []
        for d in md:
            s = d.sample(k)
            cells.append(fmt(d.comm_tot.get(s, 0)) if s is not None else
                         '<span class="hl">not in this job</span>')
        rows.append([esc(k), esc(site), "control" if ct else esc(rep or "—")] + cells)
    out.append(T("The samples of each marker, matched by name without the marker part each job adds, with "
                 "their final reads. A sample marked <i>not in this job</i> was not sequenced or not uploaded "
                 "for that marker; its site is then compared over the samples it has."))
    out.append(_tbl(["Sample", "Site", "Replicate"] + [esc(d.lab) for d in md], rows, left=(0, 1, 2)))
    low = [(k, d.lab) for d in md for s, k in d.keys.items() if 0 < d.comm_tot.get(s, 0) < 1000
           and not site_of(k)[2]]
    if low:
        out.append('<div class="note"><b>Few reads.</b> ' + "; ".join(
            f"{esc(k)} ({esc(m)})" for k, m in low[:12]) + (" and others" if len(low) > 12 else "")
            + " kept fewer than 1,000 final reads: what these samples show is less certain.</div>")

    out += _figure_pages(md)

    # ── 3. what was found at each site ────────────────────────────────────────
    out.append('<h2 class="pbreak">3. What was found at each site</h2>')
    out.append(
        "<p>A taxon counts as found at a site when at least one read of it is in one of the site&rsquo;s "
        "field replicates. Its abundance is the mean, over the site&rsquo;s replicates, of its share of each "
        "replicate&rsquo;s reads <i>for that marker</i>: shares can be compared between sites within a "
        "marker, not between markers. The number of replicates it was found in (e.g. 3/3) is the stronger "
        "evidence that it is there."
        + (" Controls are not counted in a site; the taxa they carried are marked." if controls else "") + "</p>")
    rows = []
    for site, v in sites.items():
        cells = []
        for d in md:
            found = [t for t in d.taxa if t[1] and (d.site_stats(t, v["reps"]) or (0, 0, 0))[1] > 0]
            if not any(d.sample(k) is not None for k in v["reps"]):
                cells.append('<span class="dim">n/a</span>')
                continue
            cells.append(f"{len(found)}<br><span class=\"dim\">{sum(t[1] == 'Species' for t in found)} species</span>")
        rows.append([f"<b>{esc(site)}</b>", str(len(v["reps"]))] + cells)
    out.append(T("Taxa found at each site by each marker (named taxa; of these, named to species)."))
    out.append(_tbl(["Site", "Replicates"] + [esc(d.lab) for d in md], rows))

    # site by site: the names, marker by marker
    for site, v in sites.items():
        out.append(f"<h3>{esc(site)}</h3>")
        items = []
        for d in md:
            st = [(t, d.site_stats(t, v["reps"])) for t in d.taxa if t[1]]
            st = [(t, s) for t, s in st if s and s[1] > 0]
            if not any(d.sample(k) is not None for k in v["reps"]):
                items.append(f"<li><b>{esc(d.lab)}</b>: no sample of this site in the job.</li>")
                continue
            if not st:
                items.append(f"<li><b>{esc(d.lab)}</b>: nothing named.</li>")
                continue
            st.sort(key=lambda x: (-x[1][1], -x[1][0]))
            sp = [x for x in st if x[0][1] == "Species"]
            show = st[:MAX_NAMES]
            more = len(st) - len(show)
            what = "taxon" if more == 1 else "taxa"
            items.append(
                f"<li><b>{esc(d.lab)}</b> ({_n(len(st), 'taxon', 'taxa')}, {len(sp)} to species): "
                + ", ".join(f"{_name_html(t[1], t[2])} <span class=\"dim\">{s[1]}/{s[2]}, {s[0]:.1f}%</span>"
                            for t, s in show)
                + (f", and {more} more {what}" if more > 0 else "") + ".</li>")
        out.append("<ul>" + "".join(items) + "</ul>")

    # one table per marker: taxa x sites
    out.append('<h2 class="pbreak">4. Taxa by site, marker by marker</h2>')
    for d in md:
        allt = [t for t in d.taxa]
        stats = {t: {site: d.site_stats(t, v["reps"]) for site, v in sites.items()} for t in allt}
        cts = [k for k in controls if d.sample(k) is not None]
        found = [t for t in allt if any(s and s[1] > 0 for s in stats[t].values())]
        top = lambda t: max((s[0] for s in stats[t].values() if s), default=0.0)
        gsum = {}
        for t in found:
            gsum[t[0]] = gsum.get(t[0], 0.0) + top(t)
        found.sort(key=lambda t: (t[0] == "Unassigned", -gsum[t[0]], -top(t)))
        shown = sorted(found, key=lambda t: -top(t))[:MAX_ROWS]
        keep = set(shown)
        shown = [t for t in found if t in keep]
        cut = len(found) - len(shown)
        head = (["Organism group", "Taxon"] + [esc(s) for s in sites]
                + ([f"In controls ({len(cts)})"] if cts else []))
        body, last = [], None
        for t in shown:
            ctn = sum(1 for k in cts if d.taxa[t].get(d.sample(k), 0) > 0)
            body.append("<tr>" + f"<td>{esc(t[0]) if t[0] != last else ''}</td><td>{_name_html(t[1], t[2])}</td>"
                        + "".join(f"<td>{_cell(stats[t][s])}</td>" for s in sites)
                        + (f"<td>{'<span class=hl>' + str(ctn) + '</span>' if ctn else '&ndash;'}</td>" if cts else "")
                        + "</tr>")
            last = t[0]
        out.append(f'<h3>{esc(d.lab)} &mdash; {esc(d.target)}</h3>')
        out.append(T(f"{esc(d.lab)} ({esc(d.run.label)}): every taxon by site &mdash; mean % of the site&rsquo;s "
                     "reads for this marker, and the replicates it was found in; &ndash; not found; n/a no sample of "
                     "the site in this job."
                     + (" <i>In controls</i>: number of controls that also carried it." if cts else "")
                     + (f" The {len(shown)} most abundant of {len(found)} taxa are listed; the other {cut} "
                        "are in the job&rsquo;s own report and tables." if cut else "")))
        out.append('<table class="mm"><thead><tr>' + "".join(f"<th>{h}</th>" for h in head)
                   + "</tr></thead><tbody>" + "".join(body) + "</tbody></table>")

    # ── methods and notes ─────────────────────────────────────────────────────
    out.append("<h2>5. Methods and notes</h2>")
    notes = []
    for d in md:
        p = d.run.params
        pr = " / ".join(x for x in (p.get("primer_f"), p.get("primer_r")) if x)
        notes.append(f"<b>{esc(d.lab)}</b> &mdash; job {esc(d.run.label)}; {esc(naming_of(d.run))} against "
                     f"{esc(d.run.db_name)}" + (f"; primers <span class=\"mono\">{esc(pr)}</span>" if pr else "")
                     + ". Its own report gives the full method, read tracking and diversity.")
    notes.append("Samples were matched between jobs by name, without the marker tag every name of a job shares "
                 "(e.g. &ldquo;-Riaz12S&rdquo;, &ldquo;_TAR&rdquo;, &ldquo;_COI-XT&rdquo;); the site is the "
                 "name without its replicate number, and a name ending in CT, NC, NTC, NEG or BLANK is that "
                 "site&rsquo;s control.")
    notes.append("Taxa are named at the deepest rank the marker&rsquo;s classification reached, so one organism "
                 "may appear as a species under one marker and only as a genus or family under another; the "
                 "references differ (fish mitogenomes, PR2, MIDORI2) and so do their names for the same group.")
    notes.append("Relative abundance from amplicons reflects the DNA each organism shed and how well the primers "
                 "fit it, not the number or biomass of organisms.")
    out.append("<ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")

    css = (rb._css(app, company, f"Multi-marker · {mk}") + rb._CSS_FULL.replace("@@TITLE@@", esc(title))
           + """
table.mm { font-size: 7.8pt; }
table.mmt td { vertical-align: top; line-height: 1.3; }
table.ov { font-size: 7.6pt; table-layout: fixed; }
section.side { break-before: page; }
table.side { table-layout: fixed; width: 100%; border-collapse: collapse; margin-top: 1mm; }
table.side td { vertical-align: top; text-align: left; padding: 0 2.5mm 0 0; border: none; }
table.side td + td { border-left: .5px solid #D5DEE2; padding-left: 2.5mm; }
table.side .sh { font-size: 8.6pt; } table.side .sp { font-size: 7pt; color: #6B7F88; margin-bottom: 1mm; }
table.side svg { width: 100%; height: auto; }
table.ov th, table.ov td { padding-right: 1.4mm; overflow-wrap: anywhere; }
table.ov th:nth-child(1) { width: 9%; } table.ov th:nth-child(2) { width: 15%; }
table.ov th:nth-child(3) { width: 24%; } table.ov th:nth-child(4) { width: 7%; }
table.ov th:nth-child(5), table.ov th:nth-child(6) { width: 10%; }
table.ov th:nth-child(7), table.ov th:nth-child(8), table.ov th:nth-child(9) { width: 6%; }
table.ov th:nth-child(10) { width: 7%; }
table.mm td { vertical-align: top; line-height: 1.25; }
table.mm td:nth-child(2) { text-align: left; }
table.mm thead th:nth-child(2) { text-align: left; }
""")
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f"<title>{esc(title)}</title><style>{css}</style></head><body>"
            + "".join(out) + "</body></html>")


def is_multimarker(runs) -> bool:
    """Runs of more than one marker (12S + 18S + COI ...): this report, not the reference comparison."""
    return len(runs) > 1 and len({marker_of(r) for r in runs}) > 1
