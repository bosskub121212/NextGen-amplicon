"""Result previews on the job page: what each of the eight folders of Download results will hold,
drawn from the finished run so the BI team can look before building anything.

build(run_dir) -> {"sections": [{"key", "title", "summary", "items": [item, ...]}, ...]}
item: {"kind": "svg", "title", "svg"} | {"kind": "table", "title", "head", "rows", "more"}
      | {"kind": "files", "title", "files": [{"name", "size"}]} | {"kind": "note", "text"}

The figures are the report's own (report_builder), so a preview is the picture the report and the
result's 5.Taxa / 6.Diversity PNGs carry. Every section is built on its own: one that fails says
so and the others still show. The result is cached in the run folder (.preview_cache.json) until
a file it was made from changes (metadata.csv included: groups edited after the run).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

CACHE = "packages/preview_cache.json"   # packages/ is left out of Everything (.zip), and a file there
                                        # does not make the whole-job zip look out of date
SRC = ("asv_table.csv", "taxonomy_table.csv", "taxonomy.csv", "read_tracking.csv", "metadata.csv",
       "summary.json", "run_params.json", "blast/taxonomy_blast.csv")
MAX_ROWS = 40


def _key(run_dir: Path) -> str:
    try:
        ver = json.loads((Path(__file__).resolve().parent.parent / "version.json").read_text()).get("version", "")
    except Exception:
        ver = ""
    h = hashlib.sha1(("v3|" + ver).encode())
    for f in SRC:
        p = run_dir / f
        if p.exists():
            st = p.stat()
            h.update(f"{f}|{st.st_size}|{int(st.st_mtime)}".encode())
    return h.hexdigest()


def _fmt(n) -> str:
    try:
        return f"{int(round(float(n))):,}"
    except (TypeError, ValueError):
        return "—" if n in (None, "") else str(n)


def _size(n: int) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:.0f} {u}" if u == "B" else f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} GB"


def _table(title, head, rows, note=""):
    rows = [[("" if c is None else str(c)) for c in r] for r in rows]
    return {"kind": "table", "title": title, "head": [str(h) for h in head], "rows": rows[:MAX_ROWS],
            "more": max(0, len(rows) - MAX_ROWS), "note": note}


def _svg(title, svg):
    return {"kind": "svg", "title": title, "svg": svg} if svg else None


def build(run_dir: Path, use_cache: bool = True) -> dict:
    run_dir = Path(run_dir)
    key = _key(run_dir)
    cp = run_dir / CACHE
    if use_cache and cp.exists():
        try:
            c = json.loads(cp.read_text(encoding="utf-8"))
            if c.get("key") == key:
                return c["preview"]
        except Exception:
            pass
    out = _build(run_dir)
    try:
        cp.parent.mkdir(exist_ok=True)
        cp.write_text(json.dumps({"key": key, "preview": out}), encoding="utf-8")
    except OSError:
        pass
    return out


def _build(run_dir: Path) -> dict:
    import amp_data as AD
    import report_builder as RB
    import custom_report as CR
    run = AD.load_run(run_dir)
    rrun = RB.Run(run_dir, RB.KITOME, run_dir.name)
    an = RB.Analysis(rrun)
    S = list(an.samples)
    engine = getattr(run, "engine", run.kind)
    feat = {"emu": "Taxa_table", "vsearch": "OTU"}.get(engine, "ASV")
    unit = {"emu": "taxa", "vsearch": "OTUs"}.get(engine, "ASVs")
    step3 = {"emu": "Emu", "vsearch": "VSEARCH"}.get(engine, "DADA2")
    group = {s: an.groups[s] for s in S if s in an.groups}
    many = len(S) > RB.MANY
    first = (RB._chunks(S)[0] if many else S)
    secs = []

    def section(key, title, fn):
        try:
            items, summary = fn()
        except Exception as ex:                      # one section's trouble never hides the others
            items, summary = [{"kind": "note", "text": f"Not available for this run ({type(ex).__name__}: {ex})"}], ""
        secs.append({"key": key, "title": title, "summary": summary, "items": [i for i in items if i]})

    # 1. QC
    def qc():
        cols = [c for c in run.track_cols if c != "perc"]
        names = {"raw_data": "Raw", "qc_filtered": "QC", "primer_trimmed": "Primers", "final": "Final",
                 "assigned": "Assigned", "trimmed": "Trimmed", "filtered": "Filtered", "denoisedF": "Denoised F",
                 "denoisedR": "Denoised R", "merged": "Merged", "Nonchim": "Non-chimeric"}
        have = [c for c in cols if any(run.tracking.get(s, {}).get(c) not in (None, 0) for s in S)]
        rows = [[s] + ([group.get(s, "")] if group else []) + [_fmt(run.tracking.get(s, {}).get(c)) for c in have]
                + [(f"{run.tracking[s]['perc']}%" if run.tracking.get(s, {}).get("perc") is not None else "—")]
                for s in S]
        kept = [run.tracking.get(s, {}).get("perc") for s in S if run.tracking.get(s, {}).get("perc") is not None]
        ns = f"{len(S)} sample{'s' if len(S) != 1 else ''}"
        rng = (f"{min(kept):.0f}%" if round(min(kept)) == round(max(kept)) else f"{min(kept):.0f}–{max(kept):.0f}%") if kept else ""
        summ = f"{ns} · reads kept {rng}" if kept else ns
        return [_table("Reads at each step", ["Sample"] + (["Group"] if group else []) + [names.get(c, c) for c in have]
                       + ["Kept"], rows),
                _svg("Read counts at each step", RB.svg_steps(rrun, an)),
                _svg("Reads surviving each step", RB.svg_readflow(rrun))], summ

    # 2. Cutadapt
    def cut():
        rows = []
        for s in S:
            j = run._json(f"cutadapt/{s}.cutadapt.json").get("read_counts", {}) if hasattr(run, "_json") else {}
            if j:
                inp, outp = j.get("input") or 0, j.get("output") or 0
                hit = j.get("read1_with_adapter")
                rows.append([s, _fmt(inp), (f"{100 * hit / inp:.1f}%" if inp and hit is not None else "—"), _fmt(outp),
                             (f"{100 * outp / inp:.1f}%" if inp else "—")])
        if not rows:
            return [{"kind": "note", "text": ("Nanopore runs trim the primers inside the pipeline; the reads with both "
                                               "primers are counted in 1.QC (Primers).") if run.kind == "emu" or engine == "vsearch"
                     else "No cutadapt counts in this run."}], "no cutadapt counts"
        return [_table("Primer removal (cutadapt)", ["Sample", "Reads in", "With primer", "Reads out", "Kept"], rows)], \
            f"{len(rows)} samples trimmed"

    # 3. denoising / classification
    def s3():
        items = []
        if engine == "emu":
            nt = {}
            for a in run.asvs:
                if a.get("not_taxon"):
                    for s in S:
                        nt.setdefault(s, {})[a.get("id")] = a["counts"].get(s, 0)
            rows = [[s, sum(1 for a in run.asvs if not a.get("not_taxon") and a["counts"].get(s, 0) > 0),
                     _fmt(sum(a["counts"].get(s, 0) for a in run.asvs if not a.get("not_taxon"))),
                     _fmt(nt.get(s, {}).get("unmapped", 0)), _fmt(nt.get(s, {}).get("mapped_unclassified", 0))]
                    for s in S]
            items.append(_table("Emu per sample", ["Sample", "Taxa", "Assigned reads", "Unmapped", "Mapped, unclassified"],
                                rows))
            summ = f"{sum(1 for a in run.asvs if not a.get('not_taxon'))} taxa"
        else:
            n = len(run.asvs)
            summ = f"{n} {unit}"
            items.append({"kind": "note", "text": f"{n:,} {unit} in the run" + (" (VSEARCH centroids, otus.fasta)"
                                                                                  if engine == "vsearch" else "") + "."})
            pdfs = sorted(run_dir.glob("Dada2/*.pdf"))
            if pdfs:
                items.append({"kind": "files", "title": "Error models and quality profiles (PDF)",
                              "files": [{"name": p.name, "size": _size(p.stat().st_size)} for p in pdfs]})
        return items, summ

    # 4. feature table
    def s4():
        per = [[s] + ([group.get(s, "")] if group else []) +
               [_fmt(sum(a["counts"].get(s, 0) for a in run.asvs)),
                sum(1 for a in run.asvs if a["counts"].get(s, 0) > 0 and not a.get("not_taxon"))] for s in S]
        tot = sum(a["total"] for a in run.asvs) or 1
        top = [[a["sv"], "Unassigned" if a.get("not_taxon") else AD.label(a, "Species"),
                "" if a.get("not_taxon") else AD.label(a, "Phylum"), _fmt(a["total"]), f"{100 * a['total'] / tot:.2f}%"]
               for a in run.asvs[:15]]
        return [_table(f"{unit[0].upper() + unit[1:]} per sample", ["Sample"] + (["Group"] if group else [])
                       + ["Reads", unit[0].upper() + unit[1:]], per),
                _table(f"Most abundant {unit} (IDs of the count table)", ["ID", "Name", "Phylum", "Reads", "Share"], top)], \
            f"{len(run.asvs)} rows in the count table"

    # 5. taxa
    def s5():
        items = []
        note = f" (first {len(first)} of {len(S)} samples)" if many else ""
        for rank in ("Phylum", "Genus"):
            if rank in rrun.rank_cols:
                svg = (RB.svg_stacked(an, rank, 30, samples=first, wide=True, groups=None) if many
                       else RB.svg_stacked(an, rank, 30))
                items.append(_svg(f"Top 30 {rank}{note}", svg))
        if "Genus" in rrun.rank_cols:
            items.append(_svg(f"Genus heatmap{note}", RB.svg_heatmap(an, "Genus", 30, samples=first, wide=True) if many
                              else RB.svg_heatmap(an, "Genus", 30)))
        gt = run.rank_table("Genus")
        tot = sum(gt.values()) or 1
        gen = [[k, _fmt(v), f"{100 * v / tot:.2f}%"] for k, v in list(gt.items())[:15]]
        items.append(_table("Top genera, whole run", ["Genus", "Reads", "Share"], gen))
        return items, f"{len(gt)} genera"

    # 6. diversity
    def s6():
        al = an.alpha()
        keys = ["Observed", "Chao1", "Shannon", "Simpson"]
        rows = [[s] + ([group.get(s, "")] if group else []) + [_fmt(an.comm_tot.get(s, 0))]
                + [(f"{al[s][k]:.3f}" if isinstance(al.get(s, {}).get(k), float) else al.get(s, {}).get(k, "—"))
                   for k in keys] for s in S]
        items = [_table("Alpha diversity", ["Sample"] + (["Group"] if group else []) + ["Reads"] + keys, rows),
                 _svg("Alpha diversity", RB.svg_alpha(an, al)), _svg("Rarefaction curves", RB.svg_rare(an))]
        if len(S) >= 3:
            items.append(_svg("PCoA (Bray–Curtis)", RB.svg_pcoa(an, wide=many)))
        if 2 <= len(S) <= 40:
            items.append(_svg("Bray–Curtis dissimilarity", RB.svg_bc_heat(an)))
        pm = an.permanova_groups() if group else None
        if pm:
            items.append(_table("PERMANOVA by group", ["Groups", "Samples", "R²", "F", "p"],
                                [[pm.get("k"), pm.get("n"), f"{pm.get('R2', 0):.3f}", f"{pm.get('F', 0):.2f}",
                                  f"{pm.get('p', 0):.3f}"]]))
        else:
            items.append({"kind": "note", "text": ("PCoA needs 3 samples or more. " if len(S) < 3 else "")
                          + "PERMANOVA needs 2 groups with 2 samples or more in one of them."})
        return items, f"Shannon {min(al[s]['Shannon'] for s in S):.2f}–{max(al[s]['Shannon'] for s in S):.2f}" if S else ""

    # 7. tables
    def s7():
        ranks = [r for r in AD.RANKS if r in rrun.rank_cols]
        return [{"kind": "note", "text": "Download results writes 7.Tables/<name>.xlsx: read tracking, the profiles at "
                 f"every rank ({', '.join(ranks)}) and alpha diversity, a sheet each; 5.Taxa/csv holds the same "
                 "profiles as CSV."}], f"{len(ranks)} ranks"

    # 8. reads
    def s8():
        files = []
        for s in S:
            for q in CR._filtered(run, s):
                files.append({"name": q.name, "size": _size(q.stat().st_size), "b": q.stat().st_size})
        if not files:
            return [{"kind": "note", "text": "No quality-filtered reads kept in this run folder."}], "none"
        total = sum(f.pop("b") for f in files)
        return [{"kind": "files", "title": "Quality-filtered reads (Download results → 8.Reads)", "files": files}], \
            f"{len(files)} file{'s' if len(files) != 1 else ''} · {_size(total)}"

    section("1.QC", "1.QC — reads at each step", qc)
    section("2.Cutadapt", "2.Cutadapt — primer removal", cut)
    section(f"3.{step3}", f"3.{step3} — " + {"Emu": "Emu classification", "VSEARCH": "OTU clustering"}.get(step3, "denoising"), s3)
    section(f"4.{feat}", f"4.{feat} — the count table", s4)
    section("5.Taxa", "5.Taxa — composition", s5)
    section("6.Diversity", "6.Diversity — alpha and beta", s6)
    section("7.Tables", "7.Tables — Excel", s7)
    section("8.Reads", "8.Reads — filtered reads", s8)
    return {"sections": secs, "samples": len(S), "engine": engine}
