#!/usr/bin/env python3
"""Read-based taxonomy for ONT WGS isolates (report section "Read-based Taxonomy Classification").

The QC-passed reads of a sample (reads.filtered.fastq.gz, ≥1000 bp / ≥Q10 by default) are
classified by up to three independent tools, each against its own database:

  Kraken2 + Bracken   k-mer LCA + Bayesian re-estimation      db_paths.json  kraken2_db
  Centrifuger         compressed-BWT classifier (or Centrifuge) db_paths.json  readtax_centrifuge
  GOTTCHA2 (-np)      unique-signature mapping, ONT mode       db_paths.json  readtax_gottcha2

Each tool's species-level result is turned into a share (%) of what that tool assigned at
species level, so the columns are comparable in the heatmap; the tools' own reports are kept
next to the table. A tool or database that is not installed is skipped with a note — the
table then simply has fewer columns.

Output in samples/<S>/readtax/: kraken2.report, bracken_species.tsv, centrifuger_report.tsv,
gottcha2.tsv (+ .full.tsv, .lineage.tsv), read_taxonomy.tsv, heatmap_read_taxonomy.svg,
readtax_info.json.

CLI (also used for a finished job — "Add read-based taxonomy"):
  readtax.py --results <job results dir> --work_dir <uploads>/_wgs_work --db_paths db_paths.json
             [--samples "S1,S2"] [--threads 8]
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ont_wgs_pipeline as W  # noqa: E402  (find_tool / run / log)

ORDER = ["kraken2", "centrifuge", "gottcha2"]
OUTPUTS = ("kraken2.report", "bracken_species.tsv", "bracken_species.report", "centrifuger_report.tsv",
           "centrifuger_classification.tsv", "centrifuge_report.tsv", "gottcha2.tsv", "gottcha2.full.tsv",
           "gottcha2.lineage.tsv", "read_taxonomy.tsv", "heatmap_read_taxonomy.svg", "readtax_info.json")
INFO = "readtax_info.json"
TABLE = "read_taxonomy.tsv"


def esc(s) -> str:
    return html.escape(str(s), quote=True)


def _norm(name: str) -> str:
    return re.sub(r"\s+", " ", str(name).replace("_", " ")).strip()


def _share(d: dict) -> dict:
    tot = sum(v for v in d.values() if v > 0)
    return {k: 100.0 * v / tot for k, v in d.items() if v > 0} if tot > 0 else {}


def _rows(path: Path, delim="\t"):
    try:
        with open(path, newline="") as fh:
            yield from csv.DictReader(fh, delimiter=delim)
    except (OSError, UnicodeDecodeError):
        return


# ── parsers (all return {species name: % of the tool's species-level assignments}) ──
def parse_bracken(path: Path) -> dict:
    d = {}
    for r in _rows(path):
        try:
            if r.get("taxonomy_lvl", "S") == "S":
                d[_norm(r["name"])] = d.get(_norm(r["name"]), 0) + float(r["new_est_reads"])
        except (KeyError, ValueError, TypeError):
            pass
    return _share(d)


def _kraken_lines(path: Path):
    try:
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 6:
                    yield p
    except OSError:
        return


def parse_kraken_species(path: Path) -> dict:
    """Fallback when Bracken is missing/failed: clade reads of rank-S rows."""
    d = {}
    for p in _kraken_lines(path):
        rank, name = p[-3], p[-1].strip()
        if rank == "S":
            try:
                d[_norm(name)] = d.get(_norm(name), 0) + float(p[1])
            except ValueError:
                pass
    return _share(d)


def kraken_classified_pct(path: Path) -> float | None:
    u = c = None
    for p in _kraken_lines(path):
        try:
            if p[-3] == "U":
                u = float(p[1])
            elif p[-3] == "R":
                c = float(p[1])
        except ValueError:
            pass
    if c is None:
        return None
    return 100.0 * c / ((u or 0) + c) if (u or 0) + c else None


def parse_centrifuge(path: Path) -> dict:
    """Centrifuge / Centrifuger report: reads under each species-rank node."""
    d = {}
    for r in _rows(path):
        try:
            if (r.get("taxRank") or "").strip().lower() == "species":
                n = _norm(r["name"])
                d[n] = d.get(n, 0) + float(r["numReads"])
        except (KeyError, ValueError, TypeError):
            pass
    return _share(d)


def parse_gottcha2(path: Path) -> dict:
    d = {}
    for r in _rows(path):
        try:
            if (r.get("LEVEL") or "").strip().lower() == "species":
                n = _norm(r["NAME"])
                d[n] = d.get(n, 0) + float(r["REL_ABUNDANCE"])
        except (KeyError, ValueError, TypeError):
            pass
    return _share(d)


# ── combined table + heatmap ──────────────────────────────────────────────────
def column_label(key: str, label: str | None = None) -> str:
    return label or {"kraken2": "Kraken2+Bracken", "gottcha2": "GOTTCHA2", "centrifuge": "Centrifuge"}[key]


def write_table(profiles: dict, dest: Path, labels: dict | None = None) -> Path:
    labels = labels or {}
    keys = [k for k in ORDER if profiles.get(k)]
    taxa = sorted({t for k in keys for t in profiles[k]},
                  key=lambda t: (-max(profiles[k].get(t, 0) for k in keys), t))
    dest.parent.mkdir(parents=True, exist_ok=True)

    def fmt(v):
        return "0" if not v else (f"{v:.4f}".rstrip("0").rstrip("."))
    with open(dest, "w") as fh:
        fh.write("\t".join(["Taxon"] + [f"{column_label(k, labels.get(k))} (%)" for k in keys]) + "\n")
        for t in taxa:
            fh.write("\t".join([t] + [fmt(profiles[k].get(t, 0)) for k in keys]) + "\n")
    return dest


def read_table(tsv: Path) -> tuple[list[str], list[tuple[str, list[float]]]]:
    try:
        lines = [l.rstrip("\n").split("\t") for l in Path(tsv).read_text().splitlines() if l.strip()]
    except OSError:
        return [], []
    if not lines:
        return [], []
    cols = [re.sub(r"\s*\(%\)$", "", c) for c in lines[0][1:]]
    rows = []
    for p in lines[1:]:
        try:
            rows.append((p[0], [float(x) for x in p[1:1 + len(cols)]]))
        except ValueError:
            pass
    return cols, rows


def _colour(v: float) -> str:
    """White → dark teal on a log scale, 0.01 % … 100 %."""
    if v <= 0:
        return "#ffffff"
    f = (math.log10(max(v, 0.01)) + 2) / 4          # 0.01 % → 0, 100 % → 1
    f = min(1.0, max(0.0, f))
    a, b = (247, 252, 253), (0, 88, 96)
    return "#%02x%02x%02x" % tuple(round(a[i] + (b[i] - a[i]) * f) for i in range(3))


def svg_heatmap(tsv: Path, top: int = 20) -> str:
    cols, rows = read_table(tsv)
    if not cols or not rows:
        return ""
    rows = rows[:top]
    cw, ch, lw, th = 92, 17, 230, 64
    w = lw + cw * len(cols) + 120
    h = th + ch * len(rows) + 40
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
         f'font-family="DejaVu Sans" font-size="9"><rect width="{w}" height="{h}" fill="#fff"/>']
    for j, c in enumerate(cols):
        x = lw + j * cw + cw / 2
        o.append(f'<text x="{x:.1f}" y="{th - 8}" text-anchor="middle" font-weight="bold">{esc(c)}</text>')
    for i, (name, vals) in enumerate(rows):
        y = th + i * ch
        o.append(f'<text x="{lw - 6}" y="{y + ch * 0.7:.1f}" text-anchor="end" font-style="italic">'
                 f'{esc(name)}</text>')
        for j, v in enumerate(vals):
            x = lw + j * cw
            fill = _colour(v)
            ink = "#fff" if v >= 10 else "#222"
            o.append(f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" fill="{fill}" stroke="#d0d0d0" '
                     f'stroke-width=".5"/><text x="{x + cw / 2}" y="{y + ch * 0.7:.1f}" text-anchor="middle" '
                     f'fill="{ink}">{"–" if v <= 0 else (f"{v:.2f}" if v < 10 else f"{v:.1f}")}</text>')
    # legend
    lx, ly = lw + cw * len(cols) + 24, th
    for k, v in enumerate((100, 10, 1, 0.1, 0.01)):
        o.append(f'<rect x="{lx}" y="{ly + k * 16}" width="14" height="14" fill="{_colour(v)}" stroke="#999" '
                 f'stroke-width=".5"/><text x="{lx + 20}" y="{ly + k * 16 + 10}">{v:g} %</text>')
    o.append(f'<text x="{lx}" y="{ly - 8}" font-weight="bold">log scale</text>')
    o.append(f'<text x="{lw}" y="{h - 12}" fill="#555">Share of species-level assignments per tool (%); '
             f'top {len(rows)} taxa.</text></svg>')
    return "".join(o)


# ── running the tools ─────────────────────────────────────────────────────────
def kmer_lengths(db: Path) -> list[int]:
    return sorted(int(m.group(1)) for f in db.glob("database*mers.kmer_distrib")
                  if (m := re.match(r"database(\d+)mers", f.name)))


def pick_read_len(avail: list[int], mean_len: float) -> int | None:
    if not avail:
        return None
    below = [a for a in avail if a <= mean_len]
    return max(below) if below else min(avail)


def centrifuge_index(prefix: str) -> tuple[str, str] | None:
    """('centrifuger'|'centrifuge', prefix) from an index prefix or a folder holding one."""
    if not prefix:
        return None
    p = Path(prefix)
    cands = [p] + (sorted({Path(str(f)[:-6]) for f in p.glob("*.1.cfr")} |
                          {Path(str(f)[:-5]) for f in p.glob("*.1.cf")}) if p.is_dir() else [])
    for c in cands:
        if Path(str(c) + ".1.cfr").exists():
            return "centrifuger", str(c)
        if Path(str(c) + ".1.cf").exists():
            return "centrifuge", str(c)
    return None


def gottcha_db(prefix: str) -> str | None:
    if not prefix:
        return None
    p = Path(prefix)
    if p.is_dir():
        hits = sorted(p.glob("*.fna.tax.tsv"))
        return str(hits[0])[: -len(".tax.tsv")] if len(hits) == 1 else None
    return str(p) if Path(str(p) + ".tax.tsv").exists() else None


def _db_label(path: str) -> str:
    """The name setup_wgs.sh wrote next to a database (db_name.txt), else the file/folder name."""
    p = Path(path)
    for f in (p / "db_name.txt", p.parent / "db_name.txt"):
        try:
            t = f.read_text().strip().splitlines()[0].strip()
            if t:
                return t
        except (OSError, IndexError):
            pass
    return p.name


def _ver(t, args=("--version",)) -> str:
    try:
        import subprocess
        p = subprocess.run([str(t.path), *args], capture_output=True, text=True, timeout=60, env=t.env())
        txt = (p.stdout or p.stderr or "").strip().splitlines()
        m = re.search(r"\d+\.\d+(\.\d+)*", " ".join(txt[:3]))
        return m.group(0) if m else ""
    except Exception:
        return ""


def run_sample(reads: Path, out: Path, dbs: dict, threads: int = 8, mean_len: float = 5000,
               find=None, log=None) -> dict:
    """Run every available tool on one sample's reads; returns (and saves) the info dict."""
    find = find or W.find_tool
    log = log or W.log
    out.mkdir(parents=True, exist_ok=True)
    for f in OUTPUTS:                     # a tool skipped this time must not leave its old result
        (out / f).unlink(missing_ok=True)
    logs = out / "logs"
    tools, profiles, labels = [], {}, {}

    # Kraken2 + Bracken
    k2, br = find("kraken2"), find("bracken")
    kdb = Path(dbs.get("kraken2_db") or "")
    ent = {"key": "kraken2", "label": "Kraken2 + Bracken", "db": _db_label(str(kdb)) if str(kdb) not in ("", ".") else ""}
    if not k2:
        ent.update(status="skipped", note="Kraken2 is not installed (setup_wgs.sh --readtax)")
    elif not (kdb / "hash.k2d").exists():
        ent.update(status="skipped", note="Kraken2 database not configured (db_paths.json kraken2_db)")
    else:
        rep = out / "kraken2.report"
        rc = W.run(k2, ["--db", kdb, "--threads", threads, "--report", rep, "--output", os.devnull, reads],
                   logs / "kraken2.log")
        if rc != 0 or not rep.exists():
            ent.update(status="failed", note="Kraken2 failed — see readtax/logs/kraken2.log")
        else:
            ent.update(status="ok", version=_ver(k2), classified_pct=kraken_classified_pct(rep))
            prof, refs = {}, ["kraken2"]
            rl = pick_read_len(kmer_lengths(kdb), mean_len) if br else None
            if rl:
                bo = out / "bracken_species.tsv"
                if W.run(br, ["-d", kdb, "-i", rep, "-o", bo, "-w", out / "bracken_species.report",
                              "-r", rl, "-l", "S", "-t", 10], logs / "bracken.log") == 0:
                    prof = parse_bracken(bo)
                    if prof:
                        ent["bracken_read_len"], refs = rl, ["kraken2", "bracken"]
            if not prof:
                prof = parse_kraken_species(rep)
                ent["label"] = "Kraken2"
                ent["note"] = ("Bracken not installed" if not br else
                               "no Bracken k-mer distribution / Bracken failed") + " — Kraken2 species reads used"
            ent["refs"] = refs
            labels["kraken2"] = "Kraken2+Bracken" if "bracken" in refs else "Kraken2"
            profiles["kraken2"] = prof
    tools.append(ent)

    # Centrifuger (preferred) or Centrifuge
    idx = centrifuge_index(dbs.get("readtax_centrifuge") or "")
    ent = {"key": "centrifuge", "label": "Centrifuger" if not idx or idx[0] == "centrifuger" else "Centrifuge",
           "db": _db_label(idx[1]) if idx else ""}
    if not idx:
        ent.update(status="skipped", note="Centrifuger database not installed (db_paths.json readtax_centrifuge)")
    elif idx[0] == "centrifuger":
        cf, cq = find("centrifuger"), find("centrifuger-quant")
        if not cf or not cq:
            ent.update(status="skipped", note="Centrifuger is not installed (setup_wgs.sh --readtax)")
        else:
            cls, rep = out / "centrifuger_classification.tsv", out / "centrifuger_report.tsv"
            ok = (W.run(cf, ["-x", idx[1], "-u", reads, "-t", threads], logs / "centrifuger.log",
                        stdout_to=cls) == 0
                  and W.run(cq, ["-x", idx[1], "-c", cls], logs / "centrifuger.log", stdout_to=rep) == 0)
            cls.unlink(missing_ok=True)          # one line per read — large, not a deliverable
            prof = parse_centrifuge(rep) if ok else {}
            if not ok:
                ent.update(status="failed", note="Centrifuger failed — see readtax/logs/centrifuger.log")
            else:
                ent.update(status="ok", version=_ver(cf, ("-v",)), refs=["centrifuger"])
                profiles["centrifuge"] = prof
    else:
        c1 = find("centrifuge")
        if not c1:
            ent.update(status="skipped", note="Centrifuge is not installed (setup_wgs.sh --readtax)")
        else:
            rep = out / "centrifuge_report.tsv"
            rc = W.run(c1, ["-x", idx[1], "-U", reads, "-p", threads, "--report-file", rep, "-S", os.devnull],
                       logs / "centrifuge.log")
            if rc != 0 or not rep.exists():
                ent.update(status="failed", note="Centrifuge failed — see readtax/logs/centrifuge.log")
            else:
                ent.update(status="ok", version=_ver(c1), refs=["centrifuge"])
                profiles["centrifuge"] = parse_centrifuge(rep)
    labels["centrifuge"] = ent["label"]
    tools.append(ent)

    # GOTTCHA2, nanopore mode; fast-profile when the database is the fast bundle
    g2, gdb = find("gottcha2"), gottcha_db(dbs.get("readtax_gottcha2") or "")
    ent = {"key": "gottcha2", "label": "GOTTCHA2", "db": _db_label(gdb) if gdb else ""}
    if not gdb:
        ent.update(status="skipped", note="GOTTCHA2 database not installed (db_paths.json readtax_gottcha2)")
    elif not g2:
        ent.update(status="skipped", note="GOTTCHA2 is not installed (setup_wgs.sh --readtax)")
    else:
        fast = Path(gdb + ".syldb").exists() and not Path(gdb + ".mmi").exists()
        mode = "fast-profile" if fast else "profile"
        rc = W.run(g2, [mode, "-d", gdb, "-i", reads, "-np", "-t", threads, "-o", out, "-p", "gottcha2"],
                   logs / "gottcha2.log")
        rep = out / "gottcha2.tsv"
        for b in list(out.glob("gottcha2*.bam")) + list(out.glob("gottcha2*.bai")):
            b.unlink(missing_ok=True)            # signature alignments: large, not a deliverable
        if rc != 0 or not rep.exists():
            ent.update(status="failed", note="GOTTCHA2 failed — see readtax/logs/gottcha2.log")
        else:
            ent.update(status="ok", version=_ver(g2), mode=mode, refs=["gottcha"])
            profiles["gottcha2"] = parse_gottcha2(rep)
    tools.append(ent)

    for t in tools:
        if t["status"] == "ok" and not profiles.get(t["key"]):
            t["note"] = (t.get("note", "") + "; " if t.get("note") else "") + "no species-level assignment"
    if any(profiles.values()):
        tsv = write_table(profiles, out / TABLE, labels)
        (out / "heatmap_read_taxonomy.svg").write_text(svg_heatmap(tsv))
    else:
        (out / TABLE).unlink(missing_ok=True)
    info = {"reads": reads.name, "tools": tools,
            "tools_ref": sorted({r for t in tools if t["status"] == "ok" for r in t.get("refs", [])})}
    (out / INFO).write_text(json.dumps(info, indent=1))
    for t in tools:
        log(f"  read taxonomy · {t['label']}: {t['status']}" + (f" — {t['note']}" if t.get("note") else ""))
    return info


def load_info(d: Path) -> dict:
    try:
        return json.loads((Path(d) / INFO).read_text())
    except (OSError, ValueError):
        return {"tools": [], "tools_ref": []}


def method_sentence(info: dict) -> str:
    ran = [t for t in info.get("tools", []) if t.get("status") == "ok"]
    rest = [t for t in info.get("tools", []) if t.get("status") != "ok"]
    parts = []
    for t in ran:
        s = t["label"]
        extra = []
        if t.get("key") == "gottcha2":
            extra.append(f"nanopore mode, {t.get('mode', 'profile')}")
        if t.get("bracken_read_len"):
            extra.append(f"Bracken read length {t['bracken_read_len']}")
        if t.get("db"):
            extra.append(f"database {t['db']}")
        parts.append(s + (f" ({'; '.join(extra)})" if extra else ""))
    txt = ("The quality-filtered reads were classified independently with " +
           (", ".join(parts[:-1]) + " and " + parts[-1] if len(parts) > 1 else (parts[0] if parts else "no tool")) +
           ". Each column is the share of that tool's species-level assignments (%); "
           "tools use different databases and models, so agreement between columns supports the call.")
    for t in rest:
        txt += f" {t['label']} was not run ({t.get('note') or t.get('status')})."
    return txt


# ── CLI: one job (new or finished) ────────────────────────────────────────────
def main(argv=None, find=None) -> int:
    ap = argparse.ArgumentParser(description="Read-based taxonomy for an ONT WGS job")
    ap.add_argument("--results", required=True)
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--db_paths", default="")
    ap.add_argument("--samples", default="", help="comma-separated; default: every sample")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    res, work = Path(a.results), Path(a.work_dir)
    dbs = {}
    if a.db_paths and Path(a.db_paths).exists():
        try:
            dbs = json.loads(Path(a.db_paths).read_text())
        except ValueError:
            W.warn("db_paths.json unreadable")
    want = [s.strip() for s in a.samples.split(",") if s.strip()]
    samples = sorted(p.name for p in (res / "samples").iterdir() if p.is_dir()) if (res / "samples").is_dir() else []
    if want:
        samples = [s for s in samples if s in want]
    done = 0
    for i, S in enumerate(samples):
        W.progress(100 * i / max(1, len(samples)), f"[{S}] Read-based taxonomy")
        out = res / "samples" / S / "readtax"
        reads = work / S / "reads.filtered.fastq.gz"
        if not reads.exists():
            if (out / TABLE).exists():              # an earlier result stays as it is
                W.warn(f"[{S}] filtered reads not found — the earlier read-based taxonomy is kept")
                continue
            out.mkdir(parents=True, exist_ok=True)
            (out / INFO).write_text(json.dumps({"tools": [], "tools_ref": [],
                                                "note": "filtered reads not found (the upload folder was "
                                                        "cleaned) — re-run the job to add read-based taxonomy"}))
            W.warn(f"[{S}] filtered reads not found in {work / S} — read-based taxonomy skipped")
            continue
        mean_len = 5000.0
        try:
            r = json.loads((res / "samples" / S / "result.json").read_text())
            mean_len = float((r.get("reads_filtered") or {}).get("mean_len") or mean_len)
        except (OSError, ValueError, TypeError):
            pass
        info = run_sample(reads, out, dbs, a.threads, mean_len, find=find)
        done += any(t["status"] == "ok" for t in info["tools"])
    W.progress(100, "Read-based taxonomy done")
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
