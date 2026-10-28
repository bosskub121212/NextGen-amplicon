#!/usr/bin/env python3
"""shotgun_func.py — functional profiles of shotgun reads with HUMAnN 3.9 (opt-in).

Called by shotgun_pipeline.py (--func) after the resistome, before MAGs:

  reads      QC-passed, host-removed reads; Illumina R1 + R2 in one file (HUMAnN
             treats mates as single reads); Nanopore reads cut into 300 bp pieces
             (seqkit sliding), pieces < 100 bp dropped — HUMAnN is built for short reads
  taxonomy   MetaPhlAn 4.1.1, database mpa_vJun23_CHOCOPhlAnSGB_202307 (the only
             pairing HUMAnN 3.9 accepts)
  function   HUMAnN 3.9: ChocoPhlAn pangenomes (nucleotide), then full UniRef90
             (DIAMOND, translated) → gene families, MetaCyc pathways and coverage
  tables     joined over samples, CPM, regrouped to KEGG KO and EC, split into
             community totals and per-taxon contributions; gene panels
             (data/function_panels.tsv) matched on KO first, UniRef90 names second

Everything goes to <out>/func/. A failure never fails the run: func/SKIPPED.txt says why.
Databases: db_paths.json key humann_db (setup_shotgun.sh --func).
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENV = "ngamp-func"
MPA_INDEX = "mpa_vJun23_CHOCOPhlAnSGB_202307"
MAPS = {"ko": "map_ko_uniref90.txt.gz", "ec": "map_level4ec_uniref90.txt.gz"}
NAMES = {"ko": "map_ko_name.txt.gz", "ec": "map_ec_name.txt.gz", "genefamilies": "map_uniref90_name.txt.bz2"}
ONT_PIECE, ONT_MIN = 300, 100
SPECIAL = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED"}
_SUFFIX = re.compile(r"(_Abundance(-RPKs|-CPM)?|_Coverage|-RPKs|-CPM)$")


def _bin(name: str, conda_bases: list[Path]) -> str | None:
    cands = [Path(b) / "envs" / ENV / "bin" / name for b in conda_bases]
    w = shutil.which(name)
    if w:
        cands.append(Path(w))
    for p in cands:
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    return None


def _parts(dbp: dict) -> tuple[Path | None, dict, list[str]]:
    """(db_dir, found parts, missing descriptions)."""
    d = dbp.get("humann_db") or ""
    if not d:
        return None, {}, ["HUMAnN databases are not installed (run setup_shotgun.sh --func)"]
    db = Path(d)
    if not db.is_dir():
        return db, {}, [f"HUMAnN databases not found at {db} (drive not mounted, or run "
                        "setup_shotgun.sh --func)"]
    found, miss = {}, []
    done = lambda d: (d / ".complete").exists()      # written by setup_shotgun.sh after each download
    ch = db / "chocophlan"
    if ch.is_dir() and done(ch):
        found["chocophlan"] = str(ch)
    else:
        miss.append(f"ChocoPhlAn ({ch})")
    un = db / "uniref"
    if un.is_dir() and any(un.glob("*.dmnd")) and done(un):
        found["uniref"] = str(un)
    else:
        miss.append(f"UniRef90 ({un})")
    um = db / "utility_mapping"
    if all((um / f).exists() for f in MAPS.values()) and done(um):
        found["mapping"] = str(um)
    else:
        miss.append(f"utility mapping ({um})")
    mp = db / "metaphlan"
    if any(mp.glob(f"{MPA_INDEX}*.bt2l")) and (mp / f"{MPA_INDEX}.pkl").exists() and done(mp):
        found["metaphlan_db"] = str(mp)
    else:
        miss.append(f"MetaPhlAn {MPA_INDEX} ({mp})")
    return db, found, miss


def find_install(dbp: dict, conda_bases: list[Path]) -> dict | None:
    """Paths of a complete HUMAnN install, or None when any part is missing."""
    db, found, miss = _parts(dbp)
    h, m = _bin("humann", conda_bases), _bin("metaphlan", conda_bases)
    if miss or not h or not m:
        return None
    env = os.environ.copy()
    env["PATH"] = f"{Path(h).parent}{os.pathsep}{env.get('PATH', '')}"
    return {"db_dir": str(db), **found, "metaphlan_index": MPA_INDEX,
            "humann": h, "metaphlan": m, "bin": str(Path(h).parent), "env": env}


def missing_reason(dbp: dict, conda_bases: list[Path]) -> str:
    """Why function profiles cannot run, as one sentence; "" when they can."""
    db, found, miss = _parts(dbp)
    if not dbp.get("humann_db") or (db is not None and not db.is_dir()):
        return miss[0]
    out = []
    if miss:
        out.append("missing " + ", ".join(miss))
    if not _bin("humann", conda_bases) or not _bin("metaphlan", conda_bases):
        out.append(f"humann / metaphlan not installed (env {ENV})")
    return ("HUMAnN cannot run: " + "; ".join(out) + " — setup_shotgun.sh --func") if out else ""


# ── gene panels ──────────────────────────────────────────────────────────────
PANELS_FILE = HERE / "data" / "function_panels.tsv"
PANEL_SETS = ("gut", "food", "environment")


def load_panels(path: Path, enabled: list[str]) -> list[dict]:
    """Rows of the panel file for the enabled sets; ko is a list, uniref_pattern a string."""
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("panel") or "").strip() not in enabled:
                continue
            pat = (r.get("uniref_pattern") or "").strip()
            if pat:
                try:
                    re.compile(pat)
                except re.error as e:
                    import sys
                    print(f"[warn] function_panels.tsv: gene {r.get('gene')}: bad pattern {pat!r} ({e}) — row skipped",
                          file=sys.stderr)
                    continue
            out.append({"panel": r["panel"].strip(), "group": (r.get("group") or "").strip(),
                        "gene": r["gene"].strip(),
                        "ko": [k.strip() for k in (r.get("ko") or "").split(",") if k.strip()],
                        "uniref_pattern": (r.get("uniref_pattern") or "").strip(),
                        "note": (r.get("note") or "").strip()})
    return out


def genus_of(taxon: str) -> str:
    """'g__Faecalibacterium.s__F_prausnitzii' → 'Faecalibacterium'; '' stays '' (community)."""
    if not taxon:
        return ""
    m = re.search(r"g__([^.|]+)", taxon)
    if m:
        return m.group(1).replace("_", " ")
    return "Unclassified" if taxon.lower().startswith("unclassified") else taxon


def match_panels(ko_strat: dict, uniref_strat: dict, panels: list[dict],
                 samples: list[str] | None = None) -> tuple[list[dict], list[dict]]:
    """Panel genes per sample (community CPM) and per genus.

    ko_strat / uniref_strat: {feature: {taxon or "" (community): {sample: cpm}}}; KO
    features may carry a name ("K00929: buk"). A gene with KO ids is matched on them
    only; a gene without is matched on its pattern against UniRef90 names.
    """
    ko_by_id = defaultdict(list)
    for feat, d in ko_strat.items():
        ko_by_id[feat.split(":")[0].strip()].append(d)
    if samples is None:
        samples = sorted({s for src in (ko_strat, uniref_strat) for d in src.values()
                          for v in d.values() for s in v})
    rows, taxa = [], []
    for p in panels:
        if p["ko"]:
            hits = [d for k in p["ko"] for d in ko_by_id.get(k, [])]
        elif p["uniref_pattern"]:
            rx = re.compile(p["uniref_pattern"], re.I)
            hits = [d for feat, d in uniref_strat.items() if rx.search(feat)]
        else:
            hits = []
        tot = defaultdict(float)
        by_genus = defaultdict(float)
        for d in hits:
            for taxon, vals in d.items():
                for s, v in vals.items():
                    if taxon == "":
                        tot[s] += v
                    else:
                        by_genus[(s, genus_of(taxon))] += v
        for s in samples:
            rows.append({"sample": s, "panel": p["panel"], "group": p["group"], "gene": p["gene"],
                         "cpm": round(tot.get(s, 0.0), 6)})
        for (s, g), v in sorted(by_genus.items()):
            if v > 0:
                taxa.append({"sample": s, "panel": p["panel"], "gene": p["gene"], "taxon": g,
                             "cpm": round(v, 6)})
    return rows, taxa


# ── HUMAnN tables ────────────────────────────────────────────────────────────
def _rows(path: Path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        cols = [_SUFFIX.sub("", c) for c in head[1:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            yield cols, f[0], f[1:]


def read_humann_table(path: Path) -> dict[str, dict[str, float]]:
    """{feature (stratified rows as their own feature): {sample: value}}; sample names
    without HUMAnN's _Abundance / -RPKs / -CPM / _Coverage suffixes."""
    out = {}
    for cols, feat, vals in _rows(path):
        out[feat] = {c: float(v or 0) for c, v in zip(cols, vals)}
    return out


def read_strat(path: Path, keep=None) -> dict[str, dict[str, dict[str, float]]]:
    """{feature: {taxon or "" (community total): {sample: value}}}; keep(feature) filters
    while reading, so a multi-GB gene-family table never sits in memory whole."""
    out = defaultdict(dict)
    for cols, feat, vals in _rows(path):
        f, _, taxon = feat.partition("|")
        if keep is not None and not keep(f):
            continue
        out[f][taxon] = {c: float(v or 0) for c, v in zip(cols, vals)}
    return dict(out)


def parse_humann_log(path: Path) -> dict:
    """Read shares from a HUMAnN log, as % of all input reads."""
    try:
        txt = Path(path).read_text(errors="replace")
    except OSError:
        return {}
    nuc = re.findall(r"Unaligned reads after nucleotide alignment:\s*([\d.]+)\s*%", txt)
    tra = re.findall(r"Unaligned reads after translated alignment:\s*([\d.]+)\s*%", txt)
    out = {}
    if nuc:
        un1 = float(nuc[-1])
        out["nucleotide_pct"] = round(100 - un1, 4)
        if tra:
            un2 = float(tra[-1])
            out["translated_pct"] = round(un1 - un2, 4)
            out["unmapped_pct"] = round(un2, 4)
        else:                                        # translated search bypassed
            out["translated_pct"] = 0.0
            out["unmapped_pct"] = round(un1, 4)
    return out


# ── reads ────────────────────────────────────────────────────────────────────
def _open(p: Path, mode="rt"):
    return gzip.open(p, mode) if str(p).endswith(".gz") else open(p, mode)


def _n_reads(p: Path) -> int:
    with _open(p) as fh:
        return sum(1 for _ in fh) // 4


def _chunk_py(src: Path, dest: Path):
    """seqkit sliding + seq -m, in Python (when seqkit is missing)."""
    with _open(src) as fi, gzip.open(dest, "wt") as fo:
        while True:
            h = fi.readline()
            if not h:
                break
            sq, _, ql = fi.readline().strip(), fi.readline(), fi.readline().strip()
            name = h[1:].split()[0]
            for i in range(0, len(sq), ONT_PIECE):
                piece = sq[i:i + ONT_PIECE]
                if len(piece) >= ONT_MIN:
                    fo.write(f"@{name}_sliding:{i + 1}-{i + len(piece)}\n{piece}\n+\n{ql[i:i + ONT_PIECE]}\n")


def prep_reads(q: dict, platform: str, dest_dir: Path, seqkit, run, logf: Path) -> Path | None:
    """One FASTQ for HUMAnN, or None when there are no reads.
    Illumina: R1 + R2 one after the other. Nanopore: 300 bp pieces, pieces < 100 bp dropped."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    srcs = [Path(x) for x in (q.get("r1"), q.get("r2")) if x]
    srcs = [x for x in srcs if x.exists() and x.stat().st_size > 0]
    if not srcs:
        return None
    if platform == "ont":
        dest = dest_dir / "reads.fastq.gz"
        if seqkit:
            tmp = dest_dir / "pieces.fastq.gz"
            rc = run(seqkit, ["sliding", "-g", "-W", ONT_PIECE, "-s", ONT_PIECE, srcs[0], "-o", tmp], logf)
            rc = rc or run(seqkit, ["seq", "-m", ONT_MIN, tmp, "-o", dest], logf)
            tmp.unlink(missing_ok=True)
            if rc:
                dest.unlink(missing_ok=True)
                raise RuntimeError(f"seqkit failed (exit {rc}) cutting the reads — see {logf}")
        else:
            _chunk_py(srcs[0], dest)
    else:
        gz = [str(x).endswith(".gz") for x in srcs]
        if all(gz) or not any(gz):                   # concatenated gzip members are valid gzip
            dest = dest_dir / ("reads.fastq.gz" if all(gz) else "reads.fastq")
            with open(dest, "wb") as fo:
                for x in srcs:
                    with open(x, "rb") as fi:
                        shutil.copyfileobj(fi, fo, 1 << 20)
        else:
            dest = dest_dir / "reads.fastq.gz"
            with gzip.open(dest, "wb") as fo:
                for x in srcs:
                    with _open(x, "rb") as fi:
                        shutil.copyfileobj(fi, fo, 1 << 20)
    if not dest.exists() or dest.stat().st_size == 0:
        return None
    with _open(dest) as fh:
        if not fh.readline().startswith("@"):
            return None
    return dest


# ── run ──────────────────────────────────────────────────────────────────────
KINDS = ("genefamilies", "pathabundance", "pathcoverage")
QC_COLS = ["sample", "reads_in", "nucleotide_pct", "translated_pct", "unmapped_pct",
           "unintegrated_pct", "max_rss_gb", "status"]


def _version(t) -> str:
    import subprocess
    try:
        p = subprocess.run([t[0], "--version"], capture_output=True, text=True, env=t[1], timeout=120)
        txt = (p.stdout + p.stderr).strip()
        m = re.search(r"(\d+\.\d+(?:\.\d+)?)", txt)
        return m.group(1) if m else (txt.splitlines()[-1] if txt else "?")
    except Exception:
        return "?"


def _rss_gb(path: Path) -> float | None:
    try:
        m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", path.read_text())
        return round(int(m.group(1)) / 1024 / 1024, 2) if m else None
    except OSError:
        return None


def _skip(fd: Path, reason: str, log) -> dict:
    fd.mkdir(parents=True, exist_ok=True)
    (fd / "SKIPPED.txt").write_text(reason + "\n")
    log(f"[warn] function profiles skipped: {reason}")
    return {"skipped": reason}


def _write_csv(path: Path, cols: list[str], rows: list[dict]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _one_sample(s, q, platform, fd, threads, inst, T, seqkit, run, log, vers) -> dict:
    sd = fd / s
    sd.mkdir(parents=True, exist_ok=True)
    logf = sd / "run.log"
    srcs = [Path(x) for x in (q.get("r1"), q.get("r2")) if x]
    sig = json.dumps({"reads": [[str(x), x.stat().st_size, x.stat().st_mtime_ns] for x in srcs if x.exists()],
                      "versions": vers, "db": inst["db_dir"], "platform": platform}, sort_keys=True)
    outs = [sd / f"{s}_{k}.tsv" for k in KINDS]
    sigf, qcf = sd / "inputs.sig", sd / "qc.json"
    if sigf.exists() and sigf.read_text() == sig and qcf.exists() and all(o.exists() for o in outs):
        log(f"[{s}] function: already done with these reads and databases — reused")
        return json.loads(qcf.read_text())
    for o in outs + [sigf, qcf]:
        o.unlink(missing_ok=True)
    row = {"sample": s, "status": "no_reads"}
    reads = prep_reads(q, platform, sd, seqkit, run, logf)
    if reads is None:
        log(f"[warn] [{s}] function: no reads to profile")
        return row
    row["reads_in"] = _n_reads(reads)
    try:
        mpa = sd / f"{s}_metaphlan.tsv"
        bt2 = sd / f"{s}.bowtie2.bz2"
        bt2.unlink(missing_ok=True)
        rc = run(T("metaphlan"), [reads, "--input_type", "fastq", "--bowtie2db", inst["metaphlan_db"],
                                  "--index", inst["metaphlan_index"], "--nproc", threads,
                                  "--bowtie2out", bt2, "-o", mpa], logf)
        if rc != 0 or not mpa.exists():
            log(f"[warn] [{s}] MetaPhlAn failed (exit {rc}) — see func/{s}/run.log")
            row["status"] = "failed"
            return row
        args = ["--input", reads, "--output", sd, "--output-basename", s, "--threads", threads,
                "--taxonomic-profile", mpa, "--nucleotide-database", inst["chocophlan"],
                "--protein-database", inst["uniref"], "--remove-temp-output", "--o-log", sd / "humann.log"]
        tb = Path("/usr/bin/time")
        if tb.exists() and os.access(tb, os.X_OK):
            rc = run((str(tb), inst["env"]), ["-v", "-o", sd / "time.txt", inst["humann"]] + args, logf)
        else:
            rc = run(T("humann"), args, logf)
        if rc != 0 or not all(o.exists() for o in outs):
            log(f"[warn] [{s}] HUMAnN failed (exit {rc}) — see func/{s}/run.log and humann.log")
            row["status"] = "failed"
            return row
    finally:
        reads.unlink(missing_ok=True)              # a copy of the reads: tens of GB on big runs
        (sd / f"{s}.bowtie2.bz2").unlink(missing_ok=True)
        # --remove-temp-output cleans up only when HUMAnN finishes; after a crash or an
        # out-of-memory kill its alignments (tens of GB) would stay on the disk
        shutil.rmtree(sd / f"{s}_humann_temp", ignore_errors=True)
    row.update(parse_humann_log(sd / "humann.log"))
    row["max_rss_gb"] = _rss_gb(sd / "time.txt")
    row["status"] = "ok"
    qcf.write_text(json.dumps(row))
    sigf.write_text(sig)
    return row


def uniref_names(path: Path, rx, log) -> dict[str, str]:
    """{UniRef90 id: name} for the names matching rx, streamed from the (bz2) name map once
    and cached next to it per pattern set."""
    import bz2
    import hashlib
    if not path.exists():
        log(f"[warn] {path.name} not found — panel genes without a KO cannot be matched")
        return {}
    key = hashlib.md5(rx.pattern.encode()).hexdigest()[:12]
    cache = path.parent / f".panel_names_{key}.json"
    try:
        if cache.exists() and cache.stat().st_mtime >= path.stat().st_mtime:
            return json.loads(cache.read_text())
    except (OSError, ValueError):
        pass
    out = {}
    op = bz2.open if path.suffix == ".bz2" else (gzip.open if path.suffix == ".gz" else open)
    with op(path, "rt") as fh:
        for line in fh:
            uid, _, name = line.rstrip("\n").partition("\t")
            if name and rx.search(name):
                out[uid] = name
    try:
        cache.write_text(json.dumps(out))
    except OSError:
        pass
    return out


def _features(path: Path) -> int:
    return sum(1 for f in read_humann_table(path) if f.split(":")[0] not in SPECIAL) if path.exists() else 0


def run_all(samples, qc, platform, out: Path, threads, dbp, panels, tool, run, log, progress,
            conda_bases) -> dict:
    fd = Path(out) / "func"
    fd.mkdir(parents=True, exist_ok=True)
    (fd / "SKIPPED.txt").unlink(missing_ok=True)
    inst = find_install(dbp, list(conda_bases))
    if not inst:
        return _skip(fd, missing_reason(dbp, list(conda_bases)), log)

    def T(name):
        return (str(Path(inst["bin"]) / name), inst["env"])
    vers = {"humann": _version(T("humann")), "metaphlan": _version(T("metaphlan")),
            "metaphlan_db": inst["metaphlan_index"]}
    if (Path(inst["bin"]) / "diamond").exists():
        vers["diamond"] = _version(T("diamond"))
    log(f"function: HUMAnN {vers['humann']}, MetaPhlAn {vers['metaphlan']} ({inst['metaphlan_index']}), "
        f"databases {inst['db_dir']}")
    seqkit = tool("seqkit")
    rows, ok = [], []
    for i, s in enumerate(samples):
        progress(80 + int(12 * i / max(1, len(samples))), f"[{s}] function: MetaPhlAn + HUMAnN")
        try:
            r = _one_sample(s, qc.get(s) or {}, platform, fd, threads, inst, T, seqkit, run, log, vers)
        except Exception as e:                     # one bad sample must not cost the others
            import traceback
            log(f"[warn] [{s}] function step failed: {e}")
            (fd / s).mkdir(parents=True, exist_ok=True)
            with open(fd / s / "run.log", "a") as fh:
                fh.write(traceback.format_exc())
            shutil.rmtree(fd / s / f"{s}_humann_temp", ignore_errors=True)
            for f in ("reads.fastq.gz", "reads.fastq", "pieces.fastq.gz"):
                (fd / s / f).unlink(missing_ok=True)
            r = {"sample": s, "status": "failed"}
        rows.append(r)
        if r.get("status") == "ok":
            ok.append(s)
    failed = [r["sample"] for r in rows if r.get("status") != "ok"]
    if not ok:
        _write_csv(fd / "func_qc.csv", QC_COLS, rows)
        info = _skip(fd, "HUMAnN gave no result for any sample — see func/<sample>/run.log", log)
        info["failed"] = failed
        return info

    progress(92, "function: joining tables")
    logf = fd / "tables.log"
    jd = fd / "_join"
    shutil.rmtree(jd, ignore_errors=True)
    for k in KINDS:
        (jd / k).mkdir(parents=True)
        for s in ok:
            src = fd / s / f"{s}_{k}.tsv"
            try:
                (jd / k / src.name).symlink_to(src)
            except OSError:
                shutil.copy(src, jd / k / src.name)
        run(T("humann_join_tables"), ["-i", jd / k, "-o", fd / f"{k}_joined.tsv", "--file_name", k], logf)
    um = Path(inst["mapping"])
    run(T("humann_renorm_table"), ["-i", fd / "genefamilies_joined.tsv", "-o", fd / "genefamilies_cpm.tsv",
                                   "--units", "cpm", "--update-snames"], logf)
    run(T("humann_renorm_table"), ["-i", fd / "pathabundance_joined.tsv", "-o", fd / "pathways.tsv",
                                   "--units", "cpm", "--update-snames"], logf)
    shutil.move(str(fd / "pathcoverage_joined.tsv"), str(fd / "pathcoverage.tsv"))
    for k in ("ko", "ec"):
        run(T("humann_regroup_table"), ["-i", fd / "genefamilies_cpm.tsv", "-c", um / MAPS[k],
                                        "-o", fd / f"{k}_raw.tsv"], logf)
    for k, src in (("ko", "ko_raw.tsv"), ("ec", "ec_raw.tsv")):
        names = um / NAMES[k]
        if names.exists():
            run(T("humann_rename_table"), ["-i", fd / src, "-c", names, "-o", fd / f"{k}.tsv"], logf)
        if not (fd / f"{k}.tsv").exists():
            shutil.copy(fd / src, fd / f"{k}.tsv")
    # gene families stay UniRef90 ids (renaming the whole table with 87 M UniRef90 names
    # needs many GB of RAM and disk); panel genes without a KO are matched on names below
    shutil.move(str(fd / "genefamilies_cpm.tsv"), str(fd / "genefamilies.tsv"))
    for k in ("pathways", "pathcoverage", "ko", "ec"):
        run(T("humann_split_stratified_table"), ["-i", fd / f"{k}.tsv", "-o", fd], logf)
    for f in ("genefamilies_joined.tsv", "pathabundance_joined.tsv", "ko_raw.tsv", "ec_raw.tsv"):
        (fd / f).unlink(missing_ok=True)
    shutil.rmtree(jd, ignore_errors=True)

    # QC: share of pathway abundance HUMAnN could not place in a pathway
    pu = read_humann_table(fd / "pathways_unstratified.tsv") if (fd / "pathways_unstratified.tsv").exists() else {}
    for r in rows:
        s = r["sample"]
        tot = sum(v.get(s, 0) for v in pu.values())
        if r.get("status") == "ok" and tot:
            r["unintegrated_pct"] = round(100 * pu.get("UNINTEGRATED", {}).get(s, 0) / tot, 4)
    _write_csv(fd / "func_qc.csv", QC_COLS, rows)

    # gene panels
    pan = load_panels(PANELS_FILE, list(panels)) if panels else []
    kos = {k for p in pan for k in p["ko"]}
    pats = [p["uniref_pattern"] for p in pan if not p["ko"] and p["uniref_pattern"]]
    rx = re.compile("|".join(f"(?:{x})" for x in pats), re.I) if pats else None
    ko_s = read_strat(fd / "ko.tsv", keep=lambda f: f.split(":")[0].strip() in kos) if kos else {}
    un_s = {}
    if rx:
        names = uniref_names(um / NAMES["genefamilies"], rx, log)
        raw = read_strat(fd / "genefamilies.tsv", keep=lambda f: f in names)
        un_s = {f"{f}: {names[f]}": d for f, d in raw.items()}
    prow, ptax = match_panels(ko_s, un_s, pan, samples=ok)
    _write_csv(fd / "panels.csv", ["sample", "panel", "group", "gene", "cpm"], prow)
    _write_csv(fd / "panels_taxa.csv", ["sample", "panel", "gene", "taxon", "cpm"], ptax)

    import hashlib
    info = {"samples": ok, "failed": failed, "versions": vers, "db_dir": inst["db_dir"],
            "n_pathways": _features(fd / "pathways_unstratified.tsv"),
            "n_ko": _features(fd / "ko_unstratified.tsv"), "n_ec": _features(fd / "ec_unstratified.tsv"),
            "panels": list(panels), "panels_md5": hashlib.md5(PANELS_FILE.read_bytes()).hexdigest()}
    log(f"function: {info['n_pathways']} pathways, {info['n_ko']} KO, {info['n_ec']} EC across "
        f"{len(ok)} sample(s)" + (f"; failed: {', '.join(failed)}" if failed else ""))
    return info


# ── statistics and chart tables ──────────────────────────────────────────────
TOP_CHART = 500                     # features per chart-editor table (by mean CPM)


def shannon(vals) -> float:
    import math
    tot = sum(v for v in vals if v > 0)
    return -sum((v / tot) * math.log(v / tot) for v in vals if v > 0) if tot else 0.0


def bray_matrix(table: dict, samples: list[str]) -> list[list[float]]:
    """Bray-Curtis between samples on proportions of {feature: {sample: value}}."""
    from report_builder import _bray
    feats = list(table)
    tot = {s: sum(table[f].get(s, 0) for f in feats) for s in samples}
    vec = {s: [table[f].get(s, 0) / tot[s] if tot[s] else 0.0 for f in feats] for s in samples}
    n = len(samples)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            D[i][j] = D[j][i] = _bray(vec[samples[i]], vec[samples[j]])
    return D


def permanova(D, groups: list, perm: int = 999, seed: int = 1) -> dict:
    """One-way PERMANOVA (Anderson 2001) on a distance matrix."""
    import random
    keep = [i for i, g in enumerate(groups) if groups.count(g) >= 2]   # a group of one tells nothing
    D = [[D[i][j] for j in keep] for i in keep]
    groups = [groups[i] for i in keep]
    n = len(D)
    levels = sorted(set(groups))
    sizes = {g: groups.count(g) for g in levels}
    if len(levels) < 2 or min(sizes.values()) < 2:
        return {"skipped": "needs at least two groups with two or more samples each"}
    sst = sum(D[i][j] ** 2 for i in range(n) for j in range(i + 1, n)) / n

    def f_stat(g):
        ssw = 0.0
        for lv in levels:
            idx = [i for i in range(n) if g[i] == lv]
            ssw += sum(D[a][b] ** 2 for k, a in enumerate(idx) for b in idx[k + 1:]) / len(idx)
        ssa = sst - ssw
        a, w = len(levels) - 1, n - len(levels)
        return (ssa / a) / (ssw / w) if ssw > 0 else float("inf"), ssa
    f0, ssa = f_stat(groups)
    rng = random.Random(seed)
    g = list(groups)
    hits = 0
    for _ in range(perm):
        rng.shuffle(g)
        if f_stat(g)[0] >= f0 - 1e-12:
            hits += 1
    return {"F": round(f0, 4) if f0 != float("inf") else None, "R2": round(ssa / sst, 4) if sst else 0.0,
            "p": round((hits + 1) / (perm + 1), 4), "perm": perm, "groups": sizes}


def mpa_genera(path: Path) -> dict[str, float]:
    """Genus relative abundance (%) from a MetaPhlAn profile."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            clade = f[0].split("|")
            if len(f) < 3 or not clade[-1].startswith("g__"):
                continue
            g = clade[-1][3:].replace("_", " ")
            out[g] = out.get(g, 0.0) + float(f[2] or 0)
    return out


def kraken_genera(out: Path) -> dict[str, dict[str, float]]:
    """{sample: {genus: % of genus-assigned reads}} from asv_table.csv + taxonomy.csv."""
    tax = {}
    with open(out / "taxonomy.csv") as fh:
        for r in csv.DictReader(fh):
            tax[r.get("id") or next(iter(r.values()))] = (r.get("Genus") or "").strip()
    res = defaultdict(lambda: defaultdict(float))
    with open(out / "asv_table.csv") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        for r in rd:
            g = tax.get(r[0], "")
            if not g or g.lower() in ("unclassified", "na", "nan"):
                continue
            for s, v in zip(head[1:], r[1:]):
                res[s][g] += float(v or 0)
    return {s: {g: 100 * v / (sum(d.values()) or 1) for g, v in d.items()} for s, d in res.items()}


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(a, b) -> float | None:
    if len(a) < 3:
        return None
    ra, rb = _ranks(a), _ranks(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return round(num / den, 4) if den else None


def _wide(path: Path, label: str, table: dict, samples: list[str], top: int = TOP_CHART):
    feats = [f for f in table if f.split(":")[0].strip() not in SPECIAL]
    feats.sort(key=lambda f: -sum(table[f].get(s, 0) for s in samples))
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([label] + samples)
        for f in feats[:top]:
            w.writerow([f] + [round(table[f].get(s, 0), 4) for s in samples])


def _groups(metadata: str) -> dict[str, str]:
    if not metadata or not Path(metadata).exists():
        return {}
    out = {}
    with open(metadata, newline="") as fh:
        for r in csv.DictReader(fh):
            vals = list(r.values())
            sid = str(r.get("sampleId") or r.get("sample") or (vals[0] if vals else "")).strip()
            g = str(r.get("group") or r.get("Group") or "").strip()
            ctl = str(r.get("control") or r.get("sample_type") or r.get("type") or "").strip().lower()
            if ctl in ("negative", "blank", "ntc", "neg", "mock", "positive"):
                continue
            if sid and g:
                out[sid] = g
    return out


def func_stats(out: Path, metadata: str) -> dict:
    """Chart-editor tables, functional diversity, PERMANOVA and the MetaPhlAn vs Kraken 2 check.
    Writes r_tables/func_*.csv and func/mpa_vs_kraken.csv; returns what it found (also
    saved as func/func_stats.json)."""
    out = Path(out)
    fd, rt = out / "func", out / "r_tables"
    rt.mkdir(parents=True, exist_ok=True)
    pw = read_humann_table(fd / "pathways_unstratified.tsv")
    ko = read_humann_table(fd / "ko_unstratified.tsv") if (fd / "ko_unstratified.tsv").exists() else {}
    samples = list(next(iter(pw.values())).keys()) if pw else []
    res = {"samples": samples}
    _wide(rt / "func_pathways.csv", "Pathway", pw, samples)
    _wide(rt / "func_ko.csv", "KO", ko, samples)
    shutil.copy(fd / "panels_taxa.csv", rt / "func_panels_taxa.csv") if (fd / "panels_taxa.csv").exists() else None
    pan = defaultdict(dict)
    if (fd / "panels.csv").exists():
        with open(fd / "panels.csv") as fh:
            for r in csv.DictReader(fh):
                pan[f"{r['panel']}: {r['gene']}"][r["sample"]] = float(r["cpm"] or 0)
    with open(rt / "func_panels.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gene"] + samples)
        for g, d in pan.items():
            w.writerow([g] + [round(d.get(s, 0), 4) for s in samples])
    # alpha: Shannon over pathways (UNMAPPED / UNINTEGRATED left out)
    real = {f: v for f, v in pw.items() if f.split(":")[0].strip() not in SPECIAL}
    with open(rt / "func_alpha.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sample", "Shannon_pathways", "Pathways"])
        for s in samples:
            vals = [real[f].get(s, 0) for f in real]
            w.writerow([s, round(shannon(vals), 4), sum(1 for v in vals if v > 0)])
    # beta: Bray-Curtis on KO, PCoA and PERMANOVA
    kreal = {f: v for f, v in ko.items() if f.split(":")[0].strip() not in SPECIAL}
    if len(samples) >= 3 and kreal:
        from report_builder import _pcoa
        D = bray_matrix(kreal, samples)
        coords, pct = _pcoa(D)
        with open(rt / "func_pcoa_ko.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Sample", "PC1", "PC2"])
            for s, c in zip(samples, coords):
                w.writerow([s, round(c[0], 6), round(c[1], 6)])
        res["pcoa_pct"] = [round(x, 2) for x in pct]
        grp = _groups(metadata)
        S = [s for s in samples if s in grp]
        if S:
            idx = [samples.index(s) for s in S]
            res["permanova"] = permanova([[D[i][j] for j in idx] for i in idx], [grp[s] for s in S])
        else:
            res["permanova"] = {"skipped": "no sample groups in the metadata"}
    else:
        res["permanova"] = {"skipped": "fewer than three samples"}
    # MetaPhlAn 4 vs Kraken 2 + Bracken, genus level
    rows = []
    try:
        kg = kraken_genera(out)
    except OSError:
        kg = {}
    from report_builder import _bray
    for s in samples:
        mp = fd / s / f"{s}_metaphlan.tsv"
        if not mp.exists() or s not in kg:
            continue
        m = mpa_genera(mp)
        mt = sum(m.values()) or 1
        m = {g: 100 * v / mt for g, v in m.items()}
        k = kg[s]
        gen = sorted(set(m) | set(k))
        a, b = [k.get(g, 0.0) for g in gen], [m.get(g, 0.0) for g in gen]
        top = sorted(gen, key=lambda g: -(k.get(g, 0) + m.get(g, 0)))[:10]
        rows.append({"sample": s, "spearman": spearman(a, b),
                     "bray_curtis_similarity": max(0.0, round(1 - _bray([x / 100 for x in a], [x / 100 for x in b]), 4)) + 0.0,
                     "top_genera_json": json.dumps([[g, round(k.get(g, 0), 2), round(m.get(g, 0), 2)]
                                                    for g in top])})
    _write_csv(fd / "mpa_vs_kraken.csv", ["sample", "spearman", "bray_curtis_similarity", "top_genera_json"], rows)
    (fd / "func_stats.json").write_text(json.dumps(res, indent=1))
    return res
