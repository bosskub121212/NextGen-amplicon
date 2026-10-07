#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NextGen-Amplicon — ONT whole-genome pipeline for bacterial ISOLATES (ONT-WGS)

    python ont_wgs_pipeline.py --input <upload_dir> --output <results_dir>
                               --work_dir <scratch_dir> [--threads 8] [...]

What a lab hands over for an isolate run is a MinKNOW output folder: one folder per
sample (or per barcode), each holding dozens of chunk files that MinKNOW rotates
every hour. Nothing downstream can use chunks, so the first job here is the one a
person otherwise does by hand — put every chunk of a sample back into one file,
under the sample name the customer knows, and keep a record of which files went
where so a result can always be traced back to what was delivered.

Steps, per sample (each one resumes: a re-run skips whatever already finished
with the same settings, so a failure at hour four does not cost four hours):

   1  merge chunks             → <work>/<sample>/reads.fastq.gz  + sample_map.csv
   2  read QC + filter          seqkit  (length / mean-Q filter, GC per read)
   3  de novo assembly          Flye  (--nano-hq for R10 / SUP, --nano-raw for R9)
   4  polishing                 Medaka  (bacterial methylation model when it applies)
   5  assembly QC               contig table, circularity, depth, CheckM2
   6  species identification    sourmash gather vs GTDB  →  skani ANI vs the
                                closest reference genome (NCBI, when reachable)
   7  16S rRNA                  barrnap: every full-length copy → vsearch vs the
                                16S database already used by the amplicon pipelines
   8  MLST                      mlst (PubMLST schemes)
   9  AMR / stress / virulence  AMRFinderPlus (--plus, organism-aware point mutations)
  10  virulence factors         abricate + VFDB
  11  plasmids                  MOB-suite (mob_recon / mob_typer) + PlasmidFinder
  12  annotation                Bakta
  13  report                    wgs_report.py → wgs_report.html (+ PDF from the app)

Every tool after Flye is optional: when it is not installed, or its database is
missing, the step is skipped with a warning that also lands in the report, and the
rest of the run carries on. Flye is the only hard requirement.

Tools are looked up in conda environments first (ngamp-wgs, wgs, medaka, bakta, mobsuite,
checkm2 — see setup_wgs.sh), then on PATH. Each command runs with its own
environment's bin directory first on PATH, which is what activation would do,
so Perl/Python helpers inside a tool resolve to the right interpreter.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime
from pathlib import Path

__version__ = "1.0.0"

# ══════════════════════════════════════════════════════════════════════════════
#  Logging / progress (same wire format as the other pipelines: main.py parses
#  "PROGRESS:<pct>|<label>" lines into the job's step list)
# ══════════════════════════════════════════════════════════════════════════════
WARNINGS: list[str] = []


def progress(pct: float, label: str):
    print(f"PROGRESS:{int(max(0, min(100, pct)))}|{label}", flush=True)


def log(msg: str = ""):
    print(msg, flush=True)


def warn(msg: str):
    WARNINGS.append(msg)
    print(f"[WARN] {msg}", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
#  Grouping chunk files into samples
# ══════════════════════════════════════════════════════════════════════════════
FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
BARCODE_RE = re.compile(r"(?<![A-Za-z0-9])(barcode\d+|unclassified)(?![0-9])", re.I)
CHUNK_RE = re.compile(r"[_.-]\d+$")
FAIL_RE = re.compile(r"(^|[_\-.])fail([_\-.]|$)", re.I)
# MinKNOW chunk names: <flowcell>_<pass|fail>_[barcodeNN_]<run id>_..._<n>.fastq.gz
MINKNOW_RE = re.compile(r"^[A-Z]{3}\d{4,6}_(pass|fail|skip)_", re.I)
# Folder names that say nothing about the sample — never used as a sample name.
GENERIC_DIR_RE = re.compile(
    r"^(fastq(_pass|_fail)?|pass|fail|merged|data|reads|raw|raw_data|barcod(e|ing)\d*|"
    r"unclassified|output|outputs|results?|upload|uploads|\.|)$", re.I)


def natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def _clean_name(s: str) -> str:
    s = re.sub(r"^(samples?)[_\- ]+", "", s.strip(), flags=re.I) or s.strip()
    s = re.sub(r"[\\/:*?\"<>|\s]+", "_", s).strip("._")
    return s or "sample"


def group_fastqs(files: list[str], folders: dict | None = None,
                 overrides: dict | None = None) -> list[dict]:
    """Group FASTQ chunk files into samples.

    files     : file names as stored in the job's upload dir
    folders   : {file name: folder it came from} — from a ZIP or a folder import.
                A lab's folder ("Sample_OP-SCN5") is the name the customer knows,
                so it wins over the barcode whenever the folder holds one sample.
    overrides : {group key: sample name} chosen in the UI

    The grouping key is (folder, barcode) — or (folder, file stem without the
    chunk counter) when there is no barcode in the name. Keying on the folder too
    means two runs that both used barcode09 for different samples stay apart.
    Reads from MinKNOW's fail folders are left out: they are below the run's own
    quality threshold and would only lower the assembly.
    """
    folders = folders or {}
    overrides = overrides or {}
    raw: "OrderedDict[str, dict]" = OrderedDict()
    skipped_fail = 0
    for f in sorted(files, key=natural_key):
        base = Path(f).name
        if not FASTQ_RE.search(base):
            continue
        folder = folders.get(f, "") or ""
        if FAIL_RE.search(FASTQ_RE.sub("", base)) or any(
                FAIL_RE.search(p) for p in Path(folder).parts):
            skipped_fail += 1
            continue
        stem = FASTQ_RE.sub("", base)
        m = BARCODE_RE.search(stem)
        # Only MinKNOW's own chunk names carry a counter to strip; on any other
        # file "_1"/"_2" may well be a replicate and must not be merged away.
        minknow = bool(MINKNOW_RE.match(base))
        unit = m.group(1).lower() if m else (CHUNK_RE.sub("", stem) if minknow else stem)
        key = f"{folder}|{unit}"
        g = raw.setdefault(key, {"key": key, "folder": folder,
                                 "barcode": m.group(1).lower() if m else "",
                                 "unit": unit, "files": [], "minknow": minknow,
                                 "stem": stem})
        g["files"].append(f)

    per_folder = Counter(g["folder"] for g in raw.values())
    groups, used = [], set()
    for g in raw.values():
        dname = Path(g["folder"]).name if g["folder"] else ""
        if dname and per_folder[g["folder"]] == 1 and not GENERIC_DIR_RE.match(dname):
            name = _clean_name(dname)                  # Sample_OP-SCN5 → OP-SCN5
        elif not g["minknow"] and len(g["files"]) == 1:
            name = _clean_name(g["stem"])              # OP-SCN5_barcode09.fastq.gz
        else:
            name = _clean_name(g["unit"])              # barcode09
        if g["key"] in overrides and str(overrides[g["key"]]).strip():
            name = _clean_name(str(overrides[g["key"]]))
        base, i = name, 2
        while name in used:          # same barcode in two runs: tell them apart by run folder
            top = Path(g["folder"]).parts[0] if g["folder"] else ""
            name = f"{base}_{_clean_name(top)}" if top and i == 2 else f"{base}_{i}"
            i += 1
        used.add(name)
        g["sample"] = name
        groups.append(g)
    if skipped_fail:
        for g in groups:
            g["skipped_fail"] = skipped_fail
    return groups


def load_groups(input_dir: Path, overrides: dict) -> list[dict]:
    files = sorted([p.name for p in input_dir.iterdir()
                    if p.is_file() or p.is_symlink()], key=natural_key)
    files = [f for f in files if FASTQ_RE.search(f)]

    # Manual sample↔file table from the UI wins over everything else.
    man = input_dir / "sample_manifest.json"
    if man.exists():
        try:
            rows = json.loads(man.read_text())
            by: "OrderedDict[str, dict]" = OrderedDict()
            for r in rows:
                s, f1 = str(r.get("sample", "")).strip(), str(r.get("file1", "")).strip()
                if s and f1 and (input_dir / f1).exists():
                    g = by.setdefault(s, {"key": f"manifest|{s}", "folder": "",
                                          "barcode": "", "unit": s,
                                          "sample": _clean_name(s), "files": []})
                    g["files"].append(f1)
                    m = BARCODE_RE.search(f1)
                    if m and not g["barcode"]:
                        g["barcode"] = m.group(1).lower()
            if by:
                log(f"  Using sample_manifest.json ({len(by)} sample(s))")
                return list(by.values())
        except Exception as e:
            warn(f"sample_manifest.json unreadable ({e}) — grouping by file names")

    folders = {}
    fmap = input_dir / "upload_folders.json"
    if fmap.exists():
        try:
            folders = json.loads(fmap.read_text())
        except Exception:
            folders = {}
    return group_fastqs(files, folders, overrides)


def merge_group(input_dir: Path, g: dict, dest: Path) -> None:
    """Concatenate a sample's chunks into one gzip file.

    Concatenated gzip members are a valid gzip stream, so .gz chunks are copied
    byte for byte (no decompression); plain FASTQ chunks are compressed on the
    way in. A single .gz chunk is linked rather than copied.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    srcs = [input_dir / f for f in g["files"]]
    if len(srcs) == 1 and srcs[0].name.lower().endswith(".gz"):
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        try:
            dest.symlink_to(srcs[0].resolve())
            return
        except OSError:
            pass
    with open(tmp, "wb") as out:
        for s in srcs:
            if s.name.lower().endswith(".gz"):
                with open(s, "rb") as fh:
                    shutil.copyfileobj(fh, out, 8 << 20)
            else:
                with open(s, "rb") as fh, gzip.GzipFile(fileobj=out, mode="wb",
                                                        compresslevel=1) as gz:
                    shutil.copyfileobj(fh, gz, 8 << 20)
    tmp.replace(dest)


# ══════════════════════════════════════════════════════════════════════════════
#  Tool discovery
# ══════════════════════════════════════════════════════════════════════════════
TOOL_ENVS = {
    "flye": ["ngamp-wgs", "wgs", "flye"],
    "medaka_consensus": ["medaka"],
    "medaka": ["medaka"],
    "bakta": ["bakta"],
    "mob_recon": ["mobsuite", "mob_suite"],
    "checkm2": ["checkm2"],
}
DEFAULT_ENVS = ["ngamp-wgs", "wgs"]   # setup_wgs.sh builds ngamp-wgs; "wgs" = hand-made


def conda_base() -> Path | None:
    cands = []
    if os.environ.get("NGAMP_CONDA_BASE"):
        cands.append(Path(os.environ["NGAMP_CONDA_BASE"]))
    if os.environ.get("CONDA_EXE"):
        cands.append(Path(os.environ["CONDA_EXE"]).resolve().parent.parent)
    if os.environ.get("MAMBA_ROOT_PREFIX"):
        cands.append(Path(os.environ["MAMBA_ROOT_PREFIX"]))
    home = Path.home()
    cands += [home / d for d in ("miniconda3", "anaconda3", "miniforge3", "mambaforge",
                                 "micromamba")]
    cands += [Path("/opt/conda"), Path("/opt/miniconda3")]
    for c in cands:
        if (c / "envs").is_dir():
            return c
    return None


_BASE = conda_base()


class Tool:
    def __init__(self, name: str, path: Path, bindir: Path):
        self.name, self.path, self.bindir = name, path, bindir

    def env(self) -> dict:
        e = os.environ.copy()
        e["PATH"] = f"{self.bindir}{os.pathsep}{e.get('PATH', '')}"
        prefix = self.bindir.parent
        if (prefix / "conda-meta").is_dir():
            e["CONDA_PREFIX"] = str(prefix)
        # Several tools (medaka/pytorch, sourmash) otherwise grab every core.
        e.setdefault("OMP_NUM_THREADS", "4")
        return e


_TOOL_CACHE: dict[str, Tool | None] = {}


def find_tool(name: str) -> Tool | None:
    if name in _TOOL_CACHE:
        return _TOOL_CACHE[name]
    found = None
    if _BASE:
        for env in TOOL_ENVS.get(name, []) + DEFAULT_ENVS:
            p = _BASE / "envs" / env / "bin" / name
            if p.exists() and os.access(p, os.X_OK):
                found = Tool(name, p, p.parent)
                break
    if not found:
        w = shutil.which(name)
        if w:
            found = Tool(name, Path(w), Path(w).parent)
    _TOOL_CACHE[name] = found
    return found


class StepFailed(RuntimeError):
    pass


def run(tool: Tool, args: list, logfile: Path, cwd: Path | None = None,
        timeout: int | None = None, stdout_to: Path | None = None,
        quiet_ok: bool = False) -> int:
    """Run a tool with its environment, output to its own log file.

    The pipeline log only gets the command and, on failure, the tail of the
    tool's own output — Flye and Medaka print thousands of progress lines that
    would otherwise bury everything else in the job log.
    """
    cmd = [str(tool.path)] + [str(a) for a in args]
    log(f"    $ {tool.name} {' '.join(str(a) for a in args)}"[:600])
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with open(logfile, "a") as lf:
        lf.write(f"\n### {datetime.now().isoformat(timespec='seconds')}  {' '.join(cmd)}\n")
        lf.flush()
        out = open(stdout_to, "w") if stdout_to else lf
        try:
            p = subprocess.run(cmd, stdout=out, stderr=lf, cwd=cwd, env=tool.env(),
                               timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = -9
            lf.write(f"\n### TIMEOUT after {timeout}s\n")
        finally:
            if stdout_to:
                out.close()
    if rc != 0 and not quiet_ok:
        try:
            tail = logfile.read_text(errors="replace").strip().splitlines()[-12:]
        except Exception:
            tail = []
        log(f"    ✗ {tool.name} exited {rc} — last lines of {logfile.name}:")
        for t in tail:
            log(f"      | {t[:300]}")
    return rc


def tool_version(name: str) -> str:
    t = find_tool(name)
    if not t:
        return ""
    flags = {"medaka_consensus": None, "seqkit": ["version"]}
    f = flags.get(name, ["--version"])
    if f is None:
        return ""
    try:
        p = subprocess.run([str(t.path)] + f, capture_output=True, text=True,
                           timeout=60, env=t.env())
        lines = [ln.replace(str(t.path), "") for ln in
                 ((p.stderr or "") + "\n" + (p.stdout or "")).splitlines()]
        stem = name.split("_")[0].lower()
        # lines naming the tool first: vsearch prints a citation with a DOI
        # ("10.7717/…") on stdout before its version on stderr
        for line in sorted(lines, key=lambda ln: stem not in ln.lower()):
            m = re.search(r"(?<![\w./])v?(\d+\.\d+(?:\.\d+)?(?:-b\d+)?)(?![\d./])", line)
            if m:
                return m.group(1)
        return ""
    except Exception:
        return ""


# ══════════════════════════════════════════════════════════════════════════════
#  Small helpers
# ══════════════════════════════════════════════════════════════════════════════
def read_fasta(path: Path):
    name, seq = None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path: Path, records, width: int = 80):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def gc_pct(seq: str) -> float:
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return 100.0 * (s.count("G") + s.count("C")) / acgt if acgt else 0.0


def n50(lengths) -> int:
    ls = sorted(lengths, reverse=True)
    half, acc = sum(ls) / 2, 0
    for L in ls:
        acc += L
        if acc >= half:
            return L
    return 0


def read_tsv(path: Path, comment_header: bool = False) -> list[dict]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, newline="", errors="replace") as fh:
        lines = fh.read().splitlines()
    if not lines:
        return []
    head = lines[0]
    if comment_header and head.startswith("#"):
        head = head[1:]
    cols = head.split("\t")
    out = []
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            continue
        vals = line.split("\t")
        out.append({c: (vals[i] if i < len(vals) else "") for i, c in enumerate(cols)})
    return out


def write_csv(path: Path, rows: list[dict], cols: list[str] | None = None):
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def pick(row: dict, *names, default=""):
    """First present column among several spellings (tools rename columns)."""
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return default


def fnum(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def human_bp(n) -> str:
    n = float(n or 0)
    for unit, div in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if n >= div:
            return f"{n / div:.2f} {unit}"
    return f"{int(n)} bp"


def parse_size(s: str) -> int:
    m = re.match(r"^\s*([\d.]+)\s*([kmg]?)b?\s*$", str(s), re.I)
    if not m:
        return 0
    return int(float(m.group(1)) * {"": 1, "k": 1e3, "m": 1e6, "g": 1e9}[m.group(2).lower()])


# ── step bookkeeping: lets a re-run skip finished work ───────────────────────
def _sig(obj) -> str:
    return hashlib.sha1(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:12]


def is_done(sdir: Path, step: str, sig: str) -> bool:
    p = sdir / f".done_{step}"
    return p.exists() and p.read_text().strip() == sig


def mark_done(sdir: Path, step: str, sig: str):
    (sdir / f".done_{step}").write_text(sig)


def clear_done(sdir: Path, *steps):
    for s in steps:
        p = sdir / f".done_{s}"
        if p.exists():
            p.unlink()


# ══════════════════════════════════════════════════════════════════════════════
#  Read header → run metadata (dorado SAM-style tags or MinKNOW key=value)
# ══════════════════════════════════════════════════════════════════════════════
def read_run_info(fq: Path) -> dict:
    info = {}
    try:
        with gzip.open(fq, "rt", errors="replace") if str(fq).endswith(".gz") \
                else open(fq, errors="replace") as fh:
            header = fh.readline().strip()
    except Exception:
        return info
    for tok in re.split(r"[\t ]+", header[1:]):
        if "=" in tok:
            k, v = tok.split("=", 1)
            info[k] = v
        elif re.match(r"^[A-Za-z]{2}:[A-Za-z]:", tok):
            k, _, v = tok.split(":", 2)
            info[k] = v
    model = info.get("basecall_model_version_id", "")
    if not model and "RG" in info:
        m = re.search(r"(dna_r[\d.]+_e[\d.]+_\d+bps_\w+@v[\d.]+)", info["RG"])
        model = m.group(1) if m else ""
    out = {
        "basecall_model": model,
        "flow_cell_id": info.get("flow_cell_id") or info.get("PU", ""),
        "run_start": info.get("DT") or info.get("start_time", ""),
        "barcode_tag": info.get("barcode") or info.get("SM") or info.get("al", ""),
        "protocol_group_id": info.get("protocol_group_id", ""),
    }
    mm = re.search(r"r(\d+)\.(\d+)", model)
    out["chemistry"] = f"R{mm.group(1)}.{mm.group(2)}" if mm else ""
    out["accuracy"] = next((a for a in ("sup", "hac", "fast") if f"_{a}" in model), "")
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  Step 2 — read QC
# ══════════════════════════════════════════════════════════════════════════════
LEN_BINS = [0, 500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000, 30000, 50000,
            75000, 100000, 10 ** 9]


def read_profile(seqkit: Tool, fq: Path, tsv: Path, threads: int, logf: Path) -> dict:
    """Per-read length / GC / mean-Q in one seqkit pass, summarised here.

    The GC histogram is the quickest isolate check there is: one organism gives
    one peak; a mixed culture gives two or more, long before an assembly says so.
    """
    if not (tsv.exists() and tsv.stat().st_size > 0):
        rc = run(seqkit, ["fx2tab", "-n", "-i", "-l", "-g", "-q", "-H", "-j", threads,
                          fq], logf, stdout_to=tsv)
        if rc != 0:
            raise StepFailed("seqkit fx2tab failed")
    lens, gcs, qs = [], [], []
    with open(tsv) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        il = head.index("length") if "length" in head else 1
        ig = head.index("GC") if "GC" in head else 2
        iq = head.index("AvgQual") if "AvgQual" in head else 3
        for line in fh:
            v = line.rstrip("\n").split("\t")
            try:
                lens.append(int(v[il]))
                gcs.append(float(v[ig]))
                qs.append(float(v[iq]))
            except (ValueError, IndexError):
                continue
    if not lens:
        raise StepFailed(f"no reads in {fq.name}")
    total = sum(lens)
    srt = sorted(lens)
    len_hist = [0] * (len(LEN_BINS) - 1)
    len_hist_bp = [0] * (len(LEN_BINS) - 1)
    for L in lens:
        for i in range(len(LEN_BINS) - 1):
            if LEN_BINS[i] <= L < LEN_BINS[i + 1]:
                len_hist[i] += 1
                len_hist_bp[i] += L
                break
    gc_hist = [0.0] * 101                 # base-weighted, 1 % bins
    for L, g in zip(lens, gcs):
        gc_hist[min(100, max(0, int(round(g))))] += L
    q_hist = [0] * 51
    for q in qs:
        q_hist[min(50, max(0, int(q)))] += 1
    peaks = gc_peaks(gc_hist)
    return {
        "reads": len(lens), "bases": total,
        "mean_len": round(total / len(lens), 1), "median_len": srt[len(srt) // 2],
        "n50": n50(lens), "max_len": srt[-1],
        "mean_q": round(sum(qs) / len(qs), 2),
        "pct_q10": round(100 * sum(1 for q in qs if q >= 10) / len(qs), 2),
        "pct_q20": round(100 * sum(1 for q in qs if q >= 20) / len(qs), 2),
        "gc": round(sum(L * g for L, g in zip(lens, gcs)) / total, 2),
        "len_bins": LEN_BINS[:-1], "len_hist": len_hist, "len_hist_bp": len_hist_bp,
        "gc_hist": [round(x) for x in gc_hist], "q_hist": q_hist,
        "gc_peaks": peaks,
    }


def gc_peaks(hist: list[float]) -> list[dict]:
    """Peaks of a base-weighted GC histogram, smoothed over ±2 %.

    A secondary peak only counts when it holds ≥5 % of all bases and sits ≥5 GC
    points from the main one — read-level GC has a few points of scatter even in
    a pure culture, so anything tighter is noise, not a second organism.
    """
    sm = []
    for i in range(len(hist)):
        w = hist[max(0, i - 2): i + 3]
        sm.append(sum(w) / len(w))
    tot = sum(hist) or 1
    peaks = []
    for i in range(1, len(sm) - 1):
        if sm[i] >= sm[i - 1] and sm[i] > sm[i + 1]:
            mass = sum(hist[max(0, i - 4): i + 5]) / tot
            peaks.append({"gc": i, "mass": round(mass, 3)})
    peaks.sort(key=lambda p: -p["mass"])
    keep = []
    for p in peaks:
        if p["mass"] < 0.05:
            continue
        if all(abs(p["gc"] - k["gc"]) >= 5 for k in keep):
            keep.append(p)
    return keep[:4]


# ══════════════════════════════════════════════════════════════════════════════
#  Steps 3-4 — assembly + polishing
# ══════════════════════════════════════════════════════════════════════════════
def parse_flye_info(path: Path) -> dict:
    info = {}
    for r in read_tsv(path, comment_header=True):
        name = pick(r, "seq_name")
        if name:
            info[name] = {"cov": fnum(pick(r, "cov."), 0),
                          "circular": pick(r, "circ.") == "Y",
                          "repeat": pick(r, "repeat") == "Y",
                          "mult": pick(r, "mult.")}
    return info


def run_medaka(reads: Path, draft: Path, outdir: Path, threads: int, model: str,
               bacterial: bool, logf: Path) -> tuple[Path | None, str]:
    """Medaka with the bacterial methylation model when it applies.

    `--bacteria` exists from medaka 2.0 and is refused when the reads were called
    with a model it has no bacterial counterpart for, so on failure this falls
    back to the plain model, and then to no polishing at all — an unpolished Flye
    assembly is still a usable genome, a failed run is not.
    """
    med = find_tool("medaka_consensus")
    if not med:
        warn("medaka not installed — assembly NOT polished (Flye consensus only). "
             "Run setup_wgs.sh to add it.")
        return None, "not installed"
    attempts = []
    base = ["-i", reads, "-d", draft, "-o", outdir, "-t", min(threads, 8), "-f"]
    if model and model != "auto":
        base += ["-m", model]
    if bacterial:
        attempts.append(("bacterial methylation model", base + ["--bacteria"]))
    attempts.append(("default model" if not model or model == "auto" else model, base))
    for label, args in attempts:
        if outdir.exists():
            shutil.rmtree(outdir, ignore_errors=True)
        rc = run(med, args, logf, timeout=6 * 3600)
        cons = outdir / "consensus.fasta"
        if rc == 0 and cons.exists() and cons.stat().st_size > 0:
            return cons, label
        log(f"    medaka ({label}) did not finish — trying next option")
    warn("medaka failed with every model — assembly NOT polished (Flye consensus only)")
    return None, "failed"


def classify_contigs(contigs: list[dict]) -> None:
    """chromosome / plasmid / contig labels before MOB-suite has a say.

    A closed bacterial chromosome is the largest replicon and well over 1 Mb;
    a circular sequence below that is a plasmid (or phage) candidate.
    """
    if not contigs:
        return
    big = max(contigs, key=lambda c: c["length"])
    for c in contigs:
        if c is big and c["length"] >= 1_000_000:
            c["type"] = "chromosome"
        elif c["circular"] and c["length"] < 1_000_000:
            c["type"] = "plasmid?"
        else:
            c["type"] = "contig"


# ══════════════════════════════════════════════════════════════════════════════
#  Step 6 — species identification
# ══════════════════════════════════════════════════════════════════════════════
def _strip_rank(s: str) -> str:
    return re.sub(r"^[a-z]__", "", (s or "").strip())


def load_lineages(path: Path, idents: set[str]) -> dict:
    """GTDB lineage rows for just the accessions we need (file is ~100 MB)."""
    out = {}
    if not path or not path.exists() or not idents:
        return out
    want = {i.split(".")[0]: i for i in idents}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            ident = (r.get("ident") or r.get("accession") or "").strip()
            core = re.sub(r"^(RS_|GB_)", "", ident).split(".")[0]
            if core in want:
                out[want[core]] = {k: _strip_rank(v) for k, v in r.items()
                                   if k in ("superkingdom", "domain", "phylum", "class",
                                            "order", "family", "genus", "species")}
                if len(out) == len(want):
                    break
    return out


def species_call(ani: float | None) -> str:
    if ani is None:
        return "unresolved"
    if ani >= 95.0:
        return "species"
    if ani >= 90.0:
        return "genus (ANI 90–95 %: possible novel species)"
    return "unresolved (ANI < 90 %)"


def run_sourmash(asm: Path, sdir: Path, db: str, lineages: str, threads: int,
                 logf: Path) -> list[dict]:
    sm = find_tool("sourmash")
    if not sm:
        warn("sourmash not installed — genome-based species ID skipped")
        return []
    if not db or not Path(db).exists():
        warn("GTDB sourmash database not configured (db_paths.json key 'gtdb_sourmash') "
             "— genome-based species ID skipped; run setup_wgs.sh --gtdb")
        return []
    sig = sdir / "assembly.sig.zip"
    if run(sm, ["sketch", "dna", "-p", "k=31,scaled=1000,abund", "--name", "query",
                "-o", sig, asm], logf) != 0:
        return []
    gather = sdir / "sourmash_gather.csv"
    if run(sm, ["gather", sig, db, "-k", "31", "--threshold-bp", "50000",
                "-o", gather], logf, timeout=3 * 3600) != 0 or not gather.exists():
        warn("sourmash gather failed — genome-based species ID skipped")
        return []
    with open(gather, newline="") as fh:
        rows = list(csv.DictReader(fh))
    idents = set()
    for r in rows:
        acc = (r.get("name") or "").split(" ")[0]
        r["_acc"] = re.sub(r"^(RS_|GB_)", "", acc)
        idents.add(r["_acc"])
    lin = load_lineages(Path(lineages), idents) if lineages else {}
    hits = []
    for rank, r in enumerate(rows, 1):
        ani = fnum(pick(r, "match_containment_ani", "average_containment_ani"))
        ani = ani * 100 if ani is not None and ani <= 1.0 else ani
        L = lin.get(r["_acc"], {})
        nm = (r.get("name") or "")
        ncbi_name = nm.split(" ", 1)[1] if " " in nm else ""
        hits.append({
            "rank": rank,
            "accession": r["_acc"],
            "gtdb_species": L.get("species", ""),
            "gtdb_genus": L.get("genus", ""),
            "gtdb_family": L.get("family", ""),
            "gtdb_lineage": ";".join(L.get(k, "") for k in
                                     ("domain" if "domain" in L else "superkingdom",
                                      "phylum", "class", "order", "family", "genus",
                                      "species")) if L else "",
            "reference_name": ncbi_name,
            "ani": round(ani, 2) if ani is not None else "",
            "f_query": round(fnum(pick(r, "f_unique_weighted", "f_unique_to_query"), 0)
                             * 100, 2),
            "f_match": round(fnum(pick(r, "f_match"), 0) * 100, 2),
            "intersect_bp": int(fnum(pick(r, "intersect_bp"), 0)),
        })
    return hits


def fetch_ncbi_genome(acc: str, cache: Path, timeout: int = 180) -> Path | None:
    """Reference genome FASTA from NCBI Datasets, cached across runs."""
    acc = re.sub(r"^(RS_|GB_)", "", acc)
    if not re.match(r"^GC[AF]_\d+\.\d+$", acc):
        return None
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / f"{acc}.fna.gz"
    if out.exists() and out.stat().st_size > 1000:
        return out
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
           f"?include_annotation_type=GENOME_FASTA")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "NextGen-Amplicon"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            fna = [n for n in zf.namelist() if n.endswith((".fna", ".fasta", ".fa"))]
            if not fna:
                return None
            with zf.open(fna[0]) as src, gzip.open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
        return out
    except Exception as e:
        log(f"    reference download failed for {acc}: {e}")
        return None


def run_skani(asm: Path, ref: Path, logf: Path, outf: Path) -> dict:
    sk = find_tool("skani")
    if not sk or not ref:
        return {}
    if run(sk, ["dist", "-q", asm, "-r", ref, "-o", outf], logf) != 0:
        return {}
    rows = read_tsv(outf)
    if not rows:
        return {}
    r = rows[0]
    return {"ani": fnum(pick(r, "ANI")), "af_query": fnum(pick(r, "Align_fraction_query")),
            "af_ref": fnum(pick(r, "Align_fraction_ref"))}


# ══════════════════════════════════════════════════════════════════════════════
#  Step 7 — 16S rRNA
# ══════════════════════════════════════════════════════════════════════════════
def resolve_16s_db(choice: str, dbp: dict) -> tuple[Path | None, Path | None, str]:
    """(fasta, taxonomy dir or None, label). Accepts an Emu database directory
    (species_taxid.fasta or sequences.fasta + taxonomy.tsv) or a SILVA DADA2
    trainset FASTA whose headers ARE the lineage."""
    cands = [choice] if choice else []
    cands += [dbp.get(k, "") for k in ("wgs_16s_db", "emu_db_mar2026", "emu_silva",
                                       "SILVA_16S_sp", "SILVA_16S")]
    for c in cands:
        if not c:
            continue
        p = Path(c)
        if p.is_file() and (p.parent / "taxonomy.tsv").exists():
            p = p.parent
        if p.is_dir() and (p / "taxonomy.tsv").exists():
            for fa in ("species_taxid.fasta", "sequences.fasta"):
                if (p / fa).exists():
                    return p / fa, p, p.name
        elif p.is_file():
            return p, None, p.name
    return None, None, ""


def classify_16s(seqs: list[tuple[str, str]], db_fa: Path, tax_dir: Path | None,
                 threads: int, work: Path, logf: Path) -> list[dict]:
    vs = find_tool("vsearch")
    if not vs or not seqs:
        if not vs:
            warn("vsearch not installed — 16S copies extracted but not classified")
        return []
    q = work / "16S_query.fasta"
    write_fasta(q, seqs)
    hits = work / "16S_hits.b6"
    rc = run(vs, ["--usearch_global", q, "--db", db_fa, "--id", "0.80",
                  "--maxaccepts", "64", "--maxrejects", "512", "--strand", "both",
                  "--query_cov", "0.80", "--threads", threads, "--blast6out", hits,
                  "--output_no_hits", "--notrunclabels"], logf, timeout=3 * 3600)
    if rc != 0 or not hits.exists():
        warn("vsearch 16S search failed")
        return []
    rows = []
    with open(hits) as fh:
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) >= 3 and v[1] != "*":
                rows.append((v[0].split()[0], v[1], float(v[2])))
    # target label → lineage
    seq2tax, taxmap = {}, {}
    if tax_dir:
        s2t = tax_dir / "seq2taxid.tsv"
        rows = [(q, t.split()[0], i) for q, t, i in rows]
        targets = {t for _, t, _ in rows}
        if s2t.exists():
            with open(s2t) as fh:
                for line in fh:
                    v = line.rstrip("\n").split("\t")
                    if len(v) >= 2 and v[0] in targets:
                        seq2tax[v[0]] = v[1]
        for t in targets:
            seq2tax.setdefault(t, t.split(":")[0])
        need = set(seq2tax.values())
        with open(tax_dir / "taxonomy.tsv", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                tid = str(r.get("tax_id", "")).strip()
                if tid in need:
                    taxmap[tid] = r

    def lineage(target: str) -> tuple[str, str]:
        if tax_dir:
            r = taxmap.get(seq2tax.get(target, ""), {})
            sp = (r.get("species") or "").strip()
            ge = (r.get("genus") or "").strip()
            return sp, ge
        if ";" in target:                      # SILVA trainset: header = lineage
            parts = [p.strip() for p in target.strip(";").split(";") if p.strip()]
            if len(parts) >= 7:
                sp = parts[6].replace("_", " ")
                return (sp if sp.startswith(parts[5]) else f"{parts[5]} {sp}"), parts[5]
            return "", parts[5] if len(parts) >= 6 else (parts[-1] if parts else "")
        words = target.split()                 # "ACC Genus species ..." (assignSpecies)
        if len(words) >= 3:
            return f"{words[1]} {words[2]}", words[1]
        return "", words[1] if len(words) == 2 else ""

    by_q = defaultdict(list)
    for qn, t, ident in rows:
        by_q[qn].append((ident, t))
    out = []
    for qn, _ in seqs:
        hs = sorted(by_q.get(qn.split()[0], []), reverse=True)
        if not hs:
            out.append({"copy": qn.split()[0], "best_identity": "", "best_species": "",
                        "best_genus": "", "tied_species": "", "n_tied": 0})
            continue
        best = hs[0][0]
        tied = OrderedDict()
        genus = ""
        for ident, t in hs:
            if ident < best - 0.2:
                break
            sp, ge = lineage(t)
            genus = genus or ge
            if sp:
                tied[sp] = max(tied.get(sp, 0), ident)
        sps = list(tied)
        out.append({"copy": qn.split()[0], "best_identity": round(best, 2),
                    "best_species": sps[0] if sps else "", "best_genus": genus,
                    "tied_species": "; ".join(sps[:8]), "n_tied": len(sps)})
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  Steps 8-12 — typing, AMR, virulence, plasmids, annotation
# ══════════════════════════════════════════════════════════════════════════════
def run_mlst(asm: Path, sdir: Path, threads: int, logf: Path) -> dict:
    t = find_tool("mlst")
    if not t:
        warn("mlst not installed — MLST skipped")
        return {}
    outp = sdir / "mlst.tsv"
    if run(t, ["--threads", threads, "--label", sdir.name, asm], logf,
           stdout_to=outp) != 0:
        return {}
    line = outp.read_text().strip().splitlines()
    if not line:
        return {}
    v = line[0].split("\t")
    return {"scheme": v[1] if len(v) > 1 else "", "st": v[2] if len(v) > 2 else "",
            "alleles": " ".join(v[3:])}


_AMR_ORGS: list[str] | None = None


def amr_organism(genus: str, species: str) -> str:
    global _AMR_ORGS
    t = find_tool("amrfinder")
    if not t:
        return ""
    if _AMR_ORGS is None:
        try:
            p = subprocess.run([str(t.path), "--list_organisms"], capture_output=True,
                               text=True, timeout=120, env=t.env())
            txt = (p.stdout or "") + (p.stderr or "")
            m = re.search(r"Available --organism options:\s*(.+)", txt)
            _AMR_ORGS = [o.strip() for o in m.group(1).split(",")] if m else []
        except Exception:
            _AMR_ORGS = []
    g = re.sub(r"_[A-Z]$", "", genus or "")
    sp = re.sub(r"_[A-Z]$", "", (species or "").split(" ")[-1]) if species else ""
    for cand in (f"{g}_{sp}", g, "Escherichia" if g in ("Shigella",) else ""):
        if cand and cand in _AMR_ORGS:
            return cand
    return ""


def run_amrfinder(asm: Path, sdir: Path, organism: str, threads: int, db: str,
                  logf: Path) -> list[dict]:
    t = find_tool("amrfinder")
    if not t:
        warn("AMRFinderPlus not installed — AMR / stress / virulence genes skipped")
        return []
    outp = sdir / "amrfinder.tsv"
    args = ["-n", asm, "--plus", "--threads", threads, "-o", outp, "--name", sdir.name]
    if organism:
        args += ["--organism", organism]
    if db:
        args += ["--database", db]
    if run(t, args, logf, timeout=2 * 3600) != 0:
        warn("AMRFinderPlus failed (database missing? run: amrfinder -u)")
        return []
    out = []
    for r in read_tsv(outp):
        out.append({
            "gene": pick(r, "Element symbol", "Gene symbol"),
            "name": pick(r, "Element name", "Sequence name"),
            "scope": pick(r, "Scope"),
            "type": pick(r, "Type", "Element type"),
            "subtype": pick(r, "Subtype", "Element subtype"),
            "class": pick(r, "Class"),
            "subclass": pick(r, "Subclass"),
            "method": pick(r, "Method"),
            "contig": pick(r, "Contig id"),
            "start": pick(r, "Start"), "stop": pick(r, "Stop"),
            "strand": pick(r, "Strand"),
            "coverage": pick(r, "% Coverage of reference", "% Coverage of reference sequence"),
            "identity": pick(r, "% Identity to reference", "% Identity to reference sequence"),
            "accession": pick(r, "Closest reference accession", "Accession of closest sequence"),
        })
    return out


def run_abricate(asm: Path, sdir: Path, db: str, threads: int, logf: Path,
                 minid=80, mincov=80) -> list[dict] | None:
    t = find_tool("abricate")
    if not t:
        return None
    outp = sdir / f"abricate_{db}.tsv"
    if run(t, ["--db", db, "--threads", threads, "--minid", minid, "--mincov", mincov,
               "--nopath", asm], logf, stdout_to=outp) != 0:
        return None
    out = []
    for r in read_tsv(outp, comment_header=True):
        out.append({"gene": pick(r, "GENE"), "contig": pick(r, "SEQUENCE"),
                    "start": pick(r, "START"), "end": pick(r, "END"),
                    "strand": pick(r, "STRAND"),
                    "coverage": pick(r, "%COVERAGE"), "identity": pick(r, "%IDENTITY"),
                    "product": pick(r, "PRODUCT"), "accession": pick(r, "ACCESSION"),
                    "resistance": pick(r, "RESISTANCE"), "database": db})
    return out


def run_mobsuite(asm: Path, sdir: Path, threads: int, logf: Path) -> dict:
    t = find_tool("mob_recon")
    if not t:
        warn("MOB-suite not installed — plasmid reconstruction skipped "
             "(PlasmidFinder replicons still reported)")
        return {}
    od = sdir / "mob_recon"
    if run(t, ["--infile", asm, "--outdir", od, "--num_threads", threads, "--force"],
           logf, timeout=3 * 3600) != 0:
        warn("MOB-suite failed (database not initialised? run: mob_init)")
        return {}
    contigs = {}
    for r in read_tsv(od / "contig_report.txt"):
        cid = pick(r, "contig_id").split(" ")[0]
        contigs[cid] = {"molecule": pick(r, "molecule_type"),
                        "cluster": pick(r, "primary_cluster_id")}
    types = []
    for r in read_tsv(od / "mobtyper_results.txt"):
        types.append({
            "plasmid": pick(r, "sample_id"),
            "size": pick(r, "size"),
            "gc": pick(r, "gc"),
            "rep_types": pick(r, "rep_type(s)"),
            "relaxase": pick(r, "relaxase_type(s)"),
            "mpf": pick(r, "mpf_type"),
            "orit": pick(r, "orit_type(s)"),
            "mobility": pick(r, "predicted_mobility"),
            "cluster": pick(r, "primary_cluster_id"),
            "host_range": pick(r, "predicted_host_range_overall_name"),
            "host_rank": pick(r, "predicted_host_range_overall_rank"),
        })
    return {"contigs": contigs, "plasmids": types}


def run_checkm2(asm: Path, sdir: Path, db: str, threads: int, logf: Path) -> dict:
    t = find_tool("checkm2")
    if not t:
        warn("CheckM2 not installed — genome completeness / contamination skipped")
        return {}
    if not db or not Path(db).exists():
        warn("CheckM2 database not configured (db_paths.json key 'checkm2_db') — skipped")
        return {}
    od = sdir / "checkm2"
    inp = sdir / "checkm2_in"
    inp.mkdir(exist_ok=True)
    tgt = inp / "assembly.fasta"
    if not tgt.exists():
        shutil.copy(asm, tgt)
    if run(t, ["predict", "--input", inp, "--output-directory", od, "--threads",
               threads, "--database_path", db, "-x", "fasta", "--force"],
           logf, timeout=3 * 3600) != 0:
        warn("CheckM2 failed")
        return {}
    shutil.rmtree(inp, ignore_errors=True)
    rows = read_tsv(od / "quality_report.tsv")
    if not rows:
        return {}
    r = rows[0]
    return {"completeness": fnum(pick(r, "Completeness")),
            "contamination": fnum(pick(r, "Contamination")),
            "model": pick(r, "Completeness_Model_Used"),
            "coding_density": fnum(pick(r, "Coding_Density"))}


def run_bakta(asm: Path, sdir: Path, sample: str, db: str, genus: str, species: str,
              complete: bool, threads: int, logf: Path) -> dict:
    t = find_tool("bakta")
    if not t:
        warn("Bakta not installed — genome annotation skipped")
        return {}
    if not db or not Path(db).exists():
        warn("Bakta database not configured (db_paths.json key 'bakta_db') — annotation skipped")
        return {}
    od = sdir / "bakta"
    args = ["--db", db, "--output", od, "--prefix", sample, "--threads", threads,
            "--force", "--keep-contig-headers"]
    g = re.sub(r"_[A-Z]$", "", genus or "")
    sp = re.sub(r"_[A-Z]$", "", (species or "").split(" ")[-1]) if species else ""
    if g:
        args += ["--genus", g]
    if g and sp and sp.lower() != g.lower():
        args += ["--species", sp]
    if complete:
        args += ["--complete"]
    if run(t, args + [asm], logf, timeout=4 * 3600) != 0:
        warn("Bakta failed")
        return {}
    summ = {}
    txt = od / f"{sample}.txt"
    if txt.exists():
        for line in txt.read_text(errors="replace").splitlines():
            m = re.match(r"^\s*([A-Za-z][A-Za-z ()/'-]+):\s*([\d.,]+)\s*$", line)
            if m:
                summ[m.group(1).strip()] = m.group(2).replace(",", "")
    files = {k: str(od / f"{sample}.{k}") for k in ("gbff", "gff3", "tsv", "faa", "ffn",
                                                     "png", "svg", "json", "embl")
             if (od / f"{sample}.{k}").exists()}
    return {"summary": summ, "files": files}


# ══════════════════════════════════════════════════════════════════════════════
#  Genome windows for the map figure (GC content / skew)
# ══════════════════════════════════════════════════════════════════════════════
def gc_windows(seq: str, n: int = 360) -> dict:
    L = len(seq)
    if L == 0:
        return {"gc": [], "skew": []}
    w = max(1000, L // n)
    s = seq.upper()
    gcs, skews = [], []
    for i in range(0, L, w):
        win = s[i:i + w]
        g, c = win.count("G"), win.count("C")
        acgt = g + c + win.count("A") + win.count("T")
        gcs.append(round(100 * (g + c) / acgt, 2) if acgt else 0)
        skews.append(round((g - c) / (g + c), 4) if g + c else 0)
    return {"window": w, "gc": gcs, "skew": skews}


# ══════════════════════════════════════════════════════════════════════════════
#  Per-sample driver
# ══════════════════════════════════════════════════════════════════════════════
def process_sample(g: dict, args, dbp: dict, work_root: Path, out_root: Path,
                   pct0: float, pct1: float, input_dir: Path) -> dict:
    S = g["sample"]
    span = pct1 - pct0

    def P(frac, label):
        progress(pct0 + span * frac, f"[{S}] {label}")

    wdir = work_root / S
    wdir.mkdir(parents=True, exist_ok=True)
    sdir = out_root / "samples" / S
    sdir.mkdir(parents=True, exist_ok=True)
    logs = wdir / "logs"
    logs.mkdir(exist_ok=True)
    res: dict = {"sample": S, "barcode": g.get("barcode", ""), "folder": g.get("folder", ""),
                 "n_files": len(g["files"]), "status": "ok", "warnings": []}
    w0 = len(WARNINGS)

    seqkit = find_tool("seqkit")
    flye = find_tool("flye")
    if not seqkit or not flye:
        raise StepFailed("seqkit and flye are required — run setup_wgs.sh")

    # ── 1. merge ─────────────────────────────────────────────────────────────
    P(0.00, "Merging chunk files")
    reads = wdir / "reads.fastq.gz"
    sig = _sig(sorted(g["files"]))
    if not is_done(wdir, "merge", sig) or not (reads.exists() or reads.is_symlink()):
        clear_done(wdir, "qc", "flye", "medaka", "post")
        merge_group(input_dir, g, reads)
        mark_done(wdir, "merge", sig)
    log(f"  [{S}] {len(g['files'])} file(s) → {reads.name}")
    res["run_info"] = read_run_info(reads)

    # ── 2. read QC + filter ──────────────────────────────────────────────────
    P(0.03, "Read QC")
    raw_prof = read_profile(seqkit, reads, wdir / "reads_profile.tsv", args.threads,
                            logs / "seqkit.log")
    res["reads_raw"] = raw_prof
    filt = wdir / "reads.filtered.fastq.gz"
    qsig = _sig([sig, args.min_read_len, args.min_read_q])
    if not is_done(wdir, "qc", qsig) or not filt.exists():
        clear_done(wdir, "flye", "medaka", "post")
        P(0.05, f"Filtering reads (≥{args.min_read_len} bp, ≥Q{args.min_read_q})")
        rc = run(seqkit, ["seq", "-m", args.min_read_len, "-Q", args.min_read_q,
                          "-j", args.threads, reads, "-o", filt], logs / "seqkit.log")
        if rc != 0:
            raise StepFailed("read filtering failed")
        mark_done(wdir, "qc", qsig)
        ftsv = wdir / "filtered_profile.tsv"
        if ftsv.exists():
            ftsv.unlink()
    filt_prof = read_profile(seqkit, filt, wdir / "filtered_profile.tsv", args.threads,
                             logs / "seqkit.log")
    res["reads_filtered"] = {k: v for k, v in filt_prof.items()
                             if k in ("reads", "bases", "mean_len", "median_len", "n50",
                                      "max_len", "mean_q", "pct_q20", "gc")}
    if len(raw_prof["gc_peaks"]) > 1:
        pk = ", ".join(f"{p['gc']}% ({p['mass']*100:.0f}% of bases)"
                       for p in raw_prof["gc_peaks"])
        warn(f"[{S}] read GC has more than one peak ({pk}) — the culture may be mixed "
             "or contaminated; check species ID and CheckM2 contamination")

    # ── 3. Flye ──────────────────────────────────────────────────────────────
    gsize_txt = args.genome_size if args.genome_size != "auto" else "5m"
    gsize = parse_size(gsize_txt) or 5_000_000
    cov_est = filt_prof["bases"] / gsize
    if cov_est < 20:
        warn(f"[{S}] only ~{cov_est:.0f}x coverage of a {human_bp(gsize)} genome after "
             "filtering — expect a fragmented assembly")
    ri = res["run_info"]
    read_type = args.read_type
    if read_type == "auto":
        hq = ri.get("chemistry", "").startswith("R10") or ri.get("accuracy") == "sup" \
            or raw_prof["mean_q"] >= 15
        read_type = "nano-hq" if hq else "nano-raw"
    flye_dir = wdir / "flye"
    fsig = _sig([qsig, read_type, gsize_txt, args.asm_coverage, args.flye_meta])
    res["flye_mode"] = read_type
    if not is_done(wdir, "flye", fsig) or not (flye_dir / "assembly.fasta").exists():
        clear_done(wdir, "medaka", "post")
        P(0.08, f"Assembling with Flye (--{read_type}, ~{cov_est:.0f}x)")
        if flye_dir.exists():
            shutil.rmtree(flye_dir, ignore_errors=True)
        fargs = [f"--{read_type}", filt, "--out-dir", flye_dir, "--threads", args.threads]
        if args.asm_coverage > 0 and cov_est > args.asm_coverage * 1.5:
            fargs += ["--genome-size", gsize_txt, "--asm-coverage", args.asm_coverage]
        if args.flye_meta:
            fargs += ["--meta"]
        rc = run(flye, fargs, logs / "flye.log", timeout=12 * 3600)
        if rc != 0 or not (flye_dir / "assembly.fasta").exists():
            raise StepFailed("Flye assembly failed — see logs/flye.log in the work folder")
        mark_done(wdir, "flye", fsig)
    flye_info = parse_flye_info(flye_dir / "assembly_info.txt")

    # ── 4. Medaka ────────────────────────────────────────────────────────────
    med_dir = wdir / "medaka"
    msig = _sig([fsig, args.medaka_model, args.medaka_bacteria, args.skip_medaka,
                 args.medaka_max_cov])
    polished = med_dir / "consensus.fasta"
    if args.skip_medaka:
        res["polish"] = "skipped (by setting)"
        draft_final = flye_dir / "assembly.fasta"
    elif is_done(wdir, "medaka", msig) and polished.exists():
        res["polish"] = (wdir / ".medaka_label").read_text() if \
            (wdir / ".medaka_label").exists() else "medaka"
        draft_final = polished
    else:
        clear_done(wdir, "post")
        draft_len = sum(len(s) for _, s in read_fasta(flye_dir / "assembly.fasta"))
        med_reads = filt
        cov = filt_prof["bases"] / max(1, draft_len)
        fl = find_tool("filtlong")
        if fl and args.medaka_max_cov > 0 and cov > args.medaka_max_cov * 1.2:
            P(0.45, f"Sub-sampling reads to {args.medaka_max_cov}x for polishing")
            med_reads = wdir / "reads.medaka.fastq.gz"
            tmpfq = wdir / "reads.medaka.fastq"
            rc = run(fl, ["--target_bases", int(draft_len * args.medaka_max_cov),
                          "--min_length", args.min_read_len, filt], logs / "filtlong.log",
                     stdout_to=tmpfq)
            if rc == 0 and tmpfq.exists() and tmpfq.stat().st_size > 0:
                with open(tmpfq, "rb") as src, gzip.open(med_reads, "wb", 1) as dst:
                    shutil.copyfileobj(src, dst, 8 << 20)
                tmpfq.unlink()
            else:
                med_reads = filt
        P(0.48, "Polishing with Medaka")
        cons, label = run_medaka(med_reads, flye_dir / "assembly.fasta", med_dir,
                                 args.threads, args.medaka_model, args.medaka_bacteria,
                                 logs / "medaka.log")
        if med_reads != filt and med_reads.exists():
            med_reads.unlink()
        if cons:
            res["polish"] = f"medaka ({label})"
            draft_final = cons
            (wdir / ".medaka_label").write_text(res["polish"])
            mark_done(wdir, "medaka", msig)
        else:
            res["polish"] = f"none — medaka {label}"
            draft_final = flye_dir / "assembly.fasta"

    # ── 5. final assembly + contig table ─────────────────────────────────────
    P(0.70, "Assembly QC")
    contigs = []
    recs = sorted(read_fasta(draft_final), key=lambda r: -len(r[1]))
    final_recs = []
    for i, (name, seq) in enumerate(recs, 1):
        fname = name.split()[0]
        fi = flye_info.get(fname, {})
        new = f"{S}_c{i}"
        c = {"contig": new, "flye_name": fname, "length": len(seq),
             "gc": round(gc_pct(seq), 2), "depth": fi.get("cov", ""),
             "circular": bool(fi.get("circular")), "repeat": bool(fi.get("repeat"))}
        contigs.append(c)
        topo = "circular" if c["circular"] else "linear"
        final_recs.append((f"{new} len={len(seq)} depth={c['depth']}x topology={topo}", seq))
    classify_contigs(contigs)
    asm = sdir / f"{S}_assembly.fasta"
    write_fasta(asm, final_recs)
    shutil.copy(flye_dir / "assembly_info.txt", sdir / "flye_assembly_info.txt")
    if (flye_dir / "assembly_graph.gfa").exists():
        shutil.copy(flye_dir / "assembly_graph.gfa", sdir / f"{S}_assembly_graph.gfa")
    lens = [c["length"] for c in contigs]
    chrom = next((c for c in contigs if c["type"] == "chromosome"), None)
    tot = sum(lens) or 1
    res["assembly"] = {
        "contigs": len(contigs), "total_length": sum(lens), "n50": n50(lens),
        "largest": max(lens) if lens else 0,
        "gc": round(sum(c["gc"] * c["length"] for c in contigs) / tot, 2),
        "circular": sum(1 for c in contigs if c["circular"]),
        "chromosome_closed": bool(chrom and chrom["circular"]),
        "mean_depth": round(sum((fnum(c["depth"], 0) or 0) * c["length"] for c in contigs)
                            / tot, 1),
        "polish": res.get("polish", ""),
    }
    res["contigs"] = contigs
    seqs_by_contig = {f"{S}_c{i}": s for i, (_, s) in enumerate(recs, 1)}
    res["gc_windows"] = {c["contig"]: gc_windows(seqs_by_contig[c["contig"]])
                         for c in contigs if c["length"] >= 20000}

    post_sig = _sig([msig, __version__, args.skip_annotation, args.skip_checkm2,
                     args.skip_mobsuite, args.skip_amr, args.skip_vf,
                     dbp.get("gtdb_sourmash", ""), dbp.get("bakta_db", "")])
    # (post steps are cheap relative to Flye/Medaka and always rewrite their own
    #  outputs; the signature only decides whether the cleanup below may run)

    # ── 6. species ID ────────────────────────────────────────────────────────
    P(0.72, "Species identification (sourmash / GTDB)")
    sm_hits = run_sourmash(asm, sdir, dbp.get("gtdb_sourmash", ""),
                           dbp.get("gtdb_lineages", ""), args.threads, logs / "sourmash.log")
    res["sourmash"] = sm_hits[:10]
    top = sm_hits[0] if sm_hits else {}
    ident = {"method": "", "species": "", "genus": "", "ani": None, "af": None,
             "reference": "", "accession": "", "call": "unresolved"}
    if top:
        ident.update({"method": "sourmash gather (GTDB)", "species": top["gtdb_species"]
                      or top["reference_name"], "genus": top["gtdb_genus"],
                      "ani": fnum(top["ani"]), "reference": top["reference_name"],
                      "accession": top["accession"], "lineage": top["gtdb_lineage"]})
        # gather hands the k-mers the top genome does not explain to the next
        # best genome — for a pure isolate that is usually a sister species of
        # the same genus soaking up strain-specific and error k-mers. Only a
        # different GENUS holding a real share of the assembly is evidence of a
        # second organism; CheckM2 contamination is the authoritative check.
        others = [h for h in sm_hits[1:] if h["gtdb_genus"] and top["gtdb_genus"]
                  and h["gtdb_genus"] != top["gtdb_genus"] and h["f_query"] >= 5]
        if others:
            warn(f"[{S}] {len(others)} other genus/genera carry ≥5 % of the assembly "
                 f"({', '.join(h['gtdb_species'] for h in others[:3])}) — possible "
                 "mixed culture")
        if args.verify_ani and top["accession"]:
            P(0.74, f"Verifying ANI against {top['accession']}")
            ref = fetch_ncbi_genome(top["accession"], Path(args.ref_cache))
            sk = run_skani(asm, ref, logs / "skani.log", sdir / "skani.tsv") if ref else {}
            if sk.get("ani") is not None:
                ident.update({"method": "skani ANI vs closest GTDB reference",
                              "ani": sk["ani"], "af": sk.get("af_query"),
                              "af_ref": sk.get("af_ref")})
            elif not ref:
                log("    (closest reference genome not downloadable — keeping the "
                    "sourmash ANI estimate)")
    ident["call"] = species_call(ident["ani"])

    # ── 7. 16S ───────────────────────────────────────────────────────────────
    P(0.76, "Extracting 16S rRNA genes")
    rrna = {"copies": 0, "partial": 0, "operons_23S": 0, "operons_5S": 0, "hits": []}
    br = find_tool("barrnap")
    if br:
        gff = sdir / "barrnap.gff"
        # Sequences are cut from the assembly here rather than taken from
        # barrnap --outseq: barrnap 0.9 names them "16S_rRNA::contig:a-b(+)",
        # 1.x only "contig:a-b(+)", so the name no longer says which gene it is.
        if run(br, ["--kingdom", "bac", "--threads", args.threads, asm],
               logs / "barrnap.log", stdout_to=gff) == 0:
            feats = []
            for line in gff.read_text().splitlines():
                if line.startswith("#"):
                    continue
                v = line.split("\t")
                if len(v) >= 9:
                    nm = re.search(r"Name=([^;]+)", v[8])
                    feats.append({"contig": v[0], "start": int(v[3]), "end": int(v[4]),
                                  "strand": v[6], "name": nm.group(1) if nm else "",
                                  "partial": "partial" in v[8].lower()})
            res["rrna_features"] = feats
            rrna["operons_23S"] = sum(1 for f in feats if f["name"].startswith("23S"))
            rrna["operons_5S"] = sum(1 for f in feats if f["name"].startswith("5S"))
            comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
            s16 = []
            for f in sorted((f for f in feats if f["name"].startswith("16S")),
                            key=lambda f: (natural_key(f["contig"]), f["start"])):
                cs = seqs_by_contig.get(f["contig"], "")
                seq = cs[f["start"] - 1: f["end"]]
                if f["strand"] == "-":
                    seq = seq.translate(comp)[::-1]
                loc = f"{f['contig']}:{f['start']}-{f['end']}({f['strand']})"
                if len(seq) >= 1200 and not f["partial"]:
                    s16.append((f"{S}_16S_{len(s16) + 1} {loc} len={len(seq)}", seq))
                else:
                    rrna["partial"] += 1
            rrna["copies"] = len(s16)
            if s16:
                write_fasta(sdir / f"{S}_16S.fasta", s16)
                uniq = len({s for _, s in s16})
                rrna["distinct_sequences"] = uniq
                db_fa, tax_dir, db_label = resolve_16s_db(args.db_16s, dbp)
                if db_fa:
                    P(0.77, f"Classifying 16S copies ({db_label})")
                    rrna["db"] = db_label
                    rrna["hits"] = classify_16s(s16, db_fa, tax_dir, args.threads, wdir,
                                                logs / "vsearch.log")
                else:
                    warn(f"[{S}] no 16S database found — 16S copies extracted but not "
                         "classified")
            else:
                warn(f"[{S}] no full-length 16S gene found in the assembly")
    else:
        warn("barrnap not installed — 16S extraction skipped")
    res["rrna"] = rrna
    if rrna["hits"]:
        # What 16S can honestly say: the species every copy agrees on, or — when
        # copies tie between close species (Cronobacter sakazakii/malonaticus,
        # B. cereus group…) — the genus with the tied names, never one of them
        # picked by a 0.1 % identity difference.
        sets = [set(filter(None, h["tied_species"].split("; "))) for h in rrna["hits"]
                if h["tied_species"]]
        common = set.intersection(*sets) if sets else set()
        genus16 = Counter(h["best_genus"] for h in rrna["hits"] if h["best_genus"])
        g16 = genus16.most_common(1)[0][0] if genus16 else ""
        if len(common) == 1:
            res["rrna"]["consensus_species"] = next(iter(common))
            res["rrna"]["resolved"] = True
        else:
            names = sorted(common or set.union(*sets) if sets else set())
            res["rrna"]["consensus_species"] = (f"{g16} sp." if g16 else "") + (
                f" (16S ties: {', '.join(names[:4])})" if names else "")
            res["rrna"]["resolved"] = False
        res["rrna"]["consensus_genus"] = g16
        res["rrna"]["max_tied"] = max(h["n_tied"] for h in rrna["hits"])
    if not ident["species"] and rrna.get("consensus_species"):
        ident.update({"method": "16S rRNA (genus-level reliability)",
                      "species": rrna["consensus_species"],
                      "genus": rrna.get("consensus_genus", ""), "call": "16S only"})
    res["identification"] = ident

    genus = ident.get("genus") or (ident.get("species") or "").split(" ")[0]
    species = ident.get("species") if ident.get("call") == "species" else ""

    # ── 8. MLST ──────────────────────────────────────────────────────────────
    P(0.79, "MLST")
    res["mlst"] = run_mlst(asm, sdir, args.threads, logs / "mlst.log")

    # ── 9. AMRFinderPlus ─────────────────────────────────────────────────────
    if not args.skip_amr:
        org = amr_organism(genus, species)
        P(0.81, f"AMR genes (AMRFinderPlus{', ' + org if org else ''})")
        res["amr_organism"] = org
        res["amrfinder"] = run_amrfinder(asm, sdir, org, args.threads,
                                         dbp.get("amrfinder_db", ""), logs / "amrfinder.log")
    else:
        res["amrfinder"] = []

    # ── 10. VFDB / PlasmidFinder via abricate ───────────────────────────────
    if not args.skip_vf:
        P(0.84, "Virulence factors (VFDB)")
        vf = run_abricate(asm, sdir, "vfdb", args.threads, logs / "abricate.log")
        if vf is None:
            warn("abricate not installed — VFDB virulence screen skipped")
        res["vfdb"] = vf or []
    else:
        res["vfdb"] = []
    P(0.85, "Plasmid replicons (PlasmidFinder)")
    res["plasmidfinder"] = run_abricate(asm, sdir, "plasmidfinder", args.threads,
                                        logs / "abricate.log", minid=80, mincov=60) or []

    # ── 11. MOB-suite ────────────────────────────────────────────────────────
    if not args.skip_mobsuite:
        P(0.86, "Plasmid reconstruction (MOB-suite)")
        mob = run_mobsuite(asm, sdir, args.threads, logs / "mobsuite.log")
        res["mobsuite"] = mob
        for c in contigs:
            m = mob.get("contigs", {}).get(c["contig"])
            if m:
                c["mob_molecule"] = m["molecule"]
                c["mob_cluster"] = m["cluster"]
                if m["molecule"] == "plasmid":
                    c["type"] = "plasmid"
                elif m["molecule"] == "chromosome" and c["type"] == "plasmid?":
                    c["type"] = "contig"
    for c in contigs:
        reps = sorted({h["gene"] for h in res["plasmidfinder"] if h["contig"] == c["contig"]})
        c["replicons"] = ", ".join(reps)
        if reps and c["type"] == "plasmid?":
            c["type"] = "plasmid"
    res["assembly"]["plasmids"] = sum(1 for c in contigs if c["type"] == "plasmid")

    # ── 5b. CheckM2 ──────────────────────────────────────────────────────────
    if not args.skip_checkm2:
        P(0.88, "Completeness / contamination (CheckM2)")
        res["checkm2"] = run_checkm2(asm, sdir, dbp.get("checkm2_db", ""), args.threads,
                                     logs / "checkm2.log")
        ck = res["checkm2"]
        if ck.get("contamination") is not None and ck["contamination"] > 5:
            warn(f"[{S}] CheckM2 contamination {ck['contamination']:.1f}% (>5 %) — "
                 "the isolate may not be pure")
        if ck.get("completeness") is not None and ck["completeness"] < 90:
            warn(f"[{S}] CheckM2 completeness {ck['completeness']:.1f}% (<90 %)")
    else:
        res["checkm2"] = {}

    # ── 12. Bakta ────────────────────────────────────────────────────────────
    if not args.skip_annotation:
        P(0.90, "Genome annotation (Bakta)")
        complete = bool(contigs) and all(c["circular"] for c in contigs)
        res["bakta"] = run_bakta(asm, sdir, S, dbp.get("bakta_db", ""), genus, species,
                                 complete, args.threads, logs / "bakta.log")
    else:
        res["bakta"] = {}

    # ── cleanup: keep merged reads + final products, drop the bulky middle ──
    if not args.keep_intermediate:
        for p in (flye_dir / "00-assembly", flye_dir / "10-consensus",
                  flye_dir / "20-repeat", flye_dir / "30-contigger",
                  flye_dir / "40-polishing"):
            shutil.rmtree(p, ignore_errors=True)
        for p in med_dir.glob("*.bam*") if med_dir.exists() else []:
            p.unlink()
        for p in med_dir.glob("*.hdf") if med_dir.exists() else []:
            p.unlink()
    for f in logs.glob("*.log"):
        shutil.copy(f, sdir / f"log_{f.name}")
    mark_done(wdir, "post", post_sig)
    res["warnings"] = WARNINGS[w0:]
    (sdir / "result.json").write_text(json.dumps(res, indent=1, default=str))
    P(1.0, "Done")
    return res


# ══════════════════════════════════════════════════════════════════════════════
#  Tables across samples
# ══════════════════════════════════════════════════════════════════════════════
def write_tables(results: list[dict], out: Path):
    summ, readqc, contig_rows, sp_rows, r16, mlst, amr, vf, pls, ann = ([] for _ in range(10))
    for r in results:
        S = r["sample"]
        if r.get("status") != "ok":
            summ.append({"sample": S, "barcode": r.get("barcode", ""),
                         "status": r.get("status"), "error": r.get("error", "")})
            continue
        a, ident = r.get("assembly", {}), r.get("identification", {})
        ck, ml = r.get("checkm2", {}), r.get("mlst", {})
        amr_core = [x for x in r.get("amrfinder", []) if x.get("type") == "AMR"]
        vir = [x for x in r.get("amrfinder", []) if x.get("type") == "VIRULENCE"]
        rq, fq = r.get("reads_raw", {}), r.get("reads_filtered", {})
        summ.append({
            "sample": S, "barcode": r.get("barcode", ""), "status": "ok",
            "species": ident.get("species", ""), "ani": ident.get("ani") or "",
            "identification": ident.get("call", ""), "id_method": ident.get("method", ""),
            "closest_reference": ident.get("reference", ""),
            "reference_accession": ident.get("accession", ""),
            "16S_species": r.get("rrna", {}).get("consensus_species", ""),
            "16S_copies": r.get("rrna", {}).get("copies", 0),
            "genome_size": a.get("total_length", ""), "contigs": a.get("contigs", ""),
            "circular": a.get("circular", ""),
            "chromosome_closed": a.get("chromosome_closed", ""),
            "plasmids": a.get("plasmids", ""), "n50": a.get("n50", ""),
            "gc": a.get("gc", ""), "depth": a.get("mean_depth", ""),
            "completeness": ck.get("completeness", ""),
            "contamination": ck.get("contamination", ""),
            "mlst_scheme": ml.get("scheme", ""), "st": ml.get("st", ""),
            "amr_genes": len(amr_core), "virulence_genes_amrfinder": len(vir),
            "vfdb_hits": len(r.get("vfdb", [])),
            "cds": r.get("bakta", {}).get("summary", {}).get("CDSs", ""),
            "polish": a.get("polish", ""), "flye_mode": r.get("flye_mode", ""),
            "reads": rq.get("reads", ""), "bases": rq.get("bases", ""),
            "read_n50": rq.get("n50", ""),
            "warnings": " | ".join(r.get("warnings", [])),
        })
        readqc.append({"sample": S, "stage": "raw", **{k: rq.get(k, "") for k in
                       ("reads", "bases", "mean_len", "median_len", "n50", "max_len",
                        "mean_q", "pct_q10", "pct_q20", "gc")}})
        readqc.append({"sample": S, "stage": "filtered", **{k: fq.get(k, "") for k in
                       ("reads", "bases", "mean_len", "median_len", "n50", "max_len",
                        "mean_q", "pct_q20", "gc")}})
        for c in r.get("contigs", []):
            contig_rows.append({"sample": S, **{k: c.get(k, "") for k in
                                ("contig", "type", "length", "gc", "depth", "circular",
                                 "replicons", "mob_cluster", "flye_name")}})
        for h in r.get("sourmash", []):
            sp_rows.append({"sample": S, **h})
        for h in r.get("rrna", {}).get("hits", []):
            r16.append({"sample": S, **h})
        if ml:
            mlst.append({"sample": S, **ml})
        for x in r.get("amrfinder", []):
            amr.append({"sample": S, **x})
        for x in r.get("vfdb", []):
            vf.append({"sample": S, **x})
        for p in r.get("mobsuite", {}).get("plasmids", []):
            pls.append({"sample": S, **p})
        for k, v in r.get("bakta", {}).get("summary", {}).items():
            ann.append({"sample": S, "feature": k, "count": v})
    write_csv(out / "wgs_summary.csv", summ)
    write_csv(out / "read_qc.csv", readqc)
    write_csv(out / "assembly_contigs.csv", contig_rows)
    write_csv(out / "species_id_gtdb.csv", sp_rows)
    write_csv(out / "rrna_16S_hits.csv", r16)
    write_csv(out / "mlst.csv", mlst)
    write_csv(out / "amr_stress_virulence_amrfinder.csv", amr)
    write_csv(out / "virulence_vfdb.csv", vf)
    write_csv(out / "plasmids_mobsuite.csv", pls)
    write_csv(out / "annotation_summary.csv", ann)


# ══════════════════════════════════════════════════════════════════════════════
#  main
# ══════════════════════════════════════════════════════════════════════════════
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="ONT whole-genome isolate pipeline")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--work_dir", default="")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--job_name", default="")
    p.add_argument("--db_paths", default="")
    p.add_argument("--sample_names", default="", help="JSON {group key: sample name}")
    p.add_argument("--min_read_len", type=int, default=1000)
    p.add_argument("--min_read_q", type=float, default=10)
    p.add_argument("--genome_size", default="auto")
    p.add_argument("--asm_coverage", type=int, default=100)
    p.add_argument("--read_type", default="auto",
                   choices=["auto", "nano-hq", "nano-raw", "nano-corr"])
    p.add_argument("--flye_meta", action="store_true")
    p.add_argument("--medaka_model", default="auto")
    p.add_argument("--medaka_bacteria", type=lambda s: str(s).lower() in
                   ("1", "true", "yes"), default=True)
    p.add_argument("--medaka_max_cov", type=int, default=150)
    p.add_argument("--skip_medaka", action="store_true")
    p.add_argument("--db_16s", default="")
    p.add_argument("--verify_ani", type=lambda s: str(s).lower() in
                   ("1", "true", "yes"), default=True)
    p.add_argument("--ref_cache", default=str(Path.home() / ".ngamp_wgs_refs"))
    p.add_argument("--skip_amr", action="store_true")
    p.add_argument("--skip_vf", action="store_true")
    p.add_argument("--skip_mobsuite", action="store_true")
    p.add_argument("--skip_checkm2", action="store_true")
    p.add_argument("--skip_annotation", action="store_true")
    p.add_argument("--keep_intermediate", action="store_true")
    p.add_argument("--app_version", default="")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    t0 = time.time()
    inp, out = Path(args.input), Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(args.work_dir) if args.work_dir else inp / "_wgs_work"
    work.mkdir(parents=True, exist_ok=True)
    dbp = {}
    if args.db_paths and Path(args.db_paths).exists():
        try:
            dbp = json.loads(Path(args.db_paths).read_text())
        except Exception as e:
            warn(f"db_paths.json unreadable: {e}")
    overrides = {}
    if args.sample_names:
        try:
            overrides = json.loads(args.sample_names)
        except Exception:
            p = Path(args.sample_names)
            if p.exists():
                overrides = json.loads(p.read_text())

    log("=" * 70)
    log(f"  NextGen-Amplicon ONT-WGS isolate pipeline v{__version__}")
    log(f"  conda base: {_BASE or '(none — using PATH)'}")
    log("=" * 70)
    progress(1, "Grouping chunk files into samples")
    groups = load_groups(inp, overrides)
    if not groups:
        log("ERROR: no FASTQ files found in the upload")
        return 2
    if groups[0].get("skipped_fail"):
        warn(f"{groups[0]['skipped_fail']} file(s) from MinKNOW 'fail' folders were left "
             "out (below the run's quality threshold)")
    map_rows = []
    for g in groups:
        log(f"  {g['sample']:<24} {g.get('barcode') or '-':<14} {len(g['files']):>4} file(s)"
            f"{'   from ' + g['folder'] if g.get('folder') else ''}")
        for f in g["files"]:
            map_rows.append({"sample": g["sample"], "barcode": g.get("barcode", ""),
                             "source_folder": g.get("folder", ""), "file": f})
    write_csv(out / "sample_map.csv", map_rows)
    (out / "sample_groups.json").write_text(json.dumps(
        [{k: g[k] for k in ("key", "sample", "barcode", "folder")} | {"n_files": len(g["files"])}
         for g in groups], indent=1))

    versions = {}
    for t in ("seqkit", "flye", "medaka", "sourmash", "skani", "barrnap", "vsearch", "mlst",
              "amrfinder", "abricate", "mob_recon", "checkm2", "bakta", "filtlong"):
        v = tool_version(t)
        if v:
            versions[t] = v
    log("  tools: " + ", ".join(f"{k} {v}" for k, v in versions.items()))
    if "flye" not in versions:
        log("ERROR: Flye not found — run setup_wgs.sh to install the WGS tools")
        return 3

    results = []
    n = len(groups)
    for i, g in enumerate(groups):
        a, b = 3 + 90 * i / n, 3 + 90 * (i + 1) / n
        log("")
        log(f"── Sample {i + 1}/{n}: {g['sample']} " + "─" * 40)
        try:
            results.append(process_sample(g, args, dbp, work, out, a, b, inp))
        except StepFailed as e:
            warn(f"[{g['sample']}] FAILED: {e}")
            results.append({"sample": g["sample"], "barcode": g.get("barcode", ""),
                            "status": "failed", "error": str(e)})
        except Exception as e:
            import traceback
            traceback.print_exc()
            warn(f"[{g['sample']}] FAILED: {type(e).__name__}: {e}")
            results.append({"sample": g["sample"], "barcode": g.get("barcode", ""),
                            "status": "failed", "error": f"{type(e).__name__}: {e}"})

    progress(94, "Writing tables")
    write_tables(results, out)
    ok = [r for r in results if r.get("status") == "ok"]
    summary = {
        "job_name": args.job_name or out.name, "marker": "ONT-WGS", "platform": "ONT",
        "pipeline": "ont-wgs", "pipeline_version": __version__,
        "app_version": args.app_version,
        "n_samples": len(results), "n_ok": len(ok),
        "samples": [r["sample"] for r in results],
        "settings": {k: getattr(args, k) for k in (
            "min_read_len", "min_read_q", "genome_size", "asm_coverage", "read_type",
            "medaka_model", "medaka_bacteria", "medaka_max_cov", "skip_medaka",
            "verify_ani", "skip_amr", "skip_vf", "skip_mobsuite", "skip_checkm2",
            "skip_annotation")},
        "databases": {k: dbp.get(k, "") for k in ("gtdb_sourmash", "gtdb_lineages",
                                                  "bakta_db", "checkm2_db", "amrfinder_db")},
        "tool_versions": versions, "warnings": WARNINGS,
        "runtime_min": round((time.time() - t0) / 60, 1),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "has_taxonomy": False,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "wgs_results.json").write_text(json.dumps(results, indent=1, default=str))

    progress(96, "Building report")
    try:
        sys.path.insert(0, str(Path(__file__).parent))
        import wgs_report
        html = wgs_report.build_html(out)
        (out / "wgs_report.html").write_text(html, encoding="utf-8")
        log("  Saved: wgs_report.html")
    except Exception as e:
        import traceback
        traceback.print_exc()
        warn(f"report not built: {e}")

    log("")
    log(f"  {len(ok)}/{len(results)} sample(s) assembled in {summary['runtime_min']} min")
    for r in results:
        if r.get("status") == "ok":
            idn = r["identification"]
            a = r["assembly"]
            ani = f"ANI {idn['ani']:.2f}%" if isinstance(idn.get("ani"), (int, float)) else ""
            log(f"  {r['sample']:<20} {idn.get('species') or '?':<32} {ani:<12} "
                f"{human_bp(a['total_length']):>9}  {a['contigs']} contig(s)  "
                f"ST {r.get('mlst', {}).get('st') or '-'}")
        else:
            log(f"  {r['sample']:<20} FAILED — {r.get('error', '')}")
    progress(100, "Complete")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
