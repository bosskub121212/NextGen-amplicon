"""qc_tools.py — FastQC and fastp checks for the customer report, with the report's
"Result:" sentences written from what the tools actually found.

FastQC (adapter-content module) and fastp (adapter detection for pairs) are run on the
raw reads of one sample; nothing is written over the reads.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

ENVS = ("ngamp-amp", "ngamp-nano", "ngamp-meta", "qc", "ngamp-wgs")


def find_tool(name: str) -> str | None:
    for b in ("miniconda3", "anaconda3", "miniforge3", "mambaforge"):
        for e in ENVS:
            c = Path.home() / b / "envs" / e / "bin" / name
            if c.exists():
                return str(c)
    return shutil.which(name)


def fastqc(files: list[Path], out: Path, tool: str | None, threads: int = 2) -> dict | None:
    """{"version", "adapter": {file stem: pass|warn|fail}} or None when FastQC is missing."""
    if not tool:
        return None
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run([tool, "-q", "-t", str(threads), "-o", str(out)] + [str(f) for f in files],
                       capture_output=True, text=True, timeout=3600)
    except (subprocess.TimeoutExpired, OSError):
        pass
    res = {"version": "", "adapter": {}}
    for z in sorted(out.glob("*_fastqc.zip")):
        try:
            with zipfile.ZipFile(z) as zf:
                name = next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))
                txt = zf.read(name).decode(errors="replace")
        except Exception:
            continue
        m = re.search(r"##FastQC\t(\S+)", txt)
        if m:
            res["version"] = m.group(1)
        a = re.search(r">>Adapter Content\t(\w+)", txt)
        res["adapter"][z.name[:-len("_fastqc.zip")]] = a.group(1) if a else "unknown"
    return res


def fastp_check(r1: Path, r2: Path | None, out: Path, tool: str | None, threads: int = 2) -> dict | None:
    """fastp's adapter detection (pairs: by overlap). {"version", "total_reads",
    "adapter_trimmed_reads", "adapters"} or None when fastp is missing."""
    if not tool:
        return None
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    js = out / "fastp.json"
    args = [tool, "-i", str(r1), "-j", str(js), "-h", str(out / "fastp.html"), "-w", str(threads)]
    if r2:
        args += ["-I", str(r2), "--detect_adapter_for_pe"]
    try:
        subprocess.run(args, capture_output=True, text=True, timeout=3600)
        d = json.loads(js.read_text())
    except Exception:
        return {"error": "fastp ran but gave no result"}
    ac = d.get("adapter_cutting", {}) or {}
    return {"version": d.get("summary", {}).get("fastp_version", ""),
            "total_reads": d.get("summary", {}).get("before_filtering", {}).get("total_reads", 0),
            "adapter_trimmed_reads": int(ac.get("adapter_trimmed_reads", 0) or 0),
            "adapters": [ac[k] for k in ("read1_adapter_sequence", "read2_adapter_sequence")
                         if ac.get(k) and ac.get(k) != "unspecified"]}


def sentence_fastqc(d: dict | None) -> str:
    if not d:
        return "FastQC was not available on this machine; adapter content was not checked."
    if not d.get("adapter"):
        return "FastQC ran but gave no result for these files; adapter content was not checked."
    bad = [k for k, v in d["adapter"].items() if v in ("warn", "fail")]
    if not bad:
        return "No adapter sequences in raw data"
    return (f"Adapter sequences were detected in {len(bad)} of {len(d['adapter'])} file(s) "
            f"(FastQC Adapter Content: {', '.join(sorted({d['adapter'][k] for k in bad}))}).")


def sentence_fastp(d: dict | None) -> str:
    if not d:
        return "fastp was not available on this machine; adapter trimming was not checked."
    if d.get("error"):
        return "fastp ran but gave no result for these files; adapters were not checked."
    n, tot = d.get("adapter_trimmed_reads", 0), d.get("total_reads", 0)
    if n == 0:
        return "No adapter was found in the reads; no read needed adapter trimming."
    pct = 100 * n / tot if tot else 0
    return f"Adapters were found and trimmed in {n:,} reads ({pct:.1f}% of {tot:,})."


# ── Nanopore reads ─────────────────────────────────────────────────────────────
def read_stats(files: list[Path]) -> dict | None:
    """Reads, bases, mean / median / N50 length and mean read quality of FASTQ files
    (plain or gzip), all reads counted. A read's quality is computed from its mean error
    rate, -10·log10(mean 10^(-Q/10)), as MinKNOW and NanoPlot do (the mean of the Phred
    numbers themselves overstates it badly); the mean of those per-read values is given."""
    import gzip
    import math
    from array import array
    err = [10 ** (-(c - 33) / 10) if c >= 33 else 1.0 for c in range(256)]
    lens, qsum, n = array("I"), 0.0, 0
    try:
        for f in files:
            op = gzip.open if str(f).endswith(".gz") else open
            with op(f, "rb") as fh:
                while True:
                    h = fh.readline()
                    if not h:
                        break
                    seq = fh.readline().rstrip(b"\r\n")
                    fh.readline()
                    q = fh.readline().rstrip(b"\r\n")
                    L = len(seq)
                    lens.append(L)
                    if L and q:
                        e_ = sum(map(err.__getitem__, q)) / len(q)
                        qsum += -10 * math.log10(max(e_, 1e-10))
                    n += 1
    except (OSError, EOFError, ValueError):            # unreadable or truncated FASTQ: no numbers, no crash
        return None
    if not n:
        return {"reads": 0}
    srt = sorted(lens, reverse=True)
    tot, acc, n50 = sum(lens), 0, 0
    for L in srt:
        acc += L
        if acc >= tot / 2:
            n50 = L
            break
    return {"reads": n, "bases": tot, "mean_len": tot / n, "median_len": srt[n // 2], "n50": n50,
            "mean_q": qsum / n}


def sentence_read_stats(d: dict | None) -> str:
    if not d:
        return "The raw reads could not be read; read statistics were not made."
    if not d.get("reads"):
        return "The raw FASTQ held no reads."
    return (f"{d['reads']:,} reads, {d['bases']:,} bases; mean length {d['mean_len']:,.0f} bases (median "
            f"{d['median_len']:,}, N50 {d['n50']:,}); mean read quality Q{d['mean_q']:.1f} (each read's quality "
            "from its mean error rate, as MinKNOW reports it).")


def nanoplot(files: list[Path], out: Path, tool: str | None, threads: int = 2) -> dict | None:
    """NanoPlot figures of the raw reads in `out`; {"version"} or None when NanoPlot is missing or fails."""
    if not tool or not files:
        return None
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run([tool, "--fastq"] + [str(f) for f in files] + ["-o", str(out), "-t", str(threads)],
                       capture_output=True, text=True, timeout=3600)
        v = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=60)
    except (subprocess.TimeoutExpired, OSError):
        return None
    if not any(out.glob("*.html")) and not any(out.glob("*.png")):
        return None
    ver = (v.stdout or v.stderr).strip().split()[-1] if (v.stdout or v.stderr).strip() else ""
    return {"version": ver}
