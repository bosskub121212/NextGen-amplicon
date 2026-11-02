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

ENVS = ("ngamp-amp", "ngamp-meta", "qc", "ngamp-wgs")


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
