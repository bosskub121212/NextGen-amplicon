#!/usr/bin/env python3
"""provenance.py — which software and which reference files produced a result.

Every amplicon result folder gets a provenance.json:

    {
      "app_version": "2.17.0",
      "pipeline":    "dada2",
      "recorded_at": "2026-10-17T14:02:11",
      "tools":       {"R": "4.5.1", "dada2": "1.36.0", "cutadapt": "4.9", ...},
      "databases":   [{"role": "taxonomy", "path": ".../silva_nr99_v144_toGenus_trainset.fa.gz",
                       "file": "...", "label": "SILVA 144", "bytes": 201326592,
                       "md5": "9e8a4c92...", "md5_expected": "9e8a4c92...", "md5_ok": true}]
    }

It is written by main.py right after a pipeline finishes, from the same
environment the pipeline ran in, so a question asked six months later ("which
DADA2, which SILVA?") has an answer in the folder itself. The report prints it in
the Method section.

Database checksums are cached in databases/db_manifest.json keyed on path, size
and modification time, so a 200 MB reference is hashed once, not once per run.
download_databases.sh writes the same manifest with the publisher's checksum,
which lets the report say whether the file in use is the file that was published.

Command line (for a result folder made before provenance existed):
    python3 provenance.py RESULT_DIR [--pipeline dada2] [--db PATH ...]
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess

import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DB_DIR = HERE / "databases"
MANIFEST = DB_DIR / "db_manifest.json"

# R packages that matter for each pipeline's numbers, in the order the report lists them
R_PKGS = {
    "dada2":  ["dada2", "DECIPHER", "phyloseq", "vegan", "Biostrings", "ShortRead"],
    "its":    ["dada2", "phyloseq", "vegan", "Biostrings", "ShortRead"],
    "cox1":   ["dada2", "lulu", "phyloseq", "vegan", "Biostrings", "ShortRead"],
    "pacbio": ["dada2", "phyloseq", "vegan", "Biostrings", "ShortRead"],
    "emu":    ["phyloseq", "vegan"],
    "qiime2_vsearch": ["phyloseq", "vegan"],
    "shotgun": [],
}
CLI_TOOLS = {
    "dada2":  ["cutadapt"],
    "its":    ["cutadapt", "vsearch"],
    "cox1":   ["cutadapt", "vsearch"],
    "pacbio": ["cutadapt"],
    "emu":    ["emu", "minimap2", "chopper", "cutadapt", "vsearch"],
    "qiime2_vsearch": ["vsearch", "cutadapt", "chopper"],
    "shotgun": ["kraken2", "bracken", "fastp", "hostile", "chopper", "seqkit"],
}


# ── tools ─────────────────────────────────────────────────────────────────────
def find_tool(name: str) -> str | None:
    """On PATH, else in one of the app's conda environments (same rule as the pipelines)."""
    hit = shutil.which(name)
    if hit:
        return hit
    bases = [Path(os.environ["CONDA_PREFIX"]).parent.parent] if os.environ.get("CONDA_PREFIX") else []
    bases += [Path.home() / d for d in ("miniconda3", "anaconda3", "miniforge3", "mambaforge")]
    for b in bases:
        for env in ("emu", "ngamp-wgs", "ngamp-phylo", "wgs", "picrust2", "ngamp-meta"):
            c = b / "envs" / env / "bin" / name
            if c.exists() and os.access(c, os.X_OK):
                return str(c)
    return None


def tool_version(name: str) -> str:
    exe = find_tool(name)
    if not exe:
        return ""
    try:
        p = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60,
                           env={**os.environ, "PATH": str(Path(exe).parent) + os.pathsep
                                + os.environ.get("PATH", "")})
    except Exception:
        return ""
    lines = ((p.stdout or "") + "\n" + (p.stderr or "")).splitlines()
    for ln in lines:                                  # a bare "4.9" line wins
        m = re.fullmatch(r"\s*v?(\d{1,3}\.\d+(?:\.\d+){0,2})\s*", ln)
        if m:
            return m.group(1)
    stem = name.lower()
    # vsearch prints a DOI ("10.7717/...") before its version: prefer lines naming the tool
    for ln in sorted(lines, key=lambda s: stem not in s.lower()):
        m = re.search(r"(?<![\w./])v?(\d{1,3}\.\d+(?:\.\d+){0,2}(?:-r\d+)?)(?![\d./])", ln)
        if m:
            return m.group(1)
    return ""


def r_versions(pkgs: list[str]) -> dict:
    rs = shutil.which("Rscript")
    if not rs:
        return {}
    code = ("cat('R\\t', paste(R.version$major, R.version$minor, sep='.'), '\\n', sep='');"
            f"for (p in c({', '.join(repr(x) for x in pkgs)})) "
            "cat(p, '\\t', tryCatch(as.character(utils::packageVersion(p)), "
            "error=function(e) ''), '\\n', sep='')")
    try:
        p = subprocess.run([rs, "--vanilla", "-e", code], capture_output=True, text=True,
                           timeout=120, env={**os.environ, "R_LIBS_USER": os.environ.get(
                               "R_LIBS_USER", str(Path.home() / "R" / "library"))})
    except Exception:
        return {}
    out = {}
    for ln in (p.stdout or "").splitlines():
        if "\t" in ln:
            k, v = ln.split("\t", 1)
            if v.strip():
                out[k.strip()] = v.strip()
    return out


# ── databases ─────────────────────────────────────────────────────────────────
def _load_manifest() -> dict:
    try:
        return json.loads(MANIFEST.read_text())
    except Exception:
        return {}


def _save_manifest(m: dict):
    try:
        DB_DIR.mkdir(parents=True, exist_ok=True)
        tmp = MANIFEST.with_suffix(".json.part")
        tmp.write_text(json.dumps(m, indent=2, ensure_ascii=False))
        tmp.replace(MANIFEST)
    except Exception:
        pass


def md5_file(p: Path, chunk=8 << 20) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def db_label(path: str) -> str:
    """A human name for a reference file, read from the file name."""
    n = Path(path).name
    low = n.lower()
    m = re.search(r"silva.*?v?(\d{3}(?:\.\d)?)", low)
    if "silva" in low and m:
        return f"SILVA {m.group(1)}"
    if "unite" in low or low.startswith("sh_general_release"):
        d = re.search(r"(\d{2})\.(\d{2})\.(\d{4})", n)
        return f"UNITE general release {d.group(3)}-{d.group(2)}-{d.group(1)}" if d else "UNITE"
    m = re.search(r"pr2_version_(\d+\.\d+\.\d+)", low)
    if m:
        return f"PR2 {m.group(1)}"
    m = re.search(r"midori2_.*?gb(\d+)_([a-z0-9]+)", low)
    if m:
        gene = {"co1": "COI", "srrna": "12S srRNA", "lrrna": "16S lrRNA",
                "cytb": "Cytb"}.get(m.group(2), m.group(2).upper())
        return f"MIDORI2 GenBank {m.group(1)} ({gene})"
    m = re.search(r"gtdb.*?r(\d+)", low)
    if m:
        return f"GTDB R{m.group(1)}"
    if "nemabase" in low:
        return "18S-NemaBase"
    if "emu" in low or (Path(path).is_dir() and (Path(path) / "taxonomy.tsv").exists()):
        return f"Emu database ({Path(path).name if Path(path).is_dir() else Path(path).parent.name})"
    return n


def db_record(role: str, path: str) -> dict | None:
    if not path:
        return None
    p = Path(path).expanduser()
    rec = {"role": role, "path": str(p), "file": p.name, "label": db_label(str(p))}
    if p.is_dir():
        # an Emu database is a folder; its sequences are what identify the build
        for f in ("species_taxid.fasta", "sequences.fasta"):
            if (p / f).exists():
                p = p / f
                rec["file"] = f"{rec['file']}/{f}"
                break
        else:
            return rec
    if not p.is_file():
        rec["missing"] = True
        return rec
    st = p.stat()
    rec["bytes"] = st.st_size
    man = _load_manifest()
    ent = man.get(str(p)) or {}
    if ent.get("bytes") == st.st_size and ent.get("mtime") == int(st.st_mtime) and ent.get("md5"):
        rec["md5"] = ent["md5"]
    else:
        t0 = time.time()
        rec["md5"] = md5_file(p)
        ent.update({"bytes": st.st_size, "mtime": int(st.st_mtime), "md5": rec["md5"],
                    "hashed_at": datetime.now().isoformat(timespec="seconds"),
                    "hash_seconds": round(time.time() - t0, 1)})
        man[str(p)] = ent
        _save_manifest(man)
    for k in ("label", "version", "source", "md5_expected"):
        if ent.get(k):
            rec[k if k != "label" else "label"] = ent[k]
    if rec.get("md5_expected"):
        rec["md5_ok"] = rec["md5_expected"] == rec["md5"]
    return rec


# ── main entry ────────────────────────────────────────────────────────────────
def app_version() -> str:
    for d in (HERE, HERE.parent):
        try:
            return json.loads((d / "version.json").read_text())["version"]
        except Exception:
            continue
    return ""


def write(output_dir: str | Path, pipeline: str, databases: list[tuple[str, str]],
          extra_tools: dict | None = None) -> dict:
    """Collect and write provenance.json; never raises (a run must not fail on this)."""
    out = Path(output_dir)
    rec = {"app_version": app_version(), "pipeline": pipeline,
           "recorded_at": datetime.now().isoformat(timespec="seconds"),
           "tools": {}, "databases": []}
    try:
        rec["tools"].update(r_versions(R_PKGS.get(pipeline, ["dada2"])))
    except Exception:
        pass
    for t in CLI_TOOLS.get(pipeline, []):
        try:
            v = tool_version(t)
            if v:
                rec["tools"][t] = v
        except Exception:
            pass
    if extra_tools:
        rec["tools"].update({k: v for k, v in extra_tools.items() if v})
    seen = set()
    for role, path in databases:
        if not path or path in seen:
            continue
        seen.add(path)
        try:
            r = db_record(role, path)
            if r:
                rec["databases"].append(r)
        except Exception as e:
            rec["databases"].append({"role": role, "path": path, "error": str(e)})
    try:
        (out / "provenance.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    except Exception:
        pass
    return rec


def species_file_from_log(output_dir: Path) -> str:
    """The species-assignment file the run actually used, as the pipeline logged it."""
    lp = Path(output_dir) / "pipeline.log"
    if not lp.exists():
        return ""
    t = lp.read_text(errors="replace")
    for pat in (r"Species file\s*:\s*(\S+)", r"addSpecies.*?\s(/\S+\.fa(?:sta)?(?:\.gz)?)",
                r"species (?:file|reference)[^/\n]*(/\S+\.fa(?:sta)?(?:\.gz)?)"):
        m = re.search(pat, t, re.I)
        if m and Path(m.group(1)).exists():
            return m.group(1)
    return ""


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("result_dir")
    ap.add_argument("--pipeline", default="")
    ap.add_argument("--db", action="append", default=[], help="reference file used (repeatable)")
    a = ap.parse_args()
    rd = Path(a.result_dir)
    pipe = a.pipeline
    params = {}
    try:
        params = json.loads((rd / "run_params.json").read_text())
    except Exception:
        pass
    if not pipe:
        try:
            pipe = json.loads((rd / "summary.json").read_text()).get("pipeline", "")
        except Exception:
            pipe = ""
        pipe = (pipe or "dada2").lower()
    dbs = [("taxonomy", d) for d in a.db] or [("taxonomy", params.get("dbPath") or "")]
    print(json.dumps(write(rd, pipe, dbs), indent=2))
