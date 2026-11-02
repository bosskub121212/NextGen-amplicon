#!/usr/bin/env python3
"""Result packages for download — "Download results" in the console, every pipeline.

  build(run_dir, upload_dir, kind, samples, split, dest, job_label, ...) -> manifest

kind "wgs"      each sample → the eight-folder package with Report_<S>.pdf/.docx
                (wgs_customer_report.build_sample).
kind "amplicon" each sample → the numbered package 1.QcReads, 2.Taxonomy, 3.Diversity,
                4.Tables with Report_<S>.pdf/.docx (amp_package.build_sample); with
                include_reads the quality-filtered reads go in 1.QcReads/filtered_reads.

split=True  → one zip per sample (samples of different customers never share a zip).
split=False → one zip holding the chosen samples' packages side by side.

Zips are written into `dest` on disk (results/<job>/packages/<build>/) — never built
in memory — and are listed in the returned manifest {"files": [{"name", "path", "size",
"samples"}], "warnings": [...]}. Files that are already compressed are stored, the rest
deflated.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

STORED = {".gz", ".bam", ".bai", ".zip", ".png", ".jpg", ".jpeg", ".pdf", ".docx", ".xlsx", ".bz2", ".xz",
          ".zst", ".cram"}
WHOLE_JOB_SKIP = {"packages", "customer_report", "_work"}


def safe(name: str) -> str:
    # \w misses Thai vowel / tone marks (ก้อง and ก่อง would both become ก_อง): keep the Thai block
    return re.sub(r"[^\w\u0E00-\u0E7F.\-]+", "_", name).strip("_.") or "sample"


def detect_kind(run_dir: Path) -> str | None:
    run_dir = Path(run_dir)
    if (run_dir / "wgs_results.json").exists():
        return "wgs"
    if (run_dir / "asv_table.csv").exists():
        return "amplicon"
    return None


def list_samples(run_dir: Path, kind: str) -> list[str]:
    run_dir = Path(run_dir)
    if kind == "wgs":
        sd = run_dir / "samples"
        return sorted(p.name for p in sd.iterdir() if p.is_dir() and (p / "result.json").exists()) \
            if sd.is_dir() else []
    import amp_data as AD
    return list(AD.load_run(run_dir).samples)


def zip_dir(src: Path, zpath: Path, arc_root: str | None = None, skip_top=()) -> Path:
    """Zip a folder to disk (written to .part, then renamed)."""
    src, zpath = Path(src), Path(zpath)
    zpath.parent.mkdir(parents=True, exist_ok=True)
    import uuid
    part = zpath.with_name(f"{zpath.name}.{uuid.uuid4().hex[:8]}.part")      # two builds never share it
    root = arc_root if arc_root is not None else src.name
    with zipfile.ZipFile(part, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True,
                         strict_timestamps=False) as zf:
        for dp, dns, fns in os.walk(src):
            rel_dir = Path(dp).relative_to(src)
            dns[:] = [d for d in dns if not d.startswith(".")]      # .raw (TYGS job answers), caches
            if rel_dir == Path("."):
                dns[:] = [d for d in dns if d not in skip_top]
            dns.sort()
            for fn in sorted(fns):
                f = Path(dp) / fn
                if f.suffix == ".part" or not f.is_file():
                    continue
                arc = str(Path(root) / rel_dir / fn) if root else str(rel_dir / fn)
                ct = zipfile.ZIP_STORED if f.suffix.lower() in STORED else zipfile.ZIP_DEFLATED
                zf.write(f, arc, compress_type=ct)
    part.replace(zpath)
    return zpath


def zip_whole_job(run_dir: Path, zpath: Path) -> Path:
    """Every result file of the job (internal use), written to disk; packages excluded."""
    return zip_dir(run_dir, zpath, arc_root="", skip_top=WHOLE_JOB_SKIP)


def unique_names(samples: list[str], all_samples: list[str]) -> dict:
    """File-safe name per sample, never shared: "FS 04" and "FS_04" (or names that differ only
    in characters a file name cannot hold) get _1, _2 … by their place in the run."""
    count = {}
    for s in all_samples:
        count[safe(s)] = count.get(safe(s), 0) + 1
    return {s: safe(s) if count.get(safe(s), 0) <= 1 else f"{safe(s)}_{all_samples.index(s) + 1}" for s in samples}


def build(run_dir: Path, upload_dir: Path | None, kind: str, samples: list[str], split: bool, dest: Path,
          job_label: str = "results", include_reads: bool = False, tools: dict | None = None,
          progress=None, formats=("pdf", "docx")) -> dict:
    run_dir, dest = Path(run_dir), Path(dest)
    progress = progress or (lambda i, n, label: None)
    known = list_samples(run_dir, kind)
    sel = list(dict.fromkeys(samples or known))
    bad = [s for s in sel if s not in known]
    if bad:
        raise ValueError(f"not in this run: {', '.join(bad)}")
    if not sel:
        raise ValueError("no sample chosen")
    dest.mkdir(parents=True, exist_ok=True)
    work = dest / "_work"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    files, warnings = [], []
    label = safe(job_label)

    if kind in ("wgs", "amplicon"):
        n = len(sel)
        names = unique_names(sel, known)
        if kind == "wgs":
            import wgs_customer_report as R

            def one(s):
                return R.build_sample(run_dir, s, work, formats, name=names[s])
        else:
            import amp_package as AP
            kw = {"tools": tools} if tools is not None else {}

            def one(s):
                return AP.build_sample(run_dir, upload_dir, s, work, formats, name=names[s],
                                       include_reads=include_reads, **kw)
        for i, s in enumerate(sel):
            progress(i, n, f"Packaging {s} ({i + 1}/{n})")
            m = one(s)
            for k, v in m.items():
                if k.endswith("_error"):
                    warnings.append(f"{s}: {k[:-6]} — {v}")
            warnings += [f"{s}: {x}" for x in m.get("notes", [])]
            for lg in work.glob("*.plots.log"):                 # R log: kept beside the zips, not inside
                lg.replace(dest / lg.name)
            if split:
                z = zip_dir(Path(m["dir"]), dest / f"{Path(m['dir']).name}.zip")
                files.append({"name": z.name, "path": str(z), "size": z.stat().st_size, "samples": [s]})
                shutil.rmtree(m["dir"], ignore_errors=True)       # the zip is the deliverable
        if not split:
            progress(n, n + 1, "Writing one zip")
            z = zip_dir(work, dest / f"{label}_results.zip", arc_root=f"{label}_results")
            files.append({"name": z.name, "path": str(z), "size": z.stat().st_size, "samples": sel})
    else:
        raise ValueError(f"unknown result kind {kind!r}")
    shutil.rmtree(work, ignore_errors=True)
    for f in dest.glob("*.part"):
        f.unlink(missing_ok=True)
    man = {"kind": kind, "samples": sel, "split": split, "include_reads": include_reads, "files": files,
           "warnings": warnings}
    (dest / "manifest.json").write_text(json.dumps(man, indent=1))
    progress(1, 1, "Done")
    return man
