#!/usr/bin/env python3
"""LULU for COI (installer, the pipeline's call) and the default CPU threads from the machine."""
from __future__ import annotations
import sys, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def main():
    print("LULU (COI) and default CPU threads")
    try:
        ext = (ROOT / "install_extensions.sh").read_text(encoding="utf-8")
        check("install_extensions.sh installs GitHub packages by owner/name (LULU, FUNGuildR never installed)",
              'remotes::install_github(repo' in ext and 'repo="tobiasgf/lulu"' in ext
              and 'repo="brendanf/FUNGuildR"' in ext and 'source="github", pkg=' not in ext)
        il = ROOT / "install_lulu.sh"
        t = il.read_text(encoding="utf-8") if il.exists() else ""
        check("install_lulu.sh: LULU from GitHub into the pipelines' R, BLAST+ checked",
              'install_github("tobiasgf/lulu"' in t and "makeblastdb" in t and "blastn" in t)
        r = (ROOT / "backend" / "r_scripts" / "cox1_pipeline.R").read_text(encoding="utf-8")
        check("cox1_pipeline.R: lulu() given a data.frame (a matrix fails inside lulu)",
              "as.data.frame(t(seqtab_translated))" in r)
        check("cox1_pipeline.R: the merged ASVs' reads go to their parent (LULU's curated table), log kept in LULU/",
              "curated$curated_table" in r and "stopifnot(sum(fin) == sum(seqtab_translated))" in r
              and "setwd(lulu_dir)" in r and "lulu_otu_map.csv" in r)
        check("cox1_pipeline.R: BLAST+ looked for in the app's conda envs too; what LULU did is recorded "
              "(ran / no_blast / empty_matchlist / error / missing / off)",
              'blast_bin <- function(tool)' in r and 'summary_list$lulu_status <- lulu_status' in r
              and 'lulu_status <- "no_blast"' in r and 'summary_list$lulu       <- identical(lulu_status, "ran")' in r)
        import types, report_builder as RB
        st = lambda sm: RB.Run.lulu_state.fget(types.SimpleNamespace(summary=sm, params={}, path=ROOT / "nope"))
        check("report: a LULU that did not run is not 'curated with LULU (0 merged)'",
              st({"lulu_requested": True, "lulu": False, "lulu_status": "no_blast"}) == "no_blast"
              and st({"lulu_requested": True, "lulu": True, "lulu_status": "ran"}) == "ran"
              and st({"lulu_requested": True, "lulu": True}) == "ran"
              and st({"lulu_requested": True, "lulu": False}) == "missing")
        mm = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
        check("default CPU threads = the CPUs the backend sees (/version cpus; forms and Run again)",
              '"cpus":         os.cpu_count() or 0' in mm and 'if(key === "nThreads" && MACHINE.cpus > 0) def = String(MACHINE.cpus);' in fe
              and "if(v && v.cpus > 0) MACHINE.cpus = v.cpus;" in fe and "MACHINE.cpus > 0 ? MACHINE.cpus" in fe)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
