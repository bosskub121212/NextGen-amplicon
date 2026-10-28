#!/usr/bin/env python3
"""End-to-end check of shotgun_pipeline.py --func with stand-in tools (no databases).

    python3 tests/run_func_e2e.py            # Illumina and Nanopore

Kraken 2, MetaPhlAn, HUMAnN and the humann_* utilities are tests/fake_func/fakefunc.py;
fastp, Bracken, hostile and the resistome are left out. Checks that the function step
runs inside the real pipeline, writes its tables, chart tables and summary entry, and
that the run still succeeds. PASS / FAIL per platform.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
import test_func_fast as TF                                   # noqa: E402


def one(platform: str, tmp: Path) -> list[str]:
    dbp, bases = TF.fake_world(tmp / "world")
    fbin = bases[0] / "envs" / "ngamp-func" / "bin"
    meta_bin = tmp / "metabin"
    TF._exe(meta_bin / "kraken2", f'#!/bin/sh\nexec python3 "{TF.FAKE}" kraken2 "$@"\n')
    kdb = tmp / "k2db"
    kdb.mkdir()
    (kdb / "hash.k2d").write_text("x")
    dbj = tmp / "db_paths.json"
    dbj.write_text(json.dumps(dbp))
    inp = tmp / "in"
    for s in ("S1", "S2", "S3"):
        TF.fake_reads(inp, s, n=20, paired=(platform == "illumina"))
    if platform == "ont":
        for f in inp.glob("*_R1.fastq"):
            f.rename(inp / f.name.replace("_R1", ""))
    meta = tmp / "metadata.csv"
    meta.write_text("sampleId,group\nS1,a\nS2,a\nS3,b\n")
    out = tmp / "out"
    env = dict(os.environ, NGAMP_META_BIN=str(meta_bin), PATH=f"{fbin}{os.pathsep}{os.environ['PATH']}")
    cmd = [sys.executable, str(ROOT / "backend" / "shotgun_pipeline.py"), "--input_dir", str(inp),
           "--output_dir", str(out), "--platform", platform, "--kraken_db", str(kdb), "--threads", "2",
           "--db_paths", str(dbj), "--no_amr", "--func", "--func_panels", "gut,food", "--metadata", str(meta)]
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    (tmp / "pipeline.log").write_text(p.stdout + p.stderr)
    bad = []
    if p.returncode != 0:
        bad.append(f"pipeline exit {p.returncode}: {(p.stdout + p.stderr)[-600:]}")
    for f in ("func/pathways.tsv", "func/ko_unstratified.tsv", "func/panels.csv", "func/func_qc.csv",
              "func/mpa_vs_kraken.csv", "r_tables/func_alpha.csv", "r_tables/func_pcoa_ko.csv",
              "asv_table.csv"):
        if not (out / f).exists():
            bad.append(f"missing {f}")
    try:
        sm = json.loads((out / "summary.json").read_text())
        fi = sm.get("func") or {}
        if fi.get("n_pathways") != 2 or fi.get("failed"):
            bad.append(f"summary func: {fi}")
        if "permanova" not in (fi.get("stats") or {}):
            bad.append("summary func.stats.permanova missing")
        if "humann" not in (sm.get("tools") or {}):
            bad.append("summary tools lacks humann")
    except Exception as e:
        bad.append(f"summary.json: {e}")
    return bad


def main() -> int:
    fails = 0
    for plat in ("illumina", "ont"):
        with tempfile.TemporaryDirectory() as td:
            bad = one(plat, Path(td))
            print(f"{'PASS' if not bad else 'FAIL'}  {plat}")
            for b in bad:
                print("      " + b)
            fails += bool(bad)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
