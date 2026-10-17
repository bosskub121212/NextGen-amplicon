#!/usr/bin/env python3
"""Mock-community tests: every amplicon pipeline, end to end, against a known truth.

For each pipeline this simulates reads from a mock community with a known
composition (tests/mock/simulate.py), runs the pipeline script with the same
arguments main.py gives it, records provenance, then checks the result with
tests/mock/check_mock.py (composition, missing / invented taxa, rank shifts) and
builds the HTML report from it (the report must build for every pipeline).

    python3 tests/run_mock_tests.py                 # every case whose database is installed
    python3 tests/run_mock_tests.py dada2 its       # just these
    python3 tests/run_mock_tests.py --keep /tmp/mt  # keep the work folder

Databases come from backend/databases/db_paths.json (SILVA_16S, UNITE_ITS2,
MIDORI2_COX1), or --silva/--unite/--midori. A case without its database is
SKIPPED, not failed. Exit status: 0 when nothing failed.

Takes 5-20 minutes depending on the machine (taxonomy assignment dominates).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BACKEND = ROOT / "backend"
RS = BACKEND / "r_scripts"
MOCK = ROOT / "tests" / "mock"
sys.path.insert(0, str(BACKEND))
sys.path.insert(0, str(MOCK))

V34_F, V34_R = "CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"
FL_F, FL_R = "AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT"
ITS_F, ITS_R = "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"        # ITS1F / ITS2
COI_F, COI_R = "GGWACWGGWTGAACWGTWTAYCCYCC", "TANACYTCNGGRTGNCCRAARAAYCA"  # mlCOIintF / jgHCO2198


def py():
    return sys.executable


CASES = {
    "dada2": dict(
        title="16S V3-V4, DADA2 (Illumina 2x300)", db="silva", marker="16S",
        sim=["--refs", "refs_16s.fasta", "--composition", "zymo_16s", "--primer_f", V34_F,
             "--primer_r", V34_R, "--platform", "illumina", "--read_len", "300"],
        params=dict(marker="16S", primer_f=V34_F, primer_r=V34_R, truncLen_F=260, truncLen_R=220,
                    maxEE_F=2.0, maxEE_R=2.0, pool="FALSE", minBoot=50, chimeraMethod="consensus"),
        cmd=lambda i, o, db, t: ["Rscript", str(RS / "dada2_pipeline.R"), "--input", i, "--output", o,
                                 "--marker", "16S", "--truncLen_F", "260", "--truncLen_R", "220",
                                 "--maxEE_F", "2", "--maxEE_R", "2", "--pool", "FALSE",
                                 "--chimeraMethod", "consensus", "--dbPath", db, "--minBoot", "50",
                                 "--threads", str(t), "--primer_f", V34_F, "--primer_r", V34_R,
                                 "--discard_untrimmed", "TRUE"],
        kind="dada2"),
    "vsearch": dict(
        title="16S V3-V4, VSEARCH OTUs (Illumina pairs merged)", db="silva", marker="16S",
        sim=["--refs", "refs_16s.fasta", "--composition", "zymo_16s", "--primer_f", V34_F,
             "--primer_r", V34_R, "--platform", "illumina", "--read_len", "300"],
        params=dict(marker="16S", sequencerType="qiime2_vsearch", primer_f=V34_F, primer_r=V34_R,
                    otuUnoiseAlpha=2.0, otuUnoiseMinsize=8),
        cmd=lambda i, o, db, t: [py(), str(BACKEND / "qiime2_vsearch_pipeline.py"), "--input", i,
                                 "--output", o, "--db_path", db, "--threads", str(t),
                                 "--primer_f", V34_F, "--primer_r", V34_R, "--marker", "16S",
                                 "--cutadapt_error_rate", "0.1", "--cutadapt_overlap", "10"],
        kind="qiime2_vsearch"),
    "pacbio": dict(
        title="16S full length, PacBio HiFi (both orientations)", db="silva", marker="PACBIO",
        sim=["--refs", "refs_16s.fasta", "--composition", "zymo_16s", "--primer_f", FL_F,
             "--primer_r", FL_R, "--platform", "pacbio", "--reads", "2000"],
        params=dict(marker="PACBIO", primer_f=FL_F, primer_r=FL_R, pb_min_len=1000, pb_max_len=1600,
                    pb_maxEE=3.0, pool="pseudo", minBoot=50),
        cmd=lambda i, o, db, t: ["Rscript", str(RS / "pacbio_pipeline.R"), "--input_dir", i,
                                 "--output_dir", o, "--min_length", "1000", "--max_length", "1600",
                                 "--maxEE", "3", "--threads", str(t), "--region", "V1-V9",
                                 "--pool", "pseudo", "--dbPath", db, "--minBoot", "50",
                                 "--primer_f", FL_F, "--primer_r", FL_R],
        kind="pacbio"),
    "its": dict(
        title="ITS1, DADA2 + UNITE (Illumina 2x250)", db="unite", marker="ITS1",
        sim=["--refs", "refs_its.fasta", "--composition", "fungi_its", "--primer_f", ITS_F,
             "--primer_r", ITS_R, "--platform", "illumina", "--read_len", "250"],
        params=dict(marker="ITS1", primer_f=ITS_F, primer_r=ITS_R, maxEE_F=2.0, maxEE_R=2.0, minBoot=50),
        cmd=lambda i, o, db, t: ["Rscript", str(RS / "its_pipeline.R"), "--input_dir", i,
                                 "--output_dir", o, "--its_region", "ITS1", "--threads", str(t),
                                 "--primer_f", ITS_F, "--primer_r", ITS_R, "--dbPath", db,
                                 "--minBoot", "50", "--discard_untrimmed", "TRUE"],
        kind="its"),
    "cox1": dict(
        title="COI (Leray fragment), DADA2 + MIDORI2 (Illumina 2x250)", db="midori", marker="COX1",
        sim=["--refs", "refs_coi.fasta", "--composition", "insects_coi", "--primer_f", COI_F,
             "--primer_r", COI_R, "--platform", "illumina", "--read_len", "250"],
        params=dict(marker="COX1", primer_f=COI_F, primer_r=COI_R, truncLen_cox1_f=200,
                    truncLen_cox1_r=180, maxEE_F=2.0, maxEE_R=5.0, codon_table=5, cox1_min_len=290,
                    cox1_max_len=330, run_lulu=True, minBoot=50),
        cmd=lambda i, o, db, t: ["Rscript", str(RS / "cox1_pipeline.R"), "--input_dir", i,
                                 "--output_dir", o, "--truncLen_f", "200", "--truncLen_r", "180",
                                 "--maxEE_f", "2", "--maxEE_r", "5", "--codon_table", "5",
                                 "--min_length", "290", "--max_length", "330", "--threads", str(t),
                                 "--lulu", "TRUE", "--primer_f", COI_F, "--primer_r", COI_R,
                                 "--dbPath", db, "--minBoot", "50"],
        kind="cox1"),
}


def find_dbs(a) -> dict:
    dp = {}
    try:
        dp = json.loads((BACKEND / "databases" / "db_paths.json").read_text())
    except Exception:
        pass

    def pick(explicit, *keys):
        if explicit:
            return explicit
        for k in keys:
            v = dp.get(k)
            if v and Path(v).exists():
                return v
        return ""
    return {"silva": pick(a.silva, "SILVA_16S", "SILVA_144", "SILVA_138_2"),
            "unite": pick(a.unite, "UNITE_ITS1", "UNITE_ITS2"),
            "midori": pick(a.midori, "MIDORI2_COX1")}


def run_pipeline(cmd, out: Path, log: Path, timeout: int) -> int:
    """Run like main.py does, answering a DADA2 checkpoint with 'continue'."""
    with open(log, "w") as lf:
        p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        stop = threading.Event()

        def answer():
            while not stop.is_set():
                if (out / "checkpoint.json").exists() and not (out / "checkpoint_signal").exists():
                    try:
                        if "Waiting for user decision" in log.read_text(errors="replace")[-4000:]:
                            (out / "checkpoint_signal").write_text("continue\n")
                    except Exception:
                        pass
                time.sleep(2)
        th = threading.Thread(target=answer, daemon=True)
        th.start()
        try:
            rc = p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            rc = -9
        stop.set()
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases", nargs="*", default=list(CASES))
    ap.add_argument("--silva")
    ap.add_argument("--unite")
    ap.add_argument("--midori")
    ap.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 2)))
    ap.add_argument("--keep", help="work folder to keep (default: a temporary one, removed)")
    ap.add_argument("--timeout", type=int, default=3600)
    a = ap.parse_args(argv)

    import check_mock
    import provenance
    import report_builder as rb

    dbs = find_dbs(a)
    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="ngamp_mock_"))
    work.mkdir(parents=True, exist_ok=True)
    results = []
    for name in a.cases:
        c = CASES.get(name)
        if not c:
            print(f"unknown case {name}; choose from {', '.join(CASES)}")
            return 2
        db = dbs.get(c["db"])
        print(f"\n=== {name}: {c['title']}")
        if not db:
            print(f"  SKIP — no {c['db']} database (pass --{c['db']} or install it)")
            results.append((name, "SKIP", f"no {c['db']} database"))
            continue
        inp, out = work / name / "input", work / name / "result"
        shutil.rmtree(work / name, ignore_errors=True)
        inp.mkdir(parents=True)
        out.mkdir(parents=True)
        subprocess.run([py(), str(MOCK / "simulate.py"), *c["sim"], "--out", str(inp)], check=True,
                       cwd=MOCK, stdout=subprocess.DEVNULL)
        truth = json.loads((inp / "truth.json").read_text())
        (out / "run_params.json").write_text(json.dumps({**c["params"], "dbPath": db}, indent=2))
        t0 = time.time()
        rc = run_pipeline(c["cmd"](str(inp), str(out), db, a.threads), out, out / "pipeline.log",
                          a.timeout)
        took = time.time() - t0
        if rc != 0:
            tail = (out / "pipeline.log").read_text(errors="replace").strip().splitlines()[-15:]
            print(f"  FAIL — pipeline exited {rc} after {took:.0f}s\n    " + "\n    ".join(tail))
            results.append((name, "FAIL", f"pipeline exit {rc}"))
            continue
        provenance.write(out, c["kind"], [("taxonomy", db)])
        rep = check_mock.check(out, truth)
        for s, r in rep["samples"].items():
            print(f"  {s}: Bray-Curtis family {r['bray_curtis_family']:.3f} genus "
                  f"{r['bray_curtis_genus']:.3f}  unexpected {r['unexpected_pct']:.2f}%  "
                  f"unclassified {r['unclassified_pct']:.2f}%"
                  + (f"  family only {r['family_only']}" if r["family_only"] else "")
                  + (f"  MISSING {r['missing']}" if r["missing"] else ""))
        first = next(iter(rep["samples"].values()), {})
        for g, e in first.get("expected", {}).items():
            print(f"      {g:22s} expected {e:6.2f}%  observed {first['observed'].get(g, 0):6.2f}%")
        try:
            run = rb.Run(out, set())
            html = rb.build_profile_html(run, title=f"Mock test — {c['title']}")
            (out / "report.html").write_text(html)
            rep_ok = True
        except Exception as e:
            rep_ok = False
            rep["problems"].append(f"report failed to build: {e!r}")
        status = "PASS" if rep["pass"] and rep_ok else "FAIL"
        print(f"  {status} in {took:.0f}s" + ("" if status == "PASS" else
                                              "\n    " + "\n    ".join(rep["problems"])))
        (out / "mock_check.json").write_text(json.dumps(rep, indent=2))
        results.append((name, status, "; ".join(rep["problems"]) or f"{took:.0f}s"))

    print("\n=== Summary")
    for n, st, msg in results:
        print(f"  {st:4s}  {n:8s}  {msg}")
    print(f"  work folder: {work}" + ("" if a.keep else " (removed)"))
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if any(st == "FAIL" for _n, st, _m in results) else 0


if __name__ == "__main__":
    sys.exit(main())
