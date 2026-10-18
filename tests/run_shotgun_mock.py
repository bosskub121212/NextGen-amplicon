#!/usr/bin/env python3
"""Shotgun metagenome mock test (beta pipeline): simulated Illumina and Nanopore
reads from the 8 bacteria of the ZymoBIOMICS D6300 standard, run through
shotgun_pipeline.py exactly as the app runs it, compared with the truth.

    python3 tests/run_shotgun_mock.py                  # both platforms
    python3 tests/run_shotgun_mock.py illumina --keep ~/sg_mock

Genomes: the type strains in ~/.ngamp_wgs_refs (setup_typestrains.sh); a
missing one is downloaded. Kraken 2 index: db_paths.json kraken2_db
(setup_shotgun.sh) or --kraken_db. Truth is the READ share of each species
(cells x genome length), which is what Bracken estimates.

PASS when, for each platform: genus-level Bray-Curtis to the truth <= 0.15,
every species of >= 1 % is found at >= 25 % of its expected share, and species
outside the mock make up < 5 % of the reads. Reads simulated from the very
genomes an index may contain flatter the classifier; this test guards the
pipeline (sample handling, QC, Bracken, tables, report), not Kraken 2's
accuracy on real samples.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))

ZYMO = {"Pseudomonas aeruginosa": 12, "Escherichia coli": 12, "Salmonella enterica": 12,
        "Limosilactobacillus fermentum": 12, "Enterococcus faecalis": 12,
        "Staphylococcus aureus": 12, "Listeria monocytogenes": 12, "Bacillus subtilis": 12}
# cells are even in D6300; vary them so a pipeline that returns "even" cannot pass by accident
CELLS = [30, 20, 12, 10, 10, 8, 6, 4]


def genomes(cache: Path, fetch: bool = True) -> list[dict]:
    import wgs_typestrain_db as T
    want = list(ZYMO)
    found = {}
    for f in (cache / "type_strains").glob("*.v3.json"):
        for r in json.loads(f.read_text()):
            if r.get("organism") in want and r["organism"] not in found:
                g = cache / f"{r['accession']}.fna.gz"
                if g.exists():
                    found[r["organism"]] = g
    missing = [s for s in want if s not in found]
    if missing and fetch:
        gen = sorted({m.split()[0] for m in missing})
        print(f"  downloading type strains of {', '.join(gen)}")
        T.main(["build", "--genera", ",".join(gen), "--cache", str(cache)])
        return genomes(cache, fetch=False)
    if missing:
        return []
    return [{"name": s, "taxid": i, "fasta": str(found[s]), "cells": CELLS[i]} for i, s in enumerate(want)]


def bray(a, b):
    k = set(a) | set(b)
    return sum(abs(a.get(x, 0) - b.get(x, 0)) for x in k) / max(1e-9, sum(a.values()) + sum(b.values()))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("platforms", nargs="*", default=["illumina", "ont"])
    ap.add_argument("--kraken_db")
    ap.add_argument("--cache", default=str(Path.home() / ".ngamp_wgs_refs"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep")
    ap.add_argument("--pairs", type=int, default=100000)
    ap.add_argument("--ont_reads", type=int, default=10000)
    a = ap.parse_args(argv)
    db = a.kraken_db
    if not db:
        try:
            db = json.loads((BACKEND / "databases" / "db_paths.json").read_text()).get("kraken2_db", "")
        except Exception:
            db = ""
    if not db or not (Path(db) / "hash.k2d").exists():
        print("SKIP — no Kraken 2 index (bash ~/r16s-app/setup_shotgun.sh)")
        return 0
    G = genomes(Path(a.cache).expanduser())
    if len(G) < len(ZYMO):
        print("SKIP — type-strain genomes of the mock species not available (setup_typestrains.sh)")
        return 0
    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="ngamp_sg_mock_"))
    work.mkdir(parents=True, exist_ok=True)
    (work / "spec.json").write_text(json.dumps(G))
    results = []
    import report_builder as rb
    for plat in a.platforms:
        print(f"\n=== shotgun {plat}: 8 Zymo species, uneven cells")
        inp, out = work / plat / "input", work / plat / "result"
        shutil.rmtree(work / plat, ignore_errors=True)
        inp.mkdir(parents=True)
        n = a.pairs if plat == "illumina" else a.ont_reads
        subprocess.run([sys.executable, str(ROOT / "tests" / "mock" / "simulate_shotgun.py"),
                        "--genomes", str(work / "spec.json"), "--platform", plat, "--reads", str(n),
                        "--out", str(inp), "--sample", "mock1", "--seed", "11"], check=True)
        truth_ids = json.loads((inp / "mock1_truth.json").read_text())["reads"]
        truth = {G[int(k)]["name"]: v for k, v in truth_ids.items()}
        (inp / "mock1_truth.json").rename(work / plat / "truth.json")
        t0 = time.time()
        p = subprocess.run([sys.executable, str(BACKEND / "shotgun_pipeline.py"), "--input_dir", str(inp),
                            "--output_dir", str(out), "--platform", plat, "--kraken_db", db,
                            "--threads", str(a.threads)], capture_output=True, text=True)
        (out / "pipeline.log").write_text(p.stdout + p.stderr)
        took = time.time() - t0
        if p.returncode != 0:
            print("  FAIL — pipeline exited", p.returncode, "\n   ", "\n    ".join(p.stdout.splitlines()[-10:]))
            results.append((plat, "FAIL", "pipeline"))
            continue
        tax = {r["tax_id"]: r for r in csv.DictReader(open(out / "taxonomy.csv"))}
        obs_sp, obs_g = {}, {}
        rows = list(csv.DictReader(open(out / "asv_table.csv")))
        tot = sum(float(r["mock1"]) for r in rows) or 1
        for r in rows:
            v = float(r["mock1"]) / tot
            t = tax.get(r["tax_id"], {})
            obs_sp[t.get("Species", "")] = obs_sp.get(t.get("Species", ""), 0) + v
            obs_g[t.get("Genus", "")] = obs_g.get(t.get("Genus", ""), 0) + v
        tg = {}
        for s, v in truth.items():
            tg[s.split()[0]] = tg.get(s.split()[0], 0) + v
        bc = bray(tg, obs_g)
        missing = [s for s, v in truth.items() if v >= 0.01 and obs_sp.get(s, 0) < 0.25 * v]
        extra = sum(v for s, v in obs_sp.items() if s not in truth)
        for s, v in truth.items():
            print(f"      {s:32s} expected {100 * v:6.2f}%  observed {100 * obs_sp.get(s, 0):6.2f}%")
        print(f"  genus Bray-Curtis {bc:.3f} · species outside the mock {100 * extra:.2f}%"
              + (f" · MISSING {missing}" if missing else ""))
        try:
            (out / "run_params.json").write_text(json.dumps({"marker": "SHOTGUN", "sequencerType": plat}))
            html = rb.build_profile_html(rb.Run(out, set()), title=f"Shotgun mock — {plat}")
            (out / "report.html").write_text(html)
            rep_ok = True
        except Exception as e:
            rep_ok = False
            print(f"  report failed: {e!r}")
        ok = bc <= 0.15 and not missing and extra < 0.05 and rep_ok
        print(f"  {'PASS' if ok else 'FAIL'} in {took:.0f}s")
        results.append((plat, "PASS" if ok else "FAIL", f"BC {bc:.3f}, outside {100 * extra:.1f}%"))
    print("\n=== Summary")
    for n_, st, msg in results:
        print(f"  {st:4s}  shotgun-{n_:9s} {msg}")
    print(f"  work folder: {work}" + ("" if a.keep else " (removed)"))
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if any(st == "FAIL" for _, st, _ in results) else 0


if __name__ == "__main__":
    sys.exit(main())
