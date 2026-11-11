#!/usr/bin/env python3
"""Fast checks for resuming an ONT 16S (Emu) run that was cut off part-way
(server restart, Windows Update): finished samples are reused, only the rest run,
nothing in the results folder is deleted.

Run alone (python3 tests/test_resume_fast.py) or from tests/test_fast.py.
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAILS = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILS.append(name)


FAKE_EMU = r'''#!/usr/bin/env python3
import sys, os
a = sys.argv; o = a[a.index("--output-dir") + 1]; fq = a[-1]
open(os.environ["EMULOG"], "a").write(os.path.basename(fq) + "\n")
if os.environ.get("EMUFAIL") and os.environ["EMUFAIL"] in os.path.basename(fq): sys.exit(3)
open(os.path.join(o, "x.sam"), "w").write("sam")
open(os.path.join(o, "x_rel-abundance-threshold-0.0001.tsv"), "w").write(
  "tax_id\tabundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\testimated counts\n"
  "1\t1.0\tE coli\tEscherichia\tEnterobacteriaceae\tEnterobacterales\tGammaproteobacteria\tPseudomonadota\tBacteria\t10\n")
'''
FAKE_CUTADAPT = '''#!/usr/bin/env python3
import sys, shutil
a = sys.argv; shutil.copy(a[-1], a[a.index("-o") + 1])
'''


def _emu_resume():
    print("Emu resume")
    with tempfile.TemporaryDirectory() as td:
        T = Path(td)
        b = T / "bin"; b.mkdir()
        (b / "chopper").write_text("#!/bin/sh\ncat\n")
        (b / "cutadapt").write_text(FAKE_CUTADAPT)
        (b / "emu").write_text(FAKE_EMU)
        for f in b.iterdir():
            f.chmod(0o755)
        inp = T / "in"; inp.mkdir()
        for s in "ABCD":
            with gzip.open(inp / f"{s}.fastq.gz", "wt") as fh:
                fh.write("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        db = T / "db"; db.mkdir()
        (db / "taxonomy.tsv").write_text("tax_id\tspecies\tgenus\n")
        (db / "species_taxid.fasta").write_text("")
        calls = T / "calls.txt"
        env = dict(os.environ, PATH=f"{b}{os.pathsep}{os.environ['PATH']}", EMULOG=str(calls),
                   HOME=str(T))            # no conda emu env under this HOME
        cmd = [sys.executable, str(ROOT / "backend" / "emu_pipeline.py"), "--input", str(inp),
               "--output", str(T / "out"), "--db_path", str(db), "--primer_f", "ACGT",
               "--primer_r", "ACGT", "--threads", "2"]
        subprocess.run(cmd, env=dict(env, EMUFAIL="C"), capture_output=True, text=True)
        first = calls.read_text().split()
        check("first run: Emu ran all 4 samples, C failed", sorted(first) == [f"{s}.fastq.gz" for s in "ABCD"],
              str(first))
        eo = T / "out" / "emu_output"
        (eo / "A" / "left.sam").write_text("big")
        calls.write_text("")
        r = subprocess.run(cmd + ["--resume"], env=env, capture_output=True, text=True)
        second = calls.read_text().split()
        check("resume: Emu runs only the unfinished sample", second == ["C.fastq.gz"], str(second))
        check("resume: says how many were reused", "3/4 sample(s) already finished" in r.stdout, r.stdout[-400:])
        check("resume: QC and trimming reused", "QC filter: 4 sample(s) reused" in r.stdout
              and "Primer trimming: 4 sample(s) reused" in r.stdout)
        check("resume: leftover SAM of a finished sample removed", not (eo / "A" / "left.sam").exists())
        head = (T / "out" / "asv_table.csv").read_text().splitlines()[0] if (T / "out" / "asv_table.csv").exists() else ""
        check("resume: table holds all 4 samples", head == ",A,B,C,D", head)
        check("no temporary .part- files left", not list((T / "out").rglob(".part-*")))


def _sources():
    print("Resume wiring")
    m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
    i = m.find('@app.post("/jobs/{job_id}/resume")')
    body = m[i:m.find("@app.post", i + 10)]
    check("POST /jobs/{id}/resume exists", i > 0)
    check("resume never deletes results", "rmtree" not in body and "unlink" not in body)
    check("resume only for stopped ONT 16S runs", '("error", "cancelled")' in body and "RESUMABLE_MARKERS" in body)
    check("emu gets --resume", 'if _resume:\n            cmd += ["--resume"]' in m)
    check("resumed run appends to pipeline.log", '"a" if _resume else "w"' in m)
    h = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
    check("Job Detail offers Resume run for a stopped ONT 16S run",
          'id="emuResume"' in h and "function emuResume" in h and '"/resume"' in h)


def main():
    for f in (_emu_resume, _sources):
        try:
            f()
        except Exception:
            traceback.print_exc()
            FAILS.append(f.__name__)
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
