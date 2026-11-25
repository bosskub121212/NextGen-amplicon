#!/usr/bin/env python3
"""Emu scheduling (v2.38.3): each sample aligned on ALL threads with Emu's own minimap2 command, one at a
time, and Emu given the SAM (its EM in the background while the next sample aligns) — same SAM, same
results, no idle cores during the single-core EM. Without minimap2 / the reference Emu aligns itself.
The fast alignment (-f 5000) of v2.38.0 is gone. Fake minimap2 / emu record what they were given."""
from __future__ import annotations
import json, os, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


FAKE_MM = """#!/usr/bin/env python3
import sys, json, os
open(os.environ["FAKE_LOG"], "a").write(json.dumps(["minimap2"] + sys.argv[1:]) + "\\n")
open(sys.argv[sys.argv.index("-o") + 1], "w").write("@HD\\tVN:1.6\\n")
"""
FAKE_EMU = """#!/usr/bin/env python3
import sys, json, os
a = sys.argv[1:]
open(os.environ["FAKE_LOG"], "a").write(json.dumps(["emu"] + a) + "\\n")
out = a[a.index("--output-dir") + 1]
inp = a[-1]
base = a[a.index("--output-basename") + 1] if "--output-basename" in a else os.path.basename(inp).rsplit(".", 1)[0]
open(os.path.join(out, base + "_rel-abundance.tsv"), "w").write("tax_id\\tabundance\\testimated counts\\n573\\t1.0\\t10\\n")
"""


def main():
    print("Emu: align on all threads, EM in the background (same results)")
    old_env = dict(os.environ)
    try:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            b = td / "bin"
            b.mkdir()
            for n, src in (("minimap2", FAKE_MM), ("emu", FAKE_EMU)):
                (b / n).write_text(src)
                (b / n).chmod(0o755)
            os.environ["PATH"] = f"{b}{os.pathsep}{os.environ.get('PATH', '')}"
            os.environ["HOME"] = str(td)                       # no ~/miniconda3/envs/emu here
            os.environ["FAKE_LOG"] = str(td / "calls.jsonl")
            db = td / "emu_db_mar2026"
            db.mkdir()
            (db / "species_taxid.fasta").write_text(">1\nACGT\n")
            (db / "taxonomy.tsv").write_text("tax_id\tspecies\n")
            fq = td / "T1R3.fastq.gz"
            fq.write_bytes(b"")
            import emu_pipeline as E
            fqs = {}
            for n in ("T1R3", "T4R2", "C3"):
                fqs[n] = td / f"{n}.fastq.gz"
                fqs[n].write_bytes(b"")
            (td / "calls.jsonl").write_text("")
            out = td / "run_piped"
            out.mkdir()
            r = E.run_emu(fqs, out, str(db), 0.0001, 8, 0)
            calls = [json.loads(x) for x in (td / "calls.jsonl").read_text().splitlines()]
            mm = [c for c in calls if c[0] == "minimap2"]
            em = [c for c in calls if c[0] == "emu"]
            check("minimap2 per sample on ALL threads with Emu's exact command (no -f), then Emu on the SAM",
                  len(mm) == 3 and all(c[c.index("-t") + 1] == "8" and "-f" not in c
                                        and c[1:3] == ["-ax", "map-ont"] and c[c.index("-N") + 1] == "50"
                                        and c[c.index("-p") + 1] == ".9" and c[c.index("-K") + 1] == "500000000"
                                        and c[-4] == str(db / "species_taxid.fasta") for c in mm)
                  and len(em) == 3 and all(c[-1].endswith(".sam") and "--output-basename" in c for c in em),
                  calls)
            left = [x.name for n in fqs for x in (out / "emu_output" / n).iterdir()]
            check("every sample: Emu's usual output name, SAM removed",
                  sorted(v.name for v in r.values()) == sorted(f"{n}.fastq_rel-abundance.tsv" for n in fqs)
                  and not any(x.endswith(".sam") for x in left), (r, left))
            (db / "species_taxid.fasta").unlink()               # no reference beside Emu: Emu aligns itself
            (td / "calls.jsonl").write_text("")
            out = td / "run_self"
            out.mkdir()
            r = E.run_emu({"T1R3": fq}, out, str(db), 0.0001, 4, 0)
            calls = [json.loads(x) for x in (td / "calls.jsonl").read_text().splitlines()]
            check("without minimap2 / reference: Emu aligns itself from the reads",
                  [c[0] for c in calls] == ["emu"] and calls[0][-1] == str(fq) and r.get("T1R3") is not None, calls)
            sys.path.insert(0, str(ROOT / "tests"))
            import test_report_many_fast as TM
            import report_builder as RB
            d = TM.make_run(td / "ont", 3)
            sm = json.loads((d / "summary.json").read_text())
            sm.update({"emu_align_speed": "fast", "emu_minimap2_f": 5000})
            (d / "summary.json").write_text(json.dumps(sm))
            html = RB.build_profile_html(RB.Run(d, set()))
            check("report of a v2.38.0 fast test run still says so, with its minimap2 command",
                  "-f 5000" in html and "minimap2 -ax map-ont -N 50 -p 0.9 -f 5000" in html, html.count("-f 5000"))
            mm_ = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            ep = (ROOT / "backend" / "emu_pipeline.py").read_text(encoding="utf-8")
            check("the fast alignment option is gone (form, Run again, backend, pipeline)",
                  "ont_align_speed" not in mm_ and "align_speed" not in ep and "FAST_MID_OCC" not in ep
                  and '["ont_align_speed"' not in fe and "delete params.ont_align_speed;" in fe)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")
    finally:
        os.environ.clear()
        os.environ.update(old_env)


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
