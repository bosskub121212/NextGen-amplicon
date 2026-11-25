#!/usr/bin/env python3
"""Emu alignment speed: 'fast' runs minimap2 with Emu's options plus -f 5000 and gives Emu the SAM (same
output names); 'exact' leaves the alignment to Emu. Fake minimap2 / emu record what they were given."""
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
    print("Emu alignment speed (exact / fast)")
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
            for mode in ("fast", "exact"):
                (td / "calls.jsonl").write_text("")
                out = td / f"run_{mode}"
                out.mkdir()
                r = E.run_emu({"T1R3": fq}, out, str(db), 0.0001, 4, 0, align_speed=mode)
                calls = [json.loads(x) for x in (td / "calls.jsonl").read_text().splitlines()]
                tools = [c[0] for c in calls]
                left = sorted(p.name for p in (out / "emu_output" / "T1R3").iterdir())
                if mode == "fast":
                    mm, em = calls[0], calls[-1]
                    check("fast: minimap2 with Emu's options + -f 5000, then Emu on the SAM",
                          tools == ["minimap2", "emu"] and mm[mm.index("-f") + 1] == "5000"
                          and all(x in mm for x in ("-ax", "map-ont", "-N", "50", "-p", ".9"))
                          and mm[-4] == str(db / "species_taxid.fasta") and em[-1].endswith(".sam")
                          and em[em.index("--output-basename") + 1] == "T1R3.fastq", calls)
                else:
                    check("exact: Emu aligns itself (no minimap2 call, no -f)",
                          tools == ["emu"] and calls[0][-1] == str(fq) and "--output-basename" not in calls[0], calls)
                check(f"{mode}: Emu's usual output name, SAM removed",
                      r.get("T1R3") is not None and r["T1R3"].name == "T1R3.fastq_rel-abundance.tsv"
                      and not any(x.endswith(".sam") for x in left), (r, left))
            sys.path.insert(0, str(ROOT / "tests"))
            import test_report_many_fast as TM
            import report_builder as RB
            d = TM.make_run(td / "ont", 3)
            sm = json.loads((d / "summary.json").read_text())
            sm.update({"emu_align_speed": "fast", "emu_minimap2_f": 5000})
            (d / "summary.json").write_text(json.dumps(sm))
            html = RB.build_profile_html(RB.Run(d, set()))
            check("report method: fast alignment said, with its minimap2 command",
                  "-f 5000" in html and "minimap2 -ax map-ont -N 50 -p 0.9 -f 5000" in html, html.count("-f 5000"))
            mm_ = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            check("form option ont_align_speed (default exact) reaches the pipeline",
                  'ont_align_speed: str = "exact"' in mm_ and '"--align_speed", "fast"' in mm_
                  and '["ont_align_speed",   "Emu alignment (minimap2)", "select", "exact"' in fe
                  and 'params.ont_align_speed = /^\\s*f/i.test(sp) ? "fast" : "exact";' in fe)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")
    finally:
        os.environ.clear()
        os.environ.update(old_env)


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
