#!/usr/bin/env python3
"""Job page result previews (result_preview.py, /jobs/{id}/preview) and the new job-page / form wiring."""
from __future__ import annotations
import json, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
KEYS = ["1.", "2.", "3.", "4.", "5.Taxonomy_Identification", "6.Diversity"]


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def _bad_notes(p):
    return [f"{s['title']}: {i.get('text')}" for s in p["sections"] for i in s["items"]
            if i.get("kind") == "note" and "Not available" in (i.get("text") or "")]


def one(label, d: Path, want_svg=True):
    import result_preview as RP
    p = RP.build(d, use_cache=False)
    secs = p.get("sections") or []
    ks = [s["key"] for s in secs]
    check(f"{label}: a section per folder of the result, in order (1.QC ... 6.Diversity)",
          len(ks) == 6 and all(k.startswith(w) for k, w in zip(ks, KEYS)) and ks[0] == "1.QC"
          and ks[1] in ("2.Filter", "2.Cutadapt"), ks)
    bad = _bad_notes(p)
    check(f"{label}: no section failed", not bad, bad[:3])
    kinds = {i["kind"] for s in secs for i in s["items"]}
    if want_svg:
        check(f"{label}: has figures and tables", "svg" in kinds and "table" in kinds, kinds)
    for s in secs:
        for i in s["items"]:
            if i["kind"] == "table":
                w = len(i["head"])
                if any(len(r) != w for r in i["rows"]):
                    check(f"{label}: table '{i.get('title')}' rows match its header", False)
                    return p
    json.dumps(p)                                    # the endpoint returns it as JSON
    # cache: second call reads it back, and the result is the same
    RP.build(d)
    check(f"{label}: cached and read back", (d / RP.CACHE).exists() and RP.build(d) == json.loads(json.dumps(p)))
    return p


def main():
    print("Result previews and the job page / form wiring")
    try:
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            from test_amp_fast import make_amp_run, make_emu_run
            one("DADA2 (12S)", make_amp_run(t / "amp"))
            p = one("Emu (ONT 16S)", make_emu_run(t / "emu"))
            s3 = next(s for s in p["sections"] if s["key"].startswith("3."))
            check("Emu: folders named 2.Filter, 3.Emu and 4.Taxa_table; a Krona picture in 5.",
                  s3["key"] == "3.Emu" and p["sections"][1]["key"] == "2.Filter"
                  and p["sections"][3]["key"] == "4.Taxa_table"
                  and any(i.get("title", "").startswith("Krona picture") for i in p["sections"][4]["items"]),
                  [s["key"] for s in p["sections"]])
            try:
                from test_its_ont_fast import make_vsearch_its_run
                one("VSEARCH (ITS ONT)", make_vsearch_its_run(t / "vs"))
            except ImportError as e:
                print("  skip  VSEARCH fixture:", e)
            from test_its_pacbio_pkg_fast import make_its, make_pacbio
            one("ITS Illumina", make_its(t / "its"))
            one("PacBio", make_pacbio(t / "pb", old=False))
            one("PacBio (old, one row per sample)", make_pacbio(t / "pbo", old=True))
            from test_report_many_fast import make_run
            one("many samples (60)", make_run(t / "many", 60))
        mm = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        a, b = mm.find('@app.get("/jobs/{job_id}/preview")'), mm.find('@app.get("/packages/{job_id}")')
        check("main.py: /jobs/{id}/preview, id checked, 409 while running, WGS/shotgun empty",
              0 < a < b and "_JOB_ID_RE.fullmatch(job_id)" in mm[a:b] and "status_code=409" in mm[a:b]
              and 'pipe == "shotgun"' in mm[a:b] and "rp.build(out)" in mm[a:b])
        fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
        check("frontend: pipeline card, settings beside the log, previews",
              'id="pipeBox"' in fe and "function paintPipe(" in fe and '<div class="two2' in fe
              and "function loadPreview(" in fe and '"/preview"' in fe)
        check("frontend: Run type beside the uploads; step bar for pipeline settings (amplicon + WGS)",
              'drop runtype' in fe and "function mountStepper(" in fe and 'mountStepper("advBody")' in fe
              and 'mountStepper("wgsBody")' in fe)
        check("frontend: the OTU similarity is sent with the run (it never was)", "p.otuSimilarity" in fe)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
