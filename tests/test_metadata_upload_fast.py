#!/usr/bin/env python3
"""Metadata prepared outside the app (BI's CSV / TSV / Excel) put into a finished job: the sample
column found by its values, names matched ignoring case and punctuation, a running-number column
dropped, the grouping to show first picked; the console's group manager and upload are wired."""
from __future__ import annotations
import io, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def main():
    print("Metadata upload and group manager")
    try:
        import metadata_upload as MU
        S = ["D1K1_C_H_R1", "D1K1_C_H_R2", "D1K2_W_HI_R1", "D1K5.1_G_H_R1"]
        rows = MU.read_table(("No.\tSample name\tStation\tHabitat\n1\td1k1-c-h-r1\tD1K1\tCoral\n"
                              "2\tD1K1-C-H-R2\tD1K1\tCoral\n3\td1k2_w_hi_r1\tD1K2\tWater\n4\tEXTRA_9\tD9\tCoral\n").encode(),
                             "meta.tsv")
        r = MU.check(rows, S)
        check("TSV: sample column found by its values, names matched ignoring case and '-'",
              r.get("ok") and r["id_col"] == "Sample name" and r["matched"] == 3, r)
        check("running number dropped; Station picked; unmatched both ways reported",
              r["columns"] == ["Station", "Habitat"] and r["pick"] == "Station"
              and r["unmatched_file"] == ["EXTRA_9"] and r["missing_job"] == ["D1K5.1_G_H_R1"], r)
        check("metadata.csv in the job's spelling and order", r["csv"].splitlines()[:2]
              == ["SampleID,Station,Habitat", "D1K1_C_H_R1,D1K1,Coral"] and "D1K5.1_G_H_R1,," in r["csv"], r["csv"])
        try:
            import openpyxl
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.append(["SampleID", "Group"])
            for s in S:
                ws.append([s, s.split("_")[1]])
            b = io.BytesIO()
            wb.save(b)
            rx = MU.check(MU.read_table(b.getvalue(), "m.xlsx"), S)
            check("Excel file read", rx.get("ok") and rx["matched"] == 4 and rx["columns"] == ["Group"], rx)
        except ImportError:
            print("  skip  Excel (openpyxl not installed)")
        check("no sample column -> a clear error", "No column" in MU.check([["a", "b"], ["x", "1"]], S).get("error", ""))
        check("duplicate sample -> error", "more than once" in MU.check(
            [["s", "g"], ["D1K1_C_H_R1", "a"], ["d1k1-c-h-r1", "b"]], S).get("error", ""))
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "asv_table.csv").write_text(",A_1,B_2\n123,1,2\n")
            check("job samples from asv_table.csv", MU.job_samples(td) == ["A_1", "B_2"])
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        check("endpoint: upload keeps the old metadata.csv as a backup",
              '@app.post("/results/{job_id}/metadata/upload")' in m and "metadata.csv.bak-" in m)
        h = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
        check("console: Create group, Manage (search, tick all shown), drag and drop, Upload metadata",
              'id="gNew"' in h and "Tick all shown" in h and '"dragstart"' in h and '"drop"' in h
              and "/metadata/upload" in h and "function groupByName(" in h)
        rb = (ROOT / "backend" / "report_builder.py").read_text(encoding="utf-8")
        check("report takes the grouping chosen in Customise charts (settings.json groupCol)",
              '"edit_charts" / "settings.json"' in rb and 'get("groupCol")' in rb)
    except Exception:
        check("metadata upload", False, traceback.format_exc(limit=5))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
