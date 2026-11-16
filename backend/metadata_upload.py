#!/usr/bin/env python3
"""Metadata prepared outside the app (a BI analyst's CSV / TSV / Excel, groups set by the customer)
put into a finished job as its metadata.csv — the file Customise charts, the R re-render and the
run report all read. The sample-name column is found by its values, not its header, and names are
matched to the job's own spelling ignoring case and punctuation (D1K1-C-H-R1 = D1K1_C_H_R1).
Nothing here writes: check() returns what would be saved and what did not match."""
from __future__ import annotations

import csv
import io
import re
from pathlib import Path

MAX_BYTES = 5 * 1024 * 1024


def key(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def job_samples(run_dir) -> list[str]:
    """The job's sample names, from asv_table.csv's header (Emu, VSEARCH and DADA2 layouts)."""
    p = Path(run_dir) / "asv_table.csv"
    if not p.exists():
        return []
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        head = next(csv.reader(fh), [])
    return [h for i, h in enumerate(head) if h.strip() and h.strip().lower() not in ("sequence", "asv", "otu", "#otu id")
            and not (i == 0 and h.strip().lower() in ("", "taxon", "tax_id", "id"))]


def read_table(raw: bytes, filename: str = "") -> list[list[str]]:
    """Rows of a CSV / TSV / TXT / XLSX file as strings, empty rows dropped."""
    fn = (filename or "").lower()
    if fn.endswith((".xlsx", ".xlsm")) or raw[:2] == b"PK":
        import openpyxl
        wb = openpyxl.load_workbook(io.BytesIO(raw), read_only=True, data_only=True)
        ws = wb.worksheets[0]
        rows = [["" if v is None else (str(int(v)) if isinstance(v, float) and v.is_integer() else str(v)).strip()
                 for v in r] for r in ws.iter_rows(values_only=True)]
    else:
        for enc in ("utf-8-sig", "cp874", "latin-1"):
            try:
                txt = raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        first = next((l for l in txt.splitlines() if l.strip()), "")
        delim = max(("\t", ",", ";"), key=lambda d: first.count(d))
        rows = [[c.strip() for c in r] for r in csv.reader(io.StringIO(txt), delimiter=delim)]
    return [r for r in rows if any(c for c in r)]


def check(rows: list[list[str]], samples: list[str]) -> dict:
    """{"ok", "id_col", "columns", "matched", "unmatched_file", "missing_job", "csv"} or {"error"}."""
    if len(rows) < 2:
        return {"error": "The file needs a header row and at least one sample row."}
    head = [h.strip() for h in rows[0]]
    width = max(len(r) for r in rows)
    head += [""] * (width - len(head))
    body = [r + [""] * (width - len(r)) for r in rows[1:]]
    want = {key(s): s for s in samples}
    hits = [sum(1 for r in body if key(r[c]) in want) for c in range(width)]
    idc = max(range(width), key=lambda c: hits[c]) if width else 0
    if not hits or hits[idc] == 0:
        return {"error": "No column of the file holds this job's sample names. First names in the file: "
                + ", ".join(r[0] for r in body[:4]) + "; in the job: " + ", ".join(samples[:4]) + "."}
    seen, dup = {}, []
    for r in body:
        k = key(r[idc])
        if k in want:
            if k in seen:
                dup.append(r[idc])
            seen[k] = r
    if dup:
        return {"error": "These samples appear more than once in the file: " + ", ".join(sorted(set(dup))) + "."}
    rows_in = [seen[k] for k in seen]

    def vals(c):
        return [r[c] for r in rows_in if r[c]]
    # a running number (1, 2, 3 …, one per row) is not a grouping
    cols = [c for c in range(width) if c != idc and head[c] and vals(c)
            and not (len(set(vals(c))) == len(vals(c)) == len(rows_in) and all(re.fullmatch(r"\d+", v) for v in vals(c)))]
    if not cols:
        return {"error": "The file has the sample names but no other column (a group) to add."}
    names = []
    for c in cols:
        n = re.sub(r'[,\n\r"]', " ", head[c]).strip()
        while n in names or n.lower() in ("sampleid", "sample"):
            n += "_2"
        names.append(n)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["SampleID"] + names)
    for s in samples:                                   # the job's order and spelling
        r = seen.get(key(s))
        w.writerow([s] + ([r[c] for c in cols] if r else [""] * len(cols)))
    # the grouping to show first: the first column that groups (2 values or more, fewer than the samples)
    pick = next((n for c, n in zip(cols, names) if 2 <= len(set(vals(c))) < len(rows_in)), names[0])
    return {"ok": True, "id_col": head[idc] or f"column {idc + 1}", "columns": names, "pick": pick,
            "matched": len(seen), "samples": len(samples),
            "unmatched_file": [r[idc] for r in body if r[idc] and key(r[idc]) not in want],
            "missing_job": [s for s in samples if key(s) not in seen], "csv": buf.getvalue()}
