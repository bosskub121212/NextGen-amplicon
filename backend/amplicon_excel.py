"""amplicon_excel.py — Excel workbooks of the customer report: one per sample and one
for all samples (openpyxl)."""
from __future__ import annotations

from pathlib import Path

from amp_data import RANKS, TRACK_COLS, RunData

BLAST_COLS = ["best_hit", "pident", "qcov", "source", "accession", "ncbi_best", "ncbi_pident", "ncbi_accession",
              "agreement", "check", "target"]


def _sheet(wb, title, head, rows):
    from openpyxl.styles import Font
    ws = wb.create_sheet(title)
    ws.append(head)
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in rows:
        ws.append(r)
    ws.freeze_panes = "A2"
    for col in ws.columns:
        w = max(len(str(c.value or "")) for c in col[:200])
        ws.column_dimensions[col[0].column_letter].width = min(60, max(8, w + 2))
    return ws


def _book(run: RunData, samples: list[str]):
    import openpyxl
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    _sheet(wb, "Read_tracking", ["SampleID", "Sample_name"] + TRACK_COLS,
           [[s, s] + [run.tracking[s][k] for k in TRACK_COLS] for s in samples])
    asvs = [a for a in run.asvs if any(a["counts"].get(s, 0) for s in samples)]
    head = ["SV", "sequence"] + RANKS + BLAST_COLS + [f"{s} reads" for s in samples] + [f"{s} %" for s in samples]
    tot = {s: sum(a["counts"].get(s, 0) for a in run.asvs) or 1 for s in samples}
    _sheet(wb, "ASV_taxonomy", head,
           [[a["sv"], a["sequence"]] + [a["tax"].get(r, "") for r in RANKS] + [a.get(k, "") for k in BLAST_COLS]
            + [a["counts"].get(s, 0) for s in samples]
            + [round(100 * a["counts"].get(s, 0) / tot[s], 4) for s in samples] for a in asvs])
    for r in RANKS:
        tabs = {s: run.rank_table(r, s) for s in samples}
        labs = sorted({k for t in tabs.values() for k in t}, key=lambda k: -sum(t.get(k, 0) for t in tabs.values()))
        _sheet(wb, r, [r] + samples + [f"{s} %" for s in samples],
               [[k] + [tabs[s].get(k, 0) for s in samples]
                + [round(100 * tabs[s].get(k, 0) / tot[s], 4) for s in samples] for k in labs])
    _sheet(wb, "Alpha_diversity", ["Sample_ID", "Sample_name", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson"],
           [[s, s] + [run.alpha[s][k] for k in ("Observed", "Chao1", "se.chao1", "Shannon", "Simpson")]
            for s in samples])
    return wb


def write_sample(path: Path, run: RunData, sample: str) -> Path:
    _book(run, [sample]).save(path)
    return Path(path)


def write_all(path: Path, run: RunData) -> Path:
    _book(run, list(run.samples)).save(path)
    return Path(path)
