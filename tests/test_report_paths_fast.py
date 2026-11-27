#!/usr/bin/env python3
"""Every file or folder a report names must be in what the reader gets: the Download results result
(amplicon: every pipeline, combined and one per sample; WGS: the eight-folder package), the customer
report folder, or the run folder for Build report. Also: the report's feature IDs are those of the
result's count table."""
from __future__ import annotations
import csv, html as H, re, sys, tempfile, traceback, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend")); sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
EXT = r"(?:csv|tsv|xlsx|png|svg|pdf|html|docx|fasta|fa|fna|faa|gbk|gff3?|nwk|txt|json|fastq\.gz|fq\.gz|zip|biom|tsv\.gz)"
TOK = re.compile(r"(?<![\w/.\-])((?:[A-Za-z0-9_.\-]+/)*[A-Za-z0-9_.\-]+\." + EXT + r")\b")
DIR = re.compile(r"(?<![\w/.\-])((?:\d+\.[A-Za-z0-9_\-]+|r_tables|r_plots|PICRUSt2|additional_figures|"
                 r"pipeline_figures|Krona|BLAST|NCBI)/(?:[A-Za-z0-9_.\-]+/)*)")


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def cited(doc: str, refs=()):
    """Paths a report names, outside command lines (tool inputs) and the reference list; reference
    database files (refs) are names of what was used, not files of the result."""
    doc = re.sub(r'<p class="cmd">.*?</p>', " ", doc, flags=re.S)
    doc = re.split(r"<h2[^>]*>\s*(?:\d+\.\s*)?References", doc)[0]
    t = H.unescape(re.sub(r"<[^>]+>", " ", doc))
    return ({x for x in TOK.findall(t) if x not in refs and "<" not in x}, set(DIR.findall(t)))


def missing(files, dirs, names):
    names = [n.replace("\\", "/") for n in names]
    base = {n.rsplit("/", 1)[-1] for n in names}
    alld = {"/".join(n.split("/")[:i]) + "/" for n in names for i in range(1, n.count("/") + 1)}
    out = [f for f in files if not ((f in base) if "/" not in f else any(n == f or n.endswith("/" + f) for n in names))]
    return out + [d for d in dirs if not any(x == d or x.endswith("/" + d) for x in alld)]


def main():
    print("Reports name only files the reader has")
    try:
        import packages as PK
        import report_builder as RB
        import test_report_many_fast as TM
        import test_amp_fast as TA
        import test_its_pacbio_pkg_fast as TI
        import test_its_ont_fast as TV
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            runs = {"ONT 16S (Emu)": TM.make_run(td / "emu", 4), "12S (DADA2 + BLAST)": TA.make_amp_run(td / "d2"),
                    "ITS Illumina": TI.make_its(td / "its"), "PacBio": TI.make_pacbio(td / "pb", False),
                    "PacBio (old table)": TI.make_pacbio(td / "pbo", True),
                    "ITS Nanopore (VSEARCH)": TV.make_vsearch_its_run(td / "vs")}
            S0 = (runs["ONT 16S (Emu)"] / "asv_table.csv").read_text().splitlines()[0].split(",")[1:]
            (runs["ONT 16S (Emu)"] / "metadata.csv").write_text(
                "SampleID,Type\n" + "".join(f"{x},{'AB'[i % 2]}\n" for i, x in enumerate(S0)))
            import json as _json                   # the raw-read QC section names qc/read_qc.json (v2.43.2)
            (runs["ONT 16S (Emu)"] / "qc").mkdir(exist_ok=True)
            (runs["ONT 16S (Emu)"] / "qc" / "read_qc.json").write_text(_json.dumps({
                "reads_per_file_examined": 20000, "platform": "long", "low_depth": [], "tools": {},
                "samples": {x: {"files": [{"file": f"{x}.fastq.gz", "reads_examined": 2000, "mean_len": 1400.0,
                                           "min_len": 120, "max_len": 3000, "median_len": 1500, "mean_q": 17.5,
                                           "pct_q30": 79.0, "long_reads": True, "n50": 1550,
                                           "len_hist_100bp": [0, 10, 20] + [5] * 12 + [900, 1000] + [3] * 23,
                                           "meanq_hist": [0] * 8 + [10, 50, 200, 400, 600, 400, 200, 100, 30, 10]
                                           + [0] * 23}]} for x in S0}}))
            for k, d in runs.items():
                r = RB.Run(d, set())
                refs = {Path(str(r.params.get("dbPath") or r.summary.get("db_path") or "x")).name}
                f, dd = cited(RB.build_profile_html(r), refs)
                miss = missing(f, dd, [str(p.relative_to(d)) for p in d.rglob("*") if p.is_file()])
                check(f"{k}: Build report names only files of the run folder", not miss, miss)
                for combined in (True, False):
                    m = PK.build(d, None, "amplicon", None, not combined, td / f"pk_{len(k)}_{combined}_{d.name}",
                                 job_label="J", formats=("html",), tools={}, combined=combined)
                    bad, ids_bad = [], []
                    for fz in m["files"]:
                        z = zipfile.ZipFile(fz["path"])
                        names = z.namelist()
                        for rep in [n for n in names if re.search(r"(^|/)Report_[^/]*\.html$", n)]:
                            root = rep.rsplit("/", 1)[0] + "/" if "/" in rep else ""
                            sub = [n[len(root):] for n in names if n.startswith(root)]
                            doc = z.read(rep).decode()
                            f, dd = cited(doc, refs)
                            bad += [f"{rep}: {x}" for x in missing(f, dd, sub)]
                            ct = [n for n in sub if re.match(r"4\.[^/]+/[^/]+_count_table_named\.csv$", n)]
                            if ct:
                                tid = {row[0] for row in csv.reader(z.read(root + ct[0]).decode().splitlines()[1:])}
                                shown = set(re.findall(r'<td[^>]*><span class="mono">([^<]+)</span>', doc))
                                ids_bad += [x for x in shown if re.match(r"(SV|ASV|OTU|T)\d+$", x) and x not in tid]
                    mode = "combined" if combined else "one per sample"
                    check(f"{k}: Download results ({mode}) report names only files of the result", not bad, bad[:6])
                    check(f"{k}: Download results ({mode}) report's IDs are the count table's", not ids_bad, ids_bad[:6])
            # WGS: the eight-folder package
            import test_wgs_pkg_fast as TW
            import wgs_customer_report as WR
            for full in (True, False):
                run = TW.make_wgs_run(td / f"wgs{full}", full=full)
                for s in ("FS03", "FS 04"):
                    m = WR.build_sample(run, s, td / f"wpk{full}", ["html"])
                    h = Path(m["html"])
                    t = h.read_text()
                    res = [H.unescape(x).strip() for x in re.findall(r'<p class="res">Results (.*?)</p>', t)]
                    bad = [x for x in res if not (h.parent / x).is_dir()]
                    check(f"WGS package ({'full' if full else 'without trees'}, {s}): every 'Results <folder>' is "
                          "a folder of the package", res and not bad, bad)
            # customer report (lab format): quoted paths and the 'All File' lists
            import custom_report as CR
            fq, fp = TA._fake_qc_tools(td / "qc")
            tools = {"fastqc": fq, "fastp": fp, "Rscript": TA._fake_rscript(td), "ktImportText": None}
            for k, d in (("12S", TA.make_amp_run(td / "c12")), ("ONT 16S", TA.make_emu_run(td / "cemu")),
                         ("ITS Nanopore", TV.make_vsearch_its_run(td / "cvs"))):
                out = td / f"cr_{k.replace(' ', '')}"
                CR.build(d, None, out, tools=tools, formats=("html",), job_label="JOB")
                h = next(out.rglob("*_report.html"))
                t = h.read_text()
                q = {H.unescape(x).replace(" / ", "/").replace(" /", "/").strip()
                     for x in re.findall(r"&ldquo;(.*?)&rdquo;", t)}
                q = {x for x in q if "/" in x or re.search(r"\.\w{2,5}$", x)}
                bad = [x for x in q if not (h.parent / x).exists()]
                for fold, ul in re.findall(r'All File</u></b>\s*([^<]*)</p><ul class="files">(.*?)</ul>', t):
                    bad += [x for x in (H.unescape(y).strip() for y in re.findall(r"<li>-\s*(.*?)</li>", ul))
                            if not (h.parent / H.unescape(fold).strip() / x).exists()]
                check(f"customer report ({k}): every quoted path and listed file is in its folder", q and not bad, bad[:6])
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
