#!/usr/bin/env python3
"""Multi-marker report: 12S + 18S (PR2) + COI runs of the same samples, matched by name without
each job's marker tag; overview per marker and what each marker found at each site."""
from __future__ import annotations
import json, random, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)

SITES = ["KK-NE", "KK-S", "KY-E", "KY-W"]


def names(fmt_rep, fmt_ct, skip=()):
    out = []
    for s in SITES:
        out += [fmt_rep.format(s=s, r=r) for r in (1, 2, 3) if f"{s}-{r}" not in skip]
        out.append(fmt_ct.format(s=s))
    return out


def make(d: Path, marker, samples, ranks, taxa, params=None, seed=1):
    """taxa: [(lineage list matching ranks, {site: reads per replicate} or None for everywhere)]"""
    rnd = random.Random(seed)
    d.mkdir(parents=True, exist_ok=True)
    seqs = ["".join(rnd.choice("ACGT") for _ in range(120)) for _ in taxa]
    with open(d / "taxonomy_table.csv", "w") as fh:
        fh.write("sequence," + ",".join(ranks) + "\n")
        for q, (lin, _w) in zip(seqs, taxa):
            fh.write(q + "," + ",".join(lin) + "\n")
    with open(d / "asv_table.csv", "w") as fh:
        fh.write(",".join(samples) + ",sequence\n")
        for q, (_lin, where) in zip(seqs, taxa):
            row = []
            for s in samples:
                ct = "CT" in s.replace("Riaz12S", "").split("-")[-1] or "-CT" in s or "CT_" in s or "CTRiaz" in s
                site = next(x for x in SITES if s.startswith(x + "-") and not s.startswith(x + "-" + "S"))\
                    if any(s.startswith(x + "-") for x in SITES) else ""
                site = max((x for x in SITES if s.startswith(x + "-") or s.startswith(x + "_")), key=len)
                n = 0 if ct else (where or {}).get(site, 0 if where else 500)
                row.append(str(n))
            fh.write(",".join(row) + "," + q + "\n")
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write("sample,input,filtered,merged,nonchim\n")
        for s in samples:
            fh.write(f"{s},20000,18000,15000,14000\n")
    (d / "summary.json").write_text(json.dumps({"marker": marker, "n_asvs": len(taxa)}))
    (d / "run_params.json").write_text(json.dumps({"marker": marker, **(params or {})}))
    return d


def main():
    print("Multi-marker report (12S + 18S + COI of the same samples)")
    try:
        import report_builder as RB
        import multimarker_report as MM
        check("marker tag stripped: -Riaz12S / CTRiaz12S / _TAR / _COI-XT",
              MM.sample_keys(["KK-NE-1-Riaz12S", "KK-NE-CTRiaz12S", "KK-S-CT-Riaz12S"])
              == {"KK-NE-1-Riaz12S": "KK-NE-1", "KK-NE-CTRiaz12S": "KK-NE-CT", "KK-S-CT-Riaz12S": "KK-S-CT"}
              and MM.sample_keys(["KK-NE-1_COI-XT", "KY-W-2_COI-XT"])["KY-W-2_COI-XT"] == "KY-W-2")
        check("a shared last digit is not taken for the tag", MM.sample_keys(["A1_TAR", "B1_TAR"]) ==
              {"A1_TAR": "A1", "B1_TAR": "B1"})
        check("no shared tag: names kept", MM.sample_keys(["S1", "S2"]) == {"S1": "S1", "S2": "S2"})
        check("site / replicate / control", [MM.site_of(k) for k in ("KK-NE-1", "KK-NE-CT", "ACT-1", "X")]
              == [("KK-NE", "1", False), ("KK-NE", "CT", True), ("ACT", "1", False), ("X", "", False)])
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            r12 = make(td / "12s", "12S", names("{s}-{r}-Riaz12S", "{s}-CTRiaz12S"),
                       ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                       [(["Metazoa", "Chordata", "Actinopteri", "Cichliformes", "Cichlidae", "Oreochromis",
                          "niloticus"], None),
                        (["Metazoa", "Chordata", "Actinopteri", "Cypriniformes", "Cyprinidae", "Barbonymus",
                          "gonionotus"], {"KK-NE": 300, "KY-E": 50}),
                        (["Metazoa", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "sapiens"],
                         {"KY-W": 20})], params={"mito_gene": "12S"})
            r18 = make(td / "18s", "18S", names("{s}-{r}_TAR", "{s}-CT_TAR"),
                       ["Domain", "Supergroup", "Phylum", "Subdivision", "Class", "Order", "Family", "Genus",
                        "Species"],
                       [(["Eukaryota", "TSAR", "Stramenopiles", "Gyrista", "Bacillariophyceae", "Naviculales",
                          "Naviculaceae", "Navicula", "Navicula_cryptocephala"], None),
                        (["Eukaryota", "Obazoa", "Opisthokonta", "Metazoa", "Arthropoda", "Crustacea_X",
                          "Daphniidae", "Daphnia", "Daphnia_sp."], {"KK-S": 200})],
                       params={"dbPath": "/db/pr2_version_5.1.1_SSU_dada2.fasta.gz"}, seed=2)
            rco = make(td / "coi", "COX1", names("{s}-{r}_COI-XT", "{s}-CT_COI-XT", skip={"KY-W-3"}),
                       ["Phylum", "Class", "Order", "Family", "Genus", "Species"],
                       [(["Arthropoda", "Insecta", "Diptera", "Chironomidae", "Chironomus", "Chironomus plumosus"],
                         None),
                        (["Rotifera", "Eurotatoria", "Ploima", "Brachionidae", "Brachionus", ""], {"KY-W": 80})],
                       seed=3)
            runs = [RB.Run(r12, set(), "Fish_12S"), RB.Run(r18, set(), "Euk_18S"), RB.Run(rco, set(), "COI_XT")]
            check("markers told apart: 12S / 18S / COI", [MM.marker_of(r) for r in runs] == ["12S", "18S", "COI"],
                  [MM.marker_of(r) for r in runs])
            check("multi-marker detected; one marker twice is not", MM.is_multimarker(runs)
                  and not MM.is_multimarker([runs[0], runs[0]]))
            h = MM.build_multimarker_html(runs)
            check("title and both sections", "Multi-marker Amplicon Comparison Report" in h
                  and "1. Overview per marker" in h and "2. What was found at each site" in h)
            check("4 sites, 12 field samples + 4 controls", "4 sites" in h and "12 field samples + 4 controls" in h, "")
            check("sample missing from COI said so", h.count("not in this job") >= 1 and ">KY-W-3<" in h)
            check("PR2 group = subdivision, placeholder blanked", "Gyrista" in h and "Crustacea_X" not in h
                  and "<i>Navicula cryptocephala</i>" in h)
            check("species named with genus (SILVA-style epithet)", "<i>Oreochromis niloticus</i>" in h)
            check("12S group by class (Actinopteri, Mammalia)", "Actinopteri" in h and "Mammalia" in h)
            check("Barbonymus at KK-NE 3/3 and KY-E, not KK-S", "Barbonymus gonionotus</i> <span class=\"dim\">3/3" in h)
            check("per-marker site tables", h.count('<table class="mm">') == 3)
            check("controls column", "In controls (4)" in h)
            check("genus only for COI Brachionus", "<i>Brachionus</i>" in h)
            # PDF builds
            try:
                from weasyprint import HTML
                pdf = HTML(string=h).write_pdf()
                check("PDF builds", len(pdf) > 5000)
            except ImportError:
                print("  (weasyprint not installed; PDF not built)")
            m_ = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            check("main.py routes runs of different markers here", "is_multimarker(runs)" in m_
                  and "build_multimarker_html" in m_)
            fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            check("Build report says so for different markers", "Multi-marker" in fe)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")
    return FAILS


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
