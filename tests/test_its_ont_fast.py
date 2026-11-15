#!/usr/bin/env python3
"""Nanopore ITS: routed to the VSEARCH OTU pipeline with UNITE (its_pipeline.R is DADA2
paired-end only), MinKNOW/Guppy chunk files (fastq_runid_<id>_<n>) merged into one
sample, UNITE headers parsed, fungal identity floors."""
from __future__ import annotations
import gzip, json, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)

def main():
    print("Nanopore ITS (VSEARCH + UNITE)")
    try:
        import qiime2_vsearch_pipeline as qv
        h = ("Abrothallus_subhalei|MT153946|SH1227328.10FU|reps|k__Fungi;p__Ascomycota;c__Dothideomycetes;"
             "o__Abrothallales;f__Abrothallaceae;g__Abrothallus;s__Abrothallus_subhalei")
        t = qv.parse_lineage(h)
        check("UNITE header parsed", t.get("superkingdom") == "Fungi" and t.get("genus") == "Abrothallus"
              and t.get("species") == "Abrothallus subhalei", str(t))
        t2 = qv.parse_lineage("Fungi_sp|X1|SH9|reps|k__Fungi;p__Ascomycota;c__unidentified;o__unidentified;"
                              "f__unidentified;g__unidentified;s__Fungi_sp")
        check("UNITE 'unidentified' / '_sp' left blank", t2 == {"superkingdom": "Fungi", "phylum": "Ascomycota"}, str(t2))
        check("SILVA trainset still parsed", qv.parse_lineage("Bacteria;Firmicutes;Bacilli;Lactobacillales;"
              "Lactobacillaceae;Lactobacillus;").get("genus") == "Lactobacillus")
        check("ITS identity floors (Vu et al. 2019)", qv._ITS_IDENTITY_FLOOR["species"] == 98.41
              and qv._ITS_IDENTITY_FLOOR["phylum"] == 75.0)
        with tempfile.TemporaryDirectory() as td:
            d = Path(td) / "in"; d.mkdir()
            names = [f"fastq_runid_cae8cbea-7b23_{i}.fastq.gz" for i in range(3)]
            for i, n in enumerate(names):
                with gzip.open(d / n, "wt") as fh:
                    fh.write(f"@r{i}\nACGT\n+\nIIII\n")
            (d / "sample_manifest.json").write_text(json.dumps(
                [{"sample": "PGU-ITS", "file1": n, "file2": ""} for n in names]))
            s = qv.find_fastq(d, Path(td) / "out" / "merged_reads")
            with gzip.open(s["PGU-ITS"], "rt") as fh:
                n_reads = sum(1 for l in fh) // 4
            check("chunk files of one sample merged", list(s) == ["PGU-ITS"] and n_reads == 3, f"{s} {n_reads}")
        import ont_wgs_pipeline as w
        fs = [f"fastq_runid_cae8cbea-7b23-4538-a575-f68dfe0df9cc_{i}.fastq.gz" for i in range(4)] \
            + ["fastq_runid_cae8cbea-7b23-4538-a575-f68dfe0df9cc_4_0.fastq.gz"]
        g = w.group_fastqs(fs, {f: "PGU-ITS" for f in fs})
        check("fastq_runid chunks grouped into the folder's sample", [(x["sample"], len(x["files"])) for x in g]
              == [("PGU-ITS", 5)], str([(x["sample"], len(x["files"])) for x in g]))
        g2 = w.group_fastqs(["L21_1.fastq.gz", "L21_2.fastq.gz"], {})
        check("replicates _1/_2 are not merged", len(g2) == 2)
        cmds = []
        qv.run_cmd = lambda cmd, label="": cmds.append(cmd)
        qv.find_tool = lambda n: n
        with tempfile.TemporaryDirectory() as td:
            qv.trim_primers({"S": Path(td) / "S.fastq.gz"}, Path(td), "TCCGTAGGTGAACCTGCGG",
                            "GCTGCGTTCTTCATCGATGC", 2, 0.2, 10, 100, 900)
        c = cmds[0] if cmds else []
        check("cutadapt (single-end) removes both primers of a read (-n 2) and a 7-base forward-primer stub",
              "-n" in c and c[c.index("-n") + 1] == "2" and "TCCGTAGGTGAACCTGCGG;o=7" in c, c)
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        i = m.find('if marker in ("ITS1", "ITS2", "ITS") and params.sequencerType in ("ont", "qiime2_vsearch")')
        j = m.find('elif marker in ("ITS1", "ITS2", "ITS"):')
        check("main: Nanopore ITS -> VSEARCH pipeline with UNITE before the DADA2 ITS branch",
              0 < i < j and "qiime2_vsearch_pipeline.py" in m[i:j] and "UNITE_ITS1" in m[i:j]
              and '"--marker",         region' in m[i:j])
        h_ = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
        check("form: ITS offers MinKNOW folder merging and the read-length window",
              'it.key === "ITS1" || it.key === "ITS2"' in h_ and "Minimum read length (Nanopore / VSEARCH)" in h_
              and "fastq_runid_" in h_)
    except Exception:
        check("its ont checks", False, traceback.format_exc(limit=4))
    return 1 if FAILS else 0

def make_vsearch_its_run(d: Path):
    """A finished Nanopore ITS run of the VSEARCH OTU pipeline, in its real layout."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "asv_table.csv").write_text(",PGU-ITS\nOTU_1;size=1319,1769\nOTU_2;size=684,932\nOTU_3;size=222,298\n")
    (d / "taxonomy.csv").write_text(
        "ASV,Kingdom,Phylum,Class,Order,Family,Genus,Species\n"
        "OTU_1;size=1319,Fungi,Ascomycota,Eurotiomycetes,Onygenales,Nannizziopsidaceae,Nannizziopsis,Nannizziopsis arthrosporioides\n"
        "OTU_2;size=684,Fungi,Ascomycota,Saccharomycetes,Saccharomycetales,Debaryomycetaceae,Candida,Candida albicans\n"
        "OTU_3;size=222,Fungi,Ascomycota,Eurotiomycetes,Eurotiales,Aspergillaceae,,\n")
    (d / "read_tracking.csv").write_text("sample,input,qc_filtered,primer_trimmed,final\nPGU-ITS,3200,3100,3000,2999\n")
    (d / "otus.fasta").write_text(">OTU_1;size=1319\nACGTACGT\n>OTU_2;size=684\nTTGGCCAA\n>OTU_3;size=222\nGGGTTTAA\n")
    (d / "summary.json").write_text(json.dumps({"pipeline": "qiime2_vsearch", "marker": "ITS1",
        "db_path": "/x/UNITE/sh_general_release_dynamic_19.02.2025.fasta", "otu_similarity": 0.97,
        "unoise_alpha": 2.0, "unoise_minsize": 8, "qc_min_qual": 10, "trim_min_len": 100, "trim_max_len": 900,
        "cutadapt_error_rate": 0.2, "cutadapt_overlap": 10, "tax_identity": 0.75}))
    (d / "run_params.json").write_text(json.dumps({"marker": "ITS1", "sequencerType": "ont",
        "primer_f": "TCCGTAGGTGAACCTGCGG", "primer_r": "GCTGCGTTCTTCATCGATGC",
        "dbPath": "/x/UNITE/sh_general_release_dynamic_19.02.2025.fasta"}))
    for sub in ("qc_filtered", "trimmed", "merged_reads"):
        (d / sub).mkdir()
        with gzip.open(d / sub / "PGU-ITS.fastq.gz", "wt") as fh:
            fh.write("@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
    return d


def package_checks():
    print("Nanopore ITS: Download results package")
    try:
        import amp_data as AD, amp_package as AP, custom_report as CR, packages as PK
        with tempfile.TemporaryDirectory() as td:
            d = make_vsearch_its_run(Path(td) / "run")
            r = AD.load_run(d)
            check("read as feature rows, engine vsearch", r.kind == "emu" and r.engine == "vsearch"
                  and r.samples == ["PGU-ITS"] and r.counts_are_reads, f"{r.kind} {getattr(r, 'engine', '')}")
            t = r.tracking["PGU-ITS"]
            check("tracking read", t["raw_data"] == 3200 and t["qc_filtered"] == 3100 and t["final"] == 2999, str(t))
            check("package kind amplicon", PK.detect_kind(d) == "amplicon")
            man = AP.build_sample(d, None, "PGU-ITS", Path(td) / "pk", formats=("html",), tools={})
            root = Path(man["dir"])
            html = (root / "Report_PGU-ITS.html").read_text()
            check("folders 1-4 made", all((root / f).is_dir() for f in AP.FOLDERS))
            check("OTU centroids of the sample", (root / "2.Taxonomy" / "outfiles" / "otus.fasta").exists())
            fa = (root / "2.Taxonomy" / "outfiles" / "otus.fasta").read_text().splitlines()
            check("otus.fasta named by the tables' SV numbers", fa[0].startswith(">SV001 OTU_1 reads=1769")
                  and "Nannizziopsis arthrosporioides" in fa[0] and fa[1] == "ACGTACGT"
                  and fa[2].startswith(">SV002 OTU_2"), fa[:4])
            check("report: fungal ITS title and VSEARCH method",
                  "Fungal ITS1 amplicon sequencing report" in html and "VSEARCH OTU pipeline" in html
                  and "UNOISE3" in html and "Vu et al. (2019)" in html and "UNITE general release 2025-02-19" in html)
            import re as _re
            # UNOISE2's own title names 16S; figures are embedded as base64, whose random text can hold "Emu"
            body = _re.sub(r"data:[^\"' )]+", "", html.split("<h2>References</h2>")[0])
            check("report: no Emu / 16S wording", "Emu" not in body and "16S" not in body and "minimap2" not in body,
                  [w for w in ("Emu", "16S", "minimap2") if w in body])
            check("report: primers named ITS1 / ITS2", 'ITS1="TCCGTAGGTGAACCTGCGG"' in html)
            refs = CR._refs(r, CR._tax_text(r)[1] + " " + CR._tax_text(r)[0])
            check("references: VSEARCH, UNOISE, UNITE, Vu, chopper; no Emu / DADA2 / MIDORI2",
                  any("VSEARCH" in x for x in refs) and any("UNOISE2" in x for x in refs)
                  and any("UNITE" in x for x in refs) and any("Vu D" in x for x in refs)
                  and any("NanoPack2" in x for x in refs) and not any("Emu:" in x or "DADA2" in x or "MIDORI2" in x
                                                                       for x in refs), refs)
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        check("customer package gate admits Nanopore ITS (VSEARCH)", 'pipe == "qiime2_vsearch"' in m)
    except Exception:
        check("package checks", False, traceback.format_exc(limit=5))


_main_orig = main
def main():
    rc = _main_orig()
    package_checks()
    return 1 if FAILS else rc

if __name__ == "__main__":
    sys.exit(main())
