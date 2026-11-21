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
        cmds.clear()
        with tempfile.TemporaryDirectory() as td:
            o = Path(td) / "otus.fasta"
            o.write_text(">OTU_1;size=9\nACGT\n>OTU_2;size=3\nACGA\n")
            def fake(cmd, label=""):
                cmds.append(cmd); Path(cmd[cmd.index("--centroids") + 1]).write_text(">OTU_1;size=12\nACGT\n"); return True
            qv.run_cmd = fake
            f, n0, n1 = qv.merge_otu_variants("vsearch", o, Path(td), 0.97, 2)
            c = cmds[0] if cmds else []
            check("OTU variants merged with cluster_size at the chosen identity, UNOISE set kept",
                  "--cluster_size" in c and c[c.index("--id") + 1] == "0.97" and (n0, n1) == (2, 1)
                  and (Path(td) / "otus_unoise.fasta").exists() and f == o, (c, n0, n1))
            qv.run_cmd = lambda cmd, label="": False
            o.write_text(">OTU_1;size=9\nACGT\n>OTU_2;size=3\nACGA\n")
            f, n0, n1 = qv.merge_otu_variants("vsearch", o, Path(td), 0.97, 2)
            check("failed merge falls back to the UNOISE OTUs", (n0, n1) == (2, 2) and o.read_text().count(">") == 2)
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        check("main: otuMergeId sent to both VSEARCH branches", m.count('"--otu_merge_id",   str(params.otuMergeId or 0)') == 2)
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
            check("report: no merge sentence when OTU variants were not merged", "were then merged" not in html)
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


def ncbi_checks():
    print("Nanopore ITS: NCBI check of unassigned OTUs")
    try:
        import ncbi_unassigned as NU
        lins = {"1": {"superkingdom": "Bacteria", "phylum": "Pseudomonadota", "family": "Lysobacteraceae",
                      "genus": "Stenotrophomonas"},
                "2": {"superkingdom": "Bacteria", "phylum": "Pseudomonadota", "family": "Brucellaceae",
                      "genus": "Brucella"},
                "3": {"superkingdom": "Bacteria", "phylum": "Pseudomonadota", "family": "Brucellaceae",
                      "genus": "Ochrobactrum"},
                "4": {"superkingdom": "Eukaryota", "kingdom": "Viridiplantae", "phylum": "Streptophyta",
                      "family": "Poaceae", "genus": "Dichanthium"}}
        h = lambda tid, pid, cov, bits, sp="x": {"taxid": tid, "pident": pid, "qcov": cov, "bits": bits, "species": sp}
        r = NU.name_hit([h("1", 99.3, 100, 700)], lins)
        check("named to genus at >=95% over >=80%", (r["group"], r["named_as"], r["rank"]) == ("Bacteria", "Stenotrophomonas", "genus"), r)
        r = NU.name_hit([h("2", 99.6, 98, 450), h("3", 99.6, 98, 450)], lins)
        check("tied records of two genera: family", (r["named_as"], r["rank"]) == ("Brucellaceae", "family"), r)
        r = NU.name_hit([h("1", 98.0, 58, 400)], lins)
        check("coverage below 80%: group only", (r["named_as"], r["rank"]) == ("Bacteria", "group"), r)
        r = NU.name_hit([h("1", 85.0, 100, 300)], lins)
        check("85% identity: phylum", r["named_as"] == "Pseudomonadota", r)
        check("plant record -> Plants", NU.name_hit([h("4", 100, 100, 450)], lins)["group"] == "Plants")
        check("no hit -> No NCBI match", NU.name_hit([], lins)["group"] == "No NCBI match")
        with tempfile.TemporaryDirectory() as td:
            d = make_vsearch_its_run(Path(td) / "run")
            (d / "asv_table.csv").write_text((d / "asv_table.csv").read_text() + "OTU_4;size=50,77\n")
            (d / "taxonomy.csv").write_text((d / "taxonomy.csv").read_text() + "OTU_4;size=50,,,,,,,\n")
            (d / "otus.fasta").write_text((d / "otus.fasta").read_text() + ">OTU_4;size=50\nGATTACAGATTACA\n")
            check("unassigned OTUs found (no phylum)", [x[0] for x in NU.unassigned_otus(d)] == ["OTU_4"])
            calls = []
            js = json.dumps({"BlastOutput2": [{"report": {"results": {"search": {
                "query_title": "OTU_4", "query_len": 14, "hits": [{"description": [{"sciname": "Stenotrophomonas maltophilia",
                "taxid": 1, "accession": "CP1", "title": "S. maltophilia chromosome"}],
                "hsps": [{"query_from": 1, "query_to": 14, "identity": 14, "align_len": 14, "bit_score": 28}]}]}}}}]})
            xml = ("<TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>Stenotrophomonas maltophilia</ScientificName>"
                   "<Rank>species</Rank><LineageEx><Taxon><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>"
                   "<Taxon><ScientificName>Pseudomonadota</ScientificName><Rank>phylum</Rank></Taxon>"
                   "<Taxon><ScientificName>Lysobacteraceae</ScientificName><Rank>family</Rank></Taxon>"
                   "<Taxon><ScientificName>Stenotrophomonas</ScientificName><Rank>genus</Rank></Taxon></LineageEx>"
                   "</Taxon></TaxaSet>")
            def http(url, data=None):
                calls.append(url)
                if data:
                    return b"RID = R1\nRTOE = 5\n"
                if "SearchInfo" in url:
                    return b"Status=READY"
                if "JSON2_S" in url:
                    return js.encode()
                return xml.encode()
            info = NU.run(d, http=http, sleep=lambda x: None, log=lambda *a: None)
            rows, inf2 = NU.load(d)
            check("NCBI check written: OTU_4 = Stenotrophomonas (Bacteria)", info["checked"] == 1
                  and rows.get("OTU_4", {}).get("named_as") == "Stenotrophomonas"
                  and rows["OTU_4"]["group"] == "Bacteria" and any("efetch" in c for c in calls), (info, rows))
            import amp_data as AD, amp_package as AP
            r = AD.load_run(d)
            check("run carries the NCBI row on the OTU", any(a.get("ncbi_un") for a in r.asvs))
            man = AP.build_sample(d, None, "PGU-ITS", Path(td) / "pk", formats=("html",), tools={})
            root = Path(man["dir"])
            html = (root / "Report_PGU-ITS.html").read_text()
            check("report: NCBI section with group table and the caveat",
                  "Sequences the reference could not name (NCBI check)" in html and "Stenotrophomonas" in html
                  and "what is present, not how much" in html and "16S rRNA sequencing" in html)
            check("package: 2.Taxonomy/NCBI/unassigned_otus_ncbi.tsv",
                  "Stenotrophomonas" in (root / "2.Taxonomy" / "NCBI" / "unassigned_otus_ncbi.tsv").read_text())
            import report_builder as RB
            rh = RB.build_profile_html(RB.Run(d, set()))
            check("run report: NCBI section, Bacteria share per sample",
                  "Sequences the reference could not name (NCBI check)" in rh and "Stenotrophomonas" in rh
                  and "16S rRNA sequencing" in rh)
            (d / "ncbi_unassigned.json").write_text(json.dumps({"error": "timed out", "checked": 0}))
            (d / "ncbi_unassigned.tsv").unlink()
            man = AP.build_sample(d, None, "PGU-ITS", Path(td) / "pk2", formats=("html",), tools={})
            html = (Path(man["dir"]) / "Report_PGU-ITS.html").read_text()
            check("report: NCBI unreachable is said, not hidden", "NCBI could not be reached" in html)
        rs = {f: (ROOT / "backend" / "r_scripts" / f).read_text(encoding="utf-8")
              for f in ("dada2_pipeline.R", "cox1_pipeline.R", "its_pipeline.R")}
        check("DADA2 / COX1 / ITS cutadapt remove both primers of a read-through pair (-n 2)",
              rs["dada2_pipeline.R"].count('"-n", "2"') == 2 and "-O %d -n 2" in rs["cox1_pipeline.R"]
              and "-O %d -n 2" in rs["its_pipeline.R"])
        th = (ROOT / "backend" / "r_scripts" / "tax_helpers.R").read_text(encoding="utf-8")
        check("unreadable reference: the log says why (web page / empty / no lineage), COX1 names the fix",
              "tax_reference_problem <- function" in th and "web page (HTML)" in th
              and 'iconv(hdrs, "UTF-8", "UTF-8", sub = "?")' in th
              and "tax_reference_problem(dbs$genus)" in rs["cox1_pipeline.R"]
              and "download_databases.sh --only midori" in rs["cox1_pipeline.R"])
        m = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
        q = (ROOT / "backend" / "qiime2_vsearch_pipeline.py").read_text(encoding="utf-8")
        check("main sends --ncbi_unassigned to both VSEARCH branches; pipeline runs it non-fatally",
              m.count('(["--ncbi_unassigned"] if params.ncbiUnassigned else [])') == 2
              and "NU.run(out_dir" in q and "not done:" in q)
        check("NCBI check is its own step (10 of 10), after the plots, off unless ticked",
              'f"Step 10/{NS} — NCBI check' in q and "NS = 10 if args.ncbi_unassigned else 9" in q
              and q.index("run_viz(out_dir") < q.index("Step 10/{NS}") and "ncbiUnassigned: bool = False" in m
              and "Step 9/9" not in q)
    except Exception:
        check("ncbi checks", False, traceback.format_exc(limit=5))


_main_orig = main
def main():
    rc = _main_orig()
    package_checks()
    ncbi_checks()
    return 1 if FAILS else rc

if __name__ == "__main__":
    sys.exit(main())
