#!/usr/bin/env python3
"""Fast checks for the ONT WGS customer package (folders 1–8 + lab-style report),
read-based taxonomy, TYGS import and per-sample downloads.

Run alone (python3 tests/test_wgs_pkg_fast.py) or from tests/test_fast.py.
"""
from __future__ import annotations

import json
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))

FAILS = []
TESTS = []
FOLDERS = ["1.QcReads", "2.Assembly", "3.Completeness", "4.Taxonomy", "5.ANI", "6.Annotation", "7.AMR",
           "8.readsMappingToRef"]


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILS.append(name)


def _w(p: Path, text: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def make_wgs_run(d: Path, samples=("FS03", "FS 04"), full=True):
    """A finished ONT-WGS results folder in ont_wgs_pipeline.py's layout."""
    d.mkdir(parents=True, exist_ok=True)
    _w(d / "run_params.json", json.dumps({"marker": "ONT-WGS", "min_read_len": 1000, "min_read_q": 10,
                                          "job_name": "TH-BIOENTIST-00002"}))
    _w(d / "pipeline.log", "internal log\n")
    for k, S in enumerate(samples):
        s = d / "samples" / S
        seqA = ("ACGT" * 1500)            # 6000 bp, GC 50 %
        seqB = ("AATT" * 750)             # 3000 bp, GC 0 %
        seqC = ("GGCC" * 250)             # 1000 bp, GC 100 %
        _w(s / f"{S}_assembly.fasta", f">{S}_c1\n{seqA}\n>{S}_c2\n{seqB}\n>{S}_c3\n{seqC}\n")
        _w(s / f"{S}_assembly_graph.gfa", "H\tVN:Z:1.0\n")
        _w(s / "flye_assembly_info.txt", "#seq_name\tlength\ncontig_1\t6000\n")
        _w(s / f"{S}_16S.fasta", f">{S}_16S_1\nACGT\n")
        _w(s / "barrnap.gff", "##gff-version 3\n")
        _w(s / "sourmash_gather.csv", "match_name,f_unique_weighted\nGCF_1 Bacillus subtilis,0.99\n")
        _w(s / "sourmash_prefetch.csv", "internal\n")
        _w(s / "assembly.sig.zip", "internal")
        _w(s / "skani.tsv", "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\nref\tq\t99.1\t95\n")
        _w(s / "mlst.tsv", f"{S}\tbsubtilis\t{12 + k}\n")
        _w(s / "amrfinder.tsv", "Element symbol\tClass\nvmlR\tLINCOSAMIDE\n")
        _w(s / "abricate_vfdb.tsv", "#FILE\tGENE\n")
        _w(s / "abricate_plasmidfinder.tsv", "#FILE\tGENE\n")
        _w(s / "mob_recon" / "contig_report.txt", "sample_id\tmolecule_type\n")
        _w(s / "mob_recon" / "mge.report.txt", "x\n")
        _w(s / "checkm2" / "quality_report.tsv", "Name\tCompleteness\tContamination\nasm\t100\t0.08\n")
        _w(s / "checkm2" / "diamond_output" / "DIAMOND_RESULTS.tsv", "big internal\n")
        _w(s / "checkm2" / "checkm2.log", "log\n")
        for ext in ("gff3", "gbff", "faa", "ffn", "tsv", "json", "txt", "png", "svg", "embl", "fna"):
            _w(s / "bakta" / f"{S}.{ext}", f"bakta {ext}\n")
        _w(s / "bakta" / f"{S}.log", "bakta log\n")
        if full:
            _w(s / "taxonomy" / "type_strain_ani.csv", "organism,strain,accession,anib\nBacillus subtilis,168,GCF_1,99.3\n")
            _w(s / "taxonomy" / "type_strain_16S.csv", "organism,identity\nBacillus subtilis,99.9\n")
            _w(s / "taxonomy" / f"{S}_genome_tree.nwk", f"({S}:0.1,Bsub:0.1);\n")
            _w(s / "taxonomy" / f"{S}_16S_tree.nwk", f"({S}:0.1,Bsub:0.1);\n")
            _w(s / "taxonomy" / "ani_matrix_skani.tsv", "x\n")
            _w(s / "mapping" / "reference_coverage_bins.tsv", "start\tdepth\n0\t200\n")
            _w(s / "mapping" / "reference_absent_regions.csv", "start,end\n")
            _w(s / "mapping" / "reads_vs_reference.bam", "BAM")
            _w(s / "mapping" / "reads_vs_reference.bam.bai", "BAI")
            _w(s / "mapping" / "assembly_vs_reference.paf", "paf\n")
        for t in ("flye", "medaka", "checkm2", "bakta", "seqkit"):
            _w(s / f"log_{t}.log", "log\n")
        _w(s / ".done_assembly", "sig")
        res = {"sample": S, "status": "ok",
               "reads_raw": {"reads": 254102, "bases": 1502546932, "mean_len": 5913.16, "median_len": 5000,
                             "n50": 9000, "mean_q": 18.2, "len_bins": [0, 500, 1000, 2000, 5000],
                             "len_hist": [10, 20, 30, 40, 50], "len_hist_bp": [100, 2000, 30000, 400000, 5000000]},
               "reads_filtered": {"reads": 221340, "bases": 1479063858, "mean_len": 6682.32, "median_len": 6000,
                                  "n50": 9500, "mean_q": 19.0},
               "assembly": {"contigs": 3, "total_length": 10000, "n50": 6000, "largest": 6000, "gc": 50.0,
                            "circular": 1, "mean_depth": 202.0, "polish": "medaka (bacterial methylation model)"},
               "contigs": [{"contig": f"{S}_c1", "length": 6000, "gc": 50.0, "depth": 200.0, "circular": True,
                            "type": "chromosome"}],
               "checkm2": {"completeness": 100.0, "contamination": 0.08, "model": "Neural Network (Specific Model)"},
               "identification": {"method": "ANIb vs type strain (NCBI type material)", "species": "Bacillus subtilis",
                                  "genus": "Bacillus", "ani": 99.3, "reference": "Bacillus subtilis 168 (T)",
                                  "accession": "GCF_1", "call": "species"},
               "mlst": {"scheme": "bsubtilis", "st": str(12 + k)},
               "bakta": {"summary": {"CDSs": "4390", "rRNAs": "30", "tRNAs": "86", "tmRNAs": "1"}}}
        res["sourmash"] = [{"rank": 1, "accession": "GCF_1", "gtdb_species": "s__Bacillus subtilis",
                            "gtdb_lineage": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;"
                                            "g__Bacillus;s__Bacillus subtilis", "f_query": 98.5, "ani": 99.1}]
        if full:
            res["typestrain"] = {"verdict": {"call": "known species", "level": "known", "species": "Bacillus subtilis"},
                                 "type_strains": [{"organism": "Bacillus subtilis", "strain": "168", "accession": "GCF_1",
                                                   "anib": 99.3, "anib_cov_q": 92.1, "anim": 99.5, "tetra": 0.9991,
                                                   "skani_ani": 99.2}],
                                 "genome_tree": f"({S}:0.01,(Bsub_168:0.02,Bvel:0.05):0.01);",
                                 "genome_tree_method": "skani ANI distances + FastME",
                                 "rrna_tree": f"({S}:0.001,(Bsub_168:0.002,Bvel:0.004):0.001);",
                                 "rrna_tree_method": "MAFFT + IQ-TREE"}
            res["mapping"] = {"reference": "Bacillus subtilis", "accession": "NZ_CP010052.1", "ref_length": 4215619,
                              "reads_mapped_pct": 62.99, "mean_depth": 202.01, "breadth_1x": 90.88,
                              "breadth_10x": 89.5, "snps": 1234, "indels": 56, "bin_size": 1000,
                              "coverage_bins": [["NZ_CP010052.1", i * 1000, 150 + (i % 7) * 20] for i in range(50)]}
        _w(s / "result.json", json.dumps(res))
    return d


def _tree(p: Path) -> list[str]:
    return sorted(str(x.relative_to(p)) for x in p.rglob("*") if x.is_file())


# ── Task 1: package layout ─────────────────────────────────────────────────────
def t_assembly_stats(tmp: Path):
    import wgs_package as WP
    run = make_wgs_run(tmp / "run")
    st = dict(WP.assembly_stats(run / "samples" / "FS03" / "FS03_assembly.fasta"))
    check("assembly stats: contigs / largest / total", st.get("# contigs") == 3 and st.get("Largest contig") == 6000
          and st.get("Total length") == 10000, str(st))
    check("assembly stats: length thresholds", st.get("# contigs (>= 5000 bp)") == 1
          and st.get("Total length (>= 1000 bp)") == 10000 and st.get("# contigs (>= 25000 bp)") == 0, str(st))
    check("assembly stats: N50 / N90 / L50 / L90 / auN", st.get("N50") == 6000 and st.get("N90") == 3000
          and st.get("L50") == 1 and st.get("L90") == 2 and st.get("auN") == 4600, str(st))
    check("assembly stats: GC %", abs(st.get("GC (%)", 0) - 40.0) < 0.01, str(st))


def t_package_layout(tmp: Path):
    import wgs_package as WP
    run = make_wgs_run(tmp / "run")
    m = WP.build_package(run, "FS03", tmp / "out")
    root = Path(m["dir"])
    check("package: root named after the sample", root.name == "FS03", str(root))
    check("package: the eight folders, exactly named", sorted(p.name for p in root.iterdir() if p.is_dir()) == FOLDERS,
          str(sorted(p.name for p in root.iterdir())))
    files = _tree(root)
    want = ["1.QcReads/read_statistics.tsv", "2.Assembly/FS03_assembly.fasta", "2.Assembly/assembly_statistics.tsv",
            "2.Assembly/FS03_assembly_graph.gfa", "3.Completeness/checkm2_quality_report.tsv",
            "4.Taxonomy/Contig-based Taxonomy/sourmash_gather_gtdb.csv", "4.Taxonomy/16S/FS03_16S.fasta",
            "4.Taxonomy/MLST/mlst.tsv", "4.Taxonomy/Phylogeny/FS03_genome_tree.nwk",
            "4.Taxonomy/identification.tsv", "5.ANI/type_strain_ani.csv", "5.ANI/skani_closest_genome.tsv",
            "6.Annotation/FS03.gbff", "6.Annotation/FS03.gff3", "7.AMR/amrfinderplus.tsv",
            "7.AMR/virulence_vfdb.tsv", "7.AMR/plasmidfinder.tsv", "7.AMR/MOB-suite/contig_report.txt",
            "8.readsMappingToRef/reads_vs_reference.bam", "8.readsMappingToRef/reference_coverage_bins.tsv"]
    miss = [w for w in want if w not in files]
    check("package: every step's files in its folder", not miss, str(miss))
    bad = [f for f in files if f.endswith(".log") or "/.done" in f or "diamond_output" in f or f.endswith("result.json")
           or "sig.zip" in f or "prefetch" in f]
    check("package: no tool logs or internal files", not bad, str(bad))
    check("package: nothing of the other sample", not [f for f in files if "FS 04" in f or "FS_04" in f])
    rs = (root / "1.QcReads" / "read_statistics.tsv").read_text()
    check("package: read statistics before / after filtering", "254102" in rs and "221340" in rs and "87.11" in rs, rs)
    check("package: manifest lists the files per folder", "2.Assembly/FS03_assembly.fasta" in m["files"]
          and set(m["folders"]) == set(FOLDERS), str(m.get("folders")))


def t_package_missing_steps(tmp: Path):
    import wgs_package as WP
    run = make_wgs_run(tmp / "run", full=False)
    m = WP.build_package(run, "FS 04", tmp / "out")
    root = Path(m["dir"])
    check("package: odd sample name made a safe folder name", root.name == "FS_04", root.name)
    note = root / "8.readsMappingToRef" / "NOT_AVAILABLE.txt"
    check("package: a step that did not run leaves a note in its folder", note.exists()
          and "not" in note.read_text().lower(), str(_tree(root)))
    check("package: missing steps reported in the manifest", "8.readsMappingToRef" in m["missing"], str(m["missing"]))


TESTS += [t_assembly_stats, t_package_layout, t_package_missing_steps]



# ── Task 2: lab-style report ──────────────────────────────────────────────────
WGS_HEADINGS = ["Whole genome sequencing report", "Sample: FS03", "1.QC", "Before Trimming", "After Trimming",
                "Results 1.QcReads", "2.Genome Assembly", "Results 2.Assembly",
                "3. Estimates of genome completeness and contamination", "Results 3.Completeness",
                "4. Genome-based taxonomy", "Table 4.1 Genome-based identification of strain", "Phylogenetic tree",
                "Results 4.Taxonomy/Phylogeny", "Read-based Taxonomy Classification",
                "Results 4.Taxonomy/Read-based Taxonomy", "Contig-based Community Profiling",
                "Results 4.Taxonomy/Contig-based Taxonomy", "MLST", "Results 4.Taxonomy/MLST", "5. ANI analysis",
                "Results 5.ANI", "6. Genome annotation", "Results 6.Annotation", "7.Antimicrobial Resistance Genes",
                "Results 7.AMR", "8. Reads Mapped to Reference", "Results 8.readsMappingToRef", "Reference"]


def _in_order(text, items):
    pos = 0
    for it in items:
        i = text.find(it, pos)
        if i < 0:
            return False, it
        pos = i + len(it)
    return True, ""


def _docx_ok():
    try:
        import docx  # noqa: F401
        return True
    except ImportError:
        return False


def t_report_sections(tmp: Path):
    import wgs_customer_report as WR
    run = make_wgs_run(tmp / "run")
    m = WR.build_sample(run, "FS03", tmp / "out", formats=("html", "docx") if _docx_ok() else ("html",))
    root = Path(m["dir"])
    html = (root / "Report_FS03.html").read_text()
    ok, miss = _in_order(html, WGS_HEADINGS)
    check("WGS report: lab headings in order, each section ends with its Results line", ok, f"missing: {miss}")
    check("WGS report: no cover page", "Client Name" not in html and "Order Number" not in html)
    check("WGS report: page header 'Page |'", "Page |" in html)
    check("WGS report: QC numbers as the lab writes them", "Reads #: 254102" in html
          and "Reads #: 221340 (87.11 %)" in html and "Total bases: 1479063858 (98.44 %)" in html
          and "Mean Reads Length: 6682.32" in html, "")
    check("WGS report: our assembler, real numbers", "FS03 was assembled using Flye" in html and "3 contigs" in html
          and "10,000 bp" in html and "40.00%" in html and "unicycler" not in html.lower(), "")
    check("WGS report: our tools named, the lab's other tools not claimed",
          all(k in html for k in ("Bakta", "AMRFinderPlus", "CheckM2", "ANIb")) and "PROKKA" not in html
          and "CARD" not in html and "dDDH value" not in html, "")
    check("WGS report: identification table from our ANI verdict", "belongs to known species" in html
          and "<i>Bacillus subtilis</i>" in html and "99.3" in html, "")
    check("WGS report: read-based taxonomy said not run when absent", "Read-based taxonomy was not run" in html)
    check("WGS report: dDDH said not imported when no TYGS result", "TYGS" in html and "not imported" in html)
    check("WGS report: mapping table uses only computed metrics", "62.99" in html and "1,234" in html)
    check("WGS report: written into the package root", (root / "Report_FS03.html").exists()
          and (not _docx_ok() or (root / "Report_FS03.docx").exists()), str(sorted(p.name for p in root.iterdir())))
    try:
        import cairosvg  # noqa: F401
        check("WGS report: figures also saved as PNG in their folders",
              (root / "1.QcReads" / "read_length_histogram.png").exists()
              and (root / "8.readsMappingToRef" / "coverage_plot.png").exists()
              and (root / "4.Taxonomy" / "Phylogeny" / "genome_tree.png").exists(), "")
    except ImportError:
        print("  skip  PNG figures (cairosvg not installed here)")


def t_report_missing_steps(tmp: Path):
    import wgs_customer_report as WR
    run = make_wgs_run(tmp / "run", full=False)
    m = WR.build_sample(run, "FS 04", tmp / "out", formats=("html",))
    html = (Path(m["dir"]) / "Report_FS_04.html").read_text()
    ok, miss = _in_order(html, ["8. Reads Mapped to Reference", "not performed", "Results 8.readsMappingToRef"])
    check("WGS report: a step that did not run says so, Results line kept", ok, miss)
    check("WGS report: sample name shown as given", "Sample: FS 04" in html)
    check("WGS report: no other sample in it", "FS03" not in html)


def t_report_skani_only(d):
    import wgs_customer_report as R
    res = d / "res"
    make_wgs_run(res, samples=("FS03",))
    rj = res / "samples" / "FS03" / "result.json"
    r = json.loads(rj.read_text())
    for t in r["typestrain"]["type_strains"]:
        for k in ("anib", "anib_cov_q", "anim", "tetra"):
            t.pop(k, None)
    r["identification"]["method"] = "skani ANI vs closest GTDB reference"
    rj.write_text(json.dumps(r))
    html = Path(R.build_sample(res, "FS03", d / "pk", ["html"])["html"]).read_text()
    check("WGS report: ANIb run off → no ANIb / ANIm / TETRA / JSpecies named",
          "ANIb" not in html and "ANIm" not in html and "TETRA" not in html and "JSpecies" not in html
          and "skani" in html, "")


TESTS += [t_report_sections, t_report_missing_steps, t_report_skani_only]


# ── read-based taxonomy (readtax.py) ─────────────────────────────────────────
BRACKEN = ("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\t"
           "fraction_total_reads\n"
           "Bacillus subtilis\t1423\tS\t900\t60\t960\t0.96000\n"
           "Bacillus velezensis\t492670\tS\t30\t10\t40\t0.04000\n")
K2REP = ("  2.00\t20\t20\tU\t0\tunclassified\n"
         " 98.00\t980\t0\tR\t1\troot\n"
         " 90.00\t900\t900\tS\t1423\t                    Bacillus subtilis\n"
         "  3.00\t30\t30\tS\t492670\t                    Bacillus velezensis\n")
CFREP = ("name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n"
         "Bacillus_subtilis_subsp._subtilis_str._168\t224308\tstrain\t4215606\t700\t650\t0.7\n"
         "Bacillus subtilis\t1423\tspecies\t4215606\t800\t700\t0.0\n"
         "Bacillus\t1386\tgenus\t0\t1000\t900\t0.0\n"
         "Escherichia_coli\t562\tspecies\t5000000\t200\t150\t0.0\n")
G2SUM = ("LEVEL\tNAME\tTAXID\tREAD_COUNT\tTOTAL_BP_MAPPED\tTOTAL_BP_MISMATCH\tLINEAR_LENGTH\tLINEAR_DOC\t"
         "ROLLUP_DOC\tREL_ABUNDANCE\n"
         "genus\tBacillus\t1386\t500\t1\t1\t1\t1\t1\t1.0\n"
         "species\tBacillus subtilis\t1423\t480\t1\t1\t1\t1\t1\t0.9\n"
         "species\tBacillus spizizenii\t96241\t20\t1\t1\t1\t1\t1\t0.1\n"
         "strain\tBacillus subtilis 168\t224308\t400\t1\t1\t1\t1\t1\t1.0\n")


def t_readtax_parsers(d):
    import readtax as RT
    b = RT.parse_bracken(_w(d / "b.tsv", BRACKEN))
    check("readtax: Bracken species share in %", abs(b["Bacillus subtilis"] - 96.0) < 1e-6
          and abs(b["Bacillus velezensis"] - 4.0) < 1e-6, str(b))
    k = RT.parse_kraken_species(_w(d / "k.report", K2REP))
    check("readtax: Kraken report fallback = share of species-level reads",
          abs(k["Bacillus subtilis"] - 96.774) < 0.01 and "root" not in k and "unclassified" not in k, str(k))
    check("readtax: Kraken classified % from the report", abs(RT.kraken_classified_pct(d / "k.report") - 98.0) < 1e-6)
    c = RT.parse_centrifuge(_w(d / "c.tsv", CFREP))
    check("readtax: Centrifuge(r) species rows only, reads share, underscores → spaces",
          abs(c["Bacillus subtilis"] - 80.0) < 1e-6 and abs(c["Escherichia coli"] - 20.0) < 1e-6
          and "Bacillus" not in c and not any("168" in n for n in c), str(c))
    g = RT.parse_gottcha2(_w(d / "g.tsv", G2SUM))
    check("readtax: GOTTCHA2 species rows, REL_ABUNDANCE in %",
          abs(g["Bacillus subtilis"] - 90.0) < 1e-6 and abs(g["Bacillus spizizenii"] - 10.0) < 1e-6
          and len(g) == 2, str(g))
    check("readtax: empty / missing outputs give {}", RT.parse_bracken(d / "nope") == {}
          and RT.parse_gottcha2(_w(d / "e.tsv", "")) == {})


def t_readtax_table_heatmap(d):
    import readtax as RT
    prof = {"kraken2": {"Bacillus subtilis": 96.0, "Bacillus velezensis": 4.0},
            "gottcha2": {"Bacillus subtilis": 90.0, "Bacillus spizizenii": 10.0}}
    tsv = RT.write_table(prof, d / "read_taxonomy.tsv")
    rows = [l.split("\t") for l in tsv.read_text().splitlines()]
    check("readtax: table header = Taxon + tools in fixed order",
          rows[0] == ["Taxon", "Kraken2+Bracken (%)", "GOTTCHA2 (%)"], str(rows[0]))
    check("readtax: taxa sorted by highest value, absent = 0",
          rows[1][0] == "Bacillus subtilis" and ["Bacillus spizizenii", "0", "10"] in
          [[r[0], r[1], r[2]] for r in rows], str(rows))
    sv = RT.svg_heatmap(tsv)
    check("readtax: heatmap SVG names tools and taxa", sv.startswith("<svg") and "Kraken2+Bracken" in sv
          and "GOTTCHA2" in sv and "Bacillus subtilis" in sv and "Centrifuge" not in sv)
    one = RT.write_table({"kraken2": {"A b": 100.0}}, d / "one.tsv")
    check("readtax: heatmap with a single tool column", "<rect" in RT.svg_heatmap(one))
    many = RT.write_table({"kraken2": {f"Taxon sp{i}": 100 / 40 for i in range(40)}}, d / "many.tsv")
    check("readtax: heatmap keeps the top 20 taxa", RT.svg_heatmap(many, top=20).count("Taxon sp") == 20)


def _fake_tools(d: Path) -> Path:
    """Executables that write the fixtures where the real tools would."""
    b = d / "bin"
    b.mkdir()
    scripts = {
        "kraken2": f"""import sys; a=sys.argv; open(a[a.index('--report')+1],'w').write({K2REP!r})""",
        "bracken": f"""import sys; a=sys.argv; open(a[a.index('-o')+1],'w').write({BRACKEN!r})
open(a[a.index('-w')+1],'w').write('x')""",
        "centrifuger": """import sys; print('readID\tseqID\ttaxID'); print('r1\tNC_1\t1423')""",
        "centrifuger-quant": f"""import sys; sys.stdout.write({CFREP!r})""",
        "gottcha2": f"""import sys,os; a=sys.argv; o=a[a.index('-o')+1]; p=a[a.index('-p')+1]
os.makedirs(o,exist_ok=True); open(os.path.join(o,p+'.tsv'),'w').write({G2SUM!r})""",
    }
    for n, body in scripts.items():
        f = b / n
        f.write_text(f"#!{sys.executable}\n{body}\n")
        f.chmod(0o755)
    return b


def _finder(bindir: Path, have=None):
    import ont_wgs_pipeline as W

    def find(name):
        if have is not None and name not in have:
            return None
        p = bindir / name
        return W.Tool(name, p, bindir) if p.exists() else None
    return find


def t_readtax_run_sample(d):
    import readtax as RT
    bindir = _fake_tools(d)
    reads = _w(d / "reads.filtered.fastq.gz", "@r1\nACGT\n+\nIIII\n")
    k2db = d / "k2_pluspf_16_GB_20260626"
    _w(k2db / "hash.k2d", "x")
    _w(k2db / "database300mers.kmer_distrib", "x")
    _w(k2db / "database150mers.kmer_distrib", "x")
    cf = d / "cfr" / "cfr_hpv+gbsarscov2"
    _w(Path(str(cf) + ".1.cfr"), "x")
    dbs = {"kraken2_db": str(k2db), "readtax_centrifuge": str(cf), "readtax_gottcha2": ""}
    out = d / "readtax"
    info = RT.run_sample(reads, out, dbs, threads=2, mean_len=6000, find=_finder(bindir))
    st = {t["key"]: t["status"] for t in info["tools"]}
    check("readtax run: tools with DB ran, GOTTCHA2 without DB skipped",
          st == {"kraken2": "ok", "centrifuge": "ok", "gottcha2": "skipped"}, str(info["tools"]))
    g = [t for t in info["tools"] if t["key"] == "gottcha2"][0]
    check("readtax run: skipped tool says why", "database" in g["note"].lower(), g["note"])
    k = [t for t in info["tools"] if t["key"] == "kraken2"][0]
    check("readtax run: Bracken uses the longest k-mer distribution for long reads",
          k.get("bracken_read_len") == 300 and k["db"] == "k2_pluspf_16_GB_20260626", str(k))
    c = [t for t in info["tools"] if t["key"] == "centrifuge"][0]
    check("readtax run: Centrifuger index recognised", c["label"] == "Centrifuger"
          and c["db"] == "cfr_hpv+gbsarscov2", str(c))
    tsv = (out / "read_taxonomy.tsv").read_text()
    check("readtax run: combined table has the two tools that ran",
          tsv.splitlines()[0] == "Taxon\tKraken2+Bracken (%)\tCentrifuger (%)", tsv.splitlines()[0])
    check("readtax run: per-tool outputs kept, heatmap PNG/SVG written",
          (out / "kraken2.report").exists() and (out / "bracken_species.tsv").exists()
          and (out / "centrifuger_report.tsv").exists() and (out / "heatmap_read_taxonomy.svg").exists())
    check("readtax run: classification file (per read, large) not kept",
          not (out / "centrifuger_classification.tsv").exists())
    loaded = RT.load_info(out)
    check("readtax run: info saved with references", set(loaded["tools_ref"]) == {"kraken2", "bracken", "centrifuger"},
          str(loaded.get("tools_ref")))
    ms = RT.method_sentence(loaded)
    check("readtax: method sentence names tools, DBs and the skipped one",
          "Kraken2" in ms and "k2_pluspf_16_GB_20260626" in ms and "Centrifuger" in ms
          and "GOTTCHA2" in ms and "not run" in ms, ms)
    # only Kraken2 installed
    out2 = d / "readtax2"
    info2 = RT.run_sample(reads, out2, dbs, threads=2, mean_len=6000,
                          find=_finder(bindir, have={"kraken2", "bracken"}))
    st2 = {t["key"]: t["status"] for t in info2["tools"]}
    check("readtax run: only Kraken2 installed → one column, others noted",
          st2["kraken2"] == "ok" and st2["centrifuge"] == "skipped"
          and (out2 / "read_taxonomy.tsv").read_text().splitlines()[0] == "Taxon\tKraken2+Bracken (%)",
          str(info2["tools"]))
    check("readtax run: a rerun drops the output of a tool that did not run this time",
          not (out2 / "centrifuger_report.tsv").exists())
    (out2 / "gottcha2.tsv").write_text("old")
    RT.run_sample(reads, out2, dbs, threads=2, mean_len=6000, find=_finder(bindir, have={"kraken2", "bracken"}))
    check("readtax run: stale files of skipped tools removed", not (out2 / "gottcha2.tsv").exists())
    # nothing at all
    info3 = RT.run_sample(reads, d / "rt3", {}, threads=2, mean_len=6000, find=lambda n: None)
    check("readtax run: no tools → no table, no crash", not (d / "rt3" / "read_taxonomy.tsv").exists()
          and all(t["status"] == "skipped" for t in info3["tools"]))


def t_readtax_cli_finished_job(d):
    import readtax as RT
    bindir = _fake_tools(d)
    res = d / "results"
    make_wgs_run(res, samples=("FS03", "FS 04"))
    work = d / "uploads" / "_wgs_work"
    _w(work / "FS03" / "reads.filtered.fastq.gz", "@r\nA\n+\nI\n")      # FS 04 has no reads left
    k2db = d / "k2"
    _w(k2db / "hash.k2d", "x")
    _w(k2db / "database150mers.kmer_distrib", "x")
    dbp = _w(d / "db_paths.json", json.dumps({"kraken2_db": str(k2db)}))
    rc = RT.main(["--results", str(res), "--work_dir", str(work), "--db_paths", str(dbp), "--threads", "2"],
                 find=_finder(bindir))
    check("readtax CLI: finished job, sample with reads done", rc == 0 and
          (res / "samples" / "FS03" / "readtax" / "read_taxonomy.tsv").exists(), f"rc={rc}")
    i4 = RT.load_info(res / "samples" / "FS 04" / "readtax")
    check("readtax CLI: sample whose reads are gone says so, no table",
          not (res / "samples" / "FS 04" / "readtax" / "read_taxonomy.tsv").exists()
          and "reads" in i4.get("note", "").lower(), str(i4))
    # reads gone later: an earlier result is kept, not wiped
    keep = res / "samples" / "FS03" / "readtax"
    (work / "FS03" / "reads.filtered.fastq.gz").rename(work / "FS03" / "moved.fastq.gz")
    RT.main(["--results", str(res), "--work_dir", str(work), "--db_paths", str(dbp), "--samples", "FS03"],
            find=_finder(bindir))
    check("readtax CLI: reads gone → earlier result kept", (keep / "read_taxonomy.tsv").exists()
          and RT.load_info(keep).get("tools"))
    (work / "FS03" / "moved.fastq.gz").rename(work / "FS03" / "reads.filtered.fastq.gz")
    rc2 = RT.main(["--results", str(res), "--work_dir", str(work), "--db_paths", str(dbp),
                   "--samples", "FS 04"], find=_finder(bindir))
    check("readtax CLI: --samples limits the run; nothing to do → exit 1", rc2 == 1, f"rc={rc2}")
    # report picks it up
    import wgs_customer_report as R
    r = R.build_sample(res, "FS03", d / "pkg", ["html"])
    html = Path(r["html"]).read_text()
    check("readtax → report: heatmap + method, refs for Kraken 2/Bracken",
          "Read-based taxonomy was not run" not in html and "Kraken2+Bracken" in html
          and "Improved metagenomic analysis with Kraken 2" in html and "Bracken: estimating" in html)
    check("readtax → package: files in 4.Taxonomy/Read-based Taxonomy",
          (Path(r["dir"]) / "4.Taxonomy" / "Read-based Taxonomy" / "read_taxonomy.tsv").exists())


def t_readtax_wiring(d):
    m = (ROOT / "backend" / "main.py").read_text()
    w = (ROOT / "backend" / "ont_wgs_pipeline.py").read_text()
    f = (ROOT / "frontend-beta" / "index.html").read_text()
    check("readtax wiring: WGS step + option", "--skip_readtax" in w and "readtax.run_sample(filt" in w
          and "wgs_run_readtax" in m and 'cmd += ["--skip_readtax"]' in m and '"wgs_run_readtax"' in f)
    check("readtax wiring: finished-job endpoint goes through the queue",
          '@app.post("/wgs/{job_id}/readtax")' in m and "GATE.enter(tid" in m and "_submit_task(job_id, \"readtax\"" in m)
    check("readtax wiring: only samples with reads are sent", "if s not in gone" in m)
    check("readtax wiring: console panel + queue knows tasks", "wgsReadtaxPanel(id" in f
          and "/tasks/\" + encodeURIComponent(id) + \"/cancel" in f and "TASKS.get(jid)" in m)


TESTS += [t_readtax_wiring, t_readtax_parsers, t_readtax_table_heatmap, t_readtax_run_sample, t_readtax_cli_finished_job]


# ── TYGS import (tygs.py) ─────────────────────────────────────────────────────
GUID = "29ae139a-3a41-4571-a19b-809af76a3ec5"


def _dd(q, subj, d0, d4, d6, ct="U vs. T", gc="0.12"):
    return {"comp_type": ct, "query_genome": q, "subject_genome": subj, "digital_ddh_d0": d0,
            "confidence_interval_d0": f"[{d0 - 2} - {d0 + 2}]", "digital_ddh_d4": d4,
            "confidence_interval_d4": f"[{d4 - 2} - {d4 + 2}]", "digital_ddh_d6": d6,
            "confidence_interval_d6": f"[{d6 - 2} - {d6 + 2}]", "percent_gc_diff": gc, "id": None}


def fake_tygs(finished=True, genomes=("FS03_assembly", "OTHER_customer"), status_code=200, calls=None):
    B = "https://tygs.dsmz.de"
    me = f"<b>'{genomes[0]}.fasta'</b> "
    other = f"<b>'{genomes[1]}.fasta'</b> " if len(genomes) > 1 else None
    dd = [_dd(me, "<I>Bacillus</I> <I>subtilis</I> <a href=\"x\">DSM 10</a>", 85.1, 89.3, 88.0),
          _dd(me, "<I>Bacillus</I> <I>spizizenii</I> NRRL B-23049", 40.2, 55.5, 41.0)]
    if other:
        dd += [_dd(other, "<I>Escherichia</I> <I>coli</I> DSM 30083", 90.0, 95.0, 92.0),
               _dd(me, other, 20.0, 21.0, 20.0, ct="U vs. U")]
    tips = "&#39;Bacillus subtilis DSM 10&#39;:0.01,&#39;Bacillus spizizenii NRRL B-23049&#39;:0.02"
    users = ",".join(f"&#39;{g}&#39;:0.005" for g in genomes)
    nwk = f"(({users}),({tips})88:0.01);"
    R = {f"{B}/user_requests/show/status/{GUID}":
         (status_code, "application/json", json.dumps([{"guid": GUID, "finalization_date":
                                                         "2026-11-01T10:00:00.000+00:00" if finished else None,
                                                         "genomes_list": ", ".join(g + ".fasta" for g in genomes)}])),
         f"{B}/user_results/show/digital_ddh_json/{GUID}": (200, "application/json", json.dumps(dd)),
         f"{B}/user_results/show/identification_table/{GUID}": (500, "text/html", "<html>Server error</html>"),
         f"{B}/user_results/show/strain_table_json/{GUID}": (200, "application/json", "[]"),
         f"{B}/user_results/show/bibtex/{GUID}": (200, "text/plain", "@article{x, title={Bacillus}}"),
         f"{B}/user_phylogenies/16S/newick_by_guid/{GUID}": (200, "text/html", nwk),
         f"{B}/user_phylogenies/genome/newick_by_guid/{GUID}": (200, "text/html", nwk)}

    def get(url, timeout=60):
        if calls is not None:
            calls.append(url)
        return R.get(url, (404, "text/html", "not found"))
    return get


def t_tygs_guid(d):
    import tygs as TY
    check("tygs: GUID from bare text / result link / status URL",
          TY.parse_guid(GUID) == GUID and TY.parse_guid(f"https://tygs.dsmz.de/user_results/show?guid={GUID}") == GUID
          and TY.parse_guid(" " + GUID.upper() + "\n") == GUID)
    check("tygs: junk is not a GUID", TY.parse_guid("hello") is None and TY.parse_guid("") is None
          and TY.parse_guid("29ae139a-3a41-4571-a19b") is None)


def t_tygs_pending_and_errors(d):
    import tygs as TY
    calls = []
    st = TY.import_result(GUID, "FS03", d / "tygs", get=fake_tygs(finished=False, calls=calls), now=lambda: 1000.0)
    check("tygs: unfinished job → pending, only the status was asked", st["status"] == "pending"
          and len(calls) == 1 and TY.load(d / "tygs") is None, str(st))
    calls.clear()
    st = TY.import_result(GUID, "FS03", d / "tygs", get=fake_tygs(finished=False, calls=calls), now=lambda: 1000.0 + 600)
    check("tygs: asking again within 30 min does not call TYGS", st["status"] == "pending" and not calls
          and st.get("next_check", 0) >= 1000.0 + 1800, str(st))
    st = TY.import_result(GUID, "FS03", d / "tygs", get=fake_tygs(finished=True), now=lambda: 1000.0 + 1900)
    check("tygs: after 30 min it checks again and imports", st["status"] == "done", str(st))
    bad = None
    try:
        TY.import_result(GUID, "FS03", d / "t2", get=lambda u, timeout=60: (200, "application/json", "[]"))
    except TY.TygsError as e:
        bad = str(e)
    check("tygs: unknown GUID → clear error", bad and "not found" in bad.lower(), str(bad))
    down = None

    def boom(u, timeout=60):
        raise OSError("Network is unreachable")
    try:
        TY.import_result(GUID, "FS03", d / "t3", get=boom)
    except TY.TygsError as e:
        down = str(e)
    check("tygs: network down → clear error", down and "could not reach" in down.lower(), str(down))


def t_tygs_import(d):
    import tygs as TY
    st = TY.import_result(GUID, "FS03", d / "tygs", get=fake_tygs())
    t = TY.load(d / "tygs")
    check("tygs: imported for the sample's own genome", st["status"] == "done" and st["genome"] == "FS03_assembly"
          and t is not None, str(st))
    idn = t["identification"][0]
    check("tygs: Table 4.1 from the best d4 (≥ 79 % → known species)",
          idn["strain"] == "FS03_assembly" and idn["conclusion"] == "belongs to known species"
          and idn["species"] == "Bacillus subtilis" and "89.3" in idn["remark"], str(idn))
    check("tygs: dDDH rows = this genome vs type strains, best first",
          [r["d4"] for r in t["dddh"]] == ["89.3", "55.5"] and t["dddh"][0]["subject"] == "Bacillus subtilis DSM 10"
          and t["dddh"][0]["query"] == "FS03_assembly", str(t["dddh"]))
    check("tygs: trees unescaped, own tip renamed to the sample, other job genomes removed",
          "OTHER_customer" not in t["trees"]["genome"] and "FS03:" in t["trees"]["16S"]
          and "&#39;" not in t["trees"]["16S"], t["trees"]["16S"])
    vis = [f for f in (d / "tygs").rglob("*") if f.is_file() and not any(p.startswith(".") for p in
                                                                         f.relative_to(d / "tygs").parts)]
    leak = [f.name for f in vis if "OTHER_customer" in f.read_text(errors="replace")]
    check("tygs: no other genome of the TYGS job in the deliverable files", vis and not leak, str(leak))
    check("tygs: GUID and finalisation kept", t["guid"] == GUID and t["finalization_date"].startswith("2026-11-01"))
    # subspecies / new species thresholds
    check("tygs: d4 thresholds 70 / 79 %", TY.conclusion(80) == "belongs to known species"
          and TY.conclusion(75) == "belongs to known species, potential new subspecies"
          and TY.conclusion(69.9) == "potential new species")


def t_tygs_ambiguous(d):
    import tygs as TY
    get = fake_tygs(genomes=("isolate_A", "isolate_B"))
    err = None
    try:
        TY.import_result(GUID, "FS03", d / "tygs", get=get)
    except TY.TygsAmbiguous as e:
        err = e
    check("tygs: several genomes, none named like the sample → asks which", err is not None
          and sorted(err.names) == ["isolate_A", "isolate_B"], str(err))
    st = TY.import_result(GUID, "FS03", d / "tygs", genome="isolate_A", get=get)
    t = TY.load(d / "tygs")
    check("tygs: the chosen genome is used and shown as the sample",
          st["status"] == "done" and t["identification"][0]["strain"] == "FS03_assembly"
          and "isolate_B" not in t["trees"]["genome"], str(t["identification"]))
    err1 = None
    try:
        TY.import_result(GUID, "FS03", d / "t1", get=fake_tygs(genomes=("whatever_name",)))
    except TY.TygsAmbiguous as e:
        err1 = e
    check("tygs: a one-genome job under another name is not taken silently (wrong link?)",
          err1 is not None and err1.names == ["whatever_name"])
    one = TY.import_result(GUID, "FS03", d / "t1", genome="whatever_name", get=fake_tygs(genomes=("whatever_name",)))
    check("tygs: …but is taken once the BI confirms the genome", one["status"] == "done")
    calls = []
    TY.import_result(GUID, "FS03", d / "t1", get=fake_tygs(genomes=("whatever_name",), calls=calls))
    check("tygs: importing the same finished link again does not ask TYGS", calls == [], str(calls))


def t_tygs_cache_and_relink(d):
    import tygs as TY
    calls, cache = [], d / "cache"
    get = fake_tygs(finished=False, genomes=("FS03_assembly", "FS04_assembly"), calls=calls)
    TY.import_result(GUID, "FS03", d / "a", get=get, now=lambda: 100.0, cache=cache)
    TY.import_result(GUID, "FS04", d / "b", get=get, now=lambda: 200.0, cache=cache)
    check("tygs: one status check per GUID for all samples of a run", len(calls) == 1, str(calls))
    st = TY.import_result(GUID, "FS04", d / "b", get=get, now=lambda: 300.0, cache=cache)
    check("tygs: next check counted from the shared check", st["next_check"] == 100.0 + TY.MIN_RECHECK, str(st))
    # a finished import, then a new (pending) link: the old result leaves the report
    TY.import_result(GUID, "FS03", d / "c", get=fake_tygs())
    other = "11111111-2222-3333-4444-555555555555"
    TY.import_result(other, "FS03", d / "c", get=lambda u, timeout=60: (200, "application/json", json.dumps(
        [{"guid": other, "finalization_date": None, "genomes_list": "FS03_assembly.fasta"}])))
    check("tygs: a new pending link removes the old result", TY.load(d / "c") is None
          and not (d / "c" / "dddh.tsv").exists())


def t_tygs_tree_allowlist(d):
    import tygs as TY
    nwk = ("((&#39;FS03_assembly&#39;:0.01,&#39;Mystery genome x&#39;:0.01):0.01,"
           "(&#39;Bacillus subtilis DSM 10&#39;:0.01,&#39;Bacillus_spizizenii_NRRL_B-23049&#39;:0.02):0.01);")
    allow = {TY.norm("Bacillus subtilis DSM 10"), TY.norm("Bacillus spizizenii NRRL B-23049")}
    out = TY._prune(nwk, "FS03_assembly", allow, "FS03")
    check("tygs: tree keeps only this genome + compared type strains (label spelling tolerated)",
          "Mystery" not in out and "FS03:" in out and "spizizenii" in out, out)


def t_tygs_report(d):
    import tygs as TY
    import wgs_customer_report as R
    res = d / "results"
    make_wgs_run(res, samples=("FS03",))
    TY.import_result(GUID, "FS03", res / "samples" / "FS03" / "tygs", get=fake_tygs())
    r = R.build_sample(res, "FS03", d / "pkg", ["html"])
    html = Path(r["html"]).read_text()
    check("tygs → report: dDDH table, TYGS trees, GBDP/TYGS references",
          "pairwise dDDH values" in html and "89.3" in html and "from TYGS" in html
          and "Meier-Kolthoff" in html and "not imported" not in html and "OTHER_customer" not in html)
    pk = Path(r["dir"]) / "4.Taxonomy" / "TYGS"
    names = sorted(p.name for p in pk.rglob("*") if p.is_file())
    check("tygs → package: TYGS folder has tables + trees, no raw job files",
          "dddh.tsv" in names and "genome_tree.nwk" in names and not any(n.startswith(".") for n in names)
          and not any("status" in n for n in names), str(names))


def t_tygs_wiring(d):
    m = (ROOT / "backend" / "main.py").read_text()
    f = (ROOT / "frontend-beta" / "index.html").read_text()
    check("tygs wiring: assembly download + import endpoints", '@app.get("/wgs/{job_id}/tygs_assembly/{sample}")' in m
          and '@app.post("/wgs/{job_id}/tygs/{sample}")' in m and "TygsAmbiguous" in m)
    check("tygs wiring: console panel with the DSMZ note", "wgsTygsPanel(" in f and "DSMZ" in f
          and "tygs.dsmz.de" in f)


TESTS += [t_tygs_guid, t_tygs_pending_and_errors, t_tygs_import, t_tygs_ambiguous, t_tygs_cache_and_relink,
          t_tygs_tree_allowlist, t_tygs_report, t_tygs_wiring]


# ── result packages for download (packages.py) ───────────────────────────────
def _zip_texts(z: Path) -> dict:
    import zipfile
    out = {}
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            out[n] = zf.read(n).decode("utf-8", "replace") if n.endswith((".tsv", ".csv", ".txt", ".html",
                                                                            ".json", ".fasta", ".nwk")) else ""
    return out


def t_pkg_wgs(d):
    import packages as PK
    res = d / "res"
    make_wgs_run(res, samples=("FS03", "FS 04"))
    (res / "wgs_results.json").write_text("{}")
    check("packages: kind + samples of a WGS run", PK.detect_kind(res) == "wgs"
          and PK.list_samples(res, "wgs") == ["FS 04", "FS03"])
    seen = []
    m = PK.build(res, None, "wgs", ["FS03", "FS 04"], True, d / "pk", job_label="TH-BIOENTIST-00002",
                 progress=lambda i, n, lab: seen.append((i, n, lab)), formats=("html",))
    names = sorted(f["name"] for f in m["files"])
    check("packages: one zip per sample, named by sample", names == ["FS03.zip", "FS_04.zip"], str(names))
    t3 = _zip_texts(d / "pk" / "FS03.zip")
    check("packages: WGS zip = sample folder with 1–8 + report",
          all(any(n.startswith(f"FS03/{f}/") for n in t3) for f in FOLDERS) and "FS03/Report_FS03.html" in t3,
          str(sorted(t3)[:12]))
    leak = [n for n, v in t3.items() if "FS 04" in n or "FS_04" in n or "FS 04" in v or "FS_04" in v]
    check("packages: FS03 zip holds nothing of FS 04", not leak, str(leak[:5]))
    t4 = _zip_texts(d / "pk" / "FS_04.zip")
    leak = [n for n, v in t4.items() if "FS03" in n or "FS03" in v]
    check("packages: FS 04 zip holds nothing of FS03", not leak, str(leak[:5]))
    check("packages: progress n / N reported", any(lab.startswith("Packaging FS 04 (2/2)") for _, _, lab in seen),
          str(seen))
    check("packages: no work folders / .part files left", not (d / "pk" / "_work").exists()
          and not list((d / "pk").glob("*.part")) and (d / "pk" / "manifest.json").exists())
    m1 = PK.build(res, None, "wgs", ["FS 04"], False, d / "one", job_label="TH-BIOENTIST-00002", formats=("html",))
    t1 = _zip_texts(Path(m1["files"][0]["path"]))
    check("packages: one zip with the chosen samples only", m1["files"][0]["name"] == "TH-BIOENTIST-00002_results.zip"
          and any(n.startswith("TH-BIOENTIST-00002_results/FS_04/") for n in t1)
          and not any("FS03" in n for n in t1), str(sorted(t1)[:5]))
    try:
        PK.build(res, None, "wgs", ["nope"], True, d / "x")
        check("packages: unknown sample refused", False)
    except ValueError:
        check("packages: unknown sample refused", True)
    (res / "packages").mkdir(exist_ok=True)
    (res / "packages" / "old.zip").write_text("x")
    z = PK.zip_whole_job(res, d / "all.zip")
    import zipfile
    nn = zipfile.ZipFile(z).namelist()
    check("packages: whole-job zip on disk, packages folder not inside", "samples/FS03/result.json" in nn
          and not any(n.startswith("packages/") for n in nn), str(nn[:5]))


AMP_FOLDERS = ["1.QcReads", "2.Taxonomy", "3.Diversity", "4.Tables"]
AMP_HEADINGS = ["amplicon sequencing report", "Sample: ", "1. Quality control and clean data", "Results 1.QcReads",
                "2. Taxonomic profiles", "Results 2.Taxonomy", "3. Diversity", "Results 3.Diversity", "4. Tables",
                "Results 4.Tables", "References"]


def _amp_tools(d):
    import test_amp_fast as TA
    fq, fp = TA._fake_qc_tools(d / "qc")
    return {"fastqc": fq, "fastp": fp, "Rscript": TA._fake_rscript(d), "ktImportText": None, "NanoPlot": None}


def _amp_missing(html, heads):
    pos, miss = 0, []
    for h in heads:
        k = html.find(h, pos)
        if k < 0:
            miss.append(h)
        else:
            pos = k + len(h)
    return miss


def t_amp_package_emu(d):
    import amp_package as AP
    import test_amp_fast as TA
    run = TA.make_emu_run(d / "run")
    for s_ in ("L21", "PGU-16S-full"):
        _w(run / "emu_output" / s_ / f"{s_}_rel-abundance-threshold-0.0001.tsv", f"tax_id\tabundance\n{s_}\t1\n")
        _w(run / "emu_output" / s_ / "emu.log", "log")
        _w(run / "qc_filtered" / f"{s_}.fastq.gz", "@r\nACGT\n+\nIIII\n")
    m = AP.build_sample(run, d / "up", "L21", d / "out", ("html",), tools=_amp_tools(d))
    root = Path(m["dir"])
    check("amplicon package: numbered folders + report in the root", all((root / f).is_dir() for f in AMP_FOLDERS)
          and (root / "Report_L21.html").exists(), str(sorted(p.name for p in root.iterdir())))
    html = (root / "Report_L21.html").read_text()
    miss = _amp_missing(html, AMP_HEADINGS)
    check("amplicon report: WGS-style title, sections in order, each with its Results line", not miss, str(miss))
    check("amplicon report: no cover / lab summary blocks", "Client Name" not in html and "Report Summary:" not in html
          and "Images_and_Outfiles" not in html and "16S rRNA amplicon sequencing report" in html)
    want = ["1.QcReads/read_tracking.tsv", "1.QcReads/Box_plot.png", "1.QcReads/filtered_reads/L21.fastq.gz",
            "2.Taxonomy/outfiles/expr.abundance.all.txt", "2.Taxonomy/outfiles/expr.biom",
            "2.Taxonomy/Bar charts/plot_bar_genus.png", "2.Taxonomy/Heatmap/Heat_map_top50.png",
            "2.Taxonomy/Krona/krona_static.png", "2.Taxonomy/Emu/L21_rel-abundance-threshold-0.0001.tsv",
            "3.Diversity/Rarefaction_curves.png", "3.Diversity/alpha_diversity.tsv"]
    gone = [w for w in want if not (root / w).exists()]
    check("amplicon package: files in their folders", not gone, str(gone))
    check("amplicon package: no logs / work files inside", not list(root.rglob("*.log"))
          and not [p for p in root.rglob("*") if any(x.startswith(".") or x == "_work" for x in p.relative_to(root).parts)])
    texts = {p: p.read_text(errors="replace") for p in root.rglob("*") if p.is_file() and p.suffix in
             (".html", ".tsv", ".txt", ".biom", ".json")}
    leak = [str(p.relative_to(root)) for p, v in texts.items() if "PGU" in v or "PGU" in str(p)]
    check("amplicon package: nothing of the other sample", not leak, str(leak[:5]))
    m2 = AP.build_sample(run, d / "up", "PGU-16S-full", d / "out2", ("html",), tools=_amp_tools(d),
                         include_reads=False)
    r2 = Path(m2["dir"])
    h2 = (r2 / "Report_PGU-16S-full.html").read_text()
    check("amplicon package: reads left out on request, the report does not name them",
          not (r2 / "1.QcReads" / "filtered_reads").exists() and "filtered_reads" not in h2)


def t_amp_package_dada2(d):
    import amp_package as AP
    import test_amp_fast as TA
    run = TA.make_amp_run(d / "run")
    up = d / "up"
    up.mkdir()
    for m in ("R1", "R2"):
        (up / f"S2-12S_{m}_001.fastq").write_text("@a\nACGT\n+\nIIII\n")
    m = AP.build_sample(run, up, "S2-12S", d / "out", ("html",), tools=_amp_tools(d))
    root = Path(m["dir"])
    html = (root / "Report_S2-12S.html").read_text()
    check("12S package: sections + Results lines", not _amp_missing(html, AMP_HEADINGS), str(_amp_missing(html, AMP_HEADINGS)))
    check("12S package: FastQC / fastp outputs in 1.QcReads, ASV FASTA in 2.Taxonomy",
          (root / "1.QcReads" / "FastQC").is_dir() and (root / "1.QcReads" / "fastp" / "fastp.json").exists()
          and (root / "2.Taxonomy" / "outfiles" / "expr.asv.fasta").exists()
          and (root / "1.QcReads" / "plotQualityProfile.pdf").exists())
    check("12S package: BLAST / NCBI check table of this sample's ASVs",
          (root / "2.Taxonomy" / "BLAST" / "blast_ncbi_check.tsv").exists())
    leak = [str(p.relative_to(root)) for p in root.rglob("*") if p.is_file() and p.suffix in (".html", ".tsv", ".txt",
            ".fasta", ".biom") and ("F4 CytB.1" in p.read_text(errors="replace") or "F4_CytB" in p.read_text(errors="replace"))]
    check("12S package: nothing of F4 CytB.1", not leak, str(leak[:5]))
    m3 = AP.build_sample(run, d / "noup", "F4 CytB.1", d / "out3", ("html",), tools=_amp_tools(d))
    h3 = Path(m3["html"]).read_text()
    check("12S package: raw reads missing → said so, folder still there", "were not found" in h3
          and (Path(m3["dir"]) / "1.QcReads").is_dir())


def t_pkg_amplicon(d):
    import packages as PK
    import test_amp_fast as TA
    run = TA.make_emu_run(d / "run")
    for s_ in ("L21", "PGU-16S-full"):
        _w(run / "emu_output" / s_ / f"{s_}_rel-abundance-threshold-0.0001.tsv", f"tax_id\tabundance\n{s_}\t1\n")
        _w(run / "qc_filtered" / f"{s_}.fastq.gz", "@r\nACGT\n+\nIIII\n")
    tools = _amp_tools(d)
    check("packages: kind of an amplicon run", PK.detect_kind(run) == "amplicon"
          and PK.list_samples(run, "amplicon") == ["L21", "PGU-16S-full"])
    m = PK.build(run, d / "up", "amplicon", ["L21", "PGU-16S-full"], True, d / "pk", job_label="TH-JOB",
                 include_reads=True, tools=tools, formats=("html",))
    names = sorted(f["name"] for f in m["files"])
    check("packages: amplicon one zip per sample, named by sample", names == ["L21.zip", "PGU-16S-full.zip"], str(names))
    t = _zip_texts(d / "pk" / "L21.zip")
    check("packages: amplicon zip = the numbered package", all(any(n.startswith(f"L21/{f}/") for n in t)
                                                              for f in AMP_FOLDERS) and "L21/Report_L21.html" in t,
          str(sorted(t)[:10]))
    leak = [n for n, v in t.items() if "PGU" in n or "PGU" in v or "TH-JOB" in v]
    check("packages: L21 zip holds nothing of PGU-16S-full (or the job name)", not leak, str(leak[:5]))
    m2 = PK.build(run, d / "up", "amplicon", ["PGU-16S-full"], False, d / "one", job_label="TH-JOB", tools=tools,
                  formats=("html",))
    t2 = _zip_texts(Path(m2["files"][0]["path"]))
    check("packages: one amplicon zip for a subset holds only it",
          not any("L21" in n or "L21" in v for n, v in t2.items())
          and any(n.startswith("TH-JOB_results/PGU-16S-full/2.Taxonomy/") for n in t2), str(sorted(t2)[:8]))


def t_pkg_wiring(d):
    m = (ROOT / "backend" / "main.py").read_text()
    f = (ROOT / "frontend-beta" / "index.html").read_text()
    check("packages wiring: endpoints through the queue", '@app.post("/packages/{job_id}")' in m
          and '@app.get("/packages/{job_id}/file/{name}")' in m and '_submit_task(job_id, "package"' in m)
    i = m.find('def download_zip(')
    check("packages wiring: whole-job zip written to disk, not RAM", i > 0 and "BytesIO" not in m[i:i + 2500]
          and "zip_whole_job" in m[i:i + 2500])
    check("packages wiring: download dialog for every pipeline", "openDownloadDialog(" in f
          and "One zip per sample" in f and "/packages/" in f)


def t_pkg_name_collisions(d):
    import packages as PK
    res = d / "res"
    make_wgs_run(res, samples=("FS 04", "FS_04", "ก้อง", "ก่อง"))
    (res / "wgs_results.json").write_text("{}")
    check("packages: Thai names stay distinct", PK.safe("ก้อง") != PK.safe("ก่อง") and PK.safe("ก้อง") == "ก้อง")
    m = PK.build(res, None, "wgs", ["FS 04", "FS_04", "ก้อง", "ก่อง"], True, d / "pk", formats=("html",))
    names = [f["name"] for f in m["files"]]
    check("packages: colliding names get their own zip", len(set(names)) == 4
          and all((d / "pk" / n).exists() for n in names), str(names))
    by = {f["samples"][0]: f["path"] for f in m["files"]}
    t = _zip_texts(Path(by["FS_04"]))
    rep = [n for n in t if n.endswith(".html")]
    check("packages: each colliding zip holds its own sample", rep and "Sample: FS_04" in t[rep[0]]
          and "Sample: FS 04" not in t[rep[0]], str(rep))


TESTS += [t_amp_package_emu, t_amp_package_dada2, t_pkg_wgs, t_pkg_amplicon, t_pkg_name_collisions, t_pkg_wiring]


def main() -> int:
    print("ONT WGS package / read taxonomy / TYGS / downloads")
    with tempfile.TemporaryDirectory() as td:
        for i, t in enumerate(TESTS):
            d = Path(td) / f"t{i}"
            d.mkdir()
            try:
                t(d)
            except Exception:
                check(t.__name__, False, traceback.format_exc()[-900:])
    return len(FAILS)


if __name__ == "__main__":
    n = main()
    print(f"\n{'ALL OK' if not n else f'{n} FAILED'}")
    sys.exit(1 if n else 0)
