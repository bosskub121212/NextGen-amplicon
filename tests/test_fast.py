#!/usr/bin/env python3
"""Fast checks (seconds, Python only) — run by push_github.sh before every push.

They cover the places where a past bug changed results without any error:
  * rank names read from a reference's headers (SILVA 138 / 144, PR2, MIDORI2)
    in the VSEARCH pipeline,
  * the report reading every pipeline's table layout (DADA2, ITS, COX1 old and
    new, PacBio, VSEARCH, Emu) and building a report for each,
  * COX1 results made before v2.17 having their shifted ranks put back,
  * database labels, sample column headers, paired-end file detection.

The full end-to-end mock-community test needs R and the databases and runs on
the analysis machine: python3 tests/run_mock_tests.py
"""
from __future__ import annotations

import csv
import json
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))

FAILS = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILS.append(name)


def write_csv(p: Path, header, rows):
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


SEQS = ["TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTC" + "A" * i
        for i in range(6)]


def make_result(d: Path, kind: str):
    """A tiny result folder in the layout each pipeline writes."""
    d.mkdir(parents=True)
    S = ["s1", "s2", "s3"]
    counts = [[50, 10, 0], [30, 30, 30], [5, 40, 20], [10, 0, 40], [3, 3, 3], [2, 17, 7]]
    if kind in ("dada2", "its", "pacbio", "cox1"):
        ranks = {"dada2": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"],
                 "pacbio": ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"],
                 "its": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                 "cox1": ["Phylum", "Class", "Order", "Family", "Genus", "Species"]}[kind]
        lin = {"dada2": [["Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", f"Genus{i}"]
                         for i in range(6)],
               "pacbio": [["Bacteria", "Bacillati", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae",
                           f"Genus{i}"] for i in range(6)],
               "its": [["Fungi", "Ascomycota", "Eurotiomycetes", "Eurotiales", "Aspergillaceae",
                        f"Genus{i}", f"Genus{i} niger"] for i in range(6)],
               "cox1": [["Arthropoda", "Insecta", "Diptera", "Muscidae", f"Genus{i}", f"Genus{i} domestica"]
                        for i in range(6)]}[kind]
        write_csv(d / "asv_table.csv", S + ["sequence"], [c + [s] for c, s in zip(counts, SEQS)])
        write_csv(d / "taxonomy_table.csv", [""] + ranks, [[s] + l for s, l in zip(SEQS, lin)])
        trk = {"dada2": ["input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"],
               "its": ["input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"],
               "pacbio": ["input", "filtered", "denoised", "nonchim"],
               "cox1": ["input", "filtered", "denoised_F", "denoised_R", "merged", "nonchim",
                        "length_filtered", "codon_filtered", "lulu_curated"]}[kind]
        write_csv(d / "read_tracking.csv", [""] + trk,
                  [[s] + [200 - 10 * k for k in range(len(trk))] for s in S])
        pipe = {"dada2": "", "its": "its", "pacbio": "pacbio", "cox1": "cox1"}[kind]
        (d / "summary.json").write_text(json.dumps({"pipeline": pipe} if pipe else {}))
        (d / "run_params.json").write_text(json.dumps(
            {"marker": {"dada2": "16S", "its": "ITS1", "pacbio": "PACBIO", "cox1": "COX1"}[kind]}))
    elif kind == "cox1_old":
        # v2.16 and earlier: samples as rows, ASV ids as columns, ranks shifted
        ids = [f"ASV{i+1}" for i in range(6)]
        write_csv(d / "asv_table.csv", ["sample_id"] + ids,
                  [[s] + [counts[i][j] for i in range(6)] for j, s in enumerate(S)])
        write_csv(d / "taxonomy_table.csv", ["ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"],
                  [[a, "Arthropoda_6656", "Insecta_50557", "Diptera_7147", "Muscidae_7366",
                    f"Musca{i}_7369", f"Musca{i} domestica_7370"] for i, a in enumerate(ids)])
        (d / "pipeline_summary.json").write_text("{}")
        (d / "run_params.json").write_text(json.dumps({"marker": "COX1"}))
    elif kind in ("vsearch", "emu"):
        ids = [f"OTU_{i+1}" if kind == "vsearch" else str(1000 + i) for i in range(6)]
        write_csv(d / "asv_table.csv", [""] + S, [[a] + c for a, c in zip(ids, counts)])
        write_csv(d / "taxonomy.csv", ["ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                  [[a, "Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", f"Genus{i}",
                    f"Genus{i} sp." if kind == "emu" else ""] for i, a in enumerate(ids)])
        write_csv(d / "read_tracking.csv", ["sample", "input", "qc_filtered", "primer_trimmed", "final"],
                  [[s, 200, 190, 180, sum(c[j] for c in counts)] for j, s in enumerate(S)])
        (d / "summary.json").write_text(json.dumps(
            {"pipeline": "qiime2_vsearch" if kind == "vsearch" else "emu",
             **({"emu_counts": "estimated_reads"} if kind == "emu" else {})}))
        (d / "run_params.json").write_text(json.dumps({"marker": "16S" if kind == "vsearch" else "ONT-16S"}))
    elif kind == "shotgun":
        ids = [str(2000 + i) for i in range(6)]
        write_csv(d / "asv_table.csv", ["tax_id"] + S, [[a] + c for a, c in zip(ids, counts)])
        write_csv(d / "taxonomy.csv", ["tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
                  [[a, "Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", f"Genus{i}",
                    f"Genus{i} alpha"] for i, a in enumerate(ids)])
        write_csv(d / "read_tracking.csv", ["sample", "input", "after_qc", "host_removed", "classified",
                                            "unclassified", "host_classified", "final"],
                  [[s, 200, 190, 0, 180, 10, 0, sum(c[j] for c in counts)] for j, s in enumerate(S)])
        (d / "summary.json").write_text(json.dumps({"pipeline": "shotgun", "platform": "illumina",
                                                    "tools": {"kraken2": "2.17.2"}}))
        (d / "run_params.json").write_text(json.dumps({"marker": "SHOTGUN", "sequencerType": "illumina"}))


def main():
    print("VSEARCH lineage parsing")
    import qiime2_vsearch_pipeline as q
    t = q.parse_lineage("Bacteria;Bacillati;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Limosilactobacillus;")
    check("SILVA 144: genus is the genus", t.get("genus") == "Limosilactobacillus" and t.get("phylum") == "Bacillota", t)
    t = q.parse_lineage("Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;")
    check("SILVA 138: genus is the genus", t.get("genus") == "Lactobacillus" and t.get("phylum") == "Firmicutes", t)
    t = q.parse_lineage("Arthropoda_6656;Insecta_50557;order_Muscidae_7366;Muscidae_7366;Musca_7369;Musca domestica_7370;")
    check("MIDORI2: phylum..species, taxids and placeholders removed",
          t == {"phylum": "Arthropoda", "class": "Insecta", "family": "Muscidae", "genus": "Musca",
                "species": "Musca domestica"}, t)
    t = q.parse_lineage("Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae;Gymnodinium;Gymnodinium_catenatum")
    check("PR2: 9 levels, Division as phylum", t.get("phylum") == "Alveolata" and t.get("genus") == "Gymnodinium", t)
    check("paired file names", q._PAIR_RE.match("S1_S1_L001_R2_001.fastq.gz").group("mate") == "2"
          and q._PAIR_RE.match("mock1_R1.fastq.gz").group("name") == "mock1")

    print("Database labels")
    import provenance as pv
    for n, lab in [("silva_nr99_v144_toGenus_trainset.fa.gz", "SILVA 144"),
                   ("silva_nr99_v138.1_train_set.fa.gz", "SILVA 138.1"),
                   ("sh_general_release_dynamic_19.02.2025.fasta", "UNITE general release 2025-02-19"),
                   ("pr2_version_5.1.1_SSU_dada2.fasta.gz", "PR2 5.1.1"),
                   ("MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz", "MIDORI2 GenBank 273 (COI)")]:
        check(f"label {n}", pv.db_label(n) == lab, pv.db_label(n))

    print("Report on every pipeline's layout")
    import report_builder as rb
    check("sample headers distinct", rb._col_heads(["Bacteria_Pre-R1", "Bacteria_Pre-R2", "Bacteria_R1"])
          == ["Pre-R1", "Pre-R2", "R1"])
    with tempfile.TemporaryDirectory() as td:
        for kind in ("dada2", "its", "cox1", "cox1_old", "pacbio", "vsearch", "emu"):
            d = Path(td) / kind
            make_result(d, kind)
            try:
                run = rb.Run(d, set())
                an = rb.Analysis(run)
                html = rb.build_profile_html(run)
                gen = an.rank_pct("Genus")
                exp_kind = {"cox1_old": "cox1", "vsearch": "vsearch"}.get(kind, kind)
                check(f"{kind}: read as {run.kind}", run.kind == exp_kind, run.kind)
                check(f"{kind}: 3 samples, 6 features", len(run.samples) == 3 and len(run.mat) == 6,
                      f"{run.samples} {len(run.mat)}")
                check(f"{kind}: genera are genera", all(g.startswith(("Genus", "Musca")) or g == "Unclassified"
                                                        for g in gen), list(gen)[:6])
                check(f"{kind}: report builds", "<h1>" in html and len(html) > 5000)
                if kind == "cox1_old":
                    check("cox1_old: shifted ranks restored",
                          getattr(run, "rank_shift_fixed", False) and "Species" in run.rank_cols
                          and next(iter(run.tax.values()))["Phylum"] == "Arthropoda",
                          next(iter(run.tax.values())))
                if kind == "its":
                    check("its: no reagent / archaea sections", "Reagent-associated" not in html
                          and "Archaeal community" not in html)
            except Exception:
                check(f"{kind}: report pipeline", False, traceback.format_exc(limit=3))

    print("Read QC, controls, Emu parallelism (v2.18)")
    import read_qc
    # Q36 for 250 bp, then falling to Q15: the cut sits where the median drops under Q25
    prof = [[30, 34, 36, 37, 35.0]] * 200 + [[10, 12, 15, 18, 15.0]] * 50
    P = {"median_len": 250, "per_position": prof}
    sg = read_qc.suggest_trunc(P, P, "CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC", 420)
    check("truncLen: quality cut, pairs still overlap", sg["truncLen_F"] + sg["truncLen_R"] >= 440
          and sg["truncLen_F"] <= 250 - 17 - 3, sg)
    sg2 = read_qc.suggest_trunc(P, P, "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC", None)
    check("truncLen: ITS is not truncated", sg2["truncLen_F"] == 0 and sg2["truncLen_R"] == 0, sg2)
    import controls
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "dada2"
        make_result(d, "dada2")
        run = rb.Run(d, set())
        s0, s1 = run.samples[0], run.samples[1]
        run.params["controlSamples"] = [s0]
        run.params["mockSamples"] = [s1]
        negs, mocks = controls.identify(run)
        check("controls: named in the run", negs == [s0] and mocks == [s1], (negs, mocks))
        an = rb.Analysis(run)
        ns = controls.negative_summary(an, negs, mocks)
        check("controls: blank summary", ns["negatives"][0]["sample"] == s0 and ns["genera"], ns)
        ms = controls.mock_summary(an, mocks, "zymo_d6300")
        check("controls: mock compared", len(ms["samples"]) == 1
              and 0 <= ms["samples"][0]["bray_curtis"] <= 1, ms)
        html = rb.build_profile_html(run)
        check("controls: report section", "Controls" in html)
    import emu_pipeline as ep
    par, per = ep.emu_parallelism({k: "" for k in "abcd"}, 8, 2)
    check("emu: requested parallelism honoured", par == 2 and per == 4, (par, per))
    par, per = ep.emu_parallelism({"a": ""}, 8, 0)
    check("emu: one sample runs alone", par == 1 and per == 8, (par, per))

    print("ONT-WGS offline type strains")
    import wgs_taxonomy as WT
    rep = lambda acc, name, lvl, n50: {"accession": acc, "organism": {"organism_name": name},
                                        "assembly_info": {"assembly_level": lvl},
                                        "assembly_stats": {"contig_n50": n50, "total_sequence_length": 4e6}}
    picked = WT.pick_type_strains([rep("GCA_1.1", "Bacillus subtilis DSM 10", "Contig", 50000),
                                   rep("GCF_2.1", "Bacillus subtilis DSM 10", "Complete Genome", 4e6),
                                   rep("GCF_3.1", "Bacillus sp. X", "Complete Genome", 4e6)], "Bacillus")
    check("type strains: one per species, complete genome first, no 'sp.'",
          [r["accession"] for r in picked] == ["GCF_2.1"], picked)
    import wgs_typestrain_db  # noqa: F401  (imports cleanly)

    print("Shotgun metagenome (beta)")
    import shotgun_pipeline as SG
    with tempfile.TemporaryDirectory() as td:
        rp = Path(td) / "k.report"
        rp.write_text("10.00\t10\t10\tU\t0\tunclassified\n"
                      "90.00\t90\t0\tR\t1\troot\n"
                      "90.00\t90\t0\tD\t2\t  Bacteria\n"
                      "90.00\t90\t0\tP\t1239\t    Bacillota\n"
                      "90.00\t90\t0\tG\t1386\t      Bacillus\n"
                      "60.00\t60\t60\tS\t1423\t        Bacillus subtilis\n"
                      "30.00\t30\t30\tS\t1396\t        Bacillus cereus\n")
        lin, sp, un, cl = SG.parse_report(rp)
        check("shotgun: Kraken 2 report lineage", lin[1423]["Genus"] == "Bacillus"
              and lin[1423]["Kingdom"] == "Bacteria" and sp[1396] == 30 and un == 10 and cl == 90, (lin, sp))
        check("shotgun: Bracken read length", SG.pick_read_len([100, 150, 250], 148) == 100
              and SG.pick_read_len([150, 1000], 5000) == 1000)
        d = Path(td) / "sg"
        make_result(d, "shotgun")
        run = rb.Run(d, set())
        html = rb.build_profile_html(run)
        check("shotgun: read as shotgun, own report", run.kind == "shotgun"
              and "Shotgun Metagenome Sequencing Report" in html and "Bracken" in html)
        # resistome: catalogue headers, family pooling, breadth, host from the read's taxid
        import shotgun_amr as SA
        adb = Path(td) / "amrdb"
        adb.mkdir()
        (adb / "fam.tab").write_text(
            "#node_id\tparent_node_id\tgene_symbol\thmm_id\ta\tb\tc\td\te\tf\tg\th\treportable\ttype"
            "\tsubtype\tclass\tsubclass\tfamily_name\n"
            "blaTEM\tbla-A\tblaTEM\t-\t0\t0\t0\t0\t0\t0\t0\t0\t2\tAMR\tAMR\tBETA-LACTAM\tBETA-LACTAM"
            "\tTEM family class A beta-lactamase\n")
        (adb / "AMR_CDS").write_text(">1|WP_1|NG_1|1|1|blaTEM-1|blaTEM|TEM-1 NG_1:1-1000\n" + "A" * 1000 + "\n"
                                     ">2|WP_2|NG_2|1|1|blaTEM-116|blaTEM|TEM-116 NG_2:1-1000\n" + "A" * 1000 + "\n")
        cat = SA.load_catalog(adb)
        paf = Path(td) / "a.paf"
        t1, t2 = "1|WP_1|NG_1|1|1|blaTEM-1|blaTEM|TEM-1", "2|WP_2|NG_2|1|1|blaTEM-116|blaTEM|TEM-116"
        paf.write_text("".join(f"r{i}/1\t150\t0\t150\t+\t{t1 if i % 3 else t2}\t1000\t{i * 100}\t{i * 100 + 150}"
                               f"\t149\t150\t60\tde:f:0.005\n" for i in range(8))
                       + "x/1\t150\t0\t150\t+\t" + t1 + "\t1000\t0\t150\t120\t150\t60\tde:f:0.2\n")
        rows, reads = SA.summarise(paf, cat, "illumina", 1_000_000, 0.9, 60)
        tree = SA.report_tree(rp)
        SA.assign_hosts(rows, reads, {f"r{i}": 1423 for i in range(8)}, tree)
        check("shotgun resistome: alleles pooled per family, low identity dropped, host named",
              cat[t1]["class"] == "BETA-LACTAM" and len(rows) == 1 and rows[0]["gene"] == "blaTEM"
              and rows[0]["reads"] == 8 and rows[0]["breadth"] == 85.0
              and rows[0]["host_taxon"] == "Bacillus subtilis", rows)
        SA.write_tables(d, {"s1": rows, "s2": [], "s3": rows})
        sj = json.loads((d / "summary.json").read_text())
        sj["amr"] = {"db_version": "test", "min_identity": 90, "min_breadth": 60, "gene_families": 1,
                     "amr_families": 1}
        (d / "summary.json").write_text(json.dumps(sj))
        html = rb.build_profile_html(rb.Run(d, set()))
        check("shotgun report: resistome section", "resistome" in html and "blaTEM" in html
              and "AMRFinderPlus" in html, "")
    import build_16s_typestrain_db as B16
    check("16S multicopy: alignment-defined gene ends", B16.merge([(1, 5), (100, 200), (150, 300)], 10)
          == [(1, 5), (100, 300)])

    print("Console")
    page = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
    for k in ("mockSamples", "maaslin_fixed", "ont_parallel", "run_fastqc", "/qc/", "renderSettings",
              "/license/activate", "/update/check", "/results/cleanup", "orientation-check", "mpAdd",
              "/wgs/typestrains", "tsBanner", "wgDirIn", "shotgunSpec", "sg_kraken_db", "sg_amr", "sg_mags"):
        check(f"console offers {k}", k in page)
    main_src = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
    check("console served at /, old UI at /classic", '@app.get("/classic"' in main_src
          and "frontend-beta" in main_src)

    print()
    if FAILS:
        print(f"{len(FAILS)} check(s) failed: {', '.join(FAILS)}")
        return 1
    print("all fast checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
