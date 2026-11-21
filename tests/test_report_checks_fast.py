#!/usr/bin/env python3
"""Single-run reports of the 3-targets eDNA project (12S / 18S PR2 / COI), checked one by one:
cutadapt -n 2 in the method, COX1 step order and truncLen, LULU not installed, 18S off-target
table and marker label, unnamed non-target BLAST reads, human DNA inside a vertebrate limit."""
from __future__ import annotations
import json, sys, tempfile, traceback
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend")); sys.path.insert(0, str(ROOT / "tests"))
FAILS = []
def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  {detail}"))
    if not cond:
        FAILS.append(name)


def main():
    print("Single-run report checks (12S / 18S / COI)")
    try:
        import report_builder as RB
        from test_multimarker_fast import make, names
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            # 12S, run with v2.33.1 (cutadapt -n 2), BLAST: human (target under a vertebrate limit) + unnamed non-target
            d12 = make(td / "12s", "12S", names("{s}-{r}-Riaz12S", "{s}-CTRiaz12S"),
                       ["Phylum", "Class", "Order", "Family", "Genus", "Species"],
                       [(["Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "Homo sapiens"], None),
                        (["", "", "", "", "", ""], {"KK-S": 900})],
                       params={"mito_gene": "12S", "primer_f": "TTAGATACCCCACTATGC", "primer_r": "TAGAACAGGCTCCTCTAG",
                               "tax_method": "blast", "blast_limit": "vertebrates", "discardUntrimmed": True})
            (d12 / "provenance.json").write_text(json.dumps({"app_version": "2.33.1"}))
            seqs = [l.split(",")[0] for l in (d12 / "taxonomy_table.csv").read_text().splitlines()[1:]]
            (d12 / "blast").mkdir()
            (d12 / "blast" / "taxonomy_blast.csv").write_text(
                "id,sequence,Phylum,Class,Order,Family,Genus,Species,best_hit,pident,target,agreement,check\n"
                f"ASV1,{seqs[0]},Chordata,Mammalia,Primates,Hominidae,Homo,Homo sapiens,Homo sapiens,100,yes,agree,\n"
                f"ASV2,{seqs[1]},,,,,,,Haliclona sp. (Porifera),72.5,no,ncbi_no_hit,\n")
            (d12 / "blast" / "blast_info.json").write_text(json.dumps({"limit": "vertebrates"}))
            h12 = RB.build_profile_html(RB.Run(d12, set(), "3 Targets (12S)"))
            check("12S method: the cutadapt command carries -n 2 (run made with v2.33.1)",
                  "-O 3 -n 2 -m 50" in h12 and "two rounds" in h12, "")
            check("12S: non-target reads without a name are listed with their closest hit",
                  "Not named (closest: Haliclona sp. (Porifera))" in h12)
            check("12S: human DNA inside the vertebrate limit is said", "Human DNA (<i>Homo sapiens</i>) made up" in h12
                  and "counts as vertebrate DNA" in h12)
            check("one taxon: 'the 1 most abundant phylum', not 'phyla'", "1 most abundant phyla" not in h12)
            # an older run (no provenance version): no -n 2 claimed
            (d12 / "provenance.json").write_text(json.dumps({"app_version": "2.32.0"}))
            check("run made before v2.32.1: no -n 2 claimed", "-n 2" not in RB.build_profile_html(RB.Run(d12, set())))

            # COI: truncLen_cox1, length filter before chimera removal, LULU not installed
            dco = make(td / "coi", "COX1", names("{s}-{r}_COI-XT", "{s}-CT_COI-XT"),
                       ["Phylum", "Class", "Order", "Family", "Genus", "Species"],
                       [(["Arthropoda", "Insecta", "Diptera", "Sarcophagidae", "Ravinia", "Ravinia querula"], None)],
                       params={"truncLen_F": 240, "truncLen_R": 200, "truncLen_cox1_f": 230, "truncLen_cox1_r": 200,
                               "run_lulu": True})
            S = (dco / "asv_table.csv").read_text().splitlines()[0].split(",")[:-1]
            (dco / "read_tracking.csv").write_text(
                "sample,input,filtered,merged,nonchim,length_filtered,codon_filtered,lulu_curated\n"
                + "".join(f"{x},1000,950,900,700,850,690,690\n" for x in S))
            (dco / "summary.json").write_text(json.dumps({"marker": "COX1", "pipeline": "cox1", "lulu": False,
                                                          "lulu_requested": True, "lulu_removed": 0}))
            rc = RB.Run(dco, set(), "3 Targets (COI)")
            check("COI kind and truncLen_cox1 (230 / 200; ceiling 418, not 240+200-12)",
                  rc.kind == "cox1" and rc.trunc == "230 / 200" and rc.ceiling == 418, (rc.kind, rc.trunc, rc.ceiling))
            st = [k for _l, k in rc.stages]
            check("COI steps: length filter before chimera removal; no LULU step when it did not run",
                  st.index("length_filtered") < st.index("nonchim") and "lulu_curated" not in st, st)
            hc = RB.build_profile_html(rc)
            check("COI: LULU not installed said so (method and settings)",
                  "LULU package is not installed on this machine" in hc and "skipped (LULU not installed)" in hc
                  and "were curated with LULU" not in hc)

            # 18S PR2: marker label and off-target table rows that add up
            d18 = make(td / "18s", "18S-nema", names("{s}-{r}_TAR", "{s}-CT_TAR"),
                       ["Domain", "Supergroup", "Phylum", "Subdivision", "Class", "Order", "Family", "Genus", "Species"],
                       [(["Eukaryota", "TSAR", "Alveolata", "Dinoflagellata", "Dinophyceae", "Gymnodiniales",
                          "Gymnodiniaceae", "Gyrodinium", "Gyrodinium_rubrum"], None)],
                       params={"dbPath": "/db/pr2_version_5.1.1_SSU_dada2.fasta.gz"})
            (d18 / "qc_report.json").write_text(json.dumps({"offtarget": {
                "reads": 2260, "pct": 0.19, "organelle_reads": 0, "organelle_pct": 0,
                "unassigned_reads": 497, "unassigned_pct": 0.04}}))
            h18 = RB.build_profile_html(RB.Run(d18, set(), "3 Targets (18S)"))
            check("18S PR2: labelled 18S, not 18S-nema", "18S-nema" not in h18 and "&middot; 18S &middot;" in h18)
            check("18S off-target: bacteria/archaea row (2,260 - 0), no 16S wording",
                  "Outside the target domain (bacteria and archaea)" in h18 and ">2,260<" in h18
                  and "16S primers co-amplify" not in h18)
            fe = (ROOT / "frontend-beta" / "index.html").read_text(encoding="utf-8")
            check("job page: COX1 shows its own truncLen, ITS none", "p.truncLen_cox1_f" in fe
                  and "none (ITS reads are not truncated)" in fe)
            # Download results / customer report: 18S (PR2) and COI packages build
            import amp_package as AP
            (dco / "summary.json").write_text(json.dumps({"marker": "COX1", "pipeline": "cox1", "lulu": True,
                "tax_method": "vsearch", "tax_cutoffs": {"Species": 97, "Genus": 95}, "length_range": "300-330 bp",
                "codon_table": 5, "numt_removed": 19}))
            (dco / "pipeline.log").write_text("LULU package not installed. Run install_extensions.sh first.\n")
            (dco / "run_params.json").write_text(json.dumps({"marker": "COX1", "truncLen_F": 240, "truncLen_R": 200,
                "truncLen_cox1_f": 230, "truncLen_cox1_r": 200,
                "dbPath": "/db/MIDORI2/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz"}))
            for d_, s_, want in ((d18, "KK-NE-1_TAR", ["18S rRNA (eukaryotes)", "PR2 5.1.1"]),
                                 (dco, "KK-NE-1_COI-XT", ["COI amplicon sequencing report", "MIDORI2 GenBank 273 (COI)",
                                                          "trimming the forward reads at 230", "at least 20 bases",
                                                          "Stop-codon filter", "LULU package is not installed",
                                                          "VSEARCH: a versatile"])):
                AP.build_sample(d_, None, s_, td / ("pk_" + d_.name), formats=("html",), tools={})
                hp = next((td / ("pk_" + d_.name)).rglob("Report_*.html")).read_text()
                check(f"Download results / customer report builds for {d_.name}", all(w in hp for w in want),
                      [w for w in want if w not in hp])
            m_ = (ROOT / "backend" / "main.py").read_text(encoding="utf-8")
            check("Download results allowed for 18S and COX1 jobs", '"18S-NEMA", "COX1", "COI")' in m_
                  and '"18S", "18S-NEMA", "COX1", "COI"].indexOf' in fe)
            # BLAST: a MitoFish record with no lineage (Oxyurichthys papuensis, GEOFI763-21) is a fish
            import blast_taxonomy as BT
            L = {"a": {"source": "MIDORI2;MitoFish", "Phylum": "Chordata", "Class": "Actinopteri", "Order": "Gobiiformes",
                       "Family": "Gobiidae", "Genus": "Oxyurichthys", "Species": "Oxyurichthys auchenolepis"},
                 "b": {"source": "MitoFish", "Phylum": "", "Class": "", "Order": "", "Family": "", "Genus": "",
                       "Species": "Oxyurichthys papuensis"},
                 "h": {"source": "MIDORI2", "Phylum": "Chordata", "Class": "Mammalia", "Order": "Primates",
                       "Family": "Hominidae", "Genus": "Homo", "Species": "Homo sapiens"}}
            BT.fill_lineages(L)
            th = {"Species": 98.5, "Genus": 95, "Family": 90, "Order": 85, "Class": 80}
            r = BT.assign([{"sseqid": "b", "pident": "98.0", "bitscore": "180", "qcovhsp": "100"},
                           {"sseqid": "h", "pident": "80", "bitscore": "120", "qcovhsp": "100"}], L, th, "vertebrates")
            check("BLAST: a MitoFish record without a lineage takes its genus's; 98% goby -> fish, genus Oxyurichthys",
                  r["target"] == "yes" and r["Genus"] == "Oxyurichthys" and r["Family"] == "Gobiidae", r)
            fm = (ROOT / "backend" / "fish_mito_ref.py").read_text(encoding="utf-8")
            check("fish_mito_ref fills those lineages when the reference is built", "fill_lineages(_lin)" in fm)
            rs = (ROOT / "backend" / "r_scripts" / "cox1_pipeline.R").read_text(encoding="utf-8")
            check("cox1_pipeline.R: summary says whether LULU ran", 'requireNamespace("lulu", quietly = TRUE)\nsummary_list$lulu_requested' in rs)
    except Exception:
        traceback.print_exc()
        FAILS.append("exception")
    return FAILS


if __name__ == "__main__":
    main()
    sys.exit(1 if FAILS else 0)
