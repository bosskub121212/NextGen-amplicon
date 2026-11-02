#!/usr/bin/env python3
"""Fast checks for fish eDNA amplicons: references, BLAST taxonomy, Excel, customer report.

Run alone (python3 tests/test_amp_fast.py) or from tests/test_fast.py.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))

FAILS = []
TESTS = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILS.append(name)


def skip(name, why):
    print(f"  skip  {name}  ({why})")


def _exe(p: Path, body: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body)
    p.chmod(0o755)


# ── Task 1: reference builder ─────────────────────────────────────────────────
def t_midori_latest(tmp: Path):
    import fish_mito_ref as FR
    def boom(url):
        raise OSError("offline")
    r = FR.latest_midori2(fetch=boom)
    check("MIDORI2 latest: offline -> GB273 fallback", r.get("release") == "GB273"
          and "22685701" in r.get("url", ""), str(r))
    def fake(url):
        return json.dumps({"id": 999, "metadata": {"version": "GB275"},
                           "files": [{"key": "MIDORI2_GB275_DADA2.zip"}]}).encode()
    r = FR.latest_midori2(fetch=fake)
    check("MIDORI2 latest: newest Zenodo version used", r.get("release") == "GB275"
          and r["url"].endswith("/records/999/files/MIDORI2_GB275_DADA2.zip/content"), str(r))


def t_mitofish_probe(tmp: Path):
    import datetime as dt
    import fish_mito_ref as FR
    seen = []
    def head(url):
        seen.append(url)
        return "application/x-gzip" if "/2026/09/" in url or "/2026/08/" in url else "text/html"
    m = FR.latest_mitofish(head=head, today=dt.date(2026, 10, 31))
    check("MitoFish: newest published month (current one not out yet)", m == "2026/09", f"{m} {seen}")
    m2 = FR.latest_mitofish(head=lambda u: "text/html", today=dt.date(2026, 10, 31))
    check("MitoFish: nothing found in 6 months -> None", m2 is None, str(m2))


def t_midori_header(tmp: Path):
    import fish_mito_ref as FR
    h = FR.parse_midori_header("Chordata_7711;Actinopteri_186623;Cypriniformes_7952;Cyprinidae_7953;"
                               "Cyprinus_7961;Cyprinus carpio_7962;")
    check("MIDORI2 header: ranks without taxids", h == {"Phylum": "Chordata", "Class": "Actinopteri",
          "Order": "Cypriniformes", "Family": "Cyprinidae", "Genus": "Cyprinus",
          "Species": "Cyprinus carpio"}, str(h))
    h2 = FR.parse_midori_header("Discosea_555280;Flabellinia_1485085;order_Vannellidae_95227;"
                                "Vannellidae_95227;Paravannella_1443143;Paravannella minima_1443144;")
    check("MIDORI2 header: placeholder rank left empty", h2["Order"] == "" and h2["Family"] == "Vannellidae",
          str(h2))


def _mitofish_fixture(tmp: Path):
    import pandas as pd
    t = tmp / "tables"
    t.mkdir()
    pd.DataFrame({"accession": ["A1", "A2", "A3", "A1"], "gene": ["Cytb", "Cytb", "Cytb", "12S_rRNA"],
                  "start": [3, 1, 1, 1], "end": [10, 6, 8, 4], "strand": [1, -1, 1, 1],
                  "length": [8, 6, 8, 4]}).to_parquet(t / "seq_annotation.parquet")
    pd.DataFrame({"gene": ["Cytb", "Cytb", "Cytb"], "accession": ["A1", "A3", "A3"], "chunk": [0, 0, 1],
                  "judge": ["Homospecific", "Homospecific", "Heterospecific"]}).to_parquet(
        t / "seq_heterospecific.parquet")
    pd.DataFrame({"accession": ["A1", "A2", "A3"], "create_date": ["x"] * 3, "taxon_id": ["9", "9", "9"],
                  "data_source": ["GenBank"] * 3}).to_parquet(t / "seq_taxonid.parquet")
    pd.DataFrame({"taxon_id": ["9"], "species_id": ["9"]}).to_parquet(t / "taxonid_speciesid.parquet")
    lin = [("9", 28, "species", "9"), ("9", 27, "genus", "8"), ("9", 26, "family", "7"),
           ("9", 25, "order", "6"), ("9", 20, "class", "5"), ("9", 5, "phylum", "4")]
    pd.DataFrame(lin, columns=["species_id", "rank_order", "rank_name", "lineage_id"]).to_parquet(
        t / "speciesid_lineageid.parquet")
    pd.DataFrame([("9", "species", "Cyprinus carpio"), ("8", "genus", "Cyprinus"), ("7", "family", "Cyprinidae"),
                  ("6", "order", "Cypriniformes"), ("5", "class", "Actinopteri"), ("4", "phylum", "Chordata")],
                 columns=["lineage_id", "lineage_rank", "lineage_name"]).to_parquet(t / "taxonid_name.parquet")
    fa = tmp / "mitofishdb.fa.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">A1\nAACCGGTTAACC\n>A2\nACGTAAGG\n>A3\nACGTACGTAC\n")
    return fa, t


def t_mitofish_regions(tmp: Path):
    try:
        import pandas  # noqa: F401
        import pyarrow  # noqa: F401
    except ImportError:
        return skip("MitoFish regions", "pandas/pyarrow not installed here")
    import fish_mito_ref as FR
    fa, t = _mitofish_fixture(tmp)
    recs = {a: (s, lin) for a, s, lin in FR.mitofish_regions(fa, t, "Cytb")}
    check("MitoFish: Cytb cut by annotation (1-based, inclusive)", recs.get("A1", ("",))[0] == "CCGGTTAA",
          str(recs))
    check("MitoFish: minus strand reverse-complemented", recs.get("A2", ("",))[0] == "TTACGT", str(recs))
    check("MitoFish: heterospecific record dropped", "A3" not in recs, str(recs))
    check("MitoFish: lineage from the taxonomy tables", recs.get("A1", ("", {}))[1].get("Species") == "Cyprinus carpio"
          and recs["A1"][1].get("Class") == "Actinopteri", str(recs.get("A1")))


def t_build_dedup(tmp: Path):
    try:
        import pandas  # noqa: F401
        import pyarrow  # noqa: F401
    except ImportError:
        return skip("reference build", "pandas/pyarrow not installed here")
    import fish_mito_ref as FR
    fa, t = _mitofish_fixture(tmp)
    mid = tmp / "midori.fa.gz"
    with gzip.open(mid, "wt") as fh:
        fh.write(">Chordata_7711;Actinopteri_1;Cypriniformes_2;Cyprinidae_3;Cyprinus_4;Cyprinus carpio_5;\n"
                 "CCGGTTAA\n"
                 ">Chordata_7711;Mammalia_1;Primates_2;Hominidae_3;Homo_4;Homo sapiens_5;\nGGGGCCCC\n")
    info = FR.build(tmp / "ref", "Cytb", mid, fa, t, makeblastdb=None)
    ids = (tmp / "ref" / "ids.tsv").read_text().splitlines()
    check("build: identical sequence + species from both sources kept once", info.get("n_after_dedup") == 3, str(info))
    check("build: sources joined on the merged record", any("MIDORI2;MitoFish" in l for l in ids), "\n".join(ids))
    check("build: counts recorded", info.get("n_midori2") == 2 and info.get("n_mitofish") == 2
          and info.get("n_heterospecific_removed") == 1, str(info))
    check("build: manifest written", (tmp / "ref" / "fishmito.json").exists())
    mid2 = tmp / "midori2.fa.gz"
    with gzip.open(mid2, "wt") as fh:
        fh.write(">Chordata_1;Mammalia_1;Primates_2;Hominidae_3;Homo_4;Homo sapiens_5;\n" + "ACGTNNNNNN" * 20 + "\n"
                 ">Chordata_1;Mammalia_1;Primates_2;Hominidae_3;Homo_4;Homo sapiens_5;\n" + "ACGTACGTAC" * 20 + "\n")
    info2 = FR.build(tmp / "ref2", "Cytb", mid2, None, None, makeblastdb=None)
    check("build: reference records with > 1 % ambiguous bases left out", info2.get("n_after_dedup") == 1
          and info2.get("n_ambiguous_removed") == 1, str(info2))


def t_midori_blast_format(tmp: Path):
    try:
        import pandas  # noqa: F401
        import pyarrow  # noqa: F401
    except ImportError:
        return skip("MIDORI2 BLAST format", "pandas/pyarrow not installed here")
    import fish_mito_ref as FR
    h = FR.parse_midori_header("MH910097.1.16215.17432###root_1;Eukaryota_2759;Chordata_7711;Actinopteri_1;"
                               "Cypriniformes_2;Cyprinidae_3;Cyprinus_4;Cyprinus_carpio_5")
    check("MIDORI2 BLAST header: lineage after ###", h.get("Species") == "Cyprinus carpio"
          and h.get("Phylum") == "Chordata", str(h))
    check("MIDORI2 BLAST header: accession", FR.midori_accession("MH910097.1.16215.17432###root_1") == "MH910097",
          FR.midori_accession("MH910097.1.16215.17432###root_1"))
    fa, t = _mitofish_fixture(tmp)
    mid = tmp / "m.fasta"
    mid.write_text(">A3.1.1.10###root_1;Eukaryota_2759;Chordata_7711;Actinopteri_1;Cypriniformes_2;Cyprinidae_3;"
                   "Cyprinus_4;Cyprinus_carpio_5\nACGTACGTAC\n"
                   ">Z9.1.1.8###root_1;Eukaryota_2759;Chordata_7711;Mammalia_1;Primates_2;Hominidae_3;Homo_4;"
                   "Homo_sapiens_5\nGGGGCCCC\n")
    info = FR.build(tmp / "ref", "Cytb", mid, fa, t, makeblastdb=None)
    ids = (tmp / "ref" / "ids.tsv").read_text()
    check("MIDORI2 record MitoFish judges heterospecific is left out too", "A3" not in ids
          and info.get("n_midori2_heterospecific_removed") == 1, f"{info} {ids}")
    check("MIDORI2 accession kept in ids.tsv", "Z9" in ids, ids)


def t_download_script(tmp: Path):
    sh = (ROOT / "download_databases.sh").read_text()
    for k in ("fishmito", "fish_mito_ref.py", "ngamp-amp", "MIDORI2_CYTB", "FISHMITO_CYTB", "FISHMITO_12S",
              "_BLAST.zip", "mitofishdb.fa.gz", "tables.tar"):
        check(f"download_databases.sh: {k}", k in sh)
    p = subprocess.run(["bash", "-n", str(ROOT / "download_databases.sh")], capture_output=True, text=True)
    check("download_databases.sh: bash -n", p.returncode == 0, p.stderr[-300:])


# ── Task 2: CytB under the 12S entitlement ─────────────────────────────────────
def t_cytb_wiring(tmp: Path):
    import hashlib
    m = (ROOT / "backend" / "main.py").read_text()
    check("main.py: RunParams mito_gene", 'mito_gene:          str   = "12S"' in m)
    check("main.py: CYTB prefers the MIDORI2 CytB trainset", '"CYTB": ["MIDORI2_CYTB"]' in m)
    check("main.py: DB resolved by gene for 12S runs",
          'resolve_taxonomy_db(\n        "CYTB" if (params.marker or "").upper() == "12S"' in m)
    rq = (ROOT / "backend" / "read_qc.py").read_text()
    check("read_qc: CytB insert length", '("TTCCTAGCCATACAYTAYAC", "GGTGGCKCCTCAGAAGGACATTTGKCCYCA")' in rq)
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    for k in ('label:"CytB fish  (L14912 / H15149)"', 'mito_gene:"CYTB"', '["mito_gene",',
              '/12s|mifish|srrna|cytb/i', 'setVal("f_mito_gene"'):
        check(f"console: {k}", k in ui)
    lic = hashlib.md5((ROOT / "backend" / "license.py").read_bytes()).hexdigest()[:8]
    check("license.py untouched (integrity check)", lic == "c444c850", lic)


TESTS += [t_cytb_wiring]


# ── Task 3: BLAST taxonomy ─────────────────────────────────────────────────────
def _lin(cls, order, fam, gen, sp, phy="Chordata"):
    return {"Phylum": phy, "Class": cls, "Order": order, "Family": fam, "Genus": gen, "Species": sp,
            "source": "MIDORI2", "accession": "X"}


LINS = {"f1": _lin("Actinopteri", "Cypriniformes", "Cyprinidae", "Cyprinus", "Cyprinus carpio"),
        "f2": _lin("Actinopteri", "Cypriniformes", "Cyprinidae", "Cyprinus", "Cyprinus rubrofuscus"),
        "f3": _lin("Actinopteri", "Siluriformes", "Bagridae", "Hemibagrus", "Hemibagrus nemurus"),
        "h1": _lin("Mammalia", "Primates", "Hominidae", "Homo", "Homo sapiens")}


def _hit(s, pid, bits=None):
    return {"qseqid": "q1", "sseqid": s, "pident": pid, "qcovhsp": 100.0, "bitscore": bits or pid * 3}


def t_blast_assign(tmp: Path):
    import blast_taxonomy as BT
    T = BT.DEFAULT_ID
    r = BT.assign([_hit("f1", 99.6)], LINS, T, "fish")
    check("BLAST: 99.6 % -> species", r["Species"] == "Cyprinus carpio" and r["rank_reached"] == "Species", str(r))
    r = BT.assign([_hit("f1", 97.0)], LINS, T, "fish")
    check("BLAST: 97 % -> genus only, never a species below 98.5", r["Species"] == "" and r["Genus"] == "Cyprinus",
          str(r))
    r = BT.assign([_hit("f1", 99.2), _hit("f2", 99.0)], LINS, T, "fish")
    check("BLAST: two species within 0.5 % -> lowest common rank (genus)",
          r["Species"] == "" and r["Genus"] == "Cyprinus" and r["n_ties"] == 2, str(r))
    r = BT.assign([_hit("h1", 100.0), _hit("f3", 82.0)], LINS, T, "fish")
    check("BLAST: human best hit under the fish limit -> Homo sapiens, target no",
          r["Species"] == "Homo sapiens" and r["target"] == "no", str(r))
    r = BT.assign([], LINS, T, "fish")
    check("BLAST: no hit -> Unassigned, kept", r["rank_reached"] == "none" and r["Species"] == ""
          and r["target"] == "", str(r))
    lins = dict(LINS)
    lins["f4"] = {**_lin("Actinopteri", "", "Cichlidae", "Oreochromis", "Oreochromis niloticus"), "Phylum": ""}
    lins["f5"] = _lin("Actinopteri", "Cichliformes", "Cichlidae", "Oreochromis", "Oreochromis aureus")
    r = BT.assign([_hit("f4", 98.9), _hit("f5", 98.8)], lins, T, "fish")
    check("BLAST: a tie missing a higher rank does not erase the common genus",
          r["Genus"] == "Oreochromis" and r["Family"] == "Cichlidae" and r["Species"] == "", str(r))
    r = BT.assign([_hit("f3", 91.0)], LINS, T, "fish")
    check("BLAST: 91 % -> family", r["Family"] == "Bagridae" and r["Genus"] == "" and r["rank_reached"] == "Family",
          str(r))


def t_blast_combine(tmp: Path):
    import blast_taxonomy as BT
    T = BT.DEFAULT_ID
    loc = BT.assign([_hit("f1", 99.6)], LINS, T, "fish")
    c = BT.combine(dict(loc), None, T)
    check("NCBI not reached -> not_checked, local stands", c["agreement"] == "not_checked"
          and c["Species"] == "Cyprinus carpio", str(c))
    c = BT.combine(dict(loc), [{"species": "Cyprinus carpio", "genus": "Cyprinus", "pident": 99.8,
                                "accession": "NC_1"}], T)
    check("NCBI same species -> agree", c["agreement"] == "agree" and c["check"] == "", str(c))
    c = BT.combine(dict(loc), [{"species": "Cyprinus rubrofuscus", "genus": "Cyprinus", "pident": 100.0,
                                "accession": "NC_2"}], T)
    check("NCBI other species, better -> cut to genus, flagged check",
          c["agreement"] == "differ" and c["Species"] == "" and c["Genus"] == "Cyprinus" and c["check"] == "yes"
          and c["ncbi_best"] == "Cyprinus rubrofuscus", str(c))
    un = BT.assign([], LINS, T, "fish")
    c = BT.combine(dict(un), [{"species": "Pangasius larnaudii", "genus": "Pangasius", "pident": 99.5,
                               "accession": "NC_3"}], T)
    check("no local hit, NCBI >= 98.5 -> ncbi_only, named, flagged",
          c["agreement"] == "ncbi_only" and c["Species"] == "Pangasius larnaudii" and c["check"] == "yes", str(c))


def t_blast_cli(tmp: Path):
    import blast_taxonomy as BT
    ref = tmp / "ref"
    ref.mkdir()
    (ref / "fishmito_Cytb.fasta").write_text(">f1\nACGT\n")
    (ref / "ids.tsv").write_text("id\tsource\taccession\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
                                 "f1\tMIDORI2;MitoFish\tAB1\tChordata\tActinopteri\tCypriniformes\tCyprinidae\t"
                                 "Cyprinus\tCyprinus carpio\n")
    (ref / "fishmito.json").write_text(json.dumps({"versions": {"MIDORI2": "GB273", "MitoFish": "2026/09"}}))
    fake = tmp / "blastn"
    _exe(fake, "#!/bin/sh\nfor a in \"$@\"; do [ \"$a\" = -remote ] && { echo 'network down' >&2; exit 2; }; done\n"
               "printf 'q1\\tf1\\t99.70\\t150\\t150\\t270\\t100\\t1e-80\\n'\n")
    fa = tmp / "asv.fasta"
    fa.write_text(">q1\nACGTACGT\n>q2\nTTTTGGGG\n")
    out = tmp / "out"
    os.environ["NGAMP_NCBI_URL"] = "http://127.0.0.1:9/unreachable"
    try:
        rc = BT.main(["--asv_fasta", str(fa), "--db", str(ref / "fishmito_Cytb.fasta"), "--limit", "fish",
                      "--ncbi", "on", "--blastn", str(fake), "--out", str(out)])
    finally:
        os.environ.pop("NGAMP_NCBI_URL", None)
    import csv as _csv
    rows = {r["id"]: r for r in _csv.DictReader(open(out / "taxonomy_blast.csv"))} if (out / "taxonomy_blast.csv").exists() else {}
    check("BLAST CLI: exit 0 although NCBI failed", rc == 0, str(rc))
    check("BLAST CLI: q1 named from the local reference", rows.get("q1", {}).get("Species") == "Cyprinus carpio"
          and rows["q1"].get("source") == "MIDORI2;MitoFish", str(rows.get("q1")))
    check("BLAST CLI: q2 without hit kept as Unassigned", rows.get("q2", {}).get("rank_reached") == "none", str(rows))
    check("BLAST CLI: NCBI failure recorded", rows.get("q1", {}).get("agreement") == "not_checked")
    info = json.loads((out / "blast_info.json").read_text()) if (out / "blast_info.json").exists() else {}
    check("BLAST CLI: info records reference versions and NCBI status",
          info.get("reference", {}).get("versions", {}).get("MIDORI2") == "GB273" and info.get("ncbi") == "failed",
          str(info))


def t_blast_r_hook(tmp: Path):
    r = (ROOT / "backend" / "r_scripts" / "dada2_pipeline.R").read_text()
    for k in ('make_option("--tax_method"', 'make_option("--blast_db"', 'make_option("--blast_limit"',
              'make_option("--blast_ncbi"', "blast_taxonomy.py", "taxonomy_table_dada2.csv"):
        check(f"dada2_pipeline.R: {k}", k in r)
    m = (ROOT / "backend" / "main.py").read_text()
    for k in ('tax_method:         str   = "dada2"', '"--tax_method", params.tax_method', "FISHMITO_CYTB"):
        check(f"main.py: {k}", k in m)


def t_ncbi_url_api(tmp: Path):
    import blast_taxonomy as BT
    calls = []
    js = {"BlastOutput2": [{"report": {"results": {"search": {"query_title": "q1", "query_len": 190, "hits": [
        {"description": [{"accession": "PQ1", "taxid": 64152, "sciname": "Channa striata"}],
         "hsps": [{"identity": 190, "align_len": 190, "bit_score": 351, "query_from": 1, "query_to": 190}]},
        {"description": [{"accession": "PQ2", "taxid": 1, "sciname": "Channa argus"}],
         "hsps": [{"identity": 170, "align_len": 190, "bit_score": 250, "query_from": 1, "query_to": 190}]}]}}}},
        {"report": {"results": {"search": {"query_title": "q2", "query_len": 190, "hits": []}}}}]}

    def http(url, data=None):
        calls.append((url, data))
        if data and b"CMD=Put" in data:
            return b"RID = ABC123\nRTOE = 1\n"
        if "SearchInfo" in url:
            return b"Status=READY\n"
        return json.dumps(js).encode()
    fa = tmp / "a.fasta"
    fa.write_text(">q1\nACGT\n>q2\nTTTT\n")
    r = BT.run_ncbi(fa, "txid7898[ORGN]", http=http, sleep=lambda s: None)
    check("NCBI URL API: hits parsed best first with % identity",
          r and r["q1"][0]["species"] == "Channa striata" and r["q1"][0]["pident"] == 100.0
          and r["q1"][1]["pident"] == round(100 * 170 / 190, 2) and r["q2"] == [], str(r))
    check("NCBI URL API: fish limit sent as ENTREZ_QUERY", any(d and b"ENTREZ_QUERY=txid7898" in d for _, d in calls))
    def down(url, data=None):
        raise OSError("no network")
    check("NCBI URL API: network down -> None", BT.run_ncbi(fa, "", http=down, sleep=lambda s: None) is None)


TESTS += [t_blast_assign, t_blast_combine, t_blast_cli, t_blast_r_hook, t_ncbi_url_api]


# ── Task 4: example tools (FastQC, fastp) ────────────────────────────────────
def _fake_qc_tools(tmp: Path, adapter="pass", trimmed=0):
    fq = tmp / "bin" / "fastqc"
    _exe(fq, "#!/usr/bin/env python3\nimport sys, zipfile, os\n"
             "a = sys.argv[1:]\nif '--version' in a:\n    print('FastQC v0.12.1'); sys.exit(0)\n"
             "o = a[a.index('-o') + 1]\n"
             "for f in [x for x in a if x.endswith('.fastq') or x.endswith('.gz')]:\n"
             "    b = os.path.basename(f).split('.')[0]\n"
             "    with zipfile.ZipFile(os.path.join(o, b + '_fastqc.zip'), 'w') as z:\n"
             f"        z.writestr(b + '_fastqc/fastqc_data.txt', '##FastQC\\t0.12.1\\n>>Adapter Content\\t{adapter}\\n>>END_MODULE\\n')\n")
    fp = tmp / "bin" / "fastp"
    _exe(fp, "#!/usr/bin/env python3\nimport sys, json\na = sys.argv[1:]\n"
             "if '--version' in a:\n    print('fastp 0.23.4', file=sys.stderr); sys.exit(0)\n"
             "j = a[a.index('-j') + 1]\n"
             f"json.dump({{'summary': {{'fastp_version': '0.23.4', 'before_filtering': {{'total_reads': 2000}}}},"
             f" 'adapter_cutting': {{'adapter_trimmed_reads': {trimmed}, 'read1_adapter_sequence': 'AGATCGGAAGAGC'}}}},"
             " open(j, 'w'))\n")
    return str(fq), str(fp)


def t_qc_tools(tmp: Path):
    import qc_tools as Q
    fq, fp = _fake_qc_tools(tmp)
    r1, r2 = tmp / "S1_R1.fastq", tmp / "S1_R2.fastq"
    r1.write_text("@a\nACGT\n+\nIIII\n"); r2.write_text("@a\nACGT\n+\nIIII\n")
    d = Q.fastqc([r1, r2], tmp / "fqc", fq)
    check("FastQC: version and adapter status per file", d.get("version") == "0.12.1"
          and set(d.get("adapter", {}).values()) == {"pass"}, str(d))
    check("FastQC sentence: no adapters", Q.sentence_fastqc(d) == "No adapter sequences in raw data", Q.sentence_fastqc(d))
    p = Q.fastp_check(r1, r2, tmp / "fp", fp)
    check("fastp: version and trimmed reads", p.get("version") == "0.23.4" and p.get("adapter_trimmed_reads") == 0, str(p))
    check("fastp sentence: nothing to trim", "No adapter" in Q.sentence_fastp(p), Q.sentence_fastp(p))
    fq2, fp2 = _fake_qc_tools(tmp / "b", adapter="fail", trimmed=150)
    d2 = Q.fastqc([r1], tmp / "fqc2", fq2)
    check("FastQC sentence: adapters found", "Adapter sequences were detected" in Q.sentence_fastqc(d2),
          Q.sentence_fastqc(d2))
    p2 = Q.fastp_check(r1, r2, tmp / "fp2", fp2)
    check("fastp sentence: reads trimmed with share", "150" in Q.sentence_fastp(p2) and "7.5%" in Q.sentence_fastp(p2),
          Q.sentence_fastp(p2))
    check("FastQC missing -> None, not a crash", Q.fastqc([r1], tmp / "x", None) is None)
    sh = (ROOT / "setup_amplicon_report.sh")
    s = sh.read_text() if sh.exists() else ""
    for k in ("microbiomeutilities", "venneuler", "eulerr", "openpyxl", "ngamp-amp", "krona", "--status"):
        check(f"setup_amplicon_report.sh: {k}", k in s)
    if sh.exists():
        r = subprocess.run(["bash", "-n", str(sh)], capture_output=True, text=True)
        check("setup_amplicon_report.sh: bash -n", r.returncode == 0, r.stderr[-200:])


TESTS += [t_qc_tools]


# ── Task 5: run data + Excel ───────────────────────────────────────────────────
SEQ = ["ACGT" * 40, "TTGG" * 40, "CCAA" * 40, "GATC" * 40]


def make_amp_run(d: Path, samples=("S2-12S", "F4 CytB.1")):
    """A finished 12S/CytB DADA2 run folder in the pipeline's real layout."""
    d.mkdir(parents=True, exist_ok=True)
    S = list(samples)
    counts = [[500, 0], [300, 120], [0, 80], [20, 5]]
    with open(d / "asv_table.csv", "w") as fh:
        fh.write(",".join(f'"{x}"' for x in S + ["sequence"]) + "\n")
        for c, q in zip(counts, SEQ):
            fh.write(",".join(str(x) for x in c[:len(S)]) + f',"{q}"\n')
    ranks = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    tax = [["Chordata", "Actinopteri", "Cypriniformes", "Cyprinidae", "Cyprinus", "Cyprinus carpio"],
           ["Chordata", "Actinopteri", "Siluriformes", "Bagridae", "", ""],
           ["Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "Homo sapiens"],
           ["", "", "", "", "", ""]]
    with open(d / "taxonomy_table.csv", "w") as fh:
        fh.write('"",' + ",".join(f'"{r}"' for r in ranks) + "\n")
        for q, t in zip(SEQ, tax):
            fh.write(f'"{q}",' + ",".join(f'"{x}"' if x else "NA" for x in t) + "\n")
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write('"","input","filtered","denoisedF","denoisedR","merged","nonchim"\n')
        for i, s_ in enumerate(S):
            fh.write(f'"{s_}",{1000 - i},{950 - i},{940 - i},{945 - i},{900 - i},{sum(c[i] for c in counts)}\n')
    (d / "cutadapt").mkdir(exist_ok=True)
    for s_ in S:
        (d / "cutadapt" / f"{s_}.cutadapt.json").write_text(json.dumps({
            "read_counts": {"input": 1100, "output": 1000, "read1_with_adapter": 1080, "read2_with_adapter": 1070},
            "cutadapt_version": "4.9"}))
    (d / "blast").mkdir(exist_ok=True)
    cols = ["id", "sequence"] + ranks + ["best_hit", "pident", "qcov", "source", "accession", "n_ties",
                                         "rank_reached", "target", "ncbi_best", "ncbi_pident", "ncbi_accession",
                                         "agreement", "check"]
    with open(d / "blast" / "taxonomy_blast.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        extra = [["Cyprinus carpio", "100.0", "100", "MIDORI2;MitoFish", "AB1", "1", "Species", "yes",
                  "Cyprinus carpio", "100.0", "NC1", "agree", ""],
                 ["Hemibagrus nemurus", "93.1", "100", "MIDORI2", "AB2", "1", "Family", "yes", "", "", "", "ncbi_no_hit", ""],
                 ["Homo sapiens", "100.0", "100", "MIDORI2", "AB3", "1", "Species", "no", "", "", "", "ncbi_no_hit", ""],
                 ["", "", "", "", "", "0", "none", "", "", "", "", "not_checked", ""]]
        for i, (q, t, e) in enumerate(zip(SEQ, tax, extra), 1):
            fh.write(",".join([f"q{i}", q] + t + e) + "\n")
    (d / "blast" / "blast_info.json").write_text(json.dumps({"limit": "fish", "ncbi": "done", "thresholds":
        {"Species": 98.5, "Genus": 95, "Family": 90, "Order": 85, "Class": 80},
        "reference": {"versions": {"MIDORI2": "GB273", "MitoFish": "2026/09"}}, "blastn": "blastn: 2.16.0+"}))
    (d / "summary.json").write_text(json.dumps({"marker": "12S"}))
    (d / "run_params.json").write_text(json.dumps({"marker": "12S", "mito_gene": "CYTB",
        "primer_f": "TTCCTAGCCATACAYTAYAC", "primer_r": "GGTGGCKCCTCAGAAGGACATTTGKCCYCA", "truncLen_F": 0,
        "truncLen_R": 0, "maxEE_F": 2, "maxEE_R": 5, "tax_method": "blast", "blast_limit": "fish"}))
    return d


def t_amp_data(tmp: Path):
    import amp_data as AD
    r = AD.load_run(make_amp_run(tmp / "run"))
    check("run data: samples in table order", r.samples == ["S2-12S", "F4 CytB.1"], str(r.samples))
    t = r.tracking["S2-12S"]
    check("run data: tracking with raw (cutadapt input) and trimmed (output)",
          t["raw_data"] == 1100 and t["trimmed"] == 1000 and t["Nonchim"] == 820
          and abs(t["perc"] - round(100 * 820 / 1100, 2)) < 1e-9, str(t))
    sv = {a["sv"]: a for a in r.asvs}
    check("run data: SV numbers by total abundance", sv["SV001"]["sequence"] == SEQ[0] and sv["SV002"]["sequence"] == SEQ[1],
          str([(a["sv"], a["total"]) for a in r.asvs]))
    check("run data: BLAST columns joined", sv["SV001"]["ncbi_best"] == "Cyprinus carpio"
          and sv["SV001"]["source"] == "MIDORI2;MitoFish", str(sv["SV001"]))
    lab = {a["sv"]: AD.label(a, "Species") for a in r.asvs}
    check("labels: missing rank -> '<higher> unidentified', none -> Unassigned",
          lab["SV002"] == "Bagridae unidentified" and "Unassigned" in lab.values(), str(lab))
    rt = r.rank_table("Genus", "S2-12S")
    check("rank table: reads per genus for one sample", rt.get("Cyprinus") == 500 and rt.get("Bagridae unidentified") == 300,
          str(rt))
    a = r.alpha["S2-12S"]
    check("alpha: Observed / Chao1 / Shannon / Simpson", a["Observed"] == 3 and a["Chao1"] >= 3
          and a["Shannon"] > 0 and 0 < a["Simpson"] < 1, str(a))
    check("non-target share (human) per sample", abs(r.non_target_share("F4 CytB.1") - 100 * 80 / 205) < 0.01,
          str(r.non_target_share("F4 CytB.1")))


def t_excel(tmp: Path):
    try:
        import openpyxl
    except ImportError:
        return skip("Excel", "openpyxl not installed here")
    import amp_data as AD
    import amplicon_excel as AX
    r = AD.load_run(make_amp_run(tmp / "run"))
    p = tmp / "CT_sample1.xlsx"
    AX.write_sample(p, r, "S2-12S")
    wb = openpyxl.load_workbook(p)
    check("Excel: sheets", wb.sheetnames == ["Read_tracking", "ASV_taxonomy", "Phylum", "Class", "Order", "Family",
                                              "Genus", "Species", "Alpha_diversity"], str(wb.sheetnames))
    hdr = [c.value for c in wb["Read_tracking"][1]]
    check("Excel: read-tracking columns as in the lab's report", hdr == ["SampleID", "Sample_name", "raw_data", "trimmed",
          "filtered", "denoisedF", "denoisedR", "merged", "Nonchim", "perc"], str(hdr))
    asv = [[c.value for c in row] for row in wb["ASV_taxonomy"].iter_rows()]
    check("Excel: ASV sheet holds the sequence and only this sample's ASVs",
          "sequence" in asv[0] and len(asv) == 4 and any(SEQ[0] in (x or "") for x in asv[1]), str(asv[:2]))
    p2 = tmp / "all.xlsx"
    AX.write_all(p2, r)
    wb2 = openpyxl.load_workbook(p2)
    check("Excel all samples: one column per sample in Genus", [c.value for c in wb2["Genus"][1]][:4]
          == ["Genus", "S2-12S", "F4 CytB.1", "S2-12S %"] or "F4 CytB.1" in [c.value for c in wb2["Genus"][1]],
          str([c.value for c in wb2["Genus"][1]]))


TESTS += [t_amp_data, t_excel]


# ── Task 6: customer report ───────────────────────────────────────────────────
HEADINGS = ["Report Summary:", "Clean data", "1. Quality control checks on raw sequence data by FastQC",
            "2. Adapter trimming by fastp", "3. Remove primer sequences by cutadapt", "DADA2 pipeline",
            "1. Pre-processing of reads for amplicon analysis", "1) Filter the forward and reverse reads",
            "2) Learn the Error Rates", "3) Merge paired reads", "4) Chimera Removal", "Data quality control",
            "2. Taxonomic assignment", "Explore microbiome profiles using phyloseq package in R",
            "1. Venn diagram", "2. Bar chart of taxonomic profiles", "At Phylum level", "At Class level",
            "At Order level", "At Family level", "At Genus level", "At Species level", "All File",
            "3. Top 50 ASVs Heatmap", "4. Krona plot", "5. Statistical analysis", "1) Rarefaction analysis",
            "2) Alpha diversity", "Estimate richness"]
COVER = ["Custom Metagenomic Sequencing Report Summary", "Client Name", "Order Number", "Customer Email",
         "Number of Samples to Process", "Sequencing Service Provided", "Download Link:", "Zip file name:",
         "Password zip file:", "Expiration of link download:", "Downloading Your Files",
         "Date and Time of Release of Report:", "Report by:", "Metagenomic Data Analysis Workflow Used",
         "Sources of Metagenomic Reference Sequencing", "Reference or Database Sequencing Used",
         "Data basis analysis:", "Workflow Analysis"]


def _fake_rscript(tmp: Path) -> str:
    """Stand-in for Rscript custom_report_plots.R: writes every figure as a tiny PNG."""
    png = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
    p = tmp / "fakeR" / "Rscript"
    _exe(p, "#!/usr/bin/env python3\nimport sys, csv, base64, os, json\n"
            "a = sys.argv[1:]; jobs = a[a.index('--jobs') + 1]\n"
            f"png = base64.b64decode('{png}')\n"
            "for row in csv.reader(open(jobs), delimiter='\\t'):\n"
            "    ind, out, s = row[:3]\n"
            "    for f in ['Box_plot.png', 'Taxonomic profiles/plot_bar_phylum.png', 'Taxonomic profiles/plot_bar_class.png',"
            " 'Taxonomic profiles/plot_bar_order.png', 'Taxonomic profiles/plot_bar_family.png',"
            " 'Taxonomic profiles/plot_bar_genus.png', 'Taxonomic profiles/plot_bar_species.png',"
            " 'Heatmap/Heat_map_label_ASV_number_top50ASVs.png', 'Rarefaction curves/Rarefaction_curves.png',"
            " 'Krona plot/krona_static.png']:\n"
            "        p = os.path.join(out, f); os.makedirs(os.path.dirname(p), exist_ok=True); open(p, 'wb').write(png)\n"
            "    for f in ['plotQualityProfile.pdf', 'plotQualityProfile_filtered.pdf']:\n"
            "        open(os.path.join(out, f), 'wb').write(b'%PDF-1.4')\n"
            "    json.dump({'phyloseq': '1.50', 'ggplot2': '3.5', 'heatmap_tool': 'microbiomeutilities 1.0',"
            " 'vegan': '2.6', 'dada2': '1.34'}, open(os.path.join(ind, 'plots_info.json'), 'w'))\n")
    return str(p)


def _report(tmp: Path, **kw):
    import custom_report as CR
    run = make_amp_run(tmp / "run", **kw)
    up = tmp / "up"
    up.mkdir()
    for s_ in json.loads(json.dumps(kw.get("samples", ["S2-12S", "F4 CytB.1"]))):
        for m in ("R1", "R2"):
            (up / f"{s_.replace(' ', '_')}_{m}_001.fastq").write_text("@a\nACGT\n+\nIIII\n")
    fq, fp = _fake_qc_tools(tmp / "qc")
    out = tmp / "out"
    info = CR.build(run, up, out, tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                         "ktImportText": None}, formats=("html",), job_label="JOB1")
    return info, out


def _in_order(text: str, items: list[str]) -> tuple[bool, str]:
    pos = 0
    for it in items:
        i = text.find(it, pos)
        if i < 0:
            return False, it
        pos = i + len(it)
    return True, ""


def t_custom_report(tmp: Path):
    import re
    info, out = _report(tmp)
    html = Path(info["html"]).read_text()
    ok, miss = _in_order(html, COVER)
    check("customer report: cover labels and front matter in order", ok, f"missing/out of order: {miss}")
    m = re.search(r"Client Name</td><td[^>]*>([^<]*)<", html)
    check("customer report: customer fields left blank", m is not None and m.group(1).strip() == "", m and m.group(0))
    check("customer report: no app or company name", "NextGen" not in html and "Amplicon Report" not in html)
    blocks = html.split("Report Summary:")[1:]
    check("customer report: one block per sample", len(blocks) == 2, str(len(blocks)))
    for b in blocks:
        ok, miss = _in_order("Report Summary:" + b, HEADINGS)
        check("customer report: every heading in order in a sample block", ok, f"missing/out of order: {miss}")
    caps = set(re.findall(r"&ldquo;(Sample\d+_[^&]*?)&rdquo;", html))
    missing = [c for c in caps if not (out / c.replace(" /", "/")).exists()]
    check("customer report: every caption path exists in the folder", caps and not missing, str(missing[:5]))
    check("customer report: odd sample name sanitised in folder, shown as is",
          (out / "Sample2_F4_CytB.1").is_dir() and "F4 CytB.1" in html, str(sorted(p.name for p in out.iterdir())))
    check("customer report: real settings (no truncation, maxEE 2,5, CytB primers)",
          "maxEE=2,5" in html and "TTCCTAGCCATACAYTAYAC" in html and "no truncation" in html)
    check("customer report: real taxonomy method named (MIDORI2 + MitoFish, NCBI txid 7898)",
          "MIDORI2 GB273" in html and "MitoFish 2026/09" in html and "txid 7898" in html)
    check("customer report: human / non-target share stated", "Homo sapiens" in html and "39.0%" in html)
    check("customer report: Excel per sample", (out / "Sample1_S2-12S" / "Images_and_Outfiles" / "JOB1_sample1.xlsx").exists()
          or not _has_openpyxl())
    check("customer report: outfiles written", (out / "Sample1_S2-12S" / "Images_and_Outfiles" / "Taxonomic profiles"
                                               / "outfiles" / "expr.biom").exists())
    check("customer report: zip made", Path(info.get("zip", "")).exists(), str(info))


def _has_openpyxl():
    try:
        import openpyxl  # noqa: F401
        return True
    except ImportError:
        return False


def t_custom_report_zero(tmp: Path):
    import custom_report as CR
    run = make_amp_run(tmp / "run", samples=("S1", "EMPTY"))
    lines = (run / "asv_table.csv").read_text().splitlines()
    fixed = [lines[0]] + [",".join([l.split(",")[0], "0"] + l.split(",")[2:]) for l in lines[1:]]
    (run / "asv_table.csv").write_text("\n".join(fixed) + "\n")
    fq, fp = _fake_qc_tools(tmp / "qc")
    info = CR.build(run, tmp / "noup", tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                         "ktImportText": None}, formats=("html",), job_label="J")
    html = Path(info["html"]).read_text()
    check("customer report: a sample with no reads still gets its block", html.count("Report Summary:") == 2
          and "No reads" in html, "")


def t_report_wiring(tmp: Path):
    r = (ROOT / "backend" / "r_scripts" / "custom_report_plots.R").read_text() if (ROOT / "backend" / "r_scripts" /
                                                                                     "custom_report_plots.R").exists() else ""
    for k in ("plot_taxa_heatmap", "estimate_richness", "rarecurve", "plotQualityProfile", "psmelt", "pheatmap"):
        check(f"custom_report_plots.R: {k}", k in r)
    d = (ROOT / "backend" / "r_scripts" / "dada2_pipeline.R").read_text()
    check("dada2_pipeline.R: cutadapt --json per sample when supported", "cutadapt.json" in d and "--json" in d)
    m = (ROOT / "backend" / "main.py").read_text()
    check("main.py: POST /custom_report/{job_id}", '@app.post("/custom_report/{job_id}")' in m)
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    check("console: Customer report button", "Customer report" in ui and "/custom_report/" in ui)


TESTS += [t_custom_report, t_custom_report_zero, t_report_wiring]
TESTS += [t_midori_latest, t_mitofish_probe, t_midori_header, t_mitofish_regions, t_build_dedup,
          t_midori_blast_format, t_download_script]



# ── Final-review fixes ─────────────────────────────────────────────────────────
GL = {"Cyprinus": {"Phylum": "Chordata", "Class": "Actinopteri", "Order": "Cypriniformes", "Family": "Cyprinidae",
                   "Genus": "Cyprinus"},
      "Carassius": {"Phylum": "Chordata", "Class": "Actinopteri", "Order": "Cypriniformes", "Family": "Cyprinidae",
                    "Genus": "Carassius"},
      "Pangasius": {"Phylum": "Chordata", "Class": "Actinopteri", "Order": "Siluriformes", "Family": "Pangasiidae",
                    "Genus": "Pangasius"}}


def t_fix_identical_seq_species(tmp: Path):
    import blast_taxonomy as BT
    check("BLAST: target cap >= 5000 so ties are not hidden", BT.MAX_TARGETS >= 5000, str(BT.MAX_TARGETS))
    src = (ROOT / "backend" / "blast_taxonomy.py").read_text()
    check("BLAST: -max_target_seqs uses MAX_TARGETS", "str(MAX_TARGETS)" in src)
    try:
        import pandas  # noqa: F401
    except ImportError:
        return skip("shared sequence build", "pandas not installed here")
    import fish_mito_ref as FR
    mid = tmp / "m.fa.gz"
    seq = "ACGTTGCAAC" * 15
    with gzip.open(mid, "wt") as fh:
        fh.write(">Chordata_1;Actinopteri_2;Cypriniformes_3;Cyprinidae_4;Cyprinus_5;Cyprinus carpio_6;\n" + seq + "\n"
                 ">Chordata_1;Actinopteri_2;Cypriniformes_3;Cyprinidae_4;Cyprinus_5;Cyprinus rubrofuscus_7;\n" + seq
                 + "\n")
    info = FR.build(tmp / "ref", "Cytb", mid, None, None, makeblastdb=None)
    (tmp / "ref" / "fishmito_Cytb.fasta").exists()
    lin = BT.load_ids(tmp / "ref" / "fishmito_Cytb.fasta")
    check("build: one sequence shared by two species kept once", info.get("n_after_dedup") == 1
          and info.get("n_shared_seq_merged") == 1, str(info))
    r0 = lin.get("fm0", {})
    check("build: shared record named to their common genus, species listed",
          r0.get("Genus") == "Cyprinus" and r0.get("Species") == ""
          and "Cyprinus carpio" in r0.get("species_all", "") and "Cyprinus rubrofuscus" in r0.get("species_all", ""),
          str(r0))
    r = BT.assign([_hit("fm0", 100.0)], lin, BT.DEFAULT_ID, "fish")
    check("BLAST: 100 % to a shared sequence -> genus, both species counted",
          r["Genus"] == "Cyprinus" and r["Species"] == "" and r["n_ties"] == 2 and "Cyprinus carpio" in r["best_hit"],
          str(r))


def t_fix_lca_reclassified(tmp: Path):
    import blast_taxonomy as BT
    lins = {"a": _lin("Actinopteri", "Cichliformes", "Cichlidae", "Oreochromis", "Oreochromis niloticus"),
            "b": _lin("Actinopteri", "Perciformes", "Cichlidae", "Oreochromis", "Oreochromis niloticus"),
            "c": _lin("Actinopteri", "Perciformes", "Cichlidae", "Oreochromis", "Oreochromis aureus")}
    r = BT.assign([_hit("a", 100.0), _hit("b", 99.9)], lins, BT.DEFAULT_ID, "fish")
    check("LCA: same species under two order names -> species kept, order of the top hit",
          r["Species"] == "Oreochromis niloticus" and r["Order"] == "Cichliformes", str(r))
    r = BT.assign([_hit("a", 99.2), _hit("c", 99.0)], lins, BT.DEFAULT_ID, "fish")
    check("LCA: tied species of one genus under two order names -> genus and family kept",
          r["Genus"] == "Oreochromis" and r["Family"] == "Cichlidae" and r["Species"] == ""
          and r["Order"] == "Cichliformes", str(r))


def t_fix_combine_rules(tmp: Path):
    import blast_taxonomy as BT
    T = BT.DEFAULT_ID
    loc = BT.assign([_hit("f1", 99.6)], LINS, T, "fish")
    nb = lambda sp, pid, acc="N": {"species": sp, "genus": sp.split(" ")[0], "pident": pid, "accession": acc}
    c = BT.combine(dict(loc), [nb("Cyprinus carpio", 100.0), nb("Cyprinus rubrofuscus", 100.0)], T, GL)
    check("NCBI: two species tied at the top -> species cut, flagged",
          c["agreement"] == "ncbi_tie" and c["Species"] == "" and c["Genus"] == "Cyprinus" and c["check"] == "yes",
          str(c))
    c = BT.combine(dict(loc), [nb("Cyprinus sp.", 100.0), nb("uncultured Cyprinus", 100.0),
                               nb("Cyprinus carpio", 100.0)], T, GL)
    check("NCBI: 'sp.' / uncultured names are not rival species", c["agreement"] == "agree"
          and c["Species"] == "Cyprinus carpio", str(c))
    c = BT.combine(dict(loc), [nb("Cyprinus rubrofuscus", 99.4)], T, GL)
    check("NCBI: other species within 0.5 % below -> differ, not local_better",
          c["agreement"] == "differ" and c["Species"] == "", str(c))
    c = BT.combine(dict(loc), [nb("Cyprinus rubrofuscus", 98.9)], T, GL)
    check("NCBI: other species clearly lower -> local_better", c["agreement"] == "local_better"
          and c["Species"] == "Cyprinus carpio", str(c))
    c = BT.combine(dict(loc), [nb("Carassius auratus", 100.0)], T, GL)
    check("NCBI differ, other genus of the same family -> cut to family",
          c["agreement"] == "differ" and c["Genus"] == "" and c["Family"] == "Cyprinidae", str(c))
    c = BT.combine(dict(loc), [nb("Unknownia dubia", 100.0)], T, GL)
    check("NCBI differ, genus not in the reference -> cut to order",
          c["Family"] == "" and c["Order"] == "Cypriniformes" and c["check"] == "yes", str(c))
    un = BT.assign([], LINS, T, "fish")
    c = BT.combine(dict(un), [nb("uncultured fish", 100.0), nb("Pangasius larnaudii", 99.5)], T, GL)
    check("ncbi_only: uncultured skipped, higher ranks filled from the reference",
          c["Species"] == "Pangasius larnaudii" and c["Family"] == "Pangasiidae" and c["Phylum"] == "Chordata"
          and c["agreement"] == "ncbi_only", str(c))
    c = BT.combine(dict(un), [nb("Oreochromis niloticus x Oreochromis aureus", 100.0)], T, GL)
    check("ncbi_only: a hybrid name is not given as a species", c["Species"] == "" and c["agreement"] == "ncbi_only",
          str(c))
    c = BT.combine(dict(loc), [nb("Cyprinus carpio", 100.0)], T)
    check("combine without the genus map still works", c["agreement"] == "agree", str(c))


def t_fix_ncbi_keeps_batches(tmp: Path):
    import blast_taxonomy as BT
    puts = []
    one = lambda t: {"report": {"results": {"search": {"query_title": t, "query_len": 100, "hits": [
        {"description": [{"accession": "PQ1", "sciname": "Channa striata"}],
         "hsps": [{"identity": 100, "align_len": 100, "bit_score": 180, "query_from": 1, "query_to": 100}]}]}}}}

    def http(url, data=None):
        if data and b"CMD=Put" in data:
            puts.append(data)
            if len(puts) == 2:
                raise OSError("connection reset")
            return b"RID = R1\nRTOE = 1\n"
        if "SearchInfo" in url:
            return b"Status=READY\n"
        return json.dumps({"BlastOutput2": [one(""), one("q1")]}).encode()
    fa = tmp / "a.fasta"
    fa.write_text(">q1\nACGT\n>q2\nTTTT\n")
    st = {}
    os.environ["NGAMP_NCBI_EMAIL"] = "lab@example.org"
    try:
        r = BT.run_ncbi(fa, "", http=http, sleep=lambda s: None, batch=1, status=st)
    finally:
        os.environ.pop("NGAMP_NCBI_EMAIL", None)
    check("NCBI: a failed batch keeps the finished ones", r is not None and "q1" in r and "q2" not in r, str(r))
    check("NCBI: failed batches counted", st.get("failed_batches") == 1, str(st))
    check("NCBI: EMAIL sent when set", any(b"EMAIL=lab%40example.org" in d for d in puts), str(puts[:1]))
    check("NCBI: query missing from a failed batch -> not checked (None), answered with no hit -> []",
          BT._ncbi_hits(r, "q2") is None and BT._ncbi_hits({"q1": []}, "q1") == [] and BT._ncbi_hits(None, "q1") is None)


def t_fix_r_hook_kingdom(tmp: Path):
    import blast_taxonomy as BT
    r = (ROOT / "backend" / "r_scripts" / "dada2_pipeline.R").read_text()
    for k in ('m[, k] <- ifelse(!is.na(m[, "Phylum"]), "Eukaryota"', "length(intersect(cols, RANKS_BT)) < 5",
              "shQuote(opt$blast_limit)", "shQuote(opt$blast_id)"):
        check(f"dada2_pipeline.R: {k}", k in r)
    m = (ROOT / "backend" / "main.py").read_text()
    check("main.py: blast_limit validated", "_BLAST_LIMITS" in m)
    check("main.py: blast_id validated", 're.fullmatch(r"[0-9.]+(,[0-9.]+){4}"' in m)
    i = m.find("def _db_reject_reason")
    check("main.py: fishmito BLAST reference refused as a DADA2 trainset", "fishmito_" in m[i:i + 2500])
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    check("console: fishmito files hidden from the DADA2 picker", "/^fishmito_/i" in ui)
    check("NCBI fish limit covers bony and cartilaginous fishes",
          BT.ENTREZ["fish"] == "txid7898[ORGN] OR txid7777[ORGN]", BT.ENTREZ["fish"])
    check("fish classes: Sarcopterygii (includes tetrapods) not a fish class", "Sarcopterygii" not in BT.FISH_CLASSES)
    rq = (ROOT / "backend" / "read_qc.py").read_text()
    check("read_qc: CytB insert about 240", '"GGTGGCKCCTCAGAAGGACATTTGKCCYCA"): 240' in rq)


def t_fix_dropped_sample(tmp: Path):
    import amp_data as AD
    import custom_report as CR
    run = make_amp_run(tmp / "run", samples=("S1", "GONE"))
    lines = (run / "asv_table.csv").read_text().splitlines()
    (run / "asv_table.csv").write_text("\n".join([",".join([l.split(",")[0]] + l.split(",")[2:]) for l in lines])
                                       + "\n")
    tr = (run / "read_tracking.csv").read_text().splitlines()
    (run / "read_tracking.csv").write_text("\n".join(l for l in tr if not l.startswith('"GONE"')) + "\n")
    rd = AD.load_run(run)
    check("run data: a sample the pipeline dropped is still listed", rd.samples == ["S1", "GONE"], str(rd.samples))
    check("run data: dropped sample has 0 reads, raw from cutadapt", rd.tracking.get("GONE", {}).get("Nonchim") == 0
          and rd.tracking["GONE"]["raw_data"] == 1100, str(rd.tracking.get("GONE")))
    fq, fp = _fake_qc_tools(tmp / "qc")
    info = CR.build(run, tmp / "noup", tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                         "ktImportText": None}, formats=("html",), job_label="J")
    html = Path(info["html"]).read_text()
    check("customer report: dropped sample gets its block with 'No reads'",
          html.count("Report Summary:") == 2 and "No reads" in html and "GONE" in html, "")
    check("customer report: missing raw reads said as such, not 'tool not available'",
          "raw reads of this sample were not found" in html and "not available on this machine" not in html, "")


def t_fix_truthful(tmp: Path):
    import amp_data as AD
    import custom_report as CR
    run = make_amp_run(tmp / "run")
    (run / "blast" / "taxonomy_blast.csv").unlink()
    (run / "blast" / "BLAST_SKIPPED.txt").write_text("BLAST reference not installed")
    t, ref = CR._tax_text(AD.load_run(run))
    check("report: BLAST requested but skipped is said", "BLAST naming was requested but not used" in t, t)
    run2 = make_amp_run(tmp / "run2")
    bi = json.loads((run2 / "blast" / "blast_info.json").read_text())
    bi["reference"] = {"versions": {"MIDORI2": "GB273", "MitoFish": ""}, "n_mitofish": 0}
    bi["agreement"] = {"agree": 3, "differ": 1, "local_better": 2, "ncbi_no_hit": 1, "ncbi_tie": 1, "ncbi_only": 0,
                       "not_checked": 1}
    (run2 / "blast" / "blast_info.json").write_text(json.dumps(bi))
    rp = json.loads((run2 / "run_params.json").read_text())
    rp["pool"] = "FALSE"
    (run2 / "run_params.json").write_text(json.dumps(rp))
    (run2 / "cutadapt" / "S2-12S.cutadapt.json").unlink()
    rd = AD.load_run(run2)
    t2, ref2 = CR._tax_text(rd)
    check("report: no MitoFish claim when the reference has no MitoFish record", "MitoFish" not in t2 + ref2, t2)
    for k in ("better in the local reference", "no NCBI hit", "tied in NCBI", "not checked"):
        check(f"report: NCBI sentence counts '{k}'", k in t2, t2)
    check("report: NCBI fish limit label names both txids", "txid 7898" in t2 and "txid 7777" in t2, t2)
    check("run data: raw count unknown without cutadapt json", rd.tracking["S2-12S"]["raw_data"] is None,
          str(rd.tracking["S2-12S"]))
    fq, fp = _fake_qc_tools(tmp / "qc")
    up = tmp / "up"
    up.mkdir()
    for s_ in ("S2-12S", "F4_CytB.1"):
        for m in ("R1", "R2"):
            (up / f"{s_}_{m}_001.fastq").write_text("@a\nACGT\n+\nIIII\n")
    info = CR.build(run2, up, tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                 "ktImportText": None}, formats=("html",), job_label="J")
    html = Path(info["html"]).read_text()
    check("report: fastp described as a check only", "check only" in html)
    check("report: real pool setting stated", "pool=FALSE" in html
          and "All samples of the run were denoised together" not in html)
    check("report: unknown raw count shown as unknown", "unknown" in html)
    try:
        import openpyxl
        AX = __import__("amplicon_excel")
        AX.write_sample(tmp / "x.xlsx", rd, "S2-12S")
        check("Excel: unknown raw count does not crash", (tmp / "x.xlsx").exists())
    except ImportError:
        skip("Excel unknown raw", "openpyxl not installed here")


TESTS += [t_fix_identical_seq_species, t_fix_lca_reclassified, t_fix_combine_rules, t_fix_ncbi_keeps_batches,
          t_fix_r_hook_kingdom, t_fix_dropped_sample, t_fix_truthful]



def t_fix_no_openpyxl(tmp: Path):
    import re
    saved = {k: sys.modules.get(k) for k in ("openpyxl", "amplicon_excel")}
    sys.modules["openpyxl"] = None
    sys.modules["amplicon_excel"] = None
    try:
        info, out = _report(tmp)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    html = Path(info["html"]).read_text()
    check("report without openpyxl: no caption names an Excel file that was not made",
          ".xlsx" not in html, re.findall(r"[^;>]*\.xlsx", html)[:2])
    check("report without openpyxl: says the Excel file was not made", "Excel file was not made" in html)


TESTS += [t_fix_no_openpyxl]



# ── Per-sample customer reports (v2.24) ───────────────────────────────────────
def make_emu_run(d: Path, samples=("L21", "PGU-16S-full")):
    """A finished ONT-16S Emu run folder in emu_pipeline.py's real layout."""
    d.mkdir(parents=True, exist_ok=True)
    S = list(samples)
    taxa = [("1280", ["Bacteria", "Bacillota", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus",
                      "Staphylococcus aureus"], [700, 0]),
            ("562", ["Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae",
                     "Escherichia", "Escherichia coli"], [200, 4321]),
            ("1386", ["Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus", "Bacillus sp."],
             [0, 5557]),
            ("unassigned", ["", "", "", "", "", "", ""], [33, 21])]
    with open(d / "asv_table.csv", "w") as fh:
        fh.write("," + ",".join(S) + "\n")
        for tid, _, c in taxa:
            fh.write(tid + "," + ",".join(str(x) for x in c[:len(S)]) + "\n")
    with open(d / "taxonomy.csv", "w") as fh:
        fh.write("ASV,Kingdom,Phylum,Class,Order,Family,Genus,Species\n")
        for tid, t, _ in taxa:
            fh.write(tid + "," + ",".join(t) + "\n")
    with open(d / "read_tracking.csv", "w") as fh:
        fh.write("sample,input,qc_filtered,primer_trimmed,final\n")
        for i, s_ in enumerate(S):
            fh.write(f"{s_},{1200 + 4000 * i},{1100 + 4000 * i},{1000 + 4000 * i},{1000 + 4000 * i}\n")
    (d / "summary.json").write_text(json.dumps({"marker": "ONT-16S", "platform": "ONT", "pipeline": "emu",
        "region": "V1-V9", "db_path": "/db/emu_db_mar2026", "emu_min_abundance": 0.0001,
        "emu_counts": "estimated_reads", "emu_chimera_filter": False}))
    (d / "run_params.json").write_text(json.dumps({"marker": "ONT-16S", "primer_f": "AGRGTTYGATYMTGGCTCAG",
        "primer_r": "RGYTACCTTGTTACGACTT", "ont_region": "V1-V9", "ont_min_abundance": 0.0001}))
    return d


def t_emu_data(tmp: Path):
    import amp_data as AD
    r = AD.load_run(make_emu_run(tmp / "run"))
    check("Emu run: kind and samples", r.kind == "emu" and r.samples == ["L21", "PGU-16S-full"], f"{r.kind} {r.samples}")
    t = r.tracking["L21"]
    check("Emu run: tracking raw/qc/primer/final/assigned", t["raw_data"] == 1200 and t["qc_filtered"] == 1100
          and t["primer_trimmed"] == 1000 and t["final"] == 1000 and t["assigned"] == 900
          and abs(t["perc"] - round(100 * 900 / 1200, 2)) < 1e-9, str(t))
    check("Emu run: track columns", r.track_cols == ["raw_data", "qc_filtered", "primer_trimmed", "final", "assigned",
                                                    "perc"], str(r.track_cols))
    lab = {a["id"]: AD.label(a, "Species") for a in r.asvs}
    check("Emu run: 'Genus sp.' is not a species", lab.get("1386") == "Bacillus unidentified", str(lab))
    check("Emu run: unassigned row stays Unassigned", lab.get("unassigned") == "Unassigned", str(lab))
    check("Emu run: no sequences", all(a["sequence"] == "" for a in r.asvs))


def t_only_sample(tmp: Path):
    import amp_data as AD
    r = AD.load_run(make_amp_run(tmp / "run"))
    o = r.only("F4 CytB.1")
    check("only(): one sample", o.samples == ["F4 CytB.1"] and list(o.tracking) == ["F4 CytB.1"]
          and list(o.alpha) == ["F4 CytB.1"], str(o.samples))
    check("only(): only that sample's ASVs, renumbered from its own top",
          [a["sequence"] for a in o.asvs] == [SEQ[1], SEQ[2], SEQ[3]] and o.asvs[0]["sv"] == "SV001"
          and all(list(a["counts"]) == ["F4 CytB.1"] for a in o.asvs), str([(a["sv"], a["total"]) for a in o.asvs]))
    check("only(): BLAST agreement counted on its ASVs only",
          o.blast_info.get("agreement") == {"ncbi_no_hit": 2, "not_checked": 1}, str(o.blast_info.get("agreement")))
    check("only(): parent run untouched", r.samples == ["S2-12S", "F4 CytB.1"] and len(r.asvs) == 4)
    try:
        r.only("nope")
        check("only(): unknown sample refused", False)
    except ValueError:
        check("only(): unknown sample refused", True)
    e = AD.load_run(make_emu_run(tmp / "emu")).only("L21")
    check("only(): Emu run", [a["id"] for a in e.asvs] == ["1280", "562", "unassigned"], str([a["id"] for a in e.asvs]))
    sub = r.subset(["F4 CytB.1", "S2-12S"])
    check("subset(): run order kept", sub.samples == ["S2-12S", "F4 CytB.1"], str(sub.samples))


TESTS += [t_emu_data, t_only_sample]



def _texts(d: Path) -> dict:
    """Every readable text in a package folder: html (data URIs stripped), docx / xlsx XML, text files."""
    import re
    import zipfile
    out = {}
    for f in Path(d).rglob("*"):
        if not f.is_file():
            continue
        try:
            if f.suffix in (".docx", ".xlsx", ".zip"):
                with zipfile.ZipFile(f) as z:
                    out[str(f)] = " ".join(z.read(n).decode("utf-8", "replace") for n in z.namelist()
                                           if n.endswith((".xml", ".rels", ".html", ".txt", ".tsv", ".json", ".biom",
                                                          ".fasta", ".csv")))
            elif f.suffix in (".html", ".txt", ".tsv", ".json", ".biom", ".fasta", ".csv", ".md"):
                out[str(f)] = re.sub(r"data:[a-z/+]+;base64,[A-Za-z0-9+/=]+", "", f.read_text(errors="replace"))
        except Exception as ex:
            out[str(f)] = f"<unreadable {ex}>"
        out[str(f) + "#name"] = str(f.relative_to(d))
    return out


def _leaks(d: Path, needles: list[str]) -> list[str]:
    return [f"{k}: {n}" for k, v in _texts(d).items() for n in needles if n in v]


def _split_report(tmp: Path, run: Path, samples, ups=None, fmts=("html",)):
    import custom_report as CR
    up = tmp / "up"
    up.mkdir(exist_ok=True)
    for name, body in (ups or {}).items():
        (up / name).write_text(body)
    fq, fp = _fake_qc_tools(tmp / "qc")
    return CR.build(run, up, tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                "ktImportText": None, "NanoPlot": None},
                    formats=fmts, job_label="TH-JOB", samples=list(samples), split=True)


def _docx_fmt():
    try:
        import docx  # noqa: F401
        return ("html", "docx")
    except ImportError:
        return ("html",)


def t_split_isolation(tmp: Path):
    import zipfile
    run = make_amp_run(tmp / "run")
    ups = {f"{s_}_{m}_001.fastq": "@a\nACGT\n+\nIIII\n" for s_ in ("S2-12S", "F4_CytB.1") for m in ("R1", "R2")}
    info = _split_report(tmp, run, ["S2-12S", "F4 CytB.1"], ups, _docx_fmt())
    pk = {p["sample"]: p for p in info.get("packages", [])}
    check("split: one package per sample", set(pk) == {"S2-12S", "F4 CytB.1"}, str(info))
    check("split: package zips named by sample", all(Path(p["zip"]).name == f"{n.replace(' ', '_')}_customer_report.zip"
                                                     for n, p in pk.items()), str([p.get("zip") for p in pk.values()]))
    outer = Path(info.get("zip", ""))
    names = zipfile.ZipFile(outer).namelist() if outer.exists() else []
    check("split: outer zip holds exactly the package zips", sorted(names) == sorted(Path(p["zip"]).name
                                                                                 for p in pk.values()), str(names))
    if set(pk) != {"S2-12S", "F4 CytB.1"}:
        return
    a, b = Path(pk["S2-12S"]["dir"]), Path(pk["F4 CytB.1"]["dir"])
    la = _leaks(a, ["F4 CytB.1", "F4_CytB.1", "Homo sapiens", SEQ[2], "TH-JOB"])
    check("split: S2-12S package holds nothing of F4 CytB.1", not la, str(la[:5]))
    lb = _leaks(b, ["S2-12S", "Cyprinus carpio", SEQ[0], "TH-JOB"])
    check("split: F4 CytB.1 package holds nothing of S2-12S", not lb, str(lb[:5]))
    ha = next(a.glob("*_report.html")).read_text()
    check("split: one sample block, cover says 1 sample", ha.count("Report Summary:") == 1
          and "Number of samples in this report</td><td>1<" in ha, "")
    check("split: run-wide workbook not in a package", not list(a.glob("*all_samples*")) and not list(b.glob("*all_samples*")))
    check("split: NCBI counts are this sample's own", "1 agreed at species" in ha and "0 had no NCBI hit" not in ha
          and "1 had no NCBI hit" in ha, "")
    ob = next(b.rglob("expr.abundance.all.txt")).read_text()
    check("split: SV numbers restart in each package", "SV001" in ob and "SV003" in ob and "SV004" not in ob, ob)


def t_report_subset(tmp: Path):
    import custom_report as CR
    run = make_amp_run(tmp / "run")
    fq, fp = _fake_qc_tools(tmp / "qc")
    info = CR.build(run, tmp / "noup", tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                         "ktImportText": None}, formats=("html",), job_label="J",
                    samples=["F4 CytB.1"], split=False)
    html = Path(info["html"]).read_text()
    check("one report for chosen samples: only those", html.count("Report Summary:") == 1 and "S2-12S" not in html, "")
    try:
        CR.build(run, None, tmp / "o2", tools={}, formats=("html",), samples=["nope"], split=True)
        check("unknown sample refused", False)
    except ValueError:
        check("unknown sample refused", True)


EMU_HEADINGS = ["Report Summary:", "Clean data", "1. Quality control checks on raw sequence data",
                "2. Quality and length filtering by Chopper", "3. Remove primer sequences by cutadapt", "Emu pipeline",
                "1) Chimera removal", "2) Taxonomic abundance estimation by Emu", "Data quality control",
                "2. Taxonomic assignment", "Explore microbiome profiles using phyloseq package in R", "1. Venn diagram",
                "2. Bar chart of taxonomic profiles", "At Phylum level", "At Species level", "All File",
                "3. Top 50 taxa Heatmap", "4. Krona plot", "5. Statistical analysis", "1) Rarefaction analysis",
                "2) Alpha diversity", "Estimate richness"]


def t_emu_report(tmp: Path):
    run = make_emu_run(tmp / "run")
    fq = lambda n, L: "".join(f"@r{i}\n{'A' * l}\n+\n{'5' * l}\n" for i, l in zip(range(n), L))
    ups = {"L21.fastq": fq(3, [1500, 1520, 1480]), "PGU-16S-full.fastq": fq(2, [1550, 1560])}
    info = _split_report(tmp, run, ["L21", "PGU-16S-full"], ups)
    pk = {p["sample"]: p for p in info.get("packages", [])}
    if "L21" not in pk:
        return check("Emu report: packages made", False, str(info))
    d = Path(pk["L21"]["dir"])
    html = next(d.glob("*_report.html")).read_text()
    ok, miss = _in_order(html, EMU_HEADINGS)
    check("Emu report: ONT headings in order", ok, f"missing/out of order: {miss}")
    for k in ("Q10", "2000", "emu_db_mar2026", "0.0001", "27F", "AGRGTTYGATYMTGGCTCAG", "estimated read counts",
              "3 reads", "mean length 1,500", "not applied"):
        check(f"Emu report: real setting / result '{k}'", k in html)
    check("Emu report: no DADA2 / FastQC / fastp claim", "DADA2" not in html and "FastQC" not in html
          and "fastp" not in html, "")
    check("Emu report: no ASV sequence file, and says why", "expr.asv.fasta" not in html
          and "no ASV sequences" in html and not list(d.rglob("expr.asv.fasta")), "")
    check("Emu report: tracking columns of the ONT pipeline", "qc_filtered" in html and "primer_trimmed" in html)
    check("Emu report: no non-target paragraph without BLAST", "Non-target" not in html)
    check("Emu report: cover says taxa, not ASVs, for the heatmap", "Top 50 ASVs" not in html and "TOP 50 ASVs" not in html)
    lk = _leaks(d, ["PGU-16S-full", "PGU-16S_full", "Bacillus", "Bacillaceae", "4321", "5557"])
    check("Emu report: L21 package holds nothing of PGU-16S-full", not lk, str(lk[:5]))


def t_16s_report(tmp: Path):
    import custom_report as CR
    run = make_amp_run(tmp / "run")
    rp = json.loads((run / "run_params.json").read_text())
    rp.update({"marker": "16S", "mito_gene": "", "tax_method": "dada2", "dbPath": "/db/silva_nr99_v138.2_toGenus_trainset.fa.gz",
               "primer_f": "CCTACGGGNGGCWGCAG", "primer_r": "GACTACHVGGGTATCTAATCC"})
    (run / "run_params.json").write_text(json.dumps(rp))
    import shutil as _sh
    _sh.rmtree(run / "blast")
    fq, fp = _fake_qc_tools(tmp / "qc")
    info = CR.build(run, tmp / "noup", tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                         "ktImportText": None}, formats=("html",), job_label="J",
                    samples=["S2-12S"], split=True)
    html = next(Path(info["packages"][0]["dir"]).glob("*_report.html")).read_text()
    check("16S report: DADA2 classifier and reference named", "DADA2 naive Bayesian classifier" in html
          and "silva_nr99_v138.2_toGenus_trainset" in html, "")
    check("16S report: primer names 341F / 805R", "341F" in html and "805R" in html)
    check("16S report: no BLAST / MitoFish / non-target text", "MitoFish" not in html and "Non-target" not in html
          and "BLAST naming" not in html, "")


TESTS += [t_split_isolation, t_report_subset, t_emu_report, t_16s_report]



def t_per_sample_wiring(tmp: Path):
    m = (ROOT / "backend" / "main.py").read_text()
    for k in ("class CustomReportRequest(BaseModel):", "samples: Optional[List[str]] = None", "split: bool = True",
              '@app.get("/custom_report/{job_id}/samples")', "_CUSTOM_MARKERS", 'samples=want["samples"]',
              'split=want["split"]'):
        check(f"main.py: {k}", k in m)
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    for k in ("function canCustomReport(", 'id="custSplit"', "data-cs=", '"/samples"', "split:", "samples:"):
        check(f"console: {k}", k in ui)
    check("console: customer report offered for 16S and ONT-16S", '"16S"' in ui[ui.find("function canCustomReport("):
                                                                          ui.find("function canCustomReport(") + 400]
          and '"ONT-16S"' in ui[ui.find("function canCustomReport("):ui.find("function canCustomReport(") + 400])
    em = (ROOT / "backend" / "emu_pipeline.py").read_text()
    for k in ('"qc_min_qual":', '"qc_min_len":', '"qc_max_len":', '"cutadapt_error_rate":', '"cutadapt_overlap":',
              '"qc_skipped":', '"trim_skipped":'):
        check(f"emu_pipeline summary: {k}", k in em)
    sh = (ROOT / "setup_amplicon_report.sh").read_text()
    check("setup_amplicon_report.sh: nanoplot", "nanoplot" in sh.lower())


def t_emu_trim_skipped(tmp: Path):
    run = make_emu_run(tmp / "run")
    sm = json.loads((run / "summary.json").read_text())
    sm.update({"trim_skipped": True, "qc_skipped": True})
    (run / "summary.json").write_text(json.dumps(sm))
    info = _split_report(tmp, run, ["L21"], {})
    html = next(Path(info["packages"][0]["dir"]).glob("*_report.html")).read_text()
    check("Emu report: cutadapt not run is said", "cutadapt was not available" in html and "were removed;" not in html, "")
    check("Emu report: Chopper not run is said", "Chopper was not available" in html)
    check("Emu report: raw reads not found is said", "raw reads of this sample were not found" in html)


TESTS += [t_per_sample_wiring, t_emu_trim_skipped]



# ── v2.24 review fixes ─────────────────────────────────────────────────────────
def t_fix24_name_collision(tmp: Path):
    import zipfile
    run = make_amp_run(tmp / "run", samples=("A 1", "A_1"))
    ups = {"A 1_R1_001.fastq": "@a\nAAAA\n+\nIIII\n", "A 1_R2_001.fastq": "@a\nAAAA\n+\nIIII\n",
           "A_1_R1_001.fastq": "@b\nCCCC\n+\nIIII\n", "A_1_R2_001.fastq": "@b\nCCCC\n+\nIIII\n"}
    info = _split_report(tmp, run, ["A 1", "A_1", "A_1"], ups)
    pk = info.get("packages", [])
    zips = [Path(p["zip"]).name for p in pk]
    check("collision: duplicate request entries made once, one package per sample", len(pk) == 2, str(zips))
    check("collision: package names unique", len(set(zips)) == 2, str(zips))
    outer = zipfile.ZipFile(info["zip"]).namelist()
    check("collision: outer zip entries unique", len(outer) == len(set(outer)) == 2, str(outer))
    import custom_report as CR
    import amp_data as AD
    r = AD.load_run(run)
    up = tmp / "up"
    a = CR._find_reads(up, "A 1", r.samples)[0]
    b = CR._find_reads(up, "A_1", r.samples)[0]
    check("collision: each sample finds its own reads", a is not None and b is not None and a != b
          and a.name.startswith("A 1") and b.name.startswith("A_1"), f"{a} {b}")


def t_fix24_emu_unassigned(tmp: Path):
    import amp_data as AD
    r = AD.load_run(make_emu_run(tmp / "run"))
    t = r.tracking["L21"]
    check("Emu: unassigned row not counted as assigned", t["assigned"] == 900 and t["Nonchim"] == 900, str(t))
    check("Emu: unassigned row not a taxon in alpha", r.alpha["L21"]["Observed"] == 2, str(r.alpha["L21"]))


def t_fix24_old_summary(tmp: Path):
    import amp_data as AD
    import custom_report as CR
    run = make_emu_run(tmp / "run")
    sm = json.loads((run / "summary.json").read_text())
    sm.pop("emu_counts")
    sm["emu_chimera_filter"] = True
    (run / "summary.json").write_text(json.dumps(sm))
    tr = "sample,input,qc_filtered,primer_trimmed,final\nL21,1200,1200,1200,1200\nPGU-16S-full,5200,5100,5000,5000\n"
    (run / "read_tracking.csv").write_text(tr)
    r = AD.load_run(run)
    st = CR.emu_settings(r, "L21")
    check("old summary: no Chopper output and no read lost -> Chopper not run", st["qc_skipped"] is True, str(st))
    check("old summary: no cutadapt output and no read lost -> primers not removed", st["trim_skipped"] is True, str(st))
    check("old summary: chimera filter asked but no evidence -> not claimed", st["chimera_ran"] is False, str(st))
    st2 = CR.emu_settings(r, "PGU-16S-full")
    check("old summary: reads lost at each step -> steps ran", not st2["qc_skipped"] and not st2["trim_skipped"], str(st2))
    info = _split_report(tmp, run, ["L21"], {})
    html = next(Path(info["packages"][0]["dir"]).glob("*_report.html")).read_text()
    check("old summary: chimera asked but not applied is said", "requested but not applied" in html
          and "VSEARCH" not in html and "vsearch --uchime_ref removed" not in html, "")
    run2 = make_emu_run(tmp / "run2")
    with open(run2 / "asv_table.csv", "w") as fh:
        fh.write(",L21,PGU-16S-full\n1280,600000,0\n562,400000,1000000\n")
    sm2 = json.loads((run2 / "summary.json").read_text())
    sm2.pop("emu_counts")
    (run2 / "summary.json").write_text(json.dumps(sm2))
    check("old summary: columns summing to 1e6 read as ppm", AD.load_run(run2).counts_are_reads is False)


def t_fix24_gate_and_species(tmp: Path):
    m = (ROOT / "backend" / "main.py").read_text()
    i = m.find("def _custom_check")
    check("main.py: 16S / 12S customer report only for Illumina DADA2 runs", "sequencerType" in m[i:i + 2500])
    check("main.py: a different request while building -> 409", "status_code=409" in m[m.find("def start_custom_report"):
                                                                                     m.find("def start_custom_report") + 3000])
    import custom_report as CR
    run = make_amp_run(tmp / "run")
    rp = json.loads((run / "run_params.json").read_text())
    rp.update({"marker": "16S", "mito_gene": "", "tax_method": "dada2", "dbPath": "/db/silva_nr99_v138.2_toGenus_trainset.fa.gz"})
    (run / "run_params.json").write_text(json.dumps(rp))
    import shutil as _sh
    _sh.rmtree(run / "blast")
    (run / "pipeline.log").write_text("  Species added from: silva_v138.2_assignSpecies.fa.gz \n")
    fq, fp = _fake_qc_tools(tmp / "qc")
    info = CR.build(run, None, tmp / "out", tools={"fastqc": fq, "fastp": fp, "Rscript": _fake_rscript(tmp),
                                                 "ktImportText": None}, formats=("html",), samples=["S2-12S"], split=True)
    html = next(Path(info["packages"][0]["dir"]).glob("*_report.html")).read_text()
    check("16S: addSpecies and its file stated", "addSpecies" in html and "silva_v138.2_assignSpecies.fa.gz" in html)
    check("16S: MIDORI2 reference not listed", "MIDORI2" not in html)
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    check("console: sample list can be retried", 'id="custRetry"' in ui)


TESTS += [t_fix24_name_collision, t_fix24_emu_unassigned, t_fix24_old_summary, t_fix24_gate_and_species]



def t_fix24_read_quality(tmp: Path):
    import qc_tools as Q
    f = tmp / "q.fastq"
    f.write_text("@a\nACGT\n+\n++??\n@b\nACGT\n+\n++??\n")       # Q10, Q10, Q30, Q30
    d = Q.read_stats([f])
    check("read stats: quality from the mean error rate (as MinKNOW / NanoPlot), not the mean Phred",
          d and abs(d["mean_q"] - 12.97) < 0.05, str(d))
    check("read stats: sentence says how quality is computed", "mean error rate" in Q.sentence_read_stats(d),
          Q.sentence_read_stats(d))


TESTS += [t_fix24_read_quality]


def main() -> int:
    print("fish eDNA amplicons (CytB / 12S)")
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
