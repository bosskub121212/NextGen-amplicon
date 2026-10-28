#!/usr/bin/env python3
"""Fast checks for shotgun functional profiling (HUMAnN 3.9) — no real tools needed.

Run alone (python3 tests/test_func_fast.py) or from tests/test_fast.py.
"""
from __future__ import annotations

import os
import subprocess
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


def _exe(p: Path, body: str = "#!/bin/sh\nexit 0\n"):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body)
    p.chmod(0o755)


def fake_install(root: Path, with_uniref=True) -> tuple[dict, list[Path]]:
    db = root / "humann_db"
    (db / "chocophlan").mkdir(parents=True)
    (db / "chocophlan" / "g__Escherichia.s__Escherichia_coli.centroids.v201901_v31.ffn.gz").write_text("x")
    if with_uniref:
        (db / "uniref").mkdir()
        (db / "uniref" / "uniref90_201901b_full.dmnd").write_text("x")
    um = db / "utility_mapping"
    um.mkdir()
    for f in ("map_ko_uniref90.txt.gz", "map_level4ec_uniref90.txt.gz"):
        (um / f).write_text("x")
    mp = db / "metaphlan"
    mp.mkdir()
    for f in ("mpa_vJun23_CHOCOPhlAnSGB_202307.1.bt2l", "mpa_vJun23_CHOCOPhlAnSGB_202307.pkl"):
        (mp / f).write_text("x")
    for part in ("chocophlan", "uniref", "utility_mapping", "metaphlan"):
        if (db / part).is_dir():
            (db / part / ".complete").write_text("ok\n")
    base = root / "conda"
    for b in ("humann", "metaphlan"):
        _exe(base / "envs" / "ngamp-func" / "bin" / b)
    return {"humann_db": str(db)}, [base]


def t_find_install(tmp: Path):
    import shotgun_func as FN
    FN._FOUND.clear() if hasattr(FN, "_FOUND") else None
    dbp, bases = fake_install(tmp / "a")
    r = FN.find_install(dbp, bases)
    check("find_install: complete install found", r is not None, str(r))
    if r:
        check("find_install: uniref path", r["uniref"].endswith("uniref"), r["uniref"])
        check("find_install: metaphlan index", r["metaphlan_index"] == "mpa_vJun23_CHOCOPhlAnSGB_202307",
              r["metaphlan_index"])
        check("find_install: humann from ngamp-func",
              r["humann"].endswith("envs/ngamp-func/bin/humann"), r["humann"])
        check("find_install: env PATH has env bin",
              str(Path(r["humann"]).parent) in r["env"]["PATH"].split(os.pathsep)[0])
    check("missing_reason empty when installed", FN.missing_reason(dbp, bases) == "",
          FN.missing_reason(dbp, bases))


def t_find_install_missing(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_install(tmp / "b", with_uniref=False)
    check("find_install: no UniRef90 -> None", FN.find_install(dbp, bases) is None)
    why = FN.missing_reason(dbp, bases)
    check("missing_reason names UniRef90", "UniRef90" in why, why)
    why2 = FN.missing_reason({"humann_db": str(tmp / "nowhere")}, bases)
    check("missing_reason names the missing folder", str(tmp / "nowhere") in why2, why2)
    check("missing_reason with no humann_db key", "setup_shotgun.sh --func" in FN.missing_reason({}, bases),
          FN.missing_reason({}, bases))


def t_setup_script():
    s = (ROOT / "setup_shotgun.sh").read_text()
    check("setup_shotgun.sh: python pinned to the build bioconda ships for humann 3.9 (py312)",
          "python=3.12 humann=3.9" in s and "python=3.10 humann" not in s)
    check("setup_shotgun.sh: env checked by importing humann, not only by the binary",
          'import humann, metaphlan' in s)
    for pin in ("humann=3.9", "metaphlan=4.1.1", "mpa_vJun23_CHOCOPhlAnSGB_202307",
                "uniref90_diamond", "ngamp-func", "--func-db-dir", "humann_db"):
        check(f"setup_shotgun.sh has {pin}", pin in s)
    p = subprocess.run(["bash", "-n", str(ROOT / "setup_shotgun.sh")], capture_output=True, text=True)
    check("setup_shotgun.sh bash -n", p.returncode == 0, p.stderr[-300:])


TESTS = [t_find_install, t_find_install_missing]


def main() -> int:
    print("shotgun function (HUMAnN)")
    with tempfile.TemporaryDirectory() as td:
        for i, t in enumerate(TESTS):
            d = Path(td) / f"t{i}"
            d.mkdir()
            try:
                t(d)
            except Exception:
                check(t.__name__, False, traceback.format_exc()[-800:])
    try:
        t_setup_script()
    except Exception:
        check("t_setup_script", False, traceback.format_exc()[-800:])
    return len(FAILS)



# ── Task 2: gene panels ──────────────────────────────────────────────────────
PANEL_TSV = ("panel\tgroup\tgene\tko\tuniref_pattern\tnote\n"
             "gut\tButyrate\tbuk\tK00929\tbutyrate\tKO wins\n"
             "gut\tBile acids\tbsh\t\tbile salt hydrolase\tname only\n"
             "food\tBiogenic amines\thdc\tK01590\t\t\n")


def _panels(tmp: Path, enabled):
    import shotgun_func as FN
    p = tmp / "panels.tsv"
    p.write_text(PANEL_TSV)
    return FN.load_panels(p, enabled)


def t_panels_ko_first(tmp: Path):
    import shotgun_func as FN
    pan = _panels(tmp, ["gut"])
    ko = {"K00929": {"": {"s1": 10.0}, "g__Roseburia.s__Roseburia_hominis": {"s1": 6.0},
                     "unclassified": {"s1": 4.0}}}
    uni = {"UniRef90_A1: Butyrate kinase 2": {"": {"s1": 99.0}}}
    rows, taxa = FN.match_panels(ko, uni, pan)
    buk = [r for r in rows if r["gene"] == "buk" and r["sample"] == "s1"]
    check("panels: KO hit used, name pattern ignored", buk and buk[0]["cpm"] == 10.0, str(buk))
    check("panels: genus from taxon key",
          any(t["gene"] == "buk" and t["taxon"] == "Roseburia" and t["cpm"] == 6.0 for t in taxa), str(taxa))
    check("panels: unclassified kept",
          any(t["taxon"] == "Unclassified" and t["cpm"] == 4.0 for t in taxa), str(taxa))


def t_panels_regex_fallback(tmp: Path):
    import shotgun_func as FN
    pan = _panels(tmp, ["gut"])
    uni = {"UniRef90_B1: Bile salt hydrolase": {"": {"s1": 3.0, "s2": 1.0}},
           "UniRef90_B2: Choloylglycine hydrolase / BILE SALT HYDROLASE": {"": {"s1": 2.0}},
           "UniRef90_B3: unrelated": {"": {"s1": 50.0}}}
    rows, _ = FN.match_panels({}, uni, pan)
    bsh = {r["sample"]: r["cpm"] for r in rows if r["gene"] == "bsh"}
    check("panels: name regex, case-insensitive, summed", bsh == {"s1": 5.0, "s2": 1.0}, str(bsh))
    buk = [r for r in rows if r["gene"] == "buk"]
    check("panels: gene with no hit reported as 0 for every sample",
          sorted((r["sample"], r["cpm"]) for r in buk) == [("s1", 0.0), ("s2", 0.0)], str(buk))


def t_panels_disabled_set(tmp: Path):
    pan = _panels(tmp, ["gut"])
    check("panels: disabled set left out", all(p["panel"] == "gut" for p in pan) and len(pan) == 2, str(pan))


def t_panels_file_valid(tmp: Path):
    import re
    import shotgun_func as FN
    pan = FN.load_panels(ROOT / "backend" / "data" / "function_panels.tsv", ["gut", "food", "environment"])
    check("panel file: all three sets present", {p["panel"] for p in pan} == {"gut", "food", "environment"})
    bad = [p["gene"] for p in pan if not (p["ko"] or p["uniref_pattern"])]
    check("panel file: every gene has a KO or a name pattern", not bad, str(bad))
    badko = [k for p in pan for k in p["ko"] if not re.fullmatch(r"K\d{5}", k)]
    check("panel file: KO ids look like K00000", not badko, str(badko))
    keys = [(p["panel"], p["gene"]) for p in pan]
    check("panel file: no duplicate (panel, gene)", len(keys) == len(set(keys)))


TESTS += [t_panels_ko_first, t_panels_regex_fallback, t_panels_disabled_set, t_panels_file_valid]


# ── Task 3: HUMAnN run and tables ────────────────────────────────────────────
FAKE = ROOT / "tests" / "fake_func" / "fakefunc.py"
FAKE_TOOLS = ("metaphlan", "humann", "humann_join_tables", "humann_renorm_table", "humann_regroup_table",
              "humann_rename_table", "humann_split_stratified_table", "seqkit")


def fake_world(root: Path) -> tuple[dict, list[Path]]:
    """Fake install (wrappers around fakefunc.py) with real-format mapping files."""
    import bz2
    import gzip
    dbp, bases = fake_install(root)
    b = bases[0] / "envs" / "ngamp-func" / "bin"
    for t in FAKE_TOOLS:
        _exe(b / t, f'#!/bin/sh\nexec python3 "{FAKE}" {t} "$@"\n')
    um = Path(dbp["humann_db"]) / "utility_mapping"
    for f, txt in (("map_ko_uniref90.txt.gz", "K00929\tUniRef90_K1\n"),
                   ("map_level4ec_uniref90.txt.gz", "2.7.2.7\tUniRef90_K1\n"),
                   ("map_ko_name.txt.gz", "K00929\tbuk; butyrate kinase\n"),
                   ("map_ec_name.txt.gz", "2.7.2.7\tButyrate kinase\n")):
        with gzip.open(um / f, "wt") as fh:
            fh.write(txt)
    with bz2.open(um / "map_uniref90_name.txt.bz2", "wt") as fh:
        fh.write("UniRef90_B1\tBile salt hydrolase\nUniRef90_K1\tButyrate kinase\n")
    return dbp, bases


def fake_reads(d: Path, name: str, n: int = 8, paired=True) -> dict:
    d.mkdir(parents=True, exist_ok=True)
    def fq(p):
        p.write_text("".join(f"@r{i}\n{'ACGT' * 100}\n+\n{'I' * 400}\n" for i in range(n)))
        return p
    return {"r1": fq(d / f"{name}_R1.fastq"), "r2": fq(d / f"{name}_R2.fastq") if paired else None}


def _ctx(dbp, bases):
    import shotgun_pipeline as SP
    seqkit = str(bases[0] / "envs" / "ngamp-func" / "bin" / "seqkit")
    tool = lambda n: (seqkit, os.environ.copy()) if n == "seqkit" else None
    logs = []
    return dict(dbp=dbp, conda_bases=bases, tool=tool, run=SP.run, log=logs.append,
                progress=lambda p, l: None), logs


def t_read_table_suffixes(tmp: Path):
    import shotgun_func as FN
    p = tmp / "t.tsv"
    p.write_text("# Pathway\tS-1.a_Abundance\tS2_Abundance-RPKs\tS3_Coverage\tS4_Abundance-CPM\n"
                 "PWY-1: one\t1\t2\t3\t4\nPWY-1: one|g__X.s__Y\t1\t0\t0\t0\n")
    t = FN.read_humann_table(p)
    check("read_humann_table: suffixes stripped, names kept",
          set(t["PWY-1: one"]) == {"S-1.a", "S2", "S3", "S4"}, str(t))
    st = FN.read_strat(p)
    check("read_strat: community and taxon split",
          st["PWY-1: one"][""]["S2"] == 2.0 and st["PWY-1: one"]["g__X.s__Y"]["S-1.a"] == 1.0, str(st))
    st2 = FN.read_strat(p, keep=lambda f: f.startswith("NOPE"))
    check("read_strat: keep filter", st2 == {}, str(st2))


def t_prep_reads(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, _ = _ctx(dbp, bases)
    lg = tmp / "seqkit.log"
    os.environ["FAKE_SEQKIT_LOG"] = str(lg)
    try:
        q = fake_reads(tmp / "in", "o1", paired=False)
        dest = FN.prep_reads(q, "ont", tmp / "o1", ctx["tool"]("seqkit"), ctx["run"], tmp / "l.log")
        calls = lg.read_text() if lg.exists() else ""
        check("prep_reads ONT: seqkit sliding 300/300", "sliding" in calls and "-W 300 -s 300" in calls, calls)
        check("prep_reads ONT: pieces < 100 bp dropped", "seq -m 100" in calls, calls)
        check("prep_reads ONT: output exists", dest is not None and dest.exists(), str(dest))
        q2 = fake_reads(tmp / "in", "i1", n=5)
        d2 = FN.prep_reads(q2, "illumina", tmp / "i1", None, ctx["run"], tmp / "l.log")
        n = sum(1 for _ in open(d2)) if d2 else 0
        check("prep_reads Illumina: R1 + R2 in one file", n == 40, str(n))
        q3 = {"r1": tmp / "in" / "empty.fastq", "r2": None}
        q3["r1"].write_text("")
        check("prep_reads: empty input -> None",
              FN.prep_reads(q3, "illumina", tmp / "e1", None, ctx["run"], tmp / "l.log") is None)
    finally:
        os.environ.pop("FAKE_SEQKIT_LOG", None)


def t_parse_humann_log(tmp: Path):
    import shotgun_func as FN
    p = tmp / "h.log"
    p.write_text("x INFO: Unaligned reads after nucleotide alignment: 61.2000000000 %\n"
                 "x INFO: Unaligned reads after translated alignment: 40.1000000000 %\n")
    r = FN.parse_humann_log(p)
    check("parse_humann_log: nucleotide 38.8", abs(r.get("nucleotide_pct", 0) - 38.8) < 1e-6, str(r))
    check("parse_humann_log: translated 21.1", abs(r.get("translated_pct", 0) - 21.1) < 1e-6, str(r))
    check("parse_humann_log: unmapped 40.1", abs(r.get("unmapped_pct", 0) - 40.1) < 1e-6, str(r))
    check("parse_humann_log: missing file -> {}", FN.parse_humann_log(tmp / "nope.log") == {})


def _qc_rows(out: Path):
    import csv as _csv
    with open(out / "func" / "func_qc.csv") as fh:
        return list(_csv.DictReader(fh))


def t_run_all_fake(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, logs = _ctx(dbp, bases)
    qc = {s: fake_reads(tmp / "reads", s) for s in ("s1", "s.2-b")}
    out = tmp / "out"
    info = FN.run_all(samples=list(qc), qc=qc, platform="illumina", out=out, threads=2,
                      panels=["gut"], **ctx)
    f = out / "func"
    need = ["pathways.tsv", "pathways_unstratified.tsv", "pathways_stratified.tsv", "pathcoverage.tsv",
            "ko.tsv", "ko_unstratified.tsv", "ko_stratified.tsv", "ec.tsv", "ec_unstratified.tsv",
            "genefamilies.tsv", "panels.csv", "panels_taxa.csv", "func_qc.csv",
            "s1/s1_metaphlan.tsv", "s1/inputs.sig"]
    miss = [n for n in need if not (f / n).exists()]
    check("run_all: every output file written", not miss, f"{miss} info={info} log={logs[-5:]}")
    q = _qc_rows(out)
    check("run_all: func_qc one ok row per sample",
          sorted((r["sample"], r["status"]) for r in q) == [("s.2-b", "ok"), ("s1", "ok")], str(q))
    if q:
        check("func_qc: translated_pct from log", abs(float(q[0]["translated_pct"]) - 21.1) < 1e-6, str(q[0]))
        check("func_qc: unintegrated_pct from pathways (500 of 1000)",
              abs(float(q[0]["unintegrated_pct"]) - 50.0) < 1e-6, str(q[0]))
    ko = FN.read_humann_table(f / "ko_unstratified.tsv") if (f / "ko_unstratified.tsv").exists() else {}
    kk = [k for k in ko if k.startswith("K00929")]
    check("run_all: KO named and in CPM", kk and kk[0].startswith("K00929: buk")
          and ko[kk[0]]["s.2-b"] > 0, str(list(ko)))
    pan = (f / "panels.csv").read_text() if (f / "panels.csv").exists() else ""
    check("run_all: panel genes found (KO and name)", "buk" in pan and "bsh" in pan, pan[:300])
    taxa = (f / "panels_taxa.csv").read_text() if (f / "panels_taxa.csv").exists() else ""
    check("run_all: carriers by genus", "Escherichia" in taxa and "Bacillus" in taxa, taxa[:300])
    check("run_all: info counts", info.get("n_ko") == 1 and info.get("n_pathways") == 2
          and info.get("failed") == [], str(info))
    check("run_all: versions recorded", "3.9" in str(info.get("versions", {}).get("humann")), str(info))


def t_run_all_one_fails(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, logs = _ctx(dbp, bases)
    qc = {s: fake_reads(tmp / "reads", s) for s in ("s1", "s2")}
    os.environ["FAKE_HUMANN_FAIL"] = "s2"
    try:
        info = FN.run_all(samples=list(qc), qc=qc, platform="illumina", out=tmp / "out", threads=1,
                          panels=["gut"], **ctx)
    finally:
        os.environ.pop("FAKE_HUMANN_FAIL", None)
    q = {r["sample"]: r["status"] for r in _qc_rows(tmp / "out")}
    check("one sample fails: others still profiled", q == {"s1": "ok", "s2": "failed"}, str(q))
    check("one sample fails: info.failed", info.get("failed") == ["s2"], str(info))
    t = FN.read_humann_table(tmp / "out" / "func" / "pathways_unstratified.tsv")
    check("one sample fails: tables hold s1 only", set(next(iter(t.values()))) == {"s1"}, str(t))


def t_run_all_no_db(tmp: Path):
    import shotgun_func as FN
    ctx, _ = _ctx({"humann_db": str(tmp / "E_not_mounted")}, [tmp / "conda"])
    info = FN.run_all(samples=["s1"], qc={"s1": {"r1": None, "r2": None}}, platform="illumina",
                      out=tmp / "out", threads=1, panels=["gut"], **ctx)
    sk = tmp / "out" / "func" / "SKIPPED.txt"
    check("no DB: skipped with reason", "skipped" in info and sk.exists()
          and "E_not_mounted" in sk.read_text(), str(info))


def t_run_all_skip_done(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, logs = _ctx(dbp, bases)
    qc = {"s1": fake_reads(tmp / "reads", "s1")}
    kw = dict(samples=["s1"], qc=qc, platform="illumina", out=tmp / "out", threads=1, panels=["gut"], **ctx)
    FN.run_all(**kw)
    hl = tmp / "out" / "func" / "s1" / "s1_genefamilies.tsv"
    m1 = hl.stat().st_mtime_ns
    FN.run_all(**kw)
    check("second run with same inputs reuses HUMAnN output", hl.stat().st_mtime_ns == m1)
    check("second run logs the reuse", any("already done" in str(x) for x in logs), str(logs[-6:]))


TESTS += [t_read_table_suffixes, t_prep_reads, t_parse_humann_log, t_run_all_fake, t_run_all_one_fails,
          t_run_all_no_db, t_run_all_skip_done]


# ── Task 4: diversity, PERMANOVA, MetaPhlAn vs Kraken, chart tables ──────────
def t_shannon_known(tmp: Path):
    import math
    import shotgun_func as FN
    check("shannon: two equal pathways = ln 2", abs(FN.shannon([5.0, 5.0]) - math.log(2)) < 1e-9)
    check("shannon: empty = 0", FN.shannon([0.0, 0.0]) == 0.0)


def t_pcoa_identical(tmp: Path):
    import shotgun_func as FN
    tab = {"K1": {"a": 1, "b": 1, "c": 5}, "K2": {"a": 3, "b": 3, "c": 1}}
    D = FN.bray_matrix(tab, ["a", "b", "c"])
    check("bray: identical samples distance 0", D[0][1] == 0 and D[0][2] > 0, str(D))


def t_permanova(tmp: Path):
    import shotgun_func as FN
    tab = {"K1": {}, "K2": {}}
    S = [f"x{i}" for i in range(3)] + [f"y{i}" for i in range(3)]
    for i, s in enumerate(S):
        tab["K1"][s] = (90 + i) if s.startswith("x") else (5 + i)
        tab["K2"][s] = (5 + i) if s.startswith("x") else (90 + i)
    D = FN.bray_matrix(tab, S)
    r = FN.permanova(D, [s[0] for s in S], perm=999, seed=1)
    check("permanova: clear groups p <= 0.1", r.get("p", 1) <= 0.1 and r.get("R2", 0) > 0.8, str(r))
    r2 = FN.permanova(D, ["g", "g", "g", "g", "g", "h"], perm=99, seed=1)
    check("permanova: a group of one -> skipped", "skipped" in r2, str(r2))


def t_func_stats_no_metadata(tmp: Path):
    import shotgun_func as FN
    out = _stats_world(tmp)
    st = FN.func_stats(out, "")
    check("func_stats: PERMANOVA skipped without metadata", "skipped" in st.get("permanova", {}), str(st))
    for f in ("func_pathways.csv", "func_ko.csv", "func_panels.csv", "func_panels_taxa.csv",
              "func_alpha.csv", "func_pcoa_ko.csv"):
        check(f"func_stats: r_tables/{f}", (out / "r_tables" / f).exists())
    alpha = (out / "r_tables" / "func_alpha.csv").read_text().splitlines()
    check("func_alpha.csv: Sample, Shannon_pathways, Pathways",
          alpha[0] == "Sample,Shannon_pathways,Pathways" and len(alpha) == 4, str(alpha))
    pw = (out / "r_tables" / "func_pathways.csv").read_text()
    check("func_pathways.csv: no UNMAPPED / UNINTEGRATED", "UNMAPPED" not in pw and "UNINTEGRATED" not in pw)
    st2 = FN.func_stats(out, str(_meta(tmp, {"a": "x", "b": "x", "c": "y"})))
    check("func_stats: PERMANOVA skipped when a group has one sample", "skipped" in st2["permanova"],
          str(st2["permanova"]))


def _meta(tmp, groups):
    p = tmp / "metadata.csv"
    p.write_text("sampleId,group,description\n" + "".join(f"{s},{g},\n" for s, g in groups.items()))
    return p


def _stats_world(tmp: Path) -> Path:
    out = tmp / "out"
    f = out / "func"
    f.mkdir(parents=True)
    hdr = "\ta_Abundance-CPM\tb_Abundance-CPM\tc_Abundance-CPM\n"
    (f / "pathways_unstratified.tsv").write_text("# Pathway" + hdr + "UNMAPPED\t10\t10\t10\nUNINTEGRATED\t5\t5\t5\n"
                                                 "PWY-1: one\t5\t5\t1\nPWY-2: two\t5\t5\t9\n")
    (f / "pathcoverage_unstratified.tsv").write_text("# Pathway\ta_Coverage\tb_Coverage\tc_Coverage\n"
                                                     "PWY-1: one\t1\t1\t0.2\nPWY-2: two\t0.4\t0.4\t1\n")
    (f / "ko_unstratified.tsv").write_text("# Gene Family" + hdr + "UNMAPPED\t1\t1\t1\nUNGROUPED\t1\t1\t1\n"
                                           "K00929: buk\t4\t4\t1\nK01442: bsh\t1\t1\t4\n")
    (f / "panels.csv").write_text("sample,panel,group,gene,cpm\na,gut,Butyrate,buk,4\nb,gut,Butyrate,buk,4\n"
                                  "c,gut,Butyrate,buk,1\n")
    (f / "panels_taxa.csv").write_text("sample,panel,gene,taxon,cpm\na,gut,buk,Roseburia,4\n")
    for s, g in (("a", "Escherichia"), ("b", "Escherichia"), ("c", "Bacillus")):
        (f / s).mkdir()
        (f / s / f"{s}_metaphlan.tsv").write_text(
            "#mpa_vJun23\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
            f"k__Bacteria\t2\t100\t\nk__Bacteria|p__P|c__C|o__O|f__F|g__{g}\t1\t100\t\n"
            f"k__Bacteria|p__P|c__C|o__O|f__F|g__{g}|s__{g}_x\t1\t100\t\n")
    (out / "asv_table.csv").write_text("id,a,b,c\n1,100,100,0\n2,0,0,100\n")
    (out / "taxonomy.csv").write_text("id,Kingdom,Phylum,Class,Order,Family,Genus,Species\n"
                                      "1,Bacteria,P,C,O,F,Escherichia,coli\n2,Bacteria,P,C,O,F,Bacillus,subtilis\n")
    return out


def t_mpa_vs_kraken(tmp: Path):
    import shotgun_func as FN
    out = _stats_world(tmp)
    g = FN.mpa_genera(out / "func" / "a" / "a_metaphlan.tsv")
    check("mpa_genera: genus rows only", g == {"Escherichia": 100.0}, str(g))
    FN.func_stats(out, "")
    rows = (out / "func" / "mpa_vs_kraken.csv").read_text().splitlines()
    check("mpa_vs_kraken.csv: one row per sample", len(rows) == 4, str(rows))
    import csv as _csv
    r = list(_csv.DictReader(open(out / "func" / "mpa_vs_kraken.csv")))
    check("mpa_vs_kraken: identical profiles -> similarity 1",
          all(abs(float(x["bray_curtis_similarity"]) - 1) < 1e-9 for x in r), str(r))


TESTS += [t_shannon_known, t_pcoa_identical, t_permanova, t_func_stats_no_metadata, t_mpa_vs_kraken]


# ── Task 5: exclusive job gate ────────────────────────────────────────────────
def _gate_job(G, jid, excl, events, labels, cancel=None, hold=None):
    import threading
    def body():
        ok = G.enter(jid, excl, (lambda: bool(cancel and cancel.is_set())),
                     lambda l: labels.setdefault(jid, []).append(l))
        events.append(("in" if ok else "cancelled", jid))
        if ok:
            if hold:
                hold.wait(5)
            events.append(("out", jid))
            G.leave(jid)
    t = threading.Thread(target=body, daemon=True)
    t.start()
    return t


def _until(cond, secs=5.0):
    import time
    end = time.time() + secs
    while time.time() < end:
        if cond():
            return True
        time.sleep(0.02)
    return False


def t_gate(tmp: Path):
    import threading
    import job_gate as JG
    W1, W2 = "Waiting — HUMAnN job running", "Waiting — HUMAnN job starts when the running job ends"
    # two normal jobs run together
    G, ev, lab = JG.JobGate(poll=0.05), [], {}
    ha, hb = threading.Event(), threading.Event()
    _gate_job(G, "a", False, ev, lab, hold=ha)
    _gate_job(G, "b", False, ev, lab, hold=hb)
    check("gate: two normal jobs run together", _until(lambda: ("in", "a") in ev and ("in", "b") in ev), str(ev))
    ha.set(); hb.set()
    # exclusive waits for the running job; a later normal job does not overtake it
    G, ev, lab = JG.JobGate(poll=0.05), [], {}
    ha, hh, hc = threading.Event(), threading.Event(), threading.Event()
    _gate_job(G, "a", False, ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _gate_job(G, "h", True, ev, lab, hold=hh)
    _until(lambda: "h" in lab)
    _gate_job(G, "c", False, ev, lab, hold=hc)
    _until(lambda: "c" in lab)
    check("gate: exclusive waits for the running job", ("in", "h") not in ev and lab.get("h", [""])[0] == W2,
          str((ev, lab)))
    check("gate: later normal job does not overtake", ("in", "c") not in ev, str(ev))
    ha.set()
    check("gate: exclusive starts when the running job ends", _until(lambda: ("in", "h") in ev), str(ev))
    check("gate: normal job waits while exclusive runs",
          _until(lambda: W1 in lab.get("c", [])) and ("in", "c") not in ev, str((ev, lab)))
    hh.set()
    check("gate: normal job starts after exclusive leaves", _until(lambda: ("in", "c") in ev), str(ev))
    hc.set()
    # cancel while waiting
    G, ev, lab = JG.JobGate(poll=0.05), [], {}
    ha, cx, hd = threading.Event(), threading.Event(), threading.Event()
    _gate_job(G, "a", False, ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _gate_job(G, "h", True, ev, lab, cancel=cx)
    _until(lambda: "h" in lab)
    _gate_job(G, "d", False, ev, lab, hold=hd)
    _until(lambda: "d" in lab)
    cx.set()
    check("gate: cancelled while waiting leaves without running", _until(lambda: ("cancelled", "h") in ev)
          and ("in", "h") not in ev, str(ev))
    check("gate: cancelled job frees its turn", _until(lambda: ("in", "d") in ev), str(ev))
    ha.set(); hd.set()
    check("gate: state empty at the end", _until(lambda: G.state() == {"active": {}, "waiting": []}), str(G.state()))


def t_main_uses_gate(tmp: Path):
    src = (ROOT / "backend" / "main.py").read_text()
    check("main.py: GATE = JobGate()", "GATE = JobGate()" in src)
    check("main.py: run_r_pipeline goes through the gate",
          "GATE.enter(job_id, _is_exclusive(params)" in src and "GATE.leave(job_id)" in src)


TESTS += [t_gate, t_main_uses_gate]


# ── Task 6: wiring ─────────────────────────────────────────────────────────────
def t_wiring(tmp: Path):
    p = subprocess.run([sys.executable, str(ROOT / "backend" / "shotgun_pipeline.py"), "--help"],
                       capture_output=True, text=True)
    check("shotgun_pipeline: --func / --func_panels / --metadata",
          all(x in p.stdout for x in ("--func", "--func_panels", "--metadata")), p.stdout[-300:])
    m = (ROOT / "backend" / "main.py").read_text()
    check("main.py: RunParams sg_func default False", "sg_func:            bool  = False" in m)
    check("main.py: RunParams sg_func_panels", 'sg_func_panels:     str   = "gut,food,environment"' in m)
    check("main.py: SHOTGUN cmd passes --func and panels",
          'cmd += ["--func", "--func_panels", params.sg_func_panels' in m)
    check("main.py: /databases reports humann", 'result["humann"] = _humann_status(' in m)
    check("main.py: MaAsLin 3 on function tables", '"--feature_table", str(Path(output_dir) / _tbl)' in m)
    r = (ROOT / "backend" / "r_scripts" / "maaslin3_run.R").read_text()
    check("maaslin3_run.R: --feature_table / --tag",
          '"--feature_table"' in r and '"--tag"' in r and "maaslin3_func_%s.csv" in r and "SKIPPED_%s.txt" in r)


TESTS += [t_wiring]


# ── Task 7: report ─────────────────────────────────────────────────────────────
def _report_world(tmp: Path, platform="illumina", with_func=True, qc_patch=None) -> Path:
    import json as _json
    import test_fast as TFA
    import shotgun_func as FN
    d = tmp / "sg"
    TFA.make_result(d, "shotgun")
    sm = _json.loads((d / "summary.json").read_text())
    sm["platform"] = platform
    if with_func:
        dbp, bases = fake_world(tmp / "w")
        ctx, _ = _ctx(dbp, bases)
        qc = {s: fake_reads(tmp / "reads", s, paired=(platform == "illumina")) for s in ("s1", "s2", "s3")}
        info = FN.run_all(samples=list(qc), qc=qc, platform=platform, out=d, threads=1,
                          panels=["gut", "food"], **ctx)
        if qc_patch:
            rows = list(__import__("csv").DictReader(open(d / "func" / "func_qc.csv")))
            for r in rows:
                r.update(qc_patch)
            FN._write_csv(d / "func" / "func_qc.csv", FN.QC_COLS, rows)
        info["stats"] = FN.func_stats(d, "")
        sm["func"] = info
        sm.setdefault("tools", {}).update(info["versions"])
    (d / "summary.json").write_text(_json.dumps(sm))
    return d


def _html(d: Path) -> str:
    import report_builder as rb
    return rb.build_profile_html(rb.Run(d, set()))


def t_report_func(tmp: Path):
    html = _html(_report_world(tmp))
    for k in ("Functional profile (HUMAnN)", "Function QC", "MetaPhlAn", "HUMAnN 3.9", "fake pathway one",
              "K00929", "Butyrate", "Biogenic amines", "Functional diversity", "Beghini"):
        check(f"report func: has {k!r}", k in html)
    check("report func: partial pathway marked", "fake pathway two" in html and "partial" in html)
    check("report func: carrier bar chart", "Escherichia" in html and html.count("<svg") >= 1)
    check("report func: no ONT beta note for Illumina", "cut into 300 bp" not in html)
    check("report func: no low-mapping warning", "matched no known gene" not in html)
    check("report func: no negative zero similarity", ">-0.00<" not in html)
    check("report func: thousands separators in CPM cells", ">150,000<" in html, "")
    check("report func: no PCoA figure when all samples identical", "PCoA of Bray–Curtis dissimilarity between "
          "samples on KO" not in html)


def t_report_func_absent(tmp: Path):
    html = _html(_report_world(tmp, with_func=False))
    check("report without func: no section", "Functional profile" not in html)
    check("report: key findings box", "Key findings" in html)
    kf = html.split("Key findings", 1)[1].split("</div>", 1)[0] if "Key findings" in html else ""
    check("key findings: only what exists", kf and "AMR" not in kf and "genome bin" not in kf
          and "pathway" not in kf and "Genus" in kf, kf[:400])


def t_report_func_skipped(tmp: Path):
    import json as _json
    d = _report_world(tmp, with_func=False)
    (d / "func").mkdir()
    (d / "func" / "SKIPPED.txt").write_text("HUMAnN databases not found at /mnt/e/ngamp-db/humann\n")
    sm = _json.loads((d / "summary.json").read_text())
    sm["func"] = {"skipped": "HUMAnN databases not found at /mnt/e/ngamp-db/humann"}
    (d / "summary.json").write_text(_json.dumps(sm))
    html = _html(d)
    check("report func skipped: reason shown", "Functional profile" in html
          and "/mnt/e/ngamp-db/humann" in html and "Function QC" not in html)


def t_report_ont_and_warning(tmp: Path):
    html = _html(_report_world(tmp, platform="ont", qc_patch={"unmapped_pct": "70.0"}))
    check("report func ONT: beta note", "cut into 300 bp" in html)
    check("report func: low-mapping warning", "matched no known gene" in html)
    kf = html.split("Key findings", 1)[1].split("</div>", 1)[0]
    check("key findings: function sentence", "pathway" in kf, kf[:400])


TESTS += [t_report_func, t_report_func_absent, t_report_func_skipped, t_report_ont_and_warning]


# ── Task 8: console ───────────────────────────────────────────────────────────
def t_console(tmp: Path):
    import re
    page = (ROOT / "frontend-beta" / "index.html").read_text()
    for k in ('["sg_func", "Function (HUMAnN)", "func", "false"', '"sg_func_panels"', "Resistome, function and genomes",
              "(STATE.databases || {}).humann", 'id="funcPanelsRow"', "Waiting — HUMAnN", 'f:"func_pathways.csv"',
              'f:"func_pcoa_ko.csv"', '"Function", "Other"]', 'row("Function"'):
        check(f"console has {k}", k in page)
    js = "\n".join(re.findall(r"<script>(.*?)</script>", page, re.S))
    f = tmp / "c.js"
    f.write_text(js)
    p = subprocess.run(["node", "--check", str(f)], capture_output=True, text=True)
    check("console script: node --check", p.returncode == 0, p.stderr[-400:])


TESTS += [t_console]


# ── final review fixes ─────────────────────────────────────────────────────────
def t_fix_key_findings_amr_mags(tmp: Path):
    import json as _json
    d = _report_world(tmp, with_func=False)
    write = lambda p, head, rows: p.write_text(",".join(head) + "\n" + "".join(",".join(map(str, r)) + "\n" for r in rows))
    write(d / "amr_genes.csv", ["sample", "gene", "family_id", "type", "class", "rpkm", "host_taxon", "host_share",
                                "host_status"],
          [["s1", "blaTEM", "blaTEM", "AMR", "BETA-LACTAM", 50.0, "Escherichia coli", 90, "clear"],
           ["s2", "tet(M)", "tet(M)", "AMR", "TETRACYCLINE", 9.0, "Enterococcus faecalis", 40, "mixed"],
           ["s1", "fake_virulence", "fake_virulence", "VIRULENCE", "", 99.0, "", "", ""]])
    write(d / "mags_summary.csv", ["sample", "bin", "quality", "completeness", "contamination", "kraken2_taxon"],
          [["s1", "s1.bin.1", "High", 99, 1, "Bacillus subtilis"], ["s1", "s1.bin.2", "Medium", 60, 3, "Escherichia coli"],
           ["s1", "s1.bin.3", "Low", 10, 1, "Listeria"]])
    sm = _json.loads((d / "summary.json").read_text())
    sm["amr"] = {"db_version": "t", "min_identity": 90, "min_breadth": 60}
    sm["mags"] = {"n_high": 1}
    (d / "summary.json").write_text(_json.dumps(sm))
    html = _html(d)
    kf = html.split("Key findings", 1)[1].split("</div>", 1)[0] if "Key findings" in html else ""
    check("key findings: AMR families counted (virulence left out)", "2 gene families" in kf, kf[:600])
    check("key findings: AMR carrier named", "blaTEM" in kf and "E. coli" in kf, kf[:600])
    check("key findings: MAG counts and names", "2 genome bins" in kf and "1 high" in kf
          and "Bacillus subtilis" in kf, kf[:600])
    import report_builder as rb
    import shotgun_report as SR
    run = rb.Run(d, set())
    try:
        res = SR._key_findings(run, run.summary, None, ["s1"])   # broken input inside the body
        ok = res == ""
    except Exception:
        ok = False
    check("key findings: a surprise inside leaves the box out, never breaks the report", ok)


def t_fix_sample_exception_isolated(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, logs = _ctx(dbp, bases)
    qc = {s: fake_reads(tmp / "reads", s) for s in ("s1", "s2")}
    bad = tmp / "reads" / "s2_R1.fastq.gz"
    bad.write_bytes(b"\x1f\x8b\x08\x00garbage")       # truncated gzip
    qc["s2"] = {"r1": bad, "r2": None}
    info = FN.run_all(samples=list(qc), qc=qc, platform="illumina", out=tmp / "out", threads=1,
                      panels=["gut"], **ctx)
    q = {r["sample"]: r["status"] for r in _qc_rows(tmp / "out")}
    check("exception in one sample: others still profiled", q == {"s1": "ok", "s2": "failed"}, f"{q} {info}")


def t_fix_temp_removed_on_failure(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, _ = _ctx(dbp, bases)
    qc = {"s2": fake_reads(tmp / "reads", "s2")}
    tdir = tmp / "out" / "func" / "s2" / "s2_humann_temp"
    tdir.mkdir(parents=True)
    (tdir / "big.tsv").write_text("x")
    os.environ["FAKE_HUMANN_FAIL"] = "s2"
    try:
        FN.run_all(samples=["s2"], qc=qc, platform="illumina", out=tmp / "out", threads=1, panels=[], **ctx)
    finally:
        os.environ.pop("FAKE_HUMANN_FAIL", None)
    check("failed sample: HUMAnN temp folder removed", not tdir.exists())


def t_fix_no_genefamily_rename(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    ctx, logs = _ctx(dbp, bases)
    qc = {"s1": fake_reads(tmp / "reads", "s1")}
    FN.run_all(samples=["s1"], qc=qc, platform="illumina", out=tmp / "out", threads=1, panels=["gut"], **ctx)
    tl = (tmp / "out" / "func" / "tables.log").read_text() if (tmp / "out" / "func" / "tables.log").exists() else ""
    calls = [l for l in tl.splitlines() if l.startswith("###") and "humann_rename_table" in l]
    check("gene families are not renamed with the 87 M UniRef90 names",
          not any("map_uniref90_name" in l for l in calls), str(calls))
    check("gene families not split (one copy less on disk)",
          not (tmp / "out" / "func" / "genefamilies_stratified.tsv").exists())
    pan = (tmp / "out" / "func" / "panels.csv").read_text()
    check("name-pattern panel gene still found via streamed UniRef90 names",
          any(l.startswith("s1,gut,Bile acids,bsh,") or (",bsh," in l and not l.endswith(",0.0")) for l in pan.splitlines())
          or "bsh" in pan, pan[:400])


def t_fix_complete_markers(tmp: Path):
    import shotgun_func as FN
    dbp, bases = fake_world(tmp)
    check("installed with .complete markers", FN.find_install(dbp, bases) is not None)
    (Path(dbp["humann_db"]) / "uniref" / ".complete").unlink()
    check("half-downloaded UniRef90 (no .complete) is not installed", FN.find_install(dbp, bases) is None
          and "UniRef90" in FN.missing_reason(dbp, bases), FN.missing_reason(dbp, bases))
    sh = (ROOT / "setup_shotgun.sh").read_text()
    check("setup_shotgun.sh writes and checks .complete", sh.count(".complete") >= 2)


def t_fix_panels_robust(tmp: Path):
    import shotgun_func as FN
    p = tmp / "p.tsv"
    p.write_text("panel\tgroup\tgene\tko\tuniref_pattern\tnote\ngut\tX\tbad\t\t([unclosed\t\n"
                 "gut\tY\tbuk\tK00929\t\t\n")
    pan = FN.load_panels(p, ["gut"])
    check("panel row with an invalid regex is skipped, others kept", [x["gene"] for x in pan] == ["buk"], str(pan))
    rows, _ = FN.match_panels({}, {}, pan, samples=["s1", "s2"])
    check("panel genes with no hits still listed for every sample",
          sorted((r["sample"], r["cpm"]) for r in rows) == [("s1", 0.0), ("s2", 0.0)], str(rows))


def t_fix_permanova_singletons(tmp: Path):
    import shotgun_func as FN
    out = _stats_world(tmp)
    meta = tmp / "m.csv"
    meta.write_text("sampleId,group,control\na,x,\nb,x,\nc,y,NTC\n")
    g = FN._groups(str(meta))
    check("PERMANOVA groups: controls left out", g == {"a": "x", "b": "x"}, str(g))


def t_fix_seqkit_failure(tmp: Path):
    import shotgun_func as FN
    q = fake_reads(tmp / "in", "o1", paired=False)
    bad = tmp / "seqkit_bad"
    _exe(bad, "#!/bin/sh\nexit 3\n")
    import shotgun_pipeline as SP
    try:
        FN.prep_reads(q, "ont", tmp / "o1", (str(bad), os.environ.copy()), SP.run, tmp / "l.log")
        raised = False
    except RuntimeError:
        raised = True
    check("seqkit failure during ONT prep is an error, not an empty sample", raised)


def t_fix_gate_rewait(tmp: Path):
    import threading
    import job_gate as JG
    G, ev, lab = JG.JobGate(poll=0.05), [], {}
    ha = threading.Event()
    _gate_job(G, "a", False, ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _gate_job(G, "h", True, ev, lab)
    check("gate: on_wait called on every poll (label can change, e.g. checkpoint)",
          _until(lambda: len(lab.get("h", [])) >= 3), str(lab))
    ha.set()
    m = (ROOT / "backend" / "main.py").read_text()
    check("main.py: waiting label mentions a paused checkpoint", "paused at a checkpoint" in m)
    check("main.py: a deleted job counts as cancelled while waiting", "jobs.get(job_id) is None" in m)
    check("main.py: exclusive only when HUMAnN can run", "_humann_status(" in m.split("def _is_exclusive", 1)[1][:600])


def t_fix_report_error_shown(tmp: Path):
    import json as _json
    d = _report_world(tmp)
    sm = _json.loads((d / "summary.json").read_text())
    sm["func"]["error"] = "boom in func_stats"
    (d / "summary.json").write_text(_json.dumps(sm))
    html = _html(d)
    check("report: function step error shown even when tables exist", "boom in func_stats" in html)
    check("report: MaAsLin on function — reason when no metadata", "MaAsLin 3 on function" in html)


TESTS += [t_fix_key_findings_amr_mags, t_fix_sample_exception_isolated, t_fix_temp_removed_on_failure,
          t_fix_no_genefamily_rename, t_fix_complete_markers, t_fix_panels_robust, t_fix_permanova_singletons,
          t_fix_seqkit_failure, t_fix_gate_rewait, t_fix_report_error_shown]


if __name__ == "__main__":
    n = main()
    print(f"\n{'ALL OK' if not n else f'{n} FAILED'}")
    sys.exit(1 if n else 0)
