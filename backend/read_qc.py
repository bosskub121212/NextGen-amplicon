#!/usr/bin/env python3
"""read_qc.py — read quality before analysis, and truncLen from the data.

For every sample in an upload folder this reads the first N reads of each file
(20,000 by default) and records:

  * Illumina: quality per position (10th / 25th / 50th / 75th percentile, mean)
    for R1 and R2, read length, mean quality, share of bases >= Q30;
  * long reads (ONT / PacBio): read-length distribution, N50, mean read quality
    distribution (what NanoPlot shows).

For paired Illumina amplicons it then SUGGESTS truncLen: each read is cut where
its median quality falls below Q25 for good, and the two cuts are checked
against the amplicon length (from the primer pair, or --expected_len) so the
pairs still overlap by at least 20 bases after cutting; when quality and overlap
conflict, overlap wins and the report says so.

  python3 read_qc.py INPUT_DIR [--out qc.json] [--primer_f SEQ --primer_r SEQ]
                     [--expected_len 440] [--n 20000] [--fastqc OUT_DIR]

FastQC / MultiQC are run as well when they are installed and --fastqc is given;
the numbers above never depend on them. Standard library only.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

PAIR_RE = re.compile(r"^(?P<name>.+?)(?:_S\d+)?(?:_L\d{3})?[_.](?:R)?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$")
FQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$")

# Insert length (amplicon without primers) of common primer pairs — the longest
# taxa commonly seen, so an overlap check against it is safe for the run.
INSERT_LEN = {
    ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"): 440,      # 341F / 805R  V3-V4
    ("CCTACGGGNGGCWGCAG", "GGACTACHVGGGTWTCTAAT"): 430,       # 341F / 785R  V3-V4
    ("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"): 256,     # 515F / 806R  V4
    ("GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT"): 256,     # 515F / 806R  V4 (EMP v1)
    ("AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"): 310,       # 1055F / 1392R V7-V8
    ("GTCGGTAAAACTCGTGCCAGC", "CATAGTGGGGTATCTAATCCCAGTTTG"): 185,  # MiFish-U 12S
    ("CCAGCASCYGCGGTAATTCC", "ACTTTCGTTCTTGATYRA"): 370,      # 18S V4 (TAReuk / nematode)
    ("GGWACWGGWTGAACWGTWTAYCCYCC", "TANACYTCNGGRTGNCCRAARAAYCA"): 313,  # mlCOIintF / jgHCO2198
}
ITS_PRIMERS = {"CTTGGTCATTTAGAGGAAGTAA", "TCCGTAGGTGAACCTGCGG", "GTGAATCATCGAATCTTTGAA",
               "GCTGCGTTCTTCATCGATGC", "TCCTCCGCTTATTGATATGC"}
MIN_OVERLAP = 20


def fastq_iter(path: Path, n: int):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        k = 0
        while k < n:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            if not q:
                break
            k += 1
            yield s, q


def pct(sorted_vals, p):
    if not sorted_vals:
        return 0
    i = min(len(sorted_vals) - 1, max(0, int(round(p * (len(sorted_vals) - 1)))))
    return sorted_vals[i]


def hist_pct(h, tot, p):
    """Percentile of a quality histogram (index = Q, value = count)."""
    target, acc = p * tot, 0
    for qv, c in enumerate(h):
        acc += c
        if acc >= target:
            return qv
    return len(h) - 1


def profile(path: Path, n: int) -> dict:
    """Per-file statistics from the first n reads."""
    lens, meanq, q30 = [], [], 0
    bases = 0
    pos_hist: list[list[int]] = []          # pos -> histogram over Q 0..60
    for s, q in fastq_iter(path, n):
        L = len(s)
        lens.append(L)
        qs = [ord(c) - 33 for c in q]
        if qs:
            # mean of error probabilities, as NanoPlot does, then back to Q
            pe = sum(10 ** (-x / 10) for x in qs) / len(qs)
            meanq.append(-10 * math.log10(max(pe, 1e-10)))
            q30 += sum(1 for x in qs if x >= 30)
            bases += len(qs)
        if L <= 1000:                        # per-position only for short reads
            while len(pos_hist) < L:
                pos_hist.append([0] * 61)
            for i, x in enumerate(qs):
                pos_hist[i][min(60, max(0, x))] += 1
    out = {"file": path.name, "reads_examined": len(lens)}
    if not lens:
        return out
    sl = sorted(lens)
    out.update(mean_len=round(sum(lens) / len(lens), 1), min_len=sl[0], max_len=sl[-1],
               median_len=sl[len(sl) // 2],
               mean_q=round(statistics.fmean(meanq), 1) if meanq else None,
               pct_q30=round(100 * q30 / bases, 1) if bases else None)
    long_reads = statistics.median(lens) > 600
    out["long_reads"] = long_reads
    if long_reads:
        half, acc = sum(lens) / 2, 0
        for L in sorted(lens, reverse=True):
            acc += L
            if acc >= half:
                out["n50"] = L
                break
        # histograms for plotting
        lb = [0] * 40                              # 100 bp bins to 4 kb
        for L in lens:
            lb[min(39, L // 100)] += 1
        out["len_hist_100bp"] = lb
        qb = [0] * 41
        for x in meanq:
            qb[min(40, int(x))] += 1
        out["meanq_hist"] = qb
    else:
        prof = []
        for h in pos_hist:
            tot = sum(h)
            if tot < 0.5 * len(lens):              # fewer than half the reads this long
                break
            prof.append([hist_pct(h, tot, 0.10), hist_pct(h, tot, 0.25), hist_pct(h, tot, 0.50),
                         hist_pct(h, tot, 0.75), round(sum(qv * c for qv, c in enumerate(h)) / tot, 1)])
        out["per_position"] = prof               # [p10, p25, median, p75, mean] per base
    return out


def quality_cut(prof: list, primer_len: int, floor: int = 25, window: int = 10) -> int:
    """Last position (in primer-trimmed coordinates) before the median quality
    drops below `floor` for good (a `window`-base rolling median)."""
    if not prof:
        return 0
    med = [p[2] for p in prof]
    cut = len(med)
    for i in range(len(med)):
        w = med[i:i + window]
        if len(w) == window and statistics.median(w) < floor:
            cut = i
            break
    return max(0, cut - primer_len)


def expected_insert(primer_f: str, primer_r: str) -> int | None:
    f, r = primer_f.upper(), primer_r.upper()
    if f in ITS_PRIMERS or r in ITS_PRIMERS:
        return None
    return INSERT_LEN.get((f, r))


def suggest_trunc(r1: dict, r2: dict, primer_f: str, primer_r: str, insert: int | None) -> dict:
    """truncLen F/R from the pooled R1 / R2 quality profiles."""
    pf, pr = len(primer_f or ""), len(primer_r or "")
    # filterAndTrim drops any read SHORTER than truncLen, so stay a few bases
    # under the typical trimmed length instead of exactly on it
    lenF = max(0, int(r1.get("median_len", 0)) - pf - 3)
    lenR = max(0, int(r2.get("median_len", 0)) - pr - 3)
    cf = min(lenF, quality_cut(r1.get("per_position", []), pf))
    cr = min(lenR, quality_cut(r2.get("per_position", []), pr))
    res = {"quality_cut_F": cf, "quality_cut_R": cr, "read_len_F_trimmed": lenF,
           "read_len_R_trimmed": lenR, "insert_len": insert, "min_overlap": MIN_OVERLAP}
    if primer_f.upper() in ITS_PRIMERS or primer_r.upper() in ITS_PRIMERS:
        res.update(truncLen_F=0, truncLen_R=0,
                   reason="ITS length varies between taxa, so reads are not truncated (maxEE filters them).")
        return res
    if insert is None:
        res.update(truncLen_F=cf, truncLen_R=cr,
                   reason=f"Cut where the median quality falls below Q25. The amplicon length is not known for "
                          f"these primers, so check that {cf} + {cr} bp still overlaps by {MIN_OVERLAP} bp.")
        return res
    need = insert + MIN_OVERLAP
    tf, tr = cf, cr
    if tf + tr < need:
        # give back length, worse read last, until the pairs overlap again
        short = need - (tf + tr)
        room_f, room_r = lenF - tf, lenR - tr
        add_f = min(room_f, (short + 1) // 2)
        add_r = min(room_r, short - add_f)
        add_f = min(room_f, short - add_r)
        tf, tr = tf + add_f, tr + add_r
    ok = tf + tr >= need
    res.update(truncLen_F=tf, truncLen_R=tr, overlap=tf + tr - insert, ok=ok)
    if not ok:
        res["reason"] = (f"Even untruncated the reads ({lenF} + {lenR} bp after primers) cannot overlap a "
                         f"{insert} bp insert by {MIN_OVERLAP} bp: read-pair merging will fail for most reads.")
    elif (tf, tr) == (cf, cr):
        res["reason"] = (f"Cut where the median quality falls below Q25; the reads still overlap by "
                         f"{tf + tr - insert} bp on a {insert} bp insert.")
    else:
        res["reason"] = (f"Quality alone would cut at {cf} / {cr} bp, which leaves no overlap on a {insert} bp "
                         f"insert; lengthened to {tf} / {tr} bp to keep {tf + tr - insert} bp of overlap "
                         f"(maxEE removes the worst reads).")
    return res


def find_samples(d: Path) -> dict:
    """{sample: [R1] or [R1, R2]} from the manifest, else from file names."""
    man = d / "sample_manifest.json"
    if man.exists():
        try:
            rows = json.loads(man.read_text())
            out = {}
            for r in rows:
                fs = [d / r[k] for k in ("file1", "file2") if r.get(k) and (d / r[k]).exists()]
                if r.get("sample") and fs:
                    out[r["sample"]] = fs
            if out:
                return out
        except Exception:
            pass
    files = sorted(f for f in d.rglob("*") if f.is_file() and FQ_RE.search(f.name)
                   and "_wgs_work" not in f.parts)
    pairs, single = {}, {}
    for f in files:
        m = PAIR_RE.match(f.name)
        if m:
            pairs.setdefault(m.group("name"), {})[m.group("mate")] = f
        else:
            single[FQ_RE.sub("", f.name)] = [f]
    out = {k: [v["1"], v["2"]] for k, v in pairs.items() if set(v) == {"1", "2"}}
    for k, v in pairs.items():
        if set(v) != {"1", "2"}:
            for m_, f in v.items():
                out[FQ_RE.sub("", f.name)] = [f]
    out.update(single)
    return dict(sorted(out.items()))


def pool(profiles: list[dict]) -> dict:
    """Read-weighted average of per-position profiles (for the run-level suggestion)."""
    ps = [p for p in profiles if p.get("per_position")]
    if not ps:
        return {}
    L = min(len(p["per_position"]) for p in ps)
    w = [p["reads_examined"] for p in ps]
    prof = []
    for i in range(L):
        prof.append([round(sum(p["per_position"][i][k] * wi for p, wi in zip(ps, w)) / sum(w), 1)
                     for k in range(5)])
    return {"per_position": prof,
            "median_len": int(statistics.median([p["median_len"] for p in ps])),
            "reads_examined": sum(w)}


def find_tool(name):
    hit = shutil.which(name)
    if hit:
        return hit
    for b in ("miniconda3", "anaconda3", "miniforge3", "mambaforge"):
        for e in ("emu", "ngamp-wgs", "qc"):
            c = Path.home() / b / "envs" / e / "bin" / name
            if c.exists():
                return str(c)
    return None


def run(input_dir: str, primer_f: str = "", primer_r: str = "", expected_len: int | None = None,
        n: int = 20000, fastqc_dir: str | None = None, threads: int = 2) -> dict:
    d = Path(input_dir)
    samples = find_samples(d)
    # MinKNOW folders hold hundreds of chunk files: look at no more than 60 of
    # them and at most ~400,000 reads in all, so this stays a few seconds.
    nfiles = sum(len(v) for v in samples.values())
    if nfiles > 60:
        samples = dict(list(samples.items())[:60])
        nfiles = sum(len(v) for v in samples.values())
    n = max(2000, min(n, 400000 // max(1, nfiles)))
    res = {"input_dir": str(d), "reads_per_file_examined": n, "samples": {}, "tools": {}}
    r1s, r2s = [], []
    for name, files in samples.items():
        ent = {"files": [profile(f, n) for f in files]}
        res["samples"][name] = ent
        if len(files) == 2:
            r1s.append(ent["files"][0])
            r2s.append(ent["files"][1])
    long_reads = any(f.get("long_reads") for s in res["samples"].values() for f in s["files"])
    res["platform"] = "long" if long_reads else ("paired" if r1s else "single")
    if r1s and not long_reads:
        P1, P2 = pool(r1s), pool(r2s)
        res["pooled"] = {"R1": P1, "R2": P2}
        ins = expected_len or expected_insert(primer_f, primer_r)
        res["suggestion"] = suggest_trunc(P1, P2, primer_f, primer_r, ins)
    # low-depth samples are worth flagging before a run, not after
    depth = {}
    for name, ent in res["samples"].items():
        depth[name] = ent["files"][0].get("reads_examined", 0)
    res["low_depth"] = [k for k, v in depth.items() if v < min(n, 1000)]
    if fastqc_dir:
        fq, mq = find_tool("fastqc"), find_tool("multiqc")
        if fq:
            od = Path(fastqc_dir)
            od.mkdir(parents=True, exist_ok=True)
            allf = [str(f) for fs in samples.values() for f in fs]
            subprocess.run([fq, "-q", "-t", str(threads), "-o", str(od)] + allf,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            res["tools"]["fastqc"] = str(od)
            if mq:
                subprocess.run([mq, "-q", "-f", "-o", str(od), str(od)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                if (od / "multiqc_report.html").exists():
                    res["tools"]["multiqc"] = str(od / "multiqc_report.html")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input_dir")
    ap.add_argument("--out")
    ap.add_argument("--primer_f", default="")
    ap.add_argument("--primer_r", default="")
    ap.add_argument("--expected_len", type=int)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--fastqc", help="also run FastQC + MultiQC into this folder when installed")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args(argv)
    res = run(a.input_dir, a.primer_f, a.primer_r, a.expected_len, a.n, a.fastqc, a.threads)
    txt = json.dumps(res, indent=1)
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(txt)
    s = res.get("suggestion")
    print(f"{len(res['samples'])} samples, {res['platform']} reads")
    if s:
        print(f"suggested truncLen {s['truncLen_F']} / {s['truncLen_R']} — {s['reason']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
