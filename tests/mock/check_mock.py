#!/usr/bin/env python3
"""Compare a pipeline's result folder with the known truth of a simulated mock.

    python3 check_mock.py RESULT_DIR TRUTH_JSON [--rank Genus] [--json out.json]

Reads the result through report_builder.Run — the same code the PDF/Word report
uses — so a table the report would misread fails here too. Checks, per sample:

  * every expected genus is found (at least 25 % of its true share) — or, for
    members a short amplicon cannot separate (Escherichia / Salmonella), at least
    its family,
  * reads assigned to taxa that are NOT in the mock stay under 2 %,
  * Bray-Curtis between observed and true FAMILY composition is at most 0.15
    (genus-level Bray-Curtis is reported too),
  * the Genus column holds genera — not family names ("...aceae") or binomials,
    which is what a rank shift looks like (v2.14 SILVA 144 / COX1 bugs).

Exit status 0 = pass, 1 = fail. Thresholds are deliberately loose: the reads are
simulated without PCR bias, so a real failure is far outside them.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "backend"))

import report_builder as rb  # noqa: E402

# family of each mock member's genus: two members of one family (Escherichia /
# Salmonella, Aspergillus / Penicillium) can be inseparable on a short amplicon,
# and a pipeline that stops at the family for them is right, not broken.
FAMILY = {
    "Pseudomonas": "Pseudomonadaceae", "Escherichia": "Enterobacteriaceae",
    "Salmonella": "Enterobacteriaceae", "Limosilactobacillus": "Lactobacillaceae",
    "Enterococcus": "Enterococcaceae", "Staphylococcus": "Staphylococcaceae",
    "Listeria": "Listeriaceae", "Bacillus": "Bacillaceae",
    "Saccharomyces": "Saccharomycetaceae", "Cryptococcus": "Cryptococcaceae",
    "Aspergillus": "Aspergillaceae", "Penicillium": "Aspergillaceae",
    "Candida": "Debaryomycetaceae",
    "Drosophila": "Drosophilidae", "Apis": "Apidae", "Musca": "Muscidae",
    "Bombyx": "Bombycidae", "Aedes": "Culicidae",
}
FAMILY_SYN = {"Debaryomycetaceae": {"Debaryomycetaceae", "Saccharomycetales_fam_Incertae_sedis"},
              "Cryptococcaceae": {"Cryptococcaceae", "Tremellaceae"}}

# genus names a reference may use for a mock member
SYNONYMS = {
    "Limosilactobacillus": {"Limosilactobacillus", "Lactobacillus"},
    "Escherichia": {"Escherichia", "Escherichia-Shigella", "Escherichia/Shigella", "Shigella"},
    "Salmonella": {"Salmonella"},
    "Cryptococcus": {"Cryptococcus", "Filobasidiella"},
    "Candida": {"Candida", "Candida/Lodderomyces clade", "Lodderomyces"},
}


def genus_of(member: str) -> str:
    return member.split("_")[0]


def accepted(member: str) -> set[str]:
    g = genus_of(member)
    return SYNONYMS.get(g, {g})


def bray(a: dict, b: dict) -> float:
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.values()) + sum(b.values())
    return num / den if den else 0.0


def check(result_dir: Path, truth: dict, rank: str = "Genus") -> dict:
    run = rb.Run(result_dir, set())
    an = rb.Analysis(run)
    obs_pct = an.rank_pct(rank)
    report = {"kind": run.kind, "rank_cols": run.rank_cols, "samples": {}, "problems": []}

    # rank-shift detector over the whole table
    gv = [str(r.get("Genus") or "").strip() for r in run.tax.values()]
    gv = [g for g in gv if g and not rb.BLANK_RE.match(g)]
    if gv:
        fam = sum(bool(re.search(r"aceae$", g)) for g in gv) / len(gv)
        binom = sum(" " in g and not g.startswith("Candida") for g in gv) / len(gv)
        report["genus_column"] = {"family_like": round(fam, 3), "binomial_like": round(binom, 3)}
        if fam > 0.2 or binom > 0.2:
            report["problems"].append(f"Genus column looks shifted (family-like {fam:.0%}, "
                                      f"binomial {binom:.0%})")
    else:
        report["problems"].append("no genus assignments at all")

    fam_pct = an.rank_pct("Family") if "Family" in run.rank_cols else {}
    for s, counts in truth["samples"].items():
        if s not in an.samples:
            report["problems"].append(f"sample {s} missing from the result")
            continue
        tot = sum(counts.values()) or 1
        exp, expf = {}, {}
        for m, c in counts.items():
            g = genus_of(m)
            exp[g] = exp.get(g, 0) + 100 * c / tot
            f = FAMILY.get(g, g)
            expf[f] = expf.get(f, 0) + 100 * c / tot
        obs = {}
        for t, v in obs_pct.items():
            if t == "Unclassified":
                continue
            key = next((genus_of(m) for m in counts if t in accepted(m)), t)
            obs[key] = obs.get(key, 0) + v[s]
        obsf = {}
        for t, v in fam_pct.items():
            if t == "Unclassified":
                continue
            key = next((f for f in expf if t in FAMILY_SYN.get(f, {f})), t)
            obsf[key] = obsf.get(key, 0) + v[s]
        unexpected = sum(v for k, v in obs.items() if k not in exp)
        unexpected_f = sum(v for k, v in obsf.items() if k not in expf) if obsf else 0.0
        missing = [g for g, e in exp.items() if obs.get(g, 0) < 0.25 * e]
        family_only = [g for g in missing if obsf.get(FAMILY.get(g, g), 0) >= 0.25 * exp[g]]
        missing = [g for g in missing if g not in family_only]
        unclassified = obs_pct.get("Unclassified", {}).get(s, 0.0)
        bc = bray(expf, obsf) if obsf else bray(exp, obs)
        row = {"bray_curtis_family": round(bc, 3), "bray_curtis_genus": round(bray(exp, obs), 3),
               "unexpected_pct": round(max(unexpected, unexpected_f), 2),
               "unclassified_pct": round(unclassified, 2), "missing": missing,
               "family_only": family_only,
               "expected": {k: round(v, 2) for k, v in exp.items()},
               "observed": {k: round(v, 2) for k, v in sorted(obs.items(), key=lambda kv: -kv[1])[:15]}}
        report["samples"][s] = row
        if missing:
            report["problems"].append(f"{s}: expected taxa not found even at family: {', '.join(missing)}")
        if row["unexpected_pct"] > 2.0:
            report["problems"].append(f"{s}: {row['unexpected_pct']:.1f}% of reads in taxa not in the mock")
        if bc > 0.15:
            report["problems"].append(f"{s}: Bray-Curtis (family) to the true composition {bc:.3f} > 0.15")
    report["pass"] = not report["problems"]
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("result_dir")
    ap.add_argument("truth")
    ap.add_argument("--rank", default="Genus")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    rep = check(Path(a.result_dir), json.loads(Path(a.truth).read_text()), a.rank)
    if a.json:
        Path(a.json).write_text(json.dumps(rep, indent=2))
    print(f"pipeline: {rep['kind']}   ranks: {', '.join(rep['rank_cols'])}")
    for s, r in rep["samples"].items():
        print(f"  {s}: Bray-Curtis family {r['bray_curtis_family']:.3f} / genus "
              f"{r['bray_curtis_genus']:.3f}, unexpected {r['unexpected_pct']:.2f}%, "
              f"unclassified {r['unclassified_pct']:.2f}%"
              + (f", family only: {r['family_only']}" if r["family_only"] else "")
              + (f", MISSING {r['missing']}" if r["missing"] else ""))
        for g, e in r["expected"].items():
            print(f"      {g:22s} expected {e:6.2f}%   observed {r['observed'].get(g, 0):6.2f}%")
    print("PASS" if rep["pass"] else "FAIL:\n  " + "\n  ".join(rep["problems"]))
    return 0 if rep["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
