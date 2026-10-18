#!/usr/bin/env python3
"""controls.py — negative controls and mock communities in a result, any pipeline.

Which samples are controls comes from the run itself:
  * params.controlSamples / params.mockSamples (the run form), or
  * a metadata column named control, sample_type or type with values such as
    negative / blank / NTC (negative) or mock / positive / zymo (mock).

negative_summary() shows what the blanks contain and which genera are more
abundant in the blanks than in the samples — the classic contaminant signature
(Davis et al. 2018, Salter et al. 2014) — without removing anything.
mock_summary() compares a mock sample with the theoretical composition of its
standard. Both work on report_builder.Analysis, so DADA2, ITS, COX1, PacBio, Emu
and VSEARCH results are all handled the same way.
"""
from __future__ import annotations

NEG_WORDS = {"negative", "neg", "blank", "ntc", "extraction blank", "pcr blank", "kit blank",
             "negative control", "nc", "water"}
MOCK_WORDS = {"mock", "positive", "positive control", "pc", "zymo", "standard", "mock community"}

# Theoretical compositions (genus level, % of 16S copies).
STANDARDS = {
    "zymo_d6300": {
        "label": "ZymoBIOMICS Microbial Community Standard (D6300 / D6305 / D6306), 16S theoretical",
        "genera": {"Pseudomonas": 4.2, "Escherichia": 10.1, "Salmonella": 10.4,
                   "Limosilactobacillus": 18.4, "Enterococcus": 9.9, "Staphylococcus": 15.5,
                   "Listeria": 14.1, "Bacillus": 17.4},
        "synonyms": {"Limosilactobacillus": ["Lactobacillus"],
                     "Escherichia": ["Escherichia-Shigella", "Escherichia/Shigella", "Shigella"]},
        "family": {"Escherichia": "Enterobacteriaceae", "Salmonella": "Enterobacteriaceae"},
    },
    "zymo_d6331": {
        "label": "ZymoBIOMICS Gut Microbiome Standard (D6331), 16S — genera listed, proportions not given",
        "genera": {},
    },
}


def identify(run) -> tuple[list[str], list[str]]:
    """(negative control samples, mock samples) present in this run."""
    S = set(run.samples)
    neg = [s for s in (run.params.get("controlSamples") or []) if s in S]
    mock = [s for s in (run.params.get("mockSamples") or []) if s in S]
    rows = list(run.params.get("metadata") or []) or run._rows("metadata.csv")
    for r in rows:
        keys = {str(k).strip().lower(): v for k, v in dict(r).items()}
        sid = str(keys.get("sampleid") or keys.get("sample_id") or keys.get("sample")
                  or (list(dict(r).values())[0] if r else "")).strip()
        if sid not in S:
            continue
        val = ""
        for k in ("control", "sample_type", "type"):
            if keys.get(k):
                val = str(keys[k]).strip().lower()
                break
        if val in NEG_WORDS and sid not in neg:
            neg.append(sid)
        elif val in MOCK_WORDS and sid not in mock:
            mock.append(sid)
    return neg, mock


def negative_summary(an, negs: list[str], mocks=(), top: int = 12) -> dict:
    S = [s for s in an.samples if s not in negs and s not in mocks]
    pct = an.rank_pct("Genus")
    reads = {s: an.comm_tot.get(s, 0) for s in an.samples}
    real_reads = sorted(reads[s] for s in S) or [0]
    med = real_reads[len(real_reads) // 2]
    neg_tot = sum(reads[s] for s in negs) or 1
    cnt = an.rank_counts("Genus")
    pooled = {g: 100 * sum(v[s] for s in negs) / neg_tot for g, v in cnt.items()}
    rows = []
    for g, p in sorted(pooled.items(), key=lambda kv: -kv[1])[:top]:
        if p <= 0:
            continue
        mean_s = sum(pct[g][s] for s in S) / len(S) if S else 0.0
        prev = sum(1 for s in S if pct[g][s] >= 0.1)
        rows.append({"genus": g, "pct_in_blanks": p, "mean_pct_in_samples": mean_s,
                     "samples_present": prev, "n_samples": len(S),
                     # twice as abundant in the blanks as in the samples, and not trivial
                     "flag": g != "Unclassified" and p >= 1.0 and p >= 2 * mean_s})
    return {"negatives": [{"sample": s, "reads": reads[s],
                           "ratio_to_median": (reads[s] / med) if med else None} for s in negs],
            "median_sample_reads": med, "genera": rows,
            "flagged": [r["genus"] for r in rows if r["flag"]]}


def mock_summary(an, mocks: list[str], standard: str) -> dict:
    std = STANDARDS.get(standard) or {}
    exp = std.get("genera") or {}
    if not exp:
        return {"standard": std.get("label", standard), "samples": [], "no_composition": True}
    syn = {g: {g, *std.get("synonyms", {}).get(g, [])} for g in exp}
    pct = an.rank_pct("Genus")
    fam = an.rank_pct("Family") if "Family" in an.run.rank_cols else {}
    out = []
    tot_exp = sum(exp.values())
    for s in mocks:
        obs = {}
        other = 0.0
        for t, v in pct.items():
            if t == "Unclassified":
                continue
            key = next((g for g in exp if t in syn[g]), None)
            if key:
                obs[key] = obs.get(key, 0) + v[s]
            else:
                other += v[s]
        rows, num = [], 0.0
        for g, e in exp.items():
            e_pct = 100 * e / tot_exp
            o = obs.get(g, 0.0)
            fam_only = o < 0.25 * e_pct and std.get("family", {}).get(g) in fam and \
                fam[std["family"][g]][s] >= 0.25 * e_pct
            rows.append({"genus": g, "expected": e_pct, "observed": o, "family_only": fam_only,
                         "ratio": (o / e_pct) if e_pct else None})
            num += abs(e_pct - o)
        bc = (num + other) / 200.0
        out.append({"sample": s, "bray_curtis": bc, "other_pct": other,
                    "unclassified_pct": pct.get("Unclassified", {}).get(s, 0.0), "genera": rows,
                    "missing": [r["genus"] for r in rows if r["observed"] < 0.25 * r["expected"]
                                and not r["family_only"]]})
    return {"standard": std["label"], "samples": out}
