"""amp_data.py — one finished 12S / CytB amplicon run, read once for the customer
report and the Excel files.

Reads the pipeline's own files: asv_table.csv, taxonomy_table.csv, read_tracking.csv,
cutadapt/<sample>.cutadapt.json (raw / trimmed counts, primer hits), blast/
taxonomy_blast.csv + blast_info.json (BLAST names, NCBI check), run_params.json.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
NA = {"", "NA", "nan", "None"}
TRACK_COLS = ["raw_data", "trimmed", "filtered", "denoisedF", "denoisedR", "merged", "Nonchim", "perc"]


def _v(x) -> str:
    x = (x or "").strip()
    return "" if x in NA else x


def label(asv: dict, rank: str) -> str:
    """Name at a rank, or '<deepest known higher rank> unidentified', or 'Unassigned'."""
    if _v(asv["tax"].get(rank, "")):
        return asv["tax"][rank]
    i = RANKS.index(rank)
    for r in RANKS[:i][::-1]:
        if _v(asv["tax"].get(r, "")):
            return f"{asv['tax'][r]} unidentified"
    return "Unassigned"


def alpha(counts: list[int]) -> dict:
    """phyloseq::estimate_richness values (vegan estimateR bias-corrected Chao1)."""
    c = [x for x in counts if x > 0]
    s, n = len(c), sum(c)
    f1, f2 = sum(1 for x in c if x == 1), sum(1 for x in c if x == 2)
    chao = s + f1 * (f1 - 1) / (2 * (f2 + 1))
    se = math.sqrt(f1 * (f1 - 1) / (2 * (f2 + 1)) + f1 * (2 * f1 - 1) ** 2 / (4 * (f2 + 1) ** 2)
                   + f1 ** 2 * f2 * (f1 - 1) ** 2 / (4 * (f2 + 1) ** 4)) if f1 else 0.0
    p = [x / n for x in c] if n else []
    return {"Observed": s, "Chao1": round(chao, 2), "se.chao1": round(se, 2),
            "Shannon": round(-sum(q * math.log(q) for q in p), 2) if p else 0.0,
            "Simpson": round(1 - sum(q * q for q in p), 2) if p else 0.0}


class RunData:
    def __init__(self, run_dir: Path):
        self.dir = Path(run_dir)
        self.params = self._json("run_params.json")
        self.summary = self._json("summary.json")
        self.blast_info = self._json("blast/blast_info.json")
        self._read_asvs()
        self._add_dropped()
        self._read_tracking()
        self.alpha = {s: alpha([a["counts"].get(s, 0) for a in self.asvs]) for s in self.samples}

    def _json(self, rel):
        try:
            return json.loads((self.dir / rel).read_text())
        except Exception:
            return {}

    def _read_asvs(self):
        with open(self.dir / "asv_table.csv", newline="") as fh:
            rd = csv.reader(fh)
            head = next(rd)
            si = head.index("sequence") if "sequence" in head else 0
            self.samples = [h for i, h in enumerate(head) if i != si]
            rows = [r for r in rd if r]
        tax = {}
        tp = self.dir / "taxonomy_table.csv"
        if tp.exists():
            with open(tp, newline="") as fh:
                rd = csv.reader(fh)
                th = next(rd)
                for r in rd:
                    tax[r[0]] = {h: _v(v) for h, v in zip(th[1:], r[1:])}
        bl = {}
        bp = self.dir / "blast" / "taxonomy_blast.csv"
        if bp.exists():
            with open(bp, newline="") as fh:
                bl = {r["sequence"]: r for r in csv.DictReader(fh)}
        asvs = []
        for r in rows:
            seq = r[si]
            counts = {s: int(float(v or 0)) for s, v in zip(self.samples, [x for i, x in enumerate(r) if i != si])}
            b = bl.get(seq, {})
            t = {k: _v(tax.get(seq, {}).get(k, "")) for k in RANKS}
            asvs.append({"sequence": seq, "counts": counts, "total": sum(counts.values()), "tax": t,
                         **{k: b.get(k, "") for k in ("best_hit", "pident", "qcov", "source", "accession",
                                                     "rank_reached", "target", "ncbi_best", "ncbi_pident",
                                                     "ncbi_accession", "agreement", "check")}})
        asvs.sort(key=lambda a: -a["total"])
        for i, a in enumerate(asvs, 1):
            a["sv"] = f"SV{i:03d}"
        self.asvs = asvs

    def _add_dropped(self):
        """Samples the pipeline dropped (no read left, so no column in asv_table.csv) are
        still part of the order: take them from read_tracking.csv and cutadapt/*.json."""
        known = list(self.samples)
        extra = []
        tp = self.dir / "read_tracking.csv"
        if tp.exists():
            with open(tp, newline="") as fh:
                rd = csv.reader(fh)
                next(rd, None)
                extra += [r[0] for r in rd if r]
        cd = self.dir / "cutadapt"
        if cd.is_dir():
            extra += sorted(p.name[:-len(".cutadapt.json")] for p in cd.glob("*.cutadapt.json"))
        self.dropped = []
        for s in extra:
            if s and s not in known:
                known.append(s)
                self.dropped.append(s)
        self.samples = known

    def _read_tracking(self):
        self.tracking = {}
        tp = self.dir / "read_tracking.csv"
        rows = {}
        if tp.exists():
            with open(tp, newline="") as fh:
                rd = csv.reader(fh)
                head = next(rd)
                for r in rd:
                    rows[r[0]] = dict(zip(head[1:], r[1:]))
        for s in self.samples:
            t = rows.get(s, {})
            num = lambda k: int(float(t.get(k) or 0))
            cj = self._json(f"cutadapt/{s}.cutadapt.json").get("read_counts", {})
            # raw = reads before primer removal (cutadapt input). Without cutadapt's own
            # counts it is unknown (None) — filterAndTrim's input is the trimmed count.
            raw = int(cj["input"]) if cj.get("input") is not None else None
            trimmed = int(cj.get("output") or 0) or num("input")
            nonchim = sum(a["counts"].get(s, 0) for a in self.asvs)
            self.tracking[s] = {"raw_data": raw, "trimmed": trimmed, "filtered": num("filtered"),
                                "denoisedF": num("denoisedF") or num("denoised"), "denoisedR": num("denoisedR"),
                                "merged": num("merged"), "Nonchim": nonchim,
                                "perc": (round(100 * nonchim / raw, 2) if raw else 0.0) if raw is not None else None,
                                "primer_hits": int(cj.get("read1_with_adapter") or 0),
                                "cutadapt_version": self._json(f"cutadapt/{s}.cutadapt.json").get("cutadapt_version", "")}

    def sample_asvs(self, s: str) -> list[dict]:
        return [a for a in self.asvs if a["counts"].get(s, 0) > 0]

    def rank_table(self, rank: str, s: str | None = None) -> dict:
        """{label: reads} for one sample, or for all samples summed."""
        out = {}
        for a in self.asvs:
            n = a["counts"].get(s, 0) if s else a["total"]
            if n:
                k = label(a, rank)
                out[k] = out.get(k, 0) + n
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    def non_target_share(self, s: str) -> float:
        tot = sum(a["counts"].get(s, 0) for a in self.asvs) or 1
        return 100 * sum(a["counts"].get(s, 0) for a in self.asvs if a.get("target") == "no") / tot


def load_run(run_dir: Path) -> RunData:
    return RunData(Path(run_dir))
