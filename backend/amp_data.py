"""amp_data.py — one finished amplicon run, read once for the customer report and the
Excel files.

kind "dada2" (16S, 12S / CytB): asv_table.csv, taxonomy_table.csv, read_tracking.csv,
cutadapt/<sample>.cutadapt.json (raw / trimmed counts, primer hits), blast/
taxonomy_blast.csv + blast_info.json (BLAST names, NCBI check), run_params.json.
kind "emu" (ONT 16S): asv_table.csv (tax_id rows), taxonomy.csv, read_tracking.csv
(input, qc_filtered, primer_trimmed, final), summary.json.

only(sample) / subset(samples) give a view of the run holding nothing of the other
samples — the basis of the one-package-per-sample customer report.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
NA = {"", "NA", "nan", "None"}
TRACK_COLS = ["raw_data", "trimmed", "filtered", "denoisedF", "denoisedR", "merged", "Nonchim", "perc"]
EMU_TRACK_COLS = ["raw_data", "qc_filtered", "primer_trimmed", "final", "assigned", "perc"]
# Emu rows that are not a taxon: shown as "Unassigned" in the profiles, never counted as
# assigned reads or as a taxon in richness.
EMU_NOT_TAXA = {"unassigned", "unmapped", "mapped_unclassified"}


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
        self.kind = self._kind()
        self.track_cols = EMU_TRACK_COLS if self.kind == "emu" else TRACK_COLS
        if self.kind == "emu":
            self._read_emu()
        else:
            self._read_asvs()
        self._add_dropped()
        if self.kind == "emu":
            self._read_emu_tracking()
        else:
            self._read_tracking()
        self.peer_names = list(self.samples)          # every sample of the run (name lookups only)
        self.alpha = {s: alpha([a["counts"].get(s, 0) for a in self.asvs if not a.get("not_taxon")])
                      for s in self.samples}

    def _kind(self) -> str:
        sm = self.summary
        if sm.get("pipeline") == "emu" or (sm.get("platform") == "ONT" and not (self.dir / "taxonomy_table.csv").exists()):
            return "emu"
        try:
            head = (self.dir / "asv_table.csv").open().readline()
        except OSError:
            head = ""
        if "sequence" not in head and (self.dir / "taxonomy.csv").exists():
            return "emu"
        return "dada2"

    @property
    def counts_are_reads(self) -> bool:
        if self.kind != "emu":
            return True
        kind = self.summary.get("emu_counts")
        if kind:
            return kind != "ppm"
        # summary.json from before emu_counts was recorded: ppm columns sum to one million
        sums = [sum(a["counts"].get(s, 0) for a in getattr(self, "asvs", [])) for s in self.samples]
        sums = [x for x in sums if x]
        return not (sums and all(abs(x - 1_000_000) <= 10_000 for x in sums))

    def _read_emu(self):
        with open(self.dir / "asv_table.csv", newline="") as fh:
            rd = csv.reader(fh)
            head = next(rd)
            self.samples = [h for h in head[1:]]
            rows = [r for r in rd if r]
        tax = {}
        tp = self.dir / "taxonomy.csv"
        if tp.exists():
            with open(tp, newline="") as fh:
                for r in csv.DictReader(fh):
                    t = {k: _v(r.get(k, "")) for k in RANKS}
                    # emu_pipeline writes "Genus sp." when Emu gave no species: not a species name
                    if t["Species"].endswith(" sp.") or t["Species"] == "sp.":
                        t["Species"] = ""
                    tax[r.get("ASV", "")] = t
        asvs = []
        for r in rows:
            tid = r[0]
            counts = {s: int(round(float(v or 0))) for s, v in zip(self.samples, r[1:])}
            asvs.append({"id": tid, "sequence": "", "counts": counts, "total": sum(counts.values()),
                         "not_taxon": tid.strip().lower() in EMU_NOT_TAXA,
                         "tax": tax.get(tid, {k: "" for k in RANKS}),
                         **{k: "" for k in ("best_hit", "pident", "qcov", "source", "accession", "rank_reached",
                                            "target", "ncbi_best", "ncbi_pident", "ncbi_accession", "agreement",
                                            "check")}})
        asvs.sort(key=lambda a: -a["total"])
        for i, a in enumerate(asvs, 1):
            a["sv"] = f"SV{i:03d}"
        self.asvs = asvs

    def _read_emu_tracking(self):
        rows = {}
        tp = self.dir / "read_tracking.csv"
        if tp.exists():
            with open(tp, newline="") as fh:
                for r in csv.DictReader(fh):
                    rows[r.get("sample", "")] = r
        self.tracking = {}
        for s in self.samples:
            t = rows.get(s)
            num = lambda k: int(float(t[k])) if t and t.get(k) not in (None, "") else None
            placed = sum(a["counts"].get(s, 0) for a in self.asvs if not a.get("not_taxon"))
            assigned = placed if self.counts_are_reads else None
            raw = num("input")
            self.tracking[s] = {"raw_data": raw, "qc_filtered": num("qc_filtered"),
                                "primer_trimmed": num("primer_trimmed"), "final": num("final"),
                                "assigned": assigned,
                                "Nonchim": placed,
                                "perc": (round(100 * assigned / raw, 2) if raw else 0.0)
                                if raw is not None and assigned is not None else None,
                                "primer_hits": 0, "cutadapt_version": ""}

    def subset(self, samples) -> "RunData":
        """A view holding only these samples (run order kept): their ASVs, SV numbers
        counted again from their reads, their tracking and alpha, and BLAST / NCBI
        counts recomputed from their ASVs. The run itself is not changed."""
        import copy
        want = list(samples)
        bad = [x for x in want if x not in self.samples]
        if bad:
            raise ValueError(f"not in this run: {', '.join(bad)}")
        keep = [x for x in self.samples if x in want]
        r = copy.copy(self)
        r.samples = keep
        r.dropped = [x for x in getattr(self, "dropped", []) if x in keep]
        asvs = []
        for a in self.asvs:
            c = {x: a["counts"].get(x, 0) for x in keep}
            if sum(c.values()) > 0:
                b = dict(a)
                b["counts"], b["total"] = c, sum(c.values())
                asvs.append(b)
        asvs.sort(key=lambda a: -a["total"])
        for i, a in enumerate(asvs, 1):
            a["sv"] = f"SV{i:03d}"
        r.asvs = asvs
        r.tracking = {x: self.tracking[x] for x in keep}
        r.alpha = {x: self.alpha[x] for x in keep}
        if self.blast_info:
            bi = {k: v for k, v in self.blast_info.items() if k not in ("ncbi_batches",)}
            agg = {}
            for a in asvs:
                if a.get("agreement"):
                    agg[a["agreement"]] = agg.get(a["agreement"], 0) + 1
            bi.update({"agreement": agg, "n_asv": len(asvs),
                       "n_species": sum(1 for a in asvs if a["tax"].get("Species")),
                       "n_unassigned": sum(1 for a in asvs if a.get("rank_reached") == "none"),
                       "n_non_target": sum(1 for a in asvs if a.get("target") == "no"),
                       "n_check": sum(1 for a in asvs if a.get("check") == "yes")})
            r.blast_info = bi
        return r

    def only(self, sample: str) -> "RunData":
        return self.subset([sample])

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
