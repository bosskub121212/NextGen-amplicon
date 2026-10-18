#!/usr/bin/env python3
"""Species-level comparison of 16S references on reads from genomes NOT in them.

Held-out genomes (complete RefSeq assemblies, explicitly NOT type material) of a
set of species give their own 16S copies. From those:
  * ONT 16S: full-length amplicon reads (27F..1492R) with ~1 % errors in both
    orientations, mixed in uneven proportions, classified by Emu against each
    Emu database;
  * PacBio: every distinct primer-trimmed copy taken as a perfect ASV, named by
    DADA2 assignTaxonomy + addSpecies (the app's tax_helpers.R) against each
    DADA2 reference.
All references are restricted to the same genera, so no database wins by
having fewer competitors. Species are compared by NCBI taxid (renames such as
Lactobacillus plantarum -> Lactiplantibacillus plantarum keep their taxid),
falling back to the name.

    python3 eval_16s_multicopy.py --held held/copies.json \
        --emu NAME=DIR [--emu NAME=DIR ...] --dada2 NAME=TRAINSET,SPECIES [...] \
        --genera G1,G2,... --out results/
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

FL_F, FL_R = "AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT"
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]",
         "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
         "V": "[ACG]", "N": "[ACGT]"}
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def rc(s):
    return s.translate(COMP)[::-1]


def rx(p):
    return re.compile("".join(IUPAC[c] for c in p))


def amplicon(seq: str, trim: bool) -> str | None:
    """27F..1492R region; with trim=True the primers are cut off (cutadapt-like).
    Up to 2 mismatches would be cutadapt's job; exact IUPAC match is enough for
    the 16S genes of these genera, and a copy without both sites is skipped."""
    f = rx(FL_F).search(seq)
    r = rx(rc(FL_R)).search(seq, f.end() if f else 0)
    if not f or not r:
        return None
    return seq[f.end():r.start()] if trim else seq[f.start():r.end()]


def mutate(s: str, rng: random.Random, sub=0.004, ins=0.003, dele=0.003) -> str:
    out = []
    for ch in s:
        x = rng.random()
        if x < dele:
            continue
        if x < dele + sub:
            out.append(rng.choice([b for b in "ACGT" if b != ch]))
        else:
            out.append(ch)
        if rng.random() < ins:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def restrict_emu(src: Path, dst: Path, genera: set[str]) -> int:
    dst.mkdir(parents=True, exist_ok=True)
    keep = set()
    with open(src / "taxonomy.tsv") as fh, open(dst / "taxonomy.tsv", "w") as out:
        head = fh.readline()
        out.write(head)
        cols = head.rstrip("\n").split("\t")
        gi = cols.index("genus")
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) > gi and v[gi] in genera:
                keep.add(v[0])
                out.write(line)
    n = 0
    with open(src / "species_taxid.fasta") as fh, open(dst / "species_taxid.fasta", "w") as out:
        w = False
        for line in fh:
            if line.startswith(">"):
                w = line[1:].split(":")[0] in keep
                n += w
            if w:
                out.write(line)
    return n


def restrict_dada2(train: Path, species: Path, dst: Path, genera: set[str]) -> tuple[Path, Path]:
    dst.mkdir(parents=True, exist_ok=True)
    stem = train.name.replace(".fa.gz", "")
    t_out = dst / f"{stem}.fa.gz"
    s_out = dst / species.name
    with gzip.open(train, "rt") as fh, gzip.open(t_out, "wt") as out:
        w = False
        for line in fh:
            if line.startswith(">"):
                ranks = [x for x in line[1:].strip().split(";") if x]
                w = any(r in genera for r in ranks[-2:])
            if w:
                out.write(line)
    with gzip.open(species, "rt") as fh, gzip.open(s_out, "wt") as out:
        w = False
        for line in fh:
            if line.startswith(">"):
                p = line[1:].split()
                w = len(p) >= 2 and p[1] in genera
            if w:
                out.write(line)
    return t_out, s_out


def bray(a: dict, b: dict) -> float:
    keys = set(a) | set(b)
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys) / max(1e-9, sum(a.values()) + sum(b.values()))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--held", required=True)
    ap.add_argument("--emu", action="append", default=[])
    ap.add_argument("--dada2", action="append", default=[])
    ap.add_argument("--genera", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reads", type=int, default=8000)
    ap.add_argument("--emu-bin", default="emu")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--samples", default="S1,S2", help="which simulated ONT samples to classify")
    ap.add_argument("--tax-helpers", default=str(Path(__file__).resolve().parents[2] /
                                                 "backend" / "r_scripts" / "tax_helpers.R"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    genera = {g for g in a.genera.split(",") if g}
    held = json.loads(Path(a.held).read_text())
    by_sp = defaultdict(list)
    for acc, v in held.items():
        by_sp[(v["species"], v["taxid"])].append(v["copies"])
    species = sorted(by_sp)
    res = {"species": [s[0] for s in species], "ont": {}, "pacbio": {}}

    # ── ONT 16S reads: two samples, each from a different held-out genome ─────
    rng = random.Random(7)
    weights = [rng.choice([1, 2, 4, 8, 16]) for _ in species]
    samples = {}
    for si in (0, 1):
        reads, truth = [], Counter()
        pool = []
        for (nm, tid), w in zip(species, weights):
            gens = by_sp[(nm, tid)]
            copies = [amplicon(c, False) for c in gens[si % len(gens)]]
            copies = [c for c in copies if c]
            if copies:
                pool.append((nm, tid, w, copies))
        tot = sum(p[2] for p in pool)
        for nm, tid, w, copies in pool:
            n = round(a.reads * w / tot)
            truth[tid] += n
            for _ in range(n):
                s = mutate(rng.choice(copies), rng)
                reads.append(rc(s) if rng.random() < 0.5 else s)
        rng.shuffle(reads)
        fq = out / f"ont_S{si + 1}.fastq"
        with open(fq, "w") as fh:
            for i, s in enumerate(reads):
                fh.write(f"@r{i}\n{s}\n+\n{'5' * len(s)}\n")
        samples[f"S{si + 1}"] = (fq, {k: v / sum(truth.values()) for k, v in truth.items()})
    names_by_tid = {tid: nm for nm, tid in species}

    for spec in a.emu:
        label, db = spec.split("=", 1)
        rdb = out / "db" / f"emu_{label}"
        nseq = restrict_emu(Path(db), rdb, genera)
        res["ont"][label] = {"db_sequences": nseq, "samples": {}}
        for sname, (fq, truth) in samples.items():
            if sname not in a.samples.split(","):
                continue
            od = out / f"emu_{label}_{sname}"
            if not any(od.glob("*rel-abundance.tsv")):      # an earlier run's result is reused
                subprocess.run([a.emu_bin, "abundance", str(fq), "--db", str(rdb), "--type", "map-ont",
                                "--threads", str(a.threads), "--output-dir", str(od),
                                "--min-abundance", "0.0001"], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            tsv = next(od.glob("*rel-abundance.tsv"))
            rows = [ln.rstrip("\n").split("\t") for ln in open(tsv)]
            hd = rows[0]
            obs, obs_name = Counter(), {}
            for r in rows[1:]:
                d = dict(zip(hd, r))
                try:
                    tid = int(float(d.get("tax_id", "0")))
                except ValueError:
                    continue                      # "unassigned" rows
                ab = float(d.get("abundance", 0) or 0)
                nm = d.get("species", "")
                # same species under a different taxid (merged / renamed taxa): use the name
                hit = next((t for t, n in names_by_tid.items() if n == nm), None)
                key = tid if tid in truth else (hit if hit in truth else tid)
                obs[key] += ab
                obs_name[key] = nm
            correct = [t for t in truth if obs.get(t, 0) >= 0.001]
            wrong = {t: v for t, v in obs.items() if t not in truth and v >= 0.001}
            res["ont"][label]["samples"][sname] = {
                "species_found": len(correct), "species_total": len(truth),
                "missing": [names_by_tid[t] for t in truth if t not in correct],
                "false_species": len(wrong), "false_pct": round(100 * sum(wrong.values()), 2),
                "false_list": sorted(((obs_name.get(t, str(t)), round(100 * v, 2)) for t, v in wrong.items()),
                                     key=lambda x: -x[1])[:12],
                "bray_curtis_species": round(bray(truth, dict(obs)), 3),
                "per_species": {names_by_tid[t]: [round(100 * truth[t], 2), round(100 * obs.get(t, 0), 2)]
                                for t in truth}}
            print(f"ONT {label} {sname}: {len(correct)}/{len(truth)} species, "
                  f"{len(wrong)} false ({100 * sum(wrong.values()):.1f}%), "
                  f"BC {res['ont'][label]['samples'][sname]['bray_curtis_species']}", flush=True)

    # ── PacBio: every distinct primer-trimmed copy as a perfect ASV ───────────
    asvs = {}
    for (nm, tid) in species:
        for gen in by_sp[(nm, tid)]:
            for c in gen:
                t = amplicon(c, True)
                if t:
                    asvs.setdefault(t, (nm, tid))
    fa = out / "pacbio_asvs.fasta"
    fa.write_text("".join(f">a{i}\n{s}\n" for i, s in enumerate(asvs)))
    truth_tab = out / "pacbio_truth.tsv"
    truth_tab.write_text("id\tspecies\n" + "".join(f"a{i}\t{v[0]}\n" for i, v in enumerate(asvs.values())))
    for spec in a.dada2:
        label, files = spec.split("=", 1)
        tr, sp = files.split(",")
        rtr, rsp = restrict_dada2(Path(tr), Path(sp), out / "db" / f"dada2_{label}", genera)
        rcsv = out / f"pacbio_{label}.csv"
        r = f'''
suppressPackageStartupMessages(library(dada2))
source("{a.tax_helpers}")
s <- as.character(Biostrings::readDNAStringSet("{fa}"))
lv <- tax_pick_levels("{rtr}")
tax <- if (is.null(lv)) assignTaxonomy(s, "{rtr}", minBoot=50, multithread={a.threads}) else assignTaxonomy(s, "{rtr}", minBoot=50, taxLevels=lv, multithread={a.threads})
tax <- if (exists("tax_add_species", mode="function")) tax_add_species(tax, "{rsp}") else addSpecies(tax, "{rsp}")
d <- data.frame(id=names(s), Genus=tax[, "Genus"], Species=tax[, "Species"])
write.csv(d, "{rcsv}", row.names=FALSE)
'''
        rp = subprocess.run(["Rscript", "-e", r], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if rp.returncode != 0:
            print(f"PacBio {label}: R failed\n{rp.stderr[-1500:]}", flush=True)
            continue
        rows = [ln.strip().split(",") for ln in open(rcsv)][1:]
        truth = {f"a{i}": v[0] for i, v in enumerate(asvs.values())}
        g_ok = s_ok = s_wrong = s_none = 0
        wrong_ex = Counter()
        for rid, g, s in ((x[0].strip('"'), x[1].strip('"'), x[2].strip('"')) for x in rows):
            tg, te = truth[rid].split(" ", 1)
            g_ok += g == tg
            if s in ("NA", ""):
                s_none += 1
            elif g == tg and te in s.split("/"):
                s_ok += 1
            else:
                s_wrong += 1
                wrong_ex[f"{truth[rid]} -> {g} {s}"] += 1
        n = len(rows)
        res["pacbio"][label] = {"asvs": n, "genus_correct_pct": round(100 * g_ok / n, 1),
                                "species_correct_pct": round(100 * s_ok / n, 1),
                                "species_wrong_pct": round(100 * s_wrong / n, 1),
                                "species_unassigned_pct": round(100 * s_none / n, 1),
                                "wrong_examples": wrong_ex.most_common(10)}
        print(f"PacBio {label}: genus {100 * g_ok / n:.1f}%, species correct {100 * s_ok / n:.1f}%, "
              f"wrong {100 * s_wrong / n:.1f}%, unassigned {100 * s_none / n:.1f}% (n={n})", flush=True)
    (out / "results.json").write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
