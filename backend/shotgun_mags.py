#!/usr/bin/env python3
"""shotgun_mags.py — metagenome-assembled genomes (beta, opt-in), per sample.

Called by shotgun_pipeline.py after profiling when the run asks for MAGs:

  1. assembly        Illumina: MEGAHIT · Nanopore: Flye --meta (+ Medaka, opt-in)
  2. depth           reads mapped back with minimap2 → samtools →
                     jgi_summarize_bam_contig_depths
  3. binning         MetaBAT 2 (contigs >= 1.5 kb)
  4. taxonomy        Kraken 2 on the binned contigs (NCBI names, the same as
                     the read profile; majority by length)
  5. quality         CheckM2 completeness / contamination → MIMAG-style tier
  6. refinement      bins > 10 % contaminated are split by contig genus, and
                     pieces of one species at similar depth are joined — each
                     kept only when CheckM2 scores the result better
  7. GTDB            sourmash, when the WGS GTDB database is installed
  8. AMR per genome  AMRFinderPlus on the contigs, each hit tied to its bin

Tools come from their own conda envs: ngamp-mags (megahit, metabat2 —
setup_shotgun.sh --mags), ngamp-wgs (flye, sourmash, amrfinder), checkm2.
A missing tool skips its step with a note; nothing here can fail the profile.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import shotgun_amr as AMR

ENV_PREF = {
    "megahit": ["ngamp-mags", "ngamp-meta", "megahit"],
    "metabat2": ["ngamp-mags", "ngamp-meta", "metabat2"],
    "jgi_summarize_bam_contig_depths": ["ngamp-mags", "ngamp-meta", "metabat2"],
    "flye": ["ngamp-wgs", "wgs", "flye"],
    "checkm2": ["checkm2", "ngamp-wgs", "wgs"],
    "sourmash": ["ngamp-wgs", "wgs", "sourmash"],
    "amrfinder": ["ngamp-amr", "ngamp-wgs", "wgs", "amrfinder"],
    "minimap2": ["ngamp-meta", "ngamp-mags", "ngamp-wgs", "ngamp-phylo"],
    "samtools": ["ngamp-mags", "ngamp-meta", "ngamp-wgs", "ngamp-phylo"],
    "seqkit": ["ngamp-meta", "ngamp-wgs"],
}
_FOUND: dict = {}


def find(name: str, bases: list[Path]):
    """(path, env) like shotgun_pipeline.tool(), searching the env that ships it."""
    if name in _FOUND:
        return _FOUND[name]
    cands = [b / "envs" / e / "bin" / name for b in bases for e in ENV_PREF.get(name, [])]
    w = shutil.which(name)
    if w:
        cands.append(Path(w))
    hit = None
    for p in cands:
        if p.exists() and os.access(p, os.X_OK):
            env = os.environ.copy()
            env["PATH"] = f"{p.parent}{os.pathsep}{env.get('PATH', '')}"
            env.setdefault("OMP_NUM_THREADS", "4")
            hit = (str(p), env)
            break
    _FOUND[name] = hit
    return hit


def read_fasta(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def n50(lengths: list[int]) -> int:
    tot, acc = sum(lengths), 0
    for L in sorted(lengths, reverse=True):
        acc += L
        if acc >= tot / 2:
            return L
    return 0


def tier(comp, cont) -> str:
    if comp is None or cont is None:
        return "—"
    if comp > 90 and cont < 5:
        return "High"
    if comp >= 50 and cont < 10:
        return "Medium"
    return "Low"


def find_dbs(dbp: dict) -> dict:
    """checkm2 / GTDB databases: db_paths.json first, then the WGS folder."""
    out = {k: dbp.get(k, "") for k in ("checkm2_db", "gtdb_sourmash", "gtdb_lineages")}
    wgs = Path(__file__).resolve().parent / "databases" / "wgs"
    if wgs.is_dir():
        if not (out["checkm2_db"] and Path(out["checkm2_db"]).exists()):
            dm = sorted(wgs.rglob("uniref100.KO.1.dmnd"))
            out["checkm2_db"] = str(dm[0]) if dm else ""
        if not (out["gtdb_sourmash"] and Path(out["gtdb_sourmash"]).exists()):
            sig = sorted(wgs.glob("gtdb-*reps.k31*.sig.zip"), key=lambda x: -x.stat().st_size)
            out["gtdb_sourmash"] = str(sig[0]) if sig else ""
        if not (out["gtdb_lineages"] and Path(out["gtdb_lineages"]).exists()):
            lin = sorted(wgs.glob("gtdb-*reps.lineages.csv*"))
            out["gtdb_lineages"] = str(lin[0]) if lin else ""
    return {k: v for k, v in out.items() if v and Path(v).exists()}


def mean_q(seqkit, fq: Path) -> float:
    if not seqkit:
        return 0.0
    p = subprocess.run([seqkit[0], "stats", "-a", "-T", str(fq)], capture_output=True, text=True,
                       env=seqkit[1])
    try:
        h, v = p.stdout.splitlines()[:2]
        return float(dict(zip(h.split("\t"), v.split("\t"))).get("AvgQual", 0))
    except Exception:
        return 0.0


# ── one sample ───────────────────────────────────────────────────────────────
def checkm2(t, in_dir: Path, out_dir: Path, threads, db, s: str, run, logs: Path) -> dict:
    """{sample.bin: (completeness, contamination)} for the .fa files in in_dir."""
    q = out_dir / "quality_report.tsv"
    if not q.exists():
        run(t, ["predict", "--input", in_dir, "--output-directory", out_dir, "-x", "fa",
                "--threads", threads, "--database_path", db, "--force"], logs / "checkm2.log")
    out = {}
    if q.exists():
        with open(q) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    out[f"{s}.{r['Name']}"] = (float(r["Completeness"]), float(r["Contamination"]))
                except (KeyError, ValueError):
                    pass
    return out


def bin_species(contigs: list[str], lens: dict, ctax: dict) -> str:
    """The Kraken 2 species holding most of a bin's classified length (>= 50 %), else ''."""
    sv, tot = defaultdict(int), 0
    for c in contigs:
        if c in ctax:
            tot += lens.get(c, 0)
            if ctax[c][1]:
                sv[ctax[c][1]] += lens.get(c, 0)
    if not sv:
        return ""
    sp, v = max(sv.items(), key=lambda x: x[1])
    return sp if v >= 0.5 * tot else ""


def split_by_genus(contigs: list[str], lens: dict, ctax: dict, min_bp: int = 300_000,
                   min_share: float = 0.15) -> list[tuple[str, list[str]]]:
    """[(genus, contigs)] when a bin's contigs fall into two or more genera that each hold
    >= min_bp and >= min_share of the bin; contigs of no / other genera join the largest."""
    tot = sum(lens.get(c, 0) for c in contigs) or 1
    gbp = defaultdict(int)
    for c in contigs:
        if c in ctax:
            gbp[ctax[c][0]] += lens.get(c, 0)
    big = [g for g, v in sorted(gbp.items(), key=lambda x: -x[1]) if v >= min_bp and v / tot >= min_share]
    if len(big) < 2:
        return []
    parts = {g: [] for g in big}
    for c in contigs:
        g = ctax.get(c, ("",))[0]
        parts[g if g in parts else big[0]].append(c)
    return [(g, parts[g]) for g in big]


def one_sample(s: str, r1: Path, r2: Path | None, platform: str, wd: Path, md: Path,
               threads: int, T: dict, dbs: dict, kraken, amr_db, min_contig: int,
               run, log, step=lambda label: None, polish: bool = False) -> dict:
    res = {"sample": s, "bins": [], "amr": [], "notes": []}
    # work left by an earlier run is reused only for the very same reads and settings
    sig = json.dumps([platform, min_contig] + [[str(x), x.stat().st_size] for x in (r1, r2) if x]
                     + (["medaka"] if platform == "ont" and polish else []))
    sigf = wd / "inputs.sig"
    if wd.exists() and (not sigf.exists() or sigf.read_text() != sig):
        shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True, exist_ok=True)
    sigf.write_text(sig)
    logs = wd / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    # 1. assembly
    contigs = wd / "contigs.fasta"
    if contigs.exists() and contigs.stat().st_size > 0:
        log(f"[{s}] MAGs: reusing the assembly from the last run")
    elif platform == "illumina":
        step("assembling with MEGAHIT")
        if not T.get("megahit"):
            res["notes"].append("MEGAHIT not installed (setup_shotgun.sh --mags)")
            return res
        ad = wd / "megahit"
        shutil.rmtree(ad, ignore_errors=True)
        args = (["-1", r1, "-2", r2] if r2 else ["-r", r1]) + \
            ["-o", ad, "-t", threads, "--min-contig-len", 1000, "-m", 0.5]
        if run(T["megahit"], args, logs / "megahit.log") != 0 or not (ad / "final.contigs.fa").exists():
            res["notes"].append("MEGAHIT failed — see logs/megahit.log")
            return res
        shutil.move(str(ad / "final.contigs.fa"), contigs)
        shutil.rmtree(ad, ignore_errors=True)
        res["assembler"] = "MEGAHIT"
    else:
        if not T.get("flye"):
            res["notes"].append("Flye not installed (setup_wgs.sh)")
            return res
        q = mean_q(T.get("seqkit"), r1)
        mode = "nano-hq" if q >= 15 else "nano-raw"
        step(f"assembling with Flye --meta --{mode}")
        res["flye_mode"] = mode
        ad = wd / "flye"
        shutil.rmtree(ad, ignore_errors=True)
        if run(T["flye"], ["--meta", f"--{mode}", r1, "--out-dir", ad, "--threads", threads],
               logs / "flye.log") != 0 or not (ad / "assembly.fasta").exists():
            res["notes"].append("Flye failed — see logs/flye.log")
            return res
        shutil.move(str(ad / "assembly.fasta"), contigs)
        info = ad / "assembly_info.txt"
        if info.exists():
            shutil.copy(info, wd / "assembly_info.txt")
        shutil.rmtree(ad, ignore_errors=True)
    res.setdefault("assembler", "MEGAHIT" if platform == "illumina" else "Flye --meta")

    # 1b. Nanopore: Medaka. Unpolished Nanopore contigs carry small indels that break
    #     genes, and CheckM2 (a protein-based estimate) then reads near-complete genomes
    #     as 40-70 % complete (Zymo mock: B. subtilis 3.9 Mb in 2 contigs at 58 %).
    #     Opt-in: at 8 threads Medaka took the whole WSL down on nextapp (28 GB), so it
    #     runs at 4 threads and only when asked for.
    if platform == "ont" and polish:
        pol = wd / "polish.txt"
        if pol.exists():
            res["polish"] = pol.read_text().strip()
        else:
            import ont_wgs_pipeline as W
            if W.find_tool("medaka_consensus"):
                step("polishing the contigs with Medaka")
                cons, label = W.run_medaka(r1, contigs, wd / "medaka", min(4, threads), "auto", True,
                                           logs / "medaka.log")
                if cons:
                    shutil.move(str(contigs), str(wd / "contigs_flye.fasta"))
                    shutil.move(str(cons), str(contigs))
                    res["polish"] = f"Medaka ({label})"
                    # anything computed on the unpolished contigs is stale now
                    for x in ("depth.txt", "amrfinder.tsv"):
                        (wd / x).unlink(missing_ok=True)
                    for x in ("bins", "bins_unsplit", "bins_merged", "checkm2", "checkm2_split",
                              "checkm2_merge", "sourmash"):
                        shutil.rmtree(wd / x, ignore_errors=True)
                    (wd / "merged.json").unlink(missing_ok=True)
                else:
                    res["polish"] = f"none (Medaka {label})"
                shutil.rmtree(wd / "medaka", ignore_errors=True)
            else:
                res["polish"] = "none (Medaka not installed)"
                res["notes"].append("Medaka not installed — Nanopore contigs not polished (setup_wgs.sh)")
            pol.write_text(res["polish"])
    lens = {n: len(q) for n, q in read_fasta(contigs)}
    res["assembly"] = {"contigs": len(lens), "bp": sum(lens.values()), "n50": n50(list(lens.values())),
                       "longest": max(lens.values()) if lens else 0,
                       "contigs_ge_min": sum(1 for v in lens.values() if v >= min_contig)}
    log(f"[{s}] MAGs: assembly {len(lens):,} contigs, {sum(lens.values()) / 1e6:.1f} Mb, "
        f"N50 {res['assembly']['n50']:,} bp")

    # 2. depth
    depth_txt = wd / "depth.txt"
    mm2, sam = T.get("minimap2"), T.get("samtools")
    jgi = T.get("jgi_summarize_bam_contig_depths")
    depth = {}
    if not depth_txt.exists() and mm2 and sam and jgi:
        step("mapping reads back for contig depth")
        bam = wd / "reads.bam"
        preset = "sr" if platform == "illumina" else "map-ont"
        cmd1 = [mm2[0], "-ax", preset, "-t", str(threads), str(contigs), str(r1)] + \
            ([str(r2)] if r2 else [])
        cmd2 = [sam[0], "sort", "-@", str(max(1, threads // 2)), "-m", "1G", "-o", str(bam), "-"]
        log(f"  $ minimap2 -ax {preset} contigs reads | samtools sort")
        with open(logs / "map.log", "a") as lf:
            p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=lf, env=mm2[1])
            p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stderr=lf, env=sam[1])
            p1.stdout.close()
            p2.wait()
            p1.wait()
        if p1.returncode == 0 and p2.returncode == 0:
            # jgi counts only reads >= 97 % identical to the contig by default: right for
            # Illumina, but it drops nearly every Nanopore read (90–95 %) and leaves all
            # depths at zero, so MetaBAT 2 forms no bins
            pid = ["--percentIdentity", 97 if platform == "illumina" else 80]
            run(jgi, pid + ["--outputDepth", depth_txt, bam], logs / "depth.log")
        bam.unlink(missing_ok=True)
    if depth_txt.exists():
        with open(depth_txt) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    depth[r["contigName"]] = float(r["totalAvgDepth"])
                except (KeyError, ValueError):
                    pass

    # 3. binning
    bins_dir = wd / "bins"
    bins = {}
    if T.get("metabat2") and depth_txt.exists():
        if not any(bins_dir.glob("*.fa")):
            step("binning with MetaBAT 2")
            bins_dir.mkdir(exist_ok=True)
            run(T["metabat2"], ["-i", contigs, "-a", depth_txt, "-o", bins_dir / "bin",
                                "-m", max(1500, min_contig), "-t", threads, "--seed", 1],
                logs / "metabat2.log")
        for f in sorted(bins_dir.glob("bin.*.fa"), key=lambda x: int(re.sub(r"\D", "", x.stem) or 0)):
            bins[f"{s}.{f.stem}"] = [n for n, _ in read_fasta(f)]
    else:
        res["notes"].append("MetaBAT 2 " + ("not installed (setup_shotgun.sh --mags)" if not T.get("metabat2")
                                            else "skipped: read depth could not be computed"))
    contig_bin = {c: b for b, cs in bins.items() for c in cs}
    log(f"[{s}] MAGs: {len(bins)} bin(s)")
    if not bins and depth_txt.exists():
        if depth and not any(v > 0 for v in depth.values()):
            res["notes"].append("no bins: every contig has read depth 0 (see logs/depth.log)")
        else:
            try:
                tail = (logs / "metabat2.log").read_text(errors="replace").strip().splitlines()[-1]
            except Exception:
                tail = ""
            res["notes"].append("MetaBAT 2 formed no bins" + (f" ({tail[:160]})" if tail else ""))

    # 4. Kraken 2 on the binned contigs: each contig's genus / species, used to name
    #    the bins and to split bins that mix two organisms (step 6)
    k2, kdb, mmap = kraken
    ctax = {}                                      # contig -> (genus, species)
    if bins and k2:
        step("naming bins with Kraken 2")
        allb = wd / "binned_contigs.fa"
        with open(allb, "w") as fh:
            for f in bins_dir.glob("bin.*.fa"):
                for n, q in read_fasta(f):
                    fh.write(f">{n}\n{q}\n")
        rep, kout = wd / "bins_kraken2.report", wd / "bins_kraken2.out"
        args = ["--db", kdb, "--threads", threads, "--report", rep, "--output", kout]
        if mmap:
            args.append("--memory-mapping")
        if run(k2, args + [allb], logs / "kraken2_bins.log") == 0 and rep.exists():
            tree = AMR.report_tree(rep)
            for line in open(kout):
                v = line.split("\t")
                if len(v) < 4 or v[0] != "C":
                    continue
                g = AMR._up_to(tree, int(v[2]), "G")
                if g:
                    ctax[v[1]] = (g, AMR._up_to(tree, int(v[2]), "S") or "")
        allb.unlink(missing_ok=True)
        kout.unlink(missing_ok=True)

    # 5. CheckM2
    ck = {}
    can_ck = bool(bins and T.get("checkm2") and dbs.get("checkm2_db"))
    if can_ck:
        step(f"CheckM2 on {len(bins)} bins")
        ck = checkm2(T["checkm2"], bins_dir, wd / "checkm2", threads, dbs["checkm2_db"], s, run, logs)
    elif bins:
        res["notes"].append("CheckM2 " + ("not installed" if not T.get("checkm2") else
                                          "database not configured (checkm2_db)") +
                            " — bin quality not assessed")

    # 6. Split bins that hold two organisms. With one sample MetaBAT 2 has only
    #    composition and depth to go on, so organisms of similar GC at similar depth
    #    end up together (Zymo even mock: Salmonella + E. coli, Listeria + E. faecalis,
    #    each 75-83 % "contaminated"). A bin CheckM2 calls > 10 % contaminated is cut
    #    by the Kraken 2 genus of its contigs; the cut is kept only when CheckM2 scores
    #    the parts better (completeness - 5 x contamination, summed).
    refined = {}
    if can_ck and ctax and not (wd / "bins_unsplit").exists():
        cands = {}
        for b, cs in bins.items():
            comp, cont = ck.get(b, (None, None))
            if cont is None or cont <= 10:
                continue
            parts = split_by_genus(cs, lens, ctax)
            if len(parts) >= 2:
                cands[b] = parts
        if cands:
            step(f"splitting {len(cands)} mixed bin(s) by contig taxonomy")
            sd_ = wd / "split_try"
            shutil.rmtree(sd_, ignore_errors=True)
            sd_.mkdir()
            seqs = {}
            for b in cands:
                seqs.update(dict(read_fasta(bins_dir / (b[len(s) + 1:] + ".fa"))))
            for b, parts in cands.items():
                stem = b[len(s) + 1:]
                for i, (g, cs) in enumerate(parts):
                    with open(sd_ / f"{stem}{chr(97 + i)}.fa", "w") as fh:
                        for c in cs:
                            fh.write(f">{c}\n{seqs[c]}\n")
            ck2 = checkm2(T["checkm2"], sd_, wd / "checkm2_split", threads, dbs["checkm2_db"], s, run, logs)
            score = lambda cc: max(0.0, cc[0] - 5 * cc[1]) if cc else 0.0
            keep = []
            for b, parts in cands.items():
                stem = b[len(s) + 1:]
                names = [f"{s}.{stem}{chr(97 + i)}" for i in range(len(parts))]
                if all(n in ck2 for n in names) and sum(score(ck2[n]) for n in names) > score(ck[b]):
                    keep.append((b, stem, names, parts))
            if keep:
                (wd / "bins_unsplit").mkdir(exist_ok=True)
            for b, stem, names, parts in keep:
                shutil.move(str(bins_dir / f"{stem}.fa"), str(wd / "bins_unsplit" / f"{stem}.fa"))
                del bins[b]
                ck.pop(b, None)
                for n, (g, cs) in zip(names, parts):
                    shutil.move(str(sd_ / f"{n[len(s) + 1:]}.fa"), str(bins_dir / f"{n[len(s) + 1:]}.fa"))
                    bins[n] = cs
                    ck[n] = ck2[n]
                    refined[n] = stem
                log(f"[{s}] MAGs: {stem} held {len(parts)} organisms "
                    f"({', '.join(g for g, _ in parts)}) — split by contig taxonomy")
            shutil.rmtree(sd_, ignore_errors=True)
    elif (wd / "bins_unsplit").exists():           # split in an earlier run of this folder
        for f in bins_dir.glob("bin.*[a-z].fa"):
            if not f.stem.endswith("m"):
                refined[f"{s}.{f.stem}"] = re.sub(r"[a-z]$", "", f.stem)
        if can_ck:
            ck.update(checkm2(T["checkm2"], wd / "split_try", wd / "checkm2_split", threads,
                              dbs["checkm2_db"], s, run, logs))

    # 6b. ... and join bins that are pieces of one genome: same Kraken 2 species, similar
    #     depth (within 1.5x, so two strains at different abundance stay apart), kept
    #     only when CheckM2 finds the joined bin >= 10 points more complete and still
    #     < 10 % contaminated (Zymo Nanopore: E. faecalis came out as 45 % + 39 %).
    merged = {}
    mj = wd / "merged.json"
    if can_ck and ctax and not mj.exists():
        groups = defaultdict(list)
        for b, cs in bins.items():
            sp = bin_species(cs, lens, ctax)
            if sp:
                groups[sp].append(b)
        dep = {b: sum(depth.get(c, 0) * lens.get(c, 0) for c in cs) / max(1, sum(lens.get(c, 0) for c in cs))
               for b, cs in bins.items()}
        cands = {}
        for sp, bs in groups.items():
            if len(bs) < 2:
                continue
            bs.sort(key=lambda b: -sum(lens.get(c, 0) for c in bins[b]))
            ref = dep[bs[0]] or 1e-9
            take = [b for b in bs if dep[b] and max(dep[b], ref) / min(dep[b], ref) <= 1.5]
            if len(take) >= 2:
                cands[f"{s}.{take[0][len(s) + 1:]}m"] = (sp, take)
        if cands:
            step(f"joining fragments of {len(cands)} genome(s)")
            md_ = wd / "merge_try"
            shutil.rmtree(md_, ignore_errors=True)
            md_.mkdir()
            for n, (sp, take) in cands.items():
                with open(md_ / f"{n[len(s) + 1:]}.fa", "w") as fh:
                    for b in take:
                        for c, q in read_fasta(bins_dir / f"{b[len(s) + 1:]}.fa"):
                            fh.write(f">{c}\n{q}\n")
            ck3 = checkm2(T["checkm2"], md_, wd / "checkm2_merge", threads, dbs["checkm2_db"], s, run, logs)
            for n, (sp, take) in cands.items():
                if n not in ck3 or any(b not in ck for b in take):
                    continue
                comp, cont = ck3[n]
                if comp >= max(ck[b][0] for b in take) + 10 and cont < 10 \
                        and cont <= max(ck[b][1] for b in take) + 5:
                    (wd / "bins_merged").mkdir(exist_ok=True)
                    for b in take:
                        shutil.move(str(bins_dir / f"{b[len(s) + 1:]}.fa"),
                                    str(wd / "bins_merged" / f"{b[len(s) + 1:]}.fa"))
                        bins.pop(b, None)
                        ck.pop(b, None)
                        refined.pop(b, None)
                    shutil.move(str(md_ / f"{n[len(s) + 1:]}.fa"), str(bins_dir / f"{n[len(s) + 1:]}.fa"))
                    bins[n] = [c for c, _ in read_fasta(bins_dir / f"{n[len(s) + 1:]}.fa")]
                    ck[n] = ck3[n]
                    merged[n] = " + ".join(b[len(s) + 1:] for b in take)
                    log(f"[{s}] MAGs: {merged[n]} are one {sp} genome — joined as {n[len(s) + 1:]} "
                        f"({comp:.0f}% complete, {cont:.1f}% contamination)")
            shutil.rmtree(md_, ignore_errors=True)
        mj.write_text(json.dumps(merged))
    elif mj.exists():
        merged = json.loads(mj.read_text())
        if can_ck and merged:
            ck.update(checkm2(T["checkm2"], wd / "merge_try", wd / "checkm2_merge", threads,
                              dbs["checkm2_db"], s, run, logs))
    contig_bin = {c: b for b, cs in bins.items() for c in cs}

    # bin names from the contig taxonomy (majority of the bin's length)
    ktax = {}
    for b, cs in bins.items():
        tot = sum(lens.get(c, 0) for c in cs) or 1
        gv, sv = defaultdict(float), defaultdict(float)
        for c in cs:
            if c in ctax:
                g, sp = ctax[c]
                gv[g] += lens.get(c, 0)
                if sp:
                    sv[sp] += lens.get(c, 0)
        if gv:
            g, gl = max(gv.items(), key=lambda x: x[1])
            name = g
            sp = [(k, v) for k, v in sv.items() if k.split()[0] == g]
            if sp:
                k, v = max(sp, key=lambda x: x[1])
                if v >= 0.6 * gl:
                    name = k
            ktax[b] = (name, round(100 * gl / tot, 1))

    # 7. GTDB (sourmash multigather), when the WGS database is installed
    gtdb = {}
    sm = T.get("sourmash")
    if bins and sm and dbs.get("gtdb_sourmash"):
        sd = wd / "sourmash"
        sd.mkdir(exist_ok=True)
        sigs = []
        for f in bins_dir.glob("bin.*.fa"):
            sig = sd / f"{f.stem}.sig.zip"
            if not sig.exists():
                run(sm, ["sketch", "dna", "-p", "k=31,scaled=1000,abund", "--name", f"{s}.{f.stem}",
                         "-o", sig, f], logs / "sourmash.log")
            if sig.exists():
                sigs.append(sig)
        if sigs and not any(sd.glob("*.csv")):
            step("GTDB naming (sourmash)")
            # multigather reads the whole GTDB database once per bin (half an hour for
            # 8 bins on nextapp). One prefetch pass with all bins merged keeps just the
            # genomes that share k-mers with any bin; the per-bin gather then runs on those.
            gdb = dbs["gtdb_sourmash"]
            allsig, matches = sd / "all_bins.sig.zip", sd / "gtdb_matches.zip"
            if run(sm, ["sig", "merge"] + sigs + ["-k", 31, "-o", allsig], logs / "sourmash.log") == 0 \
                    and run(sm, ["prefetch", allsig, gdb, "-k", 31, "--threshold-bp", 50000,
                                 "--save-matches", matches], logs / "sourmash.log") == 0 \
                    and matches.exists() and matches.stat().st_size > 0:
                gdb = matches
            run(sm, ["multigather", "--query"] + sigs + ["--db", gdb, "-k", 31,
                     "--threshold-bp", 50000, "--output-dir", sd], logs / "sourmash.log")
        rows = []
        for f in sd.glob("*.csv"):
            if f.name.endswith(".prefetch.csv"):
                continue
            with open(f, newline="") as fh:
                rows += list(csv.DictReader(fh))
        best = {}
        for r in rows:
            q = r.get("query_name", "")
            try:
                fq = float(r.get("f_unique_weighted") or r.get("f_unique_to_query") or 0)
            except ValueError:
                fq = 0
            if q and (q not in best or fq > best[q][0]):
                best[q] = (fq, r)
        import ont_wgs_pipeline as W
        accs = {re.sub(r"^(RS_|GB_)", "", (r.get("name") or "").split(" ")[0]) for _, r in best.values()}
        lin = W.load_lineages(Path(dbs["gtdb_lineages"]), accs) if dbs.get("gtdb_lineages") else {}
        for q, (fq, r) in best.items():
            acc = re.sub(r"^(RS_|GB_)", "", (r.get("name") or "").split(" ")[0])
            try:
                ani = float(r.get("match_containment_ani") or r.get("average_containment_ani") or 0)
            except ValueError:
                ani = 0
            if not ani:
                # multigather leaves the ANI columns empty: estimate it from the share of
                # the bin's k-mers found in the match (containment ANI = C^(1/k))
                try:
                    ani = float(r.get("f_orig_query") or 0) ** (1 / 31)
                except ValueError:
                    ani = 0
            ani = ani * 100 if ani <= 1 else ani
            L = lin.get(acc, {})
            gtdb[q] = {"gtdb_species": L.get("species", ""), "gtdb_genus": L.get("genus", ""),
                       "gtdb_accession": acc, "gtdb_ani": round(ani, 2),
                       "gtdb_call": W.species_call(ani)}

    # 8. AMRFinderPlus on the contigs
    amrf = T.get("amrfinder")
    if amrf and amr_db:
        at = wd / "amrfinder.tsv"
        if not at.exists():
            step("AMRFinderPlus on the contigs")
            big = wd / "contigs_1kb.fasta"
            with open(big, "w") as fh:
                for n, q in read_fasta(contigs):
                    if len(q) >= 1000:
                        fh.write(f">{n}\n{q}\n")
            run(amrf, ["-n", big, "--plus", "--threads", threads, "--database", amr_db, "-o", at],
                logs / "amrfinder.log")
            big.unlink(missing_ok=True)
        if at.exists():
            with open(at) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    g = lambda *k: next((r[x] for x in k if x in r and r[x] not in (None, "")), "")
                    c = g("Contig id")
                    res["amr"].append({
                        "sample": s, "bin": contig_bin.get(c, ""), "contig": c,
                        "gene": g("Element symbol", "Gene symbol"), "type": g("Type", "Element type"),
                        "subtype": g("Subtype", "Element subtype"), "class": g("Class"),
                        "subclass": g("Subclass"), "name": g("Element name", "Sequence name"),
                        "method": g("Method"), "identity": g("% Identity to reference",
                                                             "% Identity to reference sequence"),
                        "coverage": g("% Coverage of reference", "% Coverage of reference sequence")})

    # bin table
    tot_cov = sum(depth.get(c, 0) * L for c, L in lens.items() if L >= min_contig) or 1
    md.mkdir(parents=True, exist_ok=True)
    for f in sorted(bins_dir.glob("bin.*.fa"), key=lambda x: int(re.sub(r"\D", "", x.stem) or 0)):
        b = f"{s}.{f.stem}"
        cs = bins.get(b, [])
        L = [lens.get(c, 0) for c in cs]
        seq = "".join(q for _, q in read_fasta(f))
        gc = 100 * sum(seq.count(x) for x in "GCgc") / max(1, len(seq))
        cov = sum(depth.get(c, 0) * lens.get(c, 0) for c in cs)
        comp, cont = ck.get(b, (None, None))
        kt = ktax.get(b, ("", ""))
        gt = gtdb.get(b, {})
        shutil.copy(f, md / f"{b}.fa")
        res["bins"].append({
            "sample": s, "bin": b, "size_bp": sum(L), "contigs": len(cs), "n50": n50(L),
            "gc": round(gc, 1), "mean_depth": round(cov / max(1, sum(L)), 1),
            "share_of_assembly_coverage": round(100 * cov / tot_cov, 1),
            "completeness": comp, "contamination": cont, "quality": tier(comp, cont),
            "kraken2_taxon": kt[0], "kraken2_bp_share": kt[1],
            "gtdb_species": gt.get("gtdb_species", ""), "gtdb_ani": gt.get("gtdb_ani", ""),
            "gtdb_call": gt.get("gtdb_call", ""),
            "split_from": refined.get(b, ""),
            "merged_from": merged.get(b, ""),
            "amr_genes": ";".join(sorted({a["gene"] for a in res["amr"]
                                          if a["bin"] == b and a["type"] == "AMR"}))})
    res["bins"].sort(key=lambda r: ({"High": 0, "Medium": 1, "Low": 2}.get(r["quality"], 3),
                                    -(r["completeness"] or 0)))
    if contigs.exists():
        with open(contigs, "rb") as src, gzip.open(md / f"{s}.contigs.fasta.gz", "wb", 3) as dst:
            shutil.copyfileobj(src, dst)
    return res


BIN_COLS = ["sample", "bin", "quality", "completeness", "contamination", "size_bp", "contigs", "n50",
            "gc", "mean_depth", "share_of_assembly_coverage", "kraken2_taxon", "kraken2_bp_share",
            "gtdb_species", "gtdb_ani", "gtdb_call", "amr_genes", "split_from", "merged_from"]
AMR_COLS = ["sample", "bin", "contig", "gene", "type", "subtype", "class", "subclass", "name",
            "method", "identity", "coverage"]


def run_all(samples, qc, platform, out: Path, threads, dbp, kraken, amr_db, min_contig, tool, run,
            log, progress, conda_bases, polish: bool = False, amr_bin: str | None = None) -> dict:
    bases = list(conda_bases)
    T = {n: find(n, bases) for n in ENV_PREF}
    T["minimap2"] = T["minimap2"] or tool("minimap2")
    T["samtools"] = T["samtools"] or tool("samtools")
    T["seqkit"] = T["seqkit"] or tool("seqkit")
    if amr_bin:                                    # the amrfinder that matches amr_db's format
        env = os.environ.copy()
        env["PATH"] = f"{Path(amr_bin).parent}{os.pathsep}{env.get('PATH', '')}"
        T["amrfinder"] = (amr_bin, env)
    dbs = find_dbs(dbp)
    vers = {}
    used = ["megahit" if platform == "illumina" else "flye", "metabat2"] + \
        (["checkm2"] if dbs.get("checkm2_db") else []) + (["sourmash"] if dbs.get("gtdb_sourmash") else [])
    for n, flag in (("megahit", "--version"), ("metabat2", "--help"), ("flye", "--version"),
                    ("checkm2", "--version"), ("sourmash", "--version")):
        if T.get(n) and n in used:
            p = subprocess.run([T[n][0], flag], capture_output=True, text=True, env=T[n][1])
            m = re.search(r"(\d+\.\d+(?:\.\d+)?)", p.stdout + p.stderr)
            vers[n] = m.group(1) if m else "?"
    need = "megahit" if platform == "illumina" else "flye"
    log("MAGs: tools " + ", ".join(f"{k} {v}" for k, v in vers.items())
        + ("" if T.get(need) else f" — {need} NOT installed"))
    mags_root = out / "mags"
    results = []
    for i, s in enumerate(samples):
        progress(94 + int(5 * i / max(1, len(samples))), f"[{s}] MAGs: assembly and binning")
        q = qc.get(s) or {}
        if not q.get("r1"):
            continue
        wd = out / "shotgun" / s / "mags_work"
        wd.mkdir(parents=True, exist_ok=True)
        pct = 94 + int(5 * i / max(1, len(samples)))
        results.append(one_sample(s, Path(q["r1"]), Path(q["r2"]) if q.get("r2") else None, platform,
                                  wd, mags_root / s, threads, T, dbs, kraken, amr_db, min_contig,
                                  run, log, step=lambda label, s=s, pct=pct:
                                  progress(pct, f"[{s}] MAGs: {label}"), polish=polish))
    bins = [b for r in results for b in r["bins"]]
    with open(out / "mags_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BIN_COLS, extrasaction="ignore")
        w.writeheader()
        for b in bins:
            w.writerow(b)
    with open(out / "mags_amr.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=AMR_COLS, extrasaction="ignore")
        w.writeheader()
        for r in results:
            for a in r["amr"]:
                w.writerow(a)
    info = {"tools": vers, "assembler": "MEGAHIT" if platform == "illumina" else "Flye --meta",
            "binner": "MetaBAT 2" if T.get("metabat2") else "",
            "checkm2": bool(T.get("checkm2") and dbs.get("checkm2_db")),
            "gtdb": bool(T.get("sourmash") and dbs.get("gtdb_sourmash")),
            "min_contig": min_contig,
            "samples": {r["sample"]: {"assembly": r.get("assembly", {}), "notes": r["notes"],
                                      "polish": r.get("polish", ""),
                                      "flye_mode": r.get("flye_mode", ""), "bins": len(r["bins"])}
                        for r in results},
            "n_bins": len(bins),
            "n_high": sum(1 for b in bins if b["quality"] == "High"),
            "n_medium": sum(1 for b in bins if b["quality"] == "Medium"),
            "n_split": len({b["split_from"] for b in bins if b.get("split_from")}),
            "n_merged": sum(1 for b in bins if b.get("merged_from"))}
    for r in results:
        for n in r["notes"]:
            log(f"[warn] [{r['sample']}] MAGs: {n}")
    log(f"MAGs: {len(bins)} bin(s) — {info['n_high']} high, {info['n_medium']} medium quality — "
        "mags_summary.csv / mags_amr.csv, genomes in mags/")
    (out / "mags" ).mkdir(exist_ok=True)
    (out / "mags" / "mags_info.json").write_text(json.dumps(info, indent=1))
    return info
