#!/usr/bin/env python3
"""shotgun_pipeline.py — shotgun metagenome profiling (beta), Illumina or ONT.

    python3 shotgun_pipeline.py --input_dir UP --output_dir OUT --platform illumina|ont
        --kraken_db DIR [--host none|human|mouse|/path/genome.fa] [--threads 8]
        [--confidence 0.1] [--bracken_level S] [--bracken_threshold 10]
        [--min_len 500 --min_q 10]   (ONT read filter)

Per sample:
  1. read QC            Illumina: fastp (adapters, quality, length ≥ 50)
                        ONT: chopper (mean Q, length), else seqkit
  2. host removal       optional: hostile (Bowtie2 for Illumina pairs,
                        minimap2 for ONT) against human T2T+HLA, mouse, or a
                        genome given as a FASTA file
  3. classification     Kraken2 against the chosen index (k-mer, every read)
  4. abundance          Bracken re-estimates species read counts from the
                        Kraken2 report, using the k-mer distribution for the
                        read length nearest the sample's
  5. resistome          (on by default) reads aligned with minimap2 to the
                        AMRFinderPlus gene catalogue; the organism carrying
                        each gene from Kraken2's per-read calls (shotgun_amr.py)
  6. function           (opt-in) MetaPhlAn 4.1.1 + HUMAnN 3.9 against full UniRef90:
                        MetaCyc pathways, KEGG KO, EC, gene panels and the genera
                        carrying them (shotgun_func.py); needs ~25 GB RAM, so the
                        backend runs such a job alone
  7. MAGs               (opt-in) assembly, binning, CheckM2, taxonomy and AMR
                        per genome bin (shotgun_mags.py)
Output — the layout the report, diversity, MaAsLin 3 and Controls already read:
  asv_table.csv   one row per species (id = NCBI taxid), estimated reads per sample
  taxonomy.csv    id, Kingdom … Species (from the Kraken2 report lineage)
  read_tracking.csv, summary.json, shotgun/<sample>/ (fastp / Kraken2 / Bracken files)
Host (Homo sapiens / Mus musculus) reads that Kraken2 still classifies are
counted in read_tracking and kept out of the community table.

Separate from every amplicon and WGS pipeline: it changes nothing they use.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import OrderedDict, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import shotgun_amr as AMR                         # noqa: E402
FQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
PAIR_RE = re.compile(r"^(?P<name>.+?)[._-](?:R)?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$", re.I)
HOST_TAXA = {9606: "Homo sapiens", 10090: "Mus musculus"}
RANK_CODES = {"D": "Kingdom", "K": "Kingdom", "P": "Phylum", "C": "Class", "O": "Order",
              "F": "Family", "G": "Genus", "S": "Species"}
RANK_ORDER = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
ENVS = ["ngamp-meta"]
LOGF: Path | None = None


def log(msg: str):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


def progress(pct: int, label: str):
    print(f"PROGRESS:{pct}|{label}", flush=True)


# ── tools ────────────────────────────────────────────────────────────────────
def conda_bases() -> list[Path]:
    out = []
    for b in (os.environ.get("CONDA_PREFIX_1"),
              os.environ.get("CONDA_EXE") and str(Path(os.environ["CONDA_EXE"]).parent.parent),
              str(Path.home() / "miniconda3"), str(Path.home() / "miniforge3"),
              str(Path.home() / "anaconda3"), str(Path.home() / "mambaforge")):
        if b and Path(b).is_dir():
            out.append(Path(b))
    return out


_TOOLS: dict = {}


def tool(name: str):
    """(path, env) of a tool from the ngamp-meta env, else PATH; None if absent."""
    if name in _TOOLS:
        return _TOOLS[name]
    found = None
    extra = os.environ.get("NGAMP_META_BIN", "")
    cands = [Path(extra) / name] if extra else []
    for base in conda_bases():
        cands += [base / "envs" / e / "bin" / name for e in ENVS]
    w = shutil.which(name)
    if w:
        cands.append(Path(w))
    for p in cands:
        if p.exists() and os.access(p, os.X_OK):
            env = os.environ.copy()
            env["PATH"] = f"{p.parent}{os.pathsep}{env.get('PATH', '')}"
            found = (str(p), env)
            break
    _TOOLS[name] = found
    return found


def run(t, args: list, logfile: Path, stdout_to: Path | None = None) -> int:
    cmd = [t[0]] + [str(a) for a in args]
    log("  $ " + " ".join([Path(t[0]).name] + [str(a) for a in args])[:500])
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with open(logfile, "a") as lf:
        lf.write(f"\n### {' '.join(cmd)}\n")
        lf.flush()
        out = open(stdout_to, "w") if stdout_to else lf
        try:
            return subprocess.run(cmd, stdout=out, stderr=lf, env=t[1]).returncode
        finally:
            if stdout_to:
                out.close()


def count_reads(fq: Path, t_seqkit) -> int:
    if t_seqkit:
        p = subprocess.run([t_seqkit[0], "stats", "-T", str(fq)], capture_output=True, text=True,
                           env=t_seqkit[1])
        try:
            return int(p.stdout.splitlines()[1].split("\t")[3].replace(",", ""))
        except Exception:
            pass
    import gzip
    op = gzip.open if str(fq).endswith(".gz") else open
    with op(fq, "rt") as fh:
        return sum(1 for _ in fh) // 4


# ── samples ──────────────────────────────────────────────────────────────────
def find_samples(d: Path, platform: str) -> "OrderedDict[str, dict]":
    """{sample: {"r1": [files], "r2": [files]}} — the UI's manifest first."""
    out: "OrderedDict[str, dict]" = OrderedDict()
    man = d / "sample_manifest.json"
    if man.exists():
        try:
            for r in json.loads(man.read_text()):
                s = str(r.get("sample", "")).strip()
                f1, f2 = str(r.get("file1", "") or ""), str(r.get("file2", "") or "")
                if not s or not f1 or not (d / f1).exists():
                    continue
                g = out.setdefault(s, {"r1": [], "r2": []})
                g["r1"].append(d / f1)
                if f2 and (d / f2).exists():
                    g["r2"].append(d / f2)
            if out:
                return out
        except Exception as e:
            log(f"[warn] sample_manifest.json unreadable ({e}) — detecting samples from file names")
    files = sorted(f for f in d.rglob("*") if (f.is_file() or f.is_symlink()) and FQ_RE.search(f.name))
    if platform == "ont":
        sys.path.insert(0, str(HERE))
        import ont_wgs_pipeline as W                  # MinKNOW chunk grouping, fastq_fail left out
        rel = [str(f.relative_to(d)) for f in files]
        folders = {r: str(Path(r).parent) if str(Path(r).parent) != "." else "" for r in rel}
        for g in W.group_fastqs(rel, folders):
            out[g["sample"]] = {"r1": [d / f for f in g["files"]], "r2": []}
        return out
    pairs = defaultdict(dict)
    for f in files:
        m = PAIR_RE.match(f.name)
        if m:
            pairs[m.group("name")][m.group("mate")] = f
        else:
            out[FQ_RE.sub("", f.name)] = {"r1": [f], "r2": []}
    for k, v in sorted(pairs.items()):
        out[k] = {"r1": [v["1"]] if "1" in v else [], "r2": [v["2"]] if "2" in v else []}
        if not out[k]["r1"]:
            out[k] = {"r1": [v["2"]], "r2": []}
    return out


def concat(files: list[Path], dest: Path) -> Path:
    """One file per mate: gzip members concatenate into a valid gzip stream."""
    if len(files) == 1:
        return files[0]
    gz = all(str(f).endswith(".gz") for f in files)
    dest = dest.with_suffix(".fastq.gz" if gz else ".fastq")
    with open(dest, "wb") as out:
        for f in files:
            with open(f, "rb") as src:
                shutil.copyfileobj(src, out)
    return dest


# ── Kraken2 report / Bracken ─────────────────────────────────────────────────
def parse_report(path: Path) -> tuple[dict, dict, int, int]:
    """Kraken-style report → ({taxid: lineage dict}, {taxid: clade reads at S},
    unclassified reads, classified reads)."""
    lineage, species, stack = {}, {}, []
    uncl = clas = 0
    for line in open(path):
        v = line.rstrip("\n").split("\t")
        if len(v) < 6:
            continue
        clade, rank, tid, name = int(v[1]), v[3].strip(), int(v[4]), v[5]
        depth = (len(name) - len(name.lstrip(" "))) // 2
        name = name.strip()
        if rank == "U":
            uncl = clade
            continue
        if rank == "R":
            clas = clade
        while stack and stack[-1][0] >= depth:
            stack.pop()
        stack.append((depth, rank, name, tid))
        if rank == "S":
            L = {}
            for _, r, n, _t in stack:
                col = RANK_CODES.get(r)
                if col and col not in L:
                    L[col] = n
            lineage[tid] = L
            species[tid] = clade
    return lineage, species, uncl, clas


def kmer_lengths(db: Path) -> list[int]:
    return sorted(int(m.group(1)) for f in db.glob("database*mers.kmer_distrib")
                  if (m := re.match(r"database(\d+)mers", f.name)))


def pick_read_len(avail: list[int], mean_len: float) -> int | None:
    if not avail:
        return None
    below = [a for a in avail if a <= mean_len]
    return max(below) if below else min(avail)


def parse_bracken(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out[int(r["taxonomy_id"])] = (r["name"], float(r["new_est_reads"]))
            except (KeyError, ValueError):
                pass
    return out


# ── main ─────────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--input_dir", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--platform", choices=["illumina", "ont"], default="illumina")
    ap.add_argument("--kraken_db", required=True)
    ap.add_argument("--host", default="none", help="none | human | mouse | path to a host FASTA")
    ap.add_argument("--threads", type=int, default=8)
    # 0, Kraken 2's own default. With the capped (8 / 16 GB) indexes and with
    # Nanopore reads, 0.05–0.1 left most reads unclassified and dropped whole
    # species (S. aureus, B. subtilis) in the mock test; Bracken's minimum-read
    # threshold is the false-positive filter instead.
    ap.add_argument("--confidence", type=float, default=0.0)
    ap.add_argument("--bracken_level", default="S", choices=["S", "G"])
    ap.add_argument("--bracken_threshold", type=int, default=10)
    ap.add_argument("--min_len", type=int, default=500, help="ONT: minimum read length")
    ap.add_argument("--min_q", type=float, default=10, help="ONT: minimum mean read quality")
    ap.add_argument("--min_abundance", type=float, default=0.01,
                    help="%% — species below this share in every sample are left out of the "
                         "tables (k-mer spill-over onto close relatives); 0 keeps everything")
    ap.add_argument("--memory_mapping", action="store_true",
                    help="Kraken2 reads the index from disk instead of RAM (slower, less memory)")
    ap.add_argument("--db_paths", default=str(HERE / "databases" / "db_paths.json"))
    ap.add_argument("--no_amr", action="store_true", help="skip the read-based resistome")
    ap.add_argument("--amr_min_identity", type=float, default=0,
                    help="%% identity of a read to an AMR gene (default 90 Illumina, 85 Nanopore)")
    ap.add_argument("--amr_min_breadth", type=float, default=60,
                    help="%% of a gene's length the reads must cover for it to be called present")
    ap.add_argument("--func", action="store_true",
                    help="function profiles with MetaPhlAn + HUMAnN (opt-in, ~25 GB RAM, hours per sample)")
    ap.add_argument("--func_panels", default="gut,food,environment",
                    help="gene panels to report (data/function_panels.tsv): any of gut,food,environment")
    ap.add_argument("--metadata", default="", help="metadata.csv (sample groups for PERMANOVA)")
    ap.add_argument("--mags", action="store_true", help="assemble and bin genomes (slow, opt-in)")
    ap.add_argument("--mag_min_contig", type=int, default=1500)
    ap.add_argument("--mag_polish", action="store_true",
                    help="Nanopore MAGs: polish the contigs with Medaka before binning (slow, memory-hungry)")
    a = ap.parse_args(argv)

    inp, out = Path(a.input_dir), Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "shotgun"
    work.mkdir(exist_ok=True)
    logs = out / "logs"
    logs.mkdir(exist_ok=True)
    db = Path(a.kraken_db)

    log("=" * 70)
    log(f"NextGen-Amplicon shotgun metagenome pipeline (beta) — {a.platform}")
    log("=" * 70)
    k2, br = tool("kraken2"), tool("bracken")
    fastp, chopper, seqkit, hostile = tool("fastp"), tool("chopper"), tool("seqkit"), tool("hostile")
    if not k2:
        log("[ERROR] kraken2 not found — run: bash ~/r16s-app/setup_shotgun.sh")
        return 2
    if not (db / "hash.k2d").exists():
        log(f"[ERROR] {db} is not a Kraken2 index (no hash.k2d) — run: bash ~/r16s-app/setup_shotgun.sh")
        return 2
    vers = {}
    vflag = {"bracken": ["-v"], "seqkit": ["version"]}
    for n, t in (("kraken2", k2), ("bracken", br), ("fastp", fastp), ("chopper", chopper),
                 ("seqkit", seqkit), ("hostile", hostile)):
        if t:
            p = subprocess.run([t[0]] + vflag.get(n, ["--version"]), capture_output=True, text=True,
                               env=t[1])
            m = re.search(r"(?:version|v|\s)\s*(\d+\.\d+(?:\.\d+)?)", " " + p.stdout + p.stderr)
            vers[n] = m.group(1) if m else "?"
    dbp = {}
    if a.db_paths and Path(a.db_paths).exists():
        try:
            dbp = json.loads(Path(a.db_paths).read_text())
        except Exception as e:
            log(f"[warn] db_paths.json unreadable: {e}")
    mm2 = tool("minimap2")
    amr_db = amr_cat = amr_bin = None
    if not a.no_amr:
        amr_db, amr_bin = AMR.find_install(dbp, conda_bases())
        if amr_bin:
            p = subprocess.run([amr_bin, "--version"], capture_output=True, text=True)
            vers["amrfinder"] = (p.stdout or p.stderr).strip().splitlines()[-1] if (p.stdout or p.stderr) else "?"
        if not mm2 or not amr_db:
            log("[warn] resistome skipped — " + ("minimap2 not found" if not mm2 else
                "AMRFinderPlus database not found (installed by setup_wgs.sh; or set "
                "db_paths.json key amrfinder_db)"))
        else:
            amr_cat = AMR.load_catalog(amr_db)
            p = subprocess.run([mm2[0], "--version"], capture_output=True, text=True, env=mm2[1])
            vers["minimap2"] = p.stdout.strip() or "?"
    if not a.amr_min_identity:
        a.amr_min_identity = 90.0 if a.platform == "illumina" else 85.0
    log("tools: " + ", ".join(f"{k} {v}" for k, v in vers.items()))
    if amr_cat:
        log(f"resistome: AMRFinderPlus catalogue {AMR.db_version(amr_db)} ({len(amr_cat):,} genes), "
            f"identity >= {a.amr_min_identity:g}%, gene breadth >= {a.amr_min_breadth:g}%")
    if a.host != "none" and not hostile:
        log("[warn] hostile not installed — host reads are NOT removed before classification; "
            "Kraken2 still counts them as host and they are kept out of the table")

    samples = find_samples(inp, a.platform)
    if not samples:
        log("[ERROR] no FASTQ files found")
        return 2
    log(f"{len(samples)} sample(s): " + ", ".join(samples))
    if not a.memory_mapping:
        try:
            mem_avail = next(int(ln.split()[1]) * 1024 for ln in open("/proc/meminfo")
                             if ln.startswith("MemAvailable:"))
            need = (db / "hash.k2d").stat().st_size
            if need > 0.9 * mem_avail:
                a.memory_mapping = True
                log(f"[note] the index ({need / 1e9:.1f} GB) does not fit in free memory "
                    f"({mem_avail / 1e9:.1f} GB) — reading it from disk (--memory-mapping, slower)")
        except Exception:
            pass
    if a.confidence > 0:
        log(f"[note] Kraken 2 confidence {a.confidence}: fewer reads are classified; with a capped "
            "index or Nanopore reads this can drop species entirely (0 is recommended)")
    avail = kmer_lengths(db)
    if not br or not avail:
        log("[warn] Bracken " + ("not installed" if not br else "k-mer distributions missing in the index")
            + " — species reads are Kraken2's own clade counts (less accurate for close species)")

    track, per_sample, lineages, used_len = [], {}, {}, {}
    amr_rows, qc_reads = {}, {}
    n = len(samples)
    for i, (s, fs) in enumerate(samples.items()):
        sd = work / s
        sd.mkdir(exist_ok=True)
        base = 5 + int(85 * i / n)
        step = max(1, int(85 / n))
        r1 = concat(fs["r1"], sd / "merged_R1")
        r2 = concat(fs["r2"], sd / "merged_R2") if fs["r2"] else None
        paired = a.platform == "illumina" and r2 is not None
        row = {"sample": s, "input": 0}

        # 1. QC
        progress(base, f"[{s}] Read QC")
        mean_len = 150.0
        if a.platform == "illumina" and fastp:
            q1 = sd / "qc_R1.fastq.gz"
            q2 = sd / "qc_R2.fastq.gz"
            args = ["-i", r1, "-o", q1, "-w", min(16, a.threads), "-l", 50,
                    "-j", sd / "fastp.json", "-h", sd / "fastp.html"]
            if paired:
                args += ["-I", r2, "-O", q2, "--detect_adapter_for_pe"]
            if run(fastp, args, logs / "fastp.log") != 0:
                log(f"[ERROR] [{s}] fastp failed — see logs/fastp.log")
                return 1
            fj = json.loads((sd / "fastp.json").read_text())
            row["input"] = fj["summary"]["before_filtering"]["total_reads"] // (2 if paired else 1)
            row["after_qc"] = fj["summary"]["after_filtering"]["total_reads"] // (2 if paired else 1)
            mean_len = float(fj["summary"]["after_filtering"].get("read1_mean_length") or 150)
            r1, r2 = q1, (q2 if paired else None)
        elif a.platform == "ont":
            row["input"] = count_reads(r1, seqkit)
            q1 = sd / "qc.fastq.gz"
            if chopper:
                import gzip
                op = gzip.open if str(r1).endswith(".gz") else open
                with op(r1, "rb") as src, gzip.open(q1, "wb", compresslevel=3) as dst:
                    p = subprocess.Popen([chopper[0], "-q", str(a.min_q), "-l", str(a.min_len),
                                          "-t", str(min(8, a.threads))],
                                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=open(logs / "chopper.log", "a"), env=chopper[1])
                    import threading
                    th = threading.Thread(target=lambda: (shutil.copyfileobj(src, p.stdin), p.stdin.close()))
                    th.start()
                    shutil.copyfileobj(p.stdout, dst)
                    th.join()
                    p.wait()
            elif seqkit:
                run(seqkit, ["seq", "-m", a.min_len, "-Q", a.min_q, "-o", q1, r1], logs / "seqkit.log")
            else:
                q1 = r1
            r1 = q1
            row["after_qc"] = count_reads(r1, seqkit)
            if seqkit:
                p = subprocess.run([seqkit[0], "stats", "-T", str(r1)], capture_output=True, text=True,
                                   env=seqkit[1])
                try:
                    mean_len = float(p.stdout.splitlines()[1].split("\t")[6].replace(",", ""))
                except Exception:
                    mean_len = 3000.0
            else:
                mean_len = 3000.0
        else:
            row["input"] = row["after_qc"] = count_reads(r1, seqkit)

        # 2. host removal
        if a.host != "none" and hostile:
            progress(base + step // 4, f"[{s}] Removing host reads ({a.host})")
            idx = {"human": "human-t2t-hla", "mouse": "mouse-mm39"}.get(a.host, a.host)
            hd = sd / "host_removed"
            hd.mkdir(exist_ok=True)
            args = ["clean", "--fastq1", r1] + (["--fastq2", r2] if paired else []) + \
                   ["--index", idx, "--output", hd, "--threads", a.threads, "--force"]
            custom = a.host not in ("human", "mouse")
            if a.platform == "illumina":
                # standard indexes ship for Bowtie2; a host given as a FASTA file is
                # indexed on the fly by minimap2 (hostile builds no Bowtie2 index)
                args += ["--aligner", "minimap2" if custom else "bowtie2"]
            js = sd / "hostile.json"
            if run(hostile, args, logs / "hostile.log", stdout_to=js) == 0:
                cl = sorted(hd.glob("*.clean*.fastq.gz"))
                c1 = [f for f in cl if "clean_2" not in f.name]
                c2 = [f for f in cl if "clean_2" in f.name]
                if c1:
                    r1 = c1[0]
                    r2 = c2[0] if (paired and c2) else r2
                try:
                    h = json.loads(js.read_text())
                    h = h[0] if isinstance(h, list) else h
                    row["host_removed"] = int(h.get("reads_removed", 0)) // (2 if paired else 1)
                except Exception:
                    pass
            else:
                log(f"[warn] [{s}] hostile failed — see logs/hostile.log; continuing without host removal")

        qc_reads[s] = {"r1": r1, "r2": r2 if paired else None}

        # 3a. resistome: reads against the AMR gene catalogue (before Kraken 2, so
        #     that Kraken 2's per-read output can be kept for just these reads)
        gene_reads = {}
        if amr_cat:
            progress(base + step // 3, f"[{s}] AMR genes in reads")
            paf = sd / "amr.paf"
            if AMR.align(mm2, AMR._cds_file(amr_db), [x for x in (r1, r2) if x], a.platform,
                         a.threads, paf, logs / "minimap2_amr.log", run):
                # reads (a pair counts once, as in the gene counts)
                n_reads = max(1, row.get("after_qc", 0) - row.get("host_removed", 0))
                amr_rows[s], gene_reads = AMR.summarise(paf, amr_cat, a.platform, n_reads,
                                                        a.amr_min_identity / 100, a.amr_min_breadth)
                n_amr = sum(1 for x in amr_rows[s] if x["type"] == "AMR")
                log(f"[{s}] resistome: {n_amr} AMR gene famil{'y' if n_amr == 1 else 'ies'}"
                    f" ({len(amr_rows[s]) - n_amr} stress / virulence)")
            else:
                log(f"[warn] [{s}] minimap2 against the AMR catalogue failed — see logs/minimap2_amr.log")
        wanted = set().union(*gene_reads.values()) if gene_reads else set()

        # 3. Kraken2
        progress(base + step // 2, f"[{s}] Kraken2 classification")
        rep = sd / "kraken2.report"
        args = ["--db", db, "--threads", a.threads, "--confidence", a.confidence,
                "--report", rep, "--output", "-"]
        if a.memory_mapping:
            args.append("--memory-mapping")
        args += (["--paired", r1, r2] if paired else [r1])
        read_tax = {}
        if wanted:
            # "--output -" suppresses Kraken 2's per-read lines; stdout keeps them
            args[args.index("--output") + 1] = "/dev/stdout"
            cmd = [k2[0]] + [str(x) for x in args]
            log("  $ kraken2 " + " ".join(str(x) for x in args)[:500])
            with open(logs / "kraken2.log", "a") as lf:
                lf.write(f"\n### {' '.join(cmd)}\n")
                lf.flush()
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=lf, env=k2[1], text=True,
                                     bufsize=1 << 20)
                read_tax = AMR.kraken_keep(p.stdout, wanted)
                k2rc = p.wait()
        else:
            k2rc = run(k2, args, logs / "kraken2.log", stdout_to=Path(os.devnull))
        if k2rc != 0 or not rep.exists():
            log(f"[ERROR] [{s}] Kraken2 failed — see logs/kraken2.log")
            return 1
        if amr_rows.get(s):
            AMR.assign_hosts(amr_rows[s], gene_reads, read_tax, AMR.report_tree(rep))
        lin, kspecies, uncl, clas = parse_report(rep)
        lineages.update(lin)
        row["classified"] = clas
        row["unclassified"] = uncl

        # 4. Bracken
        est = None
        rl = pick_read_len(avail, mean_len) if br else None
        if rl:
            progress(base + 3 * step // 4, f"[{s}] Bracken abundance (read length {rl})")
            bo, bw = sd / "bracken_species.tsv", sd / "bracken_species.report"
            if run(br, ["-d", db, "-i", rep, "-o", bo, "-w", bw, "-r", rl, "-l", a.bracken_level,
                        "-t", a.bracken_threshold], logs / "bracken.log") == 0 and bo.exists():
                est = parse_bracken(bo)
                used_len[s] = rl
                if bw.exists():
                    lineages.update(parse_report(bw)[0])
            else:
                log(f"[warn] [{s}] Bracken failed — using Kraken2 clade counts")
        if est is None:
            est = {t: ("", float(c)) for t, c in kspecies.items()}
        host = {t: v for t, v in est.items() if t in HOST_TAXA}
        row["host_classified"] = int(round(sum(v[1] for v in host.values())))
        per_sample[s] = {t: v[1] for t, v in est.items() if t not in HOST_TAXA and v[1] > 0}
        row["final"] = int(round(sum(per_sample[s].values())))
        track.append(row)
        log(f"[{s}] {row.get('input', 0):,} reads → {row.get('after_qc', 0):,} after QC → "
            f"{clas:,} classified ({row['final']:,} to species, {row['host_classified']:,} host)")

    # ── tables in the amplicon layout ──────────────────────────────────────────
    progress(92, "Writing tables")
    names = list(samples)
    taxa = sorted({t for d in per_sample.values() for t in d},
                  key=lambda t: -sum(per_sample[s].get(t, 0) for s in names))
    # Bracken hands a few reads to dozens of close relatives of every abundant
    # species (shared k-mers): a pure B. subtilis culture came out as 90 species.
    # Keep a species when it reaches min_abundance % in at least one sample.
    dropped, dropped_reads = 0, {s: 0.0 for s in names}
    if a.min_abundance > 0:
        tot = {s: sum(per_sample[s].values()) or 1.0 for s in names}
        keep = [t for t in taxa if any(100 * per_sample[s].get(t, 0) / tot[s] >= a.min_abundance
                                       for s in names)]
        dropped = len(taxa) - len(keep)
        for t in set(taxa) - set(keep):
            for s in names:
                dropped_reads[s] += per_sample[s].get(t, 0)
        taxa = keep
        for r in track:
            r["below_min_abundance"] = int(round(dropped_reads[r["sample"]]))
            r["final"] = int(round(r["final"] - dropped_reads[r["sample"]]))
        if dropped:
            log(f"{dropped} species below {a.min_abundance}% in every sample left out "
                f"({', '.join(f'{s}: {int(dropped_reads[s]):,} reads' for s in names)})")
    with open(out / "asv_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tax_id"] + names)
        for t in taxa:
            w.writerow([t] + [int(round(per_sample[s].get(t, 0))) for s in names])
    with open(out / "taxonomy.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tax_id"] + RANK_ORDER)
        for t in taxa:
            L = lineages.get(t, {})
            w.writerow([t] + [L.get(r, "") for r in RANK_ORDER])
    cols = ["sample", "input", "after_qc", "host_removed", "classified", "unclassified",
            "host_classified", "below_min_abundance", "final"]
    with open(out / "read_tracking.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in track:
            w.writerow(r)
    dom = defaultdict(float)
    for s in names:
        for t in taxa:
            dom[lineages.get(t, {}).get("Kingdom", "Unknown")] += per_sample[s].get(t, 0)
    meta = {}
    for f in ("inspect.txt", "library_report.tsv"):
        if (db / f).exists():
            meta["index_file"] = f
    amr_info = {}
    if amr_cat:
        n_fam = AMR.write_tables(out, {s: amr_rows.get(s, []) for s in names})
        amr_info = {"db": str(amr_db), "db_version": AMR.db_version(amr_db),
                    "min_identity": a.amr_min_identity, "min_breadth": a.amr_min_breadth,
                    "gene_families": n_fam,
                    "amr_families": len({r["gene"] for v in amr_rows.values() for r in v
                                         if r["type"] == "AMR"})}
        log(f"resistome: {amr_info['amr_families']} AMR gene famil"
            f"{'y' if amr_info['amr_families'] == 1 else 'ies'} across {len(names)} sample(s) — "
            "amr_genes.csv / amr_rpkm.csv written")
    func_info = {}
    if a.func:
        import shotgun_func as FN
        progress(80, "Function: MetaPhlAn + HUMAnN")
        try:
            func_info = FN.run_all(samples=names, qc=qc_reads, platform=a.platform, out=out, threads=a.threads,
                                   dbp=dbp, panels=[x.strip() for x in a.func_panels.split(",") if x.strip()],
                                   tool=tool, run=run, log=log, progress=progress, conda_bases=conda_bases())
            for k, v in (func_info.get("versions") or {}).items():
                vers[k] = v
            if not func_info.get("skipped"):
                func_info["stats"] = FN.func_stats(out, a.metadata)
        except Exception as e:                     # function profiles are an extra: never fail the profile
            import traceback
            log(f"[warn] function step failed: {e}")
            (logs / "func_traceback.txt").write_text(traceback.format_exc())
            (out / "func").mkdir(exist_ok=True)
            (out / "func" / "SKIPPED.txt").write_text(f"the function step stopped with an error: {e}\n")
            func_info = {**func_info, "error": str(e)} if isinstance(func_info, dict) else {"error": str(e)}
    mag_info = {}
    if a.mags:
        import shotgun_mags as MAG
        progress(94, "MAGs: assembly and binning")
        try:
            mag_info = MAG.run_all(samples=names, qc=qc_reads, platform=a.platform, out=out,
                                   threads=a.threads, dbp=dbp, kraken=(k2, db, a.memory_mapping),
                                   amr_db=amr_db, amr_bin=amr_bin, min_contig=a.mag_min_contig, polish=a.mag_polish,
                                   tool=tool, run=run,
                                   log=log, progress=progress, conda_bases=conda_bases())
        except Exception as e:                     # MAGs are an extra: never fail the profile
            import traceback
            log(f"[warn] MAG step failed: {e}")
            (logs / "mags_traceback.txt").write_text(traceback.format_exc())
            mag_info = {"error": str(e)}
    (out / "summary.json").write_text(json.dumps({
        "pipeline": "shotgun", "platform": a.platform, "beta": True,
        "kraken_db": str(db), "kraken_db_name": db.name, "confidence": a.confidence,
        "bracken_level": a.bracken_level, "bracken_read_length": used_len,
        "bracken_threshold": a.bracken_threshold, "host": a.host, "tools": vers,
        "min_abundance_pct": a.min_abundance, "species_below_min_abundance": dropped,
        "emu_counts": "estimated_reads", "n_samples": len(names), "n_taxa": len(taxa),
        "domains": {k: round(v) for k, v in sorted(dom.items(), key=lambda kv: -kv[1])},
        "amr": amr_info, "func": func_info, "mags": mag_info},
        indent=1))
    progress(100, "Done")
    log(f"{len(taxa)} species across {len(names)} sample(s) — asv_table.csv / taxonomy.csv written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
