"""wgs_package.py — one ONT-WGS sample as the lab delivers it.

    build_package(result_dir, sample, out_dir) -> manifest

out_dir/<sample>/ gets eight numbered folders (the lab's layout) holding that sample's
result files, plus tables the app derives from them:

    1.QcReads            read_statistics.tsv (before / after filtering)
    2.Assembly           assembly FASTA + graph + Flye info, assembly_statistics.tsv (QUAST-style)
    3.Completeness       CheckM2 quality report
    4.Taxonomy           identification.tsv, Contig-based Taxonomy/, 16S/, MLST/, Phylogeny/,
                         Read-based Taxonomy/ (readtax.py), TYGS/ (tygs.py) when present
    5.ANI                type-strain ANIb / ANIm / TETRA, skani to the closest genome, ANI matrix
    6.Annotation         Bakta outputs
    7.AMR                AMRFinderPlus, VFDB, PlasmidFinder, MOB-suite
    8.readsMappingToRef  reads / assembly vs the closest reference genome

Tool logs and pipeline internals never go in. A folder whose step produced nothing gets
NOT_AVAILABLE.txt. Large files are hard-linked when possible, never read into memory.
"""
from __future__ import annotations

import csv
import json
import os
import re
import shutil
from pathlib import Path

FOLDERS = ["1.QcReads", "2.Assembly", "3.Completeness", "4.Taxonomy", "5.ANI", "6.Annotation", "7.AMR",
           "8.readsMappingToRef"]


def safe(name: str) -> str:
    # \w misses Thai vowel / tone marks (ก้อง and ก่อง would both become ก_อง): keep the Thai block
    return re.sub(r"[^\w\u0E00-\u0E7F.\-]+", "_", name).strip("_.") or "sample"


def _fasta_lengths(path: Path) -> tuple[list[int], int]:
    lens, gc, cur = [], 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
                continue
            s = line.strip().upper()
            cur += len(s)
            gc += s.count("G") + s.count("C")
    if cur:
        lens.append(cur)
    return lens, gc


def assembly_stats(fasta: Path) -> list[tuple[str, object]]:
    """QUAST-style 'statistics without reference' of an assembly FASTA."""
    lens, gc = _fasta_lengths(Path(fasta))
    lens.sort(reverse=True)
    tot = sum(lens)
    rows: list[tuple[str, object]] = [("# contigs", len(lens))]
    for t in (0, 1000, 5000, 10000, 25000, 50000):
        rows.append((f"# contigs (>= {t} bp)", sum(1 for L in lens if L >= t)))
    rows += [("Largest contig", lens[0] if lens else 0), ("Total length", tot)]
    for t in (0, 1000, 5000, 10000, 25000, 50000):
        rows.append((f"Total length (>= {t} bp)", sum(L for L in lens if L >= t)))

    def nx(f):
        acc = 0
        for i, L in enumerate(lens, 1):
            acc += L
            if acc >= f * tot:
                return L, i
        return 0, 0
    n50, l50 = nx(0.5)
    n90, l90 = nx(0.9)
    aun = round(sum(L * L for L in lens) / tot) if tot else 0
    rows += [("N50", n50), ("N90", n90), ("auN", aun), ("L50", l50), ("L90", l90),
             ("GC (%)", round(100 * gc / tot, 2) if tot else 0.0)]
    return rows


def _put(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _tsv(path: Path, head: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(head)
        w.writerows(rows)


def _is_internal(p: Path) -> bool:
    n = p.name
    return n.endswith(".log") or n.startswith(".done") or n.startswith("log_")


_STATE_FILES = {"state.json", "tygs.json"}


def _copy_tree(src: Path, dst: Path) -> None:
    for f in sorted(src.rglob("*")):
        rel = f.relative_to(src)
        if any(p.startswith(".") for p in rel.parts) or f.name in _STATE_FILES:
            continue                       # raw TYGS answers (other genomes of the job), app state
        if f.is_file() and not _is_internal(f):
            _put(f, dst / f.relative_to(src))


def build_package(result_dir: Path, sample: str, out_dir: Path, name: str | None = None) -> dict:
    """The sample's eight-folder package; returns {"dir", "folders", "files", "missing"}."""
    result_dir, out_dir = Path(result_dir), Path(out_dir)
    sd = result_dir / "samples" / sample
    if not sd.is_dir():
        raise FileNotFoundError(f"no results for sample {sample!r}")
    try:
        res = json.loads((sd / "result.json").read_text())
    except Exception:
        res = {}
    S = sample
    root = out_dir / (name or safe(S))
    if root.exists():
        shutil.rmtree(root)
    F = {f: root / f for f in FOLDERS}
    for p in F.values():
        p.mkdir(parents=True)

    def take(rel: str, folder: str, dest: str | None = None) -> bool:
        p = sd / rel
        if p.is_file():
            _put(p, F[folder] / (dest or p.name))
            return True
        return False

    # 1. reads
    raw, flt = res.get("reads_raw") or {}, res.get("reads_filtered") or {}
    if raw:
        pct = lambda a, b: round(100 * a / b, 2) if a is not None and b else ""
        _tsv(F["1.QcReads"] / "read_statistics.tsv", ["Metric", "Before filtering", "After filtering", "% kept"],
             [["Reads", raw.get("reads", ""), flt.get("reads", ""), pct(flt.get("reads"), raw.get("reads"))],
              ["Total bases", raw.get("bases", ""), flt.get("bases", ""), pct(flt.get("bases"), raw.get("bases"))],
              ["Mean read length", raw.get("mean_len", ""), flt.get("mean_len", ""), ""],
              ["Median read length", raw.get("median_len", ""), flt.get("median_len", ""), ""],
              ["Read N50", raw.get("n50", ""), flt.get("n50", ""), ""],
              ["Mean read quality (Q)", raw.get("mean_q", ""), flt.get("mean_q", ""), ""]])
    # 2. assembly
    asm = sd / f"{S}_assembly.fasta"
    if take(f"{S}_assembly.fasta", "2.Assembly"):
        _tsv(F["2.Assembly"] / "assembly_statistics.tsv", ["Statistics without reference", f"{S}_assembly"],
             [[k, v] for k, v in assembly_stats(asm)])
    take(f"{S}_assembly_graph.gfa", "2.Assembly")
    take("flye_assembly_info.txt", "2.Assembly")
    if res.get("contigs"):
        keys = ["contig", "length", "gc", "depth", "circular", "type"]
        _tsv(F["2.Assembly"] / "contigs.tsv", keys, [[c.get(k, "") for k in keys] for c in res["contigs"]])
    # 3. completeness
    take("checkm2/quality_report.tsv", "3.Completeness", "checkm2_quality_report.tsv")
    # 4. taxonomy
    idn = res.get("identification") or {}
    if idn:
        _tsv(F["4.Taxonomy"] / "identification.tsv", ["Strain", "Identification result", "Call", "Method", "ANI (%)",
                                                       "Reference", "Accession"],
             [[S, idn.get("species", ""), idn.get("call", ""), idn.get("method", ""), idn.get("ani", "") or "",
               idn.get("reference", ""), idn.get("accession", "")]])
    take("sourmash_gather.csv", "4.Taxonomy", "Contig-based Taxonomy/sourmash_gather_gtdb.csv")
    take(f"{S}_16S.fasta", "4.Taxonomy", f"16S/{S}_16S.fasta")
    take("barrnap.gff", "4.Taxonomy", "16S/barrnap.gff")
    take("taxonomy/type_strain_16S.csv", "4.Taxonomy", "16S/type_strain_16S.csv")
    take("mlst.tsv", "4.Taxonomy", "MLST/mlst.tsv")
    for f in sorted((sd / "taxonomy").glob("*_tree.nwk")) + sorted((sd / "taxonomy").glob("*_alignment.fasta")):
        _put(f, F["4.Taxonomy"] / "Phylogeny" / f.name)
    for sub, dest in (("readtax", "Read-based Taxonomy"), ("tygs", "TYGS")):
        if (sd / sub).is_dir():
            _copy_tree(sd / sub, F["4.Taxonomy"] / dest)
    # 5. ANI
    take("taxonomy/type_strain_ani.csv", "5.ANI")
    take("taxonomy/ani_matrix_skani.tsv", "5.ANI")
    take("skani.tsv", "5.ANI", "skani_closest_genome.tsv")
    # 6. annotation
    if (sd / "bakta").is_dir():
        _copy_tree(sd / "bakta", F["6.Annotation"])
    # 7. AMR
    take("amrfinder.tsv", "7.AMR", "amrfinderplus.tsv")
    take("abricate_vfdb.tsv", "7.AMR", "virulence_vfdb.tsv")
    take("abricate_plasmidfinder.tsv", "7.AMR", "plasmidfinder.tsv")
    if (sd / "mob_recon").is_dir():
        _copy_tree(sd / "mob_recon", F["7.AMR"] / "MOB-suite")
    # 8. mapping
    if (sd / "mapping").is_dir():
        _copy_tree(sd / "mapping", F["8.readsMappingToRef"])

    files, missing = [], {}
    for f in FOLDERS:
        got = sorted(str(p.relative_to(root)) for p in F[f].rglob("*") if p.is_file())
        if not got:
            missing[f] = "this step produced no output in this run"
            (F[f] / "NOT_AVAILABLE.txt").write_text(
                f"{f}: not available — this step did not produce output in this run "
                "(the report says why).\n")
        files += got
    return {"dir": str(root), "sample": S, "folders": FOLDERS, "files": files, "missing": missing}
