#!/usr/bin/env python3
"""build_16s_typestrain_db.py — a full-length 16S reference with EVERY rRNA copy
of every NCBI type-strain genome, for ONT 16S (Emu) and PacBio (DADA2).

Why: single-copy 16S references (one sequence per species, e.g. NCBI RefSeq 16S)
miss the copy-to-copy differences inside one genome (often 0.1–1 %), so a read
from a minority copy matches a neighbouring species as well or better. Taking
all copies from the type-strain genomes puts that variation in the reference,
under the name the type strain defines.

It reads the genome cache that setup_typestrains.sh fills (~/.ngamp_wgs_refs)
and writes a NEW database next to the existing ones — nothing is replaced:

    backend/databases/16S_typestrain_multicopy/
        emu/species_taxid.fasta, emu/taxonomy.tsv          (Emu, ONT 16S)
        dada2/ts16S_multicopy_v<date>_toGenus_trainset.fa.gz   (PacBio assignTaxonomy)
        dada2/ts16S_multicopy_v<date>_species.fa.gz            (addSpecies, exact match)
        manifest.json                                      (counts, sources, date)

    python3 build_16s_typestrain_db.py                     # everything in the cache
    python3 build_16s_typestrain_db.py --genera Bacillus,Escherichia --out /tmp/x

How copies are found: minimap2 aligns one full-length RefSeq 16S per genus
(bacteria and archaea, ~3,000 sequences) to the genome — about a second per
genome — and each locus is cut at the ends of its best full-length alignment.
barrnap (Infernal, ~20 s per genome) is only the fallback for a genome where
nothing aligns. Kept: complete genes of 1,400–1,700 bp without ambiguous bases.
Taxonomy (phylum … species) comes from NCBI Taxonomy for each genome's organism,
looked up once and cached. Re-running resumes: finished genomes are skipped.

EVALUATION (v2.19.0, tests/mock/eval_16s_multicopy.py) — read before using it.
Reads from held-out, non-type-strain genomes of 17 species, all references
restricted to the same 15 genera:
  ONT / Emu    Emu default 15/15 species, Bray-Curtis 0.036 (3.3 % false)
               this DB    14/15, BC 0.214 (21 % false: K. pneumoniae's type
                          genome is a draft without a complete 16S, its reads
                          went to K. quasivariicola)
               Emu default + this DB merged: 15/15, BC 0.038 — no gain
  PacBio       SILVA 138.2 species 34 % correct / 1 % wrong / 65 % unassigned
               this DB     14 % / 0 % / 86 % (exact matching needs the strain's
                           own copy; one type strain per species rarely has it)
Most type-strain genomes are short-read drafts (median 1 complete 16S copy
per genome), so the "multi-copy" idea mostly fails on the data that exists.
Keep the existing databases as the default; this one is an experiment.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import wgs_typestrain_db as TSDB        # noqa: E402  (NCBI helpers, same User-Agent)

MIN_LEN, MAX_LEN = 1400, 1700
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def say(*a):
    print(*a, flush=True)


# ── tools ────────────────────────────────────────────────────────────────────
def tool(name: str, override: str = "") -> tuple[str, dict] | None:
    """(path, env) — conda envs of the WGS setup first, like ont_wgs_pipeline."""
    cands = [Path(override) / name] if override else []
    for base in (os.environ.get("CONDA_PREFIX_1"), os.environ.get("CONDA_EXE") and
                 str(Path(os.environ["CONDA_EXE"]).parent.parent),
                 str(Path.home() / "miniconda3"), str(Path.home() / "miniforge3"),
                 str(Path.home() / "anaconda3"), str(Path.home() / "mambaforge")):
        if base:
            for env in ("ngamp-wgs", "ngamp-phylo", "barrnap"):
                cands.append(Path(base) / "envs" / env / "bin" / name)
    w = shutil.which(name)
    if w:
        cands.append(Path(w))
    for p in cands:
        if p.exists() and os.access(p, os.X_OK):
            env = os.environ.copy()
            env["PATH"] = f"{p.parent}{os.pathsep}{env.get('PATH', '')}"
            return str(p), env
    return None


# ── sequences ────────────────────────────────────────────────────────────────
def read_fasta(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def merge(iv: list[tuple[int, int]], gap: int = 200) -> list[tuple[int, int]]:
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + gap:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def barrnap_16s(fa: Path, kingdom: str, br, threads: int = 1) -> list[tuple[str, int, int, str, bool]]:
    """(contig, start1, end1, strand, partial) of each 16S gene barrnap reports.
    Only the fallback: barrnap 1.x runs Infernal and takes ~20 s per genome."""
    p = subprocess.run([br[0], "--kingdom", kingdom, "--threads", str(threads), "--quiet", str(fa)],
                       capture_output=True, text=True, env=br[1])
    hits = []
    for line in p.stdout.splitlines():
        if line.startswith("#"):
            continue
        v = line.split("\t")
        if len(v) >= 9 and ("16S" in v[8]):
            hits.append((v[0], int(v[3]), int(v[4]), v[6], "partial" in v[8].lower()))
    return hits


def extract_copies(acc: str, genome: Path, panel: Path | None, kingdom: str, br, mm,
                   work: Path) -> tuple[list[str], str]:
    """All complete 16S copies of one genome, 5'→3', and how they were found.

    Each locus is where full-length reference 16S genes (≥ 1,450 bp, one per
    genus from NCBI RefSeq 16S) align; the gene ends come from the best of those
    alignments, extended over any unaligned query ends — the same gene model the
    reference uses, to within a few bases. barrnap is the fallback for a genome
    where nothing aligns."""
    contigs = dict((n.split()[0], s.upper()) for n, s in read_fasta(genome))
    seqs, any_hit = [], False
    if panel is not None and mm is not None:
        p = subprocess.run([mm[0], "-x", "asm20", "-N", "200", "-p", "0.3", "--secondary=yes",
                            "-t", "1", str(genome), str(panel)],
                           capture_output=True, text=True, env=mm[1])
        hits = defaultdict(list)                    # (contig, strand) -> hits
        for line in p.stdout.splitlines():
            v = line.split("\t")
            if len(v) < 12:
                continue
            any_hit = any_hit or int(v[10]) >= 300
            ql, qs, qe, st, tn, ts, te, nm = (int(v[1]), int(v[2]), int(v[3]), v[4], v[5],
                                              int(v[7]), int(v[8]), int(v[9]))
            if ql < 1450 or (qe - qs) < 0.9 * ql:
                continue
            if st == "+":
                gs, ge = ts - qs, te + (ql - qe)
            else:
                gs, ge = ts - (ql - qe), te + qs
            hits[(tn, st)].append((nm / max(1, int(v[10])), nm, gs, ge))
        for (tn, st), hs in hits.items():
            hs.sort(key=lambda h: h[2])
            groups, cur = [], []
            for h in hs:                            # overlapping placements = one locus
                if cur and h[2] > max(x[3] for x in cur):
                    groups.append(cur)
                    cur = []
                cur.append(h)
            if cur:
                groups.append(cur)
            seq_c = contigs.get(tn, "")
            for g in groups:
                best = max(g, key=lambda h: (h[0], h[1]))
                gs, ge = max(0, best[2]), min(len(seq_c), best[3])
                s = seq_c[gs:ge]
                seqs.append(revcomp(s) if st == "-" else s)
    how = "aligned"
    # Nothing aligned at all (a lineage the panel does not reach): barrnap on the
    # genome. When pieces aligned but none was complete, the genes are cut by
    # contig ends — a draft assembly — and barrnap would only find them partial.
    if not seqs and not any_hit and br is not None:
        how = "barrnap"
        tmp = work / f"{acc}.fna"
        with open(tmp, "w") as fh:
            for c, s in contigs.items():
                fh.write(f">{c}\n{s}\n")
        for ctg, a, b, strand, partial in barrnap_16s(tmp, kingdom, br):
            if not partial:
                s = contigs[ctg][a - 1:b]
                seqs.append(revcomp(s) if strand == "-" else s)
        tmp.unlink(missing_ok=True)
    return [s for s in seqs if MIN_LEN <= len(s) <= MAX_LEN and not re.search(r"[^ACGT]", s)], how


# ── taxonomy ─────────────────────────────────────────────────────────────────
def lineages(accs: list[str], cache_file: Path) -> dict[str, dict]:
    """accession -> {domain..species: (name, taxid)} from NCBI, cached."""
    known = {}
    if cache_file.exists():
        try:
            known = json.loads(cache_file.read_text())
        except Exception:
            known = {}
    need = [a for a in accs if a not in known]
    if need:
        say(f"   NCBI taxonomy for {len(need):,} genomes")
        org = {}
        for i in range(0, len(need), 100):
            chunk = need[i:i + 100]
            d = TSDB.get_json("/genome/accession/" + urllib.parse.quote(",".join(chunk))
                              + "/dataset_report?page_size=100")
            for r in (d or {}).get("reports", []):
                tid = (r.get("organism") or {}).get("tax_id")
                for a in {r.get("accession"), r.get("current_accession")}:
                    if a and tid:
                        org[a] = int(tid)
        tids = sorted(set(org.values()))
        cls = {}
        for i in range(0, len(tids), 200):
            d = TSDB.get_json("/taxonomy/taxon/" + ",".join(map(str, tids[i:i + 200]))
                              + "/dataset_report")
            for r in (d or {}).get("reports", []):
                t = r.get("taxonomy", {})
                c = t.get("classification") or {}
                cls[int(t.get("tax_id", 0))] = {k: [c[k].get("name", ""), c[k].get("id")]
                                                for k in RANKS if c.get(k)}
            if (i // 200) % 10 == 0:
                say(f"     {min(i + 200, len(tids)):,} / {len(tids):,} taxa")
        for a in need:
            if a in org and org[a] in cls:
                known[a] = cls[org[a]]
        cache_file.write_text(json.dumps(known))
    return known


# ── build ────────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", default=str(Path.home() / ".ngamp_wgs_refs"))
    ap.add_argument("--out", default=str(HERE / "databases" / "16S_typestrain_multicopy"))
    ap.add_argument("--genera", default="", help="comma-separated genera (default: all in the cache)")
    ap.add_argument("--n16s", default="", help="NCBI RefSeq 16S FASTA used to locate copies "
                    "(default: db_paths.json ncbi_16s_type)")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--tool-dir", default="", help="folder holding barrnap and minimap2")
    ap.add_argument("--no-register", action="store_true", help="do not add to db_paths.json")
    a = ap.parse_args(argv)

    cache, out = Path(a.cache).expanduser(), Path(a.out).expanduser()
    work = out / "work"
    (work / "copies").mkdir(parents=True, exist_ok=True)
    br, mm = tool("barrnap", a.tool_dir), tool("minimap2", a.tool_dir)
    if not mm:
        say("minimap2 not found — install the ONT-WGS tools first: bash ~/r16s-app/setup_wgs.sh")
        return 2
    if not br:
        say("  ! barrnap not found — genomes where no reference 16S aligns are left out")
    dbp_file = HERE / "databases" / "db_paths.json"
    dbp = json.loads(dbp_file.read_text()) if dbp_file.exists() else {}
    n16 = Path(a.n16s or dbp.get("ncbi_16s_type", "") or
               HERE / "databases" / "wgs" / "ncbi_16s" / "bacteria.16SrRNA.fna")
    if not n16.exists():
        say(f"  ! no RefSeq 16S file at {n16} — every genome gets a whole-genome barrnap (slow)")
        n16 = None

    # 1. genomes in the cache
    wanted = {g.strip() for g in a.genera.split(",") if g.strip()}
    recs = []
    for f in sorted((cache / "type_strains").glob("*.v3.json")):
        for r in json.loads(f.read_text()):
            if wanted and r.get("genus") not in wanted:
                continue
            g = cache / f"{r['accession']}.fna.gz"
            if g.exists() and g.stat().st_size > 1000:
                recs.append(r)
    say(f"1. {len(recs):,} type-strain genomes in {cache}")
    if not recs:
        say("   nothing to do — run setup_typestrains.sh first")
        return 1

    # 2. taxonomy
    say("2. Taxonomy")
    lin = lineages([r["accession"] for r in recs], out / "lineages.json")

    # 3. one full-length 16S per genus (RefSeq 16S, bacteria + archaea) to locate genes
    panel = work / "panel_16S.fa"
    if not panel.exists():
        srcs = [n16] if n16 else []
        arc = work / "archaea.16SrRNA.fna"
        if not arc.exists():
            try:
                import urllib.request
                url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
                with urllib.request.urlopen(urllib.request.Request(url, headers=TSDB.UA),
                                            timeout=300) as r, open(str(arc) + ".gz", "wb") as fh:
                    shutil.copyfileobj(r, fh)
                with gzip.open(str(arc) + ".gz", "rb") as src, open(arc, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                Path(str(arc) + ".gz").unlink()
            except Exception as e:
                say(f"  ! archaeal RefSeq 16S not downloaded ({e}) — archaea fall back to barrnap")
        if arc.exists():
            srcs.append(arc)
        best = {}
        for src in srcs:
            for name, sq in read_fasta(src):
                w = name.split()
                if len(w) >= 3 and len(sq) >= 1450 and not re.search(r"[^ACGT]", sq.upper()):
                    if w[1] not in best or len(sq) > len(best[w[1]][1]):
                        best[w[1]] = (w[0], sq.upper())
        if best:
            panel.write_text("".join(f">{acc}\n{sq}\n" for acc, sq in best.values()))
            say(f"   {len(best):,} genera in the locating panel")
    if not (panel.exists() and panel.stat().st_size):
        panel = None

    # 4. copies, one genome at a time (resumable)
    todo = [r for r in recs if not (work / "copies" / f"{r['accession']}.fa").exists()]
    say(f"3. Extracting 16S copies: {len(todo):,} genomes to do, {a.threads} at a time")
    t0, done, empty, fallback = time.time(), 0, 0, 0
    tmpd = Path(tempfile.mkdtemp(prefix="ts16s_"))

    def one(r):
        acc = r["accession"]
        L = lin.get(acc, {})
        kingdom = "arc" if (L.get("domain") or [""])[0] == "Archaea" else "bac"
        seqs, how = extract_copies(acc, cache / f"{acc}.fna.gz", panel, kingdom, br, mm, tmpd)
        tmp = work / "copies" / f"{acc}.fa.part"
        tmp.write_text("".join(f">{acc}_{i + 1}\n{s}\n" for i, s in enumerate(seqs)))
        tmp.replace(work / "copies" / f"{acc}.fa")
        return len(seqs), how

    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        futs = [ex.submit(one, r) for r in todo]
        for fu in as_completed(futs):
            n, how = fu.result()
            done += 1
            empty += n == 0
            fallback += how == "barrnap"
            if done % 200 == 0 or done == len(todo):
                rate = done / max(1.0, time.time() - t0)
                say(f"   {done:,} / {len(todo):,} genomes · ~{(len(todo) - done) / rate / 60:.0f} min left")
    shutil.rmtree(tmpd, ignore_errors=True)
    if fallback:
        say(f"   {fallback} genome(s) needed barrnap (nothing aligned to the panel)")

    # 5. write the databases
    say("4. Writing Emu and DADA2 files")
    stamp = time.strftime("%Y%m%d")
    (out / "emu").mkdir(exist_ok=True)
    (out / "dada2").mkdir(exist_ok=True)
    for old in (out / "dada2").glob("ts16S_multicopy_v*"):
        old.unlink()
    tax_rows, n_seq, n_genomes, no_tax, no_16s = {}, 0, 0, 0, 0
    copies_per_genome = []
    seen = set()                                    # (species taxid, sequence)
    tr_name = out / "dada2" / f"ts16S_multicopy_v{stamp}_toGenus_trainset.fa.gz"
    sp_name = out / "dada2" / f"ts16S_multicopy_v{stamp}_species.fa.gz"
    with open(out / "emu" / "species_taxid.fasta", "w") as emu_fa, \
            gzip.open(tr_name, "wt") as tr, gzip.open(sp_name, "wt") as sp:
        for r in recs:
            acc = r["accession"]
            L = lin.get(acc)
            seqs = [s for _, s in read_fasta(work / "copies" / f"{acc}.fa")] \
                if (work / "copies" / f"{acc}.fa").exists() else []
            if not seqs:
                no_16s += 1
                continue
            if not L or not L.get("species") or not L.get("genus"):
                no_tax += 1
                continue
            n_genomes += 1
            copies_per_genome.append(len(seqs))
            sp_name_, sp_id = L["species"]
            nm = {k: (L.get(k) or [""])[0] for k in RANKS}
            tax_rows[sp_id] = [sp_id, sp_name_, nm["genus"], nm["family"], nm["order"], nm["class"],
                               nm["phylum"], "", nm["domain"]]
            lineage = ";".join([nm["domain"], nm["phylum"], nm["class"], nm["order"], nm["family"],
                                nm["genus"]]) + ";"
            epithet = sp_name_.split(" ", 1)[1] if " " in sp_name_ else sp_name_
            for i, s in enumerate(seqs):
                if (sp_id, s) in seen:                  # identical copies: once per species
                    continue
                seen.add((sp_id, s))
                n_seq += 1
                emu_fa.write(f">{sp_id}:ts16S:{n_seq} [{acc} copy {i + 1}]\n{s}\n")
                tr.write(f">{lineage}\n{s}\n")
                sp.write(f">{acc}_{i + 1} {nm['genus']} {epithet}\n{s}\n")
    with open(out / "emu" / "taxonomy.tsv", "w") as fh:
        fh.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\n")
        for row in sorted(tax_rows.values(), key=lambda x: x[1]):
            fh.write("\t".join(str(x) for x in row) + "\n")
    cpg = sorted(copies_per_genome)
    man = {"built": time.strftime("%Y-%m-%d %H:%M"), "source": "NCBI 'assembly from type material' "
           "genomes in " + str(cache), "genomes_used": n_genomes, "species": len(tax_rows),
           "sequences": n_seq, "genomes_without_complete_16S": no_16s,
           "genomes_without_taxonomy": no_tax,
           "copies_per_genome_median": cpg[len(cpg) // 2] if cpg else 0,
           "copies_per_genome_max": cpg[-1] if cpg else 0,
           "length_range": [MIN_LEN, MAX_LEN], "genera_filter": sorted(wanted) or "all",
           "emu_db": str(out / "emu"), "dada2_trainset": str(tr_name), "dada2_species": str(sp_name)}
    (out / "manifest.json").write_text(json.dumps(man, indent=1))
    say(f"   {n_seq:,} distinct 16S sequences, {len(tax_rows):,} species, {n_genomes:,} genomes "
        f"(median {man['copies_per_genome_median']} copies per genome)")
    if no_16s or no_tax:
        say(f"   left out: {no_16s} genome(s) without a complete 16S, {no_tax} without NCBI taxonomy")

    if not a.no_register and not wanted:
        dbp["emu_typestrain_multicopy"] = str(out / "emu")
        dbp["TS16S_multicopy"] = str(tr_name)
        dbp_file.parent.mkdir(parents=True, exist_ok=True)
        dbp_file.write_text(json.dumps(dbp, indent=2))
        say("   registered in db_paths.json: emu_typestrain_multicopy, TS16S_multicopy")
    say("Done. ONT 16S: choose the Emu database '16S_typestrain_multicopy/emu'. "
        "PacBio: choose ts16S_multicopy_…_toGenus_trainset.fa.gz.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
