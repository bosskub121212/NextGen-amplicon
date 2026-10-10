"""Type-strain taxonomy and reference comparison for the ONT-WGS pipeline.

What the DSMZ TYGS server and JSpeciesWS do for one isolate, rebuilt from
open tools so it runs offline-first on the analysis computer:

  1. 16S vs type strains   vsearch against NCBI RefSeq 16S (type material only)
  2. type-strain genomes   NCBI Datasets "assembly from type material" for the
                           candidate genera; skani picks the closest ones
  3. JSpecies-style ANI    ANIb (Goris et al. 2007: 1020-nt fragments, blastn,
                           ≥70 % coverage / ≥30 % identity, mean pident — the
                           pyANI-plus implementation), ANIm (nucmer --mum,
                           delta-filter -1) and TETRA (Teeling et al. 2004
                           tetranucleotide z-score correlation)
  4. species verdict       known species / borderline / potential novel species
  5. trees                 genome tree: FastME (BioNJ + SPR) on skani ANI
                           distances; 16S tree: MAFFT + IQ-TREE (ModelFinder,
                           1000 ultrafast bootstraps)
  6. reference mapping     minimap2 reads → closest genome (coverage, absent
                           regions) and assembly → genome (SNPs, indels,
                           isolate-specific regions, synteny dot plot)

dDDH (GBDP / GGDC) is not reproduced: the DSMZ software is not distributed,
and an imitation would be a different number with the same name.

The pipeline binds its helpers (tool discovery, logging, downloads) with
bind() so this module shares its warning list and tool cache.
"""
from __future__ import annotations

import gzip
import json
import math
import os
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

H = SimpleNamespace()          # find_tool, run, log, warn, fetch_genome, read_fasta, write_fasta

NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_16S_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
UA = {"User-Agent": "NextGen-Amplicon (ONT-WGS type-strain module)"}

ANIB_FRAG = 1020
ANIB_MIN_COV = 0.7
ANIB_MIN_ID = 0.3
SPECIES_ANI = 95.0
SPECIES_ANI_STRICT = 96.0
TETRA_SPECIES = 0.99
RRNA_SPECIES = 98.65          # Kim et al. 2014
MAP_VERSION = 1


def bind(**kw):
    for k, v in kw.items():
        setattr(H, k, v)


# ══════════════════════════════════════════════════════════════════════════════
#  NCBI Datasets
# ══════════════════════════════════════════════════════════════════════════════
def ncbi_json(path: str, timeout: int = 60):
    url = NCBI_API + path
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=timeout) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < 3:
                time.sleep(1.5 * (attempt + 1))
                continue
            H.log(f"    NCBI request failed ({e.code}): {path[:120]}")
            return None
        except Exception as e:
            if attempt < 2:
                time.sleep(1.5 * (attempt + 1))
                continue
            H.log(f"    NCBI request failed ({e}): {path[:120]}")
            return None
    return None


def genus_taxid(name: str, cache: Path) -> tuple[int | None, str]:
    """NCBI taxid of a bacterial/archaeal genus by name (GTDB '_A' suffixes dropped)."""
    name = re.sub(r"_[A-Z]+$", "", (name or "").strip())
    if not name:
        return None, ""
    cf = cache / "genus_taxids.json"
    known = {}
    if cf.exists():
        try:
            known = json.loads(cf.read_text())
        except Exception:
            known = {}
    if name in known:
        return known[name], name
    d = ncbi_json(f"/taxonomy/taxon/{urllib.parse.quote(name)}/dataset_report")
    tid = None
    for r in (d or {}).get("reports", []):
        t = r.get("taxonomy", {})
        dom = (t.get("classification", {}).get("domain") or {}).get("name", "")
        if t.get("rank") == "GENUS" and dom in ("Bacteria", "Archaea"):
            tid = t.get("tax_id")
            break
    if tid:
        known[name] = tid
        cache.mkdir(parents=True, exist_ok=True)
        cf.write_text(json.dumps(known, indent=1))
    return tid, name


LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def species_key(name: str) -> str:
    """'Cronobacter sakazakii NBRC 102416' → 'Cronobacter sakazakii' (subsp. kept)."""
    w = (name or "").split()
    # subspecies are taxa of their own; serovars / pathovars / biovars are not,
    # and the same type strain is often deposited under both names
    if len(w) >= 4 and w[2] == "subsp.":
        return " ".join(w[:4])
    return " ".join(w[:2])


def clean_strain(strain: str) -> str:
    """'KCTC3922(T)' / 'BCT-7112T' / 'DSM 1T' → strain without its type mark."""
    st = re.sub(r"\s*\((T|t)\)\s*$", "", (strain or "").strip())
    st = re.sub(r"(?<=[0-9])T$", "", st)
    st = re.sub(r"^(type strain|strain)\s*:?\s*", "", st, flags=re.I)
    return st.strip()


def ts_label(organism: str, strain: str, acc: str = "") -> str:
    st = clean_strain(strain)
    return " ".join(x for x in (organism, (st + "T") if st else "T", acc) if x) \
        if st else " ".join(x for x in (organism + " T", acc) if x)


def type_strains(genus: str, cache: Path, max_age_days: int = 90) -> list[dict]:
    """One genome per (sub)species from NCBI 'assembly from type material'.

    Cached per genus; a stale cache is still used when NCBI cannot be reached.
    """
    tid, gname = genus_taxid(genus, cache)
    if not tid:
        return []
    cf = cache / "type_strains" / f"{tid}.v3.json"
    if cf.exists() and time.time() - cf.stat().st_mtime < max_age_days * 86400:
        try:
            return json.loads(cf.read_text())
        except Exception:
            pass
    reps, token = [], ""
    for _ in range(20):
        q = (f"/genome/taxon/{tid}/dataset_report?filters.is_type_material=true"
             f"&filters.exclude_atypical=true&page_size=1000"
             + (f"&page_token={urllib.parse.quote(token)}" if token else ""))
        d = ncbi_json(q, timeout=120)
        if d is None:
            reps = None
            break
        reps += d.get("reports", [])
        token = d.get("next_page_token") or ""
        if not token:
            break
    if reps is None:
        if cf.exists():
            H.log(f"    NCBI unreachable — using cached type-strain list for {gname}")
            return json.loads(cf.read_text())
        return []
    best: dict[str, dict] = {}
    for r in reps:
        acc = r.get("accession", "")
        org = r.get("organism", {})
        name = org.get("organism_name", "")
        if not acc or not name or " sp." in name:
            continue
        st = r.get("assembly_stats", {})
        info = r.get("assembly_info", {})
        key_name = species_key(name)
        strain = clean_strain((org.get("infraspecific_names") or {}).get("strain", "") or
                              name[len(key_name):].strip())
        rec = {"accession": acc, "paired": r.get("paired_accession", ""),
               "organism": key_name, "organism_ncbi": name, "strain": strain,
               "level": info.get("assembly_level", ""),
               "length": int(st.get("total_sequence_length") or 0),
               "n50": int(st.get("contig_n50") or 0),
               "release": info.get("release_date", ""),
               "type_label": (r.get("type_material") or {}).get("type_display_text", ""),
               "genus": gname}
        key = (LEVEL_RANK.get(rec["level"], 9), 0 if acc.startswith("GCF_") else 1,
               -rec["n50"], rec["release"])
        cur = best.get(key_name)
        if cur is None or key < cur["_key"]:
            rec["_key"] = key
            best[key_name] = rec
    out = []
    for rec in best.values():
        rec.pop("_key", None)
        out.append(rec)
    out.sort(key=lambda r: r["organism"])
    cf.parent.mkdir(parents=True, exist_ok=True)
    cf.write_text(json.dumps(out, indent=1))
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  16S vs NCBI RefSeq type-strain 16S
# ══════════════════════════════════════════════════════════════════════════════
def ensure_16s_db(dbp: dict, db_root: Path) -> Path | None:
    p = dbp.get("ncbi_16s_type", "")
    if p and Path(p).exists():
        return Path(p)
    d = db_root / "ncbi_16s"
    fa = d / "bacteria.16SrRNA.fna"
    if fa.exists() and fa.stat().st_size > 1_000_000:
        return fa
    try:
        d.mkdir(parents=True, exist_ok=True)
        H.log("    downloading NCBI RefSeq 16S type-strain database (≈8 MB)")
        tmp = d / "bacteria.16SrRNA.fna.gz.part"
        with urllib.request.urlopen(urllib.request.Request(NCBI_16S_URL, headers=UA),
                                    timeout=300) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        with gzip.open(tmp, "rb") as src, open(fa, "wb") as dst:
            shutil.copyfileobj(src, dst)
        tmp.unlink()
        return fa
    except Exception as e:
        H.warn(f"NCBI 16S type-strain database not available ({e}) — run setup_wgs.sh "
               "--dbs-only; 16S comparison with type strains skipped")
        return None


def parse_16s_header(h: str) -> dict:
    """'NR_119358.1 Acinetobacter baumannii strain ATCC 19606 16S ribosomal RNA, …'"""
    acc, _, rest = h.partition(" ")
    rest = re.sub(r"\s+16S ribosomal RNA.*$", "", rest).strip()
    w = rest.split()
    if len(w) >= 4 and w[2] in ("subsp.", "pv.", "bv."):
        sp, strain = " ".join(w[:4]), " ".join(w[4:])
    else:
        sp, strain = " ".join(w[:2]), " ".join(w[2:])
    strain = clean_strain(re.sub(r"^strain\s+", "", strain))
    return {"accession": acc, "species": sp, "strain": strain}


def search_16s(seqs: list[tuple[str, str]], db: Path, work: Path, threads: int,
               logf: Path) -> list[dict]:
    """Best identity of any 16S copy to each type strain (one row per NR record)."""
    vs = H.find_tool("vsearch")
    if not vs or not seqs or not db:
        return []
    q = work / "ts16S_query.fasta"
    uniq = {}
    for n, s in seqs:
        uniq.setdefault(s.upper(), n.split()[0])
    H.write_fasta(q, [(n, s) for s, n in uniq.items()])
    out = work / "ts16S_hits.b6"
    if H.run(vs, ["--usearch_global", q, "--db", db, "--id", "0.85", "--maxaccepts", "300",
                  "--maxrejects", "3000", "--strand", "both", "--query_cov", "0.6",
                  "--threads", threads, "--blast6out", out, "--notrunclabels"],
             logf) != 0 or not out.exists():
        return []
    best: dict[str, dict] = {}
    for line in out.read_text().splitlines():
        v = line.split("\t")
        if len(v) < 4:
            continue
        ident, alen = float(v[2]), int(v[3])
        rec = parse_16s_header(v[1])
        cur = best.get(rec["accession"])
        if cur is None or ident > cur["identity"]:
            rec.update({"identity": round(ident, 2), "aln_len": alen,
                        "copy": v[0].split()[0]})
            best[rec["accession"]] = rec
    rows = sorted(best.values(), key=lambda r: (-r["identity"], -r["aln_len"]))
    return rows


def species_16s(rows: list[dict]) -> list[dict]:
    """Collapse NR records to species (best record per species)."""
    seen, out = set(), []
    for r in rows:
        if r["species"] in seen:
            continue
        seen.add(r["species"])
        out.append(r)
    return out


def fetch_16s_records(db: Path, accs: set[str]) -> dict[str, str]:
    out = {}
    if not accs:
        return out
    for name, seq in H.read_fasta(db):
        a = name.split()[0]
        if a in accs:
            out[a] = seq
            if len(out) == len(accs):
                break
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  ANIb / ANIm / TETRA
# ══════════════════════════════════════════════════════════════════════════════
def _plain_fasta(src: Path, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with (gzip.open(src, "rt") if str(src).endswith(".gz") else open(src)) as fi, \
            open(dest, "w") as fo:
        shutil.copyfileobj(fi, fo)
    return dest


def fragment(fa: Path, out: Path, size: int = ANIB_FRAG) -> int:
    n = 0
    with open(out, "w") as fh:
        for _, seq in H.read_fasta(fa):
            for i in range(0, len(seq), size):
                n += 1
                fh.write(f">frag{n:05d}\n{seq[i:i + size]}\n")
    return n


def parse_anib(tsv: Path) -> tuple[float | None, int, int]:
    """(mean pident of accepted fragments, aligned bases, accepted fragments)."""
    tot, cnt, alen_sum, prev = 0.0, 0, 0, ""
    if not tsv.exists():
        return None, 0, 0
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            q, pid, length, mism, qlen, gaps = f[0], float(f[2]), int(f[3]), int(f[4]), \
                int(f[5]), int(f[6])
            alnlen = length - gaps
            cov = alnlen / qlen
            idn = (alnlen - mism) / qlen
            if cov > ANIB_MIN_COV and idn > ANIB_MIN_ID and q != prev:
                tot += pid
                cnt += 1
                alen_sum += alnlen
                prev = q
    return (tot / cnt if cnt else None), alen_sum, cnt


def _merge_len(iv: list[tuple[int, int]]) -> int:
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    if cs is not None:
        tot += ce - cs + 1
    return tot


def parse_delta(path: Path) -> dict:
    """ANIm as in pyANI-plus: weighted identical bases / aligned bases."""
    reg_r, reg_q = defaultdict(list), defaultdict(list)
    aligned = weighted = 0
    cur_r = cur_q = None
    has = False
    for line in path.read_text().splitlines():
        v = line.split()
        if not v or v[0] == "NUCMER":
            continue
        if v[0].startswith(">"):
            cur_r, cur_q = v[0][1:], v[1]
            continue
        if len(v) == 7:
            has = True
            r = tuple(sorted((int(v[0]), int(v[1]))))
            q = tuple(sorted((int(v[2]), int(v[3]))))
            reg_r[cur_r].append(r)
            reg_q[cur_q].append(q)
            rl, ql = r[1] - r[0] + 1, q[1] - q[0] + 1
            aligned += rl + ql
            weighted += rl + ql - 2 * int(v[4])
    if not has:
        return {}
    return {"ani": 100.0 * weighted / aligned if aligned else None,
            "ref_aligned": sum(_merge_len(x) for x in reg_r.values()),
            "query_aligned": sum(_merge_len(x) for x in reg_q.values())}


_COMP = str.maketrans("ACGT", "TGCA")


def tetra_z(fa: Path) -> dict[str, float]:
    """Teeling et al. (2004) tetranucleotide z-scores on both strands (pyani)."""
    c4, c3, c2 = Counter(), Counter(), Counter()
    for _, seq in H.read_fasta(fa):
        s = seq.upper()
        for strand in (s, s.translate(_COMP)[::-1]):
            for part in re.split(r"[^ACGT]+", strand):
                n = len(part)
                if n < 4:
                    continue
                c4.update(part[i:i + 4] for i in range(n - 3))
                c3.update(part[i:i + 3] for i in range(n - 2))
                c2.update(part[i:i + 2] for i in range(n - 1))
    z = {}
    bases = "ACGT"
    for a in bases:
        for b in bases:
            for c in bases:
                for d in bases:
                    t = a + b + c + d
                    n23 = c2[b + c]
                    n123, n234 = c3[a + b + c], c3[b + c + d]
                    if not n23:
                        z[t] = 0.0
                        continue
                    exp = n123 * n234 / n23
                    var = exp * ((n23 - n123) * (n23 - n234)) / (n23 * n23)
                    z[t] = (c4[t] - exp) / math.sqrt(var) if var > 0 else 0.0
    return z


def pearson(a: dict, b: dict) -> float | None:
    keys = sorted(a)
    x = [a[k] for k in keys]
    y = [b.get(k, 0.0) for k in keys]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((i - mx) * (j - my) for i, j in zip(x, y))
    sxx = sum((i - mx) ** 2 for i in x)
    syy = sum((j - my) ** 2 for j in y)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else None


def jspecies(sample_fa: Path, refs: list[dict], work: Path, threads: int,
             logf: Path, cache_file: Path) -> dict[str, dict]:
    """ANIb (both directions), ANIm and TETRA of the sample against each reference.

    Results are cached per (assembly checksum, reference accession) so a
    Re-analyse only computes new pairs.
    """
    blastn, mkdb = H.find_tool("blastn"), H.find_tool("makeblastdb")
    nucmer, dfilt = H.find_tool("nucmer"), H.find_tool("delta-filter")
    work.mkdir(parents=True, exist_ok=True)
    import hashlib
    smd5 = hashlib.md5(sample_fa.read_bytes()).hexdigest()[:12]
    cache = {}
    if cache_file.exists():
        try:
            cache = json.loads(cache_file.read_text())
        except Exception:
            cache = {}
    out: dict[str, dict] = {}
    todo = []
    for r in refs:
        key = f"{smd5}|{r['accession']}"
        if key in cache and cache[key].get("v") == 2:
            out[r["accession"]] = cache[key]
        else:
            todo.append(r)
    if not todo:
        return out

    s_plain = _plain_fasta(sample_fa, work / "sample.fna")
    s_len = sum(len(x) for _, x in H.read_fasta(s_plain))
    s_frag = work / "sample.frag"
    fragment(s_plain, s_frag)
    s_tetra = tetra_z(s_plain)
    if blastn and mkdb:
        H.run(mkdb, ["-in", s_plain, "-dbtype", "nucl", "-out", work / "sample_db"],
              logf, quiet_ok=True)

    def one(r):
        acc = r["accession"]
        rd = work / acc
        rd.mkdir(exist_ok=True)
        rf = _plain_fasta(Path(r["path"]), rd / "ref.fna")
        r_len = sum(len(x) for _, x in H.read_fasta(rf))
        res = {"v": 2, "accession": acc, "ref_len": r_len, "query_len": s_len}
        if blastn and mkdb:
            fragment(rf, rd / "ref.frag")
            H.run(mkdb, ["-in", rf, "-dbtype", "nucl", "-out", rd / "db"], logf, quiet_ok=True)
            cols = "6 qseqid sseqid pident length mismatch qlen gaps"
            common = ["-task", "blastn", "-outfmt", cols, "-xdrop_gap_final", "150",
                      "-dust", "no", "-evalue", "1e-15", "-num_threads", "1"]
            H.run(blastn, ["-query", s_frag, "-db", rd / "db", "-out", rd / "q_vs_r.tsv"]
                  + common, logf, quiet_ok=True)
            H.run(blastn, ["-query", rd / "ref.frag", "-db", work / "sample_db",
                           "-out", rd / "r_vs_q.tsv"] + common, logf, quiet_ok=True)
            a1, l1, _ = parse_anib(rd / "q_vs_r.tsv")
            a2, l2, _ = parse_anib(rd / "r_vs_q.tsv")
            res.update({"anib_qr": a1, "anib_rq": a2,
                        "anib_cov_q": 100.0 * l1 / s_len if s_len else None,
                        "anib_cov_r": 100.0 * l2 / r_len if r_len else None})
            vals = [x for x in (a1, a2) if x is not None]
            res["anib"] = sum(vals) / len(vals) if vals else None
        if nucmer and dfilt:
            pre = rd / "nucmer"
            H.run(nucmer, ["-p", pre, "--mum", rf, s_plain], logf, quiet_ok=True)
            delta = Path(str(pre) + ".delta")
            if delta.exists():
                H.run(dfilt, ["-1", delta], logf, stdout_to=rd / "nucmer.filter",
                      quiet_ok=True)
                dm = parse_delta(rd / "nucmer.filter") if (rd / "nucmer.filter").exists() \
                    else {}
                if dm:
                    res.update({"anim": dm["ani"],
                                "anim_cov_q": 100.0 * dm["query_aligned"] / s_len,
                                "anim_cov_r": 100.0 * dm["ref_aligned"] / r_len})
        res["tetra"] = pearson(s_tetra, tetra_z(rf))
        shutil.rmtree(rd, ignore_errors=True)
        return res

    with ThreadPoolExecutor(max_workers=max(1, min(len(todo), int(threads)))) as ex:
        for res in ex.map(one, todo):
            out[res["accession"]] = res
            cache[f"{smd5}|{res['accession']}"] = res
    cache_file.write_text(json.dumps(cache, indent=1, default=str))
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  Trees: Newick, NJ fallback, midpoint root, layout
# ══════════════════════════════════════════════════════════════════════════════
class Node:
    __slots__ = ("name", "length", "children", "parent", "support")

    def __init__(self, name="", length=0.0):
        self.name, self.length, self.children, self.parent, self.support = \
            name, length, [], None, None


def parse_newick(s: str) -> Node:
    s = s.strip().rstrip(";")
    if not s:
        raise ValueError("empty Newick string")
    pos = 0

    def label():
        nonlocal pos
        st = pos
        if pos < len(s) and s[pos] == "'":
            pos += 1
            st = pos
            while s[pos] != "'":
                pos += 1
            name = s[st:pos]
            pos += 1
        else:
            while pos < len(s) and s[pos] not in ":,()":
                pos += 1
            name = s[st:pos]
        ln = 0.0
        if pos < len(s) and s[pos] == ":":
            pos += 1
            st = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            ln = max(0.0, float(s[st:pos] or 0))
        return name, ln

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
            nm, ln = label()
            if nm:
                try:
                    n.support = float(nm.split("/")[-1])
                except ValueError:
                    n.name = nm
            n.length = ln
        else:
            n.name, n.length = label()
        return n

    return node()


def to_newick(n: Node) -> str:
    def rec(x):
        if x.children:
            sup = f"{x.support:g}" if x.support is not None else ""
            return "(" + ",".join(rec(c) for c in x.children) + f"){sup}:{x.length:.6f}"
        nm = x.name if re.match(r"^[\w.\-|]+$", x.name) else "'" + x.name.replace("'", "") + "'"
        return f"{nm}:{x.length:.6f}"
    return rec(n).rsplit(":", 1)[0] + ";"


def leaves(n: Node):
    if not n.children:
        yield n
    for c in n.children:
        yield from leaves(c)


def midpoint_root(root: Node) -> Node:
    lv = list(leaves(root))
    if len(lv) < 3:
        return root
    # undirected adjacency
    adj = defaultdict(list)

    def build(n):
        for c in n.children:
            adj[id(n)].append((c, c.length))
            adj[id(c)].append((n, c.length))
            build(c)
    build(root)
    nodes = {}

    def collect(n):
        nodes[id(n)] = n
        for c in n.children:
            collect(c)
    collect(root)

    def far(start):
        best, dist, prev = (start, 0.0), {id(start): 0.0}, {id(start): None}
        stack = [start]
        while stack:
            x = stack.pop()
            for y, w in adj[id(x)]:
                if id(y) not in dist:
                    dist[id(y)] = dist[id(x)] + w
                    prev[id(y)] = x
                    stack.append(y)
                    if not y.children and dist[id(y)] > best[1]:
                        best = (y, dist[id(y)])
        return best, prev, dist
    (a, _), _, _ = far(lv[0])
    (b, dab), prev, dist = far(a)
    if dab <= 0:
        return root
    half = dab / 2
    # walk from b back towards a until passing the midpoint
    x = b
    while prev[id(x)] is not None and dist[id(prev[id(x)])] > half:
        x = prev[id(x)]
    y = prev[id(x)]                     # edge (y, x), dist[y] <= half < dist[x]
    if y is None:
        return root
    w = dist[id(x)] - dist[id(y)]
    new = Node()
    left = half - dist[id(y)]           # from y to the new root

    def reroot(cur, frm):
        nn = Node(cur.name)
        nn.support = cur.support
        for z, wz in adj[id(cur)]:
            if frm is not None and z is frm:
                continue
            ch = reroot(z, cur)
            ch.length = wz
            ch.parent = nn
            nn.children.append(ch)
        return nn
    cx = reroot(x, y)
    cx.length = w - left
    cy = reroot(y, x)
    cy.length = left
    for c in (cx, cy):
        c.parent = new
        new.children.append(c)
    # collapse unary nodes left by the old root
    def squash(n):
        for i, c in enumerate(list(n.children)):
            squash(c)
            if len(c.children) == 1:
                g = c.children[0]
                g.length += c.length
                g.parent = n
                n.children[i] = g
    squash(new)
    return new


def nj(names: list[str], D: list[list[float]]) -> Node:
    """Plain neighbour joining (fallback when FastME is not installed)."""
    nodes = [Node(n) for n in names]
    D = [row[:] for row in D]
    while len(nodes) > 2:
        n = len(nodes)
        r = [sum(row) for row in D]
        best, bi, bj = None, 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                q = (n - 2) * D[i][j] - r[i] - r[j]
                if best is None or q < best:
                    best, bi, bj = q, i, j
        dij = D[bi][bj]
        li = max(0.0, 0.5 * dij + (r[bi] - r[bj]) / (2 * (n - 2)))
        lj = max(0.0, dij - li)
        u = Node()
        for k, ln in ((bi, li), (bj, lj)):
            nodes[k].length = ln
            nodes[k].parent = u
            u.children.append(nodes[k])
        newrow = [0.5 * (D[bi][k] + D[bj][k] - dij) for k in range(n)]
        keep = [k for k in range(n) if k not in (bi, bj)]
        D = [[D[a][b] for b in keep] + [newrow[a]] for a in keep]
        D.append([newrow[k] for k in keep] + [0.0])
        nodes = [nodes[k] for k in keep] + [u]
    root = Node()
    a, b = nodes
    a.length = b.length = D[0][1] / 2
    for c in (a, b):
        c.parent = root
        root.children.append(c)
    return root


def tree_layout(root: Node) -> dict:
    """Leaves in order with x (cumulative length) / y (row); internal nodes too."""
    out = {"leaves": [], "edges": [], "max_x": 0.0}
    row = [0]

    def rec(n, x):
        if not n.children:
            y = row[0]
            row[0] += 1
            out["leaves"].append({"name": n.name, "x": x, "y": y})
            out["max_x"] = max(out["max_x"], x)
            return y
        ys = []
        for c in n.children:
            yc = rec(c, x + c.length)
            out["edges"].append({"x0": x, "x1": x + c.length, "y": yc})
            ys.append(yc)
        y = (min(ys) + max(ys)) / 2
        out["edges"].append({"x0": x, "x1": x, "y0": min(ys), "y1": max(ys),
                             "support": n.support, "yn": y})
        return y
    rec(root, 0.0)
    out["rows"] = row[0]
    return out


def genome_tree(genomes: list[tuple[str, Path]], work: Path, threads: int,
                logf: Path) -> tuple[str, list[list[float | None]]]:
    """FastME (BioNJ + SPR) on 1 − ANI/100 from skani triangle; NJ fallback."""
    sk = H.find_tool("skani")
    if not sk or len(genomes) < 3:
        return "", []
    work.mkdir(parents=True, exist_ok=True)
    lst = work / "genomes.txt"
    lst.write_text("\n".join(str(p) for _, p in genomes) + "\n")
    outf = work / "skani_triangle.tsv"
    if H.run(sk, ["triangle", "-l", lst, "-E", "-t", threads, "-o", outf], logf) != 0:
        return "", []
    idx = {str(p): i for i, (_, p) in enumerate(genomes)}
    n = len(genomes)
    ani = [[None] * n for _ in range(n)]
    for i in range(n):
        ani[i][i] = 100.0
    for line in outf.read_text().splitlines()[1:]:
        v = line.split("\t")
        if len(v) < 3:
            continue
        a, b = idx.get(v[0]), idx.get(v[1])
        if a is None or b is None:
            continue
        x = float(v[2])
        ani[a][b] = ani[b][a] = x
    D = [[(1 - (ani[i][j] if ani[i][j] is not None else 70.0) / 100.0) for j in range(n)]
         for i in range(n)]
    ids = [f"g{i}" for i in range(n)]
    fm = H.find_tool("fastme")
    tree = None
    if fm:
        phy = work / "dist.phy"
        with open(phy, "w") as fh:
            fh.write(f"{n}\n")
            for i in range(n):
                fh.write(ids[i].ljust(10) + " " + " ".join(f"{d:.6f}" for d in D[i]) + "\n")
        nwk = work / "fastme.nwk"
        # FastME refuses fewer than 4 taxa and then writes an empty file
        if n >= 4 and H.run(fm, ["-i", phy, "-o", nwk, "-m", "B", "-s", "-T", 1], logf,
                            quiet_ok=True) == 0 and nwk.exists() and nwk.read_text().strip():
            try:
                tree = parse_newick(nwk.read_text())
            except Exception:
                tree = None
    if tree is None:
        tree = nj(ids, D)
    tree = midpoint_root(tree)
    names = {ids[i]: genomes[i][0] for i in range(n)}
    for lf in leaves(tree):
        lf.name = names.get(lf.name, lf.name)
    return to_newick(tree), ani


def rrna_tree(seqs: list[tuple[str, str]], work: Path, threads: int,
              logf: Path) -> tuple[str, str]:
    """MAFFT + IQ-TREE (ModelFinder, 1000 UFBoot). Returns (newick, method)."""
    if len(seqs) < 4:
        return "", ""
    work.mkdir(parents=True, exist_ok=True)
    ids = [f"s{i}" for i in range(len(seqs))]
    fa = work / "16S_input.fasta"
    H.write_fasta(fa, [(ids[i], s) for i, (_, s) in enumerate(seqs)])
    mafft = H.find_tool("mafft")
    if not mafft:
        return "", ""
    aln = work / "16S_aligned.fasta"
    if H.run(mafft, ["--auto", "--adjustdirection", "--thread", threads, fa], logf,
             stdout_to=aln) != 0:
        return "", ""
    # MAFFT prefixes reversed sequences with _R_
    recs = [(re.sub(r"^_R_", "", n), s.upper()) for n, s in H.read_fasta(aln)]
    # trim ragged ends: keep columns where ≥ 50 % of sequences have a base
    L = len(recs[0][1]) if recs else 0
    keep = [j for j in range(L)
            if sum(1 for _, s in recs if s[j] not in "-N") >= 0.5 * len(recs)]
    trimmed = [(n, "".join(s[j] for j in keep)) for n, s in recs]
    tfa = work / "16S_trimmed.fasta"
    H.write_fasta(tfa, trimmed)
    method = ""
    nwk_txt = ""
    iq = H.find_tool("iqtree3") or H.find_tool("iqtree2") or H.find_tool("iqtree")
    if iq:
        pre = work / "iq"
        if H.run(iq, ["-s", tfa, "-m", "MFP", "-B", "1000", "-T", min(4, int(threads)),
                      "--prefix", pre, "-redo", "-quiet", "--seed", "12345"], logf,
                 timeout=3600, quiet_ok=True) == 0 and Path(str(pre) + ".treefile").exists():
            nwk_txt = Path(str(pre) + ".treefile").read_text()
            model = ""
            iqf = Path(str(pre) + ".iqtree")
            if iqf.exists():
                m = re.search(r"Best-fit model according to BIC:\s*(\S+)", iqf.read_text())
                model = m.group(1) if m else ""
            method = f"IQ-TREE (ML, {model or 'ModelFinder'}; 1000 ultrafast bootstraps)"
    if not nwk_txt:
        # p-distance NJ fallback
        n = len(trimmed)
        D = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                a, b = trimmed[i][1], trimmed[j][1]
                comp = [(x, y) for x, y in zip(a, b) if x in "ACGT" and y in "ACGT"]
                d = sum(1 for x, y in comp if x != y) / len(comp) if comp else 0.5
                D[i][j] = D[j][i] = d
        tree = nj([n_ for n_, _ in trimmed], D)
        method = "neighbour joining (p-distance; IQ-TREE not installed)"
    else:
        tree = parse_newick(nwk_txt)
    tree = midpoint_root(tree)
    names = {ids[i]: seqs[i][0] for i in range(len(seqs))}
    for lf in leaves(tree):
        lf.name = names.get(lf.name, lf.name)
    shutil.copy(aln, work / "16S_alignment.fasta")
    return to_newick(tree), method


# ══════════════════════════════════════════════════════════════════════════════
#  Species verdict
# ══════════════════════════════════════════════════════════════════════════════
def verdict(best: dict | None, best16: dict | None, others: list | None = None) -> dict:
    """Known / borderline / potential novel species from the closest type strain.

    others: the next type strains. When some of them are also ≥ 95 %, ANI alone
    does not separate those species (B. cereus / B. thuringiensis, the
    Geobacillus thermoleovorans group…) and the verdict says so.
    """
    v = _verdict(best, best16)
    def a(p):
        return p.get("anib") if p.get("anib") is not None else p.get("skani_ani")
    close = [p for p in (others or []) if a(p) is not None and a(p) >= SPECIES_ANI]
    if close and v.get("level") in ("known", "borderline"):
        names = ", ".join(f"{p['organism']} {a(p):.2f} %" for p in close[:4])
        v["close_relatives"] = [p["organism"] for p in close]
        v["text"] += (f" Also ≥ 95 % to {names}: these type strains are themselves inside "
                      f"one species boundary, so ANI ranks them but cannot separate them — "
                      f"the closest one is reported; dDDH (TYGS) or the published taxonomy "
                      f"of these names decides.")
        if v["level"] == "known":
            v["call"] = "known species (close relatives ≥ 95 %)"
    return v


def _verdict(best: dict | None, best16: dict | None) -> dict:
    if not best:
        v = {"call": "no type strain within ANI range",
             "text": "No sequenced type strain is within ~80 % ANI of this genome."}
        if best16 and best16["identity"] < RRNA_SPECIES:
            v["text"] += (f" Closest type-strain 16S: {best16['species']} "
                          f"({best16['identity']:.2f} %), below the 98.65 % species "
                          "threshold — a potential novel taxon.")
        v["level"] = "novel"
        return v
    ani = best.get("anib") if best.get("anib") is not None else best.get("skani_ani")
    tet = best.get("tetra")
    name = best["organism"]
    lab = "ANIb" if best.get("anib") is not None else "ANI (skani)"
    s = f"{lab} {ani:.2f} %" if ani is not None else "ANI n/a"
    if tet is not None:
        s += f", TETRA {tet:.4f}"
    st = clean_strain(best.get("strain") or "")
    ts = f"type strain {st} ({best['accession']})" if st else f"type strain ({best['accession']})"
    if ani is None:
        return {"call": "undetermined", "level": "unknown", "species": "",
                "text": f"ANI to {name} could not be computed."}
    if ani >= SPECIES_ANI_STRICT or (ani >= SPECIES_ANI and (tet or 0) >= TETRA_SPECIES):
        return {"call": "known species", "level": "known", "species": name,
                "text": f"{name} — {s} to the {ts}; above the 95–96 % species boundary."}
    if ani >= SPECIES_ANI:
        return {"call": "borderline (95–96 % ANI)", "level": "borderline", "species": name,
                "text": f"Most likely {name} ({s} to the {ts}) but inside the 95–96 % "
                        "transition zone — confirm with dDDH (TYGS, ≥ 70 % = same species)."}
    txt = (f"Potential novel species — closest type strain {name} ({s}), below the 95 % "
           "species boundary.")
    if best16:
        txt += (f" Best 16S match to a type strain: {best16['species']} "
                f"{best16['identity']:.2f} %"
                + (" (< 98.65 %, supports novelty)." if best16["identity"] < RRNA_SPECIES
                   else " (16S cannot separate close species)."))
    txt += (" Confirm with dDDH on TYGS (< 70 %) before describing a new species; a species "
            "whose type strain has no public genome would also look novel here.")
    return {"call": "potential novel species", "level": "novel",
            "species": f"{best['genus']} sp. (closest: {name})", "text": txt}


# ══════════════════════════════════════════════════════════════════════════════
#  Reference mapping
# ══════════════════════════════════════════════════════════════════════════════
def run_shell(cmd: str, tools: list, logf: Path, timeout: int = 4 * 3600) -> int:
    env = os.environ.copy()
    for t in tools:
        if t:
            env["PATH"] = f"{t.bindir}{os.pathsep}{env['PATH']}"
    H.log(f"    $ {cmd[:400]}")
    with open(logf, "a") as lf:
        lf.write(f"\n### {cmd}\n")
        lf.flush()
        try:
            p = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], stdout=lf, stderr=lf,
                               env=env, timeout=timeout)
            return p.returncode
        except subprocess.TimeoutExpired:
            return -9


def _q(p) -> str:
    return "'" + str(p).replace("'", "'\\''") + "'"


def parse_cs(cs: str) -> tuple[int, int, int]:
    """(SNPs, indel events, indel bases) from a minimap2 cs:Z: string."""
    snp = ins = indel_b = 0
    for op, val in re.findall(r"([:*+\-~])([0-9a-z]+)", cs):
        if op == "*":
            snp += 1
        elif op in "+-":
            ins += 1
            indel_b += len(val)
    return snp, ins, indel_b


def merge_iv(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def gaps_of(iv, length, min_len):
    out, pos = [], 0
    for s, e in merge_iv(iv):
        if s - pos >= min_len:
            out.append((pos, s))
        pos = max(pos, e)
    if length - pos >= min_len:
        out.append((pos, length))
    return out


def reference_mapping(reads: Path, asm: Path, ref: dict, outdir: Path, work: Path,
                      threads: int, logf: Path, keep_bam: bool = False) -> dict:
    mm, st = H.find_tool("minimap2"), H.find_tool("samtools")
    if not mm or not st:
        H.warn("minimap2 / samtools not installed — reference mapping skipped "
               "(run setup_wgs.sh)")
        return {}
    outdir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    refp = Path(ref["path"])
    ref_plain = _plain_fasta(refp, work / "ref.fna")
    ref_seqs = [(n.split()[0], len(s)) for n, s in H.read_fasta(ref_plain)]
    ref_len = sum(L for _, L in ref_seqs)
    res = {"reference": ref["organism"], "accession": ref["accession"],
           "strain": ref.get("strain", ""), "type_strain": bool(ref.get("type_strain")),
           "ref_length": ref_len, "ref_contigs": len(ref_seqs)}

    # ── reads → reference ───────────────────────────────────────────────────
    bam = work / "reads_vs_ref.bam"
    t = max(1, int(threads))
    rc = run_shell(f"minimap2 -ax map-ont -t {t} --secondary=no {_q(ref_plain)} {_q(reads)} "
                   f"| samtools sort -@ {min(4, t)} -m 768M -o {_q(bam)} - && "
                   f"samtools index {_q(bam)}", [st, mm], logf)
    if rc == 0 and bam.exists():
        cnt = work / "flagstat.txt"
        run_shell(f"samtools flagstat -@ {min(4, t)} {_q(bam)} > {_q(cnt)}", [st], logf)
        fs = cnt.read_text() if cnt.exists() else ""
        m_tot = re.search(r"(\d+) \+ \d+ in total", fs)
        m_map = re.search(r"(\d+) \+ \d+ primary mapped", fs) or \
            re.search(r"(\d+) \+ \d+ mapped", fs)
        m_prim = re.search(r"(\d+) \+ \d+ primary\b(?! mapped)", fs)
        tot = int(m_prim.group(1)) if m_prim else (int(m_tot.group(1)) if m_tot else 0)
        mapped = int(m_map.group(1)) if m_map else 0
        res["reads_total"] = tot
        res["reads_mapped"] = mapped
        res["reads_mapped_pct"] = round(100.0 * mapped / tot, 2) if tot else None
        # depth along the reference, binned
        win = max(1000, ref_len // 1500)
        depth_f = work / "depth.tsv"
        run_shell(f"samtools depth -a -@ {min(4, t)} {_q(bam)} > {_q(depth_f)}", [st], logf)
        bins, absent = [], []
        cov1 = cov10 = dsum = 0
        if depth_f.exists():
            cur, acc_d, acc_n, bstart = None, 0, 0, 0
            zrun_s, zero_min = None, 1000
            lastpos = 0
            with open(depth_f) as fh:
                for line in fh:
                    c, p, d = line.split("\t")
                    p, d = int(p), int(d)
                    if c != cur:
                        if cur is not None:
                            if acc_n:
                                bins.append((cur, bstart, acc_d / acc_n))
                            if zrun_s is not None and lastpos - zrun_s + 1 >= zero_min:
                                absent.append((cur, zrun_s, lastpos))
                        cur, acc_d, acc_n, bstart, zrun_s = c, 0, 0, 0, None
                    dsum += d
                    if d >= 1:
                        cov1 += 1
                    if d >= 10:
                        cov10 += 1
                    if d <= 2:
                        if zrun_s is None:
                            zrun_s = p
                    else:
                        if zrun_s is not None and p - zrun_s >= zero_min:
                            absent.append((c, zrun_s, p - 1))
                        zrun_s = None
                    acc_d += d
                    acc_n += 1
                    if acc_n >= win:
                        bins.append((c, bstart, acc_d / acc_n))
                        bstart, acc_d, acc_n = p, 0, 0
                    lastpos = p
            if cur is not None:
                if acc_n:
                    bins.append((cur, bstart, acc_d / acc_n))
                if zrun_s is not None and lastpos - zrun_s + 1 >= zero_min:
                    absent.append((cur, zrun_s, lastpos))
            depth_f.unlink()
        res.update({"mean_depth": round(dsum / ref_len, 1) if ref_len else None,
                    "breadth_1x": round(100.0 * cov1 / ref_len, 2) if ref_len else None,
                    "breadth_10x": round(100.0 * cov10 / ref_len, 2) if ref_len else None,
                    "absent_regions": len(absent),
                    "absent_bp": sum(e - s + 1 for _, s, e in absent),
                    "bin_size": win})
        res["coverage_bins"] = [[c, s, round(d, 1)] for c, s, d in bins]
        res["absent"] = [[c, s, e] for c, s, e in absent[:500]]
        with open(outdir / "reference_coverage_bins.tsv", "w") as fh:
            fh.write("contig\tstart\tmean_depth\n")
            for c, s, d in bins:
                fh.write(f"{c}\t{s}\t{d:.1f}\n")
        with open(outdir / "reference_absent_regions.csv", "w") as fh:
            fh.write("reference_contig,start,end,length\n")
            for c, s, e in absent:
                fh.write(f"{c},{s},{e},{e - s + 1}\n")
        if keep_bam:
            shutil.copy(bam, outdir / "reads_vs_reference.bam")
            shutil.copy(str(bam) + ".bai", outdir / "reads_vs_reference.bam.bai")
            shutil.copy(ref_plain, outdir / f"reference_{ref['accession']}.fasta")
    else:
        H.warn(f"read mapping to {ref['accession']} failed — see logs/mapping.log")

    # ── assembly → reference: SNPs / indels / synteny ───────────────────────
    paf = outdir / "assembly_vs_reference.paf"
    if H.run(mm, ["-cx", "asm20", "--cs", "-t", t, ref_plain, asm], logf,
             stdout_to=paf) == 0 and paf.exists():
        snp = ind = ind_b = 0
        r_iv, q_iv = defaultdict(list), defaultdict(list)
        qlens, blocks, matches, alnlen = {}, [], 0, 0
        for line in paf.read_text().splitlines():
            v = line.split("\t")
            if len(v) < 12 or "tp:A:P" not in line:
                continue
            qn, ql, qs, qe, strand, tn, tl, ts, te = v[0], int(v[1]), int(v[2]), int(v[3]), \
                v[4], v[5], int(v[6]), int(v[7]), int(v[8])
            qlens[qn] = ql
            r_iv[tn].append((ts, te))
            q_iv[qn].append((qs, qe))
            matches += int(v[9])
            alnlen += int(v[10])
            cs = next((x[5:] for x in v[12:] if x.startswith("cs:Z:")), "")
            a, b, c = parse_cs(cs)
            snp += a
            ind += b
            ind_b += c
            if te - ts >= 2000:
                blocks.append([tn, ts, te, qn, qs, qe, strand])
        ref_aln = sum(e - s for iv in r_iv.values() for s, e in merge_iv(iv))
        q_aln = sum(e - s for iv in q_iv.values() for s, e in merge_iv(iv))
        q_tot = sum(qlens.values()) or sum(len(s) for _, s in H.read_fasta(asm))
        uniq = []
        all_q = {n.split()[0]: len(s) for n, s in H.read_fasta(asm)}
        for qn, L in all_q.items():
            for s, e in gaps_of(q_iv.get(qn, []), L, 5000):
                uniq.append((qn, s + 1, e))
        res.update({"asm_ref_aligned_pct": round(100.0 * ref_aln / ref_len, 2)
                    if ref_len else None,
                    "asm_query_aligned_pct": round(100.0 * q_aln / sum(all_q.values()), 2)
                    if all_q else None,
                    "snps": snp, "indels": ind, "indel_bp": ind_b,
                    "snps_per_100kb": round(snp / (ref_aln / 1e5), 1) if ref_aln else None,
                    "alignment_identity": round(100.0 * matches / alnlen, 3) if alnlen else None,
                    "isolate_unique_regions": len(uniq),
                    "isolate_unique_bp": sum(e - s + 1 for _, s, e in uniq)})
        blocks.sort(key=lambda b: -(b[2] - b[1]))
        res["synteny_blocks"] = blocks[:3000]
        res["ref_seqs"] = ref_seqs
        res["query_seqs"] = sorted(all_q.items(), key=lambda x: -x[1])
        res["unique"] = [[c, s, e] for c, s, e in uniq[:500]]
        with open(outdir / "isolate_unique_regions.csv", "w") as fh:
            fh.write("contig,start,end,length\n")
            for c, s, e in uniq:
                fh.write(f"{c},{s},{e},{e - s + 1}\n")
    shutil.rmtree(work, ignore_errors=True)
    return res


# ══════════════════════════════════════════════════════════════════════════════
#  Orchestration for one sample
# ══════════════════════════════════════════════════════════════════════════════
def _genus_of(name: str) -> str:
    return re.sub(r"_[A-Z]+$", "", (name or "").split(" ")[0]) if name else ""


def _spname(name: str) -> str:
    """'Bacillus_A paranthracis' / 'Bacillus paranthracis strain X' → 'Bacillus paranthracis'."""
    w = re.sub(r"_[A-Z]+\b", "", name or "").split()
    return " ".join(w[:2]) if len(w) >= 2 else ""


def typestrain_analysis(S: str, asm: Path, s16: list[tuple[str, str]], gtdb: dict,
                        neighbours: list[str], args, dbp: dict, db_root: Path,
                        sdir: Path, wdir: Path, logs: Path, P) -> dict:
    """TYGS/JSpecies-style comparison of one assembly with type-strain genomes.

    gtdb       : {"genus", "species", "accession", "path"} of the closest GTDB rep
    neighbours : NCBI organism names of GTDB representatives sharing k-mers
    """
    cache = Path(args.ref_cache)
    tdir = sdir / "taxonomy"
    tdir.mkdir(parents=True, exist_ok=True)
    tw = wdir / "typestrain"
    tw.mkdir(parents=True, exist_ok=True)
    logf = logs / "typestrain.log"
    out: dict = {"type_strains": [], "rrna_type": [], "verdict": {}}

    # 1. 16S vs type-strain 16S
    db16 = ensure_16s_db(dbp, db_root)
    rows16 = search_16s(s16, db16, tw, args.threads, logf) if db16 else []
    sp16 = species_16s(rows16)
    out["rrna_type"] = sp16[:25]
    out["rrna_db"] = "NCBI RefSeq 16S (type material)" if db16 else ""

    # 2. candidate genera and type-strain genomes
    genera = []
    for g in [gtdb.get("genus", ""), _genus_of(sp16[0]["species"]) if sp16 else ""] + \
            [_genus_of(r["species"]) for r in sp16[:5]]:
        g = _genus_of(g)
        if g and g not in genera:
            genera.append(g)
    genera = genera[:2]
    P(0.745, f"Type strains of {', '.join(genera) or '—'} (NCBI)")
    ts = []
    for g in genera:
        ts += type_strains(g, cache)
    out["genera"] = genera
    if not ts and not gtdb.get("path"):
        H.warn(f"[{S}] no type-strain genomes could be listed (NCBI unreachable and no "
               "cache) — type-strain comparison skipped")
        return out
    # rank candidates when the genus is large: 16S neighbours and GTDB neighbours first
    want16 = {r["species"]: r["identity"] for r in sp16}
    wantg = {_spname(n) for n in neighbours if n}
    if gtdb.get("species"):
        wantg.add(_spname(gtdb["species"]))

    def prio(r):
        base = _spname(r["organism"])
        return (0 if base in wantg else 1, -want16.get(r["organism"], want16.get(base, 0)),
                r["organism"])
    ts.sort(key=prio)
    cand = ts[: max(5, int(args.ts_download_max))]
    P(0.75, f"Downloading {len(cand)} type-strain genome(s)")
    with ThreadPoolExecutor(max_workers=3) as ex:
        paths = list(ex.map(lambda r: H.fetch_genome(r["accession"], cache), cand))
    pool = []
    for r, p in zip(cand, paths):
        if p:
            pool.append(dict(r, path=str(p), type_strain=True))
    if gtdb.get("path") and not any(p["accession"].split(".")[0][4:] ==
                                    gtdb["accession"].split(".")[0][4:] for p in pool):
        pool.append({"accession": gtdb["accession"], "organism": gtdb.get("ncbi_name") or
                     gtdb.get("species", ""), "strain": "", "genus": _genus_of(
                         gtdb.get("species", "")), "path": gtdb["path"], "type_strain": False,
                     "level": "", "note": "GTDB species representative"})
    if not pool:
        H.warn(f"[{S}] type-strain genomes could not be downloaded — comparison skipped")
        return out

    # 3. skani screen → closest N
    sk = H.find_tool("skani")
    P(0.76, f"ANI screen against {len(pool)} genome(s) (skani)")
    if sk:
        skf = tw / "skani_screen.tsv"
        H.run(sk, ["dist", "-q", asm, "-r"] + [p["path"] for p in pool] +
              ["-t", args.threads, "-o", skf], logf)
        by_path = {p["path"]: p for p in pool}
        for line in (skf.read_text().splitlines()[1:] if skf.exists() else []):
            v = line.split("\t")
            if len(v) >= 5 and v[0] in by_path:
                by_path[v[0]].update({"skani_ani": float(v[2]), "skani_af_ref": float(v[3]),
                                      "skani_af_query": float(v[4])})
    hits = sorted([p for p in pool if p.get("skani_ani") is not None],
                  key=lambda p: (-p["skani_ani"], -p.get("skani_af_query", 0)))
    top = [p for p in hits if p.get("type_strain")][: int(args.ts_max)]
    if hits:
        # closest genome; among those within 0.3 % ANI prefer the most complete
        # assembly (a 100-contig reference makes coverage and synteny hard to read)
        near = [p for p in hits if p["skani_ani"] >= hits[0]["skani_ani"] - 0.3]
        near.sort(key=lambda p: (LEVEL_RANK.get(p.get("level", ""), 5), -p["skani_ani"]))
        out["_mapping_ref"] = near[0]          # (path inside; popped by the pipeline)

    # 4. ANIb / ANIm / TETRA
    if top and not args.skip_anib:
        P(0.77, f"ANIb / ANIm / TETRA vs {len(top)} type strain(s)")
        js = jspecies(asm, top, tw / "jspecies", args.threads, logf, wdir / "jspecies_cache.json")
        for p in top:
            p.update({k: v for k, v in js.get(p["accession"], {}).items()
                      if k not in ("v", "accession")})
        top.sort(key=lambda p: (-(p.get("anib") if p.get("anib") is not None
                                  else p.get("skani_ani", 0))))
    for p in top:
        p["rrna_identity"] = want16.get(p["organism"]) or want16.get(_spname(p["organism"]))
    out["type_strains"] = [{k: v for k, v in p.items() if k != "path"} for p in top]
    best = top[0] if top else None
    best16 = sp16[0] if sp16 else None
    out["verdict"] = verdict(best, best16, top[1:])
    out["verdict"]["rrna_best"] = best16

    # 5. trees — a failure here must not cost the ANI table above
    if not args.skip_tree:
        try:
            _trees(S, asm, s16, top, rows16, sp16, db16, out, tdir, tw, args, logf, P)
        except Exception as e:
            import traceback
            fr = traceback.extract_tb(e.__traceback__)[-1]
            H.warn(f"[{S}] tree step failed: {e} (wgs_taxonomy.py:{fr.lineno} {fr.name}) — "
                   "the type-strain table is unaffected")
    _write_ts_tables(tdir, out, sp16)
    if not args.keep_intermediate:
        shutil.rmtree(tw, ignore_errors=True)
    return out


def _trees(S, asm, s16, top, rows16, sp16, db16, out, tdir, tw, args, logf, P):
    if True:
        P(0.78, "Genome tree (skani + FastME)")
        gl = [(f"{S} (this isolate)", asm)] + [
            (ts_label(p["organism"], p.get("strain", ""), p["accession"]),
             Path(p["path"])) for p in top]
        nwk, mat = genome_tree(gl, tw / "gtree", args.threads, logf)
        if nwk:
            (tdir / f"{S}_genome_tree.nwk").write_text(nwk + "\n")
            out["genome_tree"] = nwk
            out["genome_tree_method"] = ("FastME 2 (BioNJ + SPR)" if H.find_tool("fastme")
                                         else "neighbour joining") + \
                " on skani ANI distances (1 − ANI/100), midpoint-rooted"
            with open(tdir / "ani_matrix_skani.tsv", "w") as fh:
                fh.write("genome\t" + "\t".join(n for n, _ in gl) + "\n")
                for (n, _), row in zip(gl, mat):
                    fh.write(n + "\t" + "\t".join("" if x is None else f"{x:.2f}"
                                                   for x in row) + "\n")
        # 16S tree: isolate copies (distinct) + one type-strain record per species
        P(0.79, "16S tree (MAFFT + IQ-TREE)")
        iso = []
        seen = set()
        for n, s in s16:
            if s.upper() not in seen:
                seen.add(s.upper())
                iso.append((f"{S} 16S copy {len(iso) + 1}", s))
        iso = iso[:3]
        sp_for_tree = []
        for p in top:
            r = next((x for x in rows16 if x["species"] == p["organism"] or
                      x["species"] == _spname(p["organism"])), None)
            if r and r["accession"] not in {x["accession"] for x in sp_for_tree}:
                sp_for_tree.append(r)
        top16 = sp16[0]["identity"] if sp16 else 0
        for r in sp16:
            if len(sp_for_tree) >= 16 or r["identity"] < top16 - 3.0:
                break
            if r["accession"] not in {x["accession"] for x in sp_for_tree}:
                sp_for_tree.append(r)
        recs = fetch_16s_records(db16, {r["accession"] for r in sp_for_tree}) if db16 else {}
        seqs = iso + [(ts_label(r["species"], r.get("strain", ""), r["accession"]),
                       recs[r["accession"]]) for r in sp_for_tree if r["accession"] in recs]
        if len(seqs) >= 4:
            nwk16, meth = rrna_tree(seqs, tw / "rtree", args.threads, logf)
            if nwk16:
                (tdir / f"{S}_16S_tree.nwk").write_text(nwk16 + "\n")
                aln = tw / "rtree" / "16S_alignment.fasta"
                if aln.exists():
                    names = {f"s{i}": n for i, (n, _) in enumerate(seqs)}
                    H.write_fasta(tdir / f"{S}_16S_alignment.fasta",
                                  [(names.get(re.sub(r"^_R_", "", n), n), s)
                                   for n, s in H.read_fasta(aln)])
                out["rrna_tree"] = nwk16
                out["rrna_tree_method"] = "MAFFT (--auto) alignment, " + meth + \
                    ", midpoint-rooted"


def _write_ts_tables(tdir, out, sp16):
    with open(tdir / "type_strain_ani.csv", "w") as fh:
        cols = ["organism", "strain", "accession", "level", "anib", "anib_qr", "anib_rq",
                "anib_cov_q", "anib_cov_r", "anim", "anim_cov_q", "anim_cov_r", "tetra",
                "skani_ani", "skani_af_query", "skani_af_ref", "rrna_identity"]
        fh.write(",".join(cols) + "\n")
        for p in out["type_strains"]:
            fh.write(",".join(_csv(p.get(c)) for c in cols) + "\n")
    with open(tdir / "type_strain_16S.csv", "w") as fh:
        fh.write("species,strain,accession,identity,aln_len\n")
        for r in sp16:
            fh.write(",".join(_csv(r.get(c)) for c in
                              ("species", "strain", "accession", "identity", "aln_len")) + "\n")


def _csv(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4f}".rstrip("0").rstrip(".")
    s = str(v)
    return '"' + s.replace('"', "'") + '"' if ("," in s or '"' in s) else s
