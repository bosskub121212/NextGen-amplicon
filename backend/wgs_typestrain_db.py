#!/usr/bin/env python3
"""wgs_typestrain_db.py — the type-strain genomes of NCBI, downloaded ahead of time,
so the ONT-WGS type-strain comparison (ANIb / ANIm / TETRA, genome tree,
reference mapping) runs without an internet connection.

The pipeline already caches what it downloads in ~/.ngamp_wgs_refs:
    genus_taxids.json                 genus name -> NCBI taxid
    type_strains/<taxid>.v3.json      one genome per (sub)species of that genus
    <accession>.fna.gz                the genomes
This script fills that same cache for every genus (or the genera you name),
choosing genomes exactly as the pipeline does (wgs_taxonomy.pick_type_strains),
so an offline run finds everything where it would have put it itself.

    python3 wgs_typestrain_db.py build              # all bacteria + archaea  (~30-40 GB)
    python3 wgs_typestrain_db.py build --genera Bacillus,Lactiplantibacillus
    python3 wgs_typestrain_db.py build --dry-run    # count and size only, nothing downloaded
    python3 wgs_typestrain_db.py status

Safe to stop and re-run: genomes already in the cache are skipped. Lists older
than 90 days are refreshed by the pipeline whenever NCBI is reachable; offline,
a list of any age is used.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wgs_taxonomy as WT          # noqa: E402

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
UA = {"User-Agent": "NextGen-Amplicon (type-strain offline cache)"}
DOMAINS = {"Bacteria": 2, "Archaea": 2157}


def human(b: float) -> str:
    return f"{b / 1e9:.1f} GB" if b >= 1e9 else f"{b / 1e6:.0f} MB"


def say(*a):
    print(*a, flush=True)


WT.bind(log=lambda m: say(m), warn=lambda m: say("  ! " + m))


def get_json(path: str, timeout: int = 120, data: dict | None = None):
    url = API + path
    for attempt in range(6):
        try:
            req = urllib.request.Request(
                url, headers=dict(UA, **({"Content-Type": "application/json"} if data else {})),
                data=json.dumps(data).encode() if data else None)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < 5:
                time.sleep(2 * (attempt + 1))
                continue
            raise
        except Exception:
            if attempt < 5:
                time.sleep(3 * (attempt + 1))
                continue
            raise


def list_type_material(taxon: int | str) -> list[dict]:
    """Every 'assembly from type material' report under a taxon (same filter as
    wgs_taxonomy.type_strains)."""
    reps, token, page = [], "", 0
    while True:
        q = (f"/genome/taxon/{taxon}/dataset_report?filters.is_type_material=true"
             f"&filters.exclude_atypical=true&page_size=1000"
             + (f"&page_token={urllib.parse.quote(token)}" if token else ""))
        d = get_json(q)
        reps += d.get("reports", [])
        page += 1
        token = d.get("next_page_token") or ""
        say(f"    {taxon}: page {page}, {len(reps):,} assemblies")
        if not token:
            return reps


def genus_of_taxids(taxids: list[int]) -> dict[int, tuple[int, str]]:
    """species/strain taxid -> (genus taxid, genus name), 200 per request."""
    out = {}
    ids = sorted(set(t for t in taxids if t))
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        d = get_json("/taxonomy/taxon/" + ",".join(map(str, chunk)) + "/dataset_report")
        for r in (d or {}).get("reports", []):
            t = r.get("taxonomy", {})
            g = (t.get("classification") or {}).get("genus") or {}
            if t.get("tax_id") and g.get("id"):
                out[int(t["tax_id"])] = (int(g["id"]), g.get("name", ""))
        if (i // 200) % 10 == 0:
            say(f"    genus lookup {min(i + 200, len(ids)):,} / {len(ids):,}")
    return out


def download_batch(accs: list[str], cache: Path) -> tuple[int, list[str]]:
    """One NCBI Datasets zip for up to ~50 genomes, written as <acc>.fna.gz."""
    body = {"accessions": accs, "include_annotation_type": ["GENOME_FASTA"]}
    req = urllib.request.Request(API + "/genome/download", data=json.dumps(body).encode(),
                                 headers=dict(UA, **{"Content-Type": "application/json"}))
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=900) as r:
                data = r.read()
            break
        except Exception:
            if attempt == 4:
                return 0, accs
            time.sleep(5 * (attempt + 1))
    got = set()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for n in zf.namelist():
            parts = n.split("/")
            if len(parts) >= 4 and parts[-1].endswith((".fna", ".fasta", ".fa")):
                acc = parts[-2]
                out = cache / f"{acc}.fna.gz"
                tmp = out.with_suffix(".gz.part")
                with zf.open(n) as src, gzip.open(tmp, "wb", compresslevel=6) as dst:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        dst.write(b)
                tmp.replace(out)
                got.add(acc)
    return len(got), [a for a in accs if a not in got]


def build(a):
    cache = Path(a.cache).expanduser()
    cache.mkdir(parents=True, exist_ok=True)
    wanted = [g.strip() for g in (a.genera or "").split(",") if g.strip()]

    say("1. Listing NCBI 'assembly from type material'")
    reps = []
    if wanted:
        for g in wanted:
            tid, gname = WT.genus_taxid(g, cache)
            if not tid:
                say(f"  ! genus {g} not found in NCBI taxonomy — skipped")
                continue
            for r in list_type_material(tid):
                r["_genus"] = (tid, gname)
                reps.append(r)
    else:
        for dom, tid in DOMAINS.items():
            reps += list_type_material(tid)
        say("   genus of every species (NCBI taxonomy)")
        g = genus_of_taxids([int((r.get("organism") or {}).get("tax_id") or 0) for r in reps])
        for r in reps:
            r["_genus"] = g.get(int((r.get("organism") or {}).get("tax_id") or 0))

    by_genus: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in reps:
        if r.get("_genus"):
            by_genus[r["_genus"]].append(r)

    say(f"2. Choosing one genome per species in {len(by_genus):,} genera")
    known = {}
    gf = cache / "genus_taxids.json"
    if gf.exists():
        try:
            known = json.loads(gf.read_text())
        except Exception:
            known = {}
    todo, total_len, n_genomes = [], 0, 0
    for (gtid, gname), rows in sorted(by_genus.items(), key=lambda kv: kv[0][1]):
        picked = WT.pick_type_strains(rows, gname)
        if not picked:
            continue
        n_genomes += len(picked)
        if not a.dry_run:
            known[gname] = gtid
            p = cache / "type_strains" / f"{gtid}.v3.json"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(picked, indent=1))
        for rec in picked:
            acc = rec["accession"]
            f = cache / f"{acc}.fna.gz"
            if not (f.exists() and f.stat().st_size > 1000):
                todo.append(acc)
                total_len += rec.get("length", 0)
    if not a.dry_run:
        gf.write_text(json.dumps(known, indent=1))
    est = total_len * 0.30                   # gzip of bacterial DNA ≈ 30 % of its length
    say(f"   {n_genomes:,} type-strain genomes in {len(by_genus):,} genera; "
        f"{len(todo):,} still to download (≈ {human(est)} compressed)")
    if a.dry_run or not todo:
        return 0

    say(f"3. Downloading {len(todo):,} genomes, {a.batch} per request, {a.jobs} at a time")
    batches = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
    done, failed, t0 = 0, [], time.time()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(download_batch, b, cache) for b in batches]
        for k, fu in enumerate(as_completed(futs), 1):
            n, miss = fu.result()
            done += n
            failed += miss
            rate = done / max(1.0, time.time() - t0)
            eta = (len(todo) - done - len(failed)) / rate / 60 if rate else 0
            say(f"   {done:,} / {len(todo):,} genomes"
                + (f", {len(failed)} failed" if failed else "")
                + f" · ~{eta:.0f} min left")
    (cache / "offline_db.json").write_text(json.dumps({
        "built": time.strftime("%Y-%m-%d %H:%M"), "genera": len(by_genus),
        "genomes": n_genomes, "scope": wanted or "all bacteria + archaea",
        "failed": failed[:200]}, indent=1))
    if failed:
        say(f"  ! {len(failed)} genome(s) could not be downloaded — run the same command again "
            "to retry only those")
        return 1
    say("   done — the type-strain comparison now runs without an internet connection")
    return 0


def status(a):
    cache = Path(a.cache).expanduser()
    if not cache.exists():
        say(f"no cache at {cache}")
        return 1
    lists = list((cache / "type_strains").glob("*.v3.json"))
    wanted, have = 0, 0
    for p in lists:
        try:
            for r in json.loads(p.read_text()):
                wanted += 1
                f = cache / f"{r['accession']}.fna.gz"
                have += f.exists() and f.stat().st_size > 1000
        except Exception:
            pass
    size = sum(f.stat().st_size for f in cache.glob("*.fna.gz"))
    say(f"cache          {cache}")
    say(f"genera listed  {len(lists):,}")
    say(f"genomes        {have:,} of {wanted:,} on disk ({human(size)})")
    meta = cache / "offline_db.json"
    if meta.exists():
        m = json.loads(meta.read_text())
        say(f"last build     {m.get('built')} · {m.get('scope')}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--genera", default="", help="comma-separated genus names; default: all")
    b.add_argument("--dry-run", action="store_true")
    b.add_argument("--batch", type=int, default=40)
    b.add_argument("--jobs", type=int, default=3)
    s = sub.add_parser("status")
    for p in (b, s):
        p.add_argument("--cache", default=str(Path.home() / ".ngamp_wgs_refs"))
    a = ap.parse_args(argv)
    return build(a) if a.cmd == "build" else status(a)


if __name__ == "__main__":
    sys.exit(main())
