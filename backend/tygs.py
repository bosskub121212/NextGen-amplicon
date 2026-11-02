#!/usr/bin/env python3
"""Import a TYGS (Type Strain Genome Server, DSMZ) result into one ONT-WGS sample.

TYGS has no submission API: the BI uploads <S>_assembly.fasta at https://tygs.dsmz.de
and pastes the result link / GUID here. This module reads TYGS's result API by GUID
(https://tygs.dsmz.de/application_programming_interface/show):

  status        /user_requests/show/status/<GUID>              (finalization_date set = done)
  dDDH          /user_results/show/digital_ddh_json/<GUID>
  ident. table  /user_results/show/identification_table/<GUID> (kept raw when it answers)
  strains       /user_results/show/strain_table_json/<GUID>
  literature    /user_results/show/bibtex/<GUID>
  trees         /user_phylogenies/{16S,genome}/newick_by_guid/<GUID>

TYGS asks clients to check a job's status at most once or twice an hour: an unfinished
job is re-checked only after MIN_RECHECK seconds; until then the cached state is returned.

A TYGS job may hold several user genomes (e.g. a whole run, other customers). Only this
sample's genome is kept in the deliverable files: dDDH rows of this genome vs type strains,
trees with every other user genome pruned. The raw answers stay in tygs/.raw/ (never
packaged).

Stored in samples/<S>/tygs/: tygs.json (what the report reads), dddh.tsv,
genome_tree.nwk, 16S_tree.nwk, state.json; .raw/ (raw API answers).
"""
from __future__ import annotations

import csv
import html
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

BASE = "https://tygs.dsmz.de"
GUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", re.I)
MIN_RECHECK = 1800          # s — TYGS: check a job's status at most once or twice an hour
SPECIES_D4, SUBSPECIES_D4 = 70.0, 79.0      # TYGS / GBDP delineation thresholds (formula d4)


class TygsError(RuntimeError):
    pass


class TygsAmbiguous(TygsError):
    def __init__(self, names):
        super().__init__("this TYGS job holds several genomes — pick the one of this sample: " + ", ".join(names))
        self.names = list(names)


def parse_guid(text: str) -> str | None:
    m = GUID_RE.search(text or "")
    return m.group(0).lower() if m else None


def http_get(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": "NextGen-Amplicon (TYGS result import)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.headers.get("Content-Type", ""), r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return e.code, e.headers.get("Content-Type", "") if e.headers else "", ""


def plain(s) -> str:
    """TYGS HTML snippets → text: tags out, entities decoded, filename quotes/extension off."""
    t = html.unescape(re.sub(r"<[^>]+>", "", str(s or ""))).strip()
    t = t.strip("'\" ").strip()
    return re.sub(r"\.(fna|fasta|fa|fas|fsa)$", "", t, flags=re.I)


def species_of(s) -> str:
    it = re.findall(r"<i>(.*?)</i>", str(s or ""), flags=re.I)
    return " ".join(html.unescape(x).strip() for x in it[:3] if x.strip()) or " ".join(plain(s).split()[:2])


def conclusion(d4: float) -> str:
    if d4 >= SUBSPECIES_D4:
        return "belongs to known species"
    if d4 >= SPECIES_D4:
        return "belongs to known species, potential new subspecies"
    return "potential new species"


def _f(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def _fmt(v) -> str:
    x = _f(v)
    return "" if x is None else (f"{x:.1f}" if abs(x) >= 1 or x == 0 else f"{x:g}")


# ── trees ─────────────────────────────────────────────────────────────────────
def norm(s) -> str:
    """Compare TYGS labels across its tables and trees: case, _ vs space, quotes, extension."""
    return re.sub(r"\s+", " ", plain(s).replace("_", " ")).strip().lower()


def _prune(newick: str, keep_user: str, allow: set, rename_to: str) -> str:
    """Keep only this genome and the type strains TYGS compared (allow-list): any other leaf —
    another genome of the same TYGS job, or a label we cannot place — is removed."""
    import wgs_taxonomy as WT
    root = WT.parse_newick(html.unescape(newick))
    for lf in list(WT.leaves(root)):
        nm = norm(lf.name)
        if nm == norm(keep_user):
            lf.name = rename_to
        elif nm not in allow:
            p = lf.parent
            if p is None:
                continue
            p.children.remove(lf)
            while p is not None and len(p.children) == 1 and p.parent is not None:  # collapse unary nodes
                c = p.children[0]
                c.length += p.length
                c.parent = p.parent
                p.parent.children[p.parent.children.index(p)] = c
                p = c.parent
            while root.children and len(root.children) == 1 and root.children[0].children:
                root = root.children[0]
                root.parent = None
    if sum(1 for _ in WT.leaves(root)) < 3:
        raise ValueError("too few leaves left")
    return WT.to_newick(root)


# ── import ────────────────────────────────────────────────────────────────────
def _json(get, url):
    try:
        st, ct, txt = get(url, timeout=90)
    except OSError as e:
        raise TygsError(f"could not reach TYGS ({e}) — check the internet connection and try again")
    if st != 200:
        return st, None
    try:
        return st, json.loads(txt)
    except ValueError:
        return st, None


def _text(get, url):
    try:
        st, ct, txt = get(url, timeout=90)
    except OSError:
        return None
    if st != 200 or not txt.strip() or txt.lstrip().lower().startswith(("<!doctype", "<html")):
        return None
    return txt


def _cache_get(cache: Path | None, guid: str):
    try:
        return json.loads((Path(cache) / f"{guid}.json").read_text()) if cache else None
    except (OSError, ValueError):
        return None


def _cache_put(cache: Path | None, guid: str, checked_at: float, js):
    if cache:
        Path(cache).mkdir(parents=True, exist_ok=True)
        (Path(cache) / f"{guid}.json").write_text(json.dumps({"checked_at": checked_at, "status": js}))


def _clear_result(dest: Path):
    for f in ("tygs.json", "dddh.tsv", "16S_tree.nwk", "genome_tree.nwk"):
        (dest / f).unlink(missing_ok=True)


def import_result(guid_text: str, sample: str, dest: Path, genome: str | None = None,
                  get=http_get, now=time.time, cache: Path | None = None, force: bool = False) -> dict:
    """Import (or re-check) a TYGS job for one sample. Returns the state dict:
    {"status": "pending"|"done", "guid", "genome", "checked_at", "next_check", ...}."""
    guid = parse_guid(guid_text)
    if not guid:
        raise TygsError("not a TYGS result link or GUID (it looks like 29ae139a-3a41-4571-a19b-809af76a3ec5)")
    dest = Path(dest)
    sp = dest / "state.json"
    try:
        old = json.loads(sp.read_text())
    except (OSError, ValueError):
        old = {}
    t = now()
    if (not force and old.get("guid") == guid and old.get("status") == "done" and (dest / "tygs.json").exists()
            and (genome is None or genome == old.get("genome"))):
        return dict(old, cached=True)                  # already imported: TYGS is not asked again
    if old.get("guid") and old.get("guid") != guid:
        _clear_result(dest)                            # a new link: the old result must not stay in the report
    # one status check per GUID per MIN_RECHECK, whichever sample asks (a run in one TYGS job)
    c = _cache_get(cache, guid)
    if c and c.get("status") and (c["status"][0] or {}).get("finalization_date"):
        js = c["status"]
    elif c and t - float(c.get("checked_at") or 0) < MIN_RECHECK:
        js = c["status"]
    elif (old.get("guid") == guid and old.get("status") == "pending"
            and t - float(old.get("checked_at") or 0) < MIN_RECHECK):
        return dict(old, cached=True)
    else:
        st, js = _json(get, f"{BASE}/user_requests/show/status/{guid}")
        if st == 404 or (st == 200 and not js):
            raise TygsError("TYGS job not found — check the result link / GUID")
        if st != 200 or not isinstance(js, list):
            raise TygsError(f"TYGS answered HTTP {st} — try again later")
        _cache_put(cache, guid, t, js)
    info = js[0] if js else {}
    raw = dest / ".raw"
    raw.mkdir(parents=True, exist_ok=True)
    (raw / "status.json").write_text(json.dumps(js, indent=1))
    if not info.get("finalization_date"):
        checked = float((c or {}).get("checked_at") or t) if c and js is c.get("status") else t
        state = {"status": "pending", "guid": guid, "checked_at": checked, "next_check": checked + MIN_RECHECK,
                 "message": "TYGS has not finished this job yet — it usually takes a few hours"}
        sp.write_text(json.dumps(state, indent=1))
        return state

    st, dd = _json(get, f"{BASE}/user_results/show/digital_ddh_json/{guid}")
    if not isinstance(dd, list) or not dd:
        raise TygsError("the TYGS job is finished but its dDDH table could not be read — try again later")
    (raw / "digital_ddh.json").write_text(json.dumps(dd, indent=1))

    users = []
    for r in dd:
        ct = str(r.get("comp_type", ""))
        if ct.startswith("U"):
            users.append(plain(r.get("query_genome")))
        if ct.replace(" ", "").upper().endswith("VS.U"):
            users.append(plain(r.get("subject_genome")))
    for g in re.split(r",\s*(?=['\w])", str(info.get("genomes_list") or "")):   # every genome the job holds
        g = re.sub(r"\s*\([A-Z]{2,6}\d+\)\s*$", "", g.strip())                   # "'name' (RZGQ00000000)"
        if g:
            users.append(plain(g))
    users = list(dict.fromkeys(u for u in users if u))
    import wgs_package as WP
    want = {f"{WP.safe(sample)}_assembly", WP.safe(sample), f"{sample}_assembly", sample}
    if genome:
        if genome not in users:
            raise TygsAmbiguous(users)
        me = genome
    else:
        wn = {norm(w) for w in want}
        hits = list(dict.fromkeys(u for u in users if norm(u) in wn))
        if len(hits) == 1:
            me = hits[0]
        else:                       # also a one-genome job under another name: maybe a wrong link — ask
            raise TygsAmbiguous(users)
    label = f"{sample}_assembly"

    rows = []
    for r in dd:
        q, s_ = plain(r.get("query_genome")), plain(r.get("subject_genome"))
        ct = str(r.get("comp_type", "")).replace(" ", "").upper()
        if not ct.endswith("VS.T"):
            continue
        if norm(q) != norm(me):
            continue
        rows.append({"query": label, "subject": s_, "species": species_of(r.get("subject_genome")),
                     "d0": _fmt(r.get("digital_ddh_d0")), "ci_d0": str(r.get("confidence_interval_d0") or ""),
                     "d4": _fmt(r.get("digital_ddh_d4")), "ci_d4": str(r.get("confidence_interval_d4") or ""),
                     "d6": _fmt(r.get("digital_ddh_d6")), "ci_d6": str(r.get("confidence_interval_d6") or ""),
                     "gc_diff": str(r.get("percent_gc_diff") or "")})
    rows.sort(key=lambda r: -(_f(r["d4"]) or 0))
    if rows:
        b = rows[0]
        d4 = _f(b["d4"]) or 0
        ident = [{"strain": label, "conclusion": conclusion(d4),
                  "species": b["species"] if d4 >= SPECIES_D4 else "",
                  "remark": f"closest type strain {b['subject']}: dDDH d4 {b['d4']} % {b['ci_d4']} "
                            f"(species ≥ {SPECIES_D4:g} %, subspecies ≥ {SUBSPECIES_D4:g} %)"}]
    else:
        ident = [{"strain": label, "conclusion": "no type strain compared",
                  "species": "", "remark": "TYGS reported no comparison of this genome with a type strain"}]

    allow = {norm(r.get("subject_genome")) for r in dd
             if str(r.get("comp_type", "")).replace(" ", "").upper().endswith("VS.T")}
    allow.discard("")
    trees = {}
    for key in ("16S", "genome"):
        txt = _text(get, f"{BASE}/user_phylogenies/{key}/newick_by_guid/{guid}")
        if txt and "(" in txt:
            (raw / f"{key}_tree.nwk").write_text(txt)
            try:
                trees[key] = _prune(txt.strip(), me, allow, sample)
            except Exception:
                pass
    for name, url in (("identification_table.txt", f"/user_results/show/identification_table/{guid}"),
                      ("strain_table.json", f"/user_results/show/strain_table_json/{guid}"),
                      ("literature.bib", f"/user_results/show/bibtex/{guid}")):
        txt = _text(get, BASE + url)
        if txt:
            (raw / name).write_text(txt)

    # deliverable files (this genome only)
    for f in ("dddh.tsv", "16S_tree.nwk", "genome_tree.nwk"):
        (dest / f).unlink(missing_ok=True)
    with open(dest / "dddh.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Query", "Subject", "d0", "C.I. d0", "d4", "C.I. d4", "d6", "C.I. d6", "Diff. G+C %"])
        for r in rows:
            w.writerow([r[k] for k in ("query", "subject", "d0", "ci_d0", "d4", "ci_d4", "d6", "ci_d6", "gc_diff")])
    for key, nwk in trees.items():
        (dest / f"{key}_tree.nwk").write_text(nwk + "\n")
    data = {"guid": guid, "genome": me, "sample": sample, "finalization_date": info.get("finalization_date"),
            "imported_at": t, "identification": ident, "dddh": rows, "trees": trees,
            "result_link": f"{BASE}/user_results/show?guid={guid}"}
    (dest / "tygs.json").write_text(json.dumps(data, indent=1))
    state = {"status": "done", "guid": guid, "genome": me, "checked_at": t,
             "finalization_date": info.get("finalization_date")}
    sp.write_text(json.dumps(state, indent=1))
    return state


def load(d: Path) -> dict | None:
    try:
        t = json.loads((Path(d) / "tygs.json").read_text())
        return t if t.get("identification") else None
    except (OSError, ValueError):
        return None


def state(d: Path) -> dict:
    try:
        return json.loads((Path(d) / "state.json").read_text())
    except (OSError, ValueError):
        return {}
