#!/usr/bin/env bash
# ============================================================
#  NextGen-Amplicon — Reference Database Downloader
#
#  Default set (current releases, checked October 2026):
#    1. SILVA 144          16S / 18S  — DADA2 files published by SILVA (Aug 2026)
#    2. UNITE 2025-02-19   ITS fungi  — general FASTA release v10.0
#    3. PR2 5.1.1          18S eukaryotes (DADA2 format)
#       18S-NemaBase       nematodes (GitHub WormsEtAl, rewritten as a DADA2 trainset)
#    4. MIDORI2 GB273      COI and 12S (mitochondrial srRNA), DADA2 format
#    5. Tax4Fun2           16S functional prediction
#  On request:
#    --silva138  SILVA 138.2 (Nov 2024) — for comparison with published work
#    --gtdb      GTDB R232 16S trainset, built from GTDB's own SSU files
#    --unite-euk UNITE all-eukaryotes release (ITS of non-fungi too)
#
#  Every file is checked against the publisher's MD5 where one is published,
#  registered in backend/databases/db_paths.json (other keys are kept), and
#  recorded with its version and checksum in backend/databases/db_manifest.json,
#  which provenance.py uses so every report can say which file it ran on.
#  Existing files are never deleted; re-running only fetches what is missing.
#
#  Usage:  bash ~/r16s-app/download_databases.sh [--silva138] [--gtdb] [--unite-euk]
#                                                [--only silva,unite,pr2,midori,tax4fun]
# ============================================================
set -uo pipefail   # no -e: one failed download must not abort the rest

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${DB_DIR:-$APP_DIR/backend/databases}"
PY="${PYTHON:-$APP_DIR/venv/bin/python}"
[[ -x "$PY" ]] || PY="$(command -v python3)"
mkdir -p "$DB_DIR"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "  ${GREEN}✔ $*${NC}"; }
warn() { echo -e "  ${YELLOW}⚠ $*${NC}"; }
fail() { echo -e "  ${RED}✘ $*${NC}"; }
hdr()  { echo; echo -e "${BOLD}$*${NC}"; }

WANT_SILVA138=0; WANT_GTDB=0; WANT_UNITE_EUK=0; ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --silva138)  WANT_SILVA138=1 ;;
    --gtdb)      WANT_GTDB=1 ;;
    --unite-euk) WANT_UNITE_EUK=1 ;;
    --only)      ONLY="${2:-}"; shift ;;
    -h|--help)   sed -n '2,27p' "$0"; exit 0 ;;
    *) warn "unknown option $1" ;;
  esac
  shift
done
want() { [[ -z "$ONLY" || ",$ONLY," == *",$1,"* ]]; }

echo "============================================================"
echo "  NextGen-Amplicon — Reference Database Downloader"
echo "  Target: $DB_DIR"
echo "============================================================"

# fetch <dest> <md5 or -> <url> [url2 ...]   — resumes, retries, verifies
fetch() {
  local dest="$1" md5="$2"; shift 2
  mkdir -p "$(dirname "$dest")"
  if [[ -s "$dest" ]]; then
    if [[ "$md5" == "-" ]] || [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
      ok "have $(basename "$dest")"; return 0
    fi
    warn "$(basename "$dest") is incomplete or changed — downloading again"
    mv "$dest" "$dest.bad.$(date +%s)"
  fi
  local url
  for url in "$@"; do
    echo "    $url"
    if curl -fL -# --retry 5 --retry-delay 5 --connect-timeout 30 -C - -A "Mozilla/5.0" \
            -o "$dest.part" "$url"; then
      if [[ "$md5" != "-" ]] && [[ "$(md5sum "$dest.part" | cut -d' ' -f1)" != "$md5" ]]; then
        fail "checksum mismatch for $(basename "$dest") — discarded"
        rm -f "$dest.part"; continue
      fi
      mv "$dest.part" "$dest"; ok "$(basename "$dest")"; return 0
    fi
  done
  rm -f "$dest.part"
  fail "$(basename "$dest") — every source failed"
  return 1
}

# register <db_paths key or -> <path> <label> <version> <source> <md5 or ->
register() {
  "$PY" - "$DB_DIR" "$@" <<'PYEOF'
import hashlib, json, os, sys, time
from pathlib import Path
db, key, path, label, version, source, md5x = sys.argv[1:8]
p = Path(path)
if not p.exists():
    sys.exit(0)
dp = Path(db) / "db_paths.json"
d = json.loads(dp.read_text()) if dp.exists() else {}
if key != "-":
    d[key] = str(p)
    dp.write_text(json.dumps(d, indent=2))
mp = Path(db) / "db_manifest.json"
m = json.loads(mp.read_text()) if mp.exists() else {}
st = p.stat()
e = m.get(str(p), {})
if not (e.get("bytes") == st.st_size and e.get("mtime") == int(st.st_mtime) and e.get("md5")):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    e.update(bytes=st.st_size, mtime=int(st.st_mtime), md5=h.hexdigest(),
             hashed_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
e.update(label=label, version=version, source=source)
if md5x != "-":
    e["md5_expected"] = md5x
m[str(p)] = e
mp.write_text(json.dumps(m, indent=2, ensure_ascii=False))
PYEOF
}

# ── 1. SILVA ────────────────────────────────────────────────────────────────
if want silva; then
  hdr "[1/5] SILVA 144 — 16S / 18S (~330 MB)"
  S="$DB_DIR/SILVA"
  BASE="https://www.arb-silva.de/fileadmin/silva_databases/release_144/DADA2/1.36.0/SSU"
  # The genus trainset feeds assignTaxonomy(); the species file is picked up
  # automatically from the same folder by version (v144) for addSpecies().
  fetch "$S/silva_nr99_v144_toGenus_trainset.fa.gz" 9e8a4c92c81efeeecdd091b1249f40fd \
        "$BASE/silva_nr99_v144_toGenus_trainset.fa.gz" &&
    register SILVA_16S "$S/silva_nr99_v144_toGenus_trainset.fa.gz" "SILVA 144" "144 (DADA2 1.36.0 build, 2026-08-12)" \
             "$BASE/silva_nr99_v144_toGenus_trainset.fa.gz" 9e8a4c92c81efeeecdd091b1249f40fd
  fetch "$S/silva_v144_assignSpecies.fa.gz" 444de7c0cce0b66addda7a3f8b38e012 \
        "$BASE/silva_v144_assignSpecies.fa.gz" &&
    register SILVA_16S_sp "$S/silva_v144_assignSpecies.fa.gz" "SILVA 144 species" "144" \
             "$BASE/silva_v144_assignSpecies.fa.gz" 444de7c0cce0b66addda7a3f8b38e012
  "$PY" - "$DB_DIR" <<'PYEOF'
import json, sys
from pathlib import Path
p = Path(sys.argv[1]) / "db_paths.json"
d = json.loads(p.read_text()) if p.exists() else {}
d["SILVA_144"] = d.get("SILVA_16S", "")
p.write_text(json.dumps(d, indent=2))
PYEOF
  if [[ $WANT_SILVA138 == 1 ]]; then
    Z="https://zenodo.org/records/14169026/files"
    fetch "$S/silva_nr99_v138.2_toGenus_trainset.fa.gz" 1764e2a36b4500ccb1c7d5261948a414 \
          "$Z/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1" &&
      register SILVA_138_2 "$S/silva_nr99_v138.2_toGenus_trainset.fa.gz" "SILVA 138.2" "138.2 (Nov 2024)" \
               "$Z/silva_nr99_v138.2_toGenus_trainset.fa.gz" 1764e2a36b4500ccb1c7d5261948a414
    fetch "$S/silva_v138.2_assignSpecies.fa.gz" 62fd939bd1aa9832d2089847f9f0b87d \
          "$Z/silva_v138.2_assignSpecies.fa.gz?download=1" &&
      register - "$S/silva_v138.2_assignSpecies.fa.gz" "SILVA 138.2 species" "138.2" \
               "$Z/silva_v138.2_assignSpecies.fa.gz" 62fd939bd1aa9832d2089847f9f0b87d
  fi
fi

# ── 2. UNITE ────────────────────────────────────────────────────────────────
if want unite; then
  hdr "[2/5] UNITE general FASTA release 10.0, 2025-02-19 — ITS fungi (~30 MB)"
  U="$DB_DIR/UNITE"
  UF="$U/sh_general_release_dynamic_19.02.2025.fasta"
  if [[ ! -s "$UF" ]]; then
    fetch "$U/sh_general_release_19.02.2025.tgz" 8ffa0bb111e9421e5f8ca8448ac5b7e4 \
          "https://s3.hpc.ut.ee/plutof-public/original/9489f7bc-7cc1-4e0a-84dc-c732476b9acd.tgz" &&
      tar -xzf "$U/sh_general_release_19.02.2025.tgz" -C "$U" sh_general_release_dynamic_19.02.2025.fasta
  fi
  if [[ -s "$UF" ]]; then
    ok "UNITE fungi: $(basename "$UF")"
    register UNITE_ITS1 "$UF" "UNITE general release 2025-02-19 (fungi)" "10.0" \
             "https://doi.org/10.15156/BIO/3301229" -
    register UNITE_ITS2 "$UF" "UNITE general release 2025-02-19 (fungi)" "10.0" \
             "https://doi.org/10.15156/BIO/3301229" -
  else
    fail "UNITE not installed — manual: https://doi.org/10.15156/BIO/3301229 → $U/"
  fi
  if [[ $WANT_UNITE_EUK == 1 ]]; then
    UE="$U/sh_general_release_dynamic_all_19.02.2025.fasta"
    if [[ ! -s "$UE" ]]; then
      fetch "$U/sh_general_release_all_19.02.2025.tgz" - \
            "https://s3.hpc.ut.ee/plutof-public/original/e861a3d6-54f4-42dc-882a-5f129beac39a.tgz" &&
        tar -xzf "$U/sh_general_release_all_19.02.2025.tgz" -C "$U" --wildcards '*dynamic_all*' 2>/dev/null
    fi
    [[ -s "$UE" ]] && register UNITE_ALL_EUK "$UE" "UNITE general release 2025-02-19 (all eukaryotes)" \
                               "10.0" "https://doi.org/10.15156/BIO/3301231" -
  fi
fi

# ── 3. PR2 / NemaBase ───────────────────────────────────────────────────────
if want pr2; then
  hdr "[3/5] PR2 5.1.1 — 18S eukaryotes (~55 MB)"
  N="$DB_DIR/NemaBase"
  PR2URL="https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz"
  fetch "$N/pr2_version_5.1.1_SSU_dada2.fasta.gz" - "$PR2URL" &&
    register PR2_18S "$N/pr2_version_5.1.1_SSU_dada2.fasta.gz" "PR2 5.1.1" "5.1.1 (2025-10-27)" "$PR2URL" -
  # 18S-NemaBase (Gattoni et al. 2023) is published on GitHub as a plain FASTA
  # with a 16-level WoRMS lineage after the accession; the old Zenodo link in
  # this script never existed. It is rewritten here into a DADA2 trainset with
  # the seven ranks the pipelines and report use (Kingdom..Species).
  NB="$N/18S_NemaBase_dada2.fasta.gz"
  if [[ ! -s "$NB" ]]; then
    NBURL="https://raw.githubusercontent.com/WormsEtAl/18SNemaBase/main/18S_NemaBase_fix.fasta.zip"
    fetch "$N/18S_NemaBase_fix.fasta.zip" - "$NBURL" && "$PY" - "$N/18S_NemaBase_fix.fasta.zip" "$NB" <<'PYEOF'
import gzip, io, sys, zipfile
src, dst = sys.argv[1], sys.argv[2]
z = zipfile.ZipFile(src)
name = next(n for n in z.namelist() if n.endswith((".fasta", ".fa")) and not n.startswith("__MACOSX"))
# accession_Eukaryota;Kingdom;Phylum;subphylum;superclass;Class;subclass;superorder;Order;
#   suborder;infraorder;superfamily;Family;subfamily;Genus;Genus_species
pick = [1, 2, 5, 8, 12, 14, 15]
n = 0
with z.open(name) as fh, gzip.open(dst, "wt") as out:
    keep = False
    for ln in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
        if ln.startswith(">"):
            lin = ln[1:].strip().split("_", 1)[-1].split(";")
            if len(lin) < 16:
                keep = False
                continue
            ranks = ["" if lin[i].startswith("NA_") or lin[i] == "NA" else lin[i] for i in pick]
            ranks[-1] = ranks[-1].replace("_", " ")
            out.write(">" + ";".join(ranks) + ";\n")
            keep, n = True, n + 1
        elif keep:
            out.write(ln)
print(f"  18S-NemaBase: {n} nematode sequences -> {dst}")
PYEOF
  fi
  if [[ -s "$NB" ]]; then
    register NemaBase_18S "$NB" "18S-NemaBase (Gattoni et al. 2023)" "GitHub WormsEtAl/18SNemaBase" \
             "https://github.com/WormsEtAl/18SNemaBase" -
  else
    warn "18S-NemaBase not available — PR2 is used for 18S nematodes"
  fi
  # install.sh used to save PR2 5.0 under the name NemaBase_18S.fa.gz and register
  # it as NemaBase. A NemaBase_18S entry whose headers have PR2's nine levels is
  # pointed at the real NemaBase, or at the current PR2, instead.
  "$PY" - "$DB_DIR" <<'PYEOF'
import gzip, json, sys
from pathlib import Path
p = Path(sys.argv[1]) / "db_paths.json"
d = json.loads(p.read_text()) if p.exists() else {}
cur = d.get("NemaBase_18S", "")
def depth(f):
    try:
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt") as fh:
            for ln in fh:
                if ln.startswith(">"):
                    return len([x for x in ln[1:].strip().split(";") if x])
    except Exception:
        return 0
    return 0
real = Path(sys.argv[1]) / "NemaBase" / "18S_NemaBase_dada2.fasta.gz"
if not cur or not Path(cur).exists() or depth(cur) >= 9:
    new = str(real) if real.exists() else d.get("PR2_18S", "")
    if new and new != cur:
        d["NemaBase_18S"] = new
        p.write_text(json.dumps(d, indent=2))
        print(f"  NemaBase_18S -> {Path(new).name}" + (f" (was {Path(cur).name}, a PR2 copy)" if cur else ""))
PYEOF
fi

# ── 4. MIDORI2 ──────────────────────────────────────────────────────────────
if want midori; then
  hdr "[4/5] MIDORI2 GenBank 273 (2026-09-10) — COI and 12S (~110 MB)"
  M="$DB_DIR/MIDORI2"; mkdir -p "$M"
  # MIDORI2 ships each release as one ZIP of every gene (~380 MB); only the two
  # files we use are read out of it, by HTTP range requests.
  ZIP="https://zenodo.org/api/records/22685701/files/MIDORI2_GB273_DADA2.zip/content"
  need=()
  [[ -s "$M/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz" ]] || \
    need+=(DADA2/uniq/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz "$M/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz")
  [[ -s "$M/MIDORI2_UNIQ_NUC_GB273_srRNA_DADA2.fasta.gz" ]] || \
    need+=(DADA2/uniq/MIDORI2_UNIQ_NUC_GB273_srRNA_DADA2.fasta.gz "$M/MIDORI2_UNIQ_NUC_GB273_srRNA_DADA2.fasta.gz")
  if [[ ${#need[@]} -gt 0 ]]; then
    "$PY" "$APP_DIR/backend/python_scripts/zip_member_fetch.py" "$ZIP" "${need[@]}" \
      || fail "MIDORI2 download failed — manual: https://doi.org/10.5281/zenodo.22685700 (MIDORI2_GB273_DADA2.zip)"
  fi
  [[ -s "$M/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz" ]] && ok "MIDORI2 COI" &&
    register MIDORI2_COX1 "$M/MIDORI2_UNIQ_NUC_GB273_CO1_DADA2.fasta.gz" "MIDORI2 GenBank 273 (COI)" \
             "GB273 (2026-09-10)" "https://doi.org/10.5281/zenodo.22685700" -
  [[ -s "$M/MIDORI2_UNIQ_NUC_GB273_srRNA_DADA2.fasta.gz" ]] && ok "MIDORI2 12S (srRNA)" &&
    register MIDORI2_12S "$M/MIDORI2_UNIQ_NUC_GB273_srRNA_DADA2.fasta.gz" "MIDORI2 GenBank 273 (12S srRNA)" \
             "GB273 (2026-09-10)" "https://doi.org/10.5281/zenodo.22685700" -
  # the "_sp" species-trainset name the COX1 pipeline used to look for is no
  # longer needed: MIDORI2 trainsets already end at species.
fi

# ── 5. Tax4Fun2 ─────────────────────────────────────────────────────────────
if want tax4fun; then
  hdr "[5/5] Tax4Fun2 reference data v2 — 16S functional prediction (~2 GB)"
  T="$DB_DIR/Tax4Fun2"
  if [[ -d "$T/Tax4Fun2_ReferenceData_v2" ]]; then
    ok "have Tax4Fun2 reference data"
  else
    # Tax4Fun2's own downloader points at CloudStor, which closed in 2023; the
    # authors' Zenodo deposit (10035668) carries the same reference data.
    T2URL="https://zenodo.org/api/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz/content"
    fetch "$T/Tax4Fun2_ReferenceData_v2.tar.gz" 5c2edff2f597e5f061f79a5782a27567 "$T2URL" &&
      tar -xzf "$T/Tax4Fun2_ReferenceData_v2.tar.gz" -C "$T" && rm -f "$T/Tax4Fun2_ReferenceData_v2.tar.gz" ||
      warn "Tax4Fun2 reference data failed — manual: https://doi.org/10.5281/zenodo.10035668"
  fi
  [[ -d "$T/Tax4Fun2_ReferenceData_v2" ]] && "$PY" - "$DB_DIR" "$T/Tax4Fun2_ReferenceData_v2" <<'PYEOF'
import json, sys
from pathlib import Path
p = Path(sys.argv[1]) / "db_paths.json"
d = json.loads(p.read_text()) if p.exists() else {}
d["Tax4Fun2_ref"] = sys.argv[2]
p.write_text(json.dumps(d, indent=2))
PYEOF
fi

# ── GTDB (on request) ───────────────────────────────────────────────────────
if [[ $WANT_GTDB == 1 ]]; then
  hdr "[+] GTDB R232 (2026-04-15) — 16S trainset built from GTDB's SSU files (~35 MB download)"
  G="$DB_DIR/GTDB"
  "$PY" "$APP_DIR/backend/python_scripts/build_gtdb_dada2.py" --release 232 --out "$G" &&
    register GTDB_16S "$G/gtdb_r232_ssu_toGenus_trainset.fa.gz" "GTDB R232" "R232 (2026-04-15)" \
             "https://data.gtdb.ecogenomic.org/releases/release232/" - &&
    register - "$G/gtdb_r232_ssu_assignSpecies.fa.gz" "GTDB R232 species" "R232" \
             "https://data.gtdb.ecogenomic.org/releases/release232/" - ||
    fail "GTDB build failed"
fi

# ── Summary ─────────────────────────────────────────────────────────────────
hdr "Summary — backend/databases/db_paths.json"
"$PY" - "$DB_DIR" <<'PYEOF'
import json, sys
from pathlib import Path
db = Path(sys.argv[1])
d = json.loads((db / "db_paths.json").read_text()) if (db / "db_paths.json").exists() else {}
m = json.loads((db / "db_manifest.json").read_text()) if (db / "db_manifest.json").exists() else {}
for k in ("SILVA_16S", "SILVA_16S_sp", "SILVA_138_2", "GTDB_16S", "UNITE_ITS1", "UNITE_ALL_EUK",
          "PR2_18S", "NemaBase_18S", "MIDORI2_COX1", "MIDORI2_12S", "Tax4Fun2_ref"):
    v = d.get(k)
    if not v:
        continue
    e = m.get(v, {})
    okk = Path(v).exists()
    chk = ("md5 ok" if e.get("md5_expected") == e.get("md5") else "MD5 MISMATCH") if e.get("md5_expected") else ""
    print(f"  {'✔' if okk else '✘'} {k:14s} {e.get('label') or Path(v).name:45s} {chk}")
PYEOF
echo
echo "  Old reference files are left where they are; pick between versions in the run settings."
echo "============================================================"
