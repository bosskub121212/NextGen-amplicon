#!/usr/bin/env bash
# =============================================================================
#  setup_shotgun.sh — tools and reference index for SHOTGUN METAGENOME runs (beta)
#
#    bash ~/r16s-app/setup_shotgun.sh                  # tools + Kraken2 Standard-16 + human host index
#    bash ~/r16s-app/setup_shotgun.sh --tools-only
#    bash ~/r16s-app/setup_shotgun.sh --index standard_16   # or standard_08, pluspf_16, pluspf_08, standard
#    bash ~/r16s-app/setup_shotgun.sh --no-host        # skip the human index for host removal
#    bash ~/r16s-app/setup_shotgun.sh --status
#
#  Everything goes into places of its own, nothing existing is touched:
#    conda env  ngamp-meta   python 3.11, kraken2, bracken, fastp, hostile >= 2, chopper, seqkit
#    index      backend/databases/kraken2/<name>/   (db_paths.json key: kraken2_db)
#    host index hostile's cache (~/.local/share/hostile), human T2T-CHM13 + HLA
#
#  Kraken2 indexes (Ben Langmead's builds, genome-idx.s3.amazonaws.com, newest is taken):
#    standard_16  RefSeq archaea, bacteria, viral, plasmid, human, UniVec   ~12 GB, needs ~16 GB RAM
#    pluspf_16    the same + protozoa and fungi                             ~12 GB, needs ~16 GB RAM
#    standard_08 / pluspf_08                                                 ~6 GB,  needs ~8 GB RAM
#    standard     full size (most sensitive)                                ~85 GB, needs ~90 GB RAM
#  The capped (16 / 8 GB) indexes keep a subset of k-mers: slightly less sensitive
#  for rare taxa, same species calls for abundant ones.
# =============================================================================
set -uo pipefail
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DBDIR="$APP/backend/databases/kraken2"
DBP="$APP/backend/databases/db_paths.json"
ENV=ngamp-meta
G='\033[0;32m'; Y='\033[1;33m'; R='\033[0;31m'; B='\033[1m'; N='\033[0m'
ok(){ echo -e "  ${G}✓${N} $*"; }; warn(){ echo -e "  ${Y}!${N} $*"; }; bad(){ echo -e "  ${R}✗${N} $*"; }
step(){ echo -e "\n${B}══ $* ══${N}"; }

TOOLS=1; INDEX=1; HOST=1; KIND="standard_16"; STATUS=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tools-only) INDEX=0; HOST=0 ;;
    --index) KIND="$2"; shift ;;
    --no-host) HOST=0 ;;
    --status) STATUS=1 ;;
    -h|--help) sed -n 2,25p "$0"; exit 0 ;;
    *) warn "unknown option $1" ;;
  esac
  shift
done

CONDA=""
for c in mamba conda; do command -v "$c" >/dev/null 2>&1 && { CONDA="$c"; break; }; done
[[ -z "$CONDA" && -x "$HOME/miniconda3/bin/conda" ]] && CONDA="$HOME/miniconda3/bin/conda"
# The conda base folder. `mamba info --base` (mamba 2) prints " base environment : /path"
# instead of the bare path, which broke every tool path below on nextapp — so
# take the last path-looking token of whichever command answers, and check it.
BASE=""
for c in conda mamba "$HOME/miniconda3/bin/conda"; do
  if command -v "$c" >/dev/null 2>&1 || [[ -x "$c" ]]; then
    b=$("$c" info --base 2>/dev/null | grep -o '/[^ ]*' | tail -1)
    [[ -n "$b" && -d "$b/envs" ]] && { BASE="$b"; break; }
  fi
done
if [[ -z "$BASE" ]]; then
  for b in "${CONDA_EXE:+$(dirname "$(dirname "$CONDA_EXE")")}" "$HOME/miniconda3" "$HOME/miniforge3" \
           "$HOME/anaconda3" "$HOME/mambaforge"; do
    [[ -n "$b" && -d "$b/envs" ]] && { BASE="$b"; break; }
  done
fi
[[ -n "$BASE" ]] || BASE="$HOME/miniconda3"
EBIN="$BASE/envs/$ENV/bin"

if [[ $STATUS == 1 ]]; then
  for t in kraken2 bracken fastp hostile chopper seqkit; do
    [[ -x "$EBIN/$t" ]] && ok "$t" || bad "$t missing"
  done
  python3 - "$DBP" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text()) if pathlib.Path(sys.argv[1]).exists() else {}
k = d.get("kraken2_db", "")
ok = k and (pathlib.Path(k) / "hash.k2d").exists()
print(("\033[0;32m  ✓\033[0m" if ok else "\033[1;33m  !\033[0m") + f"  kraken2_db  {k or '(not set)'}")
PY
  exit 0
fi

if [[ $TOOLS == 1 ]]; then
  step "Conda env $ENV"
  if [[ -x "$EBIN/kraken2" && -x "$EBIN/hostile" ]]; then
    ok "already installed"
  else
    # python pinned: an unpinned solve picked python 3.15 and a hostile release
    # that does not start on it
    $CONDA create -y -n $ENV -c conda-forge -c bioconda python=3.11 "hostile>=2" kraken2 bracken \
      fastp chopper seqkit || { bad "conda env failed"; exit 1; }
    ok "env $ENV created"
  fi
  "$EBIN/kraken2" --version | head -1
  "$EBIN/hostile" --version
fi

if [[ $INDEX == 1 ]]; then
  step "Kraken2 index ($KIND)"
  mkdir -p "$DBDIR"
  # newest archive of that kind, from the bucket listing (names changed from
  # k2_standard_16gb_YYYYMMDD to k2_standard_16_GB_YYYYMMDD in 2025)
  pat=$(python3 - "$KIND" <<'PY'
import sys
k = sys.argv[1]
m = {"standard_16": r"k2_standard_16(gb|_GB)_\d{8}", "standard_08": r"k2_standard_08(gb|_GB)_\d{8}",
     "pluspf_16": r"k2_pluspf_16(gb|_GB)_\d{8}", "pluspf_08": r"k2_pluspf_08(gb|_GB)_\d{8}",
     "standard": r"k2_standard_\d{8}"}
print(m.get(k, ""))
PY
)
  [[ -n "$pat" ]] || { bad "unknown index kind $KIND"; exit 1; }
  NAME=$(python3 - "$pat" <<'PY'
import re, sys, urllib.parse, urllib.request
pat = re.compile(r"kraken/(" + sys.argv[1] + r")\.tar\.gz")
names, token = [], ""
for _ in range(20):
    url = "https://genome-idx.s3.amazonaws.com/?list-type=2&prefix=kraken/k2_" + \
          (f"&continuation-token={urllib.parse.quote(token)}" if token else "")
    x = urllib.request.urlopen(url, timeout=60).read().decode()
    names += [m.group(1) for m in pat.finditer(x)]
    t = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", x)
    if not t:
        break
    token = t.group(1)
names = sorted(set(names), key=lambda n: n[-8:])
print(names[-1] if names else "")
PY
)
  [[ -n "$NAME" ]] || { bad "could not list genome-idx.s3.amazonaws.com — check the internet connection"; exit 1; }
  DEST="$DBDIR/$NAME"
  if [[ -s "$DEST/hash.k2d" ]]; then
    ok "$NAME already installed"
  else
    avail=$(df -BG --output=avail "$DBDIR" | tail -1 | tr -dc 0-9)
    echo "  downloading $NAME (newest) — ${avail} GB free"
    mkdir -p "$DEST"
    curl -fL --retry 5 -C - -o "$DEST.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/$NAME.tar.gz" \
      || { bad "download failed — re-run to resume"; exit 1; }
    tar -xzf "$DEST.tar.gz" -C "$DEST" && rm -f "$DEST.tar.gz" && ok "$NAME unpacked" \
      || { bad "unpack failed"; exit 1; }
  fi
  python3 - "$DBP" "$DEST" <<'PY'
import json, sys, pathlib
p = pathlib.Path(sys.argv[1]); d = json.loads(p.read_text()) if p.exists() else {}
d["kraken2_db"] = sys.argv[2]
p.write_text(json.dumps(d, indent=2))
PY
  ok "registered as kraken2_db"
  ls "$DEST" | grep -c kmer_distrib | xargs -I{} echo "  Bracken read lengths available: {}"
fi

if [[ $HOST == 1 ]]; then
  step "Human host index (hostile)"
  # hostile calls bowtie2 / minimap2 by name: the env's bin folder has to be on PATH
  export PATH="$EBIN:$PATH"
  "$EBIN/hostile" index fetch --name human-t2t-hla --bowtie2 && ok "Bowtie2 index (Illumina)" \
    || warn "Bowtie2 index not fetched — it downloads on the first run that needs it"
  "$EBIN/hostile" index fetch --name human-t2t-hla --minimap2 && ok "Minimap2 index (Nanopore)" \
    || warn "Minimap2 index not fetched — it downloads on the first run that needs it"
fi
echo ""
ok "Done. No backend restart needed for tools or indexes."
