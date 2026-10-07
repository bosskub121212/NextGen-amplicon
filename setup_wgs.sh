#!/usr/bin/env bash
# =============================================================================
#  NextGen-Amplicon — ONT-WGS (bacterial isolate genomes) tools + databases
#
#  Run once on the analysis machine (WSL):
#    bash ~/r16s-app/setup_wgs.sh                 # everything (≈ 15 GB download)
#    bash ~/r16s-app/setup_wgs.sh --tools-only    # conda environments only
#    bash ~/r16s-app/setup_wgs.sh --dbs-only      # databases only
#    bash ~/r16s-app/setup_wgs.sh --gtdb-small    # 0.4 GB GTDB sketch instead of 3.6 GB
#    bash ~/r16s-app/setup_wgs.sh --skip-bakta --skip-checkm2 --skip-mobsuite
#
#  Safe to re-run: environments and databases that are already in place are
#  left alone, so a run that stopped half-way just carries on.
#
#  Conda environments (kept apart because their dependencies do not mix):
#    ngamp-wgs flye seqkit barrnap sourmash skani mlst abricate
#              ncbi-amrfinderplus vsearch filtlong minimap2 samtools
#              (its own env, python 3.11: adding these to a hand-made env
#               pinned to python 3.13 left the solver running all night)
#    medaka    medaka                (polishing)
#    bakta     bakta                 (annotation)
#    mobsuite  mob_suite             (plasmid reconstruction)
#    checkm2   checkm2               (completeness / contamination)
#
#  Databases → ~/r16s-app/backend/databases/wgs/, registered in db_paths.json:
#    gtdb_sourmash   GTDB rs226 species representatives, sourmash k=31   3.6 GB
#    gtdb_lineages   GTDB rs226 lineages                                   4 MB
#    bakta_db        Bakta light database                              ≈ 1.5 GB
#    checkm2_db      CheckM2 DIAMOND database                          ≈ 3 GB
#    (AMRFinderPlus, MOB-suite and abricate/VFDB keep their own databases
#     inside their environments: amrfinder -u, mob_init)
# =============================================================================
set -uo pipefail

GREEN='\033[0;32m'; YEL='\033[1;33m'; CYAN='\033[0;36m'; RED='\033[0;31m'
BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YEL}  !${NC}  $1"; }
fail() { echo -e "${RED}  ✗${NC}  $1"; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

DO_TOOLS=1; DO_DBS=1; GTDB_SMALL=0; NO_GTDB=0
SKIP_BAKTA=0; SKIP_CHECKM2=0; SKIP_MOB=0
for a in "$@"; do
  case "$a" in
    --tools-only)    DO_DBS=0 ;;
    --dbs-only)      DO_TOOLS=0 ;;
    --gtdb-small)    GTDB_SMALL=1 ;;
    --no-gtdb)       NO_GTDB=1 ;;
    --skip-bakta)    SKIP_BAKTA=1 ;;
    --skip-checkm2)  SKIP_CHECKM2=1 ;;
    --skip-mobsuite) SKIP_MOB=1 ;;
    -h|--help) sed -n 2,30p "$0"; exit 0 ;;
    *) warn "unknown option $a" ;;
  esac
done

APP_DIR="$HOME/r16s-app"
DB_DIR="$APP_DIR/backend/databases"
WGS_DB="$DB_DIR/wgs"
DBP="$DB_DIR/db_paths.json"
mkdir -p "$WGS_DB"

# ── conda / mamba ─────────────────────────────────────────────────────────────
step "conda"
CONDA_BASE=""
for c in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge"; do
  [[ -x "$c/bin/conda" ]] && { CONDA_BASE="$c"; break; }
done
if [[ -z "$CONDA_BASE" ]] && command -v conda >/dev/null; then
  CONDA_BASE="$(conda info --base)"
fi
[[ -z "$CONDA_BASE" ]] && { fail "conda not found — install Miniconda3 first (install.sh)"; exit 1; }
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
SOLVER=(conda)
[[ -x "$CONDA_BASE/bin/mamba" ]] && SOLVER=("$CONDA_BASE/bin/mamba")
ok "conda base $CONDA_BASE (solver: ${SOLVER[0]##*/})"
CH=(-c conda-forge -c bioconda --override-channels --strict-channel-priority)
WENV="ngamp-wgs"

env_bin() { echo "$CONDA_BASE/envs/$1/bin/$2"; }
has()     { [[ -x "$(env_bin "$1" "$2")" ]]; }

ensure_env() {   # ensure_env <env> <probe-binary> <packages...>
  local env="$1" probe="$2"; shift 2
  if has "$env" "$probe"; then
    ok "$env: $probe present"
    return 0
  fi
  if [[ -d "$CONDA_BASE/envs/$env" ]]; then
    info "$env exists — adding: $*"
    "${SOLVER[@]}" install -y -n "$env" "${CH[@]}" "$@" || { fail "install into $env failed"; return 1; }
  else
    info "creating $env: $*"
    "${SOLVER[@]}" create -y -n "$env" "${CH[@]}" "$@" || { fail "create $env failed"; return 1; }
  fi
  has "$env" "$probe" && ok "$env ready" || { fail "$env: $probe still missing"; return 1; }
}

if [[ $DO_TOOLS == 1 ]]; then
  step "Conda environments"
  # One fresh environment for the light tools. Never installed INTO an
  # existing env: an env made by hand (e.g. "wgs" with flye/seqkit on python
  # 3.13) pins versions the Perl tools (mlst, abricate) cannot meet, and the
  # solver then searches for hours instead of failing.
  WGS_PKGS=(python=3.11 flye seqkit barrnap sourmash skani mlst abricate
            ncbi-amrfinderplus vsearch filtlong minimap2 samtools)
  PROBES=(flye seqkit barrnap sourmash skani mlst abricate amrfinder vsearch
          filtlong minimap2 samtools)
  MISSING=0
  for b in "${PROBES[@]}"; do has "$WENV" "$b" || MISSING=1; done
  if [[ $MISSING == 0 ]]; then
    ok "$WENV: all tools present"
  else
    [[ -d "$CONDA_BASE/envs/$WENV" ]] && { info "rebuilding incomplete $WENV"; \
      "${SOLVER[@]}" env remove -y -n "$WENV" >/dev/null 2>&1 || true; }
    info "creating $WENV (a few minutes)"
    "${SOLVER[@]}" create -y -n "$WENV" "${CH[@]}" "${WGS_PKGS[@]}" || fail "$WENV create failed"
    for b in "${PROBES[@]}"; do has "$WENV" "$b" || fail "$WENV: $b missing"; done
  fi
  ensure_env medaka medaka_consensus "medaka>=2.0"
  [[ $SKIP_BAKTA   == 0 ]] && ensure_env bakta    bakta     bakta
  [[ $SKIP_MOB     == 0 ]] && ensure_env mobsuite mob_recon mob_suite
  [[ $SKIP_CHECKM2 == 0 ]] && ensure_env checkm2  checkm2   "checkm2>=1.1"
fi

# ── databases ────────────────────────────────────────────────────────────────
register() {   # register <key> <path>
  python3 - "$DBP" "$1" "$2" <<'PY'
import json, sys, pathlib
p, k, v = pathlib.Path(sys.argv[1]), sys.argv[2], sys.argv[3]
d = json.loads(p.read_text()) if p.exists() else {}
d[k] = v
p.write_text(json.dumps(d, indent=2))
PY
  ok "db_paths.json: $1 → $2"
}

fetch() {      # fetch <url> <dest>  (resumable)
  local url="$1" dest="$2"
  if command -v wget >/dev/null; then wget -c -q --show-progress -O "$dest" "$url"
  else curl -L -C - --progress-bar -o "$dest" "$url"; fi
}

if [[ $DO_DBS == 1 ]]; then
  step "AMRFinderPlus database"
  if has "$WENV" amrfinder; then
    CONDA_PREFIX="$CONDA_BASE/envs/$WENV" PATH="$CONDA_BASE/envs/$WENV/bin:$PATH" \
      amrfinder -u >/tmp/amrfinder_u.log 2>&1 && ok "AMRFinderPlus database up to date" \
      || warn "amrfinder -u failed — see /tmp/amrfinder_u.log"
  else warn "amrfinder not installed — skipped"; fi

  if [[ $SKIP_MOB == 0 ]] && has mobsuite mob_init; then
    step "MOB-suite database"
    # mob_suite's own conda post-link step already downloads the database
    # (≈450 MB) while the env is created; running mob_init again re-downloads
    # it silently, or waits forever on a lock that download left behind.
    MOBDB=$(ls -d "$CONDA_BASE"/envs/mobsuite/lib/python3*/site-packages/mob_suite/databases 2>/dev/null | head -1)
    if [[ -n "$MOBDB" ]] && ls "$MOBDB"/ncbi_plasmid_full_seqs.fas* >/dev/null 2>&1; then
      rm -f "$MOBDB/.lock"
      ok "MOB-suite database present"
    else
      [[ -n "$MOBDB" ]] && rm -f "$MOBDB/.lock"
      info "downloading MOB-suite database (≈450 MB)"
      PATH="$CONDA_BASE/envs/mobsuite/bin:$PATH" mob_init 2>&1 | tee /tmp/mob_init.log | grep -E "%|ERROR|done|complete" \
        ; [[ ${PIPESTATUS[0]} == 0 ]] && ok "MOB-suite database ready" || warn "mob_init failed — see /tmp/mob_init.log"
    fi
  fi

  if [[ $NO_GTDB == 0 ]]; then
    step "GTDB species representatives (sourmash)"
    GBASE="https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226"
    if [[ $GTDB_SMALL == 1 ]]; then GSIG="gtdb-rs226-reps.k31-sc10k.sig.zip"
    else GSIG="gtdb-rs226-reps.k31.sig.zip"; fi
    if [[ -s "$WGS_DB/$GSIG" ]] && unzip -tq "$WGS_DB/$GSIG" >/dev/null 2>&1; then
      ok "$GSIG present"
    else
      info "downloading $GSIG"
      fetch "$GBASE/$GSIG" "$WGS_DB/$GSIG" || fail "GTDB sketch download failed"
    fi
    [[ -s "$WGS_DB/gtdb-rs226-reps.lineages.csv.gz" ]] || \
      fetch "$GBASE/gtdb-rs226-reps.lineages.csv.gz" "$WGS_DB/gtdb-rs226-reps.lineages.csv.gz"
    [[ -s "$WGS_DB/$GSIG" ]] && register gtdb_sourmash "$WGS_DB/$GSIG"
    [[ -s "$WGS_DB/gtdb-rs226-reps.lineages.csv.gz" ]] && \
      register gtdb_lineages "$WGS_DB/gtdb-rs226-reps.lineages.csv.gz"
  fi

  if [[ $SKIP_BAKTA == 0 ]] && has bakta bakta_db; then
    step "Bakta database (light)"
    BDB=$(find "$WGS_DB" -maxdepth 3 -name "version.json" -path "*bakta*" 2>/dev/null | head -1)
    if [[ -n "$BDB" ]]; then
      ok "Bakta database present"
    else
      PATH="$CONDA_BASE/envs/bakta/bin:$PATH" bakta_db download --output "$WGS_DB/bakta" \
        --type light || fail "bakta_db download failed"
      BDB=$(find "$WGS_DB" -maxdepth 3 -name "version.json" -path "*bakta*" | head -1)
    fi
    if [[ -n "$BDB" ]]; then
      register bakta_db "$(dirname "$BDB")"
      PATH="$CONDA_BASE/envs/bakta/bin:$PATH" amrfinder_update --force_update \
        --database "$(dirname "$BDB")/amrfinderplus-db" >/dev/null 2>&1 || true
    fi
  fi

  if [[ $SKIP_CHECKM2 == 0 ]] && has checkm2 checkm2; then
    step "CheckM2 database"
    CDB=$(find "$WGS_DB/checkm2" -name "*.dmnd" 2>/dev/null | head -1)
    if [[ -z "$CDB" ]]; then
      PATH="$CONDA_BASE/envs/checkm2/bin:$PATH" checkm2 database --download \
        --path "$WGS_DB/checkm2" || fail "CheckM2 database download failed"
      CDB=$(find "$WGS_DB/checkm2" -name "*.dmnd" 2>/dev/null | head -1)
    else ok "CheckM2 database present"; fi
    [[ -n "$CDB" ]] && register checkm2_db "$CDB"
  fi
fi

# ── report ───────────────────────────────────────────────────────────────────
step "Status"
for t in $WENV:flye $WENV:seqkit $WENV:sourmash $WENV:skani $WENV:barrnap $WENV:vsearch \
         $WENV:mlst $WENV:amrfinder $WENV:abricate $WENV:filtlong medaka:medaka_consensus \
         bakta:bakta mobsuite:mob_recon checkm2:checkm2; do
  e=${t%%:*}; b=${t#*:}
  has "$e" "$b" && ok "$b ($e)" || warn "$b missing ($e) — that step will be skipped"
done
python3 - "$DBP" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text()) if pathlib.Path(sys.argv[1]).exists() else {}
for k in ("gtdb_sourmash", "gtdb_lineages", "bakta_db", "checkm2_db"):
    v = d.get(k, "")
    mark = "\033[0;32m  ✓\033[0m" if v and pathlib.Path(v).exists() else "\033[1;33m  !\033[0m"
    print(f"{mark}  {k:<14} {v or '(not set)'}")
PY
echo ""
info "No backend restart needed — the pipeline finds tools and databases at run time."
