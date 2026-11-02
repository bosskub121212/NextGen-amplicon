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
#  Read-based taxonomy (separate, run once; databases are large → put them on a big drive):
#    bash ~/r16s-app/setup_wgs.sh --readtax --db-dir /mnt/e/ngamp-db/readtax
#      env ngamp-readtax: kraken2 bracken centrifuger gottcha2 sylph minimap2 samtools
#      Kraken2   reuses db_paths.json kraken2_db (setup_shotgun.sh), nothing downloaded
#      Centrifuger  cfr_hpv+gbsarscov2 (RefSeq prokaryotes + human + viruses, 2023-10)  42 GB,
#                   needs ≥ 60 GB RAM in WSL; with less RAM: Centrifuge p+h+v (2016)    8 GB
#                   (force with --centrifuge-old / --centrifuger)
#      GOTTCHA2  RefSeq-GTDB multi-domain+ "fast" bundle (fast-profile mode)            32 GB
#    Before each download the free space on the target is checked; a download that would
#    leave less than 20 GB free is refused.  --no-centrifuge / --no-gottcha skip one.
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
#    ngamp-phylo blast mummer mafft iqtree fastme minimap2 samtools
#              (type strains: ANIb / ANIm, 16S + genome trees, reference mapping)
#
#  Databases → ~/r16s-app/backend/databases/wgs/, registered in db_paths.json:
#    gtdb_sourmash   GTDB rs226 species representatives, sourmash k=31   3.6 GB
#    gtdb_lineages   GTDB rs226 lineages                                   4 MB
#    bakta_db        Bakta light database                              ≈ 1.5 GB
#    checkm2_db      CheckM2 DIAMOND database                          ≈ 3 GB
#    ncbi_16s_type   NCBI RefSeq 16S, type material only                   8 MB
#    (type-strain genomes are downloaded per genus at run time and cached in
#     ~/.ngamp_wgs_refs)
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
READTAX=0; RT_DIR=""; CF_MODE=auto; NO_CF=0; NO_G2=0; KEEP_FREE_GB=20
ARGS=("$@")
for ((i = 0; i < ${#ARGS[@]}; i++)); do
  a="${ARGS[$i]}"
  case "$a" in
    --readtax)       READTAX=1; DO_TOOLS=0; DO_DBS=0 ;;
    --db-dir)        i=$((i + 1)); RT_DIR="${ARGS[$i]:-}" ;;
    --db-dir=*)      RT_DIR="${a#--db-dir=}" ;;
    --centrifuge-old) CF_MODE=old ;;
    --centrifuger)   CF_MODE=new ;;
    --no-centrifuge) NO_CF=1 ;;
    --no-gottcha)    NO_G2=1 ;;
    --tools-only)    DO_DBS=0 ;;
    --dbs-only)      DO_TOOLS=0 ;;
    --gtdb-small)    GTDB_SMALL=1 ;;
    --no-gtdb)       NO_GTDB=1 ;;
    --skip-bakta)    SKIP_BAKTA=1 ;;
    --skip-checkm2)  SKIP_CHECKM2=1 ;;
    --skip-mobsuite) SKIP_MOB=1 ;;
    -h|--help) sed -n 2,42p "$0"; exit 0 ;;
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
  # type strains / trees / mapping — own env, created fresh (never added to ngamp-wgs)
  PHY_PROBES=(blastn makeblastdb nucmer delta-filter mafft iqtree fastme minimap2 samtools)
  PMISS=0
  for b in "${PHY_PROBES[@]}"; do has ngamp-phylo "$b" || PMISS=1; done
  if [[ $PMISS == 0 ]]; then
    ok "ngamp-phylo: all tools present"
  else
    [[ -d "$CONDA_BASE/envs/ngamp-phylo" ]] && { info "rebuilding incomplete ngamp-phylo"; \
      "${SOLVER[@]}" env remove -y -n ngamp-phylo >/dev/null 2>&1 || true; }
    info "creating ngamp-phylo (a few minutes)"
    "${SOLVER[@]}" create -y -n ngamp-phylo "${CH[@]}" python=3.12 blast mummer mafft \
      iqtree fastme minimap2 samtools || fail "ngamp-phylo create failed"
    for b in "${PHY_PROBES[@]}"; do has ngamp-phylo "$b" || fail "ngamp-phylo: $b missing"; done
  fi
  if [[ $SKIP_BAKTA == 0 ]]; then
    # Unpinned, the solver takes the newest python (3.14) and then falls back to
    # bakta 1.5 — a 2022 release that cannot read the current database and has no
    # light database. Anything older than 1.9 is rebuilt.
    # Judged by what we need rather than by a version string (an old bakta on
    # python 3.14 may not even print its version): the light database needs
    # `bakta_db download --type`, added in bakta 1.8.
    if has bakta bakta && ! "$(env_bin bakta bakta_db)" download --help 2>&1 | grep -q -- "--type"; then
      warn "bakta in the 'bakta' env is too old (no light database) — rebuilding it"
      "${SOLVER[@]}" env remove -y -n bakta >/dev/null 2>&1 || true
      rm -rf "$CONDA_BASE/envs/bakta"
    fi
    ensure_env bakta bakta "python=3.11" "bakta>=1.9"
  fi
  [[ $SKIP_MOB     == 0 ]] && ensure_env mobsuite mob_recon "python=3.11" mob_suite
  if [[ $SKIP_CHECKM2 == 0 ]]; then
    ensure_env checkm2 checkm2 "checkm2>=1.1"
    # TensorFlow's shared libraries ask for an executable stack; newer glibc
    # (Ubuntu 25.04 / glibc 2.41+) refuses to load them:
    #   "libtensorflow_cc.so.2: cannot enable executable stack ... Invalid argument"
    # Clearing that flag on the TensorFlow libraries is the standard fix.
    if has checkm2 checkm2 && ! "$(env_bin checkm2 checkm2)" --version >/dev/null 2>&1; then
      if "$(env_bin checkm2 checkm2)" --version 2>&1 | grep -q "executable stack"; then
        info "checkm2: clearing the executable-stack flag on TensorFlow libraries"
        has checkm2 patchelf || "${SOLVER[@]}" install -y -n checkm2 -c conda-forge \
          --override-channels "patchelf>=0.18" >/dev/null
        PE="$(env_bin checkm2 patchelf)"
        find "$CONDA_BASE/envs/checkm2" -type f \( -name "libtensorflow*.so*" -o \
             \( -path "*site-packages/tensorflow/*" -name "*.so*" \) \) 2>/dev/null |
        while read -r so; do
          "$PE" --print-execstack "$so" 2>/dev/null | grep -q "X" && "$PE" --clear-execstack "$so"
        done
      fi
      "$(env_bin checkm2 checkm2)" --version >/dev/null 2>&1 && ok "checkm2 runs" \
        || fail "checkm2 still does not start — run: $(env_bin checkm2 checkm2) --version"
    fi
  fi
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

# ── read-based taxonomy (--readtax) ──────────────────────────────────────────
free_gb() {    # free GB on the filesystem holding $1
  df -P -BG "$1" 2>/dev/null | awk 'NR==2 {gsub("G","",$4); print $4}'
}
room_for() {   # room_for <dir> <GB needed>  → refuses when < KEEP_FREE_GB would remain
  local f; f=$(free_gb "$1")
  [[ -z "$f" ]] && { fail "cannot read free space of $1"; return 1; }
  if (( f - $2 < KEEP_FREE_GB )); then
    fail "not enough space on $(df -P "$1" | awk 'NR==2{print $6}'): ${f} GB free, this needs $2 GB and ${KEEP_FREE_GB} GB must stay free — skipped"
    return 1
  fi
  info "space: ${f} GB free on $(df -P "$1" | awk 'NR==2{print $6}'), needs $2 GB"
}
md5_ok() { [[ "$(md5sum "$1" 2>/dev/null | cut -d' ' -f1)" == "$2" ]]; }

if [[ $READTAX == 1 ]]; then
  RT_DIR="${RT_DIR:-$DB_DIR/readtax}"
  mkdir -p "$RT_DIR" || { fail "cannot create $RT_DIR"; exit 1; }
  RT_DIR="$(cd "$RT_DIR" && pwd)"
  if [[ "$RT_DIR" == /mnt/c/* ]] || [[ "$RT_DIR" == "$HOME"* ]]; then
    warn "$RT_DIR is on C: (WSL's own disk also lives on C:) — a bigger drive is better: --db-dir /mnt/e/ngamp-db/readtax"
  fi
  step "Read-based taxonomy tools (ngamp-readtax)"
  RENV=ngamp-readtax
  ensure_env "$RENV" gottcha2 python=3.11 kraken2 bracken centrifuger "gottcha2>=2.5" sylph minimap2 samtools \
    || warn "ngamp-readtax incomplete — the missing tools are skipped in the report"
  for b in kraken2 bracken centrifuger centrifuger-quant gottcha2 sylph; do
    has "$RENV" "$b" && ok "$b" || warn "$b missing in $RENV"
  done

  step "Kraken2 database (reused)"
  K2=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1])).get('kraken2_db',''))" "$DBP" 2>/dev/null)
  if [[ -n "$K2" && -s "$K2/hash.k2d" ]]; then
    ok "kraken2_db $K2 ($(ls "$K2"/database*mers.kmer_distrib 2>/dev/null | wc -l) Bracken read lengths)"
  else warn "kraken2_db not set — run setup_shotgun.sh first (the Kraken2 column is then skipped)"; fi

  if [[ $NO_CF == 0 ]]; then
    step "Centrifuger / Centrifuge index"
    MEM=$(awk '/MemTotal/ {printf "%d", $2/1048576}' /proc/meminfo)
    MODE=$CF_MODE
    if [[ $MODE == auto ]]; then
      if (( MEM >= 60 )); then MODE=new; else MODE=old; fi
      info "WSL memory ${MEM} GB → $([[ $MODE == new ]] && echo "Centrifuger RefSeq 2023 index (needs ≈ 45 GB RAM)" || echo "Centrifuge p+h+v (needs ≈ 8 GB RAM; the 2023 index needs ≥ 60 GB)")"
    fi
    if [[ $MODE == new ]]; then
      CD="$RT_DIR/centrifuger"; mkdir -p "$CD"; PFX="$CD/cfr_hpv+gbsarscov2"
      ZB="https://zenodo.org/records/10023239/files"
      declare -A CMD5=([1]=fdcd6fdb45a2cb79be25b19c35da06c2 [2]=634d66fa33e433cdab9ed3f0d3b43193 [3]=a49a5a91b883aa98461c9cf1ed5a9588)
      if [[ -s "$PFX.1.cfr" && -f "$CD/.verified" ]]; then ok "Centrifuger index present"
      else
        NEED=46; [[ -s "$PFX.1.cfr" ]] && NEED=$(( 46 - $(du -BG "$PFX.1.cfr" | cut -f1 | tr -d G) ))
        if room_for "$CD" "$NEED"; then
          GOOD=1
          for k in 1 2 3; do
            if [[ -s "$PFX.$k.cfr" ]] && md5_ok "$PFX.$k.cfr" "${CMD5[$k]}"; then ok "cfr_hpv+gbsarscov2.$k.cfr already complete"; continue; fi
            info "downloading cfr_hpv+gbsarscov2.$k.cfr"
            fetch "$ZB/cfr_hpv+gbsarscov2.$k.cfr?download=1" "$PFX.$k.cfr"
            info "checking md5 (a few minutes for the 42 GB file)"
            if ! md5_ok "$PFX.$k.cfr" "${CMD5[$k]}"; then
              fail "md5 mismatch: cfr_hpv+gbsarscov2.$k.cfr — removed; run the command again"
              rm -f "$PFX.$k.cfr"; GOOD=0
            fi
          done
          if [[ $GOOD == 1 ]]; then
            touch "$CD/.verified"; ok "Centrifuger index ready"
            echo "cfr_hpv+gbsarscov2 (RefSeq prokaryotes, human, viruses; 2023-10)" > "$CD/db_name.txt"
          fi
        fi
      fi
      [[ -f "$CD/.verified" ]] && register readtax_centrifuge "$PFX"
    else
      has "$RENV" centrifuge || ensure_env "$RENV" centrifuge centrifuge || true
      CD="$RT_DIR/centrifuge"; mkdir -p "$CD"
      if ls "$CD"/p+h+v.1.cf >/dev/null 2>&1; then ok "Centrifuge p+h+v present"
      elif room_for "$CD" 17; then
        info "downloading Centrifuge p+h+v (8 GB)"
        fetch "https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz" "$CD/p+h+v.tar.gz"
        if tar -xzf "$CD/p+h+v.tar.gz" -C "$CD"; then
          rm -f "$CD/p+h+v.tar.gz"; ok "Centrifuge p+h+v ready"
          echo "p+h+v (RefSeq prokaryotes, human, viruses; 2016)" > "$CD/db_name.txt"
        else fail "Centrifuge index download/extract failed — run again to resume"; fi
      fi
      P1=$(ls "$CD"/*.1.cf 2>/dev/null | head -1)
      [[ -n "$P1" ]] && register readtax_centrifuge "${P1%.1.cf}"
    fi
  fi

  if [[ $NO_G2 == 0 ]]; then
    step "GOTTCHA2 database (RefSeq-GTDB multi-domain+, fast bundle)"
    GD="$RT_DIR/gottcha2"; mkdir -p "$GD"
    GTAR="$GD/gottcha_db_fast.tar"
    GP=$(ls "$GD"/*/*.fna.tax.tsv "$GD"/*.fna.tax.tsv 2>/dev/null | head -1)
    if [[ -n "$GP" && -f "$GD/.extracted" ]]; then ok "GOTTCHA2 database present"
    else
      NEED=66; [[ -s "$GTAR" ]] && NEED=$(( 66 - $(du -BG "$GTAR" | cut -f1 | tr -d G) ))
      if room_for "$GD" "$NEED"; then
        info "downloading gottcha_db_fast.tar (32 GB, resumable)"
        fetch "https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain%2B/gottcha_db_fast.tar" "$GTAR"
        if tar -xf "$GTAR" -C "$GD"; then
          rm -f "$GTAR"; touch "$GD/.extracted"; ok "GOTTCHA2 database ready"
          LM=$(curl -sIL "https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain%2B/gottcha_db_fast.tar" \
               | tr -d '\r' | awk -F': ' 'tolower($1)=="last-modified"{print $2}' | tail -1)
          LMD=$(date -d "$LM" +%Y-%m 2>/dev/null || echo "downloaded $(date +%Y-%m)")
          echo "RefSeq-GTDB multi-domain+ species signatures (fast bundle, $LMD)" > "$GD/db_name.txt"
        else fail "GOTTCHA2 download/extract failed — run again to resume"; fi
      fi
      GP=$(find "$GD" -maxdepth 3 -name "*.species.fna.tax.tsv" 2>/dev/null | head -1)
      [[ -z "$GP" ]] && GP=$(find "$GD" -maxdepth 3 -name "*.fna.tax.tsv" 2>/dev/null | head -1)
    fi
    [[ -n "$GP" && -f "$GD/.extracted" ]] && register readtax_gottcha2 "${GP%.tax.tsv}"
  fi

  step "Read-based taxonomy status"
  python3 - "$DBP" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text()) if pathlib.Path(sys.argv[1]).exists() else {}
for k, probe in (("kraken2_db", "/hash.k2d"), ("readtax_centrifuge", ""), ("readtax_gottcha2", ".tax.tsv")):
    v = d.get(k, "")
    ok = v and (pathlib.Path(v + probe).exists() if probe else any(pathlib.Path(v).parent.glob(pathlib.Path(v).name + ".1.cf*")))
    mark = "\033[0;32m  ✓\033[0m" if ok else "\033[1;33m  !\033[0m"
    print(f"{mark}  {k:<19} {v or '(not set)'}")
PY
  df -h "$RT_DIR" | tail -1 | awk '{print "  free on " $6 ": " $4}'
  info "No backend restart needed. New ONT-WGS jobs run it; finished jobs: Genome report → Add read-based taxonomy."
  exit 0
fi

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

if [[ $DO_DBS == 1 ]]; then
  step "NCBI 16S type-strain database"
  N16="$WGS_DB/ncbi_16s"; mkdir -p "$N16"
  if [[ -s "$N16/bacteria.16SrRNA.fna" ]]; then
    ok "NCBI RefSeq 16S present ($(grep -c '>' "$N16/bacteria.16SrRNA.fna") type-strain records)"
  else
    fetch "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz" \
      "$N16/bacteria.16SrRNA.fna.gz" && gunzip -f "$N16/bacteria.16SrRNA.fna.gz" \
      && ok "NCBI RefSeq 16S ready" || fail "NCBI 16S download failed"
  fi
  [[ -s "$N16/bacteria.16SrRNA.fna" ]] && register ncbi_16s_type "$N16/bacteria.16SrRNA.fna"
fi

# ── report ───────────────────────────────────────────────────────────────────
step "Status"
for t in $WENV:flye $WENV:seqkit $WENV:sourmash $WENV:skani $WENV:barrnap $WENV:vsearch \
         $WENV:mlst $WENV:amrfinder $WENV:abricate $WENV:filtlong medaka:medaka_consensus \
         bakta:bakta mobsuite:mob_recon checkm2:checkm2 ngamp-phylo:blastn \
         ngamp-phylo:nucmer ngamp-phylo:mafft ngamp-phylo:iqtree ngamp-phylo:fastme \
         ngamp-phylo:minimap2 ngamp-phylo:samtools; do
  e=${t%%:*}; b=${t#*:}
  if ! has "$e" "$b"; then warn "$b missing ($e) — that step will be skipped"
  elif [[ "$b" == checkm2 ]] && ! "$(env_bin "$e" "$b")" --version >/dev/null 2>&1; then
    warn "$b installed but does not start ($e) — that step will be skipped"
  else ok "$b ($e)"; fi
done
python3 - "$DBP" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text()) if pathlib.Path(sys.argv[1]).exists() else {}
for k in ("gtdb_sourmash", "gtdb_lineages", "bakta_db", "checkm2_db", "ncbi_16s_type"):
    v = d.get(k, "")
    mark = "\033[0;32m  ✓\033[0m" if v and pathlib.Path(v).exists() else "\033[1;33m  !\033[0m"
    print(f"{mark}  {k:<14} {v or '(not set)'}")
PY
echo ""
info "No backend restart needed — the pipeline finds tools and databases at run time."
