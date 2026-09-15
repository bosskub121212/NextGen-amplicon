#!/usr/bin/env bash
# =============================================================================
#  deploy_dev.sh — Sync the Windows workspace into the WSL runtime and restart
#
#  Run from WSL:
#    bash /mnt/c/Claude/r16s-app/deploy_dev.sh
#    bash /mnt/c/Claude/r16s-app/deploy_dev.sh --no-build    # skip npm build
#    bash /mnt/c/Claude/r16s-app/deploy_dev.sh --no-restart  # leave backend alone
#
#  This does NOT touch git. Pushing is a separate step:
#    bash ~/r16s-app/push_github.sh "your commit message"
#
#  WHY DIRECTORY SYNC INSTEAD OF A FILE LIST
#  -----------------------------------------
#  The old copy_all_extensions.sh named every file explicitly, so any NEW file
#  was silently skipped until someone remembered to add it — which bit us three
#  times (reorient_reads.py, DataPrepPanel.tsx, build_tree.R). The failure mode
#  is nasty: the app runs, but with a missing piece, and the only symptom is a
#  "[skip]" line buried in a log. Syncing whole directories with an exclude list
#  means new files are picked up automatically and only runtime data is protected.
# =============================================================================
set -euo pipefail

GREEN='\033[0;32m'; CYAN='\033[0;36m'; YEL='\033[1;33m'; RED='\033[0;31m'
BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YEL}  !${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

WIN_SRC="/mnt/c/Claude/r16s-app"
APP_DIR="$HOME/r16s-app"
DO_BUILD=1
DO_RESTART=1

for arg in "$@"; do
  case "$arg" in
    --no-build)   DO_BUILD=0 ;;
    --no-restart) DO_RESTART=0 ;;
    -h|--help)    sed -n '2,16p' "$0"; exit 0 ;;
    *) die "Unknown option: $arg" ;;
  esac
done

[[ -d "$WIN_SRC" ]] || die "Source not found: $WIN_SRC (is drive C: mounted in WSL?)"
[[ -d "$APP_DIR" ]] || die "Runtime not found: $APP_DIR (run install.sh first)"

echo ""
echo "============================================================"
echo "  NextGen-Amplicon — Deploy to WSL runtime"
echo "  Source : $WIN_SRC"
echo "  Dest   : $APP_DIR"
echo "============================================================"

# ── 1. Sync ───────────────────────────────────────────────────────────────
step "Syncing files"

# Never sync these: they are runtime state or machine-local, and clobbering
# them would destroy uploaded data, finished results, or the installed venv.
EXCL=(
  --exclude '.git/'            --exclude 'node_modules/'
  --exclude 'dist/'            --exclude 'build/'
  --exclude 'venv/'            --exclude '__pycache__/'
  --exclude '*.pyc'            --exclude '*.pyo'
  --exclude 'uploads/'         --exclude 'results/'
  --exclude 'databases/'       --exclude 'jobs_history.json'
  --exclude '.license_cache.json'
  --exclude 'logs/'            --exclude '*.log'
)

SYNC_LOG="$(mktemp)"
trap 'rm -f "$SYNC_LOG"' EXIT

if command -v rsync >/dev/null 2>&1; then
  rsync -rlt --itemize-changes "${EXCL[@]}" "$WIN_SRC/backend/"  "$APP_DIR/backend/"  >>"$SYNC_LOG"
  rsync -rlt --itemize-changes "${EXCL[@]}" "$WIN_SRC/frontend/" "$APP_DIR/frontend/" >>"$SYNC_LOG"
  rsync -rlt --itemize-changes "$WIN_SRC/version.json" "$APP_DIR/" >>"$SYNC_LOG"
  for f in "$WIN_SRC"/*.sh; do
    [[ -f "$f" ]] && rsync -rlt --itemize-changes "$f" "$APP_DIR/" >>"$SYNC_LOG"
  done
  CHANGED=$(grep -c '^[<>ch]' "$SYNC_LOG" 2>/dev/null || true)
  CHANGED=${CHANGED:-0}
  if [[ "$CHANGED" -gt 0 ]]; then
    ok "$CHANGED file(s) updated:"
    grep '^[<>ch]' "$SYNC_LOG" | awk '{print "        " $2}' | head -40
    [[ "$CHANGED" -gt 40 ]] && info "... and $((CHANGED - 40)) more"
  else
    info "Everything already up to date"
  fi
else
  warn "rsync not installed — falling back to cp (cannot report what changed)"
  warn "Install for better output:  sudo apt install rsync"
  mkdir -p "$APP_DIR/backend" "$APP_DIR/frontend"
  ( cd "$WIN_SRC" && find backend frontend -type f \
      ! -path '*/node_modules/*' ! -path '*/dist/*' ! -path '*/venv/*' \
      ! -path '*/__pycache__/*'  ! -path '*/uploads/*' ! -path '*/results/*' \
      ! -path '*/databases/*'    ! -name '*.pyc' ! -name '*.log' \
      -print0 | while IFS= read -r -d '' rel; do
        mkdir -p "$APP_DIR/$(dirname "$rel")"
        cp "$WIN_SRC/$rel" "$APP_DIR/$rel"
      done )
  cp "$WIN_SRC/version.json" "$APP_DIR/" 2>/dev/null || true
  cp "$WIN_SRC"/*.sh "$APP_DIR/" 2>/dev/null || true
  CHANGED=1
  ok "Files copied"
fi

APP_VER=$(grep -oP '"version"\s*:\s*"\K[^"]+' "$APP_DIR/version.json" 2>/dev/null || echo "?")
info "Runtime version: $APP_VER"

# ── 2. Sanity-check the pieces that fail silently ─────────────────────────
# Both of these are invoked as external processes by the R pipeline, so if
# they are missing or broken the run does not fail — it just quietly skips a
# feature. Check them here, loudly, instead.
step "Checking helper scripts"
if [[ -f "$APP_DIR/backend/python_scripts/reorient_reads.py" ]]; then
  if python3 "$APP_DIR/backend/python_scripts/reorient_reads.py" --help >/dev/null 2>&1; then
    ok "reorient_reads.py runs"
  else
    warn "reorient_reads.py present but NOT runnable — orientation repair will be skipped"
  fi
else
  warn "reorient_reads.py MISSING — orientation repair will be skipped"
fi

if [[ -f "$APP_DIR/backend/r_scripts/build_tree.R" ]]; then
  ok "build_tree.R present"
else
  warn "build_tree.R MISSING — no phylogenetic tree will be produced"
fi

# ── 3. R packages (fast check, install only what's missing) ───────────────
step "Checking R packages"
Rscript -e '
need_cran <- c("phangorn")
need_bioc <- c("decontam","DECIPHER","Biostrings")
miss_cran <- need_cran[!vapply(need_cran, requireNamespace, logical(1), quietly=TRUE)]
miss_bioc <- need_bioc[!vapply(need_bioc, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss_cran) == 0 && length(miss_bioc) == 0) {
  cat("  all present\n")
} else {
  if (length(miss_cran)) {
    cat("  installing (CRAN):", paste(miss_cran, collapse=", "), "\n")
    install.packages(miss_cran, repos="https://cloud.r-project.org", quiet=TRUE)
  }
  if (length(miss_bioc)) {
    cat("  installing (Bioconductor):", paste(miss_bioc, collapse=", "), "\n")
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(miss_bioc, ask=FALSE, update=FALSE)
  }
}' 2>&1 | grep -v "^Warning" || true
ok "R packages checked"

# ── 4. Frontend build (only when frontend sources actually changed) ───────
if [[ "$DO_BUILD" -eq 1 ]]; then
  FRONTEND_CHANGED=1
  if command -v rsync >/dev/null 2>&1; then
    grep '^[<>ch]' "$SYNC_LOG" 2>/dev/null | grep -q 'frontend' || FRONTEND_CHANGED=0
  fi
  if [[ "$FRONTEND_CHANGED" -eq 1 ]]; then
    step "Building frontend"
    cd "$APP_DIR/frontend"
    source "$HOME/.nvm/nvm.sh" 2>/dev/null || true
    nvm use 20 >/dev/null 2>&1 || true
    npm run build
    ok "Frontend built"
    cd "$APP_DIR"
  else
    step "Building frontend"
    info "No frontend changes — skipped (force with a frontend edit, or npm run build)"
  fi
else
  info "Frontend build skipped (--no-build)"
fi

# ── 5. Restart backend ────────────────────────────────────────────────────
if [[ "$DO_RESTART" -eq 1 ]]; then
  step "Restarting backend"
  # Two launchers produce different command lines — start_backend.sh gives
  # "uvicorn main:app", start.sh gives "uvicorn backend.main:app". Matching
  # only one leaves the other holding port 8000, and the new process then dies
  # with "address already in use" while the OLD code keeps serving.
  pkill -f "uvicorn.*main:app" 2>/dev/null && info "Old backend stopped" || info "No old backend running"
  sleep 2
  if command -v ss >/dev/null 2>&1 && ss -ltn 2>/dev/null | grep -q ':8000 '; then
    warn "Port 8000 is STILL in use — backend not restarted:"
    ss -ltnp 2>/dev/null | grep ':8000 ' || true
    warn "Free it, then:  bash $APP_DIR/start_backend.sh"
  else
    bash "$APP_DIR/start_backend.sh" >/dev/null 2>&1 &
    sleep 4
    if curl -fsS http://127.0.0.1:8000/ >/dev/null 2>&1; then
      ok "Backend responding on :8000"
    else
      warn "Backend launched but not responding yet — check:  bash $APP_DIR/start_backend.sh"
    fi
  fi
else
  info "Backend restart skipped (--no-restart)"
fi

echo ""
echo -e "${BOLD}${GREEN}  ✅  Deployed  (v$APP_VER)${NC}"
echo ""
echo "  App:  http://localhost:8000"
echo ""
echo "  To publish these changes to GitHub (separate step):"
echo "    bash $APP_DIR/push_github.sh \"your commit message\""
echo ""
