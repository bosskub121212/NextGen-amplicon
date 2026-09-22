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
FE_LOG="$(mktemp)"
trap 'rm -f "$SYNC_LOG" "$FE_LOG"' EXIT

if command -v rsync >/dev/null 2>&1; then
  # The frontend gets its OWN log. rsync itemizes paths relative to the transfer
  # root, so syncing "$WIN_SRC/frontend/" prints "src/App.tsx" — the word
  # "frontend" never appears in the line. The old check grepped the shared log
  # for "frontend" to decide whether to rebuild, so it matched nothing and the
  # frontend build was skipped on EVERY run: the backend updated, the browser
  # kept serving a months-old bundle, and the mismatch only showed up when new
  # backend output hit old frontend code.
  rsync -rlt --itemize-changes "${EXCL[@]}" "$WIN_SRC/backend/"  "$APP_DIR/backend/"  >>"$SYNC_LOG"
  rsync -rlt --itemize-changes "${EXCL[@]}" "$WIN_SRC/frontend/" "$APP_DIR/frontend/" >>"$FE_LOG"
  cat "$FE_LOG" >>"$SYNC_LOG"
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

# Each of these is source()d or spawned by the pipeline behind an exists() or
# file.exists() guard, so a missing one degrades quietly rather than erroring —
# exactly the failure mode that is hardest to notice. Name what is lost.
declare -A R_HELPERS=(
  [build_tree.R]="no phylogenetic tree will be produced"
  [tax_helpers.R]="SILVA 144 rank depth and '--other' genus names will be mishandled"
  [plot_helpers.R]="taxonomy bars will not match between the PDFs and the browser"
  [qc_helpers.R]="truncLen will not be checked against the reads before filtering"
  [add_species.R]="species assignment runs inline — much slower on a low-RAM machine"
)
for h in "${!R_HELPERS[@]}"; do
  if [[ -f "$APP_DIR/backend/r_scripts/$h" ]]; then
    ok "$h present"
  else
    warn "$h MISSING — ${R_HELPERS[$h]}"
  fi
done

# A helper with a syntax error is worse than a missing one: source() aborts the
# step that needed it, and the pipeline reports something unrelated.
if Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f))' \
     "$APP_DIR"/backend/r_scripts/{tax,plot,qc}_helpers.R \
     "$APP_DIR"/backend/r_scripts/{add_species,build_tree}.R >/dev/null 2>&1; then
  ok "R helper scripts parse cleanly"
else
  warn "An R helper script has a SYNTAX ERROR — run: Rscript -e 'parse(\"<file>\")'"
fi

# report_builder.py falls back to emitting HTML when weasyprint is absent, so a
# missing install shows up as "no PDF appeared" rather than as an error.
#
# Check it with the interpreter THE BACKEND ACTUALLY RUNS, not whatever python3
# resolves to in this shell. start.sh launches $APP_DIR/venv/bin/uvicorn, so a
# `pip install` typed at a conda prompt lands somewhere the backend cannot see —
# and a bare `python3 -c "import weasyprint"` here would happily report success
# while the app keeps answering "weasyprint is not installed".
VENV_PY="$APP_DIR/venv/bin/python3"
[[ -x "$VENV_PY" ]] || VENV_PY="$APP_DIR/venv/bin/python"
if [[ -x "$VENV_PY" ]]; then
  if "$VENV_PY" -c "import weasyprint" >/dev/null 2>&1; then
    ok "weasyprint present in the backend venv — PDF reports available"
  else
    warn "weasyprint missing from the backend venv — installing it there"
    if "$VENV_PY" -m pip install --quiet weasyprint >/dev/null 2>&1 \
       && "$VENV_PY" -c "import weasyprint" >/dev/null 2>&1; then
      ok "weasyprint installed into $VENV_PY"
    else
      warn "Could not install it automatically. Run:"
      warn "  $VENV_PY -m pip install weasyprint"
      warn "PDF reports will come back as HTML until then."
    fi
  fi
else
  warn "No venv at $APP_DIR/venv — cannot check weasyprint against the backend's interpreter"
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
  step "Building frontend"
  FRONTEND_CHANGED=1
  BUILD_REASON="forced"
  DIST_INDEX="$APP_DIR/frontend/dist/index.html"

  if command -v rsync >/dev/null 2>&1; then
    if grep -q '^[<>ch]' "$FE_LOG" 2>/dev/null; then
      BUILD_REASON="$(grep -c '^[<>ch]' "$FE_LOG") source file(s) changed"
    elif [[ ! -f "$DIST_INDEX" ]]; then
      BUILD_REASON="no dist/ yet"
    # Belt and braces: even with no change this run, rebuild when anything under
    # src/ is newer than the built bundle. A build that failed, was interrupted,
    # or was skipped by an earlier version of this script leaves exactly that
    # state, and it is invisible until the app misbehaves in the browser.
    elif [[ -n "$(find "$APP_DIR/frontend/src" "$APP_DIR/frontend/index.html" \
                       "$APP_DIR/frontend/package.json" \
                       -newer "$DIST_INDEX" -print -quit 2>/dev/null)" ]]; then
      BUILD_REASON="dist/ older than sources"
    else
      FRONTEND_CHANGED=0
    fi
  fi

  if [[ "$FRONTEND_CHANGED" -eq 1 ]]; then
    info "Rebuilding — $BUILD_REASON"
    cd "$APP_DIR/frontend"
    source "$HOME/.nvm/nvm.sh" 2>/dev/null || true
    nvm use 20 >/dev/null 2>&1 || true
    # Do not let a failed build pass silently: the previous bundle stays on disk
    # and the backend keeps serving it, so a broken build looks like a
    # successful deploy until the browser disagrees.
    if npm run build; then
      ok "Frontend built"
      if [[ -f "$DIST_INDEX" ]]; then
        info "Serving: $(grep -o 'assets/index-[A-Za-z0-9_-]*\.js' "$DIST_INDEX" | head -1)"
      fi
    else
      warn "The browser will keep serving the OLD bundle until this builds."
      die  "FRONTEND BUILD FAILED — backend NOT restarted. Fix the build and re-run."
    fi
    cd "$APP_DIR"
  else
    info "Up to date — dist/ is newer than every frontend source"
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
