#!/bin/bash
# =============================================================================
#  NextGen-Amplicon — Update Script
#  Run after each new release to pull code + rebuild frontend:
#
#    bash ~/r16s-app/update.sh
#
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   NextGen-Amplicon — Updating                           ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Pull latest code
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 1: Pulling latest code ─────────────────────────────"

# jobs_history.json is local-only data (job queue/history per machine).
# If git is still tracking it (was committed before .gitignore was added),
# we preserve the content, reset the tracked copy, pull, then restore.
JOBS_FILE="backend/jobs_history.json"
JOBS_BACKUP=""
if git ls-files --error-unmatch "$JOBS_FILE" &>/dev/null; then
  echo "  ℹ  Preserving local job history across update..."
  JOBS_BACKUP=$(cat "$JOBS_FILE" 2>/dev/null || echo "")
  git checkout -- "$JOBS_FILE" 2>/dev/null || true
  # Permanently stop tracking this file so this never happens again
  git rm --cached "$JOBS_FILE" 2>/dev/null || true
fi

git pull origin main
echo "  ✓ Code updated"

# Restore job history if we had one
if [ -n "$JOBS_BACKUP" ]; then
  echo "$JOBS_BACKUP" > "$JOBS_FILE"
  echo "  ✓ Job history restored"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Rebuild frontend
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 2: Rebuilding frontend ─────────────────────────────"
cd frontend
npm install --silent
npm run build --silent
cd ..
echo "  ✓ Frontend rebuilt"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Ensure start_backend.sh exists and is executable
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 3: Checking helper scripts ────────────────────────"
if [[ ! -f "$SCRIPT_DIR/start_backend.sh" ]]; then
  cat > "$SCRIPT_DIR/start_backend.sh" <<'STARTSH'
#!/usr/bin/env bash
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
STARTSH
  echo "  ✓ start_backend.sh created"
fi
chmod +x "$SCRIPT_DIR/start_backend.sh"
echo "  ✓ start_backend.sh ready"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Check if new R packages are needed
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 4: Checking R packages ─────────────────────────────"
Rscript backend/r_scripts/install_packages.R
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Update complete!                                       ║"
echo "║                                                          ║"
echo "║   Restart the backend to apply changes:                  ║"
echo "║     bash ~/r16s-app/start_backend.sh                    ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
