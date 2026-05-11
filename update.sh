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
git pull origin main
echo "  ✓ Code updated"
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
# STEP 3: Check if new R packages are needed
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 3: Checking R packages ─────────────────────────────"
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
