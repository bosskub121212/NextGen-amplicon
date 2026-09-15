#!/usr/bin/env bash
# ============================================================
#  pack.sh — Sync latest source files then create zip
#  Run from WSL: bash /mnt/c/Claude/r16s-app/pack.sh
# ============================================================
set -e

WIN_SRC="/mnt/c/Claude/r16s-app"   # Windows workspace (source of truth for src/)
WSL_APP="/home/boss/r16s-app"       # WSL copy (has full frontend project)
OUT_DIR="/mnt/c/Claude"             # Output zip goes to C:\Claude

echo "================================================="
echo "  16S/12S Amplicon App — Pack for Distribution"
echo "================================================="

# ── 1. Sync updated source files Windows → WSL ────────────────
echo ""
echo "[1/3] Syncing latest source files to WSL..."

# Backend
cp "$WIN_SRC/backend/main.py"                    "$WSL_APP/backend/main.py"
cp "$WIN_SRC/backend/requirements.txt"           "$WSL_APP/backend/requirements.txt"
cp "$WIN_SRC/backend/r_scripts/dada2_pipeline.R" "$WSL_APP/backend/r_scripts/dada2_pipeline.R"
cp "$WIN_SRC/backend/r_scripts/replot.R"         "$WSL_APP/backend/r_scripts/replot.R"

# Frontend src
cp "$WIN_SRC/frontend/src/App.tsx"               "$WSL_APP/frontend/src/App.tsx"
cp "$WIN_SRC/frontend/src/App.css"               "$WSL_APP/frontend/src/App.css"
cp "$WIN_SRC/frontend/src/components/DNAProgress.tsx"        "$WSL_APP/frontend/src/components/DNAProgress.tsx"
cp "$WIN_SRC/frontend/src/components/MetadataEditor.tsx"     "$WSL_APP/frontend/src/components/MetadataEditor.tsx"
cp "$WIN_SRC/frontend/src/components/PipelineSettings.tsx"   "$WSL_APP/frontend/src/components/PipelineSettings.tsx"
cp "$WIN_SRC/frontend/src/components/TaxonomyColorPicker.tsx" "$WSL_APP/frontend/src/components/TaxonomyColorPicker.tsx"

# Scripts + docs
cp "$WIN_SRC/setup.sh"              "$WSL_APP/setup.sh"
cp "$WIN_SRC/start.sh"              "$WSL_APP/start.sh"
cp "$WIN_SRC/start_dev.sh"          "$WSL_APP/start_dev.sh"
cp "$WIN_SRC/INSTALL.md"            "$WSL_APP/INSTALL.md"
cp "$WIN_SRC/install_r_packages.R"  "$WSL_APP/install_r_packages.R"

chmod +x "$WSL_APP/setup.sh" "$WSL_APP/start.sh" "$WSL_APP/start_dev.sh"
echo "  Sync done."

# ── 2. Build frontend (so the zip has a ready-to-run build) ───
echo ""
echo "[2/3] Building React frontend..."
cd "$WSL_APP/frontend"
npm install --silent
npm run build --silent
echo "  Build done → frontend/dist/"

# ── 3. Create zip ─────────────────────────────────────────────
echo ""
echo "[3/3] Creating zip archive..."
DATE=$(date +%Y%m%d_%H%M)
OUT="$OUT_DIR/amplicon_app_${DATE}.zip"

cd "$WSL_APP"
zip -r "$OUT" . \
  --exclude "*.git*" \
  --exclude "*node_modules*" \
  --exclude "*__pycache__*" \
  --exclude "*/venv/*" \
  --exclude "*/databases/*.fa.gz" \
  --exclude "*/databases/*.fasta.gz" \
  --exclude "*/results/*" \
  --exclude "*/uploads/*" \
  --exclude "*.pyc" \
  --exclude "*/jobs_history.json" \
  -q

SIZE=$(du -sh "$OUT" | cut -f1)
echo ""
echo "================================================="
echo "  Done! File ready at:"
echo "  C:\\Claude\\amplicon_app_${DATE}.zip  ($SIZE)"
echo "================================================="
echo ""
echo "  Contents:"
echo "   ✅ Backend (main.py, dada2_pipeline.R, replot.R)"
echo "   ✅ Frontend source + built dist/"
echo "   ✅ setup.sh, start.sh, INSTALL.md"
echo "   ⚠  Databases NOT included — recipient runs:"
echo "      wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
echo ""
