#!/usr/bin/env bash
# =============================================================
#  NextGen-Amplicon -- Copy All Extension Files to WSL Runtime
#
#  Run from WSL:
#    bash /mnt/c/Claude/r16s-app/copy_all_extensions.sh
#
#  This syncs every file edited by Claude to ~/r16s-app/
#  which is where the backend actually runs from.
# =============================================================
set -euo pipefail

SRC="/mnt/c/Claude/r16s-app"
DST="$HOME/r16s-app"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  OK${NC}  $1"; }
skip() { echo -e "${YELLOW}SKIP${NC}  $1"; }
fail() { echo -e "${RED}FAIL${NC}  $1"; }

echo ""
echo "============================================================"
echo "  NextGen-Amplicon — Sync Files to WSL Runtime"
echo "  Source : $SRC"
echo "  Dest   : $DST"
echo "============================================================"
echo ""

if [[ ! -d "$SRC" ]]; then
  echo "ERROR: Source not found: $SRC"
  echo "Make sure Windows drive C: is mounted in WSL."
  exit 1
fi

if [[ ! -d "$DST" ]]; then
  echo "ERROR: Runtime directory not found: $DST"
  echo "Run the original install script first."
  exit 1
fi

copy_file() {
  local rel="$1"          # relative path from SRC/DST
  local src_f="$SRC/$rel"
  local dst_f="$DST/$rel"
  if [[ ! -f "$src_f" ]]; then
    skip "$rel  (source missing)"
    return
  fi
  mkdir -p "$(dirname "$dst_f")"
  cp "$src_f" "$dst_f"
  ok "$rel"
}

echo "--- Backend Python ---"
copy_file "backend/main.py"
copy_file "backend/updater.py"
copy_file "backend/license.py"
copy_file "backend/picrust2_pipeline.py"

echo ""
echo "--- Tools ---"
copy_file "tools/keygen.py"
copy_file "tools/run_keygen.bat"

echo ""
echo "--- R Pipeline Scripts ---"
copy_file "backend/r_scripts/dada2_pipeline.R"
copy_file "backend/r_scripts/its_pipeline.R"
copy_file "backend/r_scripts/cox1_pipeline.R"
copy_file "backend/r_scripts/pacbio_pipeline.R"

echo ""
echo "--- Frontend Source ---"
copy_file "frontend/src/App.tsx"
copy_file "frontend/src/App.css"
copy_file "frontend/src/components/PipelineSettings.tsx"
copy_file "frontend/src/components/UpdateBanner.tsx"
copy_file "frontend/src/components/LicenseModal.tsx"
copy_file "frontend/src/components/TaxonomyColorPicker.tsx"
copy_file "frontend/src/components/MetadataEditor.tsx"
copy_file "frontend/src/components/DNAProgress.tsx"
copy_file "frontend/src/components/LicenseModal.tsx"
copy_file "frontend/src/App.css"

echo ""
echo "--- Utility / Config ---"
copy_file "version.json"
copy_file "install.sh"
copy_file "setup.sh"
copy_file "download_databases.sh"
copy_file ".dev_bypass" 2>/dev/null || true   # may not exist on all machines

echo ""
echo "============================================================"
echo "  Files synced. Now restart backend + rebuild frontend:"
echo ""
echo "  # Restart backend (if running with uvicorn --reload it's auto)"
echo "  pkill -f 'uvicorn main:app' 2>/dev/null; sleep 1"
echo "  cd ~/r16s-app && source venv/bin/activate"
echo "  uvicorn backend.main:app --reload --host 0.0.0.0 --port 8000 &"
echo ""
echo "  # Rebuild frontend (run only if TSX/CSS changed)"
echo "  cd ~/r16s-app/frontend && npm run build"
echo "  # OR for dev hot-reload:"
echo "  cd ~/r16s-app/frontend && npm run dev &"
echo "============================================================"
echo ""
