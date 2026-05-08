#!/usr/bin/env bash
# ============================================================
#  NextGen-Amplicon — Start Script
#  Usage: bash ~/r16s-app/start.sh
#
#  Opens the app at http://localhost:8000
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

UVICORN="$APP_DIR/venv/bin/uvicorn"
PYTHON="$APP_DIR/venv/bin/python3"

if [[ ! -f "$UVICORN" ]]; then
  echo "ERROR: venv not found. Run setup first:"
  echo "  bash $APP_DIR/setup.sh"
  exit 1
fi

# ── Check frontend is built ────────────────────────────────
FRONTEND_BUILT=0
for d in "$APP_DIR/frontend/dist" "$APP_DIR/frontend/build"; do
  [[ -f "$d/index.html" ]] && FRONTEND_BUILT=1 && break
done
if [[ $FRONTEND_BUILT -eq 0 ]]; then
  echo "Frontend not built — building now..."
  cd "$APP_DIR/frontend"
  npm install --silent
  npm run build
  cd "$APP_DIR"
fi

# ── Kill any existing instance ─────────────────────────────
pkill -f "uvicorn backend.main:app" 2>/dev/null || true
sleep 1

echo ""
echo "================================================="
echo "  NextGen-Amplicon"
echo "  http://localhost:8000"
echo "================================================="
echo ""

cd "$APP_DIR"
"$UVICORN" backend.main:app --host 0.0.0.0 --port 8000

# (blocking — Ctrl+C to stop)
