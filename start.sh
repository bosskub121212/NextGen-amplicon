#!/usr/bin/env bash
# ============================================================
#  16S/12S Amplicon App — Start Script
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

# ── Checks ────────────────────────────────────────────────────
if [[ ! -d "$APP_DIR/venv" ]]; then
  echo "ERROR: Python venv not found. Run setup.sh first."
  exit 1
fi
if [[ ! -d "$APP_DIR/frontend/dist" ]]; then
  echo "ERROR: Frontend not built. Run setup.sh first."
  exit 1
fi

# ── Kill any existing instances ────────────────────────────────
pkill -f "uvicorn main:app" 2>/dev/null || true
pkill -f "vite.*5173" 2>/dev/null || true
sleep 1

echo ""
echo "================================================="
echo "  Starting 16S/12S Amplicon App"
echo "================================================="

# ── Backend ───────────────────────────────────────────────────
cd "$APP_DIR/backend"
source "$APP_DIR/venv/bin/activate"
echo "  Starting backend (port 8000)..."
uvicorn main:app --host 0.0.0.0 --port 8000 &
BACKEND_PID=$!

# ── Frontend (serve built files) ──────────────────────────────
cd "$APP_DIR/frontend"
echo "  Starting frontend (port 5173)..."
npx --yes serve dist -l 5173 &
FRONTEND_PID=$!

sleep 2
echo ""
echo "  App is running!"
echo "  Open: http://localhost:5173"
echo ""
echo "  Press Ctrl+C to stop."
echo ""

# ── Wait and cleanup ─────────────────────────────────────────
trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; echo '  Stopped.'; exit 0" INT TERM
wait
