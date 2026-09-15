#!/usr/bin/env bash
# ============================================================
#  16S/12S Amplicon App — Local Machine Start Script
#  Backend : port 8000
#  Frontend: port 3000
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

PYTHON="$APP_DIR/venv/bin/python3"
UVICORN="$APP_DIR/venv/bin/uvicorn"

if [[ ! -f "$UVICORN" ]]; then
  echo "ERROR: venv not found. Run setup.sh first."
  exit 1
fi

# Kill old instances
pkill -f "uvicorn main:app" 2>/dev/null || true
pkill -f "npx.*serve"       2>/dev/null || true
pkill -f "vite"             2>/dev/null || true
sleep 1

echo ""
echo "================================================="
echo "  AmpliconApp — Local (ports 8000 / 3000)"
echo "================================================="

# ── Backend ───────────────────────────────────────────────────
echo "  Starting backend on port 8000..."
cd "$APP_DIR/backend"
"$PYTHON" -m uvicorn main:app --host 0.0.0.0 --port 8000 &
BACKEND_PID=$!

# ── Frontend ──────────────────────────────────────────────────
if [[ -d "$APP_DIR/frontend/dist" ]]; then
  echo "  Starting frontend (production) on port 3000..."
  cd "$APP_DIR/frontend"
  npx --yes serve -s dist -l 3000 &
  FRONTEND_PID=$!
else
  echo "  Starting frontend (dev/Vite) on port 3000..."
  cd "$APP_DIR/frontend"
  npm run dev -- --port 3000 &
  FRONTEND_PID=$!
fi

sleep 2
echo ""
echo "  ✅ App is running!"
echo "  Open: http://localhost:3000"
echo ""
echo "  Press Ctrl+C to stop."
echo ""

trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; echo '  Stopped.'; exit 0" INT TERM
wait
