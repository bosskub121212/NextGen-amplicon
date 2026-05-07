#!/usr/bin/env bash
# ============================================================
#  16S/12S Amplicon App — Start Script
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

# ── Find built frontend folder (Vite = dist/, CRA = build/) ──
if [[ -d "$APP_DIR/frontend/dist" ]]; then
  FRONTEND_BUILD="$APP_DIR/frontend/dist"
elif [[ -d "$APP_DIR/frontend/build" ]]; then
  FRONTEND_BUILD="$APP_DIR/frontend/build"
else
  echo "ERROR: Frontend not built. Run:"
  echo "  cd $APP_DIR/frontend && npm install && npm run build"
  exit 1
fi

# ── Use full path to venv binaries (no activate needed) ───────
UVICORN="$APP_DIR/venv/bin/uvicorn"
PYTHON="$APP_DIR/venv/bin/python3"

if [[ ! -f "$UVICORN" ]]; then
  echo "ERROR: uvicorn not found in venv. Run:"
  echo "  cd $APP_DIR && python3 -m venv venv"
  echo "  $APP_DIR/venv/bin/pip install -r backend/requirements.txt"
  exit 1
fi

# ── Kill any existing instances ────────────────────────────────
pkill -f "uvicorn main:app" 2>/dev/null || true
pkill -f "npx.*serve\|node.*serve" 2>/dev/null || true
sleep 1

echo ""
echo "================================================="
echo "  Starting 16S/12S Amplicon App"
echo "================================================="

# ── Backend (python -m uvicorn avoids permission issues) ──────
echo "  Starting backend on port 8000..."
cd "$APP_DIR/backend"
"$PYTHON" -m uvicorn main:app --host 0.0.0.0 --port 8000 &
BACKEND_PID=$!

# ── Frontend (use npx serve — no global install needed) ───────
echo "  Starting frontend on port 5173..."
cd "$APP_DIR/frontend"
npx --yes serve -s "$FRONTEND_BUILD" -l 5173 &
FRONTEND_PID=$!

sleep 2
echo ""
echo "  ✅ App is running!"
echo "  Open: http://localhost:5173"
echo ""
echo "  Press Ctrl+C to stop."
echo ""

trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; echo '  Stopped.'; exit 0" INT TERM
wait
