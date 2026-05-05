#!/usr/bin/env bash
# ============================================================
#  16S/12S App — Development Mode (hot-reload frontend)
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

pkill -f "uvicorn main:app" 2>/dev/null || true
pkill -f "vite" 2>/dev/null || true
sleep 1

echo "Starting in DEVELOPMENT mode (hot-reload)..."

# Backend
cd "$APP_DIR/backend"
source "$APP_DIR/venv/bin/activate"
uvicorn main:app --host 0.0.0.0 --port 8000 --reload &
BACKEND_PID=$!

# Frontend (Vite dev server)
cd "$APP_DIR/frontend"
npm run dev &
FRONTEND_PID=$!

echo "  Backend:  http://localhost:8000"
echo "  Frontend: http://localhost:5173"
echo "  Press Ctrl+C to stop."

trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; exit 0" INT TERM
wait
