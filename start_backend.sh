#!/usr/bin/env bash
# =============================================================================
#  NextGen-Amplicon — Start Backend
#  Run: bash ~/r16s-app/start_backend.sh
# =============================================================================
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
