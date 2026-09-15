#!/bin/bash
echo "=============================="
echo "  AmpliconApp Backend :8000"
echo "=============================="
echo ""
cd ~/r16s-app
source venv/bin/activate
uvicorn backend.main:app --reload --port 8000
echo ""
echo "Backend stopped. Press Enter to close..."
read
