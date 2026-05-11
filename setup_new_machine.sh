#!/bin/bash
# =============================================================================
#  NextGen-Amplicon — Full Setup Script for New Machine (WSL Ubuntu)
#  Run once after cloning the repo:
#
#    cd ~/r16s-app
#    bash setup_new_machine.sh
#
# =============================================================================

set -e  # stop on error
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   NextGen-Amplicon — New Machine Setup                  ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: System packages (apt)
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 1: Installing system libraries ────────────────────"
sudo apt-get update -qq
sudo apt-get install -y \
    cmake \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libuv1-dev \
    libcairo2-dev \
    libgsl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libfreetype6-dev \
    python3-pip \
    python3-venv \
    nodejs \
    npm \
    cutadapt
echo "  ✓ System libraries installed"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Python virtual environment + FastAPI backend
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 2: Python backend (FastAPI) ────────────────────────"
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo "  ✓ Created venv"
fi
source venv/bin/activate
pip install --upgrade pip -q
pip install fastapi uvicorn python-multipart aiofiles requests -q
echo "  ✓ Python packages installed"
deactivate
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Node.js frontend build
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 3: Frontend (React/Vite) ───────────────────────────"
cd frontend
npm install --silent
npm run build --silent
cd ..
echo "  ✓ Frontend built"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: R packages
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 4: R packages ──────────────────────────────────────"
if ! command -v Rscript &> /dev/null; then
    echo "  [!] R not found — installing R 4.x ..."
    sudo apt-get install -y r-base r-base-dev
fi
Rscript backend/r_scripts/install_packages.R
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Create systemd startup service (optional, auto-start on boot)
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 5: Startup helper ──────────────────────────────────"
cat > start_backend.sh << 'EOF'
#!/bin/bash
# Quick-start script for the backend
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
EOF
chmod +x start_backend.sh
echo "  ✓ Created start_backend.sh"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Setup complete!                                        ║"
echo "║                                                          ║"
echo "║   To start the backend:                                  ║"
echo "║     bash ~/r16s-app/start_backend.sh                    ║"
echo "║                                                          ║"
echo "║   To update later (after git pull):                     ║"
echo "║     bash ~/r16s-app/update.sh                           ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
