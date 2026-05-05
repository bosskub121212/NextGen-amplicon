#!/usr/bin/env bash
# ============================================================
#  16S/12S Amplicon App — Full Setup Script
#  Tested on: Ubuntu 22.04 (WSL2 or native Linux)
# ============================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$SCRIPT_DIR"

echo ""
echo "================================================="
echo "  16S/12S Amplicon Sequencing App — Setup"
echo "================================================="
echo "Install path: $APP_DIR"
echo ""

# ── Check OS ──────────────────────────────────────────────────
if [[ "$OSTYPE" != "linux-gnu"* ]]; then
  echo "ERROR: This script requires Linux (or WSL2 on Windows)."
  echo "  On Windows: open Ubuntu in WSL2, then run this script."
  exit 1
fi

# ── 1. System packages ────────────────────────────────────────
echo "[1/6] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
  curl wget git \
  python3 python3-pip python3-venv \
  r-base r-base-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libgit2-dev libgdal-dev \
  build-essential gfortran
echo "  System packages OK"

# ── 2. Node.js (for building the React frontend) ─────────────
echo ""
echo "[2/6] Checking Node.js..."
if ! command -v node &>/dev/null || [[ "$(node -v | cut -d. -f1 | tr -d 'v')" -lt 18 ]]; then
  echo "  Installing Node.js 20 LTS..."
  curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - -qq
  sudo apt-get install -y -qq nodejs
fi
echo "  Node.js $(node -v) OK"

# ── 3. Python virtual environment ────────────────────────────
echo ""
echo "[3/6] Setting up Python environment..."
cd "$APP_DIR/backend"
python3 -m venv ../venv
source ../venv/bin/activate
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt
deactivate
echo "  Python venv OK (packages: fastapi uvicorn psutil)"

# ── 4. R packages ────────────────────────────────────────────
echo ""
echo "[4/6] Installing R packages (this takes 10–30 min first time)..."

# Bioconductor + DADA2
Rscript - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── Bioconductor manager ──────────────────────────────────────
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", quiet=TRUE)
}

# ── DADA2 (Bioconductor) ──────────────────────────────────────
if (!requireNamespace("dada2", quietly=TRUE)) {
  cat("Installing dada2 (Bioconductor)...\n")
  BiocManager::install("dada2", ask=FALSE, update=FALSE, quiet=TRUE)
} else {
  cat("  dada2 already installed\n")
}

# ── Visualization packages (CRAN) ────────────────────────────
pkgs <- c("optparse","ggplot2","reshape2","vegan","ape","pheatmap","scales","jsonlite")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet=TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# ── Verify ────────────────────────────────────────────────────
cat("\n=== Package check ===\n")
all_pkgs <- c("dada2","optparse","ggplot2","reshape2","vegan","ape","pheatmap","jsonlite")
ok <- TRUE
for (pkg in all_pkgs) {
  status <- if (requireNamespace(pkg, quietly=TRUE)) "OK" else "MISSING"
  if (status == "MISSING") ok <- FALSE
  cat(sprintf("  %-12s %s\n", pkg, status))
}
if (!ok) {
  cat("\nWARNING: Some packages are missing. Re-run setup.sh to retry.\n")
  quit(status=1)
}
cat("\nAll R packages OK\n")
REOF

echo "  R packages OK"

# ── 5. Frontend build ────────────────────────────────────────
echo ""
echo "[5/6] Building React frontend..."
cd "$APP_DIR/frontend"
if [[ ! -f "package.json" ]]; then
  echo "  ERROR: frontend/package.json not found. Cannot build."
  exit 1
fi
npm install --silent
npm run build --silent
echo "  Frontend built → frontend/dist/"

# ── 6. Directories & database placeholder ────────────────────
echo ""
echo "[6/6] Creating required directories..."
mkdir -p "$APP_DIR/backend/databases"
mkdir -p "$APP_DIR/uploads"
mkdir -p "$APP_DIR/results"

# ── Database check ────────────────────────────────────────────
DB_COUNT=$(find "$APP_DIR/backend/databases" -name "*.fa.gz" -o -name "*.fasta.gz" 2>/dev/null | wc -l)
if [[ "$DB_COUNT" -eq 0 ]]; then
  echo ""
  echo "  ⚠  No taxonomy database found in backend/databases/"
  echo "  Download SILVA 138.2 (recommended for 16S):"
  echo ""
  echo "    cd $APP_DIR/backend/databases"
  echo "    wget https://zenodo.org/records/10403693/files/silva_nr99_v138.2_train_set.fa.gz"
  echo "    wget https://zenodo.org/records/10403693/files/silva_species_assignment_v138.2.fa.gz"
  echo ""
  echo "  Or for 12S (MitoFish):"
  echo "    # Place your custom database .fa.gz here"
fi

# ── Done ──────────────────────────────────────────────────────
echo ""
echo "================================================="
echo "  Setup complete!"
echo "================================================="
echo ""
echo "  Start the app:   bash $APP_DIR/start.sh"
echo "  Then open:       http://localhost:5173"
echo ""
if [[ "$DB_COUNT" -gt 0 ]]; then
  echo "  Databases found: $DB_COUNT file(s)"
fi
echo ""
