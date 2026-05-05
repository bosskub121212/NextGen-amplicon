#!/usr/bin/env bash
# =============================================================
#  NextGen-Amplicon — Full Machine Installer
#  Run this script INSIDE WSL (Ubuntu 22.04 recommended):
#
#    bash install.sh
#
#  What this script does:
#    1. Install system packages (R 4.4, Python 3, Node.js, etc.)
#    2. Clone the GitHub repo (or use existing copy)
#    3. Create Python virtual environment + install packages
#    4. Install all required R packages (DADA2, etc.)
#    5. Build the React frontend
#    6. Create start.sh / stop.sh helper scripts
#    7. (Optional) Configure systemd user auto-start
#
#  Requirements:
#    - WSL 2 with Ubuntu 22.04 or later
#    - Internet access (for package downloads)
#    - GitHub PAT with repo read access (for private repo)
# =============================================================
set -euo pipefail

# ── Colors ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YELLOW}  ⚠${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

# ── Config ────────────────────────────────────────────────────────────────────
REPO_URL="https://github.com/bosskub121212/NextGen-amplicon.git"
INSTALL_DIR="$HOME/r16s-app"
VENV_DIR="$INSTALL_DIR/venv"
NODE_VERSION="20"
R_VERSION_MAJOR="4"

# ── Banner ────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${CYAN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║       NextGen-Amplicon — Full Machine Installer          ║"
echo "  ║       16S / ITS / COX1 / 18S-nema / PacBio              ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"
echo "  Install directory : $INSTALL_DIR"
echo "  Repository        : $REPO_URL"
echo ""

# ── Sanity checks ─────────────────────────────────────────────────────────────
step "Pre-flight checks"

if [[ "$(uname -r)" != *microsoft* ]] && [[ "$(uname -r)" != *WSL* ]]; then
  warn "Not detected as WSL — this script is designed for WSL2/Ubuntu."
  warn "Continuing anyway, but some steps may need adjustment."
fi

if [[ "$EUID" -eq 0 ]]; then
  die "Do not run as root. Run as your normal user inside WSL."
fi

ok "Running as user: $(whoami)"

# ── GitHub PAT ────────────────────────────────────────────────────────────────
step "GitHub access token"

TOKEN_FILE="$HOME/.config/amplicon/github_token"
mkdir -p "$(dirname "$TOKEN_FILE")"

if [[ -f "$TOKEN_FILE" ]] && [[ -n "$(cat "$TOKEN_FILE" 2>/dev/null)" ]]; then
  GITHUB_TOKEN="$(cat "$TOKEN_FILE")"
  ok "Existing token loaded from $TOKEN_FILE"
else
  echo -e "  ${YELLOW}A GitHub Personal Access Token (PAT) is needed to:"
  echo "    - Clone the private NextGen-Amplicon repository"
  echo "    - Receive automatic updates in the future"
  echo ""
  echo -e "  Create one at: https://github.com/settings/tokens${NC}"
  echo "  Required scope: repo (read access)"
  echo ""
  read -rsp "  Paste your GitHub PAT (input hidden): " GITHUB_TOKEN
  echo ""
  if [[ -z "$GITHUB_TOKEN" ]]; then
    die "GitHub token is required."
  fi
  echo "$GITHUB_TOKEN" > "$TOKEN_FILE"
  chmod 600 "$TOKEN_FILE"
  ok "Token saved to $TOKEN_FILE"
fi

# Embed token in repo URL for authentication
REPO_AUTH_URL="https://${GITHUB_TOKEN}@github.com/bosskub121212/NextGen-amplicon.git"

# ── System packages ───────────────────────────────────────────────────────────
step "System packages (apt)"
sudo apt-get update -qq

PKGS=(
  git curl wget gnupg2 ca-certificates software-properties-common
  build-essential gfortran
  libcurl4-openssl-dev libssl-dev libxml2-dev
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
  libgsl-dev libglpk-dev libgmp3-dev zlib1g-dev
  libopenblas-dev liblapack-dev
  python3 python3-pip python3-venv python3-dev
  nodejs npm
  procps htop
)
sudo apt-get install -y "${PKGS[@]}" -qq
ok "System packages installed"

# ── R 4.4 ─────────────────────────────────────────────────────────────────────
step "R $R_VERSION_MAJOR.x (CRAN)"

R_INSTALLED=false
if command -v R &>/dev/null; then
  R_VER=$(R --version | head -1 | grep -oP '\d+\.\d+' | head -1)
  if [[ "${R_VER%%.*}" -ge "$R_VERSION_MAJOR" ]]; then
    ok "R $R_VER already installed"
    R_INSTALLED=true
  fi
fi

if [[ "$R_INSTALLED" == "false" ]]; then
  info "Installing R from CRAN..."
  # Add CRAN repo
  sudo apt-get install -y --no-install-recommends dirmngr -qq
  UBUNTU_CODENAME=$(lsb_release -cs 2>/dev/null || echo "jammy")
  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${UBUNTU_CODENAME}-cran40/" \
    | sudo tee /etc/apt/sources.list.d/cran_r.list >/dev/null
  sudo apt-get update -qq
  sudo apt-get install -y r-base r-base-dev -qq
  ok "R $(R --version | head -1 | grep -oP '\d+\.\d+\.\d+' | head -1) installed"
fi

# ── Node.js (nvm) ──────────────────────────────────────────────────────────────
step "Node.js $NODE_VERSION (via nvm)"

export NVM_DIR="$HOME/.nvm"
if [[ ! -d "$NVM_DIR" ]]; then
  info "Installing nvm..."
  curl -fsSL https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash
  ok "nvm installed"
fi

# Load nvm
# shellcheck source=/dev/null
source "$NVM_DIR/nvm.sh" 2>/dev/null || true

if ! nvm ls "$NODE_VERSION" &>/dev/null; then
  info "Installing Node.js $NODE_VERSION..."
  nvm install "$NODE_VERSION"
  nvm alias default "$NODE_VERSION"
fi
nvm use "$NODE_VERSION" >/dev/null 2>&1
ok "Node.js $(node --version) / npm $(npm --version)"

# ── Clone / update repo ────────────────────────────────────────────────────────
step "Repository"

if [[ -d "$INSTALL_DIR/.git" ]]; then
  info "Existing repo found — pulling latest..."
  cd "$INSTALL_DIR"
  git remote set-url origin "$REPO_AUTH_URL"
  git pull origin main --ff-only || warn "git pull failed — continuing with local version"
  ok "Repository up to date"
else
  info "Cloning repository to $INSTALL_DIR..."
  git clone "$REPO_AUTH_URL" "$INSTALL_DIR"
  ok "Repository cloned"
fi

cd "$INSTALL_DIR"

# Persist token for future updates
mkdir -p "$INSTALL_DIR/backend"
echo "$GITHUB_TOKEN" > "$TOKEN_FILE"
chmod 600 "$TOKEN_FILE"

# ── Python venv + packages ─────────────────────────────────────────────────────
step "Python environment"

if [[ ! -d "$VENV_DIR" ]]; then
  python3 -m venv "$VENV_DIR"
  ok "Virtual environment created"
fi

source "$VENV_DIR/bin/activate"

pip install --upgrade pip -q
pip install \
  fastapi "uvicorn[standard]" \
  pydantic python-multipart \
  requests psutil \
  pandas numpy \
  -q
ok "Core Python packages installed"

# PICRUSt2 (optional — heavy dependency, skip if conda not available)
if command -v conda &>/dev/null; then
  info "Conda found — installing PICRUSt2..."
  conda install -c bioconda -c conda-forge picrust2 -y -q 2>/dev/null || \
    warn "PICRUSt2 conda install failed — skipping (can be added later)"
  ok "PICRUSt2 installed via conda"
else
  warn "conda not found — PICRUSt2 skipped (install miniconda3 to enable)"
fi

deactivate
ok "Python environment ready"

# ── R packages ────────────────────────────────────────────────────────────────
step "R packages"

info "Installing CRAN + Bioconductor packages (this takes 15–40 min first time)..."

Rscript - <<'RSCRIPT'
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── Helper: install if missing ──────────────────────────────────────────────
install_if_missing <- function(pkgs, from = "CRAN") {
  need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(need) == 0) { cat("  All present:", paste(pkgs, collapse=", "), "\n"); return(invisible()) }
  cat("  Installing:", paste(need, collapse=", "), "\n")
  if (from == "CRAN")   install.packages(need, quiet = TRUE, dependencies = TRUE)
  if (from == "BiocMgr") BiocManager::install(need, ask = FALSE, update = FALSE)
}

# ── CRAN packages ───────────────────────────────────────────────────────────
cran_pkgs <- c(
  "optparse", "dplyr", "ggplot2", "reshape2", "gridExtra",
  "vegan", "ape", "phangorn", "scales", "RColorBrewer", "viridis",
  "jsonlite", "data.table", "stringr", "forcats", "tibble", "purrr",
  "grid", "cowplot", "ggrepel", "ggdendro",
  "seqinr", "Biostrings"
)
install_if_missing(cran_pkgs)

# ── BiocManager ────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

bioc_pkgs <- c(
  "dada2", "ShortRead", "Biostrings", "phyloseq",
  "DESeq2", "edgeR", "limma"
)
install_if_missing(bioc_pkgs, from = "BiocMgr")

# ── cutadapt (ITS primer trimming) ─────────────────────────────────────────
if (!requireNamespace("dada2", quietly = TRUE))
  warning("dada2 not installed — ITS primer trimming will fail")

# ── Tax4Fun2 (functional prediction from 16S) ──────────────────────────────
if (!requireNamespace("Tax4Fun2", quietly = TRUE)) {
  zenodo_url <- "https://zenodo.org/record/10035668/files/Tax4Fun2_1.1.5.tar.gz"
  tmp <- tempfile(fileext = ".tar.gz")
  tryCatch({
    download.file(zenodo_url, tmp, quiet = TRUE, method = "curl",
                  extra = "--retry 3 --max-time 120")
    install.packages(tmp, repos = NULL, type = "source")
    cat("  Tax4Fun2 installed from Zenodo\n")
  }, error = function(e) {
    warning(paste("Tax4Fun2 install failed:", e$message))
  })
}

cat("\n  R packages done.\n")
RSCRIPT

ok "R packages installed"

# ── Cutadapt (for ITS trimming) ────────────────────────────────────────────────
step "Cutadapt"
source "$VENV_DIR/bin/activate"
if ! command -v cutadapt &>/dev/null; then
  pip install cutadapt -q
  ok "cutadapt $(cutadapt --version) installed"
else
  ok "cutadapt already available: $(cutadapt --version)"
fi
deactivate

# ── Frontend build ─────────────────────────────────────────────────────────────
step "Frontend (React build)"

cd "$INSTALL_DIR/frontend"
source "$NVM_DIR/nvm.sh" 2>/dev/null || true
nvm use "$NODE_VERSION" >/dev/null 2>&1

info "npm install..."
npm install --legacy-peer-deps -q
info "npm build..."
npm run build
ok "Frontend built successfully"

cd "$INSTALL_DIR"

# ── Directory structure ────────────────────────────────────────────────────────
step "Directory structure"
mkdir -p "$INSTALL_DIR/backend/uploads"
mkdir -p "$INSTALL_DIR/backend/results"
mkdir -p "$INSTALL_DIR/databases/SILVA"
mkdir -p "$INSTALL_DIR/databases/UNITE"
mkdir -p "$INSTALL_DIR/databases/MIDORI2"
mkdir -p "$INSTALL_DIR/databases/NemaBase"
mkdir -p "$INSTALL_DIR/logs"

# Create db_paths.json if missing
DB_PATHS="$INSTALL_DIR/backend/databases/db_paths.json"
mkdir -p "$(dirname "$DB_PATHS")"
if [[ ! -f "$DB_PATHS" ]]; then
  cat > "$DB_PATHS" <<'JSON'
{
  "silva":          "",
  "unite":          "",
  "midori2_co1":    "",
  "midori2_co1_sp": "",
  "nemabase":       ""
}
JSON
  ok "db_paths.json created (edit paths after downloading databases)"
fi

# ── Helper scripts: start.sh / stop.sh ────────────────────────────────────────
step "Helper scripts"

cat > "$INSTALL_DIR/start.sh" <<STARTSH
#!/usr/bin/env bash
# Start NextGen-Amplicon backend (uvicorn) and frontend (vite dev server)
set -euo pipefail

INSTALL_DIR="\$HOME/r16s-app"
VENV_DIR="\$INSTALL_DIR/venv"
NVM_DIR="\$HOME/.nvm"

# Kill any old instances
pkill -f "uvicorn backend.main" 2>/dev/null || true
pkill -f "vite"                  2>/dev/null || true
sleep 1

# Start backend
source "\$VENV_DIR/bin/activate"
cd "\$INSTALL_DIR"
nohup uvicorn backend.main:app --host 0.0.0.0 --port 8000 --reload \
  > "\$INSTALL_DIR/logs/backend.log" 2>&1 &
echo "Backend started (PID \$!)"
deactivate

# Start frontend dev server
source "\$NVM_DIR/nvm.sh" 2>/dev/null || true
nvm use 20 >/dev/null 2>&1
cd "\$INSTALL_DIR/frontend"
nohup npm run dev > "\$INSTALL_DIR/logs/frontend.log" 2>&1 &
echo "Frontend started (PID \$!)"

echo ""
echo "  ✅ NextGen-Amplicon is starting..."
echo "  Backend : http://localhost:8000"
echo "  Frontend: http://localhost:5173"
echo ""
echo "  Logs:"
echo "    tail -f \$INSTALL_DIR/logs/backend.log"
echo "    tail -f \$INSTALL_DIR/logs/frontend.log"
STARTSH
chmod +x "$INSTALL_DIR/start.sh"

cat > "$INSTALL_DIR/stop.sh" <<STOPSH
#!/usr/bin/env bash
# Stop NextGen-Amplicon backend and frontend
pkill -f "uvicorn backend.main" 2>/dev/null && echo "Backend stopped" || echo "Backend was not running"
pkill -f "vite"                  2>/dev/null && echo "Frontend stopped" || echo "Frontend was not running"
STOPSH
chmod +x "$INSTALL_DIR/stop.sh"
ok "start.sh / stop.sh created"

# ── (Optional) systemd user service ───────────────────────────────────────────
step "Systemd auto-start (optional)"

SYSTEMD_DIR="$HOME/.config/systemd/user"
SERVICE_FILE="$SYSTEMD_DIR/nextgen-amplicon.service"
mkdir -p "$SYSTEMD_DIR"

cat > "$SERVICE_FILE" <<SVCEOF
[Unit]
Description=NextGen-Amplicon Amplicon Analysis Backend
After=network.target

[Service]
Type=simple
WorkingDirectory=%h/r16s-app
ExecStart=%h/r16s-app/venv/bin/uvicorn backend.main:app --host 0.0.0.0 --port 8000
Restart=on-failure
RestartSec=5
StandardOutput=append:%h/r16s-app/logs/backend.log
StandardError=append:%h/r16s-app/logs/backend.log
Environment=PATH=%h/r16s-app/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

[Install]
WantedBy=default.target
SVCEOF

if command -v systemctl &>/dev/null; then
  systemctl --user daemon-reload 2>/dev/null || true
  systemctl --user enable nextgen-amplicon.service 2>/dev/null || true
  ok "Systemd service enabled (nextgen-amplicon.service)"
  info "To start now: systemctl --user start nextgen-amplicon"
else
  warn "systemd not available in this WSL session — use start.sh manually"
fi

# ── Add shell aliases ──────────────────────────────────────────────────────────
step "Shell aliases"

ALIAS_BLOCK='
# NextGen-Amplicon shortcuts
alias ngamp-start="bash ~/r16s-app/start.sh"
alias ngamp-stop="bash ~/r16s-app/stop.sh"
alias ngamp-log="tail -f ~/r16s-app/logs/backend.log"'

for RC in ~/.bashrc ~/.zshrc; do
  if [[ -f "$RC" ]] && ! grep -q "ngamp-start" "$RC" 2>/dev/null; then
    echo "$ALIAS_BLOCK" >> "$RC"
    ok "Aliases added to $RC"
  fi
done

# ── Done ───────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║   ✅  Installation complete!                             ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"
echo "  Next steps:"
echo ""
echo "  1. Download reference databases (optional — needed for taxonomy):"
echo "     bash ~/r16s-app/download_databases.sh"
echo ""
echo "  2. Start the application:"
echo "     bash ~/r16s-app/start.sh"
echo "     OR:  ngamp-start   (after opening a new terminal)"
echo ""
echo "  3. Open in browser:"
echo "     http://localhost:5173"
echo ""
echo "  4. Enter your license key when prompted."
echo "     (7-day trial keys can be generated with tools/keygen.py)"
echo ""
