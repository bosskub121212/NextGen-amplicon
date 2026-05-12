#!/usr/bin/env bash
# =============================================================
#  NextGen-Amplicon — Full Machine Installer
#  curl -fsSL https://raw.githubusercontent.com/bosskub121212/NextGen-amplicon/main/install.sh | bash
# =============================================================
set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YELLOW}  ⚠${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

REPO_URL="https://github.com/bosskub121212/NextGen-amplicon.git"
INSTALL_DIR="$HOME/r16s-app"
VENV_DIR="$INSTALL_DIR/venv"
NODE_VERSION="20"
R_VERSION_MAJOR="4"

echo ""
echo -e "${BOLD}${CYAN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║       NextGen-Amplicon — Full Machine Installer          ║"
echo "  ║       16S / ITS / COX1 / 18S-nema / PacBio              ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

step "Pre-flight checks"
[[ "$EUID" -eq 0 ]] && die "Do not run as root."
ok "Running as user: $(whoami)"

step "System packages (apt)"
sudo apt-get update -qq
sudo apt-get install -y \
  git curl wget gnupg2 ca-certificates software-properties-common \
  build-essential gfortran cmake \
  libssl-dev libcurl4-openssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libgsl-dev libglpk-dev libgmp3-dev zlib1g-dev \
  libopenblas-dev liblapack-dev \
  libuv1-dev libcairo2-dev \
  python3 python3-pip python3-venv python3-dev \
  nodejs npm procps -qq
ok "System packages installed"

step "Node.js $NODE_VERSION"
export NVM_DIR="$HOME/.nvm"
if [[ ! -d "$NVM_DIR" ]]; then
  curl -fsSL https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash
fi
source "$NVM_DIR/nvm.sh" 2>/dev/null || true
if ! nvm ls "$NODE_VERSION" &>/dev/null; then
  nvm install "$NODE_VERSION"
  nvm alias default "$NODE_VERSION"
fi
nvm use "$NODE_VERSION" >/dev/null 2>&1
ok "Node.js $(node --version) / npm $(npm --version)"

step "R $R_VERSION_MAJOR.x"
R_INSTALLED=false
if command -v R &>/dev/null; then
  R_VER=$(R --version | head -1 | grep -oP '\d+\.\d+' | head -1)
  [[ "${R_VER%%.*}" -ge "$R_VERSION_MAJOR" ]] && R_INSTALLED=true && ok "R $R_VER already installed"
fi
if [[ "$R_INSTALLED" == "false" ]]; then
  sudo apt-get install -y --no-install-recommends dirmngr -qq
  UBUNTU_CODENAME=$(lsb_release -cs 2>/dev/null || echo "jammy")
  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${UBUNTU_CODENAME}-cran40/" \
    | sudo tee /etc/apt/sources.list.d/cran_r.list >/dev/null
  sudo apt-get update -qq
  sudo apt-get install -y r-base r-base-dev -qq
  ok "R installed"
fi

step "Repository"
if [[ -d "$INSTALL_DIR/.git" ]]; then
  info "Existing repo found — pulling latest..."
  cd "$INSTALL_DIR"
  git pull origin main --ff-only || warn "git pull failed — using local version"
  ok "Repository up to date"
else
  git clone "$REPO_URL" "$INSTALL_DIR"
  ok "Repository cloned"
fi
cd "$INSTALL_DIR"

step "Python environment (FastAPI backend)"
[[ ! -d "$VENV_DIR" ]] && python3 -m venv "$VENV_DIR" && ok "venv created"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip -q
pip install fastapi "uvicorn[standard]" pydantic python-multipart aiofiles requests psutil pandas numpy -q
ok "Python packages installed"
deactivate

step "Cutadapt (dedicated venv)"
CUTADAPT_VENV="$HOME/cutadapt-env"
if [[ ! -d "$CUTADAPT_VENV" ]]; then
  python3 -m venv "$CUTADAPT_VENV"
  ok "cutadapt-env created"
fi
source "$CUTADAPT_VENV/bin/activate"
pip install --upgrade pip -q
pip install cutadapt -q
CUTADAPT_VERSION=$(cutadapt --version)
ok "cutadapt $CUTADAPT_VERSION installed"
deactivate

# Symlink so backend can call 'cutadapt' directly without activating venv
mkdir -p "$HOME/.local/bin"
ln -sf "$CUTADAPT_VENV/bin/cutadapt" "$HOME/.local/bin/cutadapt"
# Ensure ~/.local/bin is in PATH
if ! grep -q 'HOME/.local/bin' ~/.bashrc 2>/dev/null; then
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
fi
export PATH="$HOME/.local/bin:$PATH"
ok "cutadapt symlinked → $HOME/.local/bin/cutadapt"

step "R packages (15-40 min first time)"
Rscript "$INSTALL_DIR/backend/r_scripts/install_packages.R"
ok "R packages done"

step "Frontend (React build)"
source "$NVM_DIR/nvm.sh" 2>/dev/null || true
nvm use "$NODE_VERSION" >/dev/null 2>&1
cd "$INSTALL_DIR/frontend"
npm install --legacy-peer-deps -q
npm run build
ok "Frontend built"
cd "$INSTALL_DIR"

step "Directory structure"
mkdir -p "$INSTALL_DIR"/{backend/uploads,backend/results,databases/{SILVA,UNITE,MIDORI2,NemaBase},logs}
DB_PATHS="$INSTALL_DIR/backend/databases/db_paths.json"
mkdir -p "$(dirname "$DB_PATHS")"
[[ ! -f "$DB_PATHS" ]] && echo '{"silva":"","unite":"","midori2_co1":"","midori2_co1_sp":"","nemabase":"","emu_silva":""}' > "$DB_PATHS"
# Add emu_silva key to existing db_paths.json if missing
python3 - <<'PYEOF'
import json, pathlib
p = pathlib.Path("backend/databases/db_paths.json")
if p.exists():
    d = json.loads(p.read_text())
    if "emu_silva" not in d:
        d["emu_silva"] = ""
        p.write_text(json.dumps(d, indent=2))
        print("  added emu_silva key to db_paths.json")
PYEOF
ok "Directories ready"

step "Optional tools (KronaTools + PICRUSt2)"

# KronaTools (Krona interactive charts)
if ! command -v ktImportText &>/dev/null; then
  info "Installing KronaTools..."
  sudo apt-get install -y krona -qq 2>/dev/null && \
    ok "KronaTools installed" || \
    warn "KronaTools not in apt — install manually: conda install -c bioconda krona"
else
  ok "KronaTools already installed: $(ktImportText 2>&1 | head -1 || true)"
fi

# PICRUSt2 — auto-install if conda available (KEGG, MetaCyc, COG predictions)
if command -v conda &>/dev/null; then
  if conda env list | grep -q "^picrust2"; then
    ok "PICRUSt2 conda env already exists"
  else
    echo ""
    echo -e "${BOLD}${CYAN}  PICRUSt2 enables KEGG, MetaCyc, and COG functional prediction (~3 GB, 10-20 min).${NC}"
    read -r -p "  Install PICRUSt2 now? [Y/n] " _PICRUST_ANS
    _PICRUST_ANS="${_PICRUST_ANS:-Y}"
    if [[ "$_PICRUST_ANS" =~ ^[Yy] ]]; then
      info "Installing PICRUSt2 via conda..."
      # Use mamba if available (faster), fall back to conda
      _CONDA_CMD="conda"
      command -v mamba &>/dev/null && _CONDA_CMD="mamba"
      $_CONDA_CMD create -n picrust2 -c bioconda -c conda-forge picrust2 -y 2>&1 | \
        grep -E "^(Preparing|Executing|Installing|done|ERROR)" || true
      if conda env list | grep -q "^picrust2"; then
        ok "PICRUSt2 installed — KEGG/MetaCyc/COG prediction enabled"
      else
        warn "PICRUSt2 install failed — run manually later:"
        warn "  conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y"
      fi
    else
      info "PICRUSt2 skipped — re-run install.sh to add it later"
    fi
  fi
else
  warn "conda not found — PICRUSt2 skipped"
  warn "  Install miniconda3 from https://docs.conda.io/en/latest/miniconda.html"
  warn "  Then: conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y"
fi

# FAPROTAX (Python, ~50MB)
FAPROTAX_DIR="$HOME/faprotax"
if [[ ! -d "$FAPROTAX_DIR" ]]; then
  info "Installing FAPROTAX..."
  python3 -m pip install faprotax -q 2>/dev/null && \
    ok "FAPROTAX installed" || {
      # Fallback: download from official site
      warn "pip install failed — trying direct download..."
      curl -fsSL "https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/FAPROTAX_1.2.4.zip" \
        -o /tmp/faprotax.zip 2>/dev/null && \
        unzip -q /tmp/faprotax.zip -d "$HOME/" && \
        mv "$HOME/FAPROTAX_1.2.4" "$FAPROTAX_DIR" && \
        ok "FAPROTAX downloaded to $FAPROTAX_DIR" || \
        warn "FAPROTAX install failed — install manually from http://www.loucalab.com/archive/FAPROTAX/"
    }
else
  ok "FAPROTAX already at $FAPROTAX_DIR"
fi

step "Optional tools (Emu — ONT 16S)"

# Emu — ONT full-length / sub-region 16S pipeline
if command -v conda &>/dev/null; then
  if conda run -n emu emu --version &>/dev/null 2>&1; then
    ok "Emu already installed in conda env 'emu'"
  else
    echo ""
    echo -e "${BOLD}${CYAN}  Emu enables ONT 16S analysis (V7-V8, V1-V9). Requires ~500 MB + ~4 GB database.${NC}"
    read -r -p "  Install Emu now? [Y/n] " _EMU_ANS
    _EMU_ANS="${_EMU_ANS:-Y}"
    if [[ "$_EMU_ANS" =~ ^[Yy] ]]; then
      _CONDA_CMD="conda"; command -v mamba &>/dev/null && _CONDA_CMD="mamba"
      info "Creating conda env 'emu'..."
      $_CONDA_CMD create -n emu -c bioconda -c conda-forge emu -y 2>&1 | \
        grep -E "^(Preparing|Executing|Installing|done|ERROR)" || true
      if conda run -n emu emu --version &>/dev/null 2>&1; then
        ok "Emu installed"
        echo ""
        echo -e "${BOLD}${CYAN}  Download Emu SILVA database now? (~4 GB)${NC}"
        read -r -p "  Download Emu database? [Y/n] " _EMUDB_ANS
        _EMUDB_ANS="${_EMUDB_ANS:-Y}"
        if [[ "$_EMUDB_ANS" =~ ^[Yy] ]]; then
          EMU_DB_DIR="$HOME/r16s-app/backend/databases/emu_silva"
          mkdir -p "$EMU_DB_DIR"
          info "Downloading Emu SILVA database to $EMU_DB_DIR ..."
          conda run -n emu emu download-db silva --db-dir "$EMU_DB_DIR" && \
            ok "Emu database ready at $EMU_DB_DIR" || \
            warn "Emu database download failed — run manually: emu download-db silva --db-dir ~/r16s-app/backend/databases/emu_silva"
        else
          info "Emu database skipped — download later: emu download-db silva --db-dir ~/r16s-app/backend/databases/emu_silva"
        fi
      else
        warn "Emu install failed — install manually: conda create -n emu -c bioconda emu"
      fi
    else
      info "Emu skipped — install later for ONT 16S support"
    fi
  fi
else
  warn "conda not found — Emu skipped (requires conda)"
  warn "  Install miniconda3, then: conda create -n emu -c bioconda emu"
fi

step "Helper scripts"
cat > "$INSTALL_DIR/start_backend.sh" <<'STARTSH'
#!/usr/bin/env bash
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
STARTSH
chmod +x "$INSTALL_DIR/start_backend.sh"
ok "start_backend.sh ready"

ALIAS_BLOCK='
# NextGen-Amplicon shortcuts
alias ngamp-start="bash ~/r16s-app/start_backend.sh"
alias ngamp-log="tail -f ~/r16s-app/logs/backend.log"'
for RC in ~/.bashrc ~/.zshrc; do
  [[ -f "$RC" ]] && ! grep -q "ngamp-start" "$RC" 2>/dev/null && \
    echo "$ALIAS_BLOCK" >> "$RC" && ok "Aliases added to $RC"
done

echo ""
echo -e "${BOLD}${GREEN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║   ✅  Installation complete!                             ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"
echo "  Start: bash ~/r16s-app/start_backend.sh"
echo "  Update later: bash ~/r16s-app/update.sh"
echo ""
