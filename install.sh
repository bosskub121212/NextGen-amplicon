#!/usr/bin/env bash
# =============================================================
#  NextGen-Amplicon — Full Machine Installer  v2.1.1
#  curl -fsSL https://raw.githubusercontent.com/bosskub121212/NextGen-amplicon/main/install.sh | bash
#
#  Installs:
#    • System packages, Node.js, R
#    • Python backend (FastAPI) + Cutadapt venv
#    • R packages (DADA2, phyloseq, ggplot2, vegan, etc.)
#    • QIIME2 conda environment
#    • Emu conda environment (ONT 16S)
#    • PICRUSt2 conda environment (optional)
#    • KronaTools, FAPROTAX
#    • Database downloads (SILVA, UNITE, MIDORI2, NemaBase, Emu)
#    • Emu database builder (build_emu_db.py)
#    • Frontend (React build)
#    • Helper scripts
# =============================================================
set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YELLOW}  ⚠${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }
ask()  { read -r -p "  $1 [Y/n] " _ANS; _ANS="${_ANS:-Y}"; [[ "$_ANS" =~ ^[Yy] ]]; }

REPO_URL="https://github.com/bosskub121212/NextGen-amplicon.git"
INSTALL_DIR="$HOME/r16s-app"
VENV_DIR="$INSTALL_DIR/venv"
NODE_VERSION="20"
R_VERSION_MAJOR="4"
QIIME2_ENV="qiime2-amplicon-2024.10"
QIIME2_CHANNEL="https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml"

echo ""
echo -e "${BOLD}${CYAN}"
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║       NextGen-Amplicon — Full Machine Installer  v2.1.1     ║"
echo "  ║       16S / ITS / COX1 / 18S / PacBio / ONT-16S            ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ── Pre-flight ───────────────────────────────────────────────────────────────
step "Pre-flight checks"
[[ "$EUID" -eq 0 ]] && die "Do not run as root."
ok "Running as user: $(whoami)"

# ── System packages ──────────────────────────────────────────────────────────
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
  nodejs npm procps pigz -qq
ok "System packages installed"

# ── Node.js ──────────────────────────────────────────────────────────────────
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

# ── R ────────────────────────────────────────────────────────────────────────
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

# ── Repository ───────────────────────────────────────────────────────────────
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

# ── Python backend venv ──────────────────────────────────────────────────────
step "Python environment (FastAPI backend)"
[[ ! -d "$VENV_DIR" ]] && python3 -m venv "$VENV_DIR" && ok "venv created"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip -q
pip install \
  fastapi "uvicorn[standard]" pydantic python-multipart aiofiles \
  requests psutil pandas numpy scipy biopython -q
ok "Python backend packages installed"
deactivate

# ── Cutadapt venv ────────────────────────────────────────────────────────────
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

mkdir -p "$HOME/.local/bin"
ln -sf "$CUTADAPT_VENV/bin/cutadapt" "$HOME/.local/bin/cutadapt"
if ! grep -q 'HOME/.local/bin' ~/.bashrc 2>/dev/null; then
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
fi
export PATH="$HOME/.local/bin:$PATH"
ok "cutadapt symlinked → $HOME/.local/bin/cutadapt"

# ── R packages ───────────────────────────────────────────────────────────────
step "R packages (15-40 min first time)"
Rscript "$INSTALL_DIR/backend/r_scripts/install_packages.R"
ok "R packages done"

# ── Directory structure ───────────────────────────────────────────────────────
step "Directory structure"
mkdir -p "$INSTALL_DIR"/{backend/uploads,backend/results,backend/classifiers,logs}
mkdir -p "$INSTALL_DIR/backend/databases"/{SILVA,UNITE,MIDORI2,NemaBase,emu_silva,qiime2_classifiers}
DB_PATHS="$INSTALL_DIR/backend/databases/db_paths.json"
if [[ ! -f "$DB_PATHS" ]]; then
  cat > "$DB_PATHS" <<'JSONEOF'
{
  "SILVA_16S": "",
  "SILVA_16S_sp": "",
  "PR2_18S": "",
  "NemaBase_18S": "",
  "UNITE_ITS1": "",
  "MIDORI2_COX1": "",
  "emu_silva": ""
}
JSONEOF
fi
# Ensure all expected keys exist
python3 - <<'PYEOF'
import json, pathlib
p = pathlib.Path("backend/databases/db_paths.json")
if p.exists():
    d = json.loads(p.read_text())
    keys = ["SILVA_16S","SILVA_16S_sp","PR2_18S","NemaBase_18S","UNITE_ITS1","MIDORI2_COX1","emu_silva"]
    changed = False
    for k in keys:
        if k not in d:
            d[k] = ""; changed = True
    if changed:
        p.write_text(json.dumps(d, indent=2))
        print("  ✓  db_paths.json updated with missing keys")
PYEOF
ok "Directories ready"

# ── Conda (Miniconda) ────────────────────────────────────────────────────────
step "Conda (Miniconda)"
if ! command -v conda &>/dev/null; then
  MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
  info "Conda not found — installing Miniconda3..."
  curl -fsSL "$MINICONDA_URL" -o /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
  rm /tmp/miniconda.sh
  export PATH="$HOME/miniconda3/bin:$PATH"
  conda init bash 2>/dev/null || true
  source "$HOME/.bashrc" 2>/dev/null || true
  ok "Miniconda3 installed at $HOME/miniconda3"
else
  ok "Conda already available: $(conda --version)"
fi
# Ensure conda is in PATH for this session
CONDA_PREFIX=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
export PATH="$CONDA_PREFIX/bin:$PATH"
# Use mamba if available (much faster)
_CONDA_CMD="conda"
if command -v mamba &>/dev/null; then
  _CONDA_CMD="mamba"
  ok "Using mamba for faster installs"
elif conda install -n base -c conda-forge mamba -y -q 2>/dev/null; then
  _CONDA_CMD="mamba"
  ok "mamba installed for faster conda operations"
fi

# ── QIIME2 conda env ─────────────────────────────────────────────────────────
step "QIIME2 (conda env: $QIIME2_ENV)"
if conda env list 2>/dev/null | grep -q "^${QIIME2_ENV}[[:space:]]"; then
  ok "QIIME2 env '$QIIME2_ENV' already exists"
else
  echo ""
  echo -e "${BOLD}${CYAN}  QIIME2 is required for 16S / ITS / COX1 / 18S pipelines (~5 GB, 20-40 min).${NC}"
  if ask "Install QIIME2 now?"; then
    info "Downloading QIIME2 environment definition..."
    QIIME2_YML="/tmp/qiime2-amplicon.yml"
    curl -fsSL "$QIIME2_CHANNEL" -o "$QIIME2_YML" || \
      wget -q "$QIIME2_CHANNEL" -O "$QIIME2_YML"
    info "Creating QIIME2 conda env (this may take 20–40 min)..."
    $_CONDA_CMD env create -n "$QIIME2_ENV" --file "$QIIME2_YML" 2>&1 | \
      grep -E "^(Preparing|Executing|Installing|done|ERROR|error)" || true
    rm -f "$QIIME2_YML"
    if conda env list | grep -q "^${QIIME2_ENV}"; then
      ok "QIIME2 env '$QIIME2_ENV' installed"
      # Install additional plugins
      conda run -n "$QIIME2_ENV" pip install q2-fondue 2>/dev/null && \
        ok "q2-fondue installed" || true
    else
      warn "QIIME2 install failed — try manually:"
      warn "  wget $QIIME2_CHANNEL -O qiime2.yml && conda env create -n $QIIME2_ENV --file qiime2.yml"
    fi
  else
    info "QIIME2 skipped — re-run install.sh or:"
    info "  wget $QIIME2_CHANNEL -O qiime2.yml && conda env create -n $QIIME2_ENV --file qiime2.yml"
  fi
fi

# ── Emu (ONT 16S) ────────────────────────────────────────────────────────────
step "Emu — ONT 16S (conda env: emu)"
if conda run -n emu emu --version &>/dev/null 2>&1; then
  ok "Emu already installed in conda env 'emu'"
else
  echo ""
  echo -e "${BOLD}${CYAN}  Emu enables ONT 16S analysis (V7-V8, V1-V9). Requires ~500 MB conda env.${NC}"
  if ask "Install Emu now?"; then
    info "Creating conda env 'emu'..."
    $_CONDA_CMD create -n emu -c bioconda -c conda-forge emu -y 2>&1 | \
      grep -E "^(Preparing|Executing|Installing|done|ERROR)" || true
    if conda run -n emu emu --version &>/dev/null 2>&1; then
      ok "Emu installed ($(conda run -n emu emu --version 2>/dev/null || echo 'ok'))"
    else
      warn "Emu install failed — try manually: conda create -n emu -c bioconda emu -y"
    fi
  else
    info "Emu skipped — re-run install.sh or: conda create -n emu -c bioconda emu -y"
  fi
fi

# ── PICRUSt2 (optional) ──────────────────────────────────────────────────────
step "PICRUSt2 — Functional prediction (optional)"
if conda env list 2>/dev/null | grep -q "^picrust2[[:space:]]"; then
  ok "PICRUSt2 conda env already exists"
else
  echo ""
  echo -e "${BOLD}${CYAN}  PICRUSt2 enables KEGG, MetaCyc, COG predictions (~3 GB, 10-20 min).${NC}"
  if ask "Install PICRUSt2 now?"; then
    info "Creating conda env 'picrust2'..."
    $_CONDA_CMD create -n picrust2 -c bioconda -c conda-forge picrust2 -y 2>&1 | \
      grep -E "^(Preparing|Executing|Installing|done|ERROR)" || true
    conda env list | grep -q "^picrust2" && \
      ok "PICRUSt2 installed" || \
      warn "PICRUSt2 failed — try: conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y"
  else
    info "PICRUSt2 skipped — re-run install.sh to add it later"
  fi
fi

# ── KronaTools ───────────────────────────────────────────────────────────────
step "KronaTools"
if ! command -v ktImportText &>/dev/null; then
  sudo apt-get install -y krona -qq 2>/dev/null && \
    ok "KronaTools installed" || {
    $_CONDA_CMD install -n base -c bioconda krona -y -q 2>/dev/null && \
      ok "KronaTools installed via conda" || \
      warn "KronaTools not installed — install manually: conda install -c bioconda krona"
  }
else
  ok "KronaTools already installed"
fi

# ── FAPROTAX ─────────────────────────────────────────────────────────────────
step "FAPROTAX"
if python3 -c "import faprotax" 2>/dev/null; then
  ok "FAPROTAX already installed"
elif python3 -m pip install faprotax -q 2>/dev/null; then
  ok "FAPROTAX installed via pip"
else
  warn "FAPROTAX pip install failed — install manually: pip install faprotax"
fi

# ═══════════════════════════════════════════════════════════════════════════════
#  DATABASE DOWNLOADS
# ═══════════════════════════════════════════════════════════════════════════════
step "Database downloads"
DB_DIR="$INSTALL_DIR/backend/databases"

echo ""
echo -e "${BOLD}${CYAN}  Which databases do you want to download?${NC}"
echo "  (Press Enter = Yes, type n = Skip)"
echo ""

# ── 1. SILVA 16S — DADA2 trainset (for R pipeline) ───────────────────────────
SILVA_DADA2_DIR="$DB_DIR/SILVA"
SILVA_DADA2_FILE="$SILVA_DADA2_DIR/silva_nr99_v138.1_train_set.fa.gz"
SILVA_DADA2_SP_FILE="$SILVA_DADA2_DIR/silva_nr99_v138.1_wSpecies_train_set.fa.gz"

if [[ ! -f "$SILVA_DADA2_FILE" ]]; then
  if ask "Download SILVA 138.1 — DADA2 trainset (~130 MB, genus+species)?"; then
    mkdir -p "$SILVA_DADA2_DIR"
    info "Downloading SILVA 138.1 DADA2 genus-level trainset..."
    wget -q --show-progress \
      "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" \
      -O "$SILVA_DADA2_FILE" && ok "SILVA 16S genus trainset → $SILVA_DADA2_FILE"
    info "Downloading SILVA 138.1 DADA2 species-level trainset..."
    wget -q --show-progress \
      "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz" \
      -O "$SILVA_DADA2_SP_FILE" && ok "SILVA 16S species trainset → $SILVA_DADA2_SP_FILE"
    # Update db_paths.json
    python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
d = json.loads(p.read_text())
d["SILVA_16S"]    = "$SILVA_DADA2_FILE"
d["SILVA_16S_sp"] = "$SILVA_DADA2_SP_FILE"
p.write_text(json.dumps(d, indent=2))
print("  ✓  db_paths.json updated: SILVA_16S, SILVA_16S_sp")
PYEOF
  fi
else
  ok "SILVA 16S DADA2 trainset already exists"
fi

# ── 2. SILVA QIIME2 sequences + taxonomy (for QIIME2 + custom classifier training) ──
SILVA_Q2_DIR="$DB_DIR/qiime2_classifiers"
SILVA_Q2_SEQS="$SILVA_Q2_DIR/silva-138-99-seqs.qza"
SILVA_Q2_TAX="$SILVA_Q2_DIR/silva-138-99-tax.qza"

if [[ ! -f "$SILVA_Q2_SEQS" ]]; then
  if ask "Download SILVA 138 QIIME2 sequences (~1 GB, needed for custom classifier training)?"; then
    mkdir -p "$SILVA_Q2_DIR"
    info "Downloading SILVA 138 QIIME2 reference sequences..."
    wget -q --show-progress \
      "https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza" \
      -O "$SILVA_Q2_SEQS" && ok "SILVA QIIME2 seqs → $SILVA_Q2_SEQS"
    info "Downloading SILVA 138 QIIME2 taxonomy..."
    wget -q --show-progress \
      "https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza" \
      -O "$SILVA_Q2_TAX" && ok "SILVA QIIME2 taxonomy → $SILVA_Q2_TAX"
  fi
else
  ok "SILVA QIIME2 reference sequences already exist"
fi

# ── 3. SILVA pre-trained QIIME2 classifiers (V3-V4, V4, full length) ─────────
if ask "Download pre-trained QIIME2 classifiers for SILVA (V3-V4, V4, full-length)?"; then
  mkdir -p "$SILVA_Q2_DIR"
  declare -A PREBUILT=(
    ["classifier_SILVA_V3V4.qza"]="https://data.qiime2.org/2024.10/common/silva-138-99-seqs-341-806.qza"
    ["classifier_SILVA_V4.qza"]="https://data.qiime2.org/2024.10/common/silva-138-99-seqs-515-806.qza"
    ["classifier_SILVA_full.qza"]="https://data.qiime2.org/2024.10/common/silva-138-99-nb-classifier.qza"
  )
  # Note: the V3-V4 and V4 above are seqs files; real pre-trained NB classifiers:
  declare -A PREBUILT_CLF=(
    ["classifier_SILVA_V3V4_NB.qza"]="https://data.qiime2.org/classifiers/silva-138-99-seqs-341-806.qza"
    ["classifier_SILVA_full_NB.qza"]="https://data.qiime2.org/2024.10/common/silva-138-99-nb-classifier.qza"
  )
  for FNAME in "classifier_SILVA_full_NB.qza"; do
    OUT_FILE="$INSTALL_DIR/backend/classifiers/$FNAME"
    if [[ ! -f "$OUT_FILE" ]]; then
      URL="${PREBUILT_CLF[$FNAME]}"
      info "Downloading $FNAME (~800 MB)..."
      wget -q --show-progress "$URL" -O "$OUT_FILE" 2>/dev/null && \
        ok "$FNAME ready" || warn "Download failed: $FNAME — get from https://docs.qiime2.org/2024.10/data-resources/"
    else
      ok "$FNAME already exists"
    fi
  done
fi

# ── 4. UNITE ITS ─────────────────────────────────────────────────────────────
UNITE_DIR="$DB_DIR/UNITE"
UNITE_FILE="$UNITE_DIR/unite_ITS_v10.fa.gz"
if [[ ! -f "$UNITE_FILE" ]]; then
  if ask "Download UNITE v10 — ITS Fungi (~150 MB)?"; then
    mkdir -p "$UNITE_DIR"
    info "Downloading UNITE v10 for DADA2 (dynamic representatives)..."
    wget -q --show-progress \
      "https://files.plutof.ut.ee/public/orig/98/AJ/98AJ45OYPC4M7LO975UYQZ3LDSKLJLID7VXHCOHCFRQCNJYIWFQ.gz" \
      -O "$UNITE_FILE" 2>/dev/null && ok "UNITE ITS → $UNITE_FILE" || {
      warn "UNITE download failed — download manually from https://unite.ut.ee/repository.php"
      warn "  → Save as: $UNITE_FILE"
    }
    python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
d = json.loads(p.read_text())
d["UNITE_ITS1"] = "$UNITE_FILE"
p.write_text(json.dumps(d, indent=2))
PYEOF
  fi
else
  ok "UNITE ITS already exists"
fi

# ── 5. MIDORI2 COX1 ──────────────────────────────────────────────────────────
MIDORI_DIR="$DB_DIR/MIDORI2"
MIDORI_FILE="$MIDORI_DIR/MIDORI2_COX1_dada2.fa.gz"
if [[ ! -f "$MIDORI_FILE" ]]; then
  if ask "Download MIDORI2 COX1 — Animal barcoding (~300 MB)?"; then
    mkdir -p "$MIDORI_DIR"
    info "Downloading MIDORI2 COX1 longest DADA2 release..."
    # MIDORI2 GB256 COX1 — DADA2 format
    MIDORI_URL="https://www.reference-midori.info/download/Databases/GenBank256_2023-12-15/DADA2/uniq/MIDORI2_UNIQ_NUC_GB256_CO1_DADA2.fasta.gz"
    wget -q --show-progress "$MIDORI_URL" -O "$MIDORI_FILE" 2>/dev/null && \
      ok "MIDORI2 COX1 → $MIDORI_FILE" || {
      warn "MIDORI2 download failed — download manually from https://www.reference-midori.info/"
      warn "  Format: DADA2, Gene: CO1, → Save as: $MIDORI_FILE"
    }
    python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
d = json.loads(p.read_text())
d["MIDORI2_COX1"] = "$MIDORI_FILE"
p.write_text(json.dumps(d, indent=2))
PYEOF
  fi
else
  ok "MIDORI2 COX1 already exists"
fi

# ── 6. NemaBase 18S ──────────────────────────────────────────────────────────
NEMA_DIR="$DB_DIR/NemaBase"
NEMA_FILE="$NEMA_DIR/NemaBase_18S.fa.gz"
if [[ ! -f "$NEMA_FILE" ]]; then
  if ask "Download NemaBase — Nematode 18S (~30 MB)?"; then
    mkdir -p "$NEMA_DIR"
    info "Downloading NemaBase 18S for DADA2..."
    # PR2 v5 as fallback for 18S nematode (NemaBase is typically local install)
    PR2_URL="https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz"
    wget -q --show-progress "$PR2_URL" -O "$NEMA_FILE" 2>/dev/null && \
      ok "PR2 v5 SSU (18S) → $NEMA_FILE (used as NemaBase fallback)" || {
      warn "Download failed — place your NemaBase FASTA at: $NEMA_FILE"
    }
    python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
d = json.loads(p.read_text())
d["NemaBase_18S"] = "$NEMA_FILE"
d["PR2_18S"]      = "$NEMA_FILE"
p.write_text(json.dumps(d, indent=2))
PYEOF
  fi
else
  ok "NemaBase/PR2 18S already exists"
fi

# ── 7. Emu SILVA database ─────────────────────────────────────────────────────
EMU_DB_DIR="$DB_DIR/emu_silva"
if [[ -f "$EMU_DB_DIR/species_list.tsv" ]]; then
  ok "Emu SILVA database already built at $EMU_DB_DIR"
else
  if command -v conda &>/dev/null && conda run -n emu emu --version &>/dev/null 2>&1; then
    echo ""
    echo -e "${BOLD}${CYAN}  Emu SILVA database for ONT 16S analysis (~4 GB download, ~30 min build).${NC}"
    if ask "Build Emu SILVA database now?"; then
      mkdir -p "$EMU_DB_DIR"
      if [[ -f "$INSTALL_DIR/backend/python_scripts/build_emu_db.py" ]]; then
        info "Building Emu DB via build_emu_db.py..."
        source "$VENV_DIR/bin/activate"
        python3 "$INSTALL_DIR/backend/python_scripts/build_emu_db.py" \
          --output-dir "$EMU_DB_DIR" \
          --silva-db "$SILVA_DADA2_FILE" 2>/dev/null || {
          # Fallback: use Emu's own download
          info "Trying emu download-db silva..."
          cd "$EMU_DB_DIR"
          conda run -n emu emu download-db silva 2>/dev/null && \
            ok "Emu SILVA database ready" || \
            warn "Emu DB download failed — run later: bash $INSTALL_DIR/setup_emu_db.sh"
          cd "$INSTALL_DIR"
        }
        deactivate
      else
        # Direct Emu download
        cd "$EMU_DB_DIR"
        conda run -n emu emu download-db silva && \
          ok "Emu SILVA database ready at $EMU_DB_DIR" || \
          warn "Emu DB failed — run later: bash $INSTALL_DIR/setup_emu_db.sh"
        cd "$INSTALL_DIR"
      fi
      # Update db_paths.json
      python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
d = json.loads(p.read_text())
d["emu_silva"] = "$EMU_DB_DIR"
p.write_text(json.dumps(d, indent=2))
print("  ✓  db_paths.json updated: emu_silva")
PYEOF
    else
      info "Emu database skipped — run later:"
      info "  bash $INSTALL_DIR/setup_emu_db.sh"
    fi
  else
    warn "Emu not installed — skipping Emu database setup"
    warn "  Install Emu first, then: bash $INSTALL_DIR/setup_emu_db.sh"
  fi
fi

# ── Frontend ──────────────────────────────────────────────────────────────────
step "Frontend (React build)"
source "$NVM_DIR/nvm.sh" 2>/dev/null || true
nvm use "$NODE_VERSION" >/dev/null 2>&1
cd "$INSTALL_DIR/frontend"
npm install --legacy-peer-deps -q
npm run build
ok "Frontend built"
cd "$INSTALL_DIR"

# ── Emu setup helper script ───────────────────────────────────────────────────
step "Helper scripts"

cat > "$INSTALL_DIR/setup_emu_db.sh" <<'EMUSETUP'
#!/usr/bin/env bash
# Download and build Emu SILVA database for ONT 16S
# Run this script: bash ~/r16s-app/setup_emu_db.sh
set -euo pipefail
INSTALL_DIR="$HOME/r16s-app"
EMU_DB_DIR="$INSTALL_DIR/backend/databases/emu_silva"
mkdir -p "$EMU_DB_DIR"

echo "=== Setting up Emu SILVA database ==="
echo "  Output: $EMU_DB_DIR"
echo ""

if ! conda run -n emu emu --version &>/dev/null 2>&1; then
  echo "ERROR: Emu not found. Install first:"
  echo "  conda create -n emu -c bioconda -c conda-forge emu -y"
  exit 1
fi

# Check if build_emu_db.py exists (custom builder from DADA2 SILVA trainset)
BUILDER="$INSTALL_DIR/backend/python_scripts/build_emu_db.py"
SILVA_FILE="$INSTALL_DIR/backend/databases/SILVA/silva_nr99_v138.1_train_set.fa.gz"

if [[ -f "$BUILDER" && -f "$SILVA_FILE" ]]; then
  echo "Building Emu DB from SILVA DADA2 trainset..."
  source "$INSTALL_DIR/venv/bin/activate"
  python3 "$BUILDER" --output-dir "$EMU_DB_DIR" --silva-db "$SILVA_FILE"
  deactivate
else
  echo "Downloading Emu SILVA database from NCBI (~4 GB)..."
  cd "$EMU_DB_DIR"
  conda run -n emu emu download-db silva
fi

# Update db_paths.json
python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$INSTALL_DIR/backend/databases/db_paths.json")
if p.exists():
    d = json.loads(p.read_text())
    d["emu_silva"] = "$EMU_DB_DIR"
    p.write_text(json.dumps(d, indent=2))
    print("  ✓  db_paths.json updated: emu_silva =", "$EMU_DB_DIR")
PYEOF

echo ""
echo "=== Emu SILVA database ready ==="
echo "  Restart backend to apply: pkill -f 'uvicorn main:app' && sleep 2 && bash ~/r16s-app/start_backend.sh &"
EMUSETUP
chmod +x "$INSTALL_DIR/setup_emu_db.sh"
ok "setup_emu_db.sh ready"

# ── Database download helper script ──────────────────────────────────────────
cat > "$INSTALL_DIR/download_databases.sh" <<'DLSCRIPT'
#!/usr/bin/env bash
# Download individual databases after initial install
# Usage: bash ~/r16s-app/download_databases.sh [silva|unite|midori|nema|emu|all]
set -euo pipefail
INSTALL_DIR="$HOME/r16s-app"
DB_DIR="$INSTALL_DIR/backend/databases"

update_db_path() {
  local KEY="$1" VAL="$2"
  python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_DIR/db_paths.json")
if p.exists():
    d = json.loads(p.read_text())
    d["$KEY"] = "$VAL"
    p.write_text(json.dumps(d, indent=2))
    print("  ✓  db_paths.json: $KEY =", "$VAL")
PYEOF
}

dl_silva() {
  echo "=== Downloading SILVA 138.1 DADA2 trainset ==="
  mkdir -p "$DB_DIR/SILVA"
  wget -q --show-progress "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" \
    -O "$DB_DIR/SILVA/silva_nr99_v138.1_train_set.fa.gz"
  wget -q --show-progress "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz" \
    -O "$DB_DIR/SILVA/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
  update_db_path "SILVA_16S"    "$DB_DIR/SILVA/silva_nr99_v138.1_train_set.fa.gz"
  update_db_path "SILVA_16S_sp" "$DB_DIR/SILVA/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
  echo "=== SILVA QIIME2 reference sequences (for custom classifier training) ==="
  mkdir -p "$DB_DIR/qiime2_classifiers"
  wget -q --show-progress "https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza" \
    -O "$DB_DIR/qiime2_classifiers/silva-138-99-seqs.qza"
  wget -q --show-progress "https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza" \
    -O "$DB_DIR/qiime2_classifiers/silva-138-99-tax.qza"
  echo "  ✓  SILVA done"
}

dl_unite() {
  echo "=== Downloading UNITE v10 ITS ==="
  mkdir -p "$DB_DIR/UNITE"
  wget -q --show-progress \
    "https://files.plutof.ut.ee/public/orig/98/AJ/98AJ45OYPC4M7LO975UYQZ3LDSKLJLID7VXHCOHCFRQCNJYIWFQ.gz" \
    -O "$DB_DIR/UNITE/unite_ITS_v10.fa.gz"
  update_db_path "UNITE_ITS1" "$DB_DIR/UNITE/unite_ITS_v10.fa.gz"
  echo "  ✓  UNITE done"
}

dl_midori() {
  echo "=== Downloading MIDORI2 COX1 ==="
  mkdir -p "$DB_DIR/MIDORI2"
  wget -q --show-progress \
    "https://www.reference-midori.info/download/Databases/GenBank256_2023-12-15/DADA2/uniq/MIDORI2_UNIQ_NUC_GB256_CO1_DADA2.fasta.gz" \
    -O "$DB_DIR/MIDORI2/MIDORI2_COX1_dada2.fa.gz"
  update_db_path "MIDORI2_COX1" "$DB_DIR/MIDORI2/MIDORI2_COX1_dada2.fa.gz"
  echo "  ✓  MIDORI2 done"
}

dl_nema() {
  echo "=== Downloading PR2 v5 (18S / NemaBase) ==="
  mkdir -p "$DB_DIR/NemaBase"
  wget -q --show-progress \
    "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz" \
    -O "$DB_DIR/NemaBase/NemaBase_18S.fa.gz"
  update_db_path "NemaBase_18S" "$DB_DIR/NemaBase/NemaBase_18S.fa.gz"
  update_db_path "PR2_18S"      "$DB_DIR/NemaBase/NemaBase_18S.fa.gz"
  echo "  ✓  NemaBase/PR2 done"
}

dl_emu() {
  echo "=== Building Emu SILVA database ==="
  bash "$INSTALL_DIR/setup_emu_db.sh"
}

TARGET="${1:-all}"
case "$TARGET" in
  silva)  dl_silva ;;
  unite)  dl_unite ;;
  midori) dl_midori ;;
  nema)   dl_nema ;;
  emu)    dl_emu ;;
  all)    dl_silva; dl_unite; dl_midori; dl_nema; dl_emu ;;
  *)      echo "Usage: $0 [silva|unite|midori|nema|emu|all]"; exit 1 ;;
esac
echo ""
echo "=== Database download complete ==="
echo "  Restart backend: pkill -f 'uvicorn main:app' && sleep 2 && bash ~/r16s-app/start_backend.sh &"
DLSCRIPT
chmod +x "$INSTALL_DIR/download_databases.sh"
ok "download_databases.sh ready"

# ── start_backend.sh ──────────────────────────────────────────────────────────
cat > "$INSTALL_DIR/start_backend.sh" <<'STARTSH'
#!/usr/bin/env bash
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
STARTSH
chmod +x "$INSTALL_DIR/start_backend.sh"
ok "start_backend.sh ready"

# ── Shell aliases ─────────────────────────────────────────────────────────────
ALIAS_BLOCK='
# NextGen-Amplicon shortcuts
alias ngamp-start="bash ~/r16s-app/start_backend.sh"
alias ngamp-log="tail -f ~/r16s-app/logs/backend.log"
alias ngamp-update="cd ~/r16s-app && git pull origin main && bash install.sh"
alias ngamp-db="bash ~/r16s-app/download_databases.sh"'
for RC in ~/.bashrc ~/.zshrc; do
  [[ -f "$RC" ]] && ! grep -q "ngamp-start" "$RC" 2>/dev/null && \
    echo "$ALIAS_BLOCK" >> "$RC" && ok "Aliases added to $RC"
done

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}"
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║   ✅  Installation complete!  v2.1.1                        ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"
echo "  Start backend:    bash ~/r16s-app/start_backend.sh"
echo "  Download DBs:     bash ~/r16s-app/download_databases.sh [silva|unite|midori|nema|emu|all]"
echo "  Setup Emu DB:     bash ~/r16s-app/setup_emu_db.sh"
echo "  Update app:       cd ~/r16s-app && git pull origin main && bash install.sh"
echo ""
echo "  Shortcuts (after reopening terminal):"
echo "    ngamp-start   — start backend"
echo "    ngamp-db      — download databases"
echo "    ngamp-update  — pull + reinstall"
echo ""
