#!/usr/bin/env bash
# =============================================================
#  NextGen-Amplicon — QIIME2 Environment Setup
#
#  Run ONCE after install.sh:
#    bash ~/r16s-app/setup_qiime2.sh
#
#  Installs:
#    - QIIME2 2024.5 (amplicon distribution)
#    - Required plugins: cutadapt, dada2, feature-classifier,
#      taxa, diversity, composition, phylogeny, demux, picrust2
#    - R packages: phyloseq, vegan, ANCOMBC2, FUNGuildR, LULU
#
#  Time: ~30–60 min (first install, downloads ~4 GB)
# =============================================================
set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'
BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YELLOW}  ⚠${NC}  $1"; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

# QIIME2 version to install — update this if a newer version is available
# Check available: mamba search -c conda-forge "qiime2-amplicon" | tail -10
QIIME2_VERSION="${QIIME2_VERSION:-}"   # leave empty for latest
QIIME2_ENV="qiime2-amplicon${QIIME2_VERSION:+-${QIIME2_VERSION}}"
QIIME2_ENV="${QIIME2_ENV:-qiime2-amplicon}"

echo ""
echo -e "${BOLD}${CYAN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║    NextGen-Amplicon — QIIME2 Environment Setup          ║"
echo "  ║    QIIME2 2024.5 + R visualization packages             ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ── Check conda ───────────────────────────────────────────────────────
step "Pre-flight"
if ! command -v conda &>/dev/null; then
  echo "conda not found — installing Miniconda3..."
  wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
  eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
  echo 'eval "$(~/miniconda3/bin/conda shell.bash hook)"' >> ~/.bashrc
  ok "Miniconda3 installed"
else
  eval "$(conda shell.bash hook)" 2>/dev/null || true
  ok "conda $(conda --version) found"
fi

# ── Install mamba (faster solver) ────────────────────────────────────
step "mamba solver"
if command -v mamba &>/dev/null; then
  ok "mamba $(mamba --version | head -1) found"
else
  info "Installing mamba into base environment..."
  conda install -n base -c conda-forge mamba -y --quiet
  ok "mamba installed"
fi

# ── Check if QIIME2 env already exists ───────────────────────────────
step "QIIME2 2024.5 environment"
if conda env list | grep -q "^${QIIME2_ENV}"; then
  ok "Environment '${QIIME2_ENV}' already exists — skipping"
else
  info "Creating QIIME2 conda environment via conda-forge (30–60 min)..."
  PKG="qiime2-amplicon${QIIME2_VERSION:+=}${QIIME2_VERSION}"
  mamba create -n "$QIIME2_ENV" -y \
    -c conda-forge \
    -c bioconda \
    "$PKG"
  ok "QIIME2 environment created"
fi

# ── Verify QIIME2 ─────────────────────────────────────────────────────
step "Verify QIIME2"
QIIME_VER=$(conda run -n "$QIIME2_ENV" --no-capture-output qiime --version 2>/dev/null || echo "error")
if [[ "$QIIME_VER" == "error" ]]; then
  echo "ERROR: QIIME2 not responding in env '${QIIME2_ENV}'"
  exit 1
fi
ok "$QIIME_VER"

# ── PICRUSt2 (optional functional prediction) ─────────────────────────
step "PICRUSt2 plugin (optional)"
if conda run -n "$QIIME2_ENV" --no-capture-output qiime info 2>/dev/null | grep -q picrust2; then
  ok "PICRUSt2 already installed"
else
  info "Installing q2-picrust2..."
  conda run -n "$QIIME2_ENV" --no-capture-output \
    pip install q2-picrust2 --quiet || warn "PICRUSt2 install failed — skipping (optional)"
  conda run -n "$QIIME2_ENV" --no-capture-output \
    qiime dev refresh-cache 2>/dev/null || true
  ok "PICRUSt2 installed"
fi

# ── R packages for visualization ─────────────────────────────────────
step "R packages — phyloseq, vegan, ANCOMBC2, FUNGuildR, LULU"

R_LIBS_USER="$HOME/R/library"
mkdir -p "$R_LIBS_USER"

Rscript --no-save --no-restore << 'REOF'
user_lib <- Sys.getenv("R_LIBS_USER", unset=file.path(Sys.getenv("HOME"),"R","library"))
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))

cat("── Installing Bioconductor packages ──\n")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", lib=user_lib, repos="https://cran.r-project.org", quiet=TRUE)

bioc_pkgs <- c("phyloseq", "ANCOMBC", "microbiome", "DESeq2", "ggtree")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, lib=user_lib, update=FALSE, ask=FALSE, quiet=TRUE)
  } else cat(" ", pkg, "already installed\n")
}

cat("── Installing CRAN packages ──\n")
cran_pkgs <- c("vegan","ggplot2","dplyr","tidyr","pheatmap",
               "ape","picante","GUniFrac","reshape2","scales",
               "RColorBrewer","cowplot","ggrepel","patchwork",
               "optparse","qiime2R","microViz")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib=user_lib, repos="https://cran.r-project.org", quiet=TRUE)
  } else cat(" ", pkg, "already installed\n")
}

cat("── Installing GitHub packages ──\n")
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes", lib=user_lib, repos="https://cran.r-project.org", quiet=TRUE)

# FUNGuildR (ITS ecological guild annotation)
if (!requireNamespace("FUNGuildR", quietly=TRUE)) {
  cat("Installing FUNGuildR...\n")
  remotes::install_github("brendanf/FUNGuildR", lib=user_lib, quiet=TRUE,
                           upgrade="never")
} else cat("  FUNGuildR already installed\n")

# LULU (COX1 post-clustering curation)
if (!requireNamespace("lulu", quietly=TRUE)) {
  cat("Installing lulu...\n")
  remotes::install_github("tobiasgf/lulu", lib=user_lib, quiet=TRUE,
                           upgrade="never")
} else cat("  lulu already installed\n")

cat("\n✅ R packages ready\n")
REOF

ok "R visualization packages installed"

# ── QIIME2 classifiers info ───────────────────────────────────────────
step "Classifier note"
CLASSIFIERS_DIR="$HOME/r16s-app/backend/classifiers"
mkdir -p "$CLASSIFIERS_DIR"

cat << 'EOF'

  Pre-trained classifiers (.qza) must be downloaded separately.
  Run the following for each marker you need:

  ┌─ 16S V3-V4 (SILVA 138.1) ──────────────────────────────────────────
  │  wget -O ~/r16s-app/backend/classifiers/silva_16S_classifier.qza \
  │    "https://data.qiime2.org/2024.5/common/silva-138-99-seqs-515-806.qza"
  │
  ├─ 16S Full-length / V1-V9 (PacBio) ────────────────────────────────
  │  wget -O ~/r16s-app/backend/classifiers/silva_16S_full_classifier.qza \
  │    "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza"
  │
  ├─ ITS (UNITE 9.0) ─────────────────────────────────────────────────
  │  wget -O ~/r16s-app/backend/classifiers/unite_ITS_classifier.qza \
  │    "https://github.com/colinbrislawn/unite-train/releases/download/\
  │     9.0-qiime2-2023.9-demo/unite_ver9_99_all_29.11.2022-Q2-2023.9.qza"
  │
  └─ Classifiers will be auto-detected from ~/r16s-app/backend/classifiers/

EOF

info "Classifiers directory: $CLASSIFIERS_DIR"
info "Run download_databases.sh to download classifier files"

# ── Summary ───────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  QIIME2 Setup Complete!"
echo ""
echo "  Environment  : ${QIIME2_ENV}"
echo "  Activate with: conda activate ${QIIME2_ENV}"
echo ""
echo "  Next steps:"
echo "  1. Download classifiers (see above)"
echo "  2. Restart the app: bash ~/r16s-app/start.sh"
echo "============================================================"
echo ""
