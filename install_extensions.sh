#!/usr/bin/env bash
# ============================================================
#  NextGen-Amplicon — Extension Installer
#  Installs all R packages and Python tools needed for:
#  ITS fungi, COX1, Nematode 18S, Tax4Fun2, PICRUSt2, PacBio
#
#  Run from WSL:  bash /mnt/c/Claude/r16s-app/install_extensions.sh
# ============================================================
set -euo pipefail

APP_DIR="$HOME/r16s-app"
VENV="$APP_DIR/venv"
LOG="$APP_DIR/install_extensions.log"

echo ""
echo "============================================================"
echo "  NextGen-Amplicon Extension Installer"
echo "============================================================"
echo "  Log: $LOG"
echo ""

# Activate Python venv
source "$VENV/bin/activate"

# ── 1. R packages ────────────────────────────────────────────
echo "[1/5] Installing R packages..."
Rscript - << 'REOF' 2>&1 | tee -a "$LOG"

options(repos = c(CRAN = "https://cloud.r-project.org"))
cat("=== R Package Installation ===\n")

pkg_install <- function(pkg, source="cran", ...) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat(sprintf("Installing %s ...\n", pkg))
    tryCatch({
      if (source == "cran")    install.packages(pkg, ...)
      if (source == "bioc")    BiocManager::install(pkg, ask=FALSE, update=FALSE)
      if (source == "github")  remotes::install_github(pkg, ...)
      cat(sprintf("  OK: %s\n", pkg))
    }, error=function(e) cat(sprintf("  FAILED: %s — %s\n", pkg, e$message)))
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

# Core deps
pkg_install("remotes")
pkg_install("BiocManager")

# Bioconductor packages
pkg_install("Biostrings",   source="bioc")
pkg_install("ShortRead",    source="bioc")
pkg_install("DECIPHER",     source="bioc")

# CRAN packages for new pipelines
pkg_install("optparse")
pkg_install("jsonlite")
pkg_install("ggplot2")
pkg_install("reshape2")
pkg_install("vegan")
pkg_install("ape")
pkg_install("phangorn")
pkg_install("seqinr")       # COX1 codon translation helper
pkg_install("stringr")

# LULU — post-clustering curation (COX1 NUMTs + general)
pkg_install("lulu", source="github", pkg="tobiasgf/lulu")

# FUNGuildR — fungal ecological guilds (ITS pipeline)
pkg_install("FUNGuildR", source="github", pkg="brendanf/FUNGuildR")

# Tax4Fun2 — functional prediction for 16S
if (!requireNamespace("Tax4Fun2", quietly=TRUE)) {
  cat("Installing Tax4Fun2 from GitHub...\n")
  tryCatch({
    remotes::install_github("bwemheu/Tax4Fun2")
    cat("  OK: Tax4Fun2\n")
  }, error=function(e) {
    cat(sprintf("  Tax4Fun2 GitHub failed, trying source install: %s\n", e$message))
  })
} else {
  cat("  Already installed: Tax4Fun2\n")
}

cat("\n=== R Package Installation Complete ===\n")
REOF

echo ""
echo "[2/5] Installing ITSxpress (Python, ITS primer trimming)..."
pip install itsxpress --quiet && echo "  OK: itsxpress" || echo "  FAILED: itsxpress"

echo ""
echo "[3/5] Installing vsearch (required by ITSxpress)..."
if ! command -v vsearch &>/dev/null; then
  sudo apt-get install -y vsearch 2>/dev/null || \
  conda install -y -c bioconda vsearch 2>/dev/null || \
  echo "  WARNING: vsearch not installed — install manually with: sudo apt install vsearch"
else
  echo "  Already installed: vsearch"
fi

echo ""
echo "[4/5] Installing PICRUSt2 via conda..."
if command -v conda &>/dev/null; then
  if conda env list | grep -q "picrust2"; then
    echo "  Conda env 'picrust2' already exists."
  else
    echo "  Creating conda env 'picrust2'..."
    conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y 2>&1 | tail -5 \
      && echo "  OK: picrust2 conda env created" \
      || echo "  FAILED: picrust2 install — install manually: conda create -n picrust2 -c bioconda picrust2"
  fi
else
  echo "  WARNING: conda not found."
  echo "  To install PICRUSt2 later:"
  echo "    conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y"
fi

echo ""
echo "[5/5] Installing BLAST+ (required by LULU match list generation)..."
if ! command -v makeblastdb &>/dev/null; then
  sudo apt-get install -y ncbi-blast+ 2>/dev/null \
    && echo "  OK: BLAST+" \
    || echo "  WARNING: BLAST+ not installed. Install with: sudo apt-get install ncbi-blast+"
else
  echo "  Already installed: BLAST+"
fi

echo ""
echo "============================================================"
echo "  Installation complete!"
echo "  Next step: run download_databases.sh to get reference DBs"
echo "============================================================"
echo ""
