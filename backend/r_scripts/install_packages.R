#!/usr/bin/env Rscript
# =============================================================================
#  NextGen-Amplicon — R Package Installer
#  Run once in WSL: Rscript ~/r16s-app/backend/r_scripts/install_packages.R
# =============================================================================

cat("=== NextGen-Amplicon: Installing R packages ===\n\n")

# ── User library setup ────────────────────────────────────────────────────────
user_lib <- Sys.getenv("R_LIBS_USER",
             unset=file.path(Sys.getenv("HOME"), "R", "library"))
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))
cat(sprintf("Installing to: %s\n\n", user_lib))

install_if_missing <- function(pkg, bioc=FALSE) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat(sprintf("Installing: %s ...\n", pkg))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager", repos="https://cloud.r-project.org",
                         lib=user_lib, quiet=TRUE)
      BiocManager::install(pkg, lib=user_lib, ask=FALSE, update=FALSE,
                           quiet=TRUE, force=FALSE)
    } else {
      install.packages(pkg, repos="https://cloud.r-project.org",
                       lib=user_lib, quiet=TRUE)
    }
    if (requireNamespace(pkg, quietly=TRUE))
      cat(sprintf("  ✓ %s installed\n", pkg))
    else
      cat(sprintf("  ✗ %s FAILED\n", pkg))
  } else {
    cat(sprintf("  ✓ %s already present\n", pkg))
  }
}

# ── CRAN packages ─────────────────────────────────────────────────────────────
cat("── CRAN packages ────────────────────────────────────────────\n")
cran_pkgs <- c(
  # Low-level deps (must come first — Bioconductor packages depend on these)
  "bitops",          # Rsamtools dependency
  "abind",           # SummarizedExperiment dependency
  "matrixStats",     # DelayedArray / DESeq2 dependency
  "snow",            # parallel backend
  "RcppEigen",       # many Bioc packages
  "RcppParallel",    # dada2 dependency
  "deldir",          # ShortRead dependency
  "png", "jpeg",     # ShortRead dependency
  "interp",          # ShortRead dependency
  "hwriter",         # ShortRead dependency
  "latticeExtra",    # ShortRead dependency
  # Core
  "optparse", "jsonlite", "ggplot2", "dplyr", "tidyr",
  "RColorBrewer", "ggrepel", "patchwork", "cowplot", "pheatmap",
  "vegan", "picante",
  # Phase 1 additions
  "ggVennDiagram",   # Venn diagrams
  "scales",          # axis formatting
  # Phase 2 additions
  "ggdendro",        # UPGMA dendrogram plotting
  "ggtern",          # Ternary plots (3-group composition)
  # Phase 3 additions
  "igraph",          # Co-occurrence network
  "ggraph",          # Network visualization
  "Hmisc",           # Spearman correlation (rcorr)
  "phangorn",        # Phylogenetic tree building (NJ tree for UniFrac)
  "ape"              # Tree manipulation
)
for (p in cran_pkgs) install_if_missing(p, bioc=FALSE)

# ── Bioconductor packages ─────────────────────────────────────────────────────
cat("\n── Bioconductor packages ────────────────────────────────────\n")
bioc_pkgs <- c(
  # Core (should already be installed)
  "phyloseq", "dada2", "ANCOMBC", "DESeq2",
  # Phase 1 additions
  "lefser",              # LEfSe analysis
  "SummarizedExperiment",
  "S4Vectors",
  # Phase 3b additions (UniFrac + tree viz)
  "DECIPHER",            # Multiple sequence alignment (for NJ tree)
  "Biostrings",          # DNA sequence handling
  "ggtree",              # Phylogenetic tree visualization
  "treeio",              # Tree I/O (ggtree dependency)
  # v2.6.0 additions
  "microbiomeMarker"     # LEfSe cladogram (optional, falls back gracefully)
)
for (p in bioc_pkgs) install_if_missing(p, bioc=TRUE)

# ── Summary ───────────────────────────────────────────────────────────────────
# ── Optional: FUNGuildR (ITS pipelines only — install from GitHub) ────────────
cat("\n── Optional packages (ITS only) ─────────────────────────────\n")
if (!requireNamespace("FUNGuildR", quietly=TRUE)) {
  cat("  [optional] FUNGuildR not installed.\n")
  cat("  To install (ITS pipelines only):\n")
  cat("    Rscript -e \"remotes::install_github('brendanf/FUNGuildR')\"\n")
} else {
  cat("  ✓ FUNGuildR already present\n")
}

cat("\n=== Installation complete ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
ok  <- sapply(all_pkgs, function(p) requireNamespace(p, quietly=TRUE))
cat(sprintf("  Available   : %d / %d\n", sum(ok), length(ok)))
if (any(!ok))
  cat(sprintf("  Missing     : %s\n", paste(names(ok)[!ok], collapse=", ")))
cat("\n")
