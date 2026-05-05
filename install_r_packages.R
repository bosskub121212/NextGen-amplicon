# Run this in WSL2 with: Rscript install_r_packages.R
# Installs all packages needed for the comprehensive report plots

pkgs <- c("ggplot2", "reshape2", "vegan", "ape", "pheatmap", "scales")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
    if (requireNamespace(pkg, quietly=TRUE)) {
      cat("  OK:", pkg, "\n")
    } else {
      cat("  FAILED:", pkg, "\n")
    }
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

cat("\nDone! Packages available:\n")
for (pkg in pkgs) cat(" ", pkg, "—", if(requireNamespace(pkg, quietly=TRUE)) "OK" else "MISSING", "\n")
