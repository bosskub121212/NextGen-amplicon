#!/usr/bin/env bash
# ============================================================
#  NextGen-Amplicon — LULU (COI post-clustering curation)
#  Installs the R package LULU (github.com/tobiasgf/lulu) into the R the pipelines run with, and
#  checks BLAST+ (makeblastdb / blastn), which LULU's match list needs. Safe to run again.
#
#  Run on the server, in the same shell you start the app from:
#     bash ~/r16s-app/install_lulu.sh
# ============================================================
set -uo pipefail
echo "== LULU for the COI pipeline"
RS=$(command -v Rscript || true)
if [ -z "$RS" ]; then echo "  Rscript not found on PATH — start the shell the app runs in"; exit 1; fi
echo "  R: $RS"
for t in makeblastdb blastn; do
  if command -v "$t" >/dev/null 2>&1; then echo "  $t: $(command -v $t)"; else
    echo "  $t not found — LULU needs BLAST+:  sudo apt-get install -y ncbi-blast+"; fi
done
Rscript - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))
ok <- function(p) requireNamespace(p, quietly = TRUE)
if (ok("lulu")) { cat("  lulu already installed:", as.character(packageVersion("lulu")), "\n"); quit(status = 0) }
if (!ok("remotes")) install.packages("remotes")
if (!ok("dplyr")) install.packages("dplyr")             # LULU's only dependency
tryCatch(remotes::install_github("tobiasgf/lulu", upgrade = "never"),
         error = function(e) cat("  install_github failed:", conditionMessage(e), "\n"))
if (ok("lulu")) {
  cat("  lulu installed:", as.character(packageVersion("lulu")), "\n")
} else {
  cat("  LULU is NOT installed. Library paths:", paste(.libPaths(), collapse = " | "), "\n")
  quit(status = 1)
}
REOF
rc=$?
[ $rc -eq 0 ] && echo "  Done: COI jobs run LULU from now on (Run again for jobs that skipped it)."
exit $rc
