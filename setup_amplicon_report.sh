#!/usr/bin/env bash
# =============================================================================
#  setup_amplicon_report.sh — tools for the per-sample customer report (12S / CytB)
#
#    bash ~/r16s-app/setup_amplicon_report.sh            # install what is missing
#    bash ~/r16s-app/setup_amplicon_report.sh --status   # show what is installed
#
#  The lab's customer report names these tools; each is installed here, and the
#  report states which one actually ran:
#    conda env ngamp-amp   FastQC, fastp, BLAST+, KronaTools (krona), pandas, pyarrow
#                          (also made by download_databases.sh --fishmito)
#    R  microbiomeutilities  Top 50 ASV heatmap (GitHub microsud/microbiomeutilities)
#    R  venneuler            Venn diagram (needs Java / rJava); eulerr is installed as
#                            the fallback when venneuler cannot be
#    backend venv openpyxl   Excel workbooks
#  Nothing existing is changed; re-running only adds what is missing.
# =============================================================================
set -uo pipefail
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="$APP/venv/bin/python"; [[ -x "$PY" ]] || PY="$(command -v python3)"
G='\033[0;32m'; Y='\033[1;33m'; R='\033[0;31m'; B='\033[1m'; N='\033[0m'
ok(){ echo -e "  ${G}✓${N} $*"; }; warn(){ echo -e "  ${Y}!${N} $*"; }; bad(){ echo -e "  ${R}✗${N} $*"; }
step(){ echo -e "\n${B}══ $* ══${N}"; }
STATUS=0; [[ "${1:-}" == "--status" ]] && STATUS=1

BASE=""; CONDA=""
for c in mamba conda "$HOME/miniconda3/bin/conda"; do
  if command -v "$c" >/dev/null 2>&1 || [[ -x "$c" ]]; then
    b=$("$c" info --base 2>/dev/null | grep -o '/[^ ]*' | tail -1)
    [[ -n "$b" && -d "$b/envs" ]] && { BASE="$b"; CONDA="$c"; break; }
  fi
done
AMP="${BASE:-$HOME/miniconda3}/envs/ngamp-amp/bin"
NANO="${BASE:-$HOME/miniconda3}/envs/ngamp-nano/bin"

rpkg(){ Rscript -e "cat(requireNamespace('$1', quietly=TRUE))" 2>/dev/null | grep -q TRUE; }
status(){
  for t in fastqc fastp blastn ktImportText; do
    [[ -x "$AMP/$t" ]] && ok "$t ($AMP)" || warn "$t missing"
  done
  [[ -x "$NANO/NanoPlot" ]] && ok "NanoPlot ($NANO)" || warn "NanoPlot missing (ONT 16S customer report: read statistics without NanoPlot plots)"
  for p in microbiomeutilities venneuler eulerr; do rpkg "$p" && ok "R $p" || warn "R $p missing"; done
  "$PY" -c "import openpyxl" 2>/dev/null && ok "openpyxl (backend)" || warn "openpyxl missing in $PY"
}
if [[ $STATUS == 1 ]]; then status; exit 0; fi

step "conda env ngamp-amp (FastQC, fastp, BLAST+, krona)"
if [[ -x "$AMP/fastqc" && -x "$AMP/fastp" && -x "$AMP/blastn" && -x "$AMP/ktImportText" ]]; then
  ok "already installed"
elif [[ -d "${BASE}/envs/ngamp-amp" ]]; then
  "$CONDA" install -y -n ngamp-amp -c conda-forge -c bioconda fastqc fastp blast krona pandas pyarrow || bad "install into ngamp-amp failed"
else
  "$CONDA" create -y -n ngamp-amp -c conda-forge -c bioconda python=3.12 fastqc fastp blast krona pandas pyarrow \
    || bad "conda env ngamp-amp failed"
fi

step "conda env ngamp-nano (NanoPlot, for ONT 16S reads)"
# Own env: NanoPlot pins its plotting stack, which should not move FastQC / BLAST.
if [[ -x "$NANO/NanoPlot" ]]; then
  ok "already installed"
else
  "$CONDA" create -y -n ngamp-nano -c conda-forge -c bioconda nanoplot \
    || warn "NanoPlot not installed — the ONT report still gives read statistics, without NanoPlot plots"
fi

step "R packages"
Rscript - <<'REOF'
lib <- Sys.getenv("R_LIBS_USER"); if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
repos <- "https://cloud.r-project.org"
need <- function(p) !requireNamespace(p, quietly = TRUE)
if (need("remotes")) install.packages("remotes", repos = repos, quiet = TRUE)
if (need("microbiomeutilities"))
  tryCatch(remotes::install_github("microsud/microbiomeutilities", upgrade = "never", quiet = TRUE),
           error = function(e) cat("  microbiomeutilities failed:", conditionMessage(e), "\n"))
if (need("venneuler"))
  tryCatch(install.packages("venneuler", repos = repos, quiet = TRUE),
           error = function(e) cat("  venneuler failed (needs Java):", conditionMessage(e), "\n"))
if (need("venneuler") && need("eulerr")) install.packages("eulerr", repos = repos, quiet = TRUE)
for (p in c("microbiomeutilities", "venneuler", "eulerr"))
  cat(sprintf("  %-20s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "installed" else "missing"))
REOF

step "openpyxl (Excel)"
"$PY" -c "import openpyxl" 2>/dev/null || "$PY" -m pip install --quiet openpyxl || bad "pip install openpyxl failed"

step "Status"
status
echo ""
ok "Done. No backend restart needed."
