#!/usr/bin/env Rscript
# =============================================================================
#  run_extra_viz.R — Standalone caller for dada2_extra_viz.R
#  Called by main.py's Re-run Viz for any job that has asv_table.csv
#  Usage: Rscript run_extra_viz.R <output_dir>
# =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else getwd()

# Locate dada2_extra_viz.R in same directory as this script
this_file <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sdir <- if (length(this_file) > 0) {
  dirname(normalizePath(sub("--file=", "", this_file[1])))
} else {
  dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE))
}
xscript <- file.path(sdir, "dada2_extra_viz.R")

if (!file.exists(xscript)) {
  cat("[run_extra_viz] ERROR: dada2_extra_viz.R not found at", xscript, "\n")
  quit(status = 1)
}

source(xscript)
dada2_extra_viz(out_dir)
