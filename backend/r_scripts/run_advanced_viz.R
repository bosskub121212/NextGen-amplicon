#!/usr/bin/env Rscript
# =============================================================================
#  run_advanced_viz.R — Wrapper: run all advanced visualizations for a job
#  Called by main.py's Re-run Viz (and optionally at end of pipeline)
#
#  Usage: Rscript run_advanced_viz.R <output_dir>
#
#  Runs (in order, each non-fatal):
#    1. clustering_heatmap.R  — always (if taxonomy CSV present)
#    2. venn_diagram.R        — only if metadata.csv present
#    3. phylo_tree.R          — only if tree source or bray-curtis present
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

# Locate this script's directory (where sibling scripts live)
this_file  <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(this_file) > 0) {
  dirname(normalizePath(sub("--file=", "", this_file[1])))
} else {
  dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE))
}

run_script <- function(name) {
  path <- file.path(script_dir, name)
  if (!file.exists(path)) {
    cat("[run_advanced_viz] SKIP", name, "(file not found)\n"); return(invisible(NULL))
  }
  cat("[run_advanced_viz] Running", name, "...\n")
  tryCatch(
    source(path),   # source so we inherit the out_dir set above
    error = function(e) cat("[run_advanced_viz] ERROR in", name, ":", e$message, "\n")
  )
}

# ── 1. Hierarchical Clustering Heatmap ──────────────────────────────────────
run_script("clustering_heatmap.R")

# ── 2. Venn Diagram (requires metadata.csv) ──────────────────────────────────
run_script("venn_diagram.R")

# ── 3. Phylogenetic Tree ─────────────────────────────────────────────────────
run_script("phylo_tree.R")

cat("[run_advanced_viz] All done.\n")
