#!/usr/bin/env Rscript
# =============================================================================
#  run_advanced_viz.R — Wrapper: run all advanced visualizations for a job
#  Called by main.py's Re-run Viz endpoint (and at end of pipeline optionally)
#
#  Usage: Rscript run_advanced_viz.R <output_dir>
#
#  Scripts run (in order, each non-fatal):
#    v2.5.0  clustering_heatmap.R    — pheatmap hierarchical clustering
#    v2.5.0  venn_diagram.R          — ggVennDiagram (needs metadata.csv)
#    v2.5.0  phylo_tree.R            — ggtree circular tree
#    v2.6.0  dada2_extra_viz.R       — PCoA/NMDS ellipses + Group column patch
#    v2.6.0  clust_dendro.R          — dendrogram + composition bars
#    v2.6.0  differential_analysis.R — ANOVA/KW + pairwise Wilcoxon (metastats)
#    v2.6.0  lefse_analysis.R        — LEfSe LDA scores + cladogram
#    v2.6.0  faprotax_analysis.R     — FAPROTAX ecological function prediction
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("════════════════════════════════════════════════════════════════════\n")
cat(" run_advanced_viz.R  —  out_dir:", out_dir, "\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

# ── Locate sibling scripts ─────────────────────────────────────────────────
this_file  <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(this_file) > 0) {
  dirname(normalizePath(sub("--file=", "", this_file[1])))
} else {
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
           error = function(e) dirname(normalizePath(commandArgs()[1], mustWork = FALSE)))
}
cat("  Script dir:", script_dir, "\n\n")

# ── Helper: run a sibling script as a subprocess ───────────────────────────
# Using Rscript subprocess ensures each script gets its own commandArgs()
run_script <- function(name) {
  path <- file.path(script_dir, name)
  if (!file.exists(path)) {
    cat("[skip]", name, "— not found\n\n"); return(invisible(NULL))
  }
  cat("▶", name, "\n")
  tryCatch({
    rc <- system2(
      command = "Rscript",
      args    = c("--vanilla", shQuote(path), shQuote(out_dir)),
      stdout  = "", stderr = ""
    )
    if (rc != 0) cat("  [warn]", name, "exited with code", rc, "\n")
    else         cat("  ✓", name, "done\n")
  }, error = function(e) cat("  [error]", name, ":", e$message, "\n"))
  cat("\n")
}

# ── v2.5.0 scripts ────────────────────────────────────────────────────────
run_script("clustering_heatmap.R")
run_script("venn_diagram.R")
run_script("phylo_tree.R")

# ── v2.6.0 scripts ────────────────────────────────────────────────────────
run_script("dada2_extra_viz.R")       # patches pca_scores + nmds_bray with Group + ellipse CSVs
run_script("clust_dendro.R")          # dendrogram segments + composition bars
run_script("differential_analysis.R") # ANOVA/KW results + pairwise Wilcoxon (metastats)
run_script("lefse_analysis.R")        # LEfSe LDA bar chart + optional cladogram
run_script("faprotax_analysis.R")     # FAPROTAX ecological function prediction

cat("════════════════════════════════════════════════════════════════════\n")
cat(" run_advanced_viz.R  —  all done\n")
cat("════════════════════════════════════════════════════════════════════\n")
