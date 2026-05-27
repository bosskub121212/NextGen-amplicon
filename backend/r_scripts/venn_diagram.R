#!/usr/bin/env Rscript
# =============================================================================
#  venn_diagram.R — Venn Diagram of shared/unique ASVs per sample group
#
#  Usage: Rscript venn_diagram.R <output_dir>
#
#  Reads:  asv_table.csv   (ASVs × samples, or samples × ASVs)
#          metadata.csv    (SampleID, Group, ...)
#  Writes: r_plots/venn_diagram.pdf
#          venn_counts.csv  (shared/unique counts per group combination)
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── Venn Diagram ────────────────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Input files ───────────────────────────────────────────────────────────────
asv_file  <- file.path(out_dir, "asv_table.csv")
meta_file <- file.path(out_dir, "metadata.csv")

if (!file.exists(asv_file)) {
  cat("  [skip] asv_table.csv not found\n"); quit(status = 0)
}
if (!file.exists(meta_file)) {
  cat("  [skip] metadata.csv not found — cannot group samples\n"); quit(status = 0)
}

# ── Read data ─────────────────────────────────────────────────────────────────
asv_raw <- tryCatch(
  read.csv(asv_file, check.names = FALSE),
  error = function(e) { cat("  [skip] Cannot read asv_table.csv:", e$message, "\n"); NULL }
)
if (is.null(asv_raw)) quit(status = 0)

meta <- tryCatch(
  read.csv(meta_file, stringsAsFactors = FALSE),
  error = function(e) { cat("  [skip] Cannot read metadata.csv:", e$message, "\n"); NULL }
)
if (is.null(meta)) quit(status = 0)

if (!"SampleID" %in% colnames(meta)) {
  cat("  [skip] metadata.csv missing 'SampleID' column\n"); quit(status = 0)
}
if (!"Group" %in% colnames(meta)) {
  cat("  [skip] metadata.csv missing 'Group' column\n"); quit(status = 0)
}
rownames(meta) <- meta$SampleID

# ── Orient matrix: rows=ASVs, cols=samples ────────────────────────────────────
# Drop non-numeric or sequence columns
seq_col <- which(colnames(asv_raw) == "sequence")
if (length(seq_col) > 0) asv_raw <- asv_raw[, -seq_col, drop = FALSE]

# If first column is character (ASV IDs), set as rownames
if (is.character(asv_raw[[1]]) || all(is.na(suppressWarnings(as.numeric(asv_raw[[1]]))))) {
  rownames(asv_raw) <- asv_raw[[1]]
  asv_raw           <- asv_raw[, -1, drop = FALSE]
}

mat <- as.matrix(asv_raw)
storage.mode(mat) <- "numeric"
mat[is.na(mat)] <- 0

# Detect orientation: ncols > nrows → samples-as-columns (typical DADA2 wide)
# Already expected: rows=ASVs, cols=samples
# But if samples are rows (nrows >> ncols), transpose
if (nrow(mat) < ncol(mat)) mat <- t(mat)   # now rows=ASVs, cols=samples

cat("  ASV matrix:", nrow(mat), "ASVs ×", ncol(mat), "samples\n")

# ── Group samples by metadata$Group ──────────────────────────────────────────
groups <- unique(meta$Group)
groups <- groups[!is.na(groups) & nchar(trimws(groups)) > 0]

# Limit to 7 groups for Venn (ggVennDiagram supports 2–7)
if (length(groups) > 7) {
  cat("  [warn] More than 7 groups — using first 7\n")
  groups <- groups[1:7]
}
if (length(groups) < 2) {
  cat("  [skip] Need at least 2 groups for Venn diagram\n"); quit(status = 0)
}

# For each group, collect all ASVs present (count > 0) in at least 1 sample
group_asvs <- lapply(groups, function(g) {
  samp <- meta$SampleID[meta$Group == g]
  samp <- intersect(samp, colnames(mat))
  if (length(samp) == 0) return(character(0))
  sub  <- mat[, samp, drop = FALSE]
  rownames(sub)[rowSums(sub) > 0]
})
names(group_asvs) <- groups

# Report sizes
for (g in groups)
  cat("  Group", g, "—", length(group_asvs[[g]]), "ASVs present\n")

# ── Output paths ──────────────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
pdf_out   <- file.path(plots_dir, "venn_diagram.pdf")

# ── Try ggVennDiagram first, fall back to VennDiagram ─────────────────────────
has_ggvenn  <- requireNamespace("ggVennDiagram", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2",       quietly = TRUE)
has_venn    <- requireNamespace("VennDiagram",   quietly = TRUE)

if (has_ggvenn && has_ggplot2) {
  suppressPackageStartupMessages({
    library(ggVennDiagram)
    library(ggplot2)
  })
  tryCatch({
    p <- ggVennDiagram(group_asvs, label_alpha = 0,
                       set_size = 4, label_size = 3.5) +
      scale_fill_gradient(low = "#f7fbff", high = "#2171b5",
                          name = "ASV count") +
      scale_color_manual(values = rep("grey30", length(groups))) +
      theme_void(base_size = 12) +
      labs(title = "Shared and Unique ASVs by Group",
           subtitle = paste0("Total unique ASVs: ",
                             length(unique(unlist(group_asvs)))))
    ggsave(pdf_out, p, width = 10, height = 8)
    cat("  ✓ venn_diagram.pdf  (ggVennDiagram)\n")
  }, error = function(e) cat("  [error] ggVennDiagram failed:", e$message, "\n"))

} else if (has_venn) {
  suppressPackageStartupMessages(library(VennDiagram))
  tryCatch({
    fill_cols <- c("#e41a1c","#377eb8","#4daf4a","#984ea3",
                   "#ff7f00","#a65628","#f781bf")[seq_along(groups)]
    venn.diagram(
      group_asvs,
      filename  = pdf_out,
      imagetype = "pdf",
      fill      = alpha(fill_cols, 0.5),
      col       = fill_cols,
      cat.col   = fill_cols,
      cat.cex   = 1.2, cex = 1.1,
      main      = "Shared and Unique ASVs by Group",
      main.cex  = 1.3
    )
    cat("  ✓ venn_diagram.pdf  (VennDiagram)\n")
  }, error = function(e) cat("  [error] VennDiagram failed:", e$message, "\n"))

} else {
  cat("  [skip] Neither 'ggVennDiagram' nor 'VennDiagram' installed\n")
  cat("  Install with: install.packages('ggVennDiagram')\n")
  quit(status = 0)
}

# ── Save counts CSV (for reference / interactive summary) ─────────────────────
tryCatch({
  # Build pairwise intersection matrix
  n  <- length(groups)
  counts_df <- data.frame(Group = groups, Count = sapply(group_asvs, length))
  write.csv(counts_df, file.path(out_dir, "venn_counts.csv"), row.names = FALSE)
  cat("  ✓ venn_counts.csv\n")
}, error = function(e) invisible(NULL))

cat("── Venn Diagram done ───────────────────────────────────────────────\n\n")
