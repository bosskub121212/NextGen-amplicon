#!/usr/bin/env Rscript
# =============================================================================
#  clustering_heatmap.R — Hierarchical Clustering Heatmap
#  Produces: PDF with dendrogram + annotation bar  +  CSV for interactive tab
#
#  Usage: Rscript clustering_heatmap.R <output_dir>
#
#  Reads:  taxonomy_*.csv  (species > genus > family > order > class > phylum)
#          metadata.csv    (SampleID, Group, [optional extra columns])
#  Writes: r_plots/clustering_heatmap.pdf
#          clustering_heatmap.csv   (z-score matrix, used by interactive tab)
#          clustering_heatmap_meta.csv  (sample → group mapping, for tab colors)
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── Clustering Heatmap ──────────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Required packages ────────────────────────────────────────────────────────
needed <- c("pheatmap", "RColorBrewer")
missing_pkgs <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("  [skip] Missing packages:", paste(missing_pkgs, collapse=", "), "\n")
  cat("  Install with: install.packages(c('pheatmap','RColorBrewer'))\n")
  quit(status = 0)
}
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

# ── Find best taxonomy file ──────────────────────────────────────────────────
tax_levels <- c("taxonomy_species.csv", "taxonomy_genus.csv",
                "taxonomy_family.csv",  "taxonomy_order.csv",
                "taxonomy_class.csv",   "taxonomy_phylum.csv")
tax_file   <- NULL
tax_label  <- "taxa"
for (f in tax_levels) {
  fp <- file.path(out_dir, f)
  if (file.exists(fp)) {
    tax_file  <- fp
    tax_label <- sub("taxonomy_", "", sub("\\.csv$", "", f))
    break
  }
}
if (is.null(tax_file)) {
  cat("  [skip] No taxonomy CSV found in output directory\n")
  quit(status = 0)
}
cat("  Taxonomy level:", tax_label, "\n")

# ── Read taxonomy matrix ─────────────────────────────────────────────────────
mat_raw <- tryCatch(
  read.csv(tax_file, row.names = 1, check.names = FALSE),
  error = function(e) { cat("  [skip] Cannot read", tax_file, ":", e$message, "\n"); NULL }
)
if (is.null(mat_raw) || nrow(mat_raw) == 0 || ncol(mat_raw) == 0) {
  cat("  [skip] Empty taxonomy matrix\n"); quit(status = 0)
}

# Detect orientation: if #cols >= #rows → taxa-as-columns (regular 16S) → transpose
if (ncol(mat_raw) >= nrow(mat_raw)) {
  mat_raw <- t(mat_raw)   # now: rows=taxa, cols=samples
}
# Ensure numeric
mat_num <- apply(mat_raw, 2, function(x) as.numeric(as.character(x)))
rownames(mat_num) <- rownames(mat_raw)
mat_num[is.na(mat_num)] <- 0

cat("  Matrix:", nrow(mat_num), "taxa ×", ncol(mat_num), "samples\n")

# ── Convert to relative abundance (%) ────────────────────────────────────────
col_sums <- colSums(mat_num)
col_sums[col_sums == 0] <- 1  # avoid div-by-zero
mat_rel <- sweep(mat_num, 2, col_sums, "/") * 100

# ── Select top N taxa by total abundance ─────────────────────────────────────
n_taxa   <- min(80, nrow(mat_rel))
row_tots <- rowSums(mat_rel)
top_idx  <- order(row_tots, decreasing = TRUE)[seq_len(n_taxa)]
mat_top  <- mat_rel[top_idx, , drop = FALSE]

# Trim long taxon names for readability (max 40 chars)
rownames(mat_top) <- substr(rownames(mat_top), 1, 40)

# ── Z-score scale (per taxon) ─────────────────────────────────────────────────
mat_scaled <- t(scale(t(mat_top)))
mat_scaled[is.na(mat_scaled)] <- 0   # constant rows → 0
# Clamp to ±3 for colour saturation
mat_scaled <- pmax(pmin(mat_scaled, 3), -3)
cat("  Top taxa selected:", nrow(mat_scaled), "\n")

# ── Metadata annotation ──────────────────────────────────────────────────────
meta_file      <- file.path(out_dir, "metadata.csv")
annotation_col <- NULL
ann_colors     <- NULL

if (file.exists(meta_file)) {
  meta <- tryCatch(
    read.csv(meta_file, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(meta) && "SampleID" %in% colnames(meta)) {
    rownames(meta) <- meta$SampleID
    meta$SampleID  <- NULL

    # Keep only samples present in matrix columns
    common_samp <- intersect(colnames(mat_scaled), rownames(meta))
    if (length(common_samp) > 0) {
      annotation_col <- meta[common_samp, , drop = FALSE]
      # Factor all columns for discrete colour palettes
      annotation_col[] <- lapply(annotation_col, function(x) {
        x <- as.character(x); x[is.na(x)] <- "Unknown"; factor(x)
      })

      # Build colour list: assign distinct palette per column
      palettes <- c("Set1","Set2","Set3","Paired","Dark2","Accent")
      ann_colors <- list()
      for (ci in seq_along(annotation_col)) {
        col_name <- colnames(annotation_col)[ci]
        lvls     <- levels(annotation_col[[ci]])
        n_lvls   <- length(lvls)
        pal_name <- palettes[((ci - 1) %% length(palettes)) + 1]
        max_pal  <- min(n_lvls, brewer.pal.info[pal_name, "maxcolors"])
        pal_cols <- if (n_lvls <= max_pal) {
          brewer.pal(max(3, n_lvls), pal_name)[seq_len(n_lvls)]
        } else {
          colorRampPalette(brewer.pal(8, pal_name))(n_lvls)
        }
        names(pal_cols)  <- lvls
        ann_colors[[col_name]] <- pal_cols
      }
      cat("  Annotation columns:", paste(colnames(annotation_col), collapse=", "), "\n")
    }
  }
}

# ── Save CSV for interactive tab ──────────────────────────────────────────────
csv_out  <- file.path(out_dir, "clustering_heatmap.csv")
write.csv(mat_scaled, csv_out)
cat("  ✓ clustering_heatmap.csv\n")

# Also save metadata CSV for tab colour rendering
if (!is.null(annotation_col)) {
  meta_out <- file.path(out_dir, "clustering_heatmap_meta.csv")
  out_meta          <- annotation_col
  out_meta$SampleID <- rownames(annotation_col)
  out_meta          <- out_meta[, c("SampleID", setdiff(colnames(out_meta), "SampleID")), drop=FALSE]
  write.csv(out_meta, meta_out, row.names = FALSE)
  cat("  ✓ clustering_heatmap_meta.csv\n")
}

# ── Generate PDF ──────────────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
pdf_out <- file.path(plots_dir, "clustering_heatmap.pdf")

hm_height <- max(10, min(30, nrow(mat_scaled) * 0.18 + 4))
hm_width  <- max(10, min(24, ncol(mat_scaled) * 0.55 + 5))

tryCatch({
  pheatmap(
    mat_scaled,
    color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks           = seq(-3, 3, length.out = 101),
    annotation_col   = annotation_col,
    annotation_colors= ann_colors,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    clustering_method= "ward.D2",
    show_rownames    = TRUE,
    show_colnames    = TRUE,
    fontsize_row     = max(5, min(9, 300 / nrow(mat_scaled))),
    fontsize_col     = 9,
    angle_col        = 45,
    main             = paste("Hierarchical Clustering Heatmap —", tax_label, "(z-score)"),
    filename         = pdf_out,
    width            = hm_width,
    height           = hm_height
  )
  cat("  ✓ clustering_heatmap.pdf\n")
}, error = function(e) {
  cat("  [error] PDF generation failed:", e$message, "\n")
})

cat("── Clustering Heatmap done ─────────────────────────────────────────\n\n")
