#!/usr/bin/env Rscript
# =============================================================================
#  phylo_tree.R — Circular Phylogenetic Tree
#
#  Usage: Rscript phylo_tree.R <output_dir>
#
#  Tree sources (tries in order):
#    1. phylo_tree.nwk  (copied from QIIME2 export or DADA2 FastTree)
#    2. exports/tree/tree.nwk  (QIIME2 raw export)
#    3. NJ tree from bray_curtis_distance_matrix.csv  (fallback)
#
#  Reads:  phylo_tree.nwk  OR  bray_curtis_distance_matrix.csv
#          taxonomy_*.csv  (for phylum-level coloring of tips)
#          metadata.csv    (SampleID, Group — for outer stacked bars)
#  Writes: r_plots/phylo_tree.pdf
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── Phylogenetic Tree ───────────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Required packages ────────────────────────────────────────────────────────
has_ape    <- requireNamespace("ape",       quietly = TRUE)
has_ggtree <- requireNamespace("ggtree",    quietly = TRUE)
has_ggplot <- requireNamespace("ggplot2",   quietly = TRUE)
has_dplyr  <- requireNamespace("dplyr",     quietly = TRUE)

if (!has_ape) {
  cat("  [skip] 'ape' package not installed\n")
  cat("  Install: install.packages('ape')\n"); quit(status = 0)
}
suppressPackageStartupMessages(library(ape))

# ── Locate tree file ──────────────────────────────────────────────────────────
tree_candidates <- c(
  file.path(out_dir, "phylo_tree.nwk"),
  file.path(out_dir, "exports", "tree", "tree.nwk"),
  file.path(out_dir, "exports", "tree", "tree-rooted.nwk")
)
tree_file <- NULL
for (tc in tree_candidates) {
  if (file.exists(tc)) { tree_file <- tc; break }
}

tree     <- NULL
tree_src <- "none"

if (!is.null(tree_file)) {
  tree     <- tryCatch(read.tree(tree_file), error = function(e) NULL)
  tree_src <- tree_file
  if (!is.null(tree)) cat("  Tree loaded from:", basename(tree_file), "\n")
}

# ── Fallback: NJ tree from Bray-Curtis distance matrix ───────────────────────
if (is.null(tree)) {
  bray_file <- file.path(out_dir, "bray_curtis_distance_matrix.csv")
  if (file.exists(bray_file)) {
    tryCatch({
      dm   <- read.csv(bray_file, row.names = 1, check.names = FALSE)
      dm_m <- as.matrix(dm)
      dm_m[dm_m < 0] <- 0
      tree     <- nj(as.dist(dm_m))
      tree_src <- "NJ from Bray-Curtis distance (sample-level)"
      cat("  Built NJ tree from Bray-Curtis distance matrix\n")
    }, error = function(e) cat("  [warn] NJ tree build failed:", e$message, "\n"))
  }
}

if (is.null(tree)) {
  cat("  [skip] No tree source available\n"); quit(status = 0)
}
cat("  Tree source:", tree_src, "\n")
cat("  Tips:", length(tree$tip.label), "  Nodes:", tree$Nnode, "\n")

# ── Output paths ──────────────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
pdf_out <- file.path(plots_dir, "phylo_tree.pdf")

# ── Determine tip annotation (phylum per tip) ─────────────────────────────────
tip_meta <- data.frame(label = tree$tip.label, Phylum = "Unknown",
                       stringsAsFactors = FALSE)

# Try to map ASV labels → phylum from taxonomy_table.csv
tax_file <- file.path(out_dir, "taxonomy_table.csv")
if (file.exists(tax_file)) {
  tryCatch({
    tax_raw <- read.csv(tax_file, check.names = FALSE, stringsAsFactors = FALSE)
    # Column names vary — look for "Phylum" or parse "taxonomy" string
    if ("Phylum" %in% colnames(tax_raw)) {
      asv_col <- colnames(tax_raw)[1]
      phylum_map <- setNames(tax_raw$Phylum, tax_raw[[asv_col]])
      tip_meta$Phylum <- ifelse(
        tip_meta$label %in% names(phylum_map),
        phylum_map[tip_meta$label], "Unknown")
    } else if ("taxonomy" %in% colnames(tax_raw)) {
      # QIIME2-style: "d__Bacteria; p__Firmicutes; ..."
      asv_col    <- colnames(tax_raw)[1]
      phylum_raw <- sub(".*p__([^;]+).*", "\\1", tax_raw$taxonomy)
      phylum_raw <- trimws(phylum_raw)
      phylum_map <- setNames(phylum_raw, tax_raw[[asv_col]])
      tip_meta$Phylum <- ifelse(
        tip_meta$label %in% names(phylum_map),
        phylum_map[tip_meta$label], "Unknown")
    }
    cat("  Phylum annotation: mapped",
        sum(tip_meta$Phylum != "Unknown"), "/", nrow(tip_meta), "tips\n")
  }, error = function(e) invisible(NULL))
}

# Alternatively, map from taxonomy_phylum.csv if sample names match
if (all(tip_meta$Phylum == "Unknown")) {
  phy_file <- file.path(out_dir, "taxonomy_phylum.csv")
  if (file.exists(phy_file)) {
    tryCatch({
      phy_raw  <- read.csv(phy_file, row.names = 1, check.names = FALSE)
      # This is a sample × phylum abundance table — use phylum names to color samples
      # (useful when tree tips are SAMPLES, not ASVs)
      # Find dominant phylum per sample
      if (ncol(phy_raw) >= nrow(phy_raw)) phy_raw <- t(phy_raw)
      # rows=phyla, cols=samples
      dom_phylum <- apply(phy_raw, 2, function(x) rownames(phy_raw)[which.max(x)])
      tip_meta$Phylum <- ifelse(
        tip_meta$label %in% names(dom_phylum),
        dom_phylum[tip_meta$label], "Unknown")
    }, error = function(e) invisible(NULL))
  }
}

# ── Use ggtree if available (rich circular layout) ───────────────────────────
if (has_ggtree && has_ggplot) {
  suppressPackageStartupMessages({
    library(ggtree)
    library(ggplot2)
  })

  # Phylum palette (up to 15 distinct colours)
  phylums   <- unique(tip_meta$Phylum)
  n_phy     <- length(phylums)
  pal_base  <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                 "#a65628","#f781bf","#999999","#66c2a5","#fc8d62",
                 "#8da0cb","#e78ac3","#a6d854","#ffd92f","#b3b3b3")
  phy_cols  <- setNames(
    if (n_phy <= length(pal_base)) pal_base[seq_len(n_phy)]
    else colorRampPalette(pal_base)(n_phy),
    phylums)

  # Load metadata for outer bars if available
  meta_file  <- file.path(out_dir, "metadata.csv")
  has_meta   <- FALSE
  group_cols <- NULL

  if (file.exists(meta_file)) {
    meta <- tryCatch(read.csv(meta_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(meta) && "SampleID" %in% colnames(meta) && "Group" %in% colnames(meta)) {
      has_meta <- TRUE
    }
  }

  tryCatch({
    # Root tree at midpoint for nicer layout
    tree_rooted <- tryCatch(midpoint(tree), error = function(e) tree)

    p <- ggtree(tree_rooted, layout = "circular", size = 0.3, color = "grey40") %<+%
      tip_meta +
      geom_tippoint(aes(color = Phylum), size = 1.5, alpha = 0.9) +
      geom_tiplab(aes(color = Phylum), size = 1.8, offset = 0.01, align = TRUE) +
      scale_color_manual(values = phy_cols, name = "Phylum") +
      theme_tree2() +
      theme(
        legend.position  = "right",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.4, "cm"),
        plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(title = "Circular Phylogenetic Tree",
           subtitle = paste0("Tips: ", length(tree$tip.label),
                             if (tree_src != "none") paste0("  |  Source: ", basename(tree_src)) else ""))

    h_px <- max(10, min(20, length(tree$tip.label) * 0.12 + 4))
    ggsave(pdf_out, p, width = 16, height = h_px, limitsize = FALSE)
    cat("  ✓ phylo_tree.pdf  (ggtree circular)\n")
  }, error = function(e) {
    cat("  [error] ggtree plot failed:", e$message, "\n")
    # Fallback: basic ape plot
    tryCatch({
      pdf(pdf_out, width = 14, height = 14)
      plot.phylo(tree, type = "fan", cex = 0.4, tip.color = "black",
                 main = "Circular Phylogenetic Tree")
      dev.off()
      cat("  ✓ phylo_tree.pdf  (ape fan plot — ggtree fallback)\n")
    }, error = function(e2) cat("  [error] ape plot also failed:", e2$message, "\n"))
  })

} else {
  # ── Fallback: basic ape fan plot ──────────────────────────────────────────
  cat("  ggtree not installed — using ape::plot.phylo(fan)\n")
  tryCatch({
    h_px <- max(10, min(20, length(tree$tip.label) * 0.08 + 4))
    pdf(pdf_out, width = 14, height = h_px)
    # colour tips by phylum
    phylums  <- unique(tip_meta$Phylum)
    pal_base <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                  "#a65628","#f781bf","#999999","#66c2a5","#fc8d62")
    phy_cols <- setNames(pal_base[seq_len(length(phylums))], phylums)
    tip_colors <- phy_cols[tip_meta$Phylum[match(tree$tip.label, tip_meta$label)]]
    tip_colors[is.na(tip_colors)] <- "grey50"
    plot.phylo(tree, type = "fan", cex = 0.4, tip.color = tip_colors,
               main = "Circular Phylogenetic Tree")
    legend("bottomleft", legend = names(phy_cols), fill = phy_cols,
           cex = 0.6, title = "Phylum", bty = "n")
    dev.off()
    cat("  ✓ phylo_tree.pdf  (ape fan plot)\n")
  }, error = function(e) cat("  [error]", e$message, "\n"))
}

# ── Export interactive plot data (phylo_tree_plot.csv) ───────────────────────
# Rectangular/phylogram layout via ape::plot.phylo
# Outputs: row_type (edge|tip), x0, y0, x1, y1, label, phylum
plot_csv <- file.path(out_dir, "phylo_tree_plot.csv")
tryCatch({
  n_tips  <- length(tree$tip.label)
  # Render to temp device just to extract coordinates from ape
  tmp_pdf <- tempfile(fileext = ".pdf")
  pdf(tmp_pdf, width = 10, height = max(6, n_tips * 0.25))
  pp <- plot.phylo(tree, type = "phylogram", use.edge.length = TRUE,
                   show.tip.label = FALSE, plot = TRUE)
  dev.off()
  try(file.remove(tmp_pdf), silent = TRUE)

  coords_x <- pp$xx   # length = n_tips + n_internal
  coords_y <- pp$yy

  if (length(coords_x) == 0) stop("no coordinates from plot.phylo")

  # ── Build edge segments (2 L-shaped segments per edge) ──
  edge_list <- vector("list", nrow(tree$edge) * 2)
  ei <- 0L
  for (i in seq_len(nrow(tree$edge))) {
    p_idx <- tree$edge[i, 1]; c_idx <- tree$edge[i, 2]
    px <- coords_x[p_idx]; py <- coords_y[p_idx]
    cx <- coords_x[c_idx]; cy <- coords_y[c_idx]
    ei <- ei + 1L
    edge_list[[ei]] <- data.frame(row_type="edge", x0=px, y0=py, x1=cx, y1=py,
                                   label="", phylum="", stringsAsFactors=FALSE)
    ei <- ei + 1L
    edge_list[[ei]] <- data.frame(row_type="edge", x0=cx, y0=py, x1=cx, y1=cy,
                                   label="", phylum="", stringsAsFactors=FALSE)
  }
  edges_df <- do.call(rbind, edge_list[seq_len(ei)])

  # ── Build tip node rows ──────────────────────────────────
  tip_phyla <- tip_meta$Phylum[match(tree$tip.label, tip_meta$label)]
  tip_phyla[is.na(tip_phyla)] <- "Unknown"
  tips_df <- data.frame(
    row_type = "tip",
    x0       = coords_x[seq_len(n_tips)],
    y0       = coords_y[seq_len(n_tips)],
    x1       = coords_x[seq_len(n_tips)],
    y1       = coords_y[seq_len(n_tips)],
    label    = tree$tip.label,
    phylum   = tip_phyla,
    stringsAsFactors = FALSE
  )

  plot_data <- rbind(edges_df, tips_df)
  write.csv(plot_data, plot_csv, row.names = FALSE)
  cat("  ✓ phylo_tree_plot.csv  (", n_tips, "tips,", nrow(edges_df), "segments)\n")
}, error = function(e) cat("  [warn] phylo_tree_plot.csv skipped:", e$message, "\n"))

cat("── Phylogenetic Tree done ──────────────────────────────────────────\n\n")
