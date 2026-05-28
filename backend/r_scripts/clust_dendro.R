#!/usr/bin/env Rscript
# =============================================================================
#  clust_dendro.R — Hierarchical Clustering Dendrogram + Composition Bars
#
#  Usage: Rscript clust_dendro.R <output_dir>
#
#  Reads:  bray_curtis_distance_matrix.csv  (preferred)
#          OR asv_table.csv                 (fallback — computes Bray-Curtis)
#          taxonomy_genus.csv / taxonomy_family.csv / taxonomy_phylum.csv
#          metadata.csv (SampleID, Group — for label colors)
#  Writes: r_plots/clust_dendro.pdf
#          clust_dendro_segments.csv  (dendrogram line segments)
#          clust_dendro_bars.csv      (sample × taxon abundance, ordered by hclust)
#          clust_dendro_meta.csv      (SampleID, Group, order_index)
# =============================================================================

args    <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── Clust. Dendrogram + Composition ─────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

suppressPackageStartupMessages({
  has_vegan  <- requireNamespace("vegan",    quietly=TRUE)
  has_ggplot <- requireNamespace("ggplot2",  quietly=TRUE)
  has_dplyr  <- requireNamespace("dplyr",    quietly=TRUE)
})

# ── Load distance matrix or compute Bray-Curtis ───────────────────────────────
dm_mat <- NULL
bray_f <- file.path(out_dir, "bray_curtis_distance_matrix.csv")
asv_f  <- file.path(out_dir, "asv_table.csv")

if (file.exists(bray_f)) {
  tryCatch({
    dm      <- read.csv(bray_f, row.names=1, check.names=FALSE)
    dm_mat  <- as.matrix(dm)
    cat("  Loaded Bray-Curtis matrix:", nrow(dm_mat), "samples\n")
  }, error=function(e) cat("  [warn] Could not load Bray-Curtis:", e$message, "\n"))
}

if (is.null(dm_mat) && file.exists(asv_f) && has_vegan) {
  tryCatch({
    library(vegan)
    raw <- read.csv(asv_f, check.names=FALSE)
    seq_col <- which(colnames(raw) == "sequence")
    if (length(seq_col)) raw <- raw[, -seq_col, drop=FALSE]
    mat    <- t(as.matrix(raw)); storage.mode(mat) <- "numeric"
    dm_mat <- as.matrix(vegdist(mat, method="bray"))
    cat("  Computed Bray-Curtis from asv_table.csv:", nrow(dm_mat), "samples\n")
  }, error=function(e) cat("  [warn] Could not compute Bray-Curtis:", e$message, "\n"))
}

if (is.null(dm_mat)) {
  cat("  [skip] No distance matrix available\n"); quit(status=0)
}

sample_names <- rownames(dm_mat)
n_samples    <- length(sample_names)
if (n_samples < 3) {
  cat("  [skip] Need at least 3 samples, got", n_samples, "\n"); quit(status=0)
}

# ── Hierarchical clustering ───────────────────────────────────────────────────
hc <- hclust(as.dist(dm_mat), method="ward.D2")
cat("  Clustering method: Ward.D2 (", n_samples, "samples)\n")

# Extract dendrogram segments for plotting
dend_to_segments <- function(hc) {
  dend    <- as.dendrogram(hc)
  n       <- length(hc$labels)
  heights <- hc$height
  merges  <- hc$merge

  # Use dendrapply to extract segments
  segs <- list()
  leaf_pos <- setNames(seq_len(n), hc$labels[hc$order])

  # Recursive segment extraction
  get_segments <- function(d, pos=NULL) {
    if (is.leaf(d)) return(list(x=attr(d,"nodePOS") %||% 0, y=0))
    kids  <- list(d[[1]], d[[2]])
    ht    <- attr(d, "height")
    child_info <- lapply(kids, get_segments)
    xl <- child_info[[1]]$x; yl <- child_info[[1]]$y
    xr <- child_info[[2]]$x; yr <- child_info[[2]]$y
    xm <- (xl + xr) / 2
    segs[[length(segs)+1]] <<- data.frame(x=xl, xend=xl, y=yl, yend=ht)
    segs[[length(segs)+1]] <<- data.frame(x=xr, xend=xr, y=yr, yend=ht)
    segs[[length(segs)+1]] <<- data.frame(x=xl, xend=xr, y=ht, yend=ht)
    list(x=xm, y=ht)
  }

  # Use ggdendro if available
  if (requireNamespace("ggdendro", quietly=TRUE)) {
    ddata <- ggdendro::dendro_data(dend, type="rectangle")
    segs  <- ddata$segments
    labs  <- ddata$labels
    return(list(segments=segs, labels=labs, order=hc$order))
  }

  # Fallback: manual extraction
  dend_seg <- function(d, x_from, x_to) {
    ht <- attr(d, "height")
    if (is.leaf(d)) {
      return(data.frame(x=double(), xend=double(), y=double(), yend=double()))
    }
    x_mid <- (x_from + x_to) / 2
    n_left  <- length(d[[1]])
    n_right <- length(d[[2]])
    prop_l  <- n_left  / (n_left + n_right)
    x_split <- x_from + (x_to - x_from) * prop_l
    ht_l    <- attr(d[[1]], "height")
    ht_r    <- attr(d[[2]], "height")
    this <- data.frame(
      x    = c(x_from, x_to,    x_from, x_to),
      xend = c(x_from, x_to,    x_to,   x_to),
      y    = c(ht_l,   ht_r,    ht,     ht),
      yend = c(ht,     ht,      ht,     ht)
    )
    rbind(this,
          dend_seg(d[[1]], x_from, x_split),
          dend_seg(d[[2]], x_split, x_to))
  }
  # Simple ordered segments from hclust
  merge_heights <- hc$height
  leaf_order    <- hc$order
  pos <- setNames(seq_along(leaf_order), hc$labels[leaf_order])
  seg_list <- list()
  node_pos <- rep(0, n)
  node_pos[seq_len(n)] <- seq_len(n)  # leaf positions
  node_ht  <- c(rep(0, n), merge_heights)
  cluster_pos <- pos  # leaf positions

  node_x <- rep(NA_real_, n - 1)
  for (i in seq_len(nrow(merges))) {
    l <- merges[i, 1]; r <- merges[i, 2]
    xl <- if (l < 0) cluster_pos[hc$labels[-l]] else node_x[l]
    xr <- if (r < 0) cluster_pos[hc$labels[-r]] else node_x[r]
    yl <- if (l < 0) 0 else heights[l]
    yr <- if (r < 0) 0 else heights[r]
    ht <- heights[i]
    node_x[i] <- (xl + xr) / 2
    seg_list[[length(seg_list)+1]] <- data.frame(x=xl, xend=xl, y=yl, yend=ht)
    seg_list[[length(seg_list)+1]] <- data.frame(x=xr, xend=xr, y=yr, yend=ht)
    seg_list[[length(seg_list)+1]] <- data.frame(x=xl, xend=xr, y=ht, yend=ht)
  }
  segs <- do.call(rbind, seg_list)

  label_df <- data.frame(
    x     = as.numeric(pos[hc$labels[hc$order]]),
    label = hc$labels[hc$order],
    stringsAsFactors=FALSE)

  list(segments=segs, labels=label_df, order=hc$order)
}

dend_data <- tryCatch(dend_to_segments(hc), error=function(e) NULL)
if (!is.null(dend_data) && !is.null(dend_data$segments)) {
  write.csv(dend_data$segments, file.path(out_dir, "clust_dendro_segments.csv"), row.names=FALSE)
  cat("  ✓ clust_dendro_segments.csv\n")
}

# ── Load metadata ─────────────────────────────────────────────────────────────
meta_map <- NULL
meta_f   <- file.path(out_dir, "metadata.csv")
if (file.exists(meta_f)) tryCatch({
  m <- read.csv(meta_f, stringsAsFactors=FALSE)
  sid_col <- grep("sampleid|sample_id|sample", colnames(m), ignore.case=TRUE, value=TRUE)[1]
  grp_col <- grep("^group$",                   colnames(m), ignore.case=TRUE, value=TRUE)[1]
  if (!is.na(sid_col) && !is.na(grp_col))
    meta_map <- setNames(m[[grp_col]], m[[sid_col]])
}, error=function(e) NULL)

# ── Load taxonomy (genus > family > phylum preference) ─────────────────────────
tax_df <- NULL
tax_level <- "Unknown"
for (lvl in c("genus","family","order","phylum")) {
  tf <- file.path(out_dir, paste0("taxonomy_", lvl, ".csv"))
  if (file.exists(tf)) {
    tryCatch({
      tax_df    <- read.csv(tf, check.names=FALSE, stringsAsFactors=FALSE)
      tax_level <- lvl
      cat("  Taxonomy:", lvl, "(", ncol(tax_df)-1, "taxa)\n")
      break
    }, error=function(e) NULL)
  }
}

# ── Build composition bar data in hclust order ────────────────────────────────
ordered_samples <- sample_names[hc$order]

if (!is.null(tax_df)) {
  # Detect orientation: more rows = taxa-as-rows (ONT), else taxa-as-cols (regular)
  n_rows_tax <- nrow(tax_df) - 1
  n_cols_tax <- ncol(tax_df) - 1
  taxa_as_cols <- n_cols_tax >= n_rows_tax

  if (taxa_as_cols) {
    samples  <- tax_df[[1]]
    taxa_nm  <- colnames(tax_df)[-1]
    abu_mat  <- as.matrix(tax_df[, -1]); storage.mode(abu_mat) <- "numeric"
    rownames(abu_mat) <- samples
  } else {
    taxa_nm  <- tax_df[[1]]
    abu_mat  <- t(as.matrix(tax_df[, -1])); storage.mode(abu_mat) <- "numeric"
    rownames(abu_mat) <- colnames(tax_df)[-1]
    colnames(abu_mat) <- taxa_nm
    taxa_nm  <- colnames(abu_mat)
  }

  # Top 20 taxa by total abundance
  col_totals <- colSums(abu_mat, na.rm=TRUE)
  top_taxa   <- names(sort(col_totals, decreasing=TRUE))[seq_len(min(20, ncol(abu_mat)))]

  # Relative abundance
  row_sums  <- rowSums(abu_mat, na.rm=TRUE)
  rel_mat   <- sweep(abu_mat, 1, pmax(row_sums, 1e-12), "/") * 100

  # Filter to ordered samples only
  common_samples <- intersect(ordered_samples, rownames(rel_mat))
  if (length(common_samples) < 2) {
    # Try case-insensitive match
    common_samples <- ordered_samples[
      tolower(ordered_samples) %in% tolower(rownames(rel_mat))]
  }

  if (length(common_samples) >= 2) {
    bar_mat <- rel_mat[common_samples, top_taxa, drop=FALSE]
    bar_df  <- data.frame(Sample=common_samples, bar_mat, check.names=FALSE,
                           stringsAsFactors=FALSE)
    write.csv(bar_df, file.path(out_dir, "clust_dendro_bars.csv"), row.names=FALSE)
    cat("  ✓ clust_dendro_bars.csv  (", nrow(bar_df), "samples ×",
        length(top_taxa), "taxa)\n")
  }
}

# ── Metadata order CSV ────────────────────────────────────────────────────────
meta_out <- data.frame(
  SampleID    = ordered_samples,
  Group       = if (!is.null(meta_map)) meta_map[ordered_samples] else NA,
  order_index = seq_along(ordered_samples),
  stringsAsFactors=FALSE)
meta_out$Group[is.na(meta_out$Group)] <- "Unknown"
write.csv(meta_out, file.path(out_dir, "clust_dendro_meta.csv"), row.names=FALSE)
cat("  ✓ clust_dendro_meta.csv\n")

# ── PDF output ────────────────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings=FALSE, recursive=TRUE)
pdf_out   <- file.path(plots_dir, "clust_dendro.pdf")

if (has_ggplot && !is.null(tax_df)) {
  suppressPackageStartupMessages({ library(ggplot2) })
  tryCatch({
    n_samp <- length(ordered_samples)
    h_pdf  <- max(8, min(20, n_samp * 0.35 + 4))

    # Group color palette
    grp_pal <- c("#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6",
                 "#06b6d4","#f97316","#84cc16","#ec4899","#14b8a6")
    uniq_grps <- unique(meta_out$Group)
    grp_cols  <- setNames(grp_pal[seq_along(uniq_grps)], uniq_grps)

    # Taxon color palette
    tax_pal <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854",
                 "#ffd92f","#e5c494","#b3b3b3","#1f78b4","#33a02c",
                 "#e31a1c","#ff7f00","#6a3d9a","#b15928","#a6cee3",
                 "#b2df8a","#fb9a99","#fdbf6f","#cab2d6","#ffff99")

    # Get dendrogram segments (from ggdendro or fallback)
    dend_plot  <- as.dendrogram(hc)
    seg_df     <- NULL
    if (requireNamespace("ggdendro", quietly=TRUE)) {
      ddata  <- ggdendro::dendro_data(dend_plot, type="rectangle")
      seg_df <- ddata$segments
      lab_df <- data.frame(x=seq_along(hc$order),
                            label=hc$labels[hc$order], stringsAsFactors=FALSE)
    }

    # Build bar data
    bar_src <- file.path(out_dir, "clust_dendro_bars.csv")
    if (!file.exists(bar_src) || is.null(seg_df)) {
      # Simple base-R dendrogram plot
      pdf(pdf_out, width=14, height=h_pdf)
      par(mfrow=c(1,2), mar=c(2,1,2,0))
      plot(hc, hang=-1, main="Sample Similarity", xlab="", sub="",
           labels=hc$labels[hc$order], cex=0.7)
      if (!is.null(tax_df)) {
        bar_mat <- rel_mat[ordered_samples, top_taxa, drop=FALSE]
        barplot(t(bar_mat), horiz=TRUE, las=1, col=tax_pal[seq_along(top_taxa)],
                main=paste0("Taxonomic Composition (", tax_level, ")"),
                xlab="Relative Abundance (%)", cex.names=0.5)
        legend("topright", legend=rev(top_taxa), fill=rev(tax_pal[seq_along(top_taxa)]),
               cex=0.5, bty="n")
      }
      dev.off()
      cat("  ✓ clust_dendro.pdf  (base-R)\n")
    } else {
      bar_dat <- read.csv(bar_src, check.names=FALSE, stringsAsFactors=FALSE)
      bar_long <- data.frame()
      for (tc in colnames(bar_dat)[-1]) {
        tmp <- data.frame(Sample=bar_dat$Sample, Taxon=tc,
                          Abundance=as.numeric(bar_dat[[tc]]), stringsAsFactors=FALSE)
        bar_long <- rbind(bar_long, tmp)
      }
      bar_long$Sample <- factor(bar_long$Sample, levels=rev(ordered_samples))
      bar_long$Taxon  <- factor(bar_long$Taxon,  levels=rev(top_taxa))

      n_taxa_used <- length(unique(bar_long$Taxon))
      bar_cols    <- setNames(tax_pal[seq_len(n_taxa_used)], rev(top_taxa)[seq_len(n_taxa_used)])

      sample_grp_df <- data.frame(
        Sample   = meta_out$SampleID,
        Group    = meta_out$Group,
        stringsAsFactors=FALSE)
      bar_long <- merge(bar_long, sample_grp_df, by="Sample", all.x=TRUE)
      bar_long$Group[is.na(bar_long$Group)] <- "Unknown"

      # Dendrogram plot (left)
      p_dend <- ggplot(seg_df, aes(x=y, xend=yend, y=x, yend=xend)) +
        geom_segment() +
        scale_y_continuous(breaks=seq_along(ordered_samples),
                           labels=rev(ordered_samples),
                           expand=c(0.01,0.01)) +
        scale_x_reverse() +
        labs(x="Distance (Bray-Curtis)", y=NULL, title="Similarity") +
        theme_classic(base_size=9) +
        theme(axis.text.y=element_text(
          color=grp_cols[rev(meta_out$Group)[seq_along(ordered_samples)]],
          size=8, hjust=1),
          axis.ticks.y=element_blank(),
          plot.title=element_text(size=10, face="bold", hjust=0.5),
          plot.margin=margin(2,0,2,4,"mm"))

      # Stacked bar plot (right)
      p_bar <- ggplot(bar_long, aes(x=Abundance, y=Sample, fill=Taxon)) +
        geom_bar(stat="identity", orientation="y", width=0.85) +
        scale_fill_manual(values=bar_cols, name=tools::toTitleCase(tax_level)) +
        scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
        labs(x="Relative Abundance (%)", y=NULL,
             title=paste0("Taxonomic Composition (", tools::toTitleCase(tax_level), ")")) +
        theme_classic(base_size=9) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              legend.text=element_text(size=7), legend.key.size=unit(0.35,"cm"),
              plot.title=element_text(size=10, face="bold", hjust=0.5),
              plot.margin=margin(2,4,2,0,"mm"))

      # Group legend bar (below dendrogram)
      if (requireNamespace("patchwork", quietly=TRUE)) {
        suppressPackageStartupMessages(library(patchwork))
        combined <- p_dend + p_bar + plot_layout(widths=c(1,2))
        ggsave(pdf_out, combined, width=16, height=h_pdf, limitsize=FALSE)
      } else {
        pdf(pdf_out, width=16, height=h_pdf)
        print(p_dend); print(p_bar)
        dev.off()
      }
      cat("  ✓ clust_dendro.pdf  (ggplot2)\n")
    }
  }, error=function(e) {
    cat("  [error] PDF:", e$message, "\n")
    tryCatch({
      pdf(pdf_out, width=14, height=10)
      plot(hc, hang=-1, main="Sample Similarity Clustering", cex=0.7)
      dev.off()
      cat("  ✓ clust_dendro.pdf  (hclust fallback)\n")
    }, error=function(e2) cat("  [error] fallback:", e2$message, "\n"))
  })
} else {
  # No ggplot2 — base plot
  tryCatch({
    pdf(pdf_out, width=14, height=max(8, n_samples * 0.3 + 2))
    plot(hc, hang=-1, main="Sample Similarity Clustering", cex=0.7)
    dev.off()
    cat("  ✓ clust_dendro.pdf  (base-R)\n")
  }, error=function(e) cat("  [error]", e$message, "\n"))
}

cat("── Clust. Dendrogram done ──────────────────────────────────────────\n\n")
