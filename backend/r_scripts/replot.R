# ============================================================
#  replot.R — Regenerate taxonomy plots with custom colors
#  Called by the backend when user saves colors in the UI.
#  Reads existing asv_table.csv + taxonomy_table.csv from the
#  output directory — does NOT re-run DADA2.
# ============================================================

tryCatch({
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset="")
  if (nchar(r_libs_user) > 0 && dir.exists(r_libs_user))
    .libPaths(unique(c(r_libs_user, .libPaths())))
  home_dir <- Sys.getenv("HOME", unset=path.expand("~"))
  r_root   <- file.path(home_dir, "R")
  if (dir.exists(r_root)) {
    vdirs <- list.dirs(r_root, recursive=TRUE, full.names=TRUE)
    vdirs <- vdirs[grepl("/\\d+\\.\\d+$", vdirs)]
    if (length(vdirs) > 0) .libPaths(unique(c(vdirs, .libPaths())))
  }
  simple_lib <- file.path(home_dir, "R", "library")
  if (dir.exists(simple_lib)) .libPaths(unique(c(simple_lib, .libPaths())))
}, error=function(e) NULL)

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--output",     type="character", help="Job output directory"),
  make_option("--colorsFile", type="character", default="",
              help="Path to custom_colors.json (optional)")
)
opt <- parse_args(OptionParser(option_list=option_list))

cat("=== replot.R starting ===\n")
cat("Output dir  :", opt$output, "\n")
cat("Colors file :", if (nchar(opt$colorsFile)>0) opt$colorsFile else "(none)", "\n\n")

# ── Load saved results ───────────────────────────────────────
asv_file <- file.path(opt$output, "asv_table.csv")
tax_file <- file.path(opt$output, "taxonomy_table.csv")

if (!file.exists(asv_file)) stop("asv_table.csv not found in: ", opt$output)

# Reconstruct seqtab_nochim (samples × ASVs)
asv_df  <- read.csv(asv_file, check.names=FALSE, stringsAsFactors=FALSE)
seq_col <- asv_df[["sequence"]]
cnt_mat <- asv_df[, colnames(asv_df) != "sequence", drop=FALSE]
seqtab_nochim         <- t(as.matrix(cnt_mat))
mode(seqtab_nochim)   <- "numeric"
colnames(seqtab_nochim) <- seq_col
sample_names          <- rownames(seqtab_nochim)
n_samp                <- nrow(seqtab_nochim)
cat("Loaded", n_samp, "samples,", ncol(seqtab_nochim), "ASVs\n")

# Taxonomy matrix (NULL if file missing)
tax <- NULL
if (file.exists(tax_file)) {
  tax_df <- read.csv(tax_file, check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
  tax    <- as.matrix(tax_df)
  cat("Loaded taxonomy:", nrow(tax), "ASVs,", ncol(tax), "ranks\n")
}

# ── Load custom colors ────────────────────────────────────────
# Format: { "Phylum": {"Bacillota":"#ff0000",...}, "Class":{...}, ... }
custom_colors <- list()
if (nchar(opt$colorsFile) > 0 && file.exists(opt$colorsFile)) {
  tryCatch({
    custom_colors <- jsonlite::fromJSON(opt$colorsFile, simplifyVector=FALSE)
    cat("Custom colors loaded for levels:", paste(names(custom_colors), collapse=", "), "\n")
  }, error=function(e) cat("  Could not load colors:", e$message, "\n"))
}

# ── Palettes ──────────────────────────────────────────────────
palette20 <- c("#e74c3c","#e67e22","#f1c40f","#2ecc71","#1abc9c",
               "#3498db","#9b59b6","#e91e63","#00bcd4","#8bc34a",
               "#ff5722","#607d8b","#795548","#ffc107","#03a9f4",
               "#4caf50","#673ab7","#ff9800","#009688","#9e9e9e")

make_pal <- function(n) {
  base <- c(palette20,
            "#a52a2a","#5f9ea0","#d2691e","#6495ed","#dc143c",
            "#00ced1","#ff8c00","#9400d3","#32cd32","#ff1493",
            "#1e90ff","#ffd700","#adff2f","#ff6347","#40e0d0",
            "#ee82ee","#f5deb3","#00ff7f","#87ceeb","#dda0dd")
  rep(base, length.out=n)[seq_len(n)]
}

samp_cols <- setNames(make_pal(n_samp), sample_names)

# ── Load optional packages ────────────────────────────────────
has_ggplot2  <- requireNamespace("ggplot2",  quietly=TRUE)
has_reshape2 <- requireNamespace("reshape2", quietly=TRUE)
has_pheatmap <- requireNamespace("pheatmap", quietly=TRUE)
if (has_ggplot2)  suppressPackageStartupMessages(library(ggplot2))
if (has_reshape2) suppressPackageStartupMessages(library(reshape2))
if (has_pheatmap) suppressPackageStartupMessages(library(pheatmap))

# ── Load metadata (if saved) ──────────────────────────────────
meta_df   <- NULL
meta_cols <- list()
meta_pals <- list()
meta_file <- file.path(opt$output, "metadata.csv")
if (file.exists(meta_file)) {
  tryCatch({
    meta_df <- read.csv(meta_file, stringsAsFactors=FALSE)
    sid_col <- names(meta_df)[tolower(names(meta_df)) %in%
                 c("sampleid","sample_id","sample","#sampleid")][1]
    if (!is.na(sid_col)) rownames(meta_df) <- trimws(as.character(meta_df[[sid_col]]))
    else rownames(meta_df) <- sample_names
    meta_df <- meta_df[sample_names, , drop=FALSE]
    skip_pat <- c("sampleid","sample_id","sample","#sampleid")
    for (col in names(meta_df)) {
      if (tolower(col) %in% skip_pat) next
      vals   <- trimws(as.character(meta_df[[col]]))
      vals[vals == "NA" | is.na(vals)] <- ""
      n_uniq <- length(unique(vals[nchar(vals) > 0]))
      if (n_uniq >= 2 && n_uniq <= 20 && n_uniq < n_samp) {
        gv <- setNames(vals, sample_names)
        lvs <- unique(gv[nchar(gv) > 0])
        offset <- (length(meta_cols) * 7) %% length(palette20)
        gp <- setNames(palette20[(seq_along(lvs)+offset-1) %% length(palette20)+1], lvs)
        meta_cols[[col]] <- gv
        meta_pals[[col]] <- gp
      }
    }
    cat("Metadata loaded:", nrow(meta_df), "samples,",
        length(meta_cols), "grouping column(s)\n")
  }, error=function(e) cat("  Could not load metadata:", e$message, "\n"))
}

# ── Helper: apply custom colors for a given taxonomy level ────
# Returns a named colour vector (taxon → colour), merging custom
# over the default palette. Unassigned taxa get default colours.
level_colours <- function(taxa_names, level) {
  n     <- length(taxa_names)
  cols  <- setNames(make_pal(n), taxa_names)   # defaults
  lvl_custom <- custom_colors[[level]]          # user-defined (may be NULL)
  if (!is.null(lvl_custom) && length(lvl_custom) > 0) {
    for (nm in names(lvl_custom)) {
      if (nm %in% names(cols)) cols[[nm]] <- lvl_custom[[nm]]
    }
  }
  cols
}

# ── Helper: build sample × taxon relative-abundance matrix ────
make_tax_mat <- function(level, top_n=50) {
  if (is.null(tax)) return(NULL)
  col_idx <- which(colnames(tax) == level)
  if (length(col_idx) == 0) return(NULL)
  taxa_vec <- tax[, col_idx]
  taxa_vec[is.na(taxa_vec) | taxa_vec == ""] <- "Unclassified"
  taxa_vec <- sub("^[a-z]__", "", taxa_vec)
  taxa_uniq <- unique(taxa_vec)
  mat <- sapply(taxa_uniq, function(t)
    rowSums(seqtab_nochim[, taxa_vec == t, drop=FALSE]))
  if (is.null(dim(mat)))
    mat <- matrix(mat, nrow=1, dimnames=list(sample_names, taxa_uniq))
  mat <- mat[, order(colSums(mat), decreasing=TRUE), drop=FALSE]
  if (ncol(mat) > top_n) {
    other <- rowSums(mat[,(top_n+1):ncol(mat), drop=FALSE])
    mat   <- cbind(mat[,1:top_n, drop=FALSE], Other=other)
  }
  pct <- sweep(mat, 1, rowSums(mat), "/") * 100
  pct[is.nan(pct)] <- 0
  pct
}

# ── Helper: stacked bar with custom colours ────────────────────
tax_stacked_bar <- function(pct_mat, title_str, outfile,
                             level="", x_labels=NULL, x_title="Sample") {
  if (is.null(pct_mat)) return(invisible(NULL))
  taxa_names <- colnames(pct_mat)
  cols_t     <- level_colours(taxa_names, level)
  if (is.null(x_labels)) x_labels <- rownames(pct_mat)
  n_x <- nrow(pct_mat)

  if (has_ggplot2 && has_reshape2) {
    df <- as.data.frame(pct_mat)
    df[[x_title]] <- factor(x_labels, levels=x_labels)
    long <- reshape2::melt(df, id.vars=x_title, variable.name="Taxon", value.name="Pct")
    long$Taxon <- factor(long$Taxon, levels=taxa_names)
    p <- ggplot2::ggplot(long, ggplot2::aes(x=.data[[x_title]], y=Pct, fill=Taxon)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::scale_fill_manual(values=cols_t) +
      ggplot2::labs(title=title_str, x=x_title, y="Relative Abundance (%)") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                     legend.text=ggplot2::element_text(size=7),
                     legend.key.size=ggplot2::unit(0.35,"cm"))
    w <- max(7, n_x * 1.1 + 4)
    h <- max(6, ceiling(length(taxa_names)/4)*0.4 + 5)
    ggplot2::ggsave(outfile, p, width=w, height=h, device="pdf", limitsize=FALSE)
  } else {
    pdf(outfile, width=max(9, n_x*1.1+4), height=6)
    par(mar=c(8,5,3,max(8, max(nchar(taxa_names))*0.5)), xpd=TRUE)
    barplot(t(pct_mat), beside=FALSE, col=cols_t, border=NA,
            main=title_str, ylab="Relative Abundance (%)",
            names.arg=x_labels, las=2, cex.names=0.8, ylim=c(0,100))
    legend(par("usr")[2]*1.02, par("usr")[4], legend=taxa_names,
           fill=cols_t, bty="n", cex=0.65, xpd=TRUE)
    dev.off()
  }
  cat("  Regenerated:", basename(outfile), "\n")
}

# ═══════════════════════════════════════════════════════════════
#  REGENERATE — taxonomy stacked bars (per-sample)
# ═══════════════════════════════════════════════════════════════
cat("\nRegenerating taxonomy plots...\n")
if (!is.null(tax)) {
  for (lvl in c("Phylum","Class","Order","Family","Genus","Species")) {
    tryCatch({
      m <- make_tax_mat(lvl, top_n=50)
      if (is.null(m)) next
      out_f <- file.path(opt$output, paste0("taxonomy_", tolower(lvl), ".pdf"))
      tax_stacked_bar(m,
                      title_str = paste("Relative Abundance — Top 50", lvl),
                      outfile   = out_f,
                      level     = lvl)
    }, error=function(e) cat("  [skip] taxonomy", lvl, ":", e$message, "\n"))
  }
}

# ═══════════════════════════════════════════════════════════════
#  REGENERATE — group mean stacked bars
# ═══════════════════════════════════════════════════════════════
if (length(meta_cols) > 0 && !is.null(tax)) {
  cat("Regenerating group mean plots...\n")
  for (col_name in names(meta_cols)) {
    gv   <- meta_cols[[col_name]]
    grps <- unique(gv[nchar(gv) > 0])
    if (length(grps) < 2) next
    safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)

    for (lvl in c("Phylum","Class","Order","Family","Genus")) {
      tryCatch({
        col_idx <- which(colnames(tax) == lvl)
        if (length(col_idx) == 0) next
        taxa_vec  <- sub("^[a-z]__", "", tax[, col_idx])
        taxa_vec[is.na(taxa_vec) | taxa_vec == ""] <- "Unclassified"
        rel_ab    <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/") * 100
        taxa_uniq <- unique(taxa_vec)
        mat <- sapply(taxa_uniq, function(t)
          rowSums(rel_ab[, taxa_vec == t, drop=FALSE]))
        if (is.null(dim(mat)))
          mat <- matrix(mat, nrow=1, dimnames=list(sample_names, taxa_uniq))
        grp_mat <- do.call(rbind, lapply(grps, function(g) {
          idx <- which(gv == g)
          if (length(idx) == 0) return(rep(0, ncol(mat)))
          if (length(idx) == 1) return(mat[idx,])
          colMeans(mat[idx, , drop=FALSE])
        }))
        rownames(grp_mat) <- grps
        top_n  <- min(50, ncol(grp_mat))
        ord    <- order(colMeans(grp_mat), decreasing=TRUE)
        if (ncol(grp_mat) > top_n) {
          other   <- rowSums(grp_mat[, ord[(top_n+1):ncol(grp_mat)], drop=FALSE])
          grp_mat <- cbind(grp_mat[, ord[1:top_n], drop=FALSE], Other=other)
        } else {
          grp_mat <- grp_mat[, ord, drop=FALSE]
        }
        out_f <- file.path(opt$output,
                           sprintf("group_mean_%s_%s.pdf", safe_col, tolower(lvl)))
        tax_stacked_bar(grp_mat,
                        title_str = sprintf("Mean Rel. Abundance by %s — %s", col_name, lvl),
                        outfile   = out_f,
                        level     = lvl,
                        x_title   = col_name)
      }, error=function(e) cat("  [skip] group mean", col_name, lvl, ":", e$message, "\n"))
    }
  }
}

# ═══════════════════════════════════════════════════════════════
#  REGENERATE — taxonomy heatmaps
# ═══════════════════════════════════════════════════════════════
if (!is.null(tax) && has_pheatmap) {
  cat("Regenerating heatmaps...\n")
  ann_col_hm    <- NULL
  ann_colors_hm <- list()
  if (length(meta_cols) > 0) {
    ann_col_hm <- data.frame(
      lapply(meta_cols, function(gv) factor(gv[sample_names])),
      row.names=sample_names, stringsAsFactors=FALSE
    )
    names(ann_col_hm) <- names(meta_cols)
    for (cn in names(meta_cols)) {
      gp  <- meta_pals[[cn]]
      lvs <- levels(ann_col_hm[[cn]])
      ann_colors_hm[[cn]] <- gp[lvs[lvs %in% names(gp)]]
    }
  }

  for (hm_lvl in c("Genus","Family","Phylum")) {
    tryCatch({
      col_idx_hm <- which(colnames(tax) == hm_lvl)
      if (length(col_idx_hm) == 0) next
      tv <- sub("^[a-z]__", "", tax[, col_idx_hm])
      tv[is.na(tv) | tv == ""] <- "Unclassified"
      rel_hm <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/") * 100
      tu     <- unique(tv)
      mat_hm <- sapply(tu, function(t) rowSums(rel_hm[, tv == t, drop=FALSE]))
      if (is.null(dim(mat_hm)))
        mat_hm <- matrix(mat_hm, nrow=1, dimnames=list(sample_names, tu))
      top_n_hm <- if (hm_lvl=="Phylum") min(20,ncol(mat_hm)) else min(30,ncol(mat_hm))
      top_taxa_hm <- names(sort(colMeans(mat_hm), decreasing=TRUE))[1:top_n_hm]
      heat_mat_hm <- t(mat_hm[, top_taxa_hm, drop=FALSE])

      # Apply custom row colours (per-taxon) via annotation_row colour override
      # pheatmap doesn't directly support per-row custom fill for data cells,
      # so we use the standard blue palette; custom colours affect stacked bars only.
      hm_colors <- switch(hm_lvl,
        Genus  = colorRampPalette(c("#f0f4ff","#3b82f6","#1e1b4b"))(100),
        Family = colorRampPalette(c("#fff7ed","#f97316","#431407"))(100),
        Phylum = colorRampPalette(c("#f0fff4","#22c55e","#14532d"))(100)
      )
      hm_file <- file.path(opt$output,
                           sprintf("taxonomy_heatmap_%s.pdf", tolower(hm_lvl)))
      pheatmap::pheatmap(
        heat_mat_hm,
        annotation_col    = ann_col_hm,
        annotation_colors = if (length(ann_colors_hm)>0) ann_colors_hm else NULL,
        color             = hm_colors,
        scale             = "row",
        main              = sprintf("Top %d %s — Heatmap", top_n_hm, hm_lvl),
        fontsize_row      = max(5, min(9, 200/top_n_hm)),
        fontsize_col      = 8,
        filename          = hm_file,
        width             = max(8, n_samp*0.6+4),
        height            = max(8, top_n_hm*0.35+3)
      )
      cat("  Regenerated:", basename(hm_file), "\n")
    }, error=function(e) cat("  [skip] heatmap", hm_lvl, ":", e$message, "\n"))
  }
}

# ── Update taxonomy_summary.json with new colour info ─────────
tryCatch({
  summary_file <- file.path(opt$output, "taxonomy_summary.json")
  if (file.exists(summary_file)) {
    tax_summary <- jsonlite::fromJSON(summary_file, simplifyVector=FALSE)
    # Embed custom colors into summary so frontend can read them back
    tax_summary[["custom_colors"]] <- custom_colors
    write(jsonlite::toJSON(tax_summary, auto_unbox=TRUE), summary_file)
  }
}, error=function(e) cat("  [skip] summary update:", e$message, "\n"))

cat("\n=== replot.R complete ===\n")
