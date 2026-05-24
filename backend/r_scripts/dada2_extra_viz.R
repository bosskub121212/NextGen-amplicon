#!/usr/bin/env Rscript
# =============================================================================
#  dada2_extra_viz.R — Compute interactive chart CSVs from DADA2 root outputs
#  Sourced by both dada2_pipeline.R (new jobs) and viz_pipeline.R (Re-run Viz)
#
#  Input:  bray_curtis_distance_matrix.csv, asv_table.csv, taxonomy_table.csv
#  Output: pca_scores.csv, pca_scree.csv, nmds_bray.csv,
#          jaccard_heatmap.csv, nmds_jaccard.csv, rarefaction.csv, asv_lengths.csv
# =============================================================================

dada2_extra_viz <- function(out_dir) {
  bray_file <- file.path(out_dir, "bray_curtis_distance_matrix.csv")
  asv_file  <- file.path(out_dir, "asv_table.csv")
  tax_file  <- file.path(out_dir, "taxonomy_table.csv")

  # Check: at least one usable input must exist
  have_any <- file.exists(bray_file) || file.exists(asv_file) || file.exists(tax_file)
  if (!have_any) return(invisible(NULL))

  cat("\n── DADA2 Extra Viz ─────────────────────────────────────────────\n")
  cat("  Directory:", out_dir, "\n")

  have_vegan <- suppressPackageStartupMessages(
    tryCatch({ library(vegan); TRUE }, error=function(e) {
      cat("  [skip beta/rarefaction] vegan not available\n"); FALSE }))

  # ── 1. PCoA from Bray-Curtis ────────────────────────────────────────────────
  if (!file.exists(bray_file)) {
    cat("  (bray_curtis_distance_matrix.csv missing — skip PCoA + NMDS Bray)\n")
  } else tryCatch({
    pca_out <- file.path(out_dir, "pca_scores.csv")
    if (file.exists(pca_out)) {
      cat("  (pca_scores.csv exists — skip)\n")
    } else {
      dm     <- read.csv(bray_file, row.names=1, check.names=FALSE)
      dm_mat <- as.matrix(dm)
      n_s    <- nrow(dm_mat)
      k      <- min(n_s - 1, 10)
      pcoa   <- cmdscale(dm_mat, k=k, eig=TRUE)
      eig_p  <- pcoa$eig; eig_p[eig_p < 0] <- 0
      var_p  <- if (sum(eig_p) > 0) pcoa$eig / sum(eig_p) * 100 else rep(0, length(eig_p))
      var_p[var_p < 0] <- 0
      pc_df  <- data.frame(
        Sample  = rownames(dm_mat),
        PC1     = pcoa$points[, 1],
        PC2     = pcoa$points[, 2],
        PC3     = if (k >= 3) pcoa$points[, 3] else 0,
        PC1_var = round(var_p[1], 1),
        PC2_var = round(var_p[2], 1),
        PC3_var = if (length(var_p) >= 3) round(var_p[3], 1) else 0,
        stringsAsFactors = FALSE)
      write.csv(pc_df, pca_out, row.names=FALSE)
      cat("  ✓ pca_scores.csv\n")

      n_sc  <- min(10, sum(eig_p > 0))
      sc_df <- data.frame(PC=seq_len(n_sc), Variance=round(var_p[seq_len(n_sc)], 2))
      write.csv(sc_df, file.path(out_dir, "pca_scree.csv"), row.names=FALSE)
      cat("  ✓ pca_scree.csv\n")
    }
  }, error=function(e) cat("  [skip] PCoA:", e$message, "\n"))

  # ── 2. NMDS Bray-Curtis ─────────────────────────────────────────────────────
  if (have_vegan && file.exists(bray_file)) tryCatch({
    nmds_bray_out <- file.path(out_dir, "nmds_bray.csv")
    if (file.exists(nmds_bray_out)) {
      cat("  (nmds_bray.csv exists — skip)\n")
    } else {
      dm     <- read.csv(bray_file, row.names=1, check.names=FALSE)
      dm_mat <- as.matrix(dm)
      set.seed(42)
      nm     <- metaMDS(dm_mat, distance="bray", k=2, try=20, trymax=50, trace=FALSE)
      nm_df  <- data.frame(Sample=rownames(dm_mat), NMDS1=nm$points[,1],
                            NMDS2=nm$points[,2], Stress=nm$stress,
                            stringsAsFactors=FALSE)
      write.csv(nm_df, nmds_bray_out, row.names=FALSE)
      cat("  ✓ nmds_bray.csv  (stress =", round(nm$stress, 4), ")\n")
    }
  }, error=function(e) cat("  [skip] NMDS Bray:", e$message, "\n"))

  # ── 3. Jaccard heatmap + NMDS Jaccard ───────────────────────────────────────
  asv_t_cached <- NULL
  if (have_vegan && file.exists(asv_file)) tryCatch({
    asv_raw      <- read.csv(asv_file, check.names=FALSE)
    asv_t_cached <- t(as.matrix(asv_raw))
    rownames(asv_t_cached) <- colnames(asv_raw)

    jac_out <- file.path(out_dir, "jaccard_heatmap.csv")
    if (file.exists(jac_out)) {
      cat("  (jaccard_heatmap.csv exists — skip)\n")
    } else {
      jac_mat <- as.matrix(vegdist(asv_t_cached, method="jaccard", binary=FALSE))
      jac_df  <- as.data.frame(jac_mat)
      write.csv(jac_df, jac_out)
      cat("  ✓ jaccard_heatmap.csv\n")

      set.seed(42)
      nm_j   <- metaMDS(jac_mat, distance="jaccard", k=2, try=20, trymax=50, trace=FALSE)
      nmj_df <- data.frame(Sample=rownames(jac_mat), NMDS1=nm_j$points[,1],
                            NMDS2=nm_j$points[,2], Stress=nm_j$stress,
                            stringsAsFactors=FALSE)
      write.csv(nmj_df, file.path(out_dir, "nmds_jaccard.csv"), row.names=FALSE)
      cat("  ✓ nmds_jaccard.csv  (stress =", round(nm_j$stress, 4), ")\n")
    }
  }, error=function(e) cat("  [skip] Jaccard:", e$message, "\n"))

  # ── 4. Rarefaction curves (wide format: Depth, Sample1, Sample2, ...) ───────
  if (have_vegan && file.exists(asv_file)) tryCatch({
    rar_out <- file.path(out_dir, "rarefaction.csv")
    if (file.exists(rar_out)) {
      cat("  (rarefaction.csv exists — skip)\n")
    } else {
      if (is.null(asv_t_cached)) {
        asv_raw      <- read.csv(asv_file, check.names=FALSE)
        asv_t_cached <- t(as.matrix(asv_raw))
        rownames(asv_t_cached) <- colnames(asv_raw)
      }
      row_sums  <- rowSums(asv_t_cached)
      min_reads <- min(row_sums)
      n_steps   <- min(25, min_reads)
      depths    <- unique(sort(c(
        round(seq(max(100, min_reads / n_steps), min_reads, length.out=n_steps)),
        as.integer(min_reads))))
      depths <- depths[depths >= 1 & depths <= min_reads]

      wide_df <- data.frame(Depth=depths)
      for (i in seq_len(nrow(asv_t_cached))) {
        sname <- rownames(asv_t_cached)[i]
        rich  <- suppressWarnings(rarefy(asv_t_cached[i, , drop=FALSE], sample=depths))
        wide_df[[sname]] <- round(as.numeric(rich), 1)
      }
      write.csv(wide_df, rar_out, row.names=FALSE)
      cat("  ✓ rarefaction.csv  (", nrow(asv_t_cached), "samples,",
          length(depths), "depths, min =", min_reads, "reads)\n")
    }
  }, error=function(e) cat("  [skip] Rarefaction:", e$message, "\n"))

  # ── 5. ASV length distribution ───────────────────────────────────────────────
  if (file.exists(tax_file)) tryCatch({
    len_out <- file.path(out_dir, "asv_lengths.csv")
    if (file.exists(len_out)) {
      cat("  (asv_lengths.csv exists — skip)\n")
    } else {
      tax_raw  <- read.csv(tax_file, check.names=FALSE)
      asv_seqs <- as.character(tax_raw[[1]])
      asv_seqs <- gsub('^"|"$', '', asv_seqs)
      asv_seqs <- asv_seqs[!is.na(asv_seqs) & nchar(asv_seqs) > 10]
      lengths  <- nchar(asv_seqs)
      len_tbl  <- as.data.frame(table(Length=lengths), stringsAsFactors=FALSE)
      colnames(len_tbl) <- c("Length", "Count")
      len_tbl$Length    <- as.integer(len_tbl$Length)
      len_tbl$Count     <- as.integer(len_tbl$Count)
      len_tbl <- len_tbl[order(len_tbl$Length), ]
      write.csv(len_tbl, len_out, row.names=FALSE)
      cat("  ✓ asv_lengths.csv  (range:", min(lengths), "-", max(lengths),
          "bp,", length(asv_seqs), "ASVs)\n")
    }
  }, error=function(e) cat("  [skip] ASV lengths:", e$message, "\n"))

  cat("── DADA2 Extra Viz done ─────────────────────────────────────────\n\n")
  invisible(NULL)
}
