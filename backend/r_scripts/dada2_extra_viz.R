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

  # ── Helper: load metadata.csv → named vector SampleID → Group ────────────────
  meta_groups <- local({
    mf <- file.path(out_dir, "metadata.csv")
    if (!file.exists(mf)) return(NULL)
    tryCatch({
      m <- read.csv(mf, stringsAsFactors=FALSE)
      sid_col <- grep("sampleid|sample_id|sample", colnames(m), ignore.case=TRUE, value=TRUE)[1]
      grp_col <- grep("^group$", colnames(m), ignore.case=TRUE, value=TRUE)[1]
      if (is.na(sid_col) || is.na(grp_col)) return(NULL)
      setNames(m[[grp_col]], m[[sid_col]])
    }, error=function(e) NULL)
  })

  # ── Helper: compute 95% confidence ellipse polygon for a group of points ──────
  compute_ellipse <- function(x, y, conf=0.95, n_pts=100) {
    if (length(x) < 3) return(NULL)
    tryCatch({
      cov_m  <- cov(cbind(x, y))
      center <- c(mean(x), mean(y))
      chisq  <- qchisq(conf, df=2)
      eig    <- eigen(cov_m)
      vals   <- sqrt(pmax(eig$values, 0) * chisq)
      vecs   <- eig$vectors
      theta  <- seq(0, 2*pi, length.out=n_pts)
      pts    <- t(vecs %*% (vals * rbind(cos(theta), sin(theta)))) +
                matrix(rep(center, n_pts), nrow=n_pts, byrow=TRUE)
      data.frame(x=pts[,1], y=pts[,2])
    }, error=function(e) NULL)
  }

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
      # Attach Group from metadata if available
      if (!is.null(meta_groups)) {
        pc_df$Group <- meta_groups[pc_df$Sample]
        pc_df$Group[is.na(pc_df$Group)] <- "Unknown"
      }
      write.csv(pc_df, pca_out, row.names=FALSE)
      cat("  ✓ pca_scores.csv", if (!is.null(meta_groups)) "(+Group)" else "", "\n")

      n_sc  <- min(10, sum(eig_p > 0))
      sc_df <- data.frame(PC=seq_len(n_sc), Variance=round(var_p[seq_len(n_sc)], 2))
      write.csv(sc_df, file.path(out_dir, "pca_scree.csv"), row.names=FALSE)
      cat("  ✓ pca_scree.csv\n")

      # ── PCoA ellipse coordinates per group ──────────────────────────────────
      if (!is.null(meta_groups)) tryCatch({
        ellipse_rows <- list()
        for (grp in unique(pc_df$Group)) {
          idx <- which(pc_df$Group == grp)
          ell <- compute_ellipse(pc_df$PC1[idx], pc_df$PC2[idx])
          if (!is.null(ell)) {
            ell$Group <- grp
            ellipse_rows[[grp]] <- ell
          }
        }
        if (length(ellipse_rows) > 0) {
          ell_df <- do.call(rbind, ellipse_rows)
          rownames(ell_df) <- NULL
          write.csv(ell_df, file.path(out_dir, "pca_ellipse.csv"), row.names=FALSE)
          cat("  ✓ pca_ellipse.csv  (", length(ellipse_rows), "groups)\n")
        }
      }, error=function(e) cat("  [skip] PCoA ellipse:", e$message, "\n"))
    }
  }, error=function(e) cat("  [skip] PCoA:", e$message, "\n"))

  # ── 1b. NMDS Bray: add Group column from metadata ────────────────────────────
  # (only if nmds_bray.csv was already generated; patch it with Group column)
  if (!is.null(meta_groups)) {
    nmds_bray_out <- file.path(out_dir, "nmds_bray.csv")
    if (file.exists(nmds_bray_out)) tryCatch({
      nb <- read.csv(nmds_bray_out, stringsAsFactors=FALSE)
      if (!"Group" %in% colnames(nb)) {
        nb$Group <- meta_groups[nb$Sample]
        nb$Group[is.na(nb$Group)] <- "Unknown"
        write.csv(nb, nmds_bray_out, row.names=FALSE)
        cat("  ✓ nmds_bray.csv patched with Group\n")
        # Compute NMDS ellipse
        tryCatch({
          ellipse_rows <- list()
          for (grp in unique(nb$Group)) {
            idx <- which(nb$Group == grp)
            ell <- compute_ellipse(nb$NMDS1[idx], nb$NMDS2[idx])
            if (!is.null(ell)) { ell$Group <- grp; ellipse_rows[[grp]] <- ell }
          }
          if (length(ellipse_rows) > 0) {
            ell_df <- do.call(rbind, ellipse_rows); rownames(ell_df) <- NULL
            write.csv(ell_df, file.path(out_dir, "nmds_ellipse.csv"), row.names=FALSE)
            cat("  ✓ nmds_ellipse.csv\n")
          }
        }, error=function(e) NULL)
      }
    }, error=function(e) NULL)
  }

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

  # ── Helper: load asv_table.csv → samples × ASVs matrix (drops "sequence" col)
  load_asv_matrix <- function() {
    raw <- read.csv(asv_file, check.names=FALSE)
    seq_col <- which(colnames(raw) == "sequence")
    if (length(seq_col) > 0) raw <- raw[, -seq_col, drop=FALSE]
    mat <- t(as.matrix(raw))  # samples × ASVs
    storage.mode(mat) <- "numeric"
    mat
  }

  # ── 3. Jaccard heatmap + NMDS Jaccard ───────────────────────────────────────
  asv_t_cached <- NULL
  if (have_vegan && file.exists(asv_file)) tryCatch({
    asv_t_cached <- load_asv_matrix()

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
  # Uses vegan::rarefy if available, otherwise base-R exact formula (no dependency)
  rarefy_exact <- function(counts, depth) {
    # Hurlbert (1971) exact rarefaction: E[S] = sum_i(1 - C(N-ni, n) / C(N, n))
    counts <- as.integer(counts[counts > 0])
    total  <- sum(counts)
    if (depth >= total) return(sum(counts > 0))
    # log C(N, n) = lgamma(N+1) - lgamma(n+1) - lgamma(N-n+1)
    log_CN_n  <- lgamma(total + 1) - lgamma(depth + 1) - lgamma(total - depth + 1)
    expected  <- 0
    for (ni in counts) {
      if (ni > total - depth) {
        expected <- expected + 1          # ASV always sampled
      } else {
        log_CNni_n <- lgamma(total - ni + 1) - lgamma(depth + 1) - lgamma(total - ni - depth + 1)
        expected   <- expected + (1 - exp(log_CNni_n - log_CN_n))
      }
    }
    round(expected, 1)
  }

  if (file.exists(asv_file)) tryCatch({
    rar_out <- file.path(out_dir, "rarefaction.csv")
    if (file.exists(rar_out)) {
      cat("  (rarefaction.csv exists — skip)\n")
    } else {
      mat <- if (!is.null(asv_t_cached)) asv_t_cached else load_asv_matrix()

      row_sums  <- rowSums(mat)
      min_reads <- min(row_sums)
      if (min_reads < 10) stop("min sample reads too low for rarefaction")
      n_steps   <- min(25, min_reads)
      depths    <- unique(sort(c(
        round(seq(max(100, min_reads / n_steps), min_reads, length.out=n_steps)),
        as.integer(min_reads))))
      depths <- depths[depths >= 1 & depths <= min_reads]

      wide_df <- data.frame(Depth=depths)
      for (i in seq_len(nrow(mat))) {
        sname <- rownames(mat)[i]
        if (have_vegan) {
          rich <- suppressWarnings(rarefy(mat[i, , drop=FALSE], sample=depths))
          wide_df[[sname]] <- round(as.numeric(rich), 1)
        } else {
          wide_df[[sname]] <- sapply(depths, function(d) rarefy_exact(mat[i, ], d))
        }
      }
      write.csv(wide_df, rar_out, row.names=FALSE)
      method <- if (have_vegan) "vegan" else "base-R"
      cat("  ✓ rarefaction.csv  (", nrow(mat), "samples,",
          length(depths), "depths, min =", min_reads, "reads,", method, ")\n")
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
