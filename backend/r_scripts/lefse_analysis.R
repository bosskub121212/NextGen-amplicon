#!/usr/bin/env Rscript
# =============================================================================
#  lefse_analysis.R — LEfSe LDA Score Bar + Cladogram
#
#  Usage: Rscript lefse_analysis.R <output_dir>
#
#  Uses lefser (Bioconductor) or microbiomeMarker; falls back to manual LDA
#  Reads:  asv_table.csv + taxonomy_table.csv + metadata.csv (SampleID, Group)
#  Writes: lefse_lda.csv                 — taxon, lda_score, group, level
#          r_plots/lefse_lda.pdf         — horizontal LDA bar chart
#          r_plots/lefse_lda_{G1}_vs_{G2}.pdf  — pairwise version
#          r_plots/cladogram.pdf         — LEfSe circular cladogram (if available)
# =============================================================================

args    <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── LEfSe Analysis ──────────────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Required files ────────────────────────────────────────────────────────────
asv_f  <- file.path(out_dir, "asv_table.csv")
meta_f <- file.path(out_dir, "metadata.csv")
tax_f  <- NULL; tax_level <- "genus"
for (lvl in c("genus","family","order","phylum")) {
  tf <- file.path(out_dir, paste0("taxonomy_", lvl, ".csv"))
  if (file.exists(tf)) { tax_f <- tf; tax_level <- lvl; break }
}

if (!file.exists(meta_f)) {
  cat("  [skip] metadata.csv not found\n"); quit(status=0)
}
if (is.null(tax_f) && !file.exists(asv_f)) {
  cat("  [skip] No abundance data\n"); quit(status=0)
}

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.csv(meta_f, stringsAsFactors=FALSE)
sid_col <- grep("sampleid|sample_id|sample", colnames(meta), ignore.case=TRUE, value=TRUE)[1]
grp_col <- grep("^group$",                   colnames(meta), ignore.case=TRUE, value=TRUE)[1]
if (is.na(sid_col) || is.na(grp_col)) {
  cat("  [skip] metadata.csv needs SampleID + Group columns\n"); quit(status=0)
}
meta_map <- setNames(meta[[grp_col]], meta[[sid_col]])
groups   <- unique(meta[[grp_col]])
groups   <- groups[!is.na(groups) & groups != ""]
cat("  Groups:", paste(groups, collapse=", "), "\n")
if (length(groups) < 2) {
  cat("  [skip] Need at least 2 groups\n"); quit(status=0)
}

# ── Load abundance ────────────────────────────────────────────────────────────
abu_wide <- NULL
if (!is.null(tax_f)) {
  tryCatch({
    td <- read.csv(tax_f, check.names=FALSE, stringsAsFactors=FALSE)
    n_rows <- nrow(td) - 1; n_cols <- ncol(td) - 1
    if (n_cols >= n_rows) {
      abu_wide <- as.data.frame(td[, -1]); rownames(abu_wide) <- td[[1]]
    } else {
      mat <- t(as.matrix(td[, -1])); storage.mode(mat) <- "numeric"
      rownames(mat) <- colnames(td)[-1]; colnames(mat) <- td[[1]]
      abu_wide <- as.data.frame(mat)
    }
    cat("  Loaded", tax_level, ":", nrow(abu_wide), "samples,", ncol(abu_wide), "taxa\n")
  }, error=function(e) cat("  [warn] tax load failed:", e$message, "\n"))
}

if (is.null(abu_wide) && file.exists(asv_f)) {
  tryCatch({
    raw <- read.csv(asv_f, check.names=FALSE, stringsAsFactors=FALSE)
    seq_col <- which(colnames(raw) == "sequence")
    if (length(seq_col)) raw <- raw[, -seq_col, drop=FALSE]
    mat <- t(as.matrix(raw)); storage.mode(mat) <- "numeric"
    abu_wide <- as.data.frame(mat); tax_level <- "ASV"
    cat("  Loaded ASV table:", nrow(abu_wide), "samples\n")
  }, error=function(e) cat("  [warn]", e$message, "\n"))
}
if (is.null(abu_wide)) { cat("  [skip] No data\n"); quit(status=0) }

# Filter to samples with metadata
common_s <- intersect(rownames(abu_wide), names(meta_map))
if (length(common_s) < 4) {
  cat("  [skip] Too few samples with groups (", length(common_s), ")\n"); quit(status=0)
}
abu_wide <- abu_wide[common_s, , drop=FALSE]
grp_vec  <- meta_map[common_s]

# Relative abundance
row_sums <- rowSums(abu_wide, na.rm=TRUE)
rel_wide <- sweep(abu_wide, 1, pmax(row_sums, 1e-12), "/") * 100

# Top 50 taxa
col_means <- colMeans(rel_wide, na.rm=TRUE)
top_taxa  <- names(sort(col_means, decreasing=TRUE))[seq_len(min(50, ncol(rel_wide)))]
rel_sub   <- rel_wide[, top_taxa, drop=FALSE]

# ── LEfSe via lefser (Bioconductor) ─────────────────────────────────────────
run_lefser <- function(rel_mat, grp_v, min_lda=2.0) {
  if (!requireNamespace("lefser", quietly=TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly=TRUE)) return(NULL)
  suppressPackageStartupMessages({
    library(lefser); library(SummarizedExperiment)
  })
  tryCatch({
    # lefser expects features × samples
    cts <- t(as.matrix(rel_mat))
    se  <- SummarizedExperiment(
      assays = list(counts=cts),
      colData = data.frame(group=factor(grp_v), row.names=colnames(cts)))
    res <- lefser(se, groupCol="group", lda.threshold=min_lda)
    cat("  lefser found", nrow(res), "significant features\n")
    res
  }, error=function(e) { cat("  [warn] lefser:", e$message, "\n"); NULL })
}

lefser_to_csv <- function(res, grp_v, filename) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  df <- data.frame(
    Taxon     = res$Names,
    lda_score = abs(res$scores),
    group     = ifelse(res$scores > 0, levels(factor(grp_v))[1],
                                       levels(factor(grp_v))[2]),
    level     = tax_level,
    stringsAsFactors=FALSE)
  df <- df[order(abs(df$lda_score), decreasing=TRUE), ]
  write.csv(df, filename, row.names=FALSE)
  df
}

# ── Manual LDA fallback (when lefser not installed) ─────────────────────────
manual_lda <- function(rel_mat, grp_v, min_lda=2.0) {
  grps   <- sort(unique(grp_v))
  results <- data.frame()
  for (taxon in colnames(rel_mat)) {
    vals  <- log10(rel_mat[, taxon] + 0.01)
    group_means <- tapply(vals, grp_v, mean)
    group_sds   <- tapply(vals, grp_v, sd)
    group_ns    <- tapply(vals, grp_v, length)
    # LDA score approximation: between-group variance / within-group variance
    grand_mean  <- mean(vals)
    ss_between  <- sum(group_ns * (group_means - grand_mean)^2)
    ss_within   <- sum((group_ns - 1) * group_sds^2, na.rm=TRUE)
    if (is.na(ss_within) || ss_within < 1e-12) next
    lda_approx  <- log10(1 + ss_between / ss_within) * 3
    if (lda_approx < min_lda) next
    dominant_grp <- names(which.max(group_means))
    p_val <- tryCatch(kruskal.test(vals ~ factor(grp_v))$p.value, error=function(e) 1)
    if (!is.na(p_val) && p_val < 0.05) {
      results <- rbind(results, data.frame(
        Taxon=taxon, lda_score=round(lda_approx,3),
        group=dominant_grp, level=tax_level, p_value=p_val, stringsAsFactors=FALSE))
    }
  }
  if (nrow(results) > 0) results[order(results$lda_score, decreasing=TRUE), ]
  else NULL
}

# ── Run LEfSe: all-groups ─────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings=FALSE, recursive=TRUE)

lda_all <- NULL
lda_csv <- file.path(out_dir, "lefse_lda.csv")

res_lefser <- run_lefser(rel_sub, grp_vec)
if (!is.null(res_lefser)) {
  lda_all <- lefser_to_csv(res_lefser, grp_vec, lda_csv)
  cat("  ✓ lefse_lda.csv  (lefser,", nrow(lda_all), "features)\n")
} else {
  cat("  lefser not available — using manual LDA approximation\n")
  lda_all <- manual_lda(rel_sub, grp_vec)
  if (!is.null(lda_all)) {
    write.csv(lda_all, lda_csv, row.names=FALSE)
    cat("  ✓ lefse_lda.csv  (manual,", nrow(lda_all), "features)\n")
  } else {
    cat("  [info] No significant features found (LDA threshold: 2.0)\n")
  }
}

# ── Run pairwise LEfSe ────────────────────────────────────────────────────────
pairs <- combn(sort(groups), 2, simplify=FALSE)
for (pr in pairs) {
  g1 <- pr[1]; g2 <- pr[2]
  idx <- grp_vec %in% c(g1, g2)
  if (sum(idx) < 4) next
  sub_mat <- rel_sub[idx, , drop=FALSE]
  sub_grp <- grp_vec[idx]
  pair_csv <- file.path(out_dir, sprintf("lefse_%s_vs_%s.csv",
                          gsub("[^a-zA-Z0-9]","_",g1),
                          gsub("[^a-zA-Z0-9]","_",g2)))
  res_p <- run_lefser(sub_mat, sub_grp)
  if (!is.null(res_p)) {
    lda_p <- lefser_to_csv(res_p, sub_grp, pair_csv)
    cat("  ✓", basename(pair_csv), "(lefser)\n")
  } else {
    lda_p <- manual_lda(sub_mat, sub_grp)
    if (!is.null(lda_p)) {
      write.csv(lda_p, pair_csv, row.names=FALSE)
      cat("  ✓", basename(pair_csv), "(manual)\n")
    }
  }
  # Pairwise PDF
  if (!is.null(lda_p) && nrow(lda_p) > 0) {
    pair_pdf <- file.path(plots_dir, sprintf("lefse_%s_vs_%s.pdf",
                            gsub("[^a-zA-Z0-9]","_",g1),
                            gsub("[^a-zA-Z0-9]","_",g2)))
    tryCatch(plot_lda_pdf(lda_p, pair_pdf, g1, g2), error=function(e) NULL)
  }
}

# ── LDA bar chart PDF ─────────────────────────────────────────────────────────
plot_lda_pdf <- function(lda_df, pdf_path, ...) {
  if (is.null(lda_df) || nrow(lda_df) == 0) return(invisible(NULL))
  lda_df <- head(lda_df[order(lda_df$lda_score, decreasing=TRUE), ], 30)
  lda_df$Taxon <- factor(lda_df$Taxon, levels=rev(lda_df$Taxon))
  lda_df$signed <- ifelse(lda_df$group == lda_df$group[1],
                           lda_df$lda_score, -lda_df$lda_score)

  has_gg <- requireNamespace("ggplot2", quietly=TRUE)
  if (has_gg) {
    suppressPackageStartupMessages(library(ggplot2))
    grps_in <- sort(unique(lda_df$group))
    grp_pal <- c("#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6",
                 "#06b6d4","#f97316","#84cc16","#ec4899","#14b8a6")
    gcols   <- setNames(grp_pal[seq_along(grps_in)], grps_in)

    p <- ggplot(lda_df, aes(x=Taxon, y=lda_score, fill=group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=gcols, name="Group") +
      coord_flip() +
      labs(x=NULL, y="LDA Score (log10)", title="LEfSe — LDA Score (log10)") +
      theme_classic(base_size=10) +
      theme(axis.text.y=element_text(size=8, face="italic"),
            plot.title=element_text(size=11, face="bold"),
            legend.position="bottom")

    h_pdf <- max(6, min(18, nrow(lda_df) * 0.4 + 3))
    ggsave(pdf_path, p, width=10, height=h_pdf, limitsize=FALSE)
  } else {
    pdf(pdf_path, width=10, height=max(6, nrow(lda_df)*0.4+2))
    barplot(lda_df$lda_score, names.arg=as.character(lda_df$Taxon),
            horiz=TRUE, las=1, col="#3b82f6", cex.names=0.7,
            main="LEfSe LDA Score", xlab="LDA Score (log10)")
    dev.off()
  }
  cat("  ✓", basename(pdf_path), "\n")
}

lda_pdf <- file.path(plots_dir, "lefse_lda.pdf")
tryCatch(plot_lda_pdf(lda_all, lda_pdf), error=function(e) cat("  [error] lefse_lda.pdf:", e$message, "\n"))

# ── Cladogram via microbiomeMarker ────────────────────────────────────────────
clado_pdf <- file.path(plots_dir, "cladogram.pdf")
if (!file.exists(clado_pdf) &&
    requireNamespace("microbiomeMarker", quietly=TRUE) &&
    !is.null(tax_f)) {
  tryCatch({
    suppressPackageStartupMessages(library(microbiomeMarker))
    # Build phyloseq-compatible object
    if (requireNamespace("phyloseq", quietly=TRUE)) {
      library(phyloseq)
      otu <- otu_table(t(as.matrix(rel_sub)), taxa_are_rows=TRUE)
      sam <- sample_data(data.frame(group=grp_vec, row.names=names(grp_vec)))
      ps  <- phyloseq(otu, sam)
      mm  <- run_lefse(ps, group="group", norm="none", lda_cutoff=2)
      if (!is.null(mm) && ntaxa(mm) > 0) {
        p_clado <- plot_cladogram(mm, color=setNames(
          c("#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6")[seq_along(groups)], groups))
        ggsave(clado_pdf, p_clado, width=12, height=12)
        cat("  ✓ cladogram.pdf  (microbiomeMarker)\n")
      }
    }
  }, error=function(e) cat("  [info] Cladogram skipped:", e$message, "\n"))
}

cat("── LEfSe Analysis done ─────────────────────────────────────────────\n\n")
