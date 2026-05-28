#!/usr/bin/env Rscript
# =============================================================================
#  differential_analysis.R — ANOVA (all-groups) + Metastats / Pairwise
#
#  Usage: Rscript differential_analysis.R <output_dir>
#
#  Reads:  asv_table.csv + taxonomy_table.csv + metadata.csv (SampleID, Group)
#  Writes: anova_results.csv            — all-group Kruskal-Wallis results
#          r_plots/anova_bar.pdf        — grouped bar chart with error bars
#          metastats_{G1}_vs_{G2}.csv   — pairwise Wilcoxon per pair
#          r_plots/metastats_{G1}_vs_{G2}.pdf  — forest plot per pair
# =============================================================================

args    <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── Differential Analysis ───────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Required files ────────────────────────────────────────────────────────────
asv_f  <- file.path(out_dir, "asv_table.csv")
tax_f  <- NULL
for (lvl in c("species","genus","family","order","phylum")) {
  tf <- file.path(out_dir, paste0("taxonomy_", lvl, ".csv"))
  if (file.exists(tf)) { tax_f <- tf; tax_level <- lvl; break }
}
meta_f <- file.path(out_dir, "metadata.csv")

if (!file.exists(meta_f)) {
  cat("  [skip] metadata.csv not found — need SampleID + Group columns\n")
  quit(status=0)
}
if (is.null(tax_f) && !file.exists(asv_f)) {
  cat("  [skip] No ASV/taxonomy data found\n"); quit(status=0)
}

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.csv(meta_f, stringsAsFactors=FALSE)
sid_col <- grep("sampleid|sample_id|sample", colnames(meta), ignore.case=TRUE, value=TRUE)[1]
grp_col <- grep("^group$",                   colnames(meta), ignore.case=TRUE, value=TRUE)[1]
if (is.na(sid_col) || is.na(grp_col)) {
  cat("  [skip] metadata.csv must have SampleID and Group columns\n"); quit(status=0)
}
meta_map <- setNames(meta[[grp_col]], meta[[sid_col]])
groups   <- unique(meta[[grp_col]])
groups   <- groups[!is.na(groups) & groups != ""]
cat("  Groups:", paste(groups, collapse=", "), "\n")
if (length(groups) < 2) {
  cat("  [skip] Need at least 2 groups\n"); quit(status=0)
}

# ── Load abundance data ───────────────────────────────────────────────────────
abu_wide <- NULL   # samples × taxa data.frame

if (!is.null(tax_f)) {
  tryCatch({
    td <- read.csv(tax_f, check.names=FALSE, stringsAsFactors=FALSE)
    n_rows <- nrow(td) - 1; n_cols <- ncol(td) - 1
    if (n_cols >= n_rows) {
      # taxa-as-cols (regular 16S)
      abu_wide <- as.data.frame(td[, -1]); rownames(abu_wide) <- td[[1]]
    } else {
      # taxa-as-rows (ONT-16S)
      mat <- t(as.matrix(td[, -1])); storage.mode(mat) <- "numeric"
      rownames(mat) <- colnames(td)[-1]; colnames(mat) <- td[[1]]
      abu_wide <- as.data.frame(mat)
    }
    cat("  Taxonomy (", tax_level, "):", nrow(abu_wide), "samples,", ncol(abu_wide), "taxa\n")
  }, error=function(e) cat("  [warn] Taxonomy load failed:", e$message, "\n"))
}

if (is.null(abu_wide) && file.exists(asv_f)) {
  tryCatch({
    raw <- read.csv(asv_f, check.names=FALSE, stringsAsFactors=FALSE)
    seq_col <- which(colnames(raw) == "sequence")
    if (length(seq_col)) raw <- raw[, -seq_col, drop=FALSE]
    mat <- t(as.matrix(raw)); storage.mode(mat) <- "numeric"
    abu_wide <- as.data.frame(mat)
    tax_level <- "ASV"
    cat("  ASV table:", nrow(abu_wide), "samples,", ncol(abu_wide), "ASVs\n")
  }, error=function(e) cat("  [warn] ASV table load failed:", e$message, "\n"))
}

if (is.null(abu_wide)) {
  cat("  [skip] No abundance data\n"); quit(status=0)
}

# ── Convert to relative abundance ─────────────────────────────────────────────
row_sums  <- rowSums(abu_wide, na.rm=TRUE)
rel_wide  <- sweep(abu_wide, 1, pmax(row_sums, 1e-12), "/") * 100

# Filter to samples in metadata
common_s  <- intersect(rownames(rel_wide), names(meta_map))
if (length(common_s) < 2) {
  # try trimws/case-insensitive
  common_s <- rownames(rel_wide)[
    trimws(tolower(rownames(rel_wide))) %in% trimws(tolower(names(meta_map)))]
}
if (length(common_s) < 4) {
  cat("  [skip] Too few samples with group assignments (", length(common_s), ")\n")
  quit(status=0)
}
rel_wide  <- rel_wide[common_s, , drop=FALSE]
grp_vec   <- meta_map[common_s]
cat("  Samples with group data:", length(common_s), "\n")

# Top 30 taxa by mean relative abundance
tax_means <- colMeans(rel_wide, na.rm=TRUE)
top_taxa  <- names(sort(tax_means, decreasing=TRUE))[seq_len(min(30, ncol(rel_wide)))]
rel_sub   <- rel_wide[, top_taxa, drop=FALSE]

# ── Helper: group mean ± SE ───────────────────────────────────────────────────
group_stats <- function(values, groups_vec) {
  df_out <- data.frame()
  for (g in sort(unique(groups_vec))) {
    v  <- values[groups_vec == g]
    v  <- v[!is.na(v)]
    mn <- mean(v); se <- if (length(v) > 1) sd(v) / sqrt(length(v)) else 0
    df_out <- rbind(df_out, data.frame(Group=g, Mean=mn, SE=se, N=length(v),
                                        stringsAsFactors=FALSE))
  }
  df_out
}

# ── 1. Kruskal-Wallis all-groups ─────────────────────────────────────────────
anova_rows <- list()
for (taxon in top_taxa) {
  vals   <- rel_sub[, taxon]
  gstats <- group_stats(vals, grp_vec)
  # Kruskal test
  p_val  <- tryCatch({
    kt <- kruskal.test(vals ~ factor(grp_vec))
    kt$p.value
  }, error=function(e) NA)
  row <- data.frame(Taxon=taxon, stringsAsFactors=FALSE)
  for (i in seq_len(nrow(gstats))) {
    g <- gstats$Group[i]
    row[[paste0(g, "_mean")]] <- round(gstats$Mean[i], 4)
    row[[paste0(g, "_se"  )]] <- round(gstats$SE[i],   4)
  }
  row$p_kruskal <- round(p_val, 6)
  row$sig       <- if (!is.na(p_val)) {
    if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else
    if (p_val < 0.05)  "*"   else "ns"
  } else ""
  anova_rows[[taxon]] <- row
}

anova_df <- do.call(rbind, anova_rows)
anova_df <- anova_df[order(anova_df$p_kruskal, na.last=TRUE), ]
anova_csv <- file.path(out_dir, "anova_results.csv")
write.csv(anova_df, anova_csv, row.names=FALSE)
cat("  ✓ anova_results.csv  (", nrow(anova_df), "taxa)\n")

# ── 2. Pairwise Wilcoxon (Metastats) ─────────────────────────────────────────
pairs  <- combn(sort(groups), 2, simplify=FALSE)
p_maps <- list()  # track all generated pairwise CSV names

for (pr in pairs) {
  g1 <- pr[1]; g2 <- pr[2]
  idx1 <- which(grp_vec == g1); idx2 <- which(grp_vec == g2)
  if (length(idx1) == 0 || length(idx2) == 0) next

  pair_rows <- list()
  for (taxon in top_taxa) {
    v1 <- rel_sub[idx1, taxon]; v2 <- rel_sub[idx2, taxon]
    mn1 <- mean(v1, na.rm=TRUE); se1 <- if (length(v1)>1) sd(v1,na.rm=TRUE)/sqrt(sum(!is.na(v1))) else 0
    mn2 <- mean(v2, na.rm=TRUE); se2 <- if (length(v2)>1) sd(v2,na.rm=TRUE)/sqrt(sum(!is.na(v2))) else 0
    diff_mn <- mn1 - mn2
    # 95% CI on difference (bootstrap)
    p_val <- tryCatch({
      wilcox.test(v1, v2, exact=FALSE)$p.value
    }, error=function(e) NA)
    # Approximate CI using normal approx
    se_diff  <- sqrt(se1^2 + se2^2)
    ci_lo    <- diff_mn - 1.96 * se_diff
    ci_hi    <- diff_mn + 1.96 * se_diff
    pair_rows[[taxon]] <- data.frame(
      Taxon       = taxon,
      G1_mean     = round(mn1,    4),
      G1_se       = round(se1,    4),
      G2_mean     = round(mn2,    4),
      G2_se       = round(se2,    4),
      diff_mean   = round(diff_mn,4),
      CI_low      = round(ci_lo,  4),
      CI_high     = round(ci_hi,  4),
      p_value     = round(p_val,  6),
      q_value     = NA_real_,
      G1          = g1, G2 = g2,
      stringsAsFactors=FALSE)
  }
  pair_df <- do.call(rbind, pair_rows)
  # BH q-value
  pair_df$q_value <- round(p.adjust(pair_df$p_value, method="BH"), 6)
  pair_df <- pair_df[order(pair_df$p_value, na.last=TRUE), ]

  pair_csv <- file.path(out_dir, sprintf("metastats_%s_vs_%s.csv",
                                          gsub("[^a-zA-Z0-9]","_", g1),
                                          gsub("[^a-zA-Z0-9]","_", g2)))
  write.csv(pair_df, pair_csv, row.names=FALSE)
  p_maps[[length(p_maps)+1]] <- list(g1=g1, g2=g2, file=pair_csv, df=pair_df)
  cat("  ✓", basename(pair_csv), " (", nrow(pair_df), "taxa)\n")
}

# ── PDF: ANOVA grouped bar ────────────────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings=FALSE, recursive=TRUE)

anova_pdf <- file.path(plots_dir, "anova_bar.pdf")
tryCatch({
  has_gg <- requireNamespace("ggplot2", quietly=TRUE)
  if (has_gg) {
    suppressPackageStartupMessages(library(ggplot2))
    # Top 15 significant taxa
    sig_taxa <- head(anova_df$Taxon[!is.na(anova_df$p_kruskal) & anova_df$p_kruskal < 0.05], 15)
    if (length(sig_taxa) == 0) sig_taxa <- head(anova_df$Taxon, 15)

    long_df <- data.frame()
    for (taxon in sig_taxa) {
      gstats <- group_stats(rel_sub[, taxon], grp_vec)
      gstats$Taxon <- taxon
      long_df <- rbind(long_df, gstats)
    }
    long_df$Taxon <- factor(long_df$Taxon, levels=rev(sig_taxa))
    long_df$Group <- factor(long_df$Group, levels=sort(unique(long_df$Group)))

    grp_pal <- c("#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6",
                 "#06b6d4","#f97316","#84cc16","#ec4899","#14b8a6")
    ngrp <- length(unique(long_df$Group))

    # Add significance labels
    sig_label <- anova_df[match(sig_taxa, anova_df$Taxon), c("Taxon","sig","p_kruskal")]

    p <- ggplot(long_df, aes(x=Taxon, y=Mean, fill=Group)) +
      geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.75) +
      geom_errorbar(aes(ymin=pmax(Mean-SE,0), ymax=Mean+SE),
                    position=position_dodge(width=0.8), width=0.25, size=0.4) +
      scale_fill_manual(values=setNames(grp_pal[seq_len(ngrp)],
                                         levels(long_df$Group))) +
      coord_flip() +
      scale_x_discrete(labels=function(x) {
        lbl <- sig_label$sig[match(x, sig_label$Taxon)]
        paste0(x, " ", ifelse(is.na(lbl),"",lbl))
      }) +
      labs(x=NULL, y="Relative Abundance (%)",
           title="Kruskal-Wallis: Differential Abundance (p < 0.05 taxa)",
           fill="Group") +
      theme_classic(base_size=10) +
      theme(legend.position="right",
            plot.title=element_text(size=11, face="bold"),
            axis.text.y=element_text(size=9, face="italic"))

    h_pdf <- max(8, min(20, length(sig_taxa) * 0.5 + 3))
    ggsave(anova_pdf, p, width=14, height=h_pdf, limitsize=FALSE)
    cat("  ✓ anova_bar.pdf\n")
  } else {
    # base-R fallback
    sig_taxa <- head(anova_df$Taxon[anova_df$p_kruskal < 0.05], 15)
    if (length(sig_taxa) == 0) sig_taxa <- head(anova_df$Taxon, 10)
    pdf(anova_pdf, width=14, height=max(8, length(sig_taxa) * 0.5 + 2))
    cols_fill <- rainbow(length(groups))
    barplot(
      t(sapply(sig_taxa, function(tx) {
        sapply(sort(groups), function(g) mean(rel_sub[grp_vec==g, tx], na.rm=TRUE))
      })),
      beside=TRUE, horiz=FALSE, col=cols_fill,
      names.arg=sig_taxa, las=2, cex.names=0.6,
      main="Kruskal-Wallis: Differential Abundance", ylab="Relative Abundance (%)")
    legend("topright", legend=sort(groups), fill=cols_fill, cex=0.8)
    dev.off()
    cat("  ✓ anova_bar.pdf  (base-R)\n")
  }
}, error=function(e) cat("  [error] anova_bar.pdf:", e$message, "\n"))

# ── PDF: Metastats forest plots ───────────────────────────────────────────────
for (pm in p_maps) {
  pair_pdf <- file.path(plots_dir,
    sprintf("metastats_%s_vs_%s.pdf",
            gsub("[^a-zA-Z0-9]","_", pm$g1),
            gsub("[^a-zA-Z0-9]","_", pm$g2)))
  tryCatch({
    df  <- pm$df
    sig <- df[!is.na(df$p_value) & df$p_value < 0.05, ]
    if (nrow(sig) == 0) sig <- head(df, 20)
    sig <- head(sig, 25)
    sig$Taxon <- factor(sig$Taxon, levels=rev(sig$Taxon))

    has_gg <- requireNamespace("ggplot2", quietly=TRUE)
    if (has_gg) {
      suppressPackageStartupMessages(library(ggplot2))
      # Left: bar chart of means per group
      long_s <- rbind(
        data.frame(Taxon=sig$Taxon, Group=pm$g1, Mean=sig$G1_mean, SE=sig$G1_se),
        data.frame(Taxon=sig$Taxon, Group=pm$g2, Mean=sig$G2_mean, SE=sig$G2_se))
      long_s$Taxon <- factor(long_s$Taxon, levels=levels(sig$Taxon))

      p_bar <- ggplot(long_s, aes(x=Taxon, y=Mean, fill=Group)) +
        geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.75) +
        geom_errorbar(aes(ymin=pmax(Mean-SE,0), ymax=Mean+SE),
                      position=position_dodge(width=0.8), width=0.25, size=0.4) +
        scale_fill_manual(values=c("#3b82f6","#ef4444")) +
        coord_flip() +
        labs(x=NULL, y="Relative Abundance (%)", fill="Group") +
        theme_classic(base_size=9) +
        theme(axis.text.y=element_text(size=8, face="italic"))

      # Right: CI plot
      p_ci <- ggplot(sig, aes(x=diff_mean, y=Taxon)) +
        geom_vline(xintercept=0, linetype="dashed", color="grey40") +
        geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0.25, size=0.5, color="#475569") +
        geom_point(size=2.5, color="#ef4444") +
        annotate("text", x=max(sig$CI_high, na.rm=TRUE)*1.05,
                 y=nrow(sig)/2, label="p-value", size=2.5, hjust=0) +
        geom_text(aes(x=max(sig$CI_high, na.rm=TRUE)*1.1,
                      label=sprintf("%.4f", p_value)), size=2.5, hjust=0) +
        labs(x=paste0("Difference (%)\n(", pm$g1, " − ", pm$g2, ")"), y=NULL,
             title=paste0(pm$g1, " vs ", pm$g2, " (Wilcoxon, 95% CI)")) +
        theme_classic(base_size=9) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              plot.title=element_text(size=10, face="bold"))

      if (requireNamespace("patchwork", quietly=TRUE)) {
        suppressPackageStartupMessages(library(patchwork))
        combined <- p_bar + p_ci + plot_layout(widths=c(1.5,1))
      } else { combined <- p_bar }

      h_pdf <- max(6, min(18, nrow(sig) * 0.4 + 3))
      ggsave(pair_pdf, combined, width=14, height=h_pdf, limitsize=FALSE)
      cat("  ✓", basename(pair_pdf), "\n")
    } else {
      pdf(pair_pdf, width=12, height=max(6, nrow(sig)*0.4+2))
      long_s <- rbind(
        data.frame(Taxon=as.character(sig$Taxon), Group=pm$g1, Mean=sig$G1_mean, SE=sig$G1_se),
        data.frame(Taxon=as.character(sig$Taxon), Group=pm$g2, Mean=sig$G2_mean, SE=sig$G2_se))
      with(long_s, barplot(tapply(Mean, list(Taxon, Group), mean),
                            beside=TRUE, col=c("#3b82f6","#ef4444"), las=2, cex.names=0.5,
                            main=paste0(pm$g1, " vs ", pm$g2)))
      legend("topright", legend=c(pm$g1, pm$g2), fill=c("#3b82f6","#ef4444"), cex=0.8)
      dev.off()
      cat("  ✓", basename(pair_pdf), "  (base-R)\n")
    }
  }, error=function(e) cat("  [error]", basename(pair_pdf), ":", e$message, "\n"))
}

cat("── Differential Analysis done ──────────────────────────────────────\n\n")
