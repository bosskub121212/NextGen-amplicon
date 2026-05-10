#!/usr/bin/env Rscript
# =============================================================================
#  NextGen-Amplicon — R Visualization Pipeline
#  Called by qiime2_pipeline.py after QIIME2 export
#
#  Usage:
#    Rscript viz_pipeline.R \
#      --output_dir /path/to/job_output \
#      --metadata   /path/to/metadata.tsv \   # optional
#      --group_col  treatment \               # column in metadata for grouping
#      --marker     16S                       # 16S|12S|ITS1|ITS2|COX1|18S-nema|PacBio
#
#  Expects QIIME2 exported files under output_dir/exported/:
#    feature-table/feature-table.biom
#    feature-table/feature-table.tsv
#    taxonomy/taxonomy.tsv
#    tree/tree.nwk          (16S / 12S / PacBio only)
#
#  Output: output_dir/r_plots/*.pdf  +  output_dir/r_tables/*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--output_dir", type="character", default=NULL,
              help="QIIME2 pipeline output directory"),
  make_option("--metadata",   type="character", default=NULL,
              help="QIIME2 metadata TSV file (optional)"),
  make_option("--group_col",  type="character", default="treatment",
              help="Metadata column to use for grouping"),
  make_option("--marker",     type="character", default="16S",
              help="Marker type: 16S|12S|ITS1|ITS2|COX1|18S-nema|PacBio")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$output_dir)) {
  cat("ERROR: --output_dir is required\n")
  quit(status=1)
}

OUTPUT_DIR  <- normalizePath(opt$output_dir, mustWork=FALSE)
MARKER      <- opt$marker
GROUP_COL   <- opt$group_col
METADATA_FILE <- opt$metadata

PLOTS_DIR  <- file.path(OUTPUT_DIR, "r_plots")
TABLES_DIR <- file.path(OUTPUT_DIR, "r_tables")
EXP_DIR    <- file.path(OUTPUT_DIR, "exported")

dir.create(PLOTS_DIR,  recursive=TRUE, showWarnings=FALSE)
dir.create(TABLES_DIR, recursive=TRUE, showWarnings=FALSE)

cat("\n============================================================\n")
cat("  NextGen-Amplicon R Visualization Pipeline\n")
cat(sprintf("  Marker  : %s\n", MARKER))
cat(sprintf("  Output  : %s\n", OUTPUT_DIR))
cat(sprintf("  Plots   : %s\n", PLOTS_DIR))
cat("============================================================\n\n")

# ── User library ──────────────────────────────────────────────────────────────
user_lib <- Sys.getenv("R_LIBS_USER",
             unset=file.path(Sys.getenv("HOME"), "R", "library"))
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

# ── Load packages ─────────────────────────────────────────────────────────────
load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat(sprintf("[WARN] Package '%s' not available — skipping\n", pkg))
    return(FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only=TRUE))
  return(TRUE)
}

has_phyloseq <- load_pkg("phyloseq")
has_vegan    <- load_pkg("vegan")
has_ggplot2  <- load_pkg("ggplot2")
has_dplyr    <- load_pkg("dplyr")
has_tidyr    <- load_pkg("tidyr")
has_ancombc  <- load_pkg("ANCOMBC")
has_funguild <- load_pkg("FUNGuildR")

# Extra plotting helpers (non-fatal if missing)
if (has_ggplot2) {
  load_pkg("RColorBrewer")
  load_pkg("ggrepel")
  load_pkg("patchwork")
  load_pkg("cowplot")
  load_pkg("pheatmap")
}

# ── Color palette helper ──────────────────────────────────────────────────────
make_palette <- function(n) {
  if (n <= 8 && requireNamespace("RColorBrewer", quietly=TRUE)) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else {
    scales::hue_pal()(n)
  }
}

# ── Load exported QIIME2 data ─────────────────────────────────────────────────
cat("── Loading exported QIIME2 data ──────────────────────────────\n")

biom_file  <- file.path(EXP_DIR, "feature-table", "feature-table.biom")
tax_file   <- file.path(EXP_DIR, "taxonomy",       "taxonomy.tsv")
tree_file  <- file.path(EXP_DIR, "tree",           "tree.nwk")
seq_fasta  <- file.path(EXP_DIR, "rep-seqs",       "dna-sequences.fasta")

if (!file.exists(biom_file)) {
  cat(sprintf("[ERROR] BIOM file not found: %s\n", biom_file))
  quit(status=1)
}

# Use phyloseq to import BIOM
ps <- NULL
if (has_phyloseq) {
  tryCatch({
    ps <- import_biom(biom_file)
    cat(sprintf("  Loaded BIOM: %d ASVs × %d samples\n",
                ntaxa(ps), nsamples(ps)))
  }, error=function(e) cat(sprintf("[WARN] BIOM import error: %s\n", e$message)))
}

# Load taxonomy TSV
tax_df <- NULL
if (file.exists(tax_file)) {
  tryCatch({
    tax_df <- read.table(tax_file, header=TRUE, sep="\t",
                         comment.char="", quote="", stringsAsFactors=FALSE)
    cat(sprintf("  Loaded taxonomy: %d features\n", nrow(tax_df)))

    if (has_phyloseq && !is.null(ps)) {
      # Parse taxonomy into ranks
      tax_col <- if ("Taxon" %in% colnames(tax_df)) "Taxon" else colnames(tax_df)[2]
      tax_split <- strsplit(tax_df[[tax_col]], ";\\s*")
      max_ranks <- max(sapply(tax_split, length))
      ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[seq_len(min(7, max_ranks))]
      tax_mat <- do.call(rbind, lapply(tax_split, function(x) {
        x <- sub("^[a-z]__", "", x)        # strip k__, p__, etc.
        x <- sub("^\\s+|\\s+$", "", x)
        length(x) <- length(ranks)
        x
      }))
      colnames(tax_mat) <- ranks
      rownames(tax_mat) <- tax_df[[1]]
      # Align with phyloseq ASV IDs
      common <- intersect(taxa_names(ps), rownames(tax_mat))
      if (length(common) > 0) {
        tax_table(ps) <- tax_table(tax_mat[common, , drop=FALSE])
      }
    }
  }, error=function(e) cat(sprintf("[WARN] Taxonomy load error: %s\n", e$message)))
}

# Load tree
if (file.exists(tree_file) && has_phyloseq && !is.null(ps)) {
  tryCatch({
    phy_tree(ps) <- read_tree(tree_file)
    cat("  Loaded phylogenetic tree\n")
  }, error=function(e) cat(sprintf("[WARN] Tree load error: %s\n", e$message)))
}

# Load metadata
meta_df <- NULL
if (!is.null(METADATA_FILE) && file.exists(METADATA_FILE)) {
  tryCatch({
    meta_df <- read.table(METADATA_FILE, header=TRUE, sep="\t",
                          comment.char="#", quote="", stringsAsFactors=FALSE)
    colnames(meta_df)[1] <- "SampleID"
    rownames(meta_df)    <- meta_df$SampleID
    cat(sprintf("  Loaded metadata: %d samples × %d columns\n",
                nrow(meta_df), ncol(meta_df)))
    if (has_phyloseq && !is.null(ps)) {
      common_samps <- intersect(sample_names(ps), rownames(meta_df))
      if (length(common_samps) > 0) {
        sample_data(ps) <- sample_data(meta_df[common_samps, , drop=FALSE])
        cat(sprintf("  Matched %d samples between BIOM and metadata\n",
                    length(common_samps)))
      }
    }
  }, error=function(e) cat(sprintf("[WARN] Metadata load error: %s\n", e$message)))
}

if (is.null(ps) || ntaxa(ps) == 0) {
  cat("[ERROR] Could not build phyloseq object — no ASVs loaded\n")
  quit(status=1)
}

has_meta  <- !is.null(meta_df) && GROUP_COL %in% colnames(meta_df)
has_tree  <- !is.null(phy_tree(ps)) && !inherits(tryCatch(phy_tree(ps), error=function(e) e), "error")
is_ITS    <- MARKER %in% c("ITS1", "ITS2")
is_COX1   <- MARKER == "COX1"

cat(sprintf("\n  has_metadata: %s | has_tree: %s | group_col: '%s'\n",
            has_meta, has_tree, GROUP_COL))

# ── Helper: safe PDF save ──────────────────────────────────────────────────────
save_pdf <- function(plot_obj, filename, width=10, height=7) {
  path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    if (inherits(plot_obj, "ggplot")) {
      ggsave(path, plot=plot_obj, width=width, height=height, device="pdf")
    } else {
      pdf(path, width=width, height=height)
      if (is.function(plot_obj)) plot_obj() else print(plot_obj)
      dev.off()
    }
    cat(sprintf("  ✓ Saved: %s\n", basename(path)))
  }, error=function(e) cat(sprintf("  ✗ %s: %s\n", basename(path), e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Read tracking & ASV summary
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 1: ASV Summary ────────────────────────────────────\n")
tryCatch({
  asv_counts <- sort(sample_sums(ps), decreasing=TRUE)
  asv_tbl <- data.frame(
    Sample    = names(asv_counts),
    Reads     = as.integer(asv_counts),
    n_ASV     = apply(otu_table(ps), 2, function(x) sum(x > 0))
  )
  write.csv(asv_tbl, file.path(TABLES_DIR, "asv_summary.csv"), row.names=FALSE)

  if (has_ggplot2) {
    p <- ggplot(asv_tbl, aes(x=reorder(Sample, -Reads), y=Reads)) +
      geom_col(fill="#3b82f6", alpha=0.85) +
      geom_text(aes(label=Reads), vjust=-0.4, size=2.8) +
      labs(title="Read counts per sample", x="Sample", y="Total reads") +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1, size=9))
    save_pdf(p, "01_read_counts.pdf", width=max(8, nsamples(ps)*0.5), height=5)
  }
}, error=function(e) cat(sprintf("[WARN] ASV summary: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Alpha diversity
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 2: Alpha Diversity ────────────────────────────────\n")
tryCatch({
  alpha_measures <- c("Observed","Shannon","Simpson","Chao1","ACE")
  alpha_df <- estimate_richness(ps, measures=alpha_measures)
  alpha_df$Sample <- rownames(alpha_df)

  # Add metadata grouping if available
  if (has_meta) {
    alpha_df[[GROUP_COL]] <- meta_df[alpha_df$Sample, GROUP_COL]
  }

  write.csv(alpha_df, file.path(TABLES_DIR, "alpha_diversity.csv"), row.names=FALSE)
  cat(sprintf("  Computed alpha diversity for %d samples\n", nrow(alpha_df)))

  if (has_ggplot2) {
    plot_measures <- intersect(c("Shannon","Observed","Chao1","Simpson"), colnames(alpha_df))
    for (m in plot_measures) {
      gg <- ggplot(alpha_df, aes_string(
          x=if (has_meta) GROUP_COL else "Sample",
          y=m)) +
        geom_boxplot(aes_string(fill=if (has_meta) GROUP_COL else "Sample"),
                     outlier.shape=NA, alpha=0.7) +
        geom_jitter(width=0.15, size=2, alpha=0.8) +
        labs(title=sprintf("Alpha Diversity — %s", m),
             x=if (has_meta) GROUP_COL else "Sample", y=m) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=30, hjust=1),
              legend.position="none")
      save_pdf(gg, sprintf("02_alpha_%s.pdf", tolower(m)))
    }

    # Faith's PD if tree
    if (has_tree && has_vegan) {
      tryCatch({
        load_pkg("picante")
        otu_mat <- t(as.matrix(otu_table(ps)))
        tree_obj <- phy_tree(ps)
        pd_res <- pd(otu_mat, tree_obj, include.root=TRUE)
        pd_df  <- data.frame(Sample=rownames(pd_res), PD=pd_res$PD)
        if (has_meta) pd_df[[GROUP_COL]] <- meta_df[pd_df$Sample, GROUP_COL]
        write.csv(pd_df, file.path(TABLES_DIR, "faiths_pd.csv"), row.names=FALSE)

        gg <- ggplot(pd_df, aes_string(
            x=if (has_meta) GROUP_COL else "Sample", y="PD")) +
          geom_boxplot(aes_string(fill=if (has_meta) GROUP_COL else "Sample"),
                       outlier.shape=NA, alpha=0.7) +
          geom_jitter(width=0.15, size=2, alpha=0.8) +
          labs(title="Faith's Phylogenetic Diversity", x=GROUP_COL, y="PD") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=30, hjust=1),
                legend.position="none")
        save_pdf(gg, "02_alpha_faithPD.pdf")
      }, error=function(e) cat(sprintf("  [WARN] Faith's PD: %s\n", e$message)))
    }
  }
}, error=function(e) cat(sprintf("[WARN] Alpha diversity: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Rarefaction curves
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 3: Rarefaction Curves ─────────────────────────────\n")
tryCatch({
  otu_mat <- as.matrix(otu_table(ps))
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  max_depth <- min(max(rowSums(otu_mat)), 50000)
  steps     <- unique(c(seq(100, max_depth, by=max(100, max_depth %/% 50)), max_depth))

  rarefy_row <- function(x, step) {
    sapply(steps[steps <= sum(x)], function(s) {
      rarefy(matrix(x, nrow=1), s)[1]
    })
  }

  rar_list <- lapply(seq_len(nrow(otu_mat)), function(i) {
    sname <- rownames(otu_mat)[i]
    x     <- otu_mat[i, ]
    valid <- steps[steps <= sum(x)]
    if (length(valid) == 0) return(NULL)
    data.frame(
      Sample = sname,
      Depth  = valid,
      ASVs   = sapply(valid, function(s) vegan::rarefy(matrix(x, nrow=1), s)[1])
    )
  })
  rar_df <- do.call(rbind, Filter(Negate(is.null), rar_list))

  if (!is.null(rar_df) && nrow(rar_df) > 0) {
    if (has_meta && has_ggplot2) {
      rar_df[[GROUP_COL]] <- meta_df[rar_df$Sample, GROUP_COL]
    }
    write.csv(rar_df, file.path(TABLES_DIR, "rarefaction.csv"), row.names=FALSE)

    if (has_ggplot2) {
      n_samp <- length(unique(rar_df$Sample))
      pal    <- make_palette(n_samp)
      p <- ggplot(rar_df, aes(x=Depth, y=ASVs, color=Sample)) +
        geom_line(size=0.8, alpha=0.85) +
        scale_color_manual(values=pal) +
        labs(title="Rarefaction Curves", x="Sequencing depth", y="Observed ASVs") +
        theme_bw() +
        theme(legend.text=element_text(size=8),
              legend.key.size=unit(0.5,"cm"))
      save_pdf(p, "03_rarefaction_curves.pdf", width=10, height=6)
    }
  }
}, error=function(e) cat(sprintf("[WARN] Rarefaction: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Taxonomy bar plots
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 4: Taxonomy Bar Plots ─────────────────────────────\n")
tryCatch({
  if (!has_ggplot2) stop("ggplot2 not available")
  tax_ranks_avail <- rank_names(ps)
  plot_ranks <- intersect(c("Phylum","Class","Order","Family","Genus"), tax_ranks_avail)

  for (rank in plot_ranks) {
    tryCatch({
      ps_glom <- tax_glom(ps, taxrank=rank, NArm=FALSE)
      ps_rel   <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
      melt_df  <- psmelt(ps_rel)
      melt_df[[rank]][is.na(melt_df[[rank]])] <- "Unclassified"

      # Top 15 taxa; merge rest as "Other"
      top_taxa <- melt_df %>%
        dplyr::group_by(.data[[rank]]) %>%
        dplyr::summarise(mean_abund = mean(Abundance, na.rm=TRUE)) %>%
        dplyr::arrange(dplyr::desc(mean_abund)) %>%
        dplyr::slice_head(n=15) %>%
        dplyr::pull(.data[[rank]])

      melt_df$TaxLabel <- ifelse(melt_df[[rank]] %in% top_taxa,
                                 melt_df[[rank]], "Other")
      n_col <- length(unique(melt_df$TaxLabel))
      pal   <- c(make_palette(min(15, n_col-1)), "grey80")[seq_len(n_col)]

      p <- ggplot(melt_df, aes(x=Sample, y=Abundance, fill=TaxLabel)) +
        geom_bar(stat="identity", width=0.85) +
        scale_fill_manual(values=pal, name=rank) +
        labs(title=sprintf("Relative Abundance — %s level", rank),
             x="Sample", y="Relative abundance (%)") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
              legend.text=element_text(size=8),
              legend.key.size=unit(0.4,"cm"))
      fname <- sprintf("04_taxonomy_%s.pdf", tolower(rank))
      save_pdf(p, fname, width=max(8, nsamples(ps)*0.45), height=7)

      # Save table
      tax_tbl <- melt_df %>%
        dplyr::select(Sample, !!rank, Abundance) %>%
        tidyr::pivot_wider(names_from="Sample", values_from="Abundance",
                           values_fn=mean, values_fill=0)
      write.csv(tax_tbl, file.path(TABLES_DIR, sprintf("taxonomy_%s.csv", tolower(rank))),
                row.names=FALSE)
    }, error=function(e) cat(sprintf("  [WARN] %s barplot: %s\n", rank, e$message)))
  }
}, error=function(e) cat(sprintf("[WARN] Taxonomy barplots: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Beta diversity
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 5: Beta Diversity ─────────────────────────────────\n")
tryCatch({
  if (!has_ggplot2) stop("ggplot2 not available")

  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

  # ── Bray-Curtis PCoA ──
  dist_bc   <- phyloseq::distance(ps_rel, method="bray")
  ord_pcoa  <- ordinate(ps_rel, method="PCoA", distance=dist_bc)

  eig   <- ord_pcoa$values$Eigenvalues
  var_e <- round(eig / sum(abs(eig)) * 100, 1)
  pcoa_df <- as.data.frame(ord_pcoa$vectors[, 1:2])
  colnames(pcoa_df) <- c("PC1","PC2")
  pcoa_df$Sample <- rownames(pcoa_df)
  if (has_meta) pcoa_df[[GROUP_COL]] <- meta_df[pcoa_df$Sample, GROUP_COL]

  n_grp <- if (has_meta) length(unique(pcoa_df[[GROUP_COL]])) else nsamples(ps)
  pal   <- make_palette(n_grp)

  p_pcoa <- ggplot(pcoa_df, aes_string(
      x="PC1", y="PC2",
      color=if (has_meta) GROUP_COL else "Sample",
      label="Sample")) +
    geom_point(size=3.5, alpha=0.85) +
    (if (requireNamespace("ggrepel", quietly=TRUE))
       ggrepel::geom_text_repel(size=2.8, show.legend=FALSE)
     else geom_text(vjust=-1, size=2.8, show.legend=FALSE)) +
    scale_color_manual(values=pal) +
    labs(title="Beta Diversity — Bray-Curtis PCoA",
         x=sprintf("PC1 [%.1f%%]", var_e[1]),
         y=sprintf("PC2 [%.1f%%]", var_e[2])) +
    stat_ellipse(aes_string(group=if (has_meta) GROUP_COL else "Sample"),
                 type="t", linetype=2, show.legend=FALSE) +
    theme_bw()
  save_pdf(p_pcoa, "05_beta_PCoA_BrayCurtis.pdf")

  # ── NMDS ──
  tryCatch({
    set.seed(42)
    ord_nmds <- ordinate(ps_rel, method="NMDS", distance=dist_bc)
    nmds_df  <- as.data.frame(ord_nmds$points)
    colnames(nmds_df) <- c("NMDS1","NMDS2")
    nmds_df$Sample <- rownames(nmds_df)
    if (has_meta) nmds_df[[GROUP_COL]] <- meta_df[nmds_df$Sample, GROUP_COL]

    p_nmds <- ggplot(nmds_df, aes_string(
        x="NMDS1", y="NMDS2",
        color=if (has_meta) GROUP_COL else "Sample")) +
      geom_point(size=3.5, alpha=0.85) +
      scale_color_manual(values=pal) +
      annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
               label=sprintf("stress = %.4f", ord_nmds$stress), size=3.5) +
      labs(title="Beta Diversity — NMDS (Bray-Curtis)") +
      theme_bw()
    save_pdf(p_nmds, "05_beta_NMDS.pdf")
  }, error=function(e) cat(sprintf("  [WARN] NMDS: %s\n", e$message)))

  # ── Weighted UniFrac PCoA (if tree available) ──
  if (has_tree) {
    tryCatch({
      dist_wuf <- phyloseq::distance(ps_rel, method="wunifrac")
      ord_wuf  <- ordinate(ps_rel, method="PCoA", distance=dist_wuf)
      eig2 <- ord_wuf$values$Eigenvalues
      var2 <- round(eig2 / sum(abs(eig2)) * 100, 1)
      wuf_df <- as.data.frame(ord_wuf$vectors[, 1:2])
      colnames(wuf_df) <- c("PC1","PC2")
      wuf_df$Sample <- rownames(wuf_df)
      if (has_meta) wuf_df[[GROUP_COL]] <- meta_df[wuf_df$Sample, GROUP_COL]

      p_wuf <- ggplot(wuf_df, aes_string(
          x="PC1", y="PC2",
          color=if (has_meta) GROUP_COL else "Sample")) +
        geom_point(size=3.5, alpha=0.85) +
        scale_color_manual(values=pal) +
        stat_ellipse(aes_string(group=if (has_meta) GROUP_COL else "Sample"),
                     type="t", linetype=2, show.legend=FALSE) +
        labs(title="Beta Diversity — Weighted UniFrac PCoA",
             x=sprintf("PC1 [%.1f%%]", var2[1]),
             y=sprintf("PC2 [%.1f%%]", var2[2])) +
        theme_bw()
      save_pdf(p_wuf, "05_beta_PCoA_wUniFrac.pdf")
    }, error=function(e) cat(sprintf("  [WARN] wUniFrac: %s\n", e$message)))
  }

  # ── Heatmap (sample × sample distance) ──
  tryCatch({
    dist_mat <- as.matrix(dist_bc)
    if (requireNamespace("pheatmap", quietly=TRUE)) {
      annt <- if (has_meta && GROUP_COL %in% colnames(meta_df))
        data.frame(row.names=rownames(dist_mat),
                   Group=meta_df[rownames(dist_mat), GROUP_COL])
      else NULL
      pdf(file.path(PLOTS_DIR, "05_beta_heatmap.pdf"), width=9, height=8)
      pheatmap::pheatmap(1 - dist_mat,
                         annotation_col=annt,
                         annotation_row=annt,
                         color=colorRampPalette(c("white","#2563eb"))(100),
                         main="Sample similarity (1 − Bray-Curtis)")
      dev.off()
      cat(sprintf("  ✓ Saved: 05_beta_heatmap.pdf\n"))
    }
  }, error=function(e) cat(sprintf("  [WARN] Heatmap: %s\n", e$message)))

  # ── UPGMA dendrogram ──
  tryCatch({
    upgma_tree <- hclust(dist_bc, method="average")
    pdf(file.path(PLOTS_DIR, "05_beta_UPGMA.pdf"), width=10, height=6)
    plot(upgma_tree, main="UPGMA Dendrogram (Bray-Curtis)",
         xlab="", sub="", cex=0.85)
    dev.off()
    cat(sprintf("  ✓ Saved: 05_beta_UPGMA.pdf\n"))
  }, error=function(e) cat(sprintf("  [WARN] UPGMA: %s\n", e$message)))

  # ── PERMANOVA table ──
  if (has_meta && has_vegan) {
    tryCatch({
      group_vec <- meta_df[labels(dist_bc), GROUP_COL]
      pm <- adonis2(dist_bc ~ group_vec)
      write.csv(as.data.frame(pm),
                file.path(TABLES_DIR, "permanova_BrayCurtis.csv"))
      cat("  ✓ Saved: permanova_BrayCurtis.csv\n")
    }, error=function(e) cat(sprintf("  [WARN] PERMANOVA: %s\n", e$message)))
  }

}, error=function(e) cat(sprintf("[WARN] Beta diversity: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Differential Abundance (ANCOMBC2)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 6: Differential Abundance (ANCOMBC2) ──────────────\n")
if (has_ancombc && has_meta && has_ggplot2) {
  tryCatch({
    # Collapse to genus level if rank available
    ps_da <- ps
    if ("Genus" %in% rank_names(ps)) {
      ps_da <- tax_glom(ps, taxrank="Genus", NArm=FALSE)
    }

    # Run ANCOMBC2
    res <- ANCOMBC::ancombc2(
      data        = ps_da,
      assay_name  = "counts",
      tax_level   = NULL,
      fix_formula = GROUP_COL,
      rand_formula = NULL,
      p_adj_method = "BH",
      prv_cut      = 0.10,
      lib_cut      = 1000,
      s0_perc      = 0.05,
      group        = GROUP_COL,
      struc_zero   = FALSE,
      neg_lb       = FALSE,
      alpha        = 0.05,
      n_cl         = 1,
      verbose      = FALSE
    )

    res_df  <- res$res
    sig_df  <- res_df[res_df$diff_abn == TRUE, ]

    write.csv(res_df, file.path(TABLES_DIR, "ancombc2_all.csv"),    row.names=FALSE)
    write.csv(sig_df, file.path(TABLES_DIR, "ancombc2_sig.csv"),    row.names=FALSE)
    cat(sprintf("  ANCOMBC2: %d significant taxa (p_adj < 0.05)\n", nrow(sig_df)))

    # Volcano plot
    lfc_col  <- grep("^lfc_",  colnames(res_df), value=TRUE)[1]
    padj_col <- grep("^q_",    colnames(res_df), value=TRUE)[1]

    if (!is.na(lfc_col) && !is.na(padj_col)) {
      res_df$neg_log10_q <- -log10(res_df[[padj_col]] + 1e-10)
      res_df$Significant <- res_df[[padj_col]] < 0.05

      p_volc <- ggplot(res_df, aes_string(x=lfc_col, y="neg_log10_q",
                                           color="Significant")) +
        geom_point(alpha=0.7, size=2) +
        scale_color_manual(values=c("grey60","#ef4444")) +
        geom_hline(yintercept=-log10(0.05), linetype=2, color="#6b7280") +
        geom_vline(xintercept=0, linetype=2, color="#6b7280") +
        labs(title="ANCOMBC2 — Differential Abundance",
             x="Log2 Fold Change", y="-log10(q-value)") +
        theme_bw()
      save_pdf(p_volc, "06_ancombc2_volcano.pdf")
    }

    # Bar plot — top significant taxa
    if (nrow(sig_df) > 0 && !is.na(lfc_col)) {
      top_sig <- head(sig_df[order(abs(sig_df[[lfc_col]]), decreasing=TRUE), ], 30)
      p_bar <- ggplot(top_sig, aes_string(
          x=sprintf("reorder(taxon, %s)", lfc_col),
          y=lfc_col,
          fill=sprintf("ifelse(%s > 0, 'Higher', 'Lower')", lfc_col))) +
        geom_col(alpha=0.85) +
        coord_flip() +
        scale_fill_manual(values=c("Higher"="#22c55e","Lower"="#ef4444"),
                          name="Direction") +
        labs(title="ANCOMBC2 — Top Differentially Abundant Taxa",
             x="Taxon", y="Log2 Fold Change") +
        theme_bw(base_size=10)
      save_pdf(p_bar, "06_ancombc2_barplot.pdf", width=10, height=max(6, nrow(top_sig)*0.35))
    }

  }, error=function(e) cat(sprintf("[WARN] ANCOMBC2: %s\n", e$message)))
} else {
  cat("  Skipped (needs ANCOMBC + metadata with group column)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — FUNGuildR (ITS only)
# ═══════════════════════════════════════════════════════════════════════════════
if (is_ITS) {
  cat("\n── Section 7: FUNGuildR Ecological Guild Annotation ──────────\n")
  if (has_funguild && !is.null(tax_df) && has_ggplot2) {
    tryCatch({
      # FUNGuildR needs "Taxonomy" column with full lineage string
      if ("Taxon" %in% colnames(tax_df)) {
        fg_in <- data.frame(Taxonomy = tax_df$Taxon, stringsAsFactors=FALSE)
        fg_res <- FUNGuildR::funguild_assign(fg_in)

        if (!is.null(fg_res) && nrow(fg_res) > 0) {
          fg_res$Feature.ID <- tax_df[[1]]
          write.csv(fg_res, file.path(TABLES_DIR, "funguild_annotation.csv"),
                    row.names=FALSE)
          cat(sprintf("  FUNGuildR: annotated %d / %d ASVs\n",
                      sum(!is.na(fg_res$guild)), nrow(fg_res)))

          # Guild composition pie / bar
          guild_counts <- sort(table(fg_res$guild[!is.na(fg_res$guild)]),
                               decreasing=TRUE)
          if (length(guild_counts) > 0) {
            gdf <- data.frame(Guild=names(guild_counts), n=as.integer(guild_counts))
            p_guild <- ggplot(gdf, aes(x=reorder(Guild,-n), y=n, fill=Guild)) +
              geom_col(alpha=0.85, show.legend=FALSE) +
              labs(title="FUNGuild — Ecological Guild Counts",
                   x="Guild", y="Number of ASVs") +
              theme_bw() +
              theme(axis.text.x=element_text(angle=40, hjust=1, size=9))
            save_pdf(p_guild, "07_funguild_barplot.pdf", width=10, height=6)
          }
        }
      }
    }, error=function(e) cat(sprintf("[WARN] FUNGuildR: %s\n", e$message)))
  } else {
    cat("  Skipped (FUNGuildR or ggplot2 not available)\n")
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Export final tables
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 8: Export Final Tables ───────────────────────────\n")
tryCatch({
  # ASV table (raw counts)
  otu_df <- as.data.frame(as.matrix(otu_table(ps)))
  if (taxa_are_rows(ps)) otu_df <- t(otu_df)
  write.csv(otu_df, file.path(TABLES_DIR, "asv_table.csv"))
  cat("  ✓ Saved: asv_table.csv\n")

  # Taxonomy table
  if (!is.null(access(ps, "tax_table"))) {
    tax_out <- as.data.frame(as.matrix(tax_table(ps)))
    write.csv(tax_out, file.path(TABLES_DIR, "taxonomy_table.csv"))
    cat("  ✓ Saved: taxonomy_table.csv\n")
  }
}, error=function(e) cat(sprintf("[WARN] Export tables: %s\n", e$message)))

# ─── Summary ──────────────────────────────────────────────────────────────────
plots_made <- list.files(PLOTS_DIR, pattern="\\.pdf$")
tables_made <- list.files(TABLES_DIR, pattern="\\.csv$")

cat("\n============================================================\n")
cat(sprintf("  R Visualization complete!\n"))
cat(sprintf("  PDFs created : %d  →  %s\n", length(plots_made),  PLOTS_DIR))
cat(sprintf("  CSVs created : %d  →  %s\n", length(tables_made), TABLES_DIR))
cat("============================================================\n\n")
