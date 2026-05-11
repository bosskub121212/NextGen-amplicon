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
has_lefser   <- load_pkg("lefser")

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

  # ── Pairwise Wilcoxon + BH correction ──
  if (has_meta && length(unique(alpha_df[[GROUP_COL]])) >= 2) {
    tryCatch({
      grps_w <- factor(alpha_df[[GROUP_COL]])
      wlx_list <- lapply(c("Shannon","Observed","Chao1","Simpson"), function(metric) {
        if (!metric %in% colnames(alpha_df)) return(NULL)
        vals <- setNames(alpha_df[[metric]], alpha_df$Sample)
        tryCatch({
          wt    <- pairwise.wilcox.test(vals, grps_w, p.adjust.method="BH")
          p_mat <- wt$p.value
          rows  <- rownames(p_mat); cols <- colnames(p_mat)
          do.call(rbind, lapply(seq_along(rows), function(i)
            do.call(rbind, lapply(seq_along(cols), function(j) {
              v <- p_mat[i,j]; if (is.na(v)) return(NULL)
              data.frame(metric=metric, group1=rows[i], group2=cols[j],
                         p.adj=round(v,6), stringsAsFactors=FALSE)
            }))
          ))
        }, error=function(e) NULL)
      })
      wlx_all <- do.call(rbind, Filter(Negate(is.null), wlx_list))
      if (!is.null(wlx_all) && nrow(wlx_all) > 0)
        write.csv(wlx_all, file.path(TABLES_DIR, "alpha_diversity_stats.csv"), row.names=FALSE)
      cat("  ✓ Saved: alpha_diversity_stats.csv\n")
    }, error=function(e) cat(sprintf("  [WARN] Wilcoxon alpha: %s\n", e$message)))
  }

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
# SECTION 3b — Shannon Rarefaction Curve
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 3b: Shannon Rarefaction Curve ──────────────────────\n")
if (has_vegan && has_ggplot2) {
  tryCatch({
    otu_rr    <- as.matrix(otu_table(ps))
    if (taxa_are_rows(ps)) otu_rr <- t(otu_rr)
    max_d_sh  <- min(max(rowSums(otu_rr)), 30000)
    steps_sh  <- unique(c(seq(100, max_d_sh, by=max(200, max_d_sh %/% 30)), max_d_sh))

    shan_list <- lapply(seq_len(nrow(otu_rr)), function(i) {
      sname <- rownames(otu_rr)[i]
      x     <- otu_rr[i, ]
      valid <- steps_sh[steps_sh <= sum(x)]
      if (length(valid) < 2) return(NULL)
      do.call(rbind, lapply(valid, function(s) {
        set.seed(42)
        raref <- vegan::rrarefy(matrix(x, nrow=1), s)[1, ]
        data.frame(Sample=sname, Depth=s,
                   Shannon=vegan::diversity(raref, "shannon"), stringsAsFactors=FALSE)
      }))
    })
    shan_df <- do.call(rbind, Filter(Negate(is.null), shan_list))

    if (!is.null(shan_df) && nrow(shan_df) > 0) {
      if (has_meta) shan_df[[GROUP_COL]] <- meta_df[shan_df$Sample, GROUP_COL]
      write.csv(shan_df, file.path(TABLES_DIR, "shannon_rarefaction.csv"), row.names=FALSE)
      n_s_sh <- length(unique(shan_df$Sample))
      pal_sh <- make_palette(n_s_sh)
      p_sh <- ggplot(shan_df, aes(x=Depth, y=Shannon, color=Sample)) +
        geom_line(size=0.8, alpha=0.85) +
        scale_color_manual(values=pal_sh) +
        labs(title="Shannon Diversity Rarefaction Curve",
             x="Sequencing depth", y="Shannon Index") +
        theme_bw() +
        theme(legend.text=element_text(size=8), legend.key.size=unit(0.5,"cm"))
      save_pdf(p_sh, "03b_shannon_rarefaction.pdf", width=10, height=6)
    }
  }, error=function(e) cat(sprintf("[WARN] Shannon rarefaction: %s\n", e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3c — Rank Abundance Curve
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 3c: Rank Abundance Curve ───────────────────────────\n")
if (has_ggplot2) {
  tryCatch({
    otu_ra <- as.matrix(otu_table(ps))
    if (taxa_are_rows(ps)) otu_ra <- t(otu_ra)

    ra_list <- lapply(seq_len(nrow(otu_ra)), function(i) {
      x   <- sort(otu_ra[i, otu_ra[i,] > 0], decreasing=TRUE)
      if (length(x) == 0) return(NULL)
      rel <- x / sum(x) * 100
      data.frame(Sample=rownames(otu_ra)[i], Rank=seq_along(rel),
                 RelAbundance=as.numeric(rel), stringsAsFactors=FALSE)
    })
    ra_df <- do.call(rbind, Filter(Negate(is.null), ra_list))

    if (!is.null(ra_df) && nrow(ra_df) > 0) {
      if (has_meta) ra_df[[GROUP_COL]] <- meta_df[ra_df$Sample, GROUP_COL]
      write.csv(ra_df, file.path(TABLES_DIR, "rank_abundance.csv"), row.names=FALSE)
      n_ra  <- length(unique(ra_df$Sample))
      pal_ra <- make_palette(n_ra)
      p_ra <- ggplot(ra_df, aes(x=Rank, y=RelAbundance, color=Sample)) +
        geom_line(size=0.8, alpha=0.85) +
        scale_y_log10() +
        scale_color_manual(values=pal_ra) +
        labs(title="Rank Abundance Curve",
             x="ASV Rank", y="Relative Abundance (%, log scale)") +
        theme_bw() +
        theme(legend.text=element_text(size=8), legend.key.size=unit(0.5,"cm"))
      save_pdf(p_ra, "03c_rank_abundance.pdf", width=10, height=6)
    }
  }, error=function(e) cat(sprintf("[WARN] Rank abundance: %s\n", e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3d — Species Accumulation Curve (specaccum)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 3d: Species Accumulation Curve ─────────────────────\n")
if (has_vegan && has_ggplot2 && nsamples(ps) >= 3) {
  tryCatch({
    otu_sp <- as.matrix(otu_table(ps))
    if (taxa_are_rows(ps)) otu_sp <- t(otu_sp)
    sp_acc <- vegan::specaccum(otu_sp, method="random", permutations=100)
    sp_df  <- data.frame(Sites=sp_acc$sites, Richness=sp_acc$richness, SD=sp_acc$sd)
    write.csv(sp_df, file.path(TABLES_DIR, "specaccum.csv"), row.names=FALSE)
    p_sp <- ggplot(sp_df, aes(x=Sites, y=Richness)) +
      geom_ribbon(aes(ymin=Richness-SD, ymax=Richness+SD), fill="#3b82f6", alpha=0.2) +
      geom_line(color="#3b82f6", size=1) +
      geom_point(color="#3b82f6", size=2) +
      labs(title="Species Accumulation Curve",
           x="Number of Samples", y="Cumulative ASV Richness") +
      theme_bw()
    save_pdf(p_sp, "03d_specaccum.pdf", width=8, height=5)
  }, error=function(e) cat(sprintf("[WARN] Specaccum: %s\n", e$message)))
}

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

# ── Group-level taxonomy bar charts ──────────────────────────────────────────
if (has_meta && has_ggplot2) {
  cat("  Generating group-level taxonomy bar charts...\n")
  tax_ranks_avail2 <- rank_names(ps)
  plot_ranks2 <- intersect(c("Phylum","Class","Order","Family","Genus"), tax_ranks_avail2)
  for (rank2 in plot_ranks2) {
    tryCatch({
      ps_glom2 <- tax_glom(ps, taxrank=rank2, NArm=FALSE)
      ps_rel2  <- transform_sample_counts(ps_glom2, function(x) x/sum(x)*100)
      melt2    <- psmelt(ps_rel2)
      melt2[[rank2]][is.na(melt2[[rank2]])] <- "Unclassified"
      melt2[[GROUP_COL]] <- meta_df[melt2$Sample, GROUP_COL]

      # Aggregate by group (mean)
      group_melt2 <- melt2 %>%
        dplyr::group_by(.data[[GROUP_COL]], .data[[rank2]]) %>%
        dplyr::summarise(Abundance=mean(Abundance, na.rm=TRUE), .groups="drop")

      top_t2 <- group_melt2 %>%
        dplyr::group_by(.data[[rank2]]) %>%
        dplyr::summarise(m=mean(Abundance)) %>%
        dplyr::arrange(dplyr::desc(m)) %>%
        dplyr::slice_head(n=15) %>%
        dplyr::pull(.data[[rank2]])
      group_melt2$TaxLabel <- ifelse(group_melt2[[rank2]] %in% top_t2,
                                      group_melt2[[rank2]], "Other")
      n_col2 <- length(unique(group_melt2$TaxLabel))
      pal2   <- c(make_palette(min(15, n_col2-1)), "grey80")[seq_len(n_col2)]

      p2 <- ggplot(group_melt2, aes_string(x=GROUP_COL, y="Abundance", fill="TaxLabel")) +
        geom_bar(stat="identity", width=0.7) +
        scale_fill_manual(values=pal2, name=rank2) +
        labs(title=sprintf("Group-level Relative Abundance — %s", rank2),
             x="Group", y="Mean Relative Abundance (%)") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=30, hjust=1),
              legend.text=element_text(size=8))
      save_pdf(p2, sprintf("04b_taxonomy_group_%s.pdf", tolower(rank2)),
               width=max(6, length(unique(group_melt2[[GROUP_COL]]))*1.5+3), height=7)
    }, error=function(e) cat(sprintf("  [WARN] Group bar %s: %s\n", rank2, e$message)))
  }
}

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
      group_vec_beta <- meta_df[labels(dist_bc), GROUP_COL]
      pm <- adonis2(dist_bc ~ group_vec_beta)
      write.csv(as.data.frame(pm),
                file.path(TABLES_DIR, "permanova_BrayCurtis.csv"))
      cat("  ✓ Saved: permanova_BrayCurtis.csv\n")
    }, error=function(e) cat(sprintf("  [WARN] PERMANOVA: %s\n", e$message)))
  }

  # ── Pairwise PERMANOVA ──
  if (has_meta && has_vegan) {
    tryCatch({
      gv_pw <- meta_df[sample_names(ps), GROUP_COL]
      grps_pw <- unique(gv_pw[!is.na(gv_pw) & nchar(gv_pw) > 0])
      if (length(grps_pw) >= 2) {
        pw_rows <- do.call(rbind, lapply(combn(grps_pw, 2, simplify=FALSE), function(pair) {
          subs <- names(gv_pw)[gv_pw %in% pair]
          if (length(subs) < 4) return(NULL)
          sub_d <- as.dist(as.matrix(dist_bc)[subs, subs])
          sub_g <- factor(gv_pw[subs])
          tryCatch({
            pm2 <- adonis2(sub_d ~ sub_g, permutations=999)
            data.frame(group1=pair[1], group2=pair[2],
                       R2=round(pm2$R2[1],4), F.stat=round(pm2$F[1],4),
                       p=pm2$`Pr(>F)`[1], p.adj=NA_real_, stringsAsFactors=FALSE)
          }, error=function(e) NULL)
        }))
        if (!is.null(pw_rows) && nrow(pw_rows) > 0) {
          pw_rows$p.adj <- p.adjust(pw_rows$p, method="BH")
          write.csv(pw_rows, file.path(TABLES_DIR, "beta_pairwise_permanova.csv"), row.names=FALSE)
          cat("  ✓ Saved: beta_pairwise_permanova.csv\n")
        }
      }
    }, error=function(e) cat(sprintf("  [WARN] Pairwise PERMANOVA: %s\n", e$message)))
  }

  # ── ANOSIM (global + pairwise) ──
  if (has_vegan) {
    tryCatch({
      gv_an  <- if (has_meta) meta_df[sample_names(ps), GROUP_COL] else NULL
      if (!is.null(gv_an) && length(unique(gv_an)) >= 2) {
        grp_fac_an <- factor(gv_an)
        anosi_res  <- anosim(dist_bc, grp_fac_an, permutations=999)
        sink_path  <- file.path(TABLES_DIR, "beta_anosim.txt")
        sink(sink_path)
        cat("=== ANOSIM (Bray-Curtis) ===\n\n"); print(anosi_res)
        sink()
        grps_an <- unique(gv_an[!is.na(gv_an) & nchar(gv_an) > 0])
        pw_an   <- do.call(rbind, lapply(combn(grps_an, 2, simplify=FALSE), function(pair) {
          subs <- names(gv_an)[gv_an %in% pair]
          if (length(subs) < 4) return(NULL)
          sub_d <- as.dist(as.matrix(dist_bc)[subs, subs])
          sub_g <- factor(gv_an[subs])
          tryCatch({
            an2 <- anosim(sub_d, sub_g, permutations=999)
            data.frame(group1=pair[1], group2=pair[2],
                       R=round(an2$statistic,4), p=an2$signif, stringsAsFactors=FALSE)
          }, error=function(e) NULL)
        }))
        if (!is.null(pw_an) && nrow(pw_an) > 0) {
          pw_an$p.adj <- p.adjust(pw_an$p, method="BH")
          write(capture.output(print(pw_an)), sink_path, append=TRUE)
          write.csv(pw_an, file.path(TABLES_DIR, "beta_anosim_pairwise.csv"), row.names=FALSE)
          cat("  ✓ Saved: beta_anosim_pairwise.csv\n")
        }
        cat("  ✓ Saved: beta_anosim.txt\n")
      }
    }, error=function(e) {
      if (sink.number() > 0) sink()
      cat(sprintf("  [WARN] ANOSIM: %s\n", e$message))
    })
  }

  # ── Jaccard NMDS ──
  tryCatch({
    ps_rel_j  <- transform_sample_counts(ps, function(x) x / sum(x))
    dist_jacc <- phyloseq::distance(ps_rel_j, method="jaccard", binary=TRUE)
    set.seed(42)
    ord_jacc  <- ordinate(ps_rel_j, method="NMDS", distance=dist_jacc)
    jdf       <- as.data.frame(ord_jacc$points)
    colnames(jdf) <- c("NMDS1","NMDS2")
    jdf$Sample <- rownames(jdf)
    if (has_meta) jdf[[GROUP_COL]] <- meta_df[jdf$Sample, GROUP_COL]
    if (has_ggplot2) {
      n_grp_j <- if (has_meta) length(unique(jdf[[GROUP_COL]])) else nsamples(ps)
      pal_j   <- make_palette(n_grp_j)
      p_jacc <- ggplot(jdf, aes_string(
          x="NMDS1", y="NMDS2",
          color=if (has_meta) GROUP_COL else "Sample")) +
        geom_point(size=3.5, alpha=0.85) +
        scale_color_manual(values=pal_j) +
        annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
                 label=sprintf("stress = %.4f", ord_jacc$stress), size=3.5) +
        labs(title="Beta Diversity — NMDS (Binary Jaccard)") +
        theme_bw()
      save_pdf(p_jacc, "05_beta_NMDS_Jaccard.pdf")
    }
  }, error=function(e) cat(sprintf("  [WARN] Jaccard NMDS: %s\n", e$message)))

}, error=function(e) cat(sprintf("[WARN] Beta diversity: %s\n", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5b — PCA (Principal Component Analysis, CLR-transformed)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 5b: PCA ────────────────────────────────────────────\n")
if (has_ggplot2) {
  tryCatch({
    # CLR transform (pseudocount of 1 to handle zeros)
    ps_clr <- transform_sample_counts(ps, function(x) {
      xp <- x + 1
      log(xp / exp(mean(log(xp))))
    })
    mat_clr <- as.matrix(otu_table(ps_clr))
    if (taxa_are_rows(ps_clr)) mat_clr <- t(mat_clr)

    pca_res  <- prcomp(mat_clr, scale.=FALSE, center=TRUE)
    var_exp  <- summary(pca_res)$importance[2, ] * 100
    n_pc     <- min(3, ncol(pca_res$x))
    pca_sc   <- as.data.frame(pca_res$x[, seq_len(n_pc)])
    pca_sc$Sample <- rownames(pca_sc)
    if (has_meta) pca_sc[[GROUP_COL]] <- meta_df[pca_sc$Sample, GROUP_COL]

    n_g_pca <- if (has_meta) length(unique(pca_sc[[GROUP_COL]])) else nsamples(ps)
    pal_pca <- make_palette(n_g_pca)

    make_pca_plot <- function(xcol, ycol) {
      p <- ggplot(pca_sc, aes_string(
          x=xcol, y=ycol,
          color=if (has_meta) GROUP_COL else "Sample",
          label="Sample")) +
        geom_point(size=3.5, alpha=0.85) +
        (if (requireNamespace("ggrepel", quietly=TRUE))
           ggrepel::geom_text_repel(size=2.8, show.legend=FALSE)
         else geom_text(vjust=-1, size=2.8, show.legend=FALSE)) +
        stat_ellipse(aes_string(group=if (has_meta) GROUP_COL else "Sample"),
                     type="t", linetype=2, show.legend=FALSE) +
        scale_color_manual(values=pal_pca) +
        labs(title="PCA (CLR-transformed)",
             x=sprintf("%s [%.1f%%]", xcol, var_exp[as.integer(sub("PC","",xcol))]),
             y=sprintf("%s [%.1f%%]", ycol, var_exp[as.integer(sub("PC","",ycol))])) +
        theme_bw()
      p
    }

    save_pdf(make_pca_plot("PC1","PC2"), "05b_PCA_PC1_PC2.pdf")
    if (n_pc >= 3) {
      save_pdf(make_pca_plot("PC1","PC3"), "05b_PCA_PC1_PC3.pdf")
      save_pdf(make_pca_plot("PC2","PC3"), "05b_PCA_PC2_PC3.pdf")
    }

    # Scree plot
    n_scree <- min(15, length(var_exp))
    scree_df <- data.frame(PC=seq_len(n_scree), Variance=var_exp[seq_len(n_scree)])
    p_scree <- ggplot(scree_df, aes(x=PC, y=Variance)) +
      geom_col(fill="#3b82f6", alpha=0.8) +
      geom_line(color="#1e40af", size=0.8) +
      geom_point(color="#1e40af", size=2) +
      scale_x_continuous(breaks=seq_len(n_scree)) +
      labs(title="PCA Scree Plot", x="Principal Component", y="Variance Explained (%)") +
      theme_bw()
    save_pdf(p_scree, "05b_PCA_scree.pdf", width=8, height=5)
    write.csv(pca_sc, file.path(TABLES_DIR, "pca_scores.csv"), row.names=FALSE)

  }, error=function(e) cat(sprintf("[WARN] PCA: %s\n", e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5c — Beta Distance Boxplot (within vs between groups)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 5c: Beta Distance Boxplot ──────────────────────────\n")
if (has_meta && has_ggplot2 && has_vegan) {
  tryCatch({
    ps_rel_bp  <- transform_sample_counts(ps, function(x) x/sum(x))
    dist_bc_bp <- phyloseq::distance(ps_rel_bp, method="bray")
    gv_bp      <- meta_df[sample_names(ps), GROUP_COL]
    dm_bp      <- as.matrix(dist_bc_bp)
    sn_bp      <- rownames(dm_bp)

    dist_pairs <- do.call(rbind, lapply(seq_len(nrow(dm_bp)-1), function(i) {
      do.call(rbind, lapply(seq(i+1, nrow(dm_bp)), function(j) {
        gi <- gv_bp[sn_bp[i]]; gj <- gv_bp[sn_bp[j]]
        data.frame(Distance=dm_bp[i,j],
                   Type=ifelse(gi==gj,"Within","Between"),
                   Group1=gi, Group2=gj, stringsAsFactors=FALSE)
      }))
    }))

    if (!is.null(dist_pairs) && nrow(dist_pairs) > 0) {
      write.csv(dist_pairs, file.path(TABLES_DIR, "beta_distance_pairs.csv"), row.names=FALSE)

      p_db <- ggplot(dist_pairs, aes(x=Type, y=Distance, fill=Type)) +
        geom_boxplot(alpha=0.7, outlier.shape=NA) +
        geom_jitter(width=0.12, size=1.5, alpha=0.5) +
        scale_fill_manual(values=c("Within"="#3b82f6","Between"="#ef4444")) +
        labs(title="Bray-Curtis Distance: Within vs Between Groups",
             x="", y="Bray-Curtis Distance") +
        theme_bw() + theme(legend.position="none")
      save_pdf(p_db, "05c_beta_distance_boxplot.pdf", width=5, height=6)

      # Per-comparison detail
      grp_pairs_bp <- unique(dist_pairs[dist_pairs$Type=="Between",
                                         c("Group1","Group2")])
      within_rows <- data.frame(
        Comparison=paste0("Within-", unique(gv_bp)),
        Type="Within", stringsAsFactors=FALSE
      )
      dist_pairs$Comparison <- ifelse(
        dist_pairs$Type=="Within",
        paste0("Within-", dist_pairs$Group1),
        paste(dist_pairs$Group1, "vs", dist_pairs$Group2)
      )
      p_db2 <- ggplot(dist_pairs, aes(x=Comparison, y=Distance, fill=Type)) +
        geom_boxplot(alpha=0.7, outlier.shape=NA) +
        geom_jitter(width=0.12, size=1.5, alpha=0.5) +
        scale_fill_manual(values=c("Within"="#3b82f6","Between"="#ef4444")) +
        labs(title="Bray-Curtis Distances by Group Comparison",
             x="", y="Bray-Curtis Distance") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="top")
      save_pdf(p_db2, "05c_beta_distance_groups.pdf", width=max(6, nrow(grp_pairs_bp)+3), height=6)
    }
  }, error=function(e) cat(sprintf("[WARN] Distance boxplot: %s\n", e$message)))
}

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
# SECTION 6b — DESeq2 Differential Abundance
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 6b: DESeq2 Differential Abundance ─────────────────\n")
if (has_meta && has_ggplot2 &&
    requireNamespace("DESeq2", quietly=TRUE) &&
    length(unique(meta_df[[GROUP_COL]])) >= 2 &&
    nsamples(ps) >= 4) {
  tryCatch({
    suppressPackageStartupMessages(library(DESeq2))

    ps_da <- ps
    if ("Genus" %in% rank_names(ps))
      ps_da <- tax_glom(ps, taxrank="Genus", NArm=FALSE)

    # Filter: remove features with >90% zeros
    ps_filt <- prune_taxa(
      rowSums(otu_table(ps_da) == 0) < ncol(otu_table(ps_da)) * 0.9,
      ps_da
    )
    if (ntaxa(ps_filt) < 2) stop("Too few taxa after filtering")

    cat(sprintf("  DESeq2: %d taxa × %d samples\n", ntaxa(ps_filt), nsamples(ps_filt)))
    ps_ds <- phyloseq_to_deseq2(ps_filt,
               as.formula(paste("~", GROUP_COL)))
    ds    <- DESeq2::estimateSizeFactors(ps_ds, type="poscounts")
    ds    <- DESeq2::DESeq(ds, test="Wald", fitType="parametric", quiet=TRUE)
    res   <- DESeq2::results(ds, alpha=0.05)
    res_df <- as.data.frame(res)
    res_df$taxon <- rownames(res_df)

    # Annotate labels
    if ("Genus" %in% rank_names(ps)) {
      lab <- as.character(tax_table(ps_filt)[res_df$taxon, "Genus"])
    } else {
      lab <- res_df$taxon
    }
    lab[is.na(lab) | lab == ""] <- res_df$taxon[is.na(lab) | lab == ""]
    res_df$label <- sub("^[a-z]__", "", lab)
    res_df       <- res_df[order(res_df$padj, na.last=NA), ]

    write.csv(res_df, file.path(TABLES_DIR, "deseq2_results.csv"), row.names=FALSE)
    cat("  ✓ Saved: deseq2_results.csv\n")
    cat(sprintf("  DESeq2: %d significant taxa (padj<0.05)\n",
                sum(!is.na(res_df$padj) & res_df$padj < 0.05)))

    # ── Volcano ──
    vdf <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    vdf$Significant  <- vdf$padj < 0.05
    vdf$neg_log10_p  <- -log10(vdf$padj + 1e-300)
    top_lab_vdf      <- head(vdf[vdf$Significant, ], 15)

    p_volc_ds <- ggplot(vdf, aes(x=log2FoldChange, y=neg_log10_p, colour=Significant)) +
      geom_point(alpha=0.7, size=2) +
      scale_colour_manual(values=c("grey60","#ef4444")) +
      geom_hline(yintercept=-log10(0.05), linetype=2, colour="#6b7280") +
      geom_vline(xintercept=0,            linetype=2, colour="#6b7280") +
      labs(title="DESeq2 — Differential Abundance",
           x="Log2 Fold Change", y="-log10(padj)") +
      theme_bw()
    if (requireNamespace("ggrepel", quietly=TRUE) && nrow(top_lab_vdf) > 0)
      p_volc_ds <- p_volc_ds +
        ggrepel::geom_text_repel(data=top_lab_vdf,
          aes(label=label), size=2.8, show.legend=FALSE)
    save_pdf(p_volc_ds, "06b_deseq2_volcano.pdf")

    # ── Heatmap of top 20 ──
    sig_taxa_ds <- head(res_df$taxon, 20)
    if (length(sig_taxa_ds) >= 2 && requireNamespace("pheatmap", quietly=TRUE)) {
      ps_rel_ds  <- transform_sample_counts(ps_da, function(x) x/sum(x)*100)
      ps_sig_ds  <- prune_taxa(intersect(sig_taxa_ds, taxa_names(ps_rel_ds)), ps_rel_ds)
      mat_ds     <- as.matrix(otu_table(ps_sig_ds))
      if (!taxa_are_rows(ps_sig_ds)) mat_ds <- t(mat_ds)

      if ("Genus" %in% rank_names(ps)) {
        rl <- sub("^[a-z]__", "",
                  as.character(tax_table(ps_sig_ds)[rownames(mat_ds), "Genus"]))
        rl[is.na(rl) | rl == ""] <- rownames(mat_ds)[is.na(rl) | rl == ""]
        rownames(mat_ds) <- rl
      }

      ann_col_ds <- data.frame(
        Group = meta_df[colnames(mat_ds), GROUP_COL],
        row.names = colnames(mat_ds)
      )
      grp_lvls_ds <- unique(ann_col_ds$Group)
      pal_ds      <- setNames(make_palette(length(grp_lvls_ds)), grp_lvls_ds)
      ann_colors_ds <- list(Group=pal_ds)

      if ("Phylum" %in% rank_names(ps)) {
        phy_vec_ds <- sub("^[a-z]__", "",
                          as.character(tax_table(ps_sig_ds)[, "Phylum"]))
        phy_vec_ds[is.na(phy_vec_ds)] <- "Unknown"
        phy_uniq_ds <- unique(phy_vec_ds)
        if (requireNamespace("RColorBrewer", quietly=TRUE)) {
          phy_pal_ds <- setNames(
            RColorBrewer::brewer.pal(max(3, min(length(phy_uniq_ds),12)), "Paired")[
              seq_len(length(phy_uniq_ds))],
            phy_uniq_ds)
        } else {
          phy_pal_ds <- setNames(make_palette(length(phy_uniq_ds)), phy_uniq_ds)
        }
        ann_row_ds          <- data.frame(Phylum=phy_vec_ds, row.names=rownames(mat_ds))
        ann_colors_ds$Phylum <- phy_pal_ds
      } else {
        ann_row_ds <- NULL
      }

      pdf(file.path(PLOTS_DIR, "06b_deseq2_heatmap.pdf"),
          width=max(7, nsamples(ps)*0.5+4),
          height=max(6, length(sig_taxa_ds)*0.4+3))
      pheatmap::pheatmap(
        mat_ds,
        scale             = "row",
        annotation_col    = ann_col_ds,
        annotation_row    = ann_row_ds,
        annotation_colors = ann_colors_ds,
        color             = colorRampPalette(c("#3b82f6","white","#ef4444"))(100),
        main              = "DESeq2 — Top 20 Differential Taxa (row-scaled %)",
        fontsize          = 8
      )
      dev.off()
      cat("  ✓ Saved: 06b_deseq2_heatmap.pdf\n")
    }

  }, error=function(e) cat(sprintf("[WARN] DESeq2: %s\n", e$message)))
} else {
  cat("  Skipped (needs DESeq2 + metadata + >=2 groups + >=4 samples)\n")
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

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Venn Diagram (shared ASVs between groups)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 9: Venn Diagram ────────────────────────────────────\n")
if (has_meta && has_phyloseq && length(unique(meta_df[[GROUP_COL]])) >= 2) {
  tryCatch({
    gv_venn   <- meta_df[sample_names(ps), GROUP_COL]
    grp_names <- unique(gv_venn[!is.na(gv_venn)])
    otu_mat_v <- as.matrix(otu_table(ps))
    if (!taxa_are_rows(ps)) otu_mat_v <- t(otu_mat_v)

    venn_sets <- lapply(setNames(grp_names, grp_names), function(g) {
      samps <- names(gv_venn)[!is.na(gv_venn) & gv_venn == g]
      if (length(samps) == 0) return(character(0))
      gs <- rowSums(otu_mat_v[, samps, drop=FALSE])
      names(gs)[gs > 0]
    })
    if (length(venn_sets) > 5) {
      venn_sets <- venn_sets[order(sapply(venn_sets,length), decreasing=TRUE)[1:5]]
    }
    n_v <- length(venn_sets)

    # Save membership table regardless of plot package
    all_t <- unique(unlist(venn_sets))
    if (length(all_t) > 0) {
      vmtbl <- do.call(cbind, lapply(venn_sets, function(s) as.integer(all_t %in% s)))
      rownames(vmtbl) <- all_t
      write.csv(as.data.frame(vmtbl),
                file.path(TABLES_DIR, "venn_membership.csv"))
      # Summary stats
      venn_summary <- data.frame(
        Group=names(venn_sets),
        N_unique=sapply(names(venn_sets), function(g) {
          others <- unlist(venn_sets[names(venn_sets) != g])
          sum(!(venn_sets[[g]] %in% others))
        }),
        N_total=sapply(venn_sets, length)
      )
      write.csv(venn_summary, file.path(TABLES_DIR, "venn_summary.csv"), row.names=FALSE)
    }

    if (n_v >= 2 && requireNamespace("ggVennDiagram", quietly=TRUE)) {
      suppressPackageStartupMessages(library(ggVennDiagram))
      p_venn <- ggVennDiagram(venn_sets, label_alpha=0, set_size=4) +
        scale_fill_gradient(low="white", high="#3b82f6") +
        scale_color_manual(values=rep("grey40", n_v)) +
        labs(title=sprintf("Shared ASVs between groups (%s)", GROUP_COL)) +
        theme(legend.position="right")
      save_pdf(p_venn, "09_venn_diagram.pdf", width=8, height=7)
    } else if (n_v >= 2 && has_ggplot2) {
      # Fallback: UpSet-style bar chart
      if (length(all_t) > 0) {
        cat("  [INFO] ggVennDiagram not installed — using overlap bar chart\n")
        cat("         Install with: install.packages('ggVennDiagram')\n")
        # Show unique + shared counts as a bar chart
        vc <- venn_summary
        p_vbar <- ggplot(vc, aes(x=reorder(Group,-N_unique), y=N_unique, fill=Group)) +
          geom_col(alpha=0.8) +
          labs(title="Unique ASVs per Group (Venn summary)",
               x="Group", y="# Unique ASVs") +
          theme_bw() + theme(legend.position="none")
        save_pdf(p_vbar, "09_venn_unique_bar.pdf", width=7, height=5)
      }
    }
  }, error=function(e) cat(sprintf("[WARN] Venn diagram: %s\n", e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — ANOVA per taxonomic level (+ Tukey post-hoc)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 10: ANOVA per Taxonomic Level ──────────────────────\n")
if (has_meta && has_ggplot2 && has_phyloseq &&
    length(unique(meta_df[[GROUP_COL]])) >= 2) {
  tryCatch({
    anova_ranks <- intersect(c("Phylum","Class","Order","Family","Genus"), rank_names(ps))
    for (rank_a in anova_ranks) {
      tryCatch({
        ps_glom_a <- tax_glom(ps, taxrank=rank_a, NArm=FALSE)
        ps_rel_a  <- transform_sample_counts(ps_glom_a, function(x) x/sum(x)*100)
        melt_a    <- psmelt(ps_rel_a)
        melt_a[[rank_a]][is.na(melt_a[[rank_a]])] <- "Unclassified"
        melt_a[[GROUP_COL]] <- meta_df[melt_a$Sample, GROUP_COL]

        taxa_a <- unique(melt_a[[rank_a]])
        # Keep only taxa with ≥0.1% mean in at least one group
        taxa_a <- taxa_a[sapply(taxa_a, function(t) {
          st  <- melt_a[melt_a[[rank_a]] == t, ]
          any(tapply(st$Abundance, st[[GROUP_COL]], mean, na.rm=TRUE) >= 0.1)
        })]
        if (length(taxa_a) == 0) next

        anova_rows <- do.call(rbind, lapply(taxa_a, function(t) {
          st  <- melt_a[melt_a[[rank_a]] == t, ]
          grp <- st[[GROUP_COL]]; abu <- st$Abundance
          if (length(unique(grp)) < 2 || any(table(grp) < 2)) return(NULL)
          tryCatch({
            fit   <- aov(abu ~ grp)
            p_val <- summary(fit)[[1]]$`Pr(>F)`[1]
            means <- tapply(abu, grp, mean, na.rm=TRUE)
            row   <- data.frame(taxon=t, p_value=round(p_val,6), stringsAsFactors=FALSE)
            for (g in names(means)) row[[paste0("mean_",g)]] <- round(means[g],4)
            row
          }, error=function(e) NULL)
        }))

        if (!is.null(anova_rows) && nrow(anova_rows) > 0) {
          anova_rows$p_adj <- p.adjust(anova_rows$p_value, method="BH")
          anova_rows <- anova_rows[order(anova_rows$p_adj), ]
          write.csv(anova_rows,
                    file.path(TABLES_DIR, sprintf("anova_%s.csv", tolower(rank_a))),
                    row.names=FALSE)

          sig_a <- anova_rows$taxon[!is.na(anova_rows$p_adj) & anova_rows$p_adj < 0.05]
          n_sig_a <- length(sig_a)
          cat(sprintf("  %s: %d significant taxa (padj<0.05)\n", rank_a, n_sig_a))

          if (n_sig_a >= 1) {
            top_a    <- head(sig_a, 12)
            plot_a   <- melt_a[melt_a[[rank_a]] %in% top_a, ]
            plot_a[[rank_a]] <- factor(plot_a[[rank_a]], levels=top_a)
            n_g_a   <- length(unique(plot_a[[GROUP_COL]]))
            pal_a   <- make_palette(n_g_a)

            p_an <- ggplot(plot_a, aes_string(x=GROUP_COL, y="Abundance", fill=GROUP_COL)) +
              geom_boxplot(alpha=0.7, outlier.shape=NA) +
              geom_jitter(width=0.15, size=1.5, alpha=0.7) +
              facet_wrap(as.formula(paste("~", rank_a)), scales="free_y", ncol=3) +
              scale_fill_manual(values=pal_a) +
              labs(title=sprintf("ANOVA significant %s (padj<0.05)", rank_a),
                   x=GROUP_COL, y="Relative Abundance (%)") +
              theme_bw(base_size=9) +
              theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none")
            save_pdf(p_an,
                     sprintf("10_anova_%s.pdf", tolower(rank_a)),
                     width=max(8, ceiling(length(top_a)/3)*3.5),
                     height=max(6, ceiling(length(top_a)/3)*3))
          }
        }
      }, error=function(e) cat(sprintf("  [WARN] ANOVA %s: %s\n", rank_a, e$message)))
    }
  }, error=function(e) cat(sprintf("[WARN] ANOVA section: %s\n", e$message)))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 11 — LEfSe Analysis (Linear Discriminant Analysis Effect Size)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 11: LEfSe Analysis ─────────────────────────────────\n")
if (has_meta && has_ggplot2 && has_phyloseq &&
    length(unique(meta_df[[GROUP_COL]])) >= 2 &&
    nsamples(ps) >= 4) {
  if (requireNamespace("lefser", quietly=TRUE) &&
      requireNamespace("SummarizedExperiment", quietly=TRUE)) {
    tryCatch({
      suppressPackageStartupMessages({
        library(lefser)
        library(SummarizedExperiment)
      })
      # Collapse to genus level if available
      ps_lef <- if ("Genus" %in% rank_names(ps))
                  tax_glom(ps, "Genus", NArm=FALSE) else ps

      otu_lef <- as.matrix(otu_table(ps_lef))
      if (!taxa_are_rows(ps_lef)) otu_lef <- t(otu_lef)
      # Relative abundance × 1e6 (CPM)
      otu_lef_rel <- sweep(otu_lef, 2, colSums(otu_lef), "/") * 1e6

      grp_lef <- factor(meta_df[colnames(otu_lef_rel), GROUP_COL])
      se_lef  <- SummarizedExperiment(
        assays  = list(counts=otu_lef_rel),
        colData = S4Vectors::DataFrame(group=grp_lef)
      )

      res_lef <- lefser(se_lef, groupCol="group", lda.threshold=2)

      if (!is.null(res_lef) && nrow(res_lef) > 0) {
        # Annotate with genus names
        if ("Genus" %in% rank_names(ps_lef)) {
          gmap <- sub("^[a-z]__","", as.character(tax_table(ps_lef)[,"Genus"]))
          names(gmap) <- taxa_names(ps_lef)
          res_lef$Label <- ifelse(res_lef$Names %in% names(gmap),
                                   gmap[res_lef$Names], res_lef$Names)
          res_lef$Label[is.na(res_lef$Label)|res_lef$Label==""] <-
            res_lef$Names[is.na(res_lef$Label)|res_lef$Label==""]
        } else {
          res_lef$Label <- res_lef$Names
        }

        grp_levels <- levels(grp_lef)
        res_lef$Direction <- ifelse(res_lef$scores > 0, grp_levels[2], grp_levels[1])
        res_lef <- res_lef[order(res_lef$scores), ]
        res_lef$Label <- factor(res_lef$Label, levels=res_lef$Label)
        write.csv(res_lef, file.path(TABLES_DIR, "lefse_results.csv"), row.names=FALSE)

        n_dir <- length(unique(res_lef$Direction))
        pal_lef <- setNames(make_palette(n_dir), unique(res_lef$Direction))

        p_lef <- ggplot(res_lef, aes(x=Label, y=scores, fill=Direction)) +
          geom_col(alpha=0.85) +
          coord_flip() +
          scale_fill_manual(values=pal_lef, name="Enriched in") +
          geom_hline(yintercept=0, color="grey50") +
          labs(title=sprintf("LEfSe — Significant Biomarkers (LDA≥2)\nGrouped by: %s", GROUP_COL),
               x="", y="LDA Score (log10)") +
          theme_bw(base_size=10) +
          theme(legend.position="top")
        save_pdf(p_lef, "11_lefse_lda_barplot.pdf",
                 width=9, height=max(5, nrow(res_lef)*0.32+2))
        cat(sprintf("  LEfSe: %d significant biomarkers\n", nrow(res_lef)))
      } else {
        cat("  LEfSe: no significant biomarkers found (try lowering lda.threshold)\n")
      }
    }, error=function(e) cat(sprintf("[WARN] LEfSe: %s\n", e$message)))
  } else {
    cat("  Skipped (install lefser: BiocManager::install('lefser'))\n")
  }
}

# ─── Summary ──────────────────────────────────────────────────────────────────
plots_made <- list.files(PLOTS_DIR, pattern="\\.pdf$")
tables_made <- list.files(TABLES_DIR, pattern="\\.csv$")

cat("\n============================================================\n")
cat(sprintf("  R Visualization complete!\n"))
cat(sprintf("  PDFs created : %d  →  %s\n", length(plots_made),  PLOTS_DIR))
cat(sprintf("  CSVs created : %d  →  %s\n", length(tables_made), TABLES_DIR))
cat("============================================================\n\n")
