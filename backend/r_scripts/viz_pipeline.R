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
              help="Pipeline output directory"),
  make_option("--input_dir",  type="character", default=NULL,
              help="Alias for --output_dir (Emu pipeline compatibility)"),
  make_option("--metadata",   type="character", default=NULL,
              help="Metadata file (TSV or CSV, optional)"),
  make_option("--group_col",  type="character", default="treatment",
              help="Metadata column to use for grouping"),
  make_option("--marker",     type="character", default="16S",
              help="Marker type: 16S|12S|ITS1|ITS2|COX1|18S-nema|PacBio|ONT-16S"),
  make_option("--topN",       type="integer",   default=30,
              help="Top N taxa to display [default: 30]"),
  make_option("--threads",    type="integer",   default=4,
              help="Threads (passed through)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# --input_dir is an alias for --output_dir (emu_pipeline.py compatibility)
if (is.null(opt$output_dir) && !is.null(opt$input_dir)) {
  opt$output_dir <- opt$input_dir
}

if (is.null(opt$output_dir)) {
  cat("ERROR: --output_dir is required\n")
  quit(status=1)
}

OUTPUT_DIR  <- normalizePath(opt$output_dir, mustWork=FALSE)
MARKER      <- opt$marker
GROUP_COL   <- opt$group_col
METADATA_FILE <- opt$metadata
TOP_N       <- if (!is.null(opt$topN) && opt$topN > 0) opt$topN else 30

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
}
has_pheatmap <- load_pkg("pheatmap")

# ── Color palette helper ──────────────────────────────────────────────────────
make_palette <- function(n) {
  if (n <= 8 && requireNamespace("RColorBrewer", quietly=TRUE)) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else {
    scales::hue_pal()(n)
  }
}

# ── Detect pipeline mode: Emu (CSV) vs QIIME2/DADA2 (BIOM) ───────────────────
emu_asv_file <- file.path(OUTPUT_DIR, "asv_table.csv")
emu_tax_file <- file.path(OUTPUT_DIR, "taxonomy.csv")
IS_EMU_MODE  <- file.exists(emu_asv_file) && file.exists(emu_tax_file) &&
                MARKER %in% c("ONT-16S", "ONT16S", "ONT")

biom_file  <- file.path(EXP_DIR, "feature-table", "feature-table.biom")
tax_file   <- file.path(EXP_DIR, "taxonomy",       "taxonomy.tsv")
tree_file  <- file.path(EXP_DIR, "tree",           "tree.nwk")
seq_fasta  <- file.path(EXP_DIR, "rep-seqs",       "dna-sequences.fasta")

ps      <- NULL
tax_df  <- NULL

if (IS_EMU_MODE) {
  # ── Emu mode: load asv_table.csv + taxonomy.csv ────────────────────────────
  cat("── Loading Emu output (CSV format) ───────────────────────────\n")

  tryCatch({
    # ASV table: rows = tax_id, cols = samples
    asv_raw <- read.csv(emu_asv_file, row.names=1, check.names=FALSE,
                        stringsAsFactors=FALSE)
    asv_mat <- as.matrix(asv_raw)
    storage.mode(asv_mat) <- "numeric"
    asv_mat[is.na(asv_mat)] <- 0
    # Remove all-zero rows
    asv_mat <- asv_mat[rowSums(asv_mat) > 0, , drop=FALSE]
    cat(sprintf("  Loaded ASV table: %d taxa × %d samples\n",
                nrow(asv_mat), ncol(asv_mat)))

    # Taxonomy table: rows = tax_id, cols = Kingdom..Species
    tax_raw <- read.csv(emu_tax_file, row.names=1, check.names=FALSE,
                        stringsAsFactors=FALSE)
    tax_cols <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                          colnames(tax_raw))
    tax_mat <- as.matrix(tax_raw[, tax_cols, drop=FALSE])
    tax_mat[is.na(tax_mat)] <- ""
    # Align rows
    common_ids <- intersect(rownames(asv_mat), rownames(tax_mat))
    asv_mat  <- asv_mat[common_ids, , drop=FALSE]
    tax_mat  <- tax_mat[common_ids, , drop=FALSE]
    cat(sprintf("  Loaded taxonomy: %d taxa\n", nrow(tax_mat)))

    if (has_phyloseq) {
      ps <- phyloseq(
        otu_table(asv_mat, taxa_are_rows=TRUE),
        tax_table(tax_mat)
      )
      cat(sprintf("  phyloseq object: %d ASVs × %d samples\n",
                  ntaxa(ps), nsamples(ps)))
    }
  }, error=function(e) cat(sprintf("[WARN] Emu CSV load error: %s\n", e$message)))

} else {
  # ── QIIME2/DADA2 mode: load BIOM ──────────────────────────────────────────
  cat("── Loading exported QIIME2 data ──────────────────────────────\n")

  if (!file.exists(biom_file)) {
    cat(sprintf("[ERROR] BIOM file not found: %s\n", biom_file))
    quit(status=1)
  }

  if (has_phyloseq) {
    tryCatch({
      ps <- import_biom(biom_file)
      cat(sprintf("  Loaded BIOM: %d ASVs × %d samples\n",
                  ntaxa(ps), nsamples(ps)))
    }, error=function(e) cat(sprintf("[WARN] BIOM import error: %s\n", e$message)))
  }

  # Load taxonomy TSV
  if (file.exists(tax_file)) {
    tryCatch({
      tax_df <- read.table(tax_file, header=TRUE, sep="\t",
                           comment.char="", quote="", stringsAsFactors=FALSE)
      cat(sprintf("  Loaded taxonomy: %d features\n", nrow(tax_df)))

      if (has_phyloseq && !is.null(ps)) {
        tax_col    <- if ("Taxon" %in% colnames(tax_df)) "Taxon" else colnames(tax_df)[2]
        tax_split  <- strsplit(tax_df[[tax_col]], ";\\s*")
        max_ranks  <- max(sapply(tax_split, length))
        ranks      <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[seq_len(min(7, max_ranks))]
        tax_mat2   <- do.call(rbind, lapply(tax_split, function(x) {
          x <- sub("^[a-z]__", "", x)
          x <- sub("^\\s+|\\s+$", "", x)
          length(x) <- length(ranks); x
        }))
        colnames(tax_mat2) <- ranks
        rownames(tax_mat2) <- tax_df[[1]]
        common <- intersect(taxa_names(ps), rownames(tax_mat2))
        if (length(common) > 0)
          tax_table(ps) <- tax_table(tax_mat2[common, , drop=FALSE])
      }
    }, error=function(e) cat(sprintf("[WARN] Taxonomy load error: %s\n", e$message)))
  }
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
has_tree  <- tryCatch(!is.null(phy_tree(ps)), error=function(e) FALSE)
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
    # Write WIDE format for Preview page: Depth | Sample1 | Sample2 | ...
    rar_wide <- tidyr::pivot_wider(rar_df[, c("Sample","Depth","ASVs")],
                                   names_from  = "Sample",
                                   values_from = "ASVs",
                                   values_fill = NA)
    rar_wide <- rar_wide[order(rar_wide$Depth), ]
    write.csv(rar_wide, file.path(TABLES_DIR, "rarefaction.csv"), row.names=FALSE)
    cat("  ✓ Saved: rarefaction.csv (wide format)\n")
    if (has_meta && has_ggplot2) {
      rar_df[[GROUP_COL]] <- meta_df[rar_df$Sample, GROUP_COL]
    }

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

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4b — Taxonomy Heatmaps (taxa × samples, clustered, top taxa)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 4b: Taxonomy Heatmaps ──────────────────────────────\n")
if (has_pheatmap && has_ggplot2 && has_dplyr && has_tidyr) {
  tax_ranks_avail3 <- rank_names(ps)
  hm_ranks <- intersect(c("Phylum","Family","Genus"), tax_ranks_avail3)

  ann_col_hm    <- NULL
  ann_colors_hm <- list()
  if (has_meta) {
    ann_col_hm <- data.frame(Group = factor(meta_df[sample_names(ps), GROUP_COL]),
                              row.names = sample_names(ps))
    colnames(ann_col_hm) <- GROUP_COL
    grp_levels_hm <- levels(ann_col_hm[[GROUP_COL]])
    grp_pal_hm    <- make_palette(length(grp_levels_hm))
    names(grp_pal_hm) <- grp_levels_hm
    ann_colors_hm[[GROUP_COL]] <- grp_pal_hm
  }

  for (hm_lvl in hm_ranks) {
    tryCatch({
      ps_glom_hm <- tax_glom(ps, taxrank=hm_lvl, NArm=FALSE)
      ps_rel_hm  <- transform_sample_counts(ps_glom_hm, function(x) x / sum(x) * 100)
      melt_hm    <- psmelt(ps_rel_hm)
      melt_hm[[hm_lvl]][is.na(melt_hm[[hm_lvl]])] <- "Unclassified"

      top_n_hm    <- min(20, length(unique(melt_hm[[hm_lvl]])))
      top_taxa_hm <- melt_hm %>%
        dplyr::group_by(.data[[hm_lvl]]) %>%
        dplyr::summarise(mean_abund = mean(Abundance, na.rm=TRUE), .groups="drop") %>%
        dplyr::arrange(dplyr::desc(mean_abund)) %>%
        dplyr::slice_head(n=top_n_hm) %>%
        dplyr::pull(.data[[hm_lvl]])

      wide_hm <- melt_hm %>%
        dplyr::filter(.data[[hm_lvl]] %in% top_taxa_hm) %>%
        dplyr::select(Sample, !!hm_lvl, Abundance) %>%
        tidyr::pivot_wider(names_from="Sample", values_from="Abundance",
                            values_fn=mean, values_fill=0)

      mat_hm <- as.matrix(wide_hm[, -1, drop=FALSE])
      rownames(mat_hm) <- wide_hm[[hm_lvl]]
      mat_hm <- mat_hm[order(rowMeans(mat_hm), decreasing=TRUE), , drop=FALSE]
      # Drop zero-variance rows (all-equal, e.g. all-zero) — row scaling would divide by 0
      mat_hm <- mat_hm[apply(mat_hm, 1, function(r) sd(r) > 0), , drop=FALSE]
      if (nrow(mat_hm) < 2) stop("fewer than 2 taxa with non-zero variance")

      hm_colors <- switch(hm_lvl,
        Genus  = colorRampPalette(c("#f0f4ff","#3b82f6","#1e1b4b"))(100),
        Family = colorRampPalette(c("#fff7ed","#f97316","#431407"))(100),
        Phylum = colorRampPalette(c("#f0fff4","#22c55e","#14532d"))(100)
      )
      hm_title <- sprintf("Top %d %s — Relative Abundance Heatmap", nrow(mat_hm), hm_lvl)
      hm_file  <- file.path(PLOTS_DIR, sprintf("04b_taxonomy_heatmap_%s.pdf", tolower(hm_lvl)))
      hm_w <- max(8, ncol(mat_hm) * 0.6 + 4)
      hm_h <- max(6, nrow(mat_hm) * 0.35 + 3)

      pheatmap::pheatmap(
        mat_hm,
        annotation_col    = ann_col_hm,
        annotation_colors = if (length(ann_colors_hm) > 0) ann_colors_hm else NULL,
        color             = hm_colors,
        scale             = "row",
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main              = hm_title,
        fontsize_row      = max(5, min(9, 200/nrow(mat_hm))),
        fontsize_col      = 8,
        filename          = hm_file,
        width             = hm_w,
        height            = hm_h
      )
      cat(sprintf("  ✓ Saved: %s\n", basename(hm_file)))

      write.csv(wide_hm, file.path(TABLES_DIR, sprintf("taxonomy_heatmap_%s.csv", tolower(hm_lvl))),
                row.names=FALSE)
    }, error=function(e) cat(sprintf("  [WARN] %s heatmap: %s\n", hm_lvl, e$message)))
  }
} else {
  cat("  Skipped (requires pheatmap + ggplot2 + dplyr + tidyr)\n")
}

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

  # Save PCoA scores CSV for Preview page
  pca_out <- pcoa_df[, c("Sample","PC1","PC2")]
  if (has_meta && GROUP_COL %in% colnames(pcoa_df)) pca_out$Group <- pcoa_df[[GROUP_COL]]
  pca_out$PC1_var <- round(var_e[1], 2)
  pca_out$PC2_var <- round(var_e[2], 2)
  write.csv(pca_out, file.path(TABLES_DIR, "pca_scores.csv"), row.names=FALSE)
  cat("  ✓ Saved: pca_scores.csv\n")

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
    # Save NMDS CSV for Edit Charts interactive view
    nmds_df$Stress <- round(ord_nmds$stress, 6)
    write.csv(nmds_df, file.path(TABLES_DIR, "nmds_bray.csv"), row.names=FALSE)
    cat("  ✓ Saved: nmds_bray.csv\n")
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
    # Save similarity matrix CSV for Preview page
    sim_mat <- 1 - dist_mat
    sim_df  <- as.data.frame(sim_mat)
    sim_df  <- cbind(Sample=rownames(sim_df), sim_df)
    write.csv(sim_df, file.path(TABLES_DIR, "beta_heatmap.csv"), row.names=FALSE)
    cat("  ✓ Saved: beta_heatmap.csv\n")

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
    # Save Jaccard NMDS CSV for Edit Charts interactive view
    jdf$Stress <- round(ord_jacc$stress, 6)
    write.csv(jdf, file.path(TABLES_DIR, "nmds_jaccard.csv"), row.names=FALSE)
    cat("  ✓ Saved: nmds_jaccard.csv\n")
    # Save Jaccard similarity matrix CSV (1 − Jaccard distance)
    jac_sim <- 1 - as.matrix(dist_jacc)
    jac_sim_df <- cbind(Sample=rownames(jac_sim), as.data.frame(jac_sim))
    write.csv(jac_sim_df, file.path(TABLES_DIR, "jaccard_heatmap.csv"), row.names=FALSE)
    cat("  ✓ Saved: jaccard_heatmap.csv\n")
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

    # Save pca_scores.csv with PC3 + variance % columns for Edit Charts
    pca_out <- pca_sc[, c("Sample", "PC1", "PC2")]
    if (n_pc >= 3 && "PC3" %in% colnames(pca_sc)) pca_out$PC3 <- pca_sc$PC3
    if (has_meta && GROUP_COL %in% colnames(pca_sc)) pca_out$Group <- pca_sc[[GROUP_COL]]
    pca_out$PC1_var <- round(var_exp[1], 2)
    pca_out$PC2_var <- round(var_exp[2], 2)
    if (n_pc >= 3) pca_out$PC3_var <- round(var_exp[3], 2)
    write.csv(pca_out, file.path(TABLES_DIR, "pca_scores.csv"), row.names=FALSE)
    cat("  ✓ Saved: pca_scores.csv\n")

    # Save PCA scree CSV
    write.csv(scree_df, file.path(TABLES_DIR, "pca_scree.csv"), row.names=FALSE)
    cat("  ✓ Saved: pca_scree.csv\n")

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

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — UPGMA Hierarchical Clustering Tree
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 12: UPGMA Clustering Tree ──────────────────────────\n")
if (has_phyloseq && has_vegan && has_ggplot2) {
  tryCatch({
    # Relative abundance matrix
    ps_upg  <- transform_sample_counts(ps, function(x) x / sum(x))
    otu_upg <- as.matrix(otu_table(ps_upg))
    if (taxa_are_rows(ps_upg)) otu_upg <- t(otu_upg)

    # Distance metrics to cluster
    dist_methods <- c("bray", "jaccard")
    for (dm in dist_methods) {
      tryCatch({
        d_upg  <- vegan::vegdist(otu_upg, method=dm)
        hc_upg <- hclust(d_upg, method="average")   # UPGMA = average linkage

        # Colour bars by group if metadata available
        if (has_meta) {
          grp_upg <- meta_df[hc_upg$labels, GROUP_COL]
          grp_upg[is.na(grp_upg)] <- "Unknown"
          pal_upg  <- make_palette(length(unique(grp_upg)))
          col_upg  <- setNames(pal_upg, unique(grp_upg))
          label_col <- col_upg[grp_upg]
        }

        fname <- sprintf("12_upgma_%s.pdf", dm)
        pdf(file.path(PLOTS_DIR, fname), width=max(8, nsamples(ps)*0.6), height=6)
        par(mar=c(5,4,3,2))
        plot(hc_upg,
             main=sprintf("UPGMA Clustering — %s distance", dm),
             xlab="", sub="", cex=0.85)
        if (has_meta) {
          rect.hclust(hc_upg, k=min(length(unique(grp_upg)), nsamples(ps)-1),
                      border=pal_upg)
        }
        dev.off()
        cat(sprintf("  ✓ Saved: %s\n", fname))
      }, error=function(e) cat(sprintf("  [WARN] UPGMA %s: %s\n", dm, e$message)))
    }
  }, error=function(e) cat(sprintf("[WARN] UPGMA section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + vegan + ggplot2)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 12b — ClusterTree + Bar (UPGMA dendrogram + Taxonomy stacked bar)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 12b: ClusterTree + Bar ─────────────────────────────\n")
if (has_phyloseq && has_vegan && has_ggplot2 &&
    requireNamespace("ggdendro", quietly=TRUE) &&
    requireNamespace("patchwork", quietly=TRUE) &&
    requireNamespace("tidyr",    quietly=TRUE)) {
  tryCatch({
    suppressPackageStartupMessages({
      library(ggdendro)
      library(patchwork)
      library(tidyr)
    })

    ps_ct  <- transform_sample_counts(ps, function(x) x / sum(x))
    otu_ct <- as.matrix(otu_table(ps_ct))
    if (taxa_are_rows(ps_ct)) otu_ct <- t(otu_ct)

    d_ct  <- vegan::vegdist(otu_ct, method="bray")
    hc_ct <- hclust(d_ct, method="average")
    sample_order <- hc_ct$labels[hc_ct$order]

    # ── Dendrogram ──────────────────────────────────────────────────────────
    dend_data <- ggdendro::dendro_data(hc_ct, type="rectangle")

    p_dend <- ggplot() +
      ggdendro::geom_segment(
        data = ggdendro::segment(dend_data),
        aes(x=x, y=y, xend=xend, yend=yend),
        linewidth=0.5, colour="#475569") +
      ggdendro::geom_text(
        data = ggdendro::label(dend_data),
        aes(x=x, y=y, label=label),
        hjust=1, size=2.8, colour="#1e293b") +
      coord_flip() +
      scale_y_reverse(expand=c(0.15, 0)) +
      theme_void(base_size=9) +
      theme(plot.margin=margin(4,0,4,4))

    # ── Stacked bar (Phylum) ────────────────────────────────────────────────
    tax_rank_ct <- if ("Phylum" %in% rank_names(ps)) "Phylum" else rank_names(ps)[min(2, length(rank_names(ps)))]
    ps_phy <- tryCatch(tax_glom(ps_ct, taxrank=tax_rank_ct, NArm=FALSE), error=function(e) NULL)

    if (!is.null(ps_phy)) {
      otu_bar <- as.data.frame(as.matrix(otu_table(ps_phy)))
      if (taxa_are_rows(ps_phy)) otu_bar <- as.data.frame(t(otu_bar))
      tax_names_ct <- sub("^[a-z]__", "",
                          as.character(tax_table(ps_phy)[, tax_rank_ct]))
      tax_names_ct[is.na(tax_names_ct) | tax_names_ct == ""] <- "Unknown"
      colnames(otu_bar) <- make.unique(tax_names_ct)
      otu_bar$Sample <- rownames(otu_bar)
      otu_bar$Sample <- factor(otu_bar$Sample, levels=sample_order)

      top10 <- names(sort(colSums(otu_bar[, -ncol(otu_bar)]), decreasing=TRUE))[
                 seq_len(min(10, ncol(otu_bar)-1))]
      df_melt <- tidyr::pivot_longer(otu_bar, cols=-Sample,
                                     names_to="Taxon", values_to="Abund")
      df_melt$Taxon <- ifelse(df_melt$Taxon %in% top10, df_melt$Taxon, "Other")
      df_melt$Taxon <- factor(df_melt$Taxon, levels=c(top10, "Other"))

      pal_ct <- c(make_palette(length(top10)), "#94a3b8")
      names(pal_ct) <- c(top10, "Other")

      p_bar <- ggplot(df_melt, aes(x=Sample, y=Abund, fill=Taxon)) +
        geom_bar(stat="identity", width=0.82) +
        scale_fill_manual(values=pal_ct) +
        scale_y_continuous(labels=scales::percent_format(accuracy=1),
                           expand=c(0,0)) +
        coord_flip() +
        labs(y="Relative Abundance", x=NULL,
             fill=tax_rank_ct,
             title=sprintf("ClusterTree — Bray-Curtis UPGMA + %s", tax_rank_ct)) +
        theme_bw(base_size=9) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.y=element_blank(),
              panel.grid.minor=element_blank(),
              legend.position="right",
              legend.key.size=unit(0.35,"cm"),
              legend.text=element_text(size=7.5),
              plot.title=element_text(size=9, face="bold"),
              plot.margin=margin(4,4,4,0))

      combined_ct <- p_dend + p_bar +
        patchwork::plot_layout(widths=c(1, 2.2))

      w_ct <- max(10, nsamples(ps) * 0.45 + 5)
      h_ct <- max(5,  nsamples(ps) * 0.28 + 2)
      fname_ct <- "12b_clustertree_bar.pdf"
      ggsave(file.path(PLOTS_DIR, fname_ct), plot=combined_ct,
             width=w_ct, height=h_ct, device="pdf")
      cat(sprintf("  ✓ Saved: %s\n", fname_ct))
    }
  }, error=function(e) cat(sprintf("  [WARN] ClusterTree+Bar: %s\n", e$message)))
} else {
  cat("  Skipped (requires ggdendro + patchwork + tidyr)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 13 — Metastats (T-test based differential abundance)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 13: Metastats (T-test) ─────────────────────────────\n")
if (has_phyloseq && has_meta && has_ggplot2) {
  tryCatch({
    grps_ms <- unique(meta_df[[GROUP_COL]])
    grp_pairs_ms <- combn(grps_ms, 2, simplify=FALSE)

    for (ranks_ms in c("Genus", "Family", "Phylum")) {
      if (!(ranks_ms %in% rank_names(ps))) next
      tryCatch({
        ps_ms  <- tax_glom(ps, taxrank=ranks_ms, NArm=FALSE)
        ps_ms  <- transform_sample_counts(ps_ms, function(x) x / sum(x) * 100)
        otu_ms <- as.matrix(otu_table(ps_ms))
        if (taxa_are_rows(ps_ms)) otu_ms <- t(otu_ms)

        # Get taxon names
        tax_ms <- sub("^[a-z]__", "",
                      as.character(tax_table(ps_ms)[, ranks_ms]))
        tax_ms[is.na(tax_ms) | tax_ms == ""] <- taxa_names(ps_ms)[is.na(tax_ms) | tax_ms == ""]
        colnames(otu_ms) <- tax_ms

        for (pair_ms in grp_pairs_ms) {
          g1 <- pair_ms[1]; g2 <- pair_ms[2]
          s1 <- rownames(meta_df)[meta_df[[GROUP_COL]] == g1]
          s2 <- rownames(meta_df)[meta_df[[GROUP_COL]] == g2]
          s1 <- intersect(s1, rownames(otu_ms))
          s2 <- intersect(s2, rownames(otu_ms))
          if (length(s1) < 2 || length(s2) < 2) next

          mat1 <- otu_ms[s1, , drop=FALSE]
          mat2 <- otu_ms[s2, , drop=FALSE]

          meta_res <- do.call(rbind, lapply(colnames(otu_ms), function(tx) {
            tryCatch({
              tt <- t.test(mat1[, tx], mat2[, tx], var.equal=FALSE)
              data.frame(
                taxon   = tx,
                mean_g1 = round(mean(mat1[, tx]), 4),
                mean_g2 = round(mean(mat2[, tx]), 4),
                p_value = tt$p.value,
                stringsAsFactors = FALSE
              )
            }, error=function(e) NULL)
          }))

          if (is.null(meta_res) || nrow(meta_res) == 0) next
          meta_res$p_adj <- p.adjust(meta_res$p_value, method="BH")
          meta_res       <- meta_res[order(meta_res$p_adj), ]
          colnames(meta_res)[2:3] <- c(paste0("mean_", g1), paste0("mean_", g2))

          csv_name <- sprintf("metastats_%s_%s_vs_%s.csv",
                              tolower(ranks_ms), g1, g2)
          write.csv(meta_res, file.path(TABLES_DIR, csv_name), row.names=FALSE)

          # Bar plot — top 20 significant (q≤0.05)
          sig_ms <- meta_res[!is.na(meta_res$p_adj) & meta_res$p_adj <= 0.05, ]
          if (nrow(sig_ms) == 0) {
            cat(sprintf("  %s %s vs %s: no significant taxa (q≤0.05)\n",
                        ranks_ms, g1, g2))
            next
          }
          top_ms <- head(sig_ms, 20)
          m1col  <- paste0("mean_", g1)
          m2col  <- paste0("mean_", g2)
          plot_ms <- rbind(
            data.frame(taxon=top_ms$taxon, mean=top_ms[[m1col]], group=g1,
                       stringsAsFactors=FALSE),
            data.frame(taxon=top_ms$taxon, mean=top_ms[[m2col]], group=g2,
                       stringsAsFactors=FALSE)
          )
          plot_ms$taxon <- factor(plot_ms$taxon,
                                  levels=rev(unique(top_ms$taxon)))

          p_ms <- ggplot(plot_ms, aes(x=taxon, y=mean, fill=group)) +
            geom_col(position="dodge", alpha=0.85) +
            coord_flip() +
            scale_fill_manual(values=make_palette(2)) +
            labs(title=sprintf("Metastats — %s: %s vs %s (q≤0.05)", ranks_ms, g1, g2),
                 x="", y="Mean Relative Abundance (%)", fill="Group") +
            theme_bw(base_size=10) +
            theme(legend.position="top")

          pdf_name <- sprintf("13_metastats_%s_%s_vs_%s.pdf",
                              tolower(ranks_ms), g1, g2)
          save_pdf(p_ms, pdf_name,
                   width=9, height=max(5, nrow(top_ms)*0.35+2))
          cat(sprintf("  %s %s vs %s: %d significant taxa\n",
                      ranks_ms, g1, g2, nrow(sig_ms)))
        }
      }, error=function(e) cat(sprintf("  [WARN] Metastats %s: %s\n", ranks_ms, e$message)))
    }
  }, error=function(e) cat(sprintf("[WARN] Metastats section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + metadata + ggplot2)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 14 — Ternary Plots (≥3 groups)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 14: Ternary Plots ───────────────────────────────────\n")
if (has_phyloseq && has_meta && has_ggplot2 &&
    requireNamespace("ggtern", quietly=TRUE)) {
  tryCatch({
    grps_tern <- unique(meta_df[[GROUP_COL]])
    if (length(grps_tern) < 3) {
      cat("  Skipped (need ≥3 groups for ternary plots)\n")
    } else {
      suppressPackageStartupMessages(library(ggtern))

      # Collapse to genus
      rank_tern <- if ("Genus" %in% rank_names(ps)) "Genus" else
                   if ("Family" %in% rank_names(ps)) "Family" else NULL

      if (!is.null(rank_tern)) {
        ps_tern  <- tax_glom(ps, taxrank=rank_tern, NArm=FALSE)
        ps_tern  <- transform_sample_counts(ps_tern, function(x) x / sum(x) * 100)
        otu_tern <- as.matrix(otu_table(ps_tern))
        if (taxa_are_rows(ps_tern)) otu_tern <- t(otu_tern)

        tax_tern <- sub("^[a-z]__", "",
                        as.character(tax_table(ps_tern)[, rank_tern]))
        tax_tern[is.na(tax_tern) | tax_tern == ""] <-
          taxa_names(ps_tern)[is.na(tax_tern) | tax_tern == ""]
        colnames(otu_tern) <- tax_tern

        # Compute mean per group
        grp_means_tern <- do.call(rbind, lapply(grps_tern, function(g) {
          samps_g <- intersect(rownames(meta_df)[meta_df[[GROUP_COL]] == g],
                               rownames(otu_tern))
          if (length(samps_g) == 0) return(NULL)
          colMeans(otu_tern[samps_g, , drop=FALSE])
        }))
        rownames(grp_means_tern) <- grps_tern

        # Top 20 taxa by total abundance
        top_tern <- names(sort(colSums(grp_means_tern), decreasing=TRUE))[1:min(20, ncol(grp_means_tern))]

        # All combinations of 3 groups
        group_combos <- combn(grps_tern, 3, simplify=FALSE)
        for (combo in group_combos[1:min(3, length(group_combos))]) {
          g1t <- combo[1]; g2t <- combo[2]; g3t <- combo[3]
          df_tern <- data.frame(
            taxon = top_tern,
            g1    = as.numeric(grp_means_tern[g1t, top_tern]),
            g2    = as.numeric(grp_means_tern[g2t, top_tern]),
            g3    = as.numeric(grp_means_tern[g3t, top_tern]),
            stringsAsFactors = FALSE
          )
          colnames(df_tern)[2:4] <- c(g1t, g2t, g3t)
          # Normalize rows to sum to 100
          row_sums_t <- rowSums(df_tern[,2:4])
          row_sums_t[row_sums_t == 0] <- 1
          df_tern[,2:4] <- df_tern[,2:4] / row_sums_t * 100

          p_tern <- ggtern(df_tern,
                           aes_string(x=g1t, y=g2t, z=g3t)) +
            geom_point(aes(colour=taxon), size=3, alpha=0.8) +
            theme_bw() +
            theme(legend.position="right",
                  legend.text=element_text(size=7)) +
            labs(title=sprintf("Ternary — %s/%s/%s (%s level)",
                               g1t, g2t, g3t, rank_tern),
                 colour=rank_tern)

          if (requireNamespace("ggrepel", quietly=TRUE)) {
            p_tern <- p_tern +
              ggrepel::geom_text_repel(aes(label=taxon), size=2.5,
                                       max.overlaps=8, show.legend=FALSE)
          }

          fname_tern <- sprintf("14_ternary_%s_%s_%s.pdf", g1t, g2t, g3t)
          save_pdf(p_tern, fname_tern, width=9, height=7)
        }
        cat(sprintf("  Ternary plots done (%d group combo(s))\n",
                    min(3, length(group_combos))))
      }
    }
  }, error=function(e) cat(sprintf("[WARN] Ternary section: %s\n", e$message)))
} else if (!requireNamespace("ggtern", quietly=TRUE)) {
  cat("  Skipped (install ggtern: install.packages('ggtern'))\n")
} else {
  cat("  Skipped (requires phyloseq + metadata + ≥3 groups)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 15 — Co-occurrence Network + ZiPi Plot
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 15: Co-occurrence Network + ZiPi ────────────────────\n")
if (has_phyloseq && has_ggplot2 &&
    requireNamespace("igraph",  quietly=TRUE) &&
    requireNamespace("Hmisc",   quietly=TRUE)) {
  tryCatch({
    suppressPackageStartupMessages({
      library(igraph)
      library(Hmisc)
    })

    # Collapse to genus, relative abundance
    ps_net <- if ("Genus" %in% rank_names(ps))
                tax_glom(ps, "Genus", NArm=FALSE) else ps
    ps_net <- transform_sample_counts(ps_net, function(x) x / sum(x) * 100)

    # Keep taxa present in ≥30% of samples, mean rel abund ≥0.1%
    ps_net <- filter_taxa(ps_net,
                          function(x) sum(x > 0) / length(x) >= 0.3 &
                                      mean(x) >= 0.1, TRUE)

    otu_net <- as.matrix(otu_table(ps_net))
    if (taxa_are_rows(ps_net)) otu_net <- t(otu_net)

    cat(sprintf("  Network taxa (filtered): %d\n", ncol(otu_net)))

    if (ncol(otu_net) >= 5 && nrow(otu_net) >= 5) {
      # Spearman correlation
      cor_res  <- Hmisc::rcorr(otu_net, type="spearman")
      r_mat    <- cor_res$r
      p_mat    <- cor_res$P
      diag(r_mat) <- 0

      # Threshold: |r| ≥ 0.6, p < 0.05
      adj_mat <- (abs(r_mat) >= 0.6) & (p_mat < 0.05) & !is.na(p_mat)
      diag(adj_mat) <- FALSE

      if (sum(adj_mat) > 0) {
        g_net <- igraph::graph_from_adjacency_matrix(
          adj_mat, mode="undirected", diag=FALSE)
        igraph::E(g_net)$weight <- r_mat[adj_mat]
        igraph::E(g_net)$color  <- ifelse(
          r_mat[adj_mat] > 0, "#3b82f6", "#ef4444")

        # Node properties
        deg_net   <- igraph::degree(g_net)
        btwn_net  <- igraph::betweenness(g_net, normalized=TRUE)

        # Taxon labels
        tax_lab_net <- sub("^[a-z]__", "",
          as.character(tax_table(ps_net)[igraph::V(g_net)$name,
            if ("Genus" %in% rank_names(ps_net)) "Genus" else rank_names(ps_net)[1]]))
        tax_lab_net[is.na(tax_lab_net) | tax_lab_net == ""] <-
          igraph::V(g_net)$name[is.na(tax_lab_net) | tax_lab_net == ""]
        igraph::V(g_net)$label <- tax_lab_net

        # Group colour by metadata if available
        if (has_meta && requireNamespace("ggraph", quietly=TRUE)) {
          suppressPackageStartupMessages(library(ggraph))
          set.seed(42)
          p_net <- ggraph(g_net, layout="fr") +
            geom_edge_link(aes(colour=color), alpha=0.5, width=0.7) +
            geom_node_point(aes(size=deg_net), colour="#1e40af", alpha=0.8) +
            geom_node_text(aes(label=label), repel=TRUE, size=2.5,
                           max.overlaps=15) +
            scale_edge_colour_identity() +
            scale_size_continuous(range=c(2, 8), name="Degree") +
            labs(title=sprintf("Co-occurrence Network (|r|≥0.6, p<0.05)\n%d nodes, %d edges",
                               igraph::vcount(g_net), igraph::ecount(g_net))) +
            theme_void(base_size=10) +
            theme(plot.title=element_text(hjust=0.5))
          save_pdf(p_net, "15a_cooccurrence_network.pdf", width=10, height=9)
        } else {
          # Fallback: base R plot
          pdf(file.path(PLOTS_DIR, "15a_cooccurrence_network.pdf"), width=10, height=9)
          igraph::plot.igraph(g_net,
            vertex.label      = tax_lab_net,
            vertex.label.cex  = 0.6,
            vertex.size       = sqrt(deg_net + 1) * 4,
            vertex.color      = "#93c5fd",
            edge.color        = igraph::E(g_net)$color,
            edge.width        = 1.5,
            layout            = igraph::layout_with_fr(g_net),
            main              = sprintf("Co-occurrence Network (%d nodes, %d edges)",
                                        igraph::vcount(g_net), igraph::ecount(g_net)))
          legend("bottomleft",
                 legend=c("Positive (r≥0.6)", "Negative (r≤-0.6)"),
                 col=c("#3b82f6","#ef4444"), lwd=2, bty="n", cex=0.8)
          dev.off()
          cat("  ✓ Saved: 15a_cooccurrence_network.pdf\n")
        }

        # Export edge table
        el_net <- igraph::as_data_frame(g_net, what="edges")
        el_net$r <- round(r_mat[cbind(el_net$from, el_net$to)], 4)
        write.csv(el_net, file.path(TABLES_DIR, "network_edges.csv"),
                  row.names=FALSE)

        # ── ZiPi plot ──────────────────────────────────────────────────────────
        # Module detection via fast_greedy
        comm_net <- igraph::cluster_fast_greedy(
          igraph::as.undirected(g_net))
        modules  <- igraph::membership(comm_net)

        node_df <- data.frame(
          taxon  = names(modules),
          module = as.integer(modules),
          degree = deg_net[names(modules)],
          stringsAsFactors = FALSE
        )

        # Within-module z-score (Zi)
        node_df$Zi <- sapply(seq_len(nrow(node_df)), function(i) {
          mod_nodes <- node_df$taxon[node_df$module == node_df$module[i]]
          mod_degs  <- node_df$degree[node_df$module == node_df$module[i]]
          if (length(mod_degs) < 2) return(0)
          (node_df$degree[i] - mean(mod_degs)) / sd(mod_degs)
        })
        node_df$Zi[is.nan(node_df$Zi) | is.na(node_df$Zi)] <- 0

        # Participation coefficient (Pi)
        node_df$Pi <- sapply(seq_len(nrow(node_df)), function(i) {
          nd   <- node_df$degree[i]
          if (nd == 0) return(0)
          nbrs <- names(igraph::neighbors(g_net, node_df$taxon[i]))
          km_s <- table(node_df$module[match(nbrs, node_df$taxon)])
          1 - sum((km_s / nd)^2)
        })

        # Classify roles
        node_df$Role <- with(node_df, ifelse(
          Zi >= 2.5 & Pi <  0.62, "Module hub",
          ifelse(Zi >= 2.5 & Pi >= 0.62, "Network hub",
          ifelse(Zi <  2.5 & Pi >= 0.62, "Connector",
                                          "Peripheral"))))

        write.csv(node_df, file.path(TABLES_DIR, "network_zipi.csv"),
                  row.names=FALSE)

        role_cols <- c("Peripheral"="#94a3b8","Connector"="#f59e0b",
                       "Module hub"="#3b82f6","Network hub"="#ef4444")

        p_zipi <- ggplot(node_df, aes(x=Pi, y=Zi, colour=Role)) +
          geom_point(size=3, alpha=0.85) +
          geom_vline(xintercept=0.62, linetype=2, colour="grey50") +
          geom_hline(yintercept=2.5,  linetype=2, colour="grey50") +
          scale_colour_manual(values=role_cols) +
          labs(title="ZiPi Plot — Node Ecological Roles",
               x="Participation Coefficient (Pi)",
               y="Within-module Connectivity (Zi)") +
          theme_bw(base_size=10)
        if (requireNamespace("ggrepel", quietly=TRUE)) {
          hub_nodes <- node_df[node_df$Role != "Peripheral", ]
          if (nrow(hub_nodes) > 0)
            p_zipi <- p_zipi +
              ggrepel::geom_text_repel(data=hub_nodes,
                aes(label=taxon), size=2.8, show.legend=FALSE)
        }
        save_pdf(p_zipi, "15b_zipi_plot.pdf", width=8, height=6)
        cat(sprintf("  Network: %d nodes, %d edges | hubs: %d\n",
                    igraph::vcount(g_net), igraph::ecount(g_net),
                    sum(node_df$Role %in% c("Module hub","Network hub"))))
      } else {
        cat("  No edges passed threshold (|r|≥0.6, p<0.05) — try more samples\n")
      }
    } else {
      cat("  Too few taxa or samples for network analysis\n")
    }
  }, error=function(e) cat(sprintf("[WARN] Network section: %s\n", e$message)))
} else {
  cat("  Skipped (install igraph + Hmisc: install.packages(c('igraph','Hmisc')))\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 16 — RDA / CCA (microbiome vs environmental factors)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 16: RDA / CCA ───────────────────────────────────────\n")
if (has_phyloseq && has_vegan && has_meta && has_ggplot2) {
  tryCatch({
    # Find numeric environmental columns (excluding group col and SampleID)
    skip_cols <- c("SampleID", GROUP_COL)
    env_cols  <- setdiff(colnames(meta_df), skip_cols)
    num_cols  <- env_cols[sapply(meta_df[, env_cols, drop=FALSE],
                                 function(x) is.numeric(x) || all(!is.na(suppressWarnings(as.numeric(x)))))]

    if (length(num_cols) == 0) {
      cat("  Skipped — no numeric environmental variables found in metadata\n")
      cat("  (add columns like pH, temperature, depth etc. to metadata TSV)\n")
    } else {
      cat(sprintf("  Environmental variables: %s\n", paste(num_cols, collapse=", ")))

      # Collapse to genus, CLR transform
      ps_rda <- if ("Genus" %in% rank_names(ps))
                  tax_glom(ps, "Genus", NArm=FALSE) else ps
      otu_rda <- as.matrix(otu_table(ps_rda))
      if (taxa_are_rows(ps_rda)) otu_rda <- t(otu_rda)

      # CLR transform
      otu_clr <- t(apply(otu_rda, 1, function(x) {
        xp <- x + 1; log(xp / exp(mean(log(xp))))
      }))

      # Align samples
      common_rda <- intersect(rownames(otu_clr), rownames(meta_df))
      otu_clr    <- otu_clr[common_rda, , drop=FALSE]
      env_df     <- meta_df[common_rda, num_cols, drop=FALSE]
      env_df     <- as.data.frame(lapply(env_df, as.numeric))

      # Remove columns with NA or zero variance
      env_df <- env_df[, sapply(env_df, function(x)
        !all(is.na(x)) && var(x, na.rm=TRUE) > 0), drop=FALSE]

      if (ncol(env_df) == 0 || nrow(env_df) < 5) {
        cat("  Skipped — not enough valid env data after filtering\n")
      } else {
        # RDA (linear, CLR data)
        rda_res <- vegan::rda(otu_clr ~ ., data=env_df)
        rda_sum <- summary(rda_res)

        # Extract scores
        site_sc  <- as.data.frame(scores(rda_res, display="sites",  choices=1:2))
        biplot_sc <- as.data.frame(scores(rda_res, display="bp",    choices=1:2))
        pct       <- round(rda_sum$cont$importance[2, 1:2] * 100, 1)

        site_sc$Sample <- rownames(site_sc)
        if (has_meta)
          site_sc$Group <- meta_df[site_sc$Sample, GROUP_COL]

        biplot_sc$Variable <- rownames(biplot_sc)

        n_grp_rda <- length(unique(site_sc$Group))
        pal_rda   <- make_palette(n_grp_rda)

        p_rda <- ggplot(site_sc, aes(x=RDA1, y=RDA2)) +
          geom_point(aes(colour=Group), size=3, alpha=0.85) +
          scale_colour_manual(values=pal_rda) +
          geom_segment(data=biplot_sc,
                       aes(x=0, y=0, xend=RDA1*0.8, yend=RDA2*0.8),
                       arrow=arrow(length=unit(0.25,"cm")),
                       colour="#374151", linewidth=0.7) +
          ggrepel::geom_text_repel(data=biplot_sc,
                       aes(x=RDA1*0.85, y=RDA2*0.85, label=Variable),
                       size=3.2, colour="#374151") +
          geom_hline(yintercept=0, linetype=2, colour="grey70") +
          geom_vline(xintercept=0, linetype=2, colour="grey70") +
          labs(title=sprintf("RDA — Microbiome vs Environmental Factors\nGrouped by: %s", GROUP_COL),
               x=sprintf("RDA1 (%.1f%%)", pct[1]),
               y=sprintf("RDA2 (%.1f%%)", pct[2]),
               colour=GROUP_COL) +
          theme_bw(base_size=10)
        save_pdf(p_rda, "16_rda_ordination.pdf", width=9, height=7)

        # PERMANOVA against env variables
        perm_rda <- vegan::adonis2(otu_clr ~ ., data=env_df,
                                   permutations=999, method="euclidean")
        write.csv(as.data.frame(perm_rda),
                  file.path(TABLES_DIR, "rda_permanova.csv"))
        cat(sprintf("  RDA: RDA1=%.1f%%, RDA2=%.1f%%\n", pct[1], pct[2]))
      }
    }
  }, error=function(e) cat(sprintf("[WARN] RDA section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + vegan + metadata)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 17 — Beta Distance Heatmap (sample × sample)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 17: Beta Distance Heatmap ──────────────────────────\n")
if (has_phyloseq && has_vegan && requireNamespace("pheatmap", quietly=TRUE)) {
  tryCatch({
    ps_bh  <- transform_sample_counts(ps, function(x) x / sum(x))
    otu_bh <- as.matrix(otu_table(ps_bh))
    if (taxa_are_rows(ps_bh)) otu_bh <- t(otu_bh)

    dist_methods_bh <- c("bray", "jaccard")
    for (dm_bh in dist_methods_bh) {
      tryCatch({
        d_bh  <- as.matrix(vegan::vegdist(otu_bh, method=dm_bh))

        # Annotation bar by group
        ann_col_bh <- NULL
        ann_col_colors <- NULL
        if (has_meta) {
          grp_bh <- meta_df[rownames(d_bh), GROUP_COL]
          grp_bh[is.na(grp_bh)] <- "Unknown"
          ann_col_bh <- data.frame(Group=grp_bh, row.names=rownames(d_bh))
          pal_bh <- make_palette(length(unique(grp_bh)))
          ann_col_colors <- list(Group=setNames(pal_bh, unique(grp_bh)))
        }

        fname_bh <- sprintf("17_beta_heatmap_%s.pdf", dm_bh)
        pdf(file.path(PLOTS_DIR, fname_bh),
            width=max(7, ncol(d_bh)*0.5+2),
            height=max(6, nrow(d_bh)*0.5+2))
        pheatmap::pheatmap(
          d_bh,
          color            = colorRampPalette(c("#eff6ff","#3b82f6","#1e3a8a"))(50),
          clustering_distance_rows = as.dist(d_bh),
          clustering_distance_cols = as.dist(d_bh),
          clustering_method = "average",
          annotation_col    = ann_col_bh,
          annotation_row    = ann_col_bh,
          annotation_colors = ann_col_colors,
          main    = sprintf("Beta Diversity Heatmap — %s distance", dm_bh),
          fontsize = 9,
          border_color = NA
        )
        dev.off()
        cat(sprintf("  ✓ Saved: %s\n", fname_bh))
      }, error=function(e) cat(sprintf("  [WARN] Heatmap %s: %s\n", dm_bh, e$message)))
    }
  }, error=function(e) cat(sprintf("[WARN] Beta heatmap section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + vegan + pheatmap)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 18 — OTU/ASV Distribution Chart
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 18: OTU/ASV Distribution ───────────────────────────\n")
if (has_phyloseq && has_ggplot2) {
  tryCatch({
    # Per-sample read counts and ASV richness
    samp_reads <- sample_sums(ps)
    samp_rich  <- apply(otu_table(ps), ifelse(taxa_are_rows(ps), 2, 1),
                        function(x) sum(x > 0))

    df_dist <- data.frame(
      Sample  = names(samp_reads),
      Reads   = as.integer(samp_reads),
      Richness = as.integer(samp_rich),
      stringsAsFactors = FALSE
    )
    if (has_meta)
      df_dist$Group <- meta_df[df_dist$Sample, GROUP_COL]
    else
      df_dist$Group <- "All"

    df_dist$Sample <- factor(df_dist$Sample,
                             levels=df_dist$Sample[order(df_dist$Reads, decreasing=TRUE)])
    n_grp_od <- length(unique(df_dist$Group))
    pal_od   <- make_palette(n_grp_od)

    # Panel A: read count bar
    p_reads <- ggplot(df_dist, aes(x=Sample, y=Reads, fill=Group)) +
      geom_col(alpha=0.85) +
      scale_fill_manual(values=pal_od) +
      geom_hline(yintercept=mean(df_dist$Reads), linetype=2, colour="grey40") +
      labs(title="Read Count per Sample",
           x="", y="Total Reads", fill=GROUP_COL) +
      theme_bw(base_size=10) +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")

    # Panel B: ASV richness bar
    p_rich <- ggplot(df_dist, aes(x=Sample, y=Richness, fill=Group)) +
      geom_col(alpha=0.85) +
      scale_fill_manual(values=pal_od) +
      geom_hline(yintercept=mean(df_dist$Richness), linetype=2, colour="grey40") +
      labs(title="ASV Richness per Sample",
           x="", y="Number of ASVs", fill=GROUP_COL) +
      theme_bw(base_size=10) +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="none")

    # Panel C: ASV abundance histogram (how many ASVs have each read count)
    asv_sums <- rowSums(otu_table(ps))
    if (!taxa_are_rows(ps)) asv_sums <- colSums(otu_table(ps))
    df_asv_hist <- data.frame(TotalReads=as.integer(asv_sums))

    p_hist <- ggplot(df_asv_hist, aes(x=TotalReads)) +
      geom_histogram(bins=40, fill="#3b82f6", colour="white", alpha=0.85) +
      scale_x_log10(labels=scales::comma) +
      labs(title="ASV Abundance Distribution",
           x="Total Reads per ASV (log10)", y="Count") +
      theme_bw(base_size=10)

    if (requireNamespace("patchwork", quietly=TRUE)) {
      p_od <- patchwork::wrap_plots(p_reads, p_rich, p_hist, ncol=1)
      save_pdf(p_od, "18_otu_distribution.pdf",
               width=max(8, nsamples(ps)*0.5+2), height=14)
    } else {
      save_pdf(p_reads, "18a_read_counts.pdf",
               width=max(8, nsamples(ps)*0.5+2), height=5)
      save_pdf(p_rich,  "18b_asv_richness.pdf",
               width=max(8, nsamples(ps)*0.5+2), height=5)
      save_pdf(p_hist,  "18c_asv_histogram.pdf", width=7, height=5)
    }
    write.csv(df_dist, file.path(TABLES_DIR, "otu_distribution.csv"),
              row.names=FALSE)
  }, error=function(e) cat(sprintf("[WARN] OTU distribution: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + ggplot2)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 19 — STAMP-style Extended Error Bar Plots
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 19: STAMP-style Error Bar Plots ─────────────────────\n")
if (has_phyloseq && has_meta && has_ggplot2) {
  tryCatch({
    grps_st   <- unique(meta_df[[GROUP_COL]])
    pairs_st  <- combn(grps_st, 2, simplify=FALSE)

    for (rank_st in c("Genus", "Family", "Phylum")) {
      if (!(rank_st %in% rank_names(ps))) next
      tryCatch({
        ps_st  <- tax_glom(ps, taxrank=rank_st, NArm=FALSE)
        ps_st  <- transform_sample_counts(ps_st, function(x) x / sum(x) * 100)
        otu_st <- as.matrix(otu_table(ps_st))
        if (taxa_are_rows(ps_st)) otu_st <- t(otu_st)

        tax_st <- sub("^[a-z]__", "",
                      as.character(tax_table(ps_st)[, rank_st]))
        tax_st[is.na(tax_st) | tax_st == ""] <-
          taxa_names(ps_st)[is.na(tax_st) | tax_st == ""]
        colnames(otu_st) <- tax_st

        for (pair_st in pairs_st) {
          g1s <- pair_st[1]; g2s <- pair_st[2]
          s1s <- intersect(rownames(meta_df)[meta_df[[GROUP_COL]] == g1s],
                           rownames(otu_st))
          s2s <- intersect(rownames(meta_df)[meta_df[[GROUP_COL]] == g2s],
                           rownames(otu_st))
          if (length(s1s) < 2 || length(s2s) < 2) next

          # Mean ± SE per group, t-test p-value
          stamp_rows <- do.call(rbind, lapply(colnames(otu_st), function(tx) {
            x1 <- otu_st[s1s, tx]; x2 <- otu_st[s2s, tx]
            p  <- tryCatch(t.test(x1, x2)$p.value, error=function(e) NA)
            data.frame(
              taxon    = tx,
              mean1    = mean(x1), se1 = sd(x1)/sqrt(length(x1)),
              mean2    = mean(x2), se2 = sd(x2)/sqrt(length(x2)),
              p_value  = p,
              stringsAsFactors = FALSE
            )
          }))
          stamp_rows$p_adj <- p.adjust(stamp_rows$p_value, method="BH")
          stamp_rows       <- stamp_rows[order(stamp_rows$p_adj), ]

          # Top 25 by mean abundance (regardless of significance)
          stamp_rows$mean_total <- stamp_rows$mean1 + stamp_rows$mean2
          top_st <- head(stamp_rows[order(stamp_rows$mean_total, decreasing=TRUE), ], 25)
          top_st$taxon <- factor(top_st$taxon, levels=rev(top_st$taxon))
          top_st$sig   <- ifelse(!is.na(top_st$p_adj) & top_st$p_adj < 0.05, "*", "")

          pal2_st <- make_palette(2)

          # Build long format for error bars
          df_long_st <- rbind(
            data.frame(taxon=top_st$taxon, mean=top_st$mean1,
                       se=top_st$se1, group=g1s, sig=top_st$sig,
                       stringsAsFactors=FALSE),
            data.frame(taxon=top_st$taxon, mean=top_st$mean2,
                       se=top_st$se2, group=g2s, sig=top_st$sig,
                       stringsAsFactors=FALSE)
          )

          p_st <- ggplot(df_long_st,
                         aes(x=taxon, y=mean, colour=group)) +
            geom_point(position=position_dodge(0.5), size=3) +
            geom_errorbar(aes(ymin=pmax(mean-se*1.96, 0), ymax=mean+se*1.96),
                          position=position_dodge(0.5), width=0.3, linewidth=0.7) +
            geom_text(data=subset(df_long_st, group==g1s & sig=="*"),
                      aes(label=sig, y=mean+se*1.96+0.3),
                      colour="black", size=5, show.legend=FALSE) +
            coord_flip() +
            scale_colour_manual(values=setNames(pal2_st, c(g1s, g2s))) +
            labs(title=sprintf("STAMP — %s: %s vs %s\n(top 25 by abundance, * = q<0.05)",
                               rank_st, g1s, g2s),
                 x="", y="Mean Relative Abundance (%)", colour="Group") +
            theme_bw(base_size=10) +
            theme(legend.position="top")

          fname_st <- sprintf("19_stamp_%s_%s_vs_%s.pdf",
                              tolower(rank_st), g1s, g2s)
          save_pdf(p_st, fname_st,
                   width=9, height=max(6, nrow(top_st)*0.38+2))
          cat(sprintf("  %s %s vs %s: done (%d sig taxa)\n",
                      rank_st, g1s, g2s,
                      sum(!is.na(top_st$p_adj) & top_st$p_adj < 0.05)))
        }
      }, error=function(e) cat(sprintf("  [WARN] STAMP %s: %s\n", rank_st, e$message)))
    }
  }, error=function(e) cat(sprintf("[WARN] STAMP section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + metadata + ggplot2)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 20 — UniFrac PCoA + NMDS (Weighted + Unweighted)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 20: UniFrac PCoA + NMDS ─────────────────────────────\n")
if (has_phyloseq && has_ggplot2) {
  tryCatch({
    ps_uf    <- ps
    tree_ok  <- has_tree

    # Auto-build NJ tree if not present and rep-seqs FASTA available
    if (!tree_ok && file.exists(seq_fasta) &&
        requireNamespace("DECIPHER",    quietly=TRUE) &&
        requireNamespace("phangorn",    quietly=TRUE) &&
        requireNamespace("Biostrings",  quietly=TRUE)) {
      tryCatch({
        suppressPackageStartupMessages({ library(DECIPHER); library(phangorn) })
        cat("  Building NJ tree from rep-seqs (may take a few minutes)...\n")
        seqs_uf <- Biostrings::readDNAStringSet(seq_fasta)
        if (length(seqs_uf) > 300) {
          top300  <- names(sort(taxa_sums(ps_uf), decreasing=TRUE))[1:300]
          seqs_uf <- seqs_uf[intersect(names(seqs_uf), top300)]
        }
        cat(sprintf("  Aligning %d sequences...\n", length(seqs_uf)))
        aln_uf      <- DECIPHER::AlignSeqs(seqs_uf, verbose=FALSE)
        phang_uf    <- phangorn::phyDat(as(aln_uf, "matrix"), type="DNA")
        dm_uf       <- phangorn::dist.ml(phang_uf)
        tree_uf_obj <- phangorn::midpoint(phangorn::NJ(dm_uf))
        common_uf   <- intersect(taxa_names(ps_uf), tree_uf_obj$tip.label)
        ps_uf       <- prune_taxa(common_uf, ps_uf)
        phy_tree(ps_uf) <- tree_uf_obj
        tree_ok     <- TRUE
        cat(sprintf("  ✓ NJ tree: %d taxa\n", length(common_uf)))
      }, error=function(e) cat(sprintf("  [WARN] Tree build: %s\n", e$message)))
    }

    if (!tree_ok) {
      cat("  Skipped — no tree (provide tree/tree.nwk or install DECIPHER+phangorn)\n")
    } else {
      n_grp_uf <- if (has_meta) length(unique(meta_df[[GROUP_COL]])) else 1
      pal_uf   <- make_palette(n_grp_uf)

      uf_list <- list(
        list(d="unifrac",  m="PCoA", lbl="Unweighted UniFrac PCoA",  f="20a_unifrac_unwt_pcoa.pdf"),
        list(d="wunifrac", m="PCoA", lbl="Weighted UniFrac PCoA",    f="20b_unifrac_wt_pcoa.pdf"),
        list(d="unifrac",  m="NMDS", lbl="Unweighted UniFrac NMDS",  f="20c_unifrac_unwt_nmds.pdf"),
        list(d="wunifrac", m="NMDS", lbl="Weighted UniFrac NMDS",    f="20d_unifrac_wt_nmds.pdf")
      )
      for (uf in uf_list) {
        tryCatch({
          ord_uf <- phyloseq::ordinate(ps_uf, method=uf$m, distance=uf$d)
          p_uf <- if (has_meta) {
            base_uf <- phyloseq::plot_ordination(ps_uf, ord_uf, color=GROUP_COL) +
              geom_point(size=3.5, alpha=0.85) +
              scale_colour_manual(values=pal_uf) +
              labs(title=uf$lbl,
                   subtitle=sprintf("n=%d samples | grouped by: %s",
                                    nsamples(ps_uf), GROUP_COL)) +
              theme_bw(base_size=11)
            # Add ellipse only if enough samples per group
            grp_counts_uf <- table(meta_df[[GROUP_COL]])
            if (all(grp_counts_uf >= 3))
              base_uf <- base_uf +
                stat_ellipse(aes_string(group=GROUP_COL), level=0.95,
                             linetype=2, linewidth=0.6)
            base_uf
          } else {
            phyloseq::plot_ordination(ps_uf, ord_uf) +
              geom_point(size=3.5, colour="#3b82f6", alpha=0.85) +
              labs(title=uf$lbl) + theme_bw(base_size=11)
          }
          save_pdf(p_uf, uf$f, width=8, height=6)
        }, error=function(e) cat(sprintf("  [WARN] %s: %s\n", uf$lbl, e$message)))
      }
    }
  }, error=function(e) cat(sprintf("[WARN] UniFrac section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + ggplot2)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 21 — Phylogenetic Tree Visualization (ggtree)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 21: Phylogenetic Tree (ggtree) ──────────────────────\n")
if (has_phyloseq && requireNamespace("ggtree", quietly=TRUE)) {
  tryCatch({
    suppressPackageStartupMessages(library(ggtree))

    tree_viz <- tryCatch(phy_tree(ps), error=function(e) NULL)
    if (is.null(tree_viz)) {
      cat("  Skipped — no tree in phyloseq object\n")
    } else {
      # Prune to top 100 taxa by abundance for readability
      ps_tv <- ps
      if (ntaxa(ps) > 100) {
        top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
        ps_tv  <- prune_taxa(top100, ps)
        cat("  Pruned to top 100 taxa for readability\n")
      }
      tree_tv <- phy_tree(ps_tv)

      # Build tip annotation data frame
      tip_df_tv <- data.frame(label=tree_tv$tip.label, stringsAsFactors=FALSE)
      tip_df_tv$Abundance <- log10(taxa_sums(ps_tv)[tip_df_tv$label] + 1)

      if ("Phylum" %in% rank_names(ps_tv)) {
        tip_df_tv$Phylum <- sub("^[a-z]__", "",
          as.character(tax_table(ps_tv)[tip_df_tv$label, "Phylum"]))
      }
      if ("Genus" %in% rank_names(ps_tv)) {
        tip_df_tv$Genus <- sub("^[a-z]__", "",
          as.character(tax_table(ps_tv)[tip_df_tv$label, "Genus"]))
        tip_df_tv$Genus[is.na(tip_df_tv$Genus) | tip_df_tv$Genus==""] <- ""
      }

      p_tree <- ggtree(tree_tv, layout="circular", linewidth=0.3,
                        colour="grey60") %<+% tip_df_tv

      if ("Phylum" %in% colnames(tip_df_tv)) {
        n_phy_tv <- length(unique(na.omit(tip_df_tv$Phylum)))
        pal_tv   <- make_palette(min(n_phy_tv, 12))
        p_tree <- p_tree +
          geom_tippoint(aes(colour=Phylum, size=Abundance), alpha=0.85) +
          scale_colour_manual(values=pal_tv, na.value="grey70",
                              name="Phylum") +
          scale_size_continuous(range=c(1.5, 5), name="log10(reads)") +
          labs(title=sprintf("Phylogenetic Tree — top %d ASVs", ntaxa(ps_tv)),
               subtitle="Tip colour = Phylum | size = abundance") +
          theme(legend.text=element_text(size=7),
                plot.title=element_text(hjust=0.5))
      } else {
        p_tree <- p_tree +
          geom_tippoint(aes(size=Abundance), colour="#3b82f6", alpha=0.8) +
          scale_size_continuous(range=c(1.5, 5), name="log10(reads)") +
          labs(title="Phylogenetic Tree")
      }
      save_pdf(p_tree, "21_phylogenetic_tree.pdf", width=11, height=11)
    }
  }, error=function(e) cat(sprintf("[WARN] ggtree section: %s\n", e$message)))
} else if (!requireNamespace("ggtree", quietly=TRUE)) {
  cat("  Skipped (install: BiocManager::install('ggtree'))\n")
} else {
  cat("  Skipped (requires phyloseq)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 22 — Krona Charts (per-sample taxonomy HTML)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 22: Krona Charts ────────────────────────────────────\n")
if (has_phyloseq) {
  tryCatch({
    krona_dir <- file.path(PLOTS_DIR, "krona")
    dir.create(krona_dir, showWarnings=FALSE)

    rank_order_kr <- intersect(
      c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
      rank_names(ps))

    krona_files <- c()
    for (samp_kr in sample_names(ps)) {
      tryCatch({
        counts_kr <- if (taxa_are_rows(ps)) otu_table(ps)[, samp_kr] else
                     otu_table(ps)[samp_kr, ]
        counts_kr <- counts_kr[counts_kr > 0]
        if (length(counts_kr) == 0) next

        tax_kr  <- tax_table(ps)[names(counts_kr), rank_order_kr, drop=FALSE]
        lines_kr <- sapply(seq_along(counts_kr), function(i) {
          path_kr <- sub("^[a-z]__", "", as.character(tax_kr[i, ]))
          path_kr <- path_kr[!is.na(path_kr) & nchar(trimws(path_kr)) > 0]
          paste(c(as.integer(counts_kr[i]), path_kr), collapse="\t")
        })
        fname_kr <- file.path(krona_dir, paste0(samp_kr, ".krona.txt"))
        writeLines(lines_kr, fname_kr)
        krona_files <- c(krona_files, fname_kr)
      }, error=function(e) NULL)
    }
    cat(sprintf("  Written %d Krona text files\n", length(krona_files)))

    # Call ktImportText if KronaTools is installed
    kt_bin <- Sys.which("ktImportText")
    if (nchar(kt_bin) > 0 && length(krona_files) > 0) {
      html_kr <- file.path(PLOTS_DIR, "22_krona_all_samples.html")
      cmd_kr  <- paste(c(kt_bin, "-o", shQuote(html_kr),
                          shQuote(krona_files)), collapse=" ")
      ret_kr  <- system(cmd_kr, ignore.stderr=TRUE)
      if (ret_kr == 0)
        cat("  ✓ Saved: 22_krona_all_samples.html\n")
      else
        cat("  [WARN] ktImportText failed — text files saved for manual use\n")
    } else {
      cat("  KronaTools not found — to generate HTML:\n")
      cat("    sudo apt-get install -y krona\n")
      cat(sprintf("    ktImportText -o krona.html %s/*.krona.txt\n", krona_dir))
    }
  }, error=function(e) cat(sprintf("[WARN] Krona section: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 23 — LEfSe Cladogram (taxonomy tree colored by enrichment)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Section 23: LEfSe Cladogram ─────────────────────────────────\n")
if (has_phyloseq && has_ggplot2 && requireNamespace("igraph", quietly=TRUE)) {
  tryCatch({
    lefse_csv <- file.path(TABLES_DIR, "lefse_results.csv")
    if (!file.exists(lefse_csv)) {
      cat("  Skipped — lefse_results.csv not found (run LEfSe first)\n")
    } else {
      lefse_df_cl <- read.csv(lefse_csv, stringsAsFactors=FALSE)
      if (nrow(lefse_df_cl) == 0) {
        cat("  Skipped — no significant LEfSe biomarkers\n")
      } else {
        suppressPackageStartupMessages(library(igraph))

        rank_clad <- intersect(
          c("Phylum","Class","Order","Family","Genus"), rank_names(ps))
        tax_clad  <- as.data.frame(tax_table(ps)[, rank_clad, drop=FALSE],
                                    stringsAsFactors=FALSE)
        tax_clad  <- as.data.frame(lapply(tax_clad, function(x)
          sub("^[a-z]__", "", trimws(x))), stringsAsFactors=FALSE)

        # Build edge list: Root → Phylum → Class → ...
        edges_cl <- data.frame(from="Root",
                                to=unique(na.omit(tax_clad[,1][nchar(tax_clad[,1])>0])),
                                stringsAsFactors=FALSE)
        for (r in seq_len(length(rank_clad)-1)) {
          pairs_cl <- unique(tax_clad[, c(r, r+1), drop=FALSE])
          pairs_cl <- pairs_cl[!is.na(pairs_cl[,1]) & !is.na(pairs_cl[,2]) &
                                  nchar(pairs_cl[,1])>0 & nchar(pairs_cl[,2])>0, ]
          if (nrow(pairs_cl) > 0)
            edges_cl <- rbind(edges_cl,
              data.frame(from=pairs_cl[,1], to=pairs_cl[,2],
                         stringsAsFactors=FALSE))
        }
        edges_cl <- unique(edges_cl)
        edges_cl <- edges_cl[edges_cl$from != edges_cl$to, ]

        g_cl <- igraph::graph_from_data_frame(edges_cl, directed=TRUE)
        node_names_cl <- igraph::V(g_cl)$name

        # LEfSe enrichment lookup
        grp_levels_cl <- unique(lefse_df_cl$Direction)
        n_grp_cl      <- length(grp_levels_cl)
        pal_cl_grp    <- make_palette(n_grp_cl)
        col_lookup_cl <- c(setNames(pal_cl_grp, grp_levels_cl), "none"="grey85")

        lef_dir_cl   <- setNames(lefse_df_cl$Direction, lefse_df_cl$Label)
        lef_score_cl <- setNames(abs(lefse_df_cl$scores), lefse_df_cl$Label)

        node_col_cl  <- sapply(node_names_cl, function(n)
          if (n %in% names(lef_dir_cl)) col_lookup_cl[lef_dir_cl[n]] else "grey85")
        node_size_cl <- sapply(node_names_cl, function(n)
          if (n %in% names(lef_score_cl)) 4 + lef_score_cl[n]*1.5 else 3)

        if (requireNamespace("ggraph", quietly=TRUE)) {
          suppressPackageStartupMessages(library(ggraph))
          p_clad <- ggraph(g_cl, layout="dendrogram", circular=TRUE) +
            geom_edge_diagonal(colour="grey75", alpha=0.5, linewidth=0.35) +
            geom_node_point(aes(size=node_size_cl), colour=node_col_cl,
                            alpha=0.9) +
            geom_node_text(
              aes(label=ifelse(name %in% names(lef_dir_cl), name, "")),
              repel=TRUE, size=2.2, colour=node_col_cl, max.overlaps=20) +
            scale_size_continuous(range=c(1.5, 8), guide="none") +
            labs(title=sprintf("LEfSe Cladogram — %d significant biomarkers",
                               nrow(lefse_df_cl)),
                 subtitle=paste(grp_levels_cl, collapse=" vs ")) +
            theme_void(base_size=9) +
            theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
                  plot.subtitle=element_text(hjust=0.5))
          save_pdf(p_clad, "23_lefse_cladogram.pdf", width=12, height=12)
        } else {
          pdf(file.path(PLOTS_DIR, "23_lefse_cladogram.pdf"), width=12, height=12)
          igraph::plot.igraph(g_cl,
            layout           = igraph::layout_as_tree(g_cl, circular=TRUE),
            vertex.color     = node_col_cl,
            vertex.size      = pmin(node_size_cl, 15),
            vertex.label     = ifelse(node_names_cl %in% names(lef_dir_cl),
                                      node_names_cl, ""),
            vertex.label.cex = 0.45,
            edge.arrow.size  = 0.2,
            main = sprintf("LEfSe Cladogram (%d biomarkers)", nrow(lefse_df_cl)))
          legend("bottomright", legend=c(grp_levels_cl, "Not significant"),
                 fill=c(pal_cl_grp, "grey85"), bty="n", cex=0.8)
          dev.off()
          cat("  ✓ Saved: 23_lefse_cladogram.pdf\n")
        }
        cat(sprintf("  LEfSe cladogram: %d biomarkers, %d nodes\n",
                    nrow(lefse_df_cl), igraph::vcount(g_cl)))
      }
    }
  }, error=function(e) cat(sprintf("[WARN] LEfSe cladogram: %s\n", e$message)))
} else {
  cat("  Skipped (requires phyloseq + ggplot2 + igraph)\n")
}

# ─── Summary ──────────────────────────────────────────────────────────────────
plots_made <- list.files(PLOTS_DIR, pattern="\\.pdf$")
tables_made <- list.files(TABLES_DIR, pattern="\\.csv$")

cat("\n============================================================\n")
cat(sprintf("  R Visualization complete!\n"))
cat(sprintf("  PDFs created : %d  →  %s\n", length(plots_made),  PLOTS_DIR))
cat(sprintf("  CSVs created : %d  →  %s\n", length(tables_made), TABLES_DIR))
cat("============================================================\n\n")


# ── DADA2 Extra Viz (PCoA, NMDS, Jaccard, Rarefaction, ASV Lengths) ──────────
# Source external helper → defines dada2_extra_viz() in global env.
# Falls back to inline definition if file not found or source fails.
tryCatch({
  this_file  <- grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
  sdir       <- if (length(this_file) > 0)
                  dirname(normalizePath(sub("--file=", "", this_file[1])))
                else getwd()
  xscript    <- file.path(sdir, "dada2_extra_viz.R")
  if (file.exists(xscript)) source(xscript)   # local=FALSE (default) → global env
}, error=function(e) invisible(NULL))

# Inline definition used only if source above did not define the function
if (!exists("dada2_extra_viz", mode="function")) {
  dada2_extra_viz <- function(out_dir) {
    bray_file <- file.path(out_dir, "bray_curtis_distance_matrix.csv")
    asv_file  <- file.path(out_dir, "asv_table.csv")
    tax_file  <- file.path(out_dir, "taxonomy_table.csv")
    have_any  <- file.exists(bray_file) || file.exists(asv_file) || file.exists(tax_file)
    if (!have_any) return(invisible(NULL))
    cat("\n── DADA2 Extra Viz (inline) ────────────────────────────────\n")
    cat("  Directory:", out_dir, "\n")

    have_vegan <- suppressPackageStartupMessages(
      tryCatch({ library(vegan); TRUE },
               error=function(e) { cat("  vegan not available\n"); FALSE }))

    load_asv_mat <- function() {
      raw <- read.csv(asv_file, check.names=FALSE)
      sc  <- which(colnames(raw) == "sequence")
      if (length(sc) > 0) raw <- raw[, -sc, drop=FALSE]
      mat <- t(as.matrix(raw)); storage.mode(mat) <- "numeric"; mat
    }

    rarefy_exact <- function(counts, depth) {
      counts <- as.integer(counts[counts > 0]); total <- sum(counts)
      if (depth >= total) return(sum(counts > 0))
      log_CN  <- lgamma(total+1) - lgamma(depth+1) - lgamma(total-depth+1)
      expected <- 0
      for (ni in counts) {
        if (ni > total - depth) { expected <- expected + 1 } else {
          log_CNni <- lgamma(total-ni+1) - lgamma(depth+1) - lgamma(total-ni-depth+1)
          expected <- expected + (1 - exp(log_CNni - log_CN))
        }
      }
      round(expected, 1)
    }

    # PCoA
    if (file.exists(bray_file)) tryCatch({
      pca_out <- file.path(out_dir, "pca_scores.csv")
      if (!file.exists(pca_out)) {
        dm <- as.matrix(read.csv(bray_file, row.names=1, check.names=FALSE))
        k  <- min(nrow(dm)-1, 10); pc <- cmdscale(dm, k=k, eig=TRUE)
        ep <- pc$eig; ep[ep<0] <- 0
        vp <- if(sum(ep)>0) pc$eig/sum(ep)*100 else rep(0,length(ep)); vp[vp<0]<-0
        write.csv(data.frame(Sample=rownames(dm), PC1=pc$points[,1], PC2=pc$points[,2],
          PC3=if(k>=3)pc$points[,3]else 0, PC1_var=round(vp[1],1), PC2_var=round(vp[2],1),
          PC3_var=if(length(vp)>=3)round(vp[3],1)else 0), pca_out, row.names=FALSE)
        cat("  ✓ pca_scores.csv\n")
        ns <- min(10,sum(ep>0))
        write.csv(data.frame(PC=seq_len(ns),Variance=round(vp[seq_len(ns)],2)),
                  file.path(out_dir,"pca_scree.csv"), row.names=FALSE)
        cat("  ✓ pca_scree.csv\n")
      }
    }, error=function(e) cat("  [skip] PCoA:", e$message, "\n"))

    # NMDS Bray
    if (have_vegan && file.exists(bray_file)) tryCatch({
      nb_out <- file.path(out_dir, "nmds_bray.csv")
      if (!file.exists(nb_out)) {
        dm <- as.matrix(read.csv(bray_file, row.names=1, check.names=FALSE))
        set.seed(42); nm <- metaMDS(dm, distance="bray", k=2, try=20, trymax=50, trace=FALSE)
        write.csv(data.frame(Sample=rownames(dm), NMDS1=nm$points[,1],
                             NMDS2=nm$points[,2], Stress=nm$stress), nb_out, row.names=FALSE)
        cat("  ✓ nmds_bray.csv\n")
      }
    }, error=function(e) cat("  [skip] NMDS Bray:", e$message, "\n"))

    # Jaccard + NMDS Jaccard
    asv_cached <- NULL
    if (have_vegan && file.exists(asv_file)) tryCatch({
      asv_cached <- load_asv_mat()
      jac_out <- file.path(out_dir, "jaccard_heatmap.csv")
      if (!file.exists(jac_out)) {
        jm <- as.matrix(vegdist(asv_cached, method="jaccard", binary=FALSE))
        write.csv(as.data.frame(jm), jac_out)
        cat("  ✓ jaccard_heatmap.csv\n")
        set.seed(42); nmj <- metaMDS(jm, distance="jaccard", k=2, try=20, trymax=50, trace=FALSE)
        write.csv(data.frame(Sample=rownames(jm), NMDS1=nmj$points[,1],
                             NMDS2=nmj$points[,2], Stress=nmj$stress),
                  file.path(out_dir,"nmds_jaccard.csv"), row.names=FALSE)
        cat("  ✓ nmds_jaccard.csv\n")
      }
    }, error=function(e) cat("  [skip] Jaccard:", e$message, "\n"))

    # Rarefaction (no vegan needed — uses base-R rarefy_exact)
    if (file.exists(asv_file)) tryCatch({
      rar_out <- file.path(out_dir, "rarefaction.csv")
      if (!file.exists(rar_out)) {
        mat <- if (!is.null(asv_cached)) asv_cached else load_asv_mat()
        rs  <- rowSums(mat); minr <- min(rs)
        if (minr < 10) stop("min reads too low")
        ns2 <- min(25, minr)
        dep <- unique(sort(c(round(seq(max(100,minr/ns2), minr, length.out=ns2)),
                             as.integer(minr))))
        dep <- dep[dep>=1 & dep<=minr]
        wide <- data.frame(Depth=dep)
        for (i in seq_len(nrow(mat))) {
          sn <- rownames(mat)[i]
          if (have_vegan) {
            wide[[sn]] <- round(as.numeric(suppressWarnings(
              rarefy(mat[i,,drop=FALSE], sample=dep))), 1)
          } else {
            wide[[sn]] <- sapply(dep, function(d) rarefy_exact(mat[i,], d))
          }
        }
        write.csv(wide, rar_out, row.names=FALSE)
        cat("  ✓ rarefaction.csv (", nrow(mat), "samples,", length(dep),
            "depths, min =", minr, "reads )\n")
      }
    }, error=function(e) cat("  [skip] Rarefaction:", e$message, "\n"))

    # ASV lengths
    if (file.exists(tax_file)) tryCatch({
      len_out <- file.path(out_dir, "asv_lengths.csv")
      if (!file.exists(len_out)) {
        seqs <- as.character(read.csv(tax_file, check.names=FALSE)[[1]])
        seqs <- gsub('^"|"$','',seqs); seqs <- seqs[!is.na(seqs) & nchar(seqs)>10]
        lt   <- as.data.frame(table(Length=nchar(seqs)), stringsAsFactors=FALSE)
        colnames(lt) <- c("Length","Count")
        lt$Length <- as.integer(lt$Length); lt$Count <- as.integer(lt$Count)
        write.csv(lt[order(lt$Length),], len_out, row.names=FALSE)
        cat("  ✓ asv_lengths.csv\n")
      }
    }, error=function(e) cat("  [skip] ASV lengths:", e$message, "\n"))

    cat("── DADA2 Extra Viz done ─────────────────────────────────────\n\n")
    invisible(NULL)
  }
}

tryCatch(dada2_extra_viz(OUTPUT_DIR),
         error=function(e) cat("  [skip] dada2_extra_viz:", e$message, "\n"))
