#!/usr/bin/env Rscript
# =============================================================================
#  NextGen-Amplicon — BugBase Phenotype Prediction (16S only)
#  Predicts: aerobic/anaerobic, Gram+/-, biofilm, pathogenicity, etc.
#
#  Usage:
#    Rscript bugbase_run.R \
#      --otu_table  /path/to/job/exported/feature-table/feature-table.tsv \
#      --taxonomy   /path/to/job/exported/taxonomy/taxonomy.tsv \
#      --output_dir /path/to/job/bugbase \
#      --group_col  treatment \
#      --metadata   /path/to/metadata.tsv
#
#  Note: BugBase uses GreenGenes 13.5 taxonomy. Results are approximate
#        when using SILVA-based taxonomy (default in our pipeline).
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--otu_table",  type="character", help="Feature table TSV"),
  make_option("--taxonomy",   type="character", help="Taxonomy TSV"),
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--metadata",   type="character", default=NULL),
  make_option("--group_col",  type="character", default="treatment")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$otu_table) || is.null(opt$output_dir)) {
  cat("ERROR: --otu_table and --output_dir are required\n"); quit(status=1)
}

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)
cat("=== BugBase Phenotype Prediction ===\n\n")

# ── Install BugBase if needed ─────────────────────────────────────────────────
if (!requireNamespace("bugbase", quietly=TRUE)) {
  cat("Installing BugBase from GitHub...\n")
  if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cloud.r-project.org", quiet=TRUE)
  tryCatch(
    remotes::install_github("danknights/bugbase", quiet=TRUE),
    error=function(e) {
      cat(sprintf("[ERROR] BugBase install failed: %s\n", e$message))
      cat("  Try manually: remotes::install_github('danknights/bugbase')\n")
      quit(status=1)
    }
  )
}

suppressPackageStartupMessages(library(bugbase))

# ── Load OTU table ────────────────────────────────────────────────────────────
cat("Loading OTU table...\n")
otu_df <- tryCatch({
  df <- read.table(opt$otu_table, header=TRUE, sep="\t",
                   comment.char="#", row.names=1, check.names=FALSE)
  cat(sprintf("  %d taxa × %d samples\n", nrow(df), ncol(df)))
  df
}, error=function(e) {
  cat(sprintf("[ERROR] OTU table: %s\n", e$message)); quit(status=1)
})

# ── Load metadata ─────────────────────────────────────────────────────────────
meta_df <- NULL
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  meta_df <- read.table(opt$metadata, header=TRUE, sep="\t",
                         comment.char="#", quote="", stringsAsFactors=FALSE)
  colnames(meta_df)[1] <- "SampleID"
  rownames(meta_df)    <- meta_df$SampleID
  cat(sprintf("  Metadata: %d samples\n", nrow(meta_df)))
}

# ── Run BugBase ───────────────────────────────────────────────────────────────
cat("Running BugBase prediction...\n")
tryCatch({
  # BugBase expects OTU table with samples as columns, taxa as rows
  otu_mat <- as.matrix(otu_df)

  # Run prediction
  pred <- bugbase::predict.phenotypes(otu_mat)

  # Save full results
  write.csv(pred, file.path(opt$output_dir, "bugbase_predictions.csv"))
  cat(sprintf("  Predicted phenotypes: %s\n",
              paste(colnames(pred), collapse=", ")))

  # ── Plots per phenotype ─────────────────────────────────────────────────
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    cat("  [WARN] ggplot2 not available — skipping plots\n")
  } else {
    suppressPackageStartupMessages(library(ggplot2))

    for (phenotype in colnames(pred)) {
      tryCatch({
        df_plot <- data.frame(
          Sample     = rownames(pred),
          Proportion = pred[, phenotype],
          stringsAsFactors = FALSE
        )
        if (!is.null(meta_df)) {
          df_plot$Group <- meta_df[df_plot$Sample, opt$group_col]
          df_plot$Group[is.na(df_plot$Group)] <- "Unknown"
        } else {
          df_plot$Group <- "All"
        }

        n_grp_bb <- length(unique(df_plot$Group))
        pal_bb   <- if (n_grp_bb <= 8 && requireNamespace("RColorBrewer", quietly=TRUE))
                      RColorBrewer::brewer.pal(max(3, n_grp_bb), "Set2")[1:n_grp_bb]
                    else scales::hue_pal()(n_grp_bb)

        p_bb <- ggplot(df_plot, aes(x=Group, y=Proportion, fill=Group)) +
          geom_boxplot(alpha=0.7, outlier.shape=NA) +
          geom_jitter(width=0.15, size=2, alpha=0.7) +
          scale_fill_manual(values=pal_bb) +
          scale_y_continuous(labels=scales::percent) +
          labs(title=sprintf("BugBase — %s", phenotype),
               x=opt$group_col, y="Predicted Proportion") +
          theme_bw(base_size=11) +
          theme(legend.position="none")

        fname_bb <- file.path(opt$output_dir,
                               sprintf("bugbase_%s.pdf",
                                       gsub("[^a-zA-Z0-9]", "_", phenotype)))
        ggsave(fname_bb, plot=p_bb, width=6, height=5, device="pdf")
        cat(sprintf("  ✓ %s\n", basename(fname_bb)))
      }, error=function(e) cat(sprintf("  [WARN] %s: %s\n", phenotype, e$message)))
    }
  }
  cat("\n✓ BugBase complete\n")
}, error=function(e) {
  cat(sprintf("[ERROR] BugBase prediction: %s\n", e$message))
  cat("  Note: BugBase uses GreenGenes taxonomy — results may be inaccurate with SILVA\n")
})
