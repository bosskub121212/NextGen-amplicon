#!/usr/bin/env Rscript
# =============================================================
#  NextGen-Amplicon -- PacBio CCS Long-Read 16S Pipeline
#  Supports: V1-V9 full-length 16S (~1500 bp), V1-V3, V3-V5, etc.
#
#  Key DADA2 settings for PacBio CCS:
#    - errorEstimationFunction = PacBioErrfun  (CCS-specific)
#    - BAND_SIZE = 32  (wider band for long reads)
#    - Single-end denoising (CCS reads are already full-length)
#    - No merging step required
#    - Chimera removal with bimera consensus
# =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(Biostrings)
  library(ggplot2)
  library(jsonlite)
})

# ------------------------------------------------------------------
# CLI options
# ------------------------------------------------------------------
option_list <- list(
  make_option("--input_dir",    type="character", default=NULL,
              help="Directory containing CCS FASTQ files (single-end)"),
  make_option("--output_dir",   type="character", default=NULL,
              help="Directory for all output files"),
  make_option("--primer_f",     type="character", default="AGRGTTYGATYMTGGCTCAG",
              help="Forward primer (default: 27F for V1-V9)"),
  make_option("--primer_r",     type="character", default="RGYTACCTTGTTACGACTT",
              help="Reverse primer (default: 1492R for V1-V9)"),
  make_option("--min_length",   type="integer",   default=1000,
              help="Minimum read length to keep (bp)"),
  make_option("--max_length",   type="integer",   default=1600,
              help="Maximum read length to keep (bp)"),
  make_option("--maxEE",        type="numeric",   default=3,
              help="Maximum expected errors (CCS can tolerate 2-5)"),
  make_option("--minQ",         type="integer",   default=3,
              help="Minimum Phred quality score per base"),
  make_option("--threads",      type="integer",   default=4,
              help="Number of CPU threads"),
  make_option("--dbPath",       type="character", default="",
              help="Path to SILVA or other 16S taxonomy database"),
  make_option("--db_paths_json",type="character", default="",
              help="Path to db_paths.json for auto-detecting databases"),
  make_option("--minBoot",      type="integer",   default=50,
              help="Minimum bootstrap confidence for taxonomy"),
  make_option("--topN",         type="integer",   default=10,
              help="Top N taxa for bar charts"),
  make_option("--job_name",     type="character", default="PacBio16S",
              help="Job name"),
  make_option("--metadata",     type="character", default="",
              help="Optional: path to metadata CSV/TSV (sample_id first column)"),
  make_option("--pool",         type="character", default="pseudo",
              help="Pooling: FALSE, TRUE, or pseudo"),
  make_option("--region",       type="character", default="V1-V9",
              help="16S region amplified (for record-keeping; e.g., V1-V9, V1-V3)"),
  make_option("--checkpoint_file", type="character", default=NULL,
              help="Write checkpoint progress here for backend polling")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input_dir))  stop("--input_dir is required")
if (is.null(opt$output_dir)) stop("--output_dir is required")

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)
# Note: do NOT use sink() here — backend reads stdout via proc.stdout.
# Writing to a log file via sink() makes Python receive nothing.

cat("=== NextGen-Amplicon PacBio CCS Pipeline ===\n")
cat("Job     :", opt$job_name, "\n")
cat("Region  :", opt$region, "\n")
cat("Input   :", opt$input_dir, "\n")
cat("Output  :", opt$output_dir, "\n")
cat("Threads :", opt$threads, "\n")
cat("Length  :", opt$min_length, "-", opt$max_length, "bp\n")
cat("Started :", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------
# Helper: checkpoint
# ------------------------------------------------------------------
write_checkpoint <- function(step, pct, msg="") {
  if (!is.null(opt$checkpoint_file)) {
    cp <- list(step=step, pct=pct, msg=msg, ts=format(Sys.time()))
    write(toJSON(cp, auto_unbox=TRUE), opt$checkpoint_file)
  }
  # PROGRESS: format is parsed by the backend progress poller
  label <- if (nchar(msg) > 0) paste0(step, ": ", msg) else step
  cat(sprintf("PROGRESS:%.0f|%s\n", pct, label))
  flush.console()
}

write_checkpoint("init", 2, "PacBio CCS pipeline starting")

# ------------------------------------------------------------------
# Step 1: Detect FASTQ files (single-end CCS)
# ------------------------------------------------------------------
write_checkpoint("detect_files", 4, "Detecting CCS FASTQ files")

fq_all <- list.files(opt$input_dir,
                     pattern="\\.(fastq|fastq\\.gz|fq|fq\\.gz)$",
                     full.names=TRUE, recursive=FALSE)

# PacBio CCS files are single-end — filter out any R2 files
fq_all <- fq_all[!grepl("_R2|_2\\.f(q|astq)", fq_all)]

if (length(fq_all) == 0) stop("No CCS FASTQ files found in: ", opt$input_dir)

sample_names <- sub("_R1.*|_1\\.f(q|astq).*|\\.f(q|astq)(\\.gz)?$", "", basename(fq_all))
sample_names <- make.unique(sample_names, sep="_")

cat("Samples detected:", length(sample_names), "\n")
cat(paste(" ", sample_names, collapse="\n"), "\n\n")

# ------------------------------------------------------------------
# Step 2: Quality profiles
# ------------------------------------------------------------------
write_checkpoint("quality_plots", 6, "Plotting quality profiles")

qc_dir <- file.path(opt$output_dir, "QC_plots")
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

tryCatch({
  n_show <- min(length(fq_all), 4)
  pdf(file.path(qc_dir, "quality_profile_CCS.pdf"), width=12, height=6)
  print(plotQualityProfile(fq_all[1:n_show]))
  dev.off()
}, error=function(e) cat("Quality plot warning:", e$message, "\n"))

# ------------------------------------------------------------------
# Step 3: Primer trimming (cutadapt)
# ------------------------------------------------------------------
write_checkpoint("primer_trim", 10, "Trimming primers")

trim_dir <- file.path(opt$output_dir, "Trimmed")
dir.create(trim_dir, recursive=TRUE, showWarnings=FALSE)

cutadapt_ok <- (system("cutadapt --version", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0)

fq_use <- fq_all  # will be updated to trimmed files if cutadapt succeeds

if (cutadapt_ok) {
  cat("cutadapt found. Trimming primers (linked adapter mode for PacBio)...\n")
  fq_trim <- file.path(trim_dir, basename(fq_all))

  RC <- function(seq) as.character(reverseComplement(DNAString(seq)))
  primer_f_rc <- RC(opt$primer_f)
  primer_r_rc <- RC(opt$primer_r)

  for (i in seq_along(fq_all)) {
    # For PacBio CCS, reads can be in either orientation — use linked adapters
    cmd <- sprintf(
      "cutadapt -g \"%s...%s\" --discard-untrimmed -m %d -M %d --max-n 0 -j %d -o %s %s > /dev/null 2>&1",
      opt$primer_f, primer_r_rc,
      opt$min_length, opt$max_length,
      opt$threads,
      fq_trim[i], fq_all[i]
    )
    rc_cmd <- sprintf(
      "cutadapt -g \"%s...%s\" --discard-untrimmed -m %d -M %d --max-n 0 -j %d -o %s.rc_temp %s > /dev/null 2>&1",
      opt$primer_r, primer_f_rc,
      opt$min_length, opt$max_length,
      opt$threads,
      fq_trim[i], fq_all[i]
    )
    system(cmd)
    system(rc_cmd)
    # Merge both orientations if RC file exists and has reads
    rc_temp <- paste0(fq_trim[i], ".rc_temp")
    if (file.exists(rc_temp) && file.info(rc_temp)$size > 100) {
      system(sprintf("cat %s >> %s && rm %s", rc_temp, fq_trim[i], rc_temp))
    } else if (file.exists(rc_temp)) {
      file.remove(rc_temp)
    }
  }

  # Keep only files that have reads
  exists_trim <- file.exists(fq_trim) & file.info(fq_trim)$size > 100
  if (any(exists_trim)) {
    fq_use <- fq_trim[exists_trim]
    sample_names <- sample_names[exists_trim]
    cat(sprintf("Primer trimming: %d/%d samples passed\n\n", sum(exists_trim), length(exists_trim)))
  } else {
    cat("WARNING: No reads survived primer trimming. Proceeding with untrimmed reads.\n")
    cat("Check that your primer sequences match the data.\n\n")
  }
} else {
  cat("WARNING: cutadapt not found. Proceeding without primer trimming.\n")
  cat("Install: pip install cutadapt\n\n")
}

# ------------------------------------------------------------------
# Step 4: Filter and trim (CCS-specific settings)
# ------------------------------------------------------------------
write_checkpoint("filter_trim", 18, "Quality filtering CCS reads")

filt_dir <- file.path(opt$output_dir, "Filtered")
dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)

fq_filt <- file.path(filt_dir, paste0(sample_names, "_filt.fastq.gz"))
names(fq_filt) <- sample_names

# PacBio CCS: filterAndTrim single-end, no truncLen, use minLen/maxLen
out <- filterAndTrim(
  fq_use, fq_filt,
  minQ   = opt$minQ,
  maxEE  = opt$maxEE,
  minLen = opt$min_length,
  maxLen = opt$max_length,
  maxN   = 0,
  compress = TRUE,
  multithread = opt$threads
)
rownames(out) <- sample_names
cat("Filter results:\n"); print(out); cat("\n")

# Remove samples with 0 reads after filtering
exists_filt <- file.exists(fq_filt) & file.info(fq_filt)$size > 100
fq_filt_ok   <- fq_filt[exists_filt]
snames_ok    <- sample_names[exists_filt]

if (length(fq_filt_ok) == 0) stop("No samples passed quality filtering.")
cat(sprintf("%d/%d samples passed quality filter.\n\n",
            length(fq_filt_ok), length(sample_names)))

# ------------------------------------------------------------------
# Step 5: Learn error rates — PacBio CCS error function
# ------------------------------------------------------------------
write_checkpoint("error_model", 28, "Learning PacBio CCS error rates")

cat("Learning error rates with PacBioErrfun...\n")
cat("(BAND_SIZE=32 for long reads — this may take longer than short-read mode)\n")

err <- learnErrors(
  fq_filt_ok,
  errorEstimationFunction = PacBioErrfun,
  BAND_SIZE   = 32,
  multithread = opt$threads
)

err_dir <- file.path(opt$output_dir, "Dada2")
dir.create(err_dir, recursive=TRUE, showWarnings=FALSE)

pdf(file.path(err_dir, "ErrorModel_PacBio.pdf"), width=10, height=8)
print(plotErrors(err, nominalQ=TRUE))
dev.off()
cat("Error model saved.\n\n")

# ------------------------------------------------------------------
# Step 6: Dereplicate + DADA2 denoise (single-end, CCS)
# ------------------------------------------------------------------
write_checkpoint("dada2_denoise", 45, "DADA2 denoising (single-end)")

pool_val <- switch(opt$pool,
  "TRUE"   = TRUE,
  "FALSE"  = FALSE,
  "pseudo" = "pseudo",
  FALSE
)

derep <- derepFastq(fq_filt_ok)
names(derep) <- snames_ok

# Denoise with PacBio CCS settings
dada_out <- dada(
  derep,
  err        = err,
  BAND_SIZE  = 32,
  pool       = pool_val,
  multithread = opt$threads
)

# Ensure list format even with a single sample
if (inherits(dada_out, "dada")) {
  dada_out <- setNames(list(dada_out), snames_ok)
}

cat("Denoising complete.\n")
for (s in snames_ok) {
  cat(sprintf("  %s: %d ASVs\n", s, length(dada_out[[s]]$sequence)))
}
cat("\n")

# ------------------------------------------------------------------
# Step 7: Build ASV table (no merging for single-end CCS)
# ------------------------------------------------------------------
write_checkpoint("asv_table", 60, "Building ASV table")

seqtab <- makeSequenceTable(dada_out)
cat("ASV table:", nrow(seqtab), "samples x", ncol(seqtab), "ASVs\n")

seq_lengths <- nchar(colnames(seqtab))
cat("ASV length distribution:\n"); print(table(seq_lengths)); cat("\n")

# Length filter (should mostly already be done but double-check)
keep_len <- seq_lengths >= opt$min_length & seq_lengths <= opt$max_length
seqtab_len <- seqtab[, keep_len, drop=FALSE]
cat(sprintf("After length filter [%d-%d bp]: %d ASVs remain\n\n",
            opt$min_length, opt$max_length, ncol(seqtab_len)))

pdf(file.path(err_dir, "ASV_length_distribution.pdf"), width=8, height=5)
df_len <- data.frame(length=seq_lengths)
p <- ggplot(df_len, aes(x=length)) +
  geom_histogram(binwidth=5, fill="#6366F1", color="white", alpha=0.85) +
  geom_vline(xintercept=c(opt$min_length, opt$max_length),
             color="red", linetype="dashed") +
  labs(title=sprintf("PacBio CCS ASV Length Distribution (%s)", opt$region),
       x="ASV length (bp)", y="Count") +
  theme_minimal(base_size=13)
print(p)
dev.off()

# ------------------------------------------------------------------
# Step 8: Remove chimeras
# ------------------------------------------------------------------
write_checkpoint("chimera_removal", 65, "Removing chimeras")

seqtab_nochim <- removeBimeraDenovo(
  seqtab_len,
  method      = "consensus",
  multithread = opt$threads,
  verbose     = TRUE
)
cat(sprintf("After chimera removal: %d ASVs, %.1f%% reads retained\n\n",
            ncol(seqtab_nochim),
            100 * sum(seqtab_nochim) / max(sum(seqtab_len), 1)))

# ------------------------------------------------------------------
# Step 9: Taxonomy
# ------------------------------------------------------------------
write_checkpoint("taxonomy", 72, "Assigning taxonomy")

tax <- NULL
db_path <- opt$dbPath

# Auto-detect if not specified
if (nchar(db_path) == 0 || !file.exists(db_path)) {
  json_path <- opt$db_paths_json
  if (nchar(json_path) == 0)
    json_path <- path.expand("~/r16s-app/backend/databases/db_paths.json")

  if (file.exists(json_path)) {
    cat("Auto-detecting database from db_paths.json...\n")
  }

  # Fall back to scanning default db directory
  db_dir <- path.expand("~/r16s-app/backend/databases")
  db_files <- list.files(db_dir,
                          pattern="\\.(fa|fasta)(\\.gz)?$",
                          full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
  silva_files <- db_files[grepl("silva", db_files, ignore.case=TRUE) &
                           !grepl("species|Species", db_files)]
  if (length(silva_files) > 0) {
    db_path <- silva_files[1]
    cat("Auto-selected:", db_path, "\n")
  }
}

if (nchar(db_path) > 0 && file.exists(db_path)) {
  cat("Assigning taxonomy from:", basename(db_path), "\n")
  tryCatch({
    tax <- assignTaxonomy(seqtab_nochim, db_path,
                          minBoot=opt$minBoot,
                          multithread=opt$threads,
                          verbose=FALSE)
    cat("Taxonomy assigned.\n\n")

    # Try species
    db_dir2 <- dirname(db_path)
    sp_cands <- list.files(db_dir2, pattern="(?i)assignspecies.*\\.fa(\\.gz)?$",
                           full.names=TRUE, perl=TRUE)
    if (length(sp_cands) > 0) {
      tryCatch({
        tax <- addSpecies(tax, sp_cands[1])
        cat("Species added from:", basename(sp_cands[1]), "\n\n")
      }, error=function(e) cat("Species skipped:", e$message, "\n"))
    }
  }, error=function(e) cat("Taxonomy error:", e$message, "\n\n"))
} else {
  cat("No taxonomy database found. Skipping.\n\n")
}

# ------------------------------------------------------------------
# Step 10: Read tracking
# ------------------------------------------------------------------
write_checkpoint("read_tracking", 85, "Building read tracking table")

getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input    = out[, "reads.in"],
  filtered = out[, "reads.out"],
  row.names = sample_names
)
for (s in snames_ok) {
  track[s, "denoised"] <- getN(dada_out[[s]])
}
track[snames_ok, "nonchim"] <- rowSums(seqtab_nochim)
track[is.na(track)] <- 0

write.csv(track, file.path(opt$output_dir, "read_tracking.csv"))
cat("Read tracking:\n"); print(track); cat("\n")

pdf(file.path(qc_dir, "read_tracking_plot.pdf"), width=10, height=6)
steps_show <- intersect(c("input","filtered","denoised","nonchim"), colnames(track))
track_long <- reshape(track[, steps_show, drop=FALSE],
                      varying=steps_show, v.names="reads",
                      timevar="step", times=steps_show, direction="long")
track_long$step <- factor(track_long$step, levels=steps_show)
p <- ggplot(track_long, aes(x=step, y=reads, group=id, color=id)) +
  geom_line(alpha=0.7) + geom_point(size=2.5) +
  labs(title=sprintf("PacBio CCS Read Tracking (%s)", opt$region),
       x="Step", y="Read Count") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1))
print(p); dev.off()

# ------------------------------------------------------------------
# Step 11: Save tables and FASTA
# ------------------------------------------------------------------
write_checkpoint("save_tables", 88, "Saving ASV and taxonomy tables")

asv_ids <- paste0("ASV", seq_len(ncol(seqtab_nochim)))
colnames(seqtab_nochim) <- asv_ids

asv_df <- as.data.frame(seqtab_nochim)
asv_df <- cbind(sample_id=rownames(asv_df), asv_df)
write.csv(asv_df, file.path(opt$output_dir, "asv_table.csv"), row.names=FALSE)

if (!is.null(tax)) {
  rownames(tax) <- asv_ids
  tax_df <- as.data.frame(tax, stringsAsFactors=FALSE)
  tax_df <- cbind(ASV_ID=rownames(tax_df), tax_df)
  write.csv(tax_df, file.path(opt$output_dir, "taxonomy_table.csv"), row.names=FALSE)

  # taxonomy_summary.json for frontend
  asv_counts <- colSums(seqtab_nochim)
  tax_summary_list <- list()
  for (lvl in c("Phylum","Class","Order","Family","Genus")) {
    col_idx <- which(colnames(tax) == lvl)
    if (length(col_idx) == 0) next
    taxon_vec <- tax[, col_idx]
    not_na <- !is.na(taxon_vec) & taxon_vec != ""
    if (sum(not_na) == 0) next
    tbl   <- sort(tapply(asv_counts[not_na], taxon_vec[not_na], sum), decreasing=TRUE)
    total <- sum(tbl)
    tax_summary_list[[lvl]] <- lapply(names(tbl), function(n)
      list(name=n, abundance=round(as.numeric(tbl[n])/total*100, 2)))
  }
  write(toJSON(tax_summary_list, auto_unbox=TRUE),
        file.path(opt$output_dir, "taxonomy_summary.json"))
}

# ASV FASTA
asv_seqs_orig <- colnames(seqtab_len)[keep_len]
if (length(asv_seqs_orig) != length(asv_ids)) {
  asv_seqs_orig <- colnames(seqtab_len)
  if (length(asv_seqs_orig) > length(asv_ids))
    asv_seqs_orig <- asv_seqs_orig[seq_along(asv_ids)]
}
fa_lines <- unlist(mapply(function(id, seq) c(paste0(">", id), seq),
                           asv_ids, asv_seqs_orig))
writeLines(fa_lines, file.path(opt$output_dir, "asvs.fasta"))

# ------------------------------------------------------------------
# Step 12: Taxonomy plots
# ------------------------------------------------------------------
write_checkpoint("taxonomy_plots", 92, "Generating taxonomy plots")

tax_plot_dir <- file.path(opt$output_dir, "Taxonomy_plots")
dir.create(tax_plot_dir, recursive=TRUE, showWarnings=FALSE)

pal20 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
           "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
           "#A6D854","#FFD92F","#E5C494","#B3B3B3","#1B9E77","#D95F02",
           "#7570B3","#E7298A")

make_tax_plot <- function(seqtab, taxonomy, level, topN=opt$topN) {
  if (is.null(taxonomy)) return(NULL)
  col_idx <- which(colnames(taxonomy) == level)
  if (length(col_idx) == 0) return(NULL)

  taxa_vec <- taxonomy[, col_idx]
  taxa_vec[is.na(taxa_vec)] <- "Unclassified"
  taxa_vec <- gsub("^[a-z]__", "", taxa_vec)

  count_mat <- matrix(0, nrow=nrow(seqtab), ncol=length(unique(taxa_vec)))
  rownames(count_mat) <- rownames(seqtab)
  colnames(count_mat) <- unique(taxa_vec)
  for (t in unique(taxa_vec)) {
    idx <- which(taxa_vec == t)
    count_mat[, t] <- rowSums(seqtab[, idx, drop=FALSE])
  }

  total_per_taxa <- colSums(count_mat)
  top_taxa <- names(sort(total_per_taxa, decreasing=TRUE))[1:min(topN, length(total_per_taxa))]
  other    <- rowSums(count_mat[, !colnames(count_mat) %in% top_taxa, drop=FALSE])
  count_top <- count_mat[, top_taxa, drop=FALSE]
  if (any(other > 0)) count_top <- cbind(count_top, Other=other)

  rel <- sweep(count_top, 1, rowSums(count_top), "/") * 100
  rel[is.nan(rel)] <- 0

  df <- data.frame(
    sample = rep(rownames(rel), ncol(rel)),
    taxon  = rep(colnames(rel), each=nrow(rel)),
    rel_ab = as.vector(rel)
  )
  all_taxa <- c(top_taxa, if ("Other" %in% colnames(count_top)) "Other")
  df$taxon <- factor(df$taxon, levels=all_taxa)
  colours  <- setNames(c(pal20[seq_len(min(length(top_taxa), 20))],
                          if ("Other" %in% colnames(count_top)) "#AAAAAA"),
                        c(top_taxa[seq_len(min(length(top_taxa), 20))],
                          if ("Other" %in% colnames(count_top)) "Other"))

  ggplot(df, aes(x=sample, y=rel_ab, fill=taxon)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=colours, na.value="#DDDDDD") +
    labs(title=sprintf("PacBio 16S Relative Abundance — %s (%s)", level, opt$region),
         x="Sample", y="Relative Abundance (%)", fill=level) +
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

for (lvl in c("Phylum","Class","Order","Family","Genus")) {
  tryCatch({
    p <- make_tax_plot(seqtab_nochim, tax, lvl)
    if (!is.null(p)) {
      pdf(file.path(tax_plot_dir, sprintf("taxonomy_%s.pdf", tolower(lvl))),
          width=12, height=6)
      print(p); dev.off()
      cat("Saved:", lvl, "plot\n")
    }
  }, error=function(e) cat("Plot error", lvl, ":", e$message, "\n"))
}

# ------------------------------------------------------------------
# Step 13: Alpha diversity
# ------------------------------------------------------------------
write_checkpoint("alpha_diversity", 96, "Alpha diversity")

alpha_dir <- file.path(opt$output_dir, "Alpha_diversity")
dir.create(alpha_dir, recursive=TRUE, showWarnings=FALSE)

tryCatch({
  obs     <- rowSums(seqtab_nochim > 0)
  shannon <- vegan::diversity(seqtab_nochim, "shannon")
  simpson <- vegan::diversity(seqtab_nochim, "simpson")
  chao1   <- apply(seqtab_nochim, 1, function(x) {
    sobs <- sum(x > 0); f1 <- sum(x == 1); f2 <- sum(x == 2)
    if (f2 > 0) sobs + f1^2/(2*f2) else sobs + f1*(f1-1)/2
  })
  alpha_df <- data.frame(sample_id=names(obs), observed=obs, shannon=shannon,
                          simpson=simpson, chao1=chao1)
  write.csv(alpha_df, file.path(alpha_dir, "alpha_diversity.csv"), row.names=FALSE)

  alpha_long <- reshape(alpha_df,
                         varying=c("observed","shannon","simpson","chao1"),
                         v.names="value", timevar="metric",
                         times=c("observed","shannon","simpson","chao1"),
                         direction="long")
  pdf(file.path(alpha_dir, "alpha_diversity.pdf"), width=10, height=6)
  p <- ggplot(alpha_long, aes(x=metric, y=value)) +
    geom_boxplot(fill="#6366F1", alpha=0.5, outlier.shape=NA) +
    geom_jitter(width=0.15, size=2, alpha=0.8) +
    facet_wrap(~metric, scales="free_y", nrow=1) +
    labs(title=sprintf("PacBio CCS Alpha Diversity (%s)", opt$region), x="", y="Value") +
    theme_minimal(base_size=12)
  print(p); dev.off()
}, error=function(e) cat("Alpha diversity warning:", e$message, "\n"))

# ------------------------------------------------------------------
# Step 14: Summary JSON
# ------------------------------------------------------------------
write_checkpoint("summary", 99, "Writing summary")

summary_list <- list(
  pipeline        = "PacBio_CCS_16S",
  job_name        = opt$job_name,
  region          = opt$region,
  completed_at    = format(Sys.time()),
  n_samples       = nrow(seqtab_nochim),
  n_asvs          = ncol(seqtab_nochim),
  total_reads     = sum(seqtab_nochim),
  length_range    = paste0(opt$min_length, "-", opt$max_length, " bp"),
  error_fn        = "PacBioErrfun",
  band_size       = 32,
  has_taxonomy    = !is.null(tax),
  status          = "completed",
  output_files    = list.files(opt$output_dir, full.names=FALSE, recursive=TRUE)
)
write(toJSON(summary_list, auto_unbox=TRUE, pretty=TRUE),
      file.path(opt$output_dir, "summary.json"))

cat("\n=== PacBio CCS Pipeline Complete ===\n")
cat("Region  :", opt$region, "\n")
cat("Samples :", nrow(seqtab_nochim), "\n")
cat("ASVs    :", ncol(seqtab_nochim), "\n")
cat("Reads   :", sum(seqtab_nochim), "\n")
cat("Output  :", opt$output_dir, "\n")
cat("Finished:", format(Sys.time()), "\n")

write_checkpoint("done", 100, "PacBio CCS pipeline complete")

sink(type="output"); sink(type="message"); close(con)
