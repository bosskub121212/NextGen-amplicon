#!/usr/bin/env Rscript
# =============================================================
#  NextGen-Amplicon -- COX1 Metabarcoding Pipeline
#  Marker: Cytochrome c Oxidase Subunit I (COX1 / CO1)
#
#  Workflow:
#    1. Primer / adapter trimming (cutadapt via system call)
#    2. DADA2 quality filter, learn errors, denoise
#    3. Merge paired-end reads
#    4. Chimera removal (consensus method)
#    5. Codon translation filter  -- remove NUMTs / pseudogenes
#    6. LULU post-clustering curation  -- remove artefact ASVs
#    7. MIDORI2 taxonomy assignment (DADA2 assignTaxonomy)
#    8. Standard output files (ASV table, taxonomy table, plots)
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
              help="Directory containing paired-end FASTQ files"),
  make_option("--output_dir",   type="character", default=NULL,
              help="Directory for all output files"),
  make_option("--db_paths",     type="character", default=NULL,
              help="Path to db_paths.json (optional; auto-detected if absent)"),
  make_option("--primer_f",     type="character", default="GGWACWGGWTGAACWGTWTAYCCYCC",
              help="Forward primer sequence (default: mlCOIintF)"),
  make_option("--primer_r",     type="character", default="TANACYTCNGGRTGNCCRAARAAYCA",
              help="Reverse primer sequence (default: jgHCO2198)"),
  make_option("--truncLen_f",   type="integer",   default=230,
              help="Truncate forward reads at this length (0 = no truncation)"),
  make_option("--truncLen_r",   type="integer",   default=200,
              help="Truncate reverse reads at this length (0 = no truncation)"),
  make_option("--maxEE_f",      type="numeric",   default=2,
              help="Maximum expected errors (forward)"),
  make_option("--maxEE_r",      type="numeric",   default=5,
              help="Maximum expected errors (reverse); COX1 reverse often lower quality"),
  make_option("--min_overlap",  type="integer",   default=20,
              help="Minimum overlap for paired-end merging"),
  make_option("--threads",      type="integer",   default=4,
              help="Number of CPU threads"),
  make_option("--job_name",     type="character", default="COX1_job",
              help="Job / sample set name"),
  make_option("--codon_table",  type="integer",   default=5,
              help="Genetic code for translation check (5=invertebrate mt; 2=vertebrate mt)"),
  make_option("--min_length",   type="integer",   default=300,
              help="Minimum ASV length to keep after filtering"),
  make_option("--max_length",   type="integer",   default=330,
              help="Maximum ASV length (expected COX1 amplicon ~313 bp for mlCOIintF/jgHCO2198)"),
  make_option("--topN",         type="integer",   default=10,
              help="Top N taxa for bar charts"),
  make_option("--metadata",     type="character", default=NULL,
              help="Optional: path to metadata TSV (sample_id in first column)"),
  make_option("--lulu",         type="logical",   default=TRUE,
              help="Run LULU post-clustering curation (requires BLAST)"),
  make_option("--checkpoint_file", type="character", default=NULL,
              help="Write checkpoint progress here for backend polling")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input_dir))  stop("--input_dir is required")
if (is.null(opt$output_dir)) stop("--output_dir is required")

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)
log_file <- file.path(opt$output_dir, "pipeline.log")
con <- file(log_file, open="wt")
sink(con, type="output")
sink(con, type="message")

cat("=== NextGen-Amplicon COX1 Pipeline ===\n")
cat("Job:", opt$job_name, "\n")
cat("Input:", opt$input_dir, "\n")
cat("Output:", opt$output_dir, "\n")
cat("Threads:", opt$threads, "\n")
cat("Primer F:", opt$primer_f, "\n")
cat("Primer R:", opt$primer_r, "\n")
cat("Codon table:", opt$codon_table, "(5=invertebrate mt, 2=vertebrate mt)\n")
cat("Started:", format(Sys.time()), "\n\n")

# ------------------------------------------------------------------
# Helper: checkpoint
# ------------------------------------------------------------------
write_checkpoint <- function(step, pct, msg="") {
  if (!is.null(opt$checkpoint_file)) {
    cp <- list(step=step, pct=pct, msg=msg, ts=format(Sys.time()))
    write(toJSON(cp, auto_unbox=TRUE), opt$checkpoint_file)
  }
  cat(sprintf("[CHECKPOINT] %s (%.0f%%) %s\n", step, pct, msg))
}

write_checkpoint("init", 2, "COX1 pipeline starting")

# ------------------------------------------------------------------
# Helper: locate MIDORI2 database
# ------------------------------------------------------------------
find_midori2_db <- function(db_paths_file=NULL) {
  db_genus <- NULL; db_sp <- NULL

  if (!is.null(db_paths_file) && file.exists(db_paths_file)) {
    paths <- fromJSON(db_paths_file)
    if (!is.null(paths$MIDORI2_COX1) && file.exists(paths$MIDORI2_COX1))
      db_genus <- paths$MIDORI2_COX1
    if (!is.null(paths$MIDORI2_COX1_sp) && file.exists(paths$MIDORI2_COX1_sp))
      db_sp <- paths$MIDORI2_COX1_sp
  }

  auto_search <- c(
    "~/r16s-app/backend/databases/MIDORI2/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz",
    "/home/databases/MIDORI2/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz",
    "~/databases/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz"
  )
  if (is.null(db_genus)) {
    for (p in auto_search) {
      ep <- path.expand(p)
      if (file.exists(ep)) { db_genus <- ep; break }
    }
  }

  sp_search <- gsub("DADA2\\.fasta\\.gz", "DADA2_sp.fasta.gz", auto_search)
  if (is.null(db_sp)) {
    for (p in sp_search) {
      ep <- path.expand(p)
      if (file.exists(ep)) { db_sp <- ep; break }
    }
  }

  list(genus=db_genus, species=db_sp)
}

dbs <- find_midori2_db(opt$db_paths)
cat("MIDORI2 genus DB :", ifelse(is.null(dbs$genus), "NOT FOUND", dbs$genus), "\n")
cat("MIDORI2 species DB:", ifelse(is.null(dbs$sp),   "NOT FOUND (addSpecies step will be skipped)", dbs$sp), "\n\n")

# ------------------------------------------------------------------
# Step 1: Detect FASTQ files
# ------------------------------------------------------------------
write_checkpoint("detect_files", 4, "Detecting FASTQ files")

fq_all <- list.files(opt$input_dir, pattern="\\.(fastq|fastq\\.gz|fq|fq\\.gz)$",
                     full.names=TRUE, recursive=FALSE)

fnFs <- sort(fq_all[grepl("_R1", fq_all) | grepl("_1\\.f", fq_all)])
fnRs <- sort(fq_all[grepl("_R2", fq_all) | grepl("_2\\.f", fq_all)])

if (length(fnFs) == 0) stop("No R1 FASTQ files found in: ", opt$input_dir)
if (length(fnFs) != length(fnRs)) stop("Mismatched R1/R2 counts")

sample_names <- sub("_R1.*", "", basename(fnFs))
sample_names <- sub("_1\\.f.*",  "", sample_names)
cat("Samples detected:", length(sample_names), "\n")
cat(paste(" ", sample_names, collapse="\n"), "\n\n")

# ------------------------------------------------------------------
# Step 2: Quality profiles
# ------------------------------------------------------------------
write_checkpoint("quality_plots", 6, "Plotting quality profiles")

qc_dir <- file.path(opt$output_dir, "QC_plots")
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

tryCatch({
  n_show <- min(length(fnFs), 4)
  pdf(file.path(qc_dir, "quality_profile_R1.pdf"), width=10, height=6)
  print(plotQualityProfile(fnFs[1:n_show]))
  dev.off()
  pdf(file.path(qc_dir, "quality_profile_R2.pdf"), width=10, height=6)
  print(plotQualityProfile(fnRs[1:n_show]))
  dev.off()
}, error=function(e) cat("Quality plot warning:", e$message, "\n"))

# ------------------------------------------------------------------
# Step 3: Primer trimming with cutadapt (if available)
# ------------------------------------------------------------------
write_checkpoint("primer_trim", 10, "Trimming primers with cutadapt")

trim_dir <- file.path(opt$output_dir, "Trimmed")
dir.create(trim_dir, recursive=TRUE, showWarnings=FALSE)

cutadapt_ok <- (system("cutadapt --version", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0)

if (cutadapt_ok) {
  cat("cutadapt found, trimming primers...\n")
  fnFs_trim <- file.path(trim_dir, basename(fnFs))
  fnRs_trim <- file.path(trim_dir, basename(fnRs))

  RC <- function(seq) as.character(reverseComplement(DNAString(seq)))
  primer_f_rc <- RC(opt$primer_f)
  primer_r_rc <- RC(opt$primer_r)

  for (i in seq_along(fnFs)) {
    cmd <- sprintf(
      "cutadapt -g %s -a %s -G %s -A %s --discard-untrimmed -m 50 -j %d -o %s -p %s %s %s > /dev/null 2>&1",
      opt$primer_f, primer_r_rc,
      opt$primer_r, primer_f_rc,
      opt$threads,
      fnFs_trim[i], fnRs_trim[i],
      fnFs[i], fnRs[i]
    )
    system(cmd)
  }
  # Use trimmed files going forward
  fnFs <- fnFs_trim
  fnRs <- fnRs_trim
  cat("Primer trimming done.\n\n")
} else {
  cat("WARNING: cutadapt not found. Proceeding without explicit primer trimming.\n")
  cat("Install with: pip install cutadapt\n\n")
}

# ------------------------------------------------------------------
# Step 4: Filter and trim
# ------------------------------------------------------------------
write_checkpoint("filter_trim", 18, "Quality filtering reads")

filt_dir <- file.path(opt$output_dir, "Filtered")
dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)

fnFs_filt <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
fnRs_filt <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(fnFs_filt) <- sample_names
names(fnRs_filt) <- sample_names

# COX1 uses truncLen to keep consistent amplicon size; 0 = no truncation
truncLen_vec <- c(opt$truncLen_f, opt$truncLen_r)
if (all(truncLen_vec == 0)) {
  out <- filterAndTrim(fnFs, fnFs_filt, fnRs, fnRs_filt,
                       maxEE=c(opt$maxEE_f, opt$maxEE_r),
                       truncQ=2, minLen=opt$min_length, rm.phix=TRUE,
                       compress=TRUE, multithread=opt$threads)
} else {
  out <- filterAndTrim(fnFs, fnFs_filt, fnRs, fnRs_filt,
                       truncLen=truncLen_vec,
                       maxEE=c(opt$maxEE_f, opt$maxEE_r),
                       truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=opt$threads)
}

cat("Filter results:\n")
print(out)
cat("\n")

# ------------------------------------------------------------------
# Step 5: Learn error rates
# ------------------------------------------------------------------
write_checkpoint("error_model", 30, "Learning error rates")

errF <- learnErrors(fnFs_filt[file.exists(fnFs_filt)], multithread=opt$threads)
errR <- learnErrors(fnRs_filt[file.exists(fnRs_filt)], multithread=opt$threads)

err_dir <- file.path(opt$output_dir, "Dada2")
dir.create(err_dir, recursive=TRUE, showWarnings=FALSE)

pdf(file.path(err_dir, "ErrorModel_R1.pdf"), width=10, height=8)
print(plotErrors(errF, nominalQ=TRUE))
dev.off()
pdf(file.path(err_dir, "ErrorModel_R2.pdf"), width=10, height=8)
print(plotErrors(errR, nominalQ=TRUE))
dev.off()

# ------------------------------------------------------------------
# Step 6: Dereplicate + DADA2 sample inference
# ------------------------------------------------------------------
write_checkpoint("dada2_denoise", 45, "DADA2 denoising")

exists_filt <- file.exists(fnFs_filt) & file.exists(fnRs_filt)
fnFs_ok <- fnFs_filt[exists_filt]
fnRs_ok <- fnRs_filt[exists_filt]
snames_ok <- sample_names[exists_filt]

derepF <- derepFastq(fnFs_ok)
derepR <- derepFastq(fnRs_ok)
names(derepF) <- snames_ok
names(derepR) <- snames_ok

dadaF <- dada(derepF, err=errF, multithread=opt$threads)
dadaR <- dada(derepR, err=errR, multithread=opt$threads)

# ------------------------------------------------------------------
# Step 7: Merge paired reads
# ------------------------------------------------------------------
write_checkpoint("merge_pairs", 55, "Merging paired-end reads")

mergers <- mergePairs(dadaF, derepF, dadaR, derepR,
                      minOverlap=opt$min_overlap,
                      maxMismatch=0,
                      verbose=FALSE)

# ------------------------------------------------------------------
# Step 8: Construct ASV table
# ------------------------------------------------------------------
write_checkpoint("asv_table", 60, "Building ASV table")

seqtab <- makeSequenceTable(mergers)
cat("ASV table dimensions:", dim(seqtab), "\n")

# Length distribution
seq_lengths <- nchar(colnames(seqtab))
cat("ASV length distribution:\n")
print(table(seq_lengths))

# Filter by expected amplicon length range
keep_len <- seq_lengths >= opt$min_length & seq_lengths <= opt$max_length
seqtab_len <- seqtab[, keep_len, drop=FALSE]
cat(sprintf("After length filter [%d-%d bp]: %d ASVs\n",
            opt$min_length, opt$max_length, ncol(seqtab_len)))

pdf(file.path(err_dir, "ASV_length_distribution.pdf"), width=8, height=5)
df_len <- data.frame(length=seq_lengths)
p <- ggplot(df_len, aes(x=length)) +
  geom_histogram(binwidth=1, fill="#E07B39", color="white", alpha=0.85) +
  geom_vline(xintercept=c(opt$min_length, opt$max_length), color="red", linetype="dashed") +
  labs(title="COX1 ASV Length Distribution", x="ASV length (bp)", y="Count") +
  theme_minimal(base_size=13)
print(p)
dev.off()

# ------------------------------------------------------------------
# Step 9: Remove chimeras
# ------------------------------------------------------------------
write_checkpoint("chimera_removal", 65, "Removing chimeras")

seqtab_nochim <- removeBimeraDenovo(seqtab_len, method="consensus",
                                    multithread=opt$threads, verbose=TRUE)
cat(sprintf("After chimera removal: %d ASVs (%.1f%% reads retained)\n",
            ncol(seqtab_nochim),
            100 * sum(seqtab_nochim) / sum(seqtab_len)))

# ------------------------------------------------------------------
# Step 10: Codon translation filter (remove NUMTs / pseudogenes)
# ------------------------------------------------------------------
write_checkpoint("codon_filter", 70, "Codon translation filter (NUMT removal)")

cat("\n--- Codon translation filter (genetic code:", opt$codon_table, ") ---\n")

translate_check <- function(seq, gc=opt$codon_table) {
  # Check all 3 reading frames; return TRUE if at least 1 has no stop codons
  dna <- DNAString(seq)
  for (frame in 1:3) {
    sub_seq <- subseq(dna, start=frame)
    # Trim to multiple of 3
    trim_len <- floor(nchar(sub_seq) / 3) * 3
    if (trim_len < 150) next  # too short
    sub_seq <- subseq(sub_seq, end=trim_len)
    tryCatch({
      aa <- translate(sub_seq, genetic.code=getGeneticCode(as.character(gc)))
      if (!any(aa == "*")) return(TRUE)
    }, error=function(e) NULL)
  }
  return(FALSE)
}

asvs <- colnames(seqtab_nochim)
cat("Checking", length(asvs), "ASVs for premature stop codons...\n")

# Process in batches for speed reporting
valid_codon <- logical(length(asvs))
for (i in seq_along(asvs)) {
  valid_codon[i] <- translate_check(asvs[i])
  if (i %% 100 == 0) cat(sprintf("  Checked %d / %d ASVs\r", i, length(asvs)))
}
cat("\n")

n_pass <- sum(valid_codon)
n_fail <- sum(!valid_codon)
cat(sprintf("Translation filter: %d pass, %d fail (likely NUMTs)\n", n_pass, n_fail))

seqtab_translated <- seqtab_nochim[, valid_codon, drop=FALSE]

if (ncol(seqtab_translated) == 0) {
  cat("WARNING: All ASVs failed translation filter. Skipping filter (check codon table).\n")
  seqtab_translated <- seqtab_nochim
}

# ------------------------------------------------------------------
# Step 11: LULU post-clustering curation
# ------------------------------------------------------------------
write_checkpoint("lulu_curation", 76, "LULU post-clustering curation")

seqtab_final <- seqtab_translated

if (opt$lulu && requireNamespace("lulu", quietly=TRUE)) {
  blast_ok <- (system("makeblastdb -version", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0)
  if (!blast_ok) {
    cat("WARNING: BLAST not found. Skipping LULU curation.\n")
    cat("Install with: sudo apt-get install ncbi-blast+\n")
  } else {
    cat("Running LULU curation...\n")
    asv_seqs <- colnames(seqtab_translated)
    asv_ids   <- paste0("ASV", seq_along(asv_seqs))

    # Write FASTA for BLAST
    lulu_dir   <- file.path(opt$output_dir, "LULU")
    dir.create(lulu_dir, recursive=TRUE, showWarnings=FALSE)
    fasta_path <- file.path(lulu_dir, "asvs.fasta")
    db_path    <- file.path(lulu_dir, "asv_blastdb")
    match_path <- file.path(lulu_dir, "match_list.txt")

    fa_lines <- unlist(mapply(function(id, seq) c(paste0(">", id), seq),
                              asv_ids, asv_seqs))
    writeLines(fa_lines, fasta_path)

    system(sprintf("makeblastdb -in %s -dbtype nucl -out %s > /dev/null 2>&1",
                   fasta_path, db_path))
    system(sprintf(
      "blastn -db %s -query %s -out %s -perc_identity 84 -qcov_hsp_perc 80 -num_alignments 10 -outfmt '6 qseqid sseqid pident' -num_threads %d > /dev/null 2>&1",
      db_path, fasta_path, match_path, opt$threads
    ))

    if (file.exists(match_path) && file.info(match_path)$size > 0) {
      match_list <- read.table(match_path, header=FALSE, col.names=c("child","father","pident"),
                               stringsAsFactors=FALSE)

      otu_table_lulu <- t(seqtab_translated)
      rownames(otu_table_lulu) <- asv_ids

      tryCatch({
        curated <- lulu::lulu(otu_table_lulu, match_list)
        kept_ids <- curated$curated_table |> rownames()
        kept_idx <- match(kept_ids, asv_ids)
        seqtab_final <- seqtab_translated[, kept_idx, drop=FALSE]
        cat(sprintf("LULU: %d -> %d ASVs after curation (removed %d)\n",
                    ncol(seqtab_translated), ncol(seqtab_final),
                    ncol(seqtab_translated) - ncol(seqtab_final)))
      }, error=function(e) {
        cat("LULU error:", e$message, "— skipping\n")
        seqtab_final <<- seqtab_translated
      })
    } else {
      cat("LULU: BLAST match list empty — skipping curation\n")
    }
  }
} else if (opt$lulu) {
  cat("LULU package not installed. Run install_extensions.sh first.\n")
}

# ------------------------------------------------------------------
# Step 12: Taxonomy with MIDORI2
# ------------------------------------------------------------------
write_checkpoint("taxonomy", 82, "Assigning taxonomy with MIDORI2")

if (!is.null(dbs$genus) && file.exists(dbs$genus)) {
  cat("Assigning taxonomy with MIDORI2 (genus level)...\n")
  # Note: taxLevels is omitted so DADA2 auto-detects levels from the MIDORI2
  # FASTA headers (which contain 7 levels: Kingdom→Species). Supplying an
  # explicit 6-level taxLevels vector causes an internal "subscript out of
  # bounds" error in older DADA2 versions when the DB has a different count.
  tax <- tryCatch({
    assignTaxonomy(seqtab_final, dbs$genus,
                   multithread = opt$threads,
                   minBoot     = 50,
                   verbose     = FALSE)
  }, error = function(e) {
    cat("WARNING: assignTaxonomy error:", e$message, "\n")
    cat("Retrying without minBoot...\n")
    tryCatch({
      assignTaxonomy(seqtab_final, dbs$genus,
                     multithread = opt$threads,
                     verbose     = FALSE)
    }, error = function(e2) {
      cat("ERROR: assignTaxonomy failed completely:", e2$message, "\n")
      cat("Creating empty taxonomy table and continuing.\n")
      matrix(NA_character_,
             nrow = ncol(seqtab_final), ncol = 7,
             dimnames = list(
               paste0("ASV", seq_len(ncol(seqtab_final))),
               c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
             ))
    })
  })

  if (!is.null(dbs$species) && file.exists(dbs$species)) {
    cat("Adding species-level assignments...\n")
    tryCatch({
      tax <- addSpecies(tax, dbs$species, allowMultiple=FALSE)
    }, error = function(e) {
      cat("WARNING: addSpecies failed:", e$message, "— skipping\n")
    })
  }

  rownames(tax) <- paste0("ASV", seq_len(nrow(tax)))
} else {
  cat("WARNING: MIDORI2 database not found. Creating empty taxonomy.\n")
  cat("Run download_databases.sh to get the database.\n")
  tax <- matrix(NA, nrow=ncol(seqtab_final),
                ncol=7,
                dimnames=list(paste0("ASV", seq_len(ncol(seqtab_final))),
                              c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
}

# Rename ASVs in seqtab
colnames(seqtab_final) <- paste0("ASV", seq_len(ncol(seqtab_final)))

# ------------------------------------------------------------------
# Step 13: Read tracking table
# ------------------------------------------------------------------
write_checkpoint("read_tracking", 87, "Building read tracking table")

getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input   = out[, 1],
  filtered= out[, 2],
  row.names = rownames(out)
)
# Add denoised, merged
for (s in snames_ok) {
  track[s, "denoised_F"] <- getN(dadaF[[s]])
  track[s, "denoised_R"] <- getN(dadaR[[s]])
  track[s, "merged"]     <- getN(mergers[[s]])
}
track[snames_ok, "nonchim"] <- rowSums(seqtab_nochim)
track[snames_ok, "length_filtered"] <- rowSums(seqtab_len)
track[snames_ok, "codon_filtered"] <- rowSums(seqtab_translated)
track[snames_ok, "lulu_curated"] <- rowSums(seqtab_final)
track[is.na(track)] <- 0

write.csv(track, file.path(opt$output_dir, "read_tracking.csv"))
cat("Read tracking:\n"); print(track); cat("\n")

# Track plot
pdf(file.path(qc_dir, "read_tracking_plot.pdf"), width=10, height=6)
steps_show <- intersect(c("input","filtered","merged","lulu_curated"), colnames(track))
track_long <- reshape(track[, steps_show, drop=FALSE],
                      varying=steps_show, v.names="reads",
                      timevar="step", times=steps_show, direction="long")
track_long$step <- factor(track_long$step, levels=steps_show)
p <- ggplot(track_long, aes(x=step, y=reads, group=id, color=id)) +
  geom_line(alpha=0.7) + geom_point(size=2) +
  labs(title="COX1 Read Tracking", x="Pipeline Step", y="Read Count") +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
print(p)
dev.off()

# ------------------------------------------------------------------
# Step 14: Save core tables
# ------------------------------------------------------------------
write_checkpoint("save_tables", 90, "Saving ASV and taxonomy tables")

# ASV table: samples x ASVs
asv_df <- as.data.frame(seqtab_final)
asv_df <- cbind(sample_id=rownames(asv_df), asv_df)
write.csv(asv_df, file.path(opt$output_dir, "asv_table.csv"), row.names=FALSE)

# Taxonomy table
tax_df <- as.data.frame(tax, stringsAsFactors=FALSE)
tax_df <- cbind(ASV_ID=rownames(tax_df), tax_df)
write.csv(tax_df, file.path(opt$output_dir, "taxonomy_table.csv"), row.names=FALSE)

# ASV FASTA (for external tools)
asv_seqs_final <- colnames(seqtab_translated)[valid_codon | !opt$lulu]
if (length(asv_seqs_final) != nrow(tax_df)) {
  asv_seqs_final <- colnames(seqtab_translated)
}
fa_out <- file.path(opt$output_dir, "asvs.fasta")
asv_ids_out <- rownames(tax_df)
fa_lines <- unlist(mapply(function(id, seq) c(paste0(">", id), seq),
                          asv_ids_out, colnames(seqtab_final)))
# Note: colnames of seqtab_final are now ASV1, ASV2... but we need actual sequences
# Get sequences from rownames of tax (original seqtab columns)
asv_sequences <- colnames(seqtab_translated)
if (length(asv_sequences) > length(asv_ids_out)) {
  asv_sequences <- asv_sequences[valid_codon]
}
if (length(asv_sequences) == length(asv_ids_out)) {
  fa_out_lines <- unlist(mapply(function(id, seq) c(paste0(">", id), seq),
                                asv_ids_out, asv_sequences))
  writeLines(fa_out_lines, fa_out)
  cat("ASV FASTA saved:", fa_out, "\n")
}

# ------------------------------------------------------------------
# Step 15: Taxonomy plots
# ------------------------------------------------------------------
write_checkpoint("taxonomy_plots", 93, "Generating taxonomy plots")

tax_plot_dir <- file.path(opt$output_dir, "Taxonomy_plots")
dir.create(tax_plot_dir, recursive=TRUE, showWarnings=FALSE)

# Colour palette (20 distinct colours)
pal20 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
           "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
           "#A6D854","#FFD92F","#E5C494","#B3B3B3","#1B9E77","#D95F02",
           "#7570B3","#E7298A")

make_tax_plot <- function(seqtab, taxonomy, level, topN=opt$topN, title_suffix="") {
  col_idx <- which(colnames(taxonomy) == level)
  if (length(col_idx) == 0) return(NULL)

  taxa_vec <- taxonomy[, col_idx]
  taxa_vec[is.na(taxa_vec)] <- "Unclassified"
  taxa_vec <- gsub("^[a-z]__", "", taxa_vec)

  # Aggregate by taxa
  count_mat <- matrix(0, nrow=nrow(seqtab), ncol=length(unique(taxa_vec)))
  rownames(count_mat) <- rownames(seqtab)
  colnames(count_mat) <- unique(taxa_vec)
  for (t in unique(taxa_vec)) {
    idx <- which(taxa_vec == t)
    count_mat[, t] <- rowSums(seqtab[, idx, drop=FALSE])
  }

  # Top N
  total_per_taxa <- colSums(count_mat)
  top_taxa <- names(sort(total_per_taxa, decreasing=TRUE))[1:min(topN, length(total_per_taxa))]
  other <- rowSums(count_mat[, !colnames(count_mat) %in% top_taxa, drop=FALSE])
  count_top <- count_mat[, top_taxa, drop=FALSE]
  if (any(other > 0)) {
    count_top <- cbind(count_top, Other=other)
  }

  # Relative abundance
  rel <- sweep(count_top, 1, rowSums(count_top), "/") * 100
  rel[is.nan(rel)] <- 0

  df <- data.frame(
    sample = rep(rownames(rel), ncol(rel)),
    taxon  = rep(colnames(rel), each=nrow(rel)),
    rel_ab = as.vector(rel)
  )
  df$taxon <- factor(df$taxon, levels=c(top_taxa, if("Other" %in% colnames(count_top)) "Other"))

  colours <- setNames(
    c(pal20[seq_len(min(length(top_taxa), 20))],
      if("Other" %in% colnames(count_top)) "#AAAAAA"),
    c(top_taxa[seq_len(min(length(top_taxa), 20))],
      if("Other" %in% colnames(count_top)) "Other")
  )

  p <- ggplot(df, aes(x=sample, y=rel_ab, fill=taxon)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=colours, na.value="#DDDDDD") +
    labs(title=paste("COX1 Relative Abundance —", level, title_suffix),
         x="Sample", y="Relative Abundance (%)", fill=level) +
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="right",
          legend.key.size=unit(0.4, "cm"))
  p
}

for (lvl in c("Phylum","Class","Order","Family","Genus")) {
  tryCatch({
    p <- make_tax_plot(seqtab_final, tax, lvl)
    if (!is.null(p)) {
      pdf(file.path(tax_plot_dir, sprintf("taxonomy_%s.pdf", tolower(lvl))),
          width=12, height=6)
      print(p)
      dev.off()
      cat("Saved:", lvl, "taxonomy plot\n")
    }
  }, error=function(e) cat("Plot error for", lvl, ":", e$message, "\n"))
}

# ------------------------------------------------------------------
# Step 16: Alpha diversity
# ------------------------------------------------------------------
write_checkpoint("alpha_diversity", 96, "Alpha diversity")

alpha_dir <- file.path(opt$output_dir, "Alpha_diversity")
dir.create(alpha_dir, recursive=TRUE, showWarnings=FALSE)

tryCatch({
  obs      <- rowSums(seqtab_final > 0)
  shannon  <- vegan::diversity(seqtab_final, "shannon")
  simpson  <- vegan::diversity(seqtab_final, "simpson")
  chao1    <- apply(seqtab_final, 1, function(x) {
    sobs <- sum(x > 0); f1 <- sum(x == 1); f2 <- sum(x == 2)
    if (f2 > 0) sobs + f1^2 / (2 * f2) else sobs + f1 * (f1 - 1) / 2
  })

  alpha_df <- data.frame(sample_id=names(obs), observed=obs, shannon=shannon,
                         simpson=simpson, chao1=chao1)
  write.csv(alpha_df, file.path(alpha_dir, "alpha_diversity.csv"), row.names=FALSE)

  # Box/dot plot
  alpha_long <- reshape(alpha_df, varying=c("observed","shannon","simpson","chao1"),
                        v.names="value", timevar="metric",
                        times=c("observed","shannon","simpson","chao1"), direction="long")
  alpha_long$metric <- factor(alpha_long$metric, levels=c("observed","shannon","simpson","chao1"))

  pdf(file.path(alpha_dir, "alpha_diversity.pdf"), width=10, height=6)
  p <- ggplot(alpha_long, aes(x=metric, y=value)) +
    geom_boxplot(fill="#E07B39", alpha=0.5, outlier.shape=NA) +
    geom_jitter(width=0.15, size=2, alpha=0.8) +
    facet_wrap(~metric, scales="free_y", nrow=1) +
    labs(title="COX1 Alpha Diversity", x="", y="Value") +
    theme_minimal(base_size=12)
  print(p)
  dev.off()
}, error=function(e) cat("Alpha diversity warning:", e$message, "\n"))

# ------------------------------------------------------------------
# Step 17: Summary JSON (for backend)
# ------------------------------------------------------------------
write_checkpoint("summary", 99, "Writing pipeline summary")

n_samples  <- nrow(seqtab_final)
n_asvs     <- ncol(seqtab_final)
total_reads <- sum(seqtab_final)
n_numt     <- n_fail
n_lulu_rm  <- ncol(seqtab_translated) - ncol(seqtab_final)

summary_list <- list(
  pipeline       = "COX1",
  job_name       = opt$job_name,
  completed_at   = format(Sys.time()),
  n_samples      = n_samples,
  n_asvs         = n_asvs,
  total_reads    = total_reads,
  codon_table    = opt$codon_table,
  numt_removed   = n_numt,
  lulu_removed   = max(0, n_lulu_rm),
  taxonomy_db    = ifelse(is.null(dbs$genus), "not_found", basename(dbs$genus)),
  primer_f       = opt$primer_f,
  primer_r       = opt$primer_r,
  length_range   = paste0(opt$min_length, "-", opt$max_length, " bp")
)
write(toJSON(summary_list, auto_unbox=TRUE, pretty=TRUE),
      file.path(opt$output_dir, "pipeline_summary.json"))

cat("\n=== COX1 Pipeline Complete ===\n")
cat("Samples:            ", n_samples, "\n")
cat("Final ASVs:         ", n_asvs, "\n")
cat("Total reads (final):", total_reads, "\n")
cat("NUMTs removed:      ", n_numt, "\n")
cat("LULU removed:       ", max(0, n_lulu_rm), "\n")
cat("Output dir:         ", opt$output_dir, "\n")
cat("Finished:", format(Sys.time()), "\n")

write_checkpoint("done", 100, "COX1 pipeline complete")

sink(type="output"); sink(type="message"); close(con)
