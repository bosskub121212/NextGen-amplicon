# ============================================================
#  its_pipeline.R  —  ITS1 / ITS2 Fungal Amplicon Pipeline
#  Based on DADA2 ITS workflow (benjjneb.github.io/dada2/ITS_workflow.html)
#
#  Key differences from 16S pipeline:
#   - NO truncLen (ITS has variable length 200-600 bp)
#   - ITSxpress for precise ITS1/ITS2 boundary trimming
#   - UNITE database for taxonomy
#   - FUNGuildR for ecological guild annotation
# ============================================================

tryCatch({
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset="")
  if (nchar(r_libs_user)>0 && dir.exists(r_libs_user))
    .libPaths(unique(c(r_libs_user, .libPaths())))
  home_dir <- Sys.getenv("HOME", unset=path.expand("~"))
  for (rp in c(file.path(home_dir,"R","library"),
               file.path(home_dir,"R","x86_64-pc-linux-gnu-library",
                         paste(R.version$major, substr(R.version$minor,1,1), sep=".")))) {
    if (dir.exists(rp)) .libPaths(unique(c(rp, .libPaths())))
  }
}, error=function(e) NULL)

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  # ── Names match what backend/main.py passes ──────────────────────────────
  make_option("--input_dir",   type="character", help="Input directory with FASTQ files"),
  make_option("--output_dir",  type="character", help="Output directory"),
  make_option("--its_region",  type="character", default="ITS2",
              help="ITS region: ITS1 or ITS2 [default: ITS2]"),
  make_option("--db_paths",    type="character", default="",
              help="Path to db_paths.json (auto-detected if absent)"),
  make_option("--maxEE_f",     type="double",   default=2.0),
  make_option("--maxEE_r",     type="double",   default=2.0),
  make_option("--threads",     type="integer",  default=4),
  make_option("--job_name",    type="character", default="ITS_job"),
  make_option("--primer_f",          type="character", default="",
              help="Forward primer sequence for cutadapt trimming"),
  make_option("--primer_r",          type="character", default="",
              help="Reverse primer sequence for cutadapt trimming"),
  make_option("--error_rate",        type="double",    default=0.1,
              help="cutadapt -e: fraction of mismatches allowed"),
  make_option("--min_overlap",       type="integer",   default=3,
              help="cutadapt -O: minimum overlap bases"),
  make_option("--discard_untrimmed", type="logical",   default=FALSE,
              help="Discard reads where primer was not found"),
  # ── Kept for backwards-compat / direct CLI use ───────────────────────────
  make_option("--taxDatabase", type="character", default="UNITE"),
  make_option("--dbPath",      type="character", default=""),
  make_option("--minLen",      type="integer",  default=50),
  make_option("--minBoot",     type="integer",  default=50),
  make_option("--topN",        type="integer",  default=30),
  make_option("--metadata",    type="character", default=""),
  make_option("--marker",      type="character", default="ITS2"),
  make_option("--nbases",      type="double",   default=1e8),
  make_option("--pool",        type="character", default="FALSE")
)
opt <- parse_args(OptionParser(option_list=option_list))

input_dir  <- opt$input_dir
output_dir <- opt$output_dir
region     <- toupper(opt$its_region)
if (!region %in% c("ITS1","ITS2")) region <- "ITS2"

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
log_file <- file.path(output_dir, "pipeline.log")

prog <- function(pct, label) {
  cat(sprintf("PROGRESS:%d|%s\n", as.integer(pct), label))
  flush.console()
}

cat("=== ITS Pipeline Starting ===\n")
cat("Region     :", region, "\n")
cat("Input dir  :", input_dir, "\n")
cat("Output dir :", output_dir, "\n\n")

# ── R Library Paths ───────────────────────────────────────────
cat("=== R Library Paths ===\n")
for (lp in .libPaths()) cat(lp, "\n")
cat("\n")

# ── Find FASTQ files ──────────────────────────────────────────
prog(2, "Step 1/6 — Finding FASTQ files")

# Accept both .fastq/.fastq.gz and .fq/.fq.gz, and both _R1/_R2 and _1/_2 naming
fq_all  <- list.files(input_dir,
                      pattern="\\.(fastq|fq)(\\.gz)?$",
                      full.names=TRUE, recursive=FALSE)
fnFs_raw <- sort(fq_all[grepl("_R1|_1\\.(fq|fastq)", fq_all)])
fnRs_raw <- sort(fq_all[grepl("_R2|_2\\.(fq|fastq)", fq_all)])

cat("Files found:", length(fq_all), "total,",
    length(fnFs_raw), "R1,", length(fnRs_raw), "R2\n")

if (length(fnFs_raw) == 0) stop("No R1 FASTQ files found in: ", input_dir)
if (length(fnRs_raw) == 0) stop("No R2 FASTQ files found in: ", input_dir)
if (length(fnFs_raw) != length(fnRs_raw)) stop("Unequal number of R1 and R2 files")

sample_names <- sub("_R1.*","", basename(fnFs_raw))
cat("Samples found:", length(sample_names), "\n")
cat(paste(sample_names, collapse=", "), "\n\n")

# ── Step 1a: Primer trimming with cutadapt ────────────────────
prog(5, "Step 1/6 — Trimming primers with cutadapt")
cat("Step 1/6a: Primer Trimming (cutadapt)...\n")

# Default primers by region if not specified
default_primer_f <- if (region == "ITS1") "CTTGGTCATTTAGAGGAAGTAA" else "GTGAATCATCGAATCTTTGAA"
default_primer_r <- if (region == "ITS1") "GCTGCGTTCTTCATCGATGC"   else "TCCTCCGCTTATTGATATGC"

primer_f <- if (nchar(opt$primer_f) > 0) opt$primer_f else default_primer_f
primer_r <- if (nchar(opt$primer_r) > 0) opt$primer_r else default_primer_r

cat("Primer F:", primer_f, "\n")
cat("Primer R:", primer_r, "\n")

# Reverse complement helper
RC <- function(seq) {
  chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(seq, "")[[1]]), collapse=""))
}
primer_f_rc <- RC(primer_f)
primer_r_rc <- RC(primer_r)

# system() only inherits a minimal PATH that may not include ~/.local/bin,
# where `pip install --user cutadapt` puts the binary — without this fallback,
# cutadapt_ok silently comes back FALSE and trimming is skipped with no clear
# error even though cutadapt is actually installed.
resolve_cutadapt <- function() {
  candidates <- c("cutadapt", path.expand("~/.local/bin/cutadapt"), "/usr/local/bin/cutadapt")
  for (c in candidates) {
    ok <- suppressWarnings(system(paste(shQuote(c), "--version"),
                                  ignore.stdout=TRUE, ignore.stderr=TRUE) == 0)
    if (ok) return(c)
  }
  NULL
}
CUTADAPT_BIN <- resolve_cutadapt()
cutadapt_ok <- !is.null(CUTADAPT_BIN)
trim_dir <- file.path(output_dir, "trimmed")

if (cutadapt_ok) {
  cat("Using cutadapt at:", CUTADAPT_BIN, "\n")
  dir.create(trim_dir, recursive=TRUE, showWarnings=FALSE)
  trimFs <- file.path(trim_dir, paste0(sample_names, "_F_trim.fastq.gz"))
  trimRs <- file.path(trim_dir, paste0(sample_names, "_R_trim.fastq.gz"))

  discard_flag <- if (isTRUE(opt$discard_untrimmed)) "--discard-untrimmed" else ""
  n_trimmed <- 0
  for (i in seq_along(fnFs_raw)) {
    cmd <- sprintf(
      "%s -g %s -a %s -G %s -A %s -e %s -O %d -m 50 -j %d %s -o %s -p %s %s %s > /dev/null 2>&1",
      CUTADAPT_BIN,
      primer_f, primer_r_rc,
      primer_r, primer_f_rc,
      opt$error_rate, opt$min_overlap,
      opt$threads, discard_flag,
      trimFs[i], trimRs[i], fnFs_raw[i], fnRs_raw[i]
    )
    ret <- system(cmd)
    if (ret == 0 && file.exists(trimFs[i]) && file.size(trimFs[i]) > 0) n_trimmed <- n_trimmed + 1
  }
  cat(sprintf("cutadapt: %d/%d samples trimmed OK\n", n_trimmed, length(fnFs_raw)))

  # Use trimmed files if we got any; else fall back to raw
  valid <- file.exists(trimFs) & file.size(trimFs) > 0
  if (sum(valid) > 0) {
    fnFs_raw <- trimFs[valid]
    fnRs_raw <- trimRs[valid]
    sample_names <- sample_names[valid]
    cat("Using cutadapt-trimmed reads for downstream steps\n\n")
  } else {
    cat("WARNING: cutadapt produced no output. Using untrimmed reads.\n\n")
  }
} else {
  cat("WARNING: cutadapt not found — primers NOT trimmed.\n")
  cat("Install: pip install cutadapt  (or re-run setup.sh)\n\n")
}

# ── Step 1b: Filter (NO truncLen for ITS!) ────────────────────
prog(10, "Step 1/6 — Quality filtering (no length truncation)")
cat("Step 1/6b: Quality Filtering (ITS — no truncLen)...\n")

filt_dir <- file.path(output_dir, "filtered")
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# ITS: NO truncLen; use minLen to remove tiny fragments
out <- filterAndTrim(
  fnFs_raw, filtFs, fnRs_raw, filtRs,
  maxN=0, maxEE=c(opt$maxEE_f, opt$maxEE_r),
  truncQ=2, rm.phix=TRUE,
  minLen=opt$minLen,        # key ITS setting — remove fragments < 50 bp
  compress=TRUE, multithread=opt$threads
)
cat("Filter results:\n")
print(out)

# Remove samples with 0 reads
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
sample_names <- names(filtFs)
if (length(sample_names) == 0) stop("All samples had 0 reads after filtering")
cat("\nSamples passing filter:", length(sample_names), "\n\n")

# ── Step 2: Learn Error Rates ─────────────────────────────────
prog(25, "Step 2/6 — Learning error rates")
cat("Step 2/6: Learning Error Rates...\n")

# WORKAROUND: dada2's learnErrors(nbases=...) does not reliably honor the
# nbases cutoff — a well-documented, long-standing upstream bug (dada2 GitHub
# issues #954 and #2054). Confirmed on our own 16S/12S data: changing nbases
# had zero effect on how much was actually read. Since we can't fix dada2
# itself, pre-select our own subset of files so it physically can't read more
# than intended, regardless of whether its internal cutoff logic works.
select_files_for_nbases <- function(files, reads_out, read_len, nbases) {
  # `files` and `reads_out` must already be positionally aligned. Apply one
  # boolean mask to both together — never filter+reindex them separately, or
  # a zero-read sample in the middle of the list will silently misalign
  # every entry after it.
  keep      <- file.exists(files)
  files     <- files[keep]
  reads_out <- reads_out[keep]
  if (length(files) == 0) return(files)
  reads_out[is.na(reads_out)] <- 0
  cum <- cumsum(pmax(reads_out, 0) * max(read_len, 1))
  n_needed <- which(cum >= nbases)[1]
  if (is.na(n_needed)) n_needed <- length(files)  # target never reached — use everything available
  files[seq_len(n_needed)]
}
# ITS has no fixed truncLen (variable-length reads); use minLen as a conservative
# per-read base estimate (reads are guaranteed to be >= minLen, so this slightly
# over-selects files rather than under-selecting — the safe direction).
# NOTE: this must run against the *pre-filtering* filtFs/filtRs (before the
# file.exists() subsetting a few lines above reassigned them), so pair it with
# the original `out` matrix — but since filtFs/filtRs were already
# reassigned to the filtered subset above, and `out` still has one row per
# *original* input file, only rely on this when the shapes line up.
its_readLen_est <- max(opt$minLen, 200)
its_reads_out <- if (!is.null(out) && ncol(out) >= 2 && nrow(out) == length(filtFs)) {
  out[, 2]
} else {
  rep(NA, length(filtFs))  # shape mismatch — fall back to "use everything" below
}
filtFs_for_err <- select_files_for_nbases(filtFs, its_reads_out, its_readLen_est, as.numeric(opt$nbases))
filtRs_for_err <- select_files_for_nbases(filtRs, its_reads_out, its_readLen_est, as.numeric(opt$nbases))
cat(sprintf("  Using %d/%d sample file(s) to reach the ~%.0f bases target for error learning\n",
            length(filtFs_for_err), length(filtFs), as.numeric(opt$nbases)))

errF <- learnErrors(filtFs_for_err, nbases=as.numeric(opt$nbases), multithread=opt$threads, verbose=FALSE)
errR <- learnErrors(filtRs_for_err, nbases=as.numeric(opt$nbases), multithread=opt$threads, verbose=FALSE)
cat("  Error rates learned.\n\n")

# ── Step 3: Dereplicate ───────────────────────────────────────
prog(40, "Step 3/6 — Dereplicating reads")
cat("Step 3/6: Dereplicating...\n")
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# ── Step 4: DADA2 Denoising ───────────────────────────────────
prog(50, "Step 4/6 — DADA2 denoising")
cat("Step 4/6: DADA2 Denoising...\n")

pool_setting <- switch(opt$pool,
  "TRUE"  = TRUE,
  "FALSE" = FALSE,
  "pseudo" = "pseudo",
  FALSE)

dadaFs <- dada(derepFs, err=errF, pool=pool_setting, multithread=opt$threads)
dadaRs <- dada(derepRs, err=errR, pool=pool_setting, multithread=opt$threads)
if (length(sample_names)==1) {
  dadaFs <- list(dadaFs); names(dadaFs) <- sample_names
  dadaRs <- list(dadaRs); names(dadaRs) <- sample_names
}
cat("  Denoising complete.\n\n")

# ── Step 5: Merge + Chimera removal ──────────────────────────
prog(60, "Step 5/6 — Merging and chimera removal")
cat("Step 5/6: Merging paired reads...\n")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      minOverlap=12, verbose=FALSE)  # lower minOverlap for ITS
if (is.data.frame(mergers)) mergers <- setNames(list(mergers), sample_names)

seqtab       <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=opt$threads)
cat("  ASVs before chimera removal:", ncol(seqtab), "\n")
cat("  ASVs after chimera removal:",  ncol(seqtab_nochim), "\n\n")

# ── Read tracking ─────────────────────────────────────────────
track <- cbind(
  out[sample_names,, drop=FALSE],
  sapply(dadaFs, function(d) sum(d$denoised)),
  sapply(dadaRs, function(d) sum(d$denoised)),
  sapply(mergers, function(m) sum(m$accept)),
  rowSums(seqtab_nochim)
)
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample_names
write.csv(track, file.path(output_dir, "read_tracking.csv"))

# ── Save ASV table ────────────────────────────────────────────
asv_df <- data.frame(sequence=colnames(seqtab_nochim), t(seqtab_nochim),
                     check.names=FALSE, stringsAsFactors=FALSE)
write.csv(asv_df, file.path(output_dir, "asv_table.csv"), row.names=FALSE)

# ── Step 6: Taxonomy ──────────────────────────────────────────
prog(70, "Step 6/6 — Taxonomy (UNITE)")
cat("Step 6/6: Taxonomic Assignment (UNITE)...\n")

# Find UNITE database
# Prefer --db_paths JSON (passed by main.py), fall back to --dbPath or auto-detect
unite_db <- if (nchar(opt$dbPath) > 0) opt$dbPath else ""

if (nchar(unite_db)==0 || !file.exists(unite_db)) {
  # Try db_paths JSON — accept both --db_paths (new) and auto-location (old)
  db_paths_candidates <- c(
    opt$db_paths,
    file.path(dirname(dirname(input_dir)), "databases", "db_paths.json")
  )
  for (db_paths_file in db_paths_candidates) {
    if (nchar(db_paths_file) > 0 && file.exists(db_paths_file)) {
      db_paths <- fromJSON(db_paths_file)
      key <- if (region=="ITS1") "UNITE_ITS1" else "UNITE_ITS2"
      candidate <- db_paths[[key]]
      if (!is.null(candidate) && nchar(candidate) > 0 && file.exists(candidate)) {
        unite_db <- candidate
        break
      }
    }
  }
}
if (nchar(unite_db)==0 || !file.exists(unite_db)) {
  # Search databases/ folder
  db_dir <- file.path(path.expand("~"), "r16s-app", "backend", "databases", "UNITE")
  candidates <- c(
    list.files(db_dir, pattern="sh_general.*\\.fasta(\\.gz)?$", full.names=TRUE),
    list.files(db_dir, pattern="UNITE.*\\.fasta(\\.gz)?$", full.names=TRUE)
  )
  if (length(candidates) > 0) unite_db <- candidates[1]
}

tax <- NULL
if (nchar(unite_db)>0 && file.exists(unite_db)) {
  cat("  Using UNITE database:", unite_db, "\n")
  tryCatch({
    tax <- assignTaxonomy(
      seqtab_nochim, unite_db,
      taxLevels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
      multithread=opt$threads, minBoot=opt$minBoot, tryRC=TRUE
    )
    cat("  Taxonomy assigned to", sum(!is.na(tax[,"Genus"])), "ASVs at genus level\n")
  }, error=function(e) cat("  Taxonomy error:", e$message, "\n"))
} else {
  cat("  WARNING: UNITE database not found. Run download_databases.sh first.\n")
  cat("  Skipping taxonomy...\n")
}

# ── Save taxonomy table ───────────────────────────────────────
if (!is.null(tax)) {
  tax_df <- as.data.frame(tax)
  tax_df$sequence <- rownames(tax)
  write.csv(tax_df, file.path(output_dir, "taxonomy_table.csv"))
}

# ── QC plot ───────────────────────────────────────────────────
has_ggplot2 <- requireNamespace("ggplot2", quietly=TRUE)
has_reshape2 <- requireNamespace("reshape2", quietly=TRUE)

prog(82, "Generating QC and taxonomy plots")

# Read tracking plot
tryCatch({
  track_df <- as.data.frame(track)
  track_df$sample <- rownames(track_df)
  if (has_ggplot2 && has_reshape2) {
    long <- reshape2::melt(track_df, id.vars="sample",
                           variable.name="Step", value.name="Reads")
    long$Step <- factor(long$Step, levels=colnames(track))
    p <- ggplot2::ggplot(long, ggplot2::aes(x=Step, y=Reads,
                         color=sample, group=sample)) +
      ggplot2::geom_line(size=1) + ggplot2::geom_point(size=2) +
      ggplot2::scale_y_continuous(labels=scales::comma) +
      ggplot2::labs(title="Read Counts Through Pipeline Steps",
                    x="Pipeline Step", y="Number of Reads") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
    ggplot2::ggsave(file.path(output_dir,"read_tracking_plot.pdf"),
                    p, width=10, height=6)
  }
}, error=function(e) cat("  read tracking plot error:", e$message, "\n"))

# ── Taxonomy plots ────────────────────────────────────────────
if (!is.null(tax)) {
  palette30 <- c("#e74c3c","#e67e22","#f1c40f","#2ecc71","#1abc9c",
                 "#3498db","#9b59b6","#e91e63","#00bcd4","#8bc34a",
                 "#ff5722","#607d8b","#795548","#ffc107","#03a9f4",
                 "#4caf50","#673ab7","#ff9800","#009688","#9e9e9e",
                 "#f44336","#2196f3","#4caf50","#ff9800","#9c27b0",
                 "#00bcd4","#ffeb3b","#795548","#607d8b","#e91e63")

  make_pal <- function(n) {
    if (n <= length(palette30)) palette30[1:n]
    else rep_len(palette30, n)
  }

  make_tax_mat <- function(level, top_n=opt$topN) {
    col_idx <- which(colnames(tax) == level)
    if (length(col_idx)==0) return(NULL)
    taxa_vec <- tax[,col_idx]
    taxa_vec[is.na(taxa_vec)|taxa_vec==""] <- "Unclassified"
    taxa_vec <- sub("^[a-z]__","",taxa_vec)
    taxa_uniq <- unique(taxa_vec)
    mat <- sapply(taxa_uniq, function(t)
      rowSums(seqtab_nochim[,taxa_vec==t,drop=FALSE]))
    if (is.null(dim(mat)))
      mat <- matrix(mat,nrow=1,dimnames=list(sample_names,taxa_uniq))
    mat <- mat[,order(colSums(mat),decreasing=TRUE),drop=FALSE]
    if (ncol(mat)>top_n) {
      other <- rowSums(mat[,(top_n+1):ncol(mat),drop=FALSE])
      mat <- cbind(mat[,1:top_n,drop=FALSE],Other=other)
    }
    pct <- sweep(mat,1,rowSums(mat),"/")*100
    pct[is.nan(pct)] <- 0
    pct
  }

  for (lvl in c("Phylum","Class","Order","Family","Genus")) {
    tryCatch({
      m <- make_tax_mat(lvl)
      if (is.null(m)) next
      cols <- make_pal(ncol(m))
      names(cols) <- colnames(m)
      out_f <- file.path(output_dir, paste0("taxonomy_",tolower(lvl),".pdf"))
      if (has_ggplot2 && has_reshape2) {
        df <- as.data.frame(m)
        df$Sample <- factor(rownames(m),levels=rownames(m))
        long2 <- reshape2::melt(df,id.vars="Sample",
                                variable.name="Taxon",value.name="Pct")
        long2$Taxon <- factor(long2$Taxon,levels=colnames(m))
        p2 <- ggplot2::ggplot(long2,ggplot2::aes(x=Sample,y=Pct,fill=Taxon))+
          ggplot2::geom_bar(stat="identity")+
          ggplot2::scale_fill_manual(values=cols)+
          ggplot2::labs(title=paste("Relative Abundance — Top",opt$topN,lvl),
                        y="Relative Abundance (%)")+
          ggplot2::theme_bw(base_size=11)+
          ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),
                         legend.text=ggplot2::element_text(size=7),
                         legend.key.size=ggplot2::unit(0.35,"cm"))
        w <- max(7, length(sample_names)*1.2+4)
        ggplot2::ggsave(out_f,p2,width=w,height=7,device="pdf",limitsize=FALSE)
      }
      cat("  Generated:", basename(out_f), "\n")
    }, error=function(e) cat("  taxonomy", lvl, "error:", e$message, "\n"))
  }
}

# ── FUNGuildR — Ecological guild annotation ───────────────────
prog(90, "FUNGuildR — ecological guild annotation")
cat("\nFUNGuildR: Assigning ecological guilds...\n")

if (!is.null(tax) && requireNamespace("FUNGuildR", quietly=TRUE)) {
  tryCatch({
    tax_df2 <- as.data.frame(tax, stringsAsFactors=FALSE)
    # Build semicolon-delimited taxonomy string required by FUNGuildR
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    ranks_present <- ranks[ranks %in% colnames(tax_df2)]
    tax_df2$Taxonomy <- apply(tax_df2[,ranks_present,drop=FALSE], 1,
                               function(x) paste(na.omit(x), collapse=";"))
    guild_db <- FUNGuildR::fung_db()
    guilds <- FUNGuildR::funguild_assign(tax_df2, db=guild_db)
    write.csv(guilds, file.path(output_dir, "funguild_guilds.csv"))
    n_assigned <- sum(!is.na(guilds$guild) & guilds$guild != "", na.rm=TRUE)
    cat(sprintf("  Guild annotation: %d / %d ASVs assigned\n",
                n_assigned, nrow(guilds)))

    # Summary table
    if (n_assigned > 0) {
      guild_summary <- table(guilds$trophicMode[!is.na(guilds$trophicMode)])
      guild_df <- as.data.frame(guild_summary)
      colnames(guild_df) <- c("TrophicMode","Count")
      write.csv(guild_df, file.path(output_dir,"funguild_summary.csv"),row.names=FALSE)
    }
  }, error=function(e) cat("  FUNGuildR error:", e$message, "\n"))
} else if (is.null(tax)) {
  cat("  Skipped (no taxonomy)\n")
} else {
  cat("  FUNGuildR not installed — run install_extensions.sh\n")
}

# ── Alpha diversity ───────────────────────────────────────────
prog(93, "Alpha diversity plots")
tryCatch({
  if (requireNamespace("vegan", quietly=TRUE) && has_ggplot2) {
    shannon <- vegan::diversity(seqtab_nochim, index="shannon")
    simpson <- vegan::diversity(seqtab_nochim, index="simpson")
    richness <- rowSums(seqtab_nochim > 0)
    alpha_df <- data.frame(
      Sample=sample_names,
      Shannon=shannon, Simpson=simpson, Observed=richness
    )
    write.csv(alpha_df, file.path(output_dir,"alpha_diversity.csv"),row.names=FALSE)

    long3 <- reshape2::melt(alpha_df,id.vars="Sample",
                             variable.name="Metric",value.name="Value")
    p3 <- ggplot2::ggplot(long3,ggplot2::aes(x=Sample,y=Value,fill=Sample))+
      ggplot2::geom_bar(stat="identity",show.legend=FALSE)+
      ggplot2::facet_wrap(~Metric,scales="free_y")+
      ggplot2::labs(title="Alpha Diversity (ITS)",x="Sample")+
      ggplot2::theme_bw(base_size=11)+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
    ggplot2::ggsave(file.path(output_dir,"alpha_diversity.pdf"),
                    p3,width=10,height=5,device="pdf")
    cat("  Alpha diversity: OK\n")
  }
}, error=function(e) cat("  alpha diversity error:", e$message, "\n"))

# ── Write ASV sequences as FASTA ──────────────────────────────
fasta_out <- file.path(output_dir, "asv_sequences.fasta")
seqs <- colnames(seqtab_nochim)
asv_ids <- paste0("ASV", seq_along(seqs))
writeLines(c(rbind(paste0(">",asv_ids), seqs)), fasta_out)
cat("  ASV sequences saved:", fasta_out, "\n")

prog(100, "ITS pipeline complete!")
cat("\n=== ITS Pipeline Complete ===\n")
cat("ASVs:", ncol(seqtab_nochim), "\n")
cat("Samples:", nrow(seqtab_nochim), "\n")
cat("Output:", output_dir, "\n\n")
