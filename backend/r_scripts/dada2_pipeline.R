# ============================================================
#  DADA2 Pipeline — 16S / 12S Amplicon Analysis
# ============================================================

# ── Fix R library path (works for any user / machine) ────────
tryCatch({
  # 1. Add R_LIBS_USER if set (handles non-sudo installs)
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset="")
  if (nchar(r_libs_user) > 0 && dir.exists(r_libs_user))
    .libPaths(unique(c(r_libs_user, .libPaths())))

  # 2. Scan ~/R/ for any versioned library directories
  home_dir <- Sys.getenv("HOME", unset=path.expand("~"))
  r_root   <- file.path(home_dir, "R")
  if (dir.exists(r_root)) {
    all_ver_dirs <- list.dirs(r_root, recursive=TRUE, full.names=TRUE)
    all_ver_dirs <- all_ver_dirs[grepl("/\\d+\\.\\d+$", all_ver_dirs)]
    if (length(all_ver_dirs) > 0)
      .libPaths(unique(c(all_ver_dirs, .libPaths())))
  }

  # 3. Also add ~/R/library as a fallback (simple user lib)
  simple_lib <- file.path(home_dir, "R", "library")
  if (dir.exists(simple_lib))
    .libPaths(unique(c(simple_lib, .libPaths())))

}, error=function(e) NULL)

cat("=== R Library Paths ===\n")
cat(paste(.libPaths(), collapse = "\n"), "\n\n")

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(jsonlite)
})

# ── Load optional visualization packages ───────────────────────
# NOTE: moved up here (was previously declared much later, around the
# "Generating plots" step) because earlier code — the error-model PDF plots
# right after learnErrors() — also references has_ggplot2, but ran BEFORE
# this used to be defined, causing a silent "object 'has_ggplot2' not found"
# error (caught by tryCatch, so it just skipped the plots with no real
# warning). Defining it here makes it available to every section that needs
# it, in the correct order.
has_ggplot2  <- requireNamespace("ggplot2",  quietly=TRUE)
has_reshape2 <- requireNamespace("reshape2", quietly=TRUE)
has_vegan    <- requireNamespace("vegan",    quietly=TRUE)
has_ape      <- requireNamespace("ape",      quietly=TRUE)
has_pheatmap <- requireNamespace("pheatmap", quietly=TRUE)
has_scales   <- requireNamespace("scales",   quietly=TRUE)
if (has_ggplot2)  suppressPackageStartupMessages(library(ggplot2))
if (has_reshape2) suppressPackageStartupMessages(library(reshape2))
if (has_vegan)    suppressPackageStartupMessages(library(vegan))
if (has_ape)      suppressPackageStartupMessages(library(ape))
if (has_pheatmap) suppressPackageStartupMessages(library(pheatmap))

# ── Arguments ────────────────────────────────────────────────
option_list <- list(
  make_option("--input",         type="character"),
  make_option("--output",        type="character"),
  make_option("--marker",        type="character", default="16S"),
  make_option("--truncLen_F",    type="integer",   default=240),
  make_option("--truncLen_R",    type="integer",   default=200),
  make_option("--maxEE_F",       type="double",    default=2),
  make_option("--maxEE_R",       type="double",    default=2),
  make_option("--trimLeft_F",    type="integer",   default=0),
  make_option("--trimLeft_R",    type="integer",   default=0),
  make_option("--nbases",        type="integer",   default=100000000),
  make_option("--pool",          type="character", default="FALSE"),
  make_option("--chimeraMethod", type="character", default="consensus"),
  make_option("--taxDatabase",   type="character", default="SILVA_138"),
  make_option("--dbPath",        type="character", default=""),
  make_option("--minBoot",       type="integer",   default=50),
  make_option("--topN",          type="integer",   default=30),
  make_option("--metadata",      type="character", default=""),
  make_option("--tax4fun",       type="logical",   default=FALSE,
              help="Run Tax4Fun2 functional prediction (16S only)"),
  make_option("--tax4fun_ref",   type="character", default="",
              help="Path to Tax4Fun2_ReferenceData_v2 directory"),
  make_option("--nema_db",       type="character", default="",
              help="Path to 18S-NemaBase FASTA for nematode 18S mode"),
  make_option("--db_paths_json", type="character", default="",
              help="Path to db_paths.json for auto-detecting databases"),
  make_option("--primer_f",          type="character", default="",
              help="Forward primer sequence for cutadapt trimming (blank = skip)"),
  make_option("--primer_r",          type="character", default="",
              help="Reverse primer sequence for cutadapt trimming (blank = skip)"),
  make_option("--error_rate",        type="double",    default=0.1,
              help="cutadapt -e: fraction of mismatches allowed (default 0.1)"),
  make_option("--min_overlap",       type="integer",   default=3,
              help="cutadapt -O: minimum overlap bases (default 3)"),
  make_option("--discard_untrimmed", type="logical",   default=FALSE,
              help="Discard reads where primer was not found"),
  make_option("--single_end",        type="logical",   default=FALSE,
              help="Single-end mode: ONT R10.4+ reads (no R2, no merge step)"),
  make_option("--ont_minLen",        type="integer",   default=0,
              help="ONT mode: minimum read length to keep after filtering (bp)"),
  make_option("--ont_maxLen",        type="integer",   default=0,
              help="ONT mode: maximum read length to keep after filtering (bp)"),
  make_option("--threads",           type="integer",   default=4,
              help="Number of CPU threads for DADA2 steps (default 4; TRUE=all cores)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ── Thread count ──────────────────────────────────────────────
THREADS <- if (!is.null(opt$threads) && opt$threads > 0) opt$threads else 4

# ── Single-end flag (ONT mode) ────────────────────────────────
is_single <- isTRUE(opt$single_end)

cat("=== DADA2 Pipeline Starting ===\n")
cat("Marker   :", opt$marker, "\n")
cat("Database :", opt$taxDatabase, "\n")
cat("Mode     :", if (is_single) "single-end (ONT)" else "paired-end (Illumina)", "(initial — may be overridden by manifest)\n\n")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)

# Helper: write progress for the frontend to read
prog <- function(pct, label) {
  cat(sprintf("PROGRESS:%d|%s\n", pct, label))
  flush.console()
}

prog(2, "Pipeline initialized — loading libraries")

# ── Step 1: Find FASTQ files ─────────────────────────────────
findFASTQ <- function(dir, pats) {
  for (p in pats) {
    f <- sort(list.files(dir, pattern=p, full.names=TRUE))
    if (length(f) > 0) return(f)
  }
  return(character(0))
}

manifest_file <- file.path(opt$input, "sample_manifest.json")
use_manifest  <- file.exists(manifest_file)

if (use_manifest) {
  # ── Manual sample↔file assignment (from frontend pairing UI) ─────────────
  # manifest.json: [{ "sample": "AC", "file1": "AC_1.fastq.gz", "file2": "AC_2.fastq.gz" }, ...]
  # file2 == "" (or absent) for that row → treated as single-end for that row.
  # Overrides --single_end / filename auto-detection entirely.
  cat("── Using manual sample/file manifest (sample_manifest.json) ──\n")
  manifest <- fromJSON(manifest_file, simplifyDataFrame = TRUE)
  if (!is.data.frame(manifest)) manifest <- as.data.frame(manifest, stringsAsFactors = FALSE)
  if (nrow(manifest) == 0) stop("sample_manifest.json is empty")

  sample_names <- as.character(manifest$sample)
  # Sanitize: sample names become filenames (e.g. <sample>_F_filt.fastq.gz),
  # so strip anything that could be interpreted as a path separator.
  sample_names_raw <- sample_names
  sample_names <- gsub("[\\\\/]+", "_", trimws(sample_names))
  changed <- sample_names_raw != sample_names
  if (any(changed))
    cat("  [WARN] sanitized sample name(s) with path separators: ",
        paste(sample_names_raw[changed], "->", sample_names[changed], collapse=", "), "\n")
  file1        <- as.character(manifest$file1)
  file2        <- if ("file2" %in% colnames(manifest)) as.character(manifest$file2) else rep("", length(file1))
  file2[is.na(file2)] <- ""

  dup_names <- unique(sample_names[duplicated(sample_names)])
  if (length(dup_names) > 0)
    stop("Duplicate sample names in manifest: ", paste(dup_names, collapse=", "))

  fnFs <- file.path(opt$input, file1)
  missing_f <- !file.exists(fnFs)
  if (any(missing_f))
    stop("Manifest file1 not found: ", paste(file1[missing_f], collapse=", "))

  has_pairs <- any(file2 != "")
  if (has_pairs) {
    if (!all(file2 != ""))
      stop("Manifest has a mix of paired and unpaired samples — either give every ",
           "sample a file2, or none. Samples missing file2: ",
           paste(sample_names[file2 == ""], collapse=", "))
    fnRs <- file.path(opt$input, file2)
    missing_r <- !file.exists(fnRs)
    if (any(missing_r))
      stop("Manifest file2 not found: ", paste(file2[missing_r], collapse=", "))
    is_single <- FALSE
    cat(sprintf("  %d sample(s), paired-end merge (manual manifest)\n\n", length(sample_names)))
  } else {
    fnRs <- character(0)
    is_single <- TRUE
    cat(sprintf("  %d sample(s), single-end / long-read (manual manifest)\n\n", length(sample_names)))
  }
} else if (is_single) {
  # ── Single-end (ONT): accept any .fastq / .fastq.gz file ───
  fnFs <- sort(list.files(opt$input,
    pattern="\\.f(q|astq)(\\.gz)?$", full.names=TRUE))
  if (length(fnFs) == 0) stop("No FASTQ files found in: ", opt$input)
  fnRs <- character(0)
  # Extract sample names — handles ONT Guppy/Dorado barcode naming:
  #   FBE77026_pass_barcode08_877a528d_...fastq.gz → barcode08
  #   barcode08.fastq.gz                           → barcode08
  raw_names <- sub("\\.fastq\\.gz$|\\.fastq$|\\.fq\\.gz$|\\.fq$", "",
                   basename(fnFs))
  sample_names <- gsub(".*_pass_(barcode\\d+)_.*", "\\1", raw_names)
  # Fallback: strip trailing _digits if still long
  still_long <- !grepl("^barcode", sample_names)
  sample_names[still_long] <- sub("_[0-9]+$", "", raw_names[still_long])
  cat("Found (single-end):", length(fnFs), "files\n")
  cat("Sample names:", paste(sample_names, collapse=", "), "\n\n")

  # ── Guard: duplicate sample names usually mean this is actually
  #    paired-end data (R1/R2 pairs) mistakenly run in ONT/single-end mode ──
  dup_names <- unique(sample_names[duplicated(sample_names)])
  if (length(dup_names) > 0) {
    cat("ERROR: Duplicate sample names after parsing:", paste(dup_names, collapse=", "), "\n")
    cat("This usually means the uploaded files are PAIRED-END (R1/R2 pairs,\n")
    cat("e.g. SampleA_1.fastq.gz + SampleA_2.fastq.gz) but the job was run in\n")
    cat("ONT / single-end mode, which strips the trailing _1/_2 as if it were\n")
    cat("noise and collapses both mates into one sample name.\n\n")
    cat("Fix: use the manual sample/file pairing option on the upload screen\n")
    cat("to assign File 1 + File 2 per sample and merge as paired-end, or\n")
    cat("create a new job WITHOUT the ONT/single-end option for this data.\n")
    stop("Duplicate sample names in single-end mode — see message above.")
  }
} else {
  # ── Paired-end (Illumina): require R1/R2 ───────────────────
  patterns <- list(
    R1 = c("_R1.*\\.fastq$", "_R1.*\\.fastq\\.gz$", "_1\\.f(q|astq)(\\.gz)?$"),
    R2 = c("_R2.*\\.fastq$", "_R2.*\\.fastq\\.gz$", "_2\\.f(q|astq)(\\.gz)?$")
  )
  fnFs <- findFASTQ(opt$input, patterns$R1)
  fnRs <- findFASTQ(opt$input, patterns$R2)
  if (length(fnFs) == 0) stop("No FASTQ files found in: ", opt$input)
  cat("Found R1:", length(fnFs), "files\n")
  cat("Found R2:", length(fnRs), "files\n\n")
  sample_names <- sub("_R1.*|_1\\.(fq|fastq).*", "", basename(fnFs))
}
cat("Mode (resolved):", if (is_single) "single-end" else "paired-end", "\n\n")
prog(8, sprintf("Found %d sample(s) — ready to process", length(fnFs)))

# ── Step 1b (optional): Cutadapt Primer Trimming ─────────────
RC <- function(seq) {
  chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(seq, "")[[1]]), collapse=""))
}

# system() only inherits a minimal PATH that may not include ~/.local/bin,
# where `pip install --user cutadapt` puts the binary. Without checking that
# fallback, cutadapt_ok silently comes back FALSE and primer trimming gets
# skipped with no clear error — even though cutadapt is actually installed
# and works fine for every other pipeline (Python scripts already check this
# fallback via shutil.which).
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

if (nchar(opt$primer_f) > 0 && nchar(opt$primer_r) > 0) {
  prog(10, "Step 1/8 — Cutadapt primer trimming")
  cat("Step 1b: Primer Trimming with cutadapt...\n")
  cat("  Primer F:", opt$primer_f, "\n")
  cat("  Primer R:", opt$primer_r, "\n")
  if (is_single) cat("  Mode: single-end\n")

  cutadapt_ok <- !is.null(CUTADAPT_BIN)
  if (!cutadapt_ok) cat("  [WARN] cutadapt not found on PATH or in ~/.local/bin — skipping primer trimming\n")
  else cat("  Using cutadapt at:", CUTADAPT_BIN, "\n")

  if (cutadapt_ok) {
    primer_f_rc <- RC(opt$primer_f)
    primer_r_rc <- RC(opt$primer_r)
    cut_dir     <- file.path(opt$output, "cutadapt")
    dir.create(cut_dir, recursive=TRUE, showWarnings=FALSE)

    discard_flag <- if (isTRUE(opt$discard_untrimmed)) "--discard-untrimmed" else ""
    n_cut <- 0

    if (is_single) {
      # ── Single-end cutadapt ──────────────────────────────────
      cutFs <- file.path(cut_dir, paste0(sample_names, "_cut.fastq.gz"))
      for (i in seq_along(fnFs)) {
        cmd <- paste(
          CUTADAPT_BIN,
          "-g", opt$primer_f, "-a", primer_r_rc,
          "-e", opt$error_rate,
          "-O", opt$min_overlap,
          sprintf("-m 50 --cores=%d", THREADS),
          discard_flag,
          "-o", cutFs[i],
          fnFs[i],
          "> /dev/null 2>&1"
        )
        ret <- system(cmd)
        if (ret == 0 && file.exists(cutFs[i]) && file.size(cutFs[i]) > 0) n_cut <- n_cut + 1
      }
      cat(sprintf("  cutadapt: %d/%d samples trimmed OK\n", n_cut, length(fnFs)))
      valid_cut <- file.exists(cutFs) & file.size(cutFs) > 0
      if (sum(valid_cut) > 0) {
        fnFs         <- cutFs[valid_cut]
        sample_names <- sample_names[valid_cut]
        cat("  Using cutadapt-trimmed reads for downstream steps\n\n")
      } else {
        cat("  WARNING: cutadapt produced no output — using untrimmed reads\n\n")
      }
    } else {
      # ── Paired-end cutadapt ──────────────────────────────────
      cutFs <- file.path(cut_dir, paste0(sample_names, "_F_cut.fastq.gz"))
      cutRs <- file.path(cut_dir, paste0(sample_names, "_R_cut.fastq.gz"))
      for (i in seq_along(fnFs)) {
        cmd <- paste(
          CUTADAPT_BIN,
          "-g", opt$primer_f,   "-a", primer_r_rc,
          "-G", opt$primer_r,   "-A", primer_f_rc,
          "-e", opt$error_rate,
          "-O", opt$min_overlap,
          sprintf("-m 50 --cores=%d", THREADS),
          discard_flag,
          "-o", cutFs[i], "-p", cutRs[i],
          fnFs[i], fnRs[i],
          "> /dev/null 2>&1"
        )
        ret <- system(cmd)
        if (ret == 0 && file.exists(cutFs[i]) && file.size(cutFs[i]) > 0) n_cut <- n_cut + 1
      }
      cat(sprintf("  cutadapt: %d/%d samples trimmed OK\n", n_cut, length(fnFs)))
      valid_cut <- file.exists(cutFs) & file.size(cutFs) > 0
      if (sum(valid_cut) > 0) {
        fnFs         <- cutFs[valid_cut]
        fnRs         <- cutRs[valid_cut]
        sample_names <- sample_names[valid_cut]
        cat("  Using cutadapt-trimmed reads for downstream steps\n\n")
      } else {
        cat("  WARNING: cutadapt produced no output — using untrimmed reads\n\n")
      }
    }
  } else {
    cat("  WARNING: cutadapt not found — primers NOT trimmed\n")
    cat("  Install: pip install cutadapt  (or run setup.sh)\n\n")
  }
} else {
  cat("Step 1b: No primers specified — skipping cutadapt\n\n")
}

# ── Step 2: Filter & Trim ─────────────────────────────────────
prog(12, "Step 1/8 — Filtering and trimming reads")
cat("Step 1/5: Filter & Trim...\n")
filtFs <- file.path(opt$output, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
names(filtFs) <- sample_names

if (is_single) {
  # ── Single-end filterAndTrim (ONT long-read mode) ─────────────
  # ONT amplicon reads vary in length and carry far higher (mostly
  # indel-driven) error rates than Illumina, so a hard truncLen cutoff
  # combined with a strict Illumina-style maxEE throws away almost all
  # reads. Instead: filter by expected length range (minLen/maxLen) and
  # do NOT truncate to a fixed position. Indel-heavy errors are handled
  # later at the dada() step via BAND_SIZE/HOMOPOLYMER_GAP_PENALTY.
  ont_minLen <- if (!is.null(opt$ont_minLen) && opt$ont_minLen > 0) opt$ont_minLen else 50
  ont_maxLen <- if (!is.null(opt$ont_maxLen) && opt$ont_maxLen > 0) opt$ont_maxLen else Inf
  cat(sprintf("  [ONT mode] length filter: %d-%d bp (no truncLen), maxEE=%.1f\n",
              ont_minLen, ont_maxLen, opt$maxEE_F))
  out <- filterAndTrim(
    fnFs, filtFs,
    minLen = ont_minLen, maxLen = ont_maxLen,
    trimLeft = opt$trimLeft_F,
    maxN=0, maxEE=opt$maxEE_F,
    truncQ=2, rm.phix=FALSE, compress=TRUE, multithread=THREADS
  )
  filtRs <- character(0)
} else {
  # ── Paired-end filterAndTrim ─────────────────────────────────
  filtRs <- file.path(opt$output, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
  names(filtRs) <- sample_names
  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen   = c(opt$truncLen_F, opt$truncLen_R),
    trimLeft   = c(opt$trimLeft_F, opt$trimLeft_R),
    maxN=0, maxEE=c(opt$maxEE_F, opt$maxEE_R),
    truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=THREADS
  )
}
cat("  Done. Reads passing filter:", sum(out[,2]), "\n\n")
prog(28, sprintf("Filter & Trim done — %d reads passed", sum(out[,2])))

# ── Save trimmed FastQ files → output_dir/Trim_seq/ ───────────
tryCatch({
  trim_dir <- file.path(opt$output, "Trim_seq")
  dir.create(trim_dir, showWarnings=FALSE, recursive=TRUE)
  for (i in seq_along(filtFs)) {
    if (file.exists(filtFs[i])) {
      file.copy(filtFs[i],
        file.path(trim_dir, paste0(sample_names[i], "_R1.trim.fastq.gz")),
        overwrite=TRUE)
    }
  }
  if (!is_single) {
    for (i in seq_along(filtRs)) {
      if (file.exists(filtRs[i])) {
        file.copy(filtRs[i],
          file.path(trim_dir, paste0(sample_names[i], "_R2.trim.fastq.gz")),
          overwrite=TRUE)
      }
    }
  }
  cat("  Trimmed FastQ saved to Trim_seq/\n")
}, error=function(e) cat("  [skip] Trim_seq copy:", e$message, "\n"))

# ── Step 3: Learn Error Rates ─────────────────────────────────
prog(32, "Step 2/8 — Learning error rates (this takes a while...)")
cat("Step 2/5: Learning Error Rates...\n")

# WORKAROUND: dada2's learnErrors(nbases=...) does not reliably honor the
# nbases cutoff — this is a well-documented, long-standing upstream bug
# (dada2 GitHub issues #954 and #2054: users report the exact same huge
# base/sample count no matter what nbases is set to, sometimes billions of
# bases from a single file). We confirmed this on our own data: setting
# "Bases to Learn From" to 10M vs 100M made literally zero difference — same
# 26 samples / 115,166,400 bases used both times. Since we can't fix dada2
# itself, we work around it by choosing our own subset of *files* to hand to
# learnErrors() in the first place, so it physically cannot read more than
# intended — regardless of whether its own internal cutoff logic works.
select_files_for_nbases <- function(files, reads_out, read_len, nbases, outlier_mult=20) {
  # `files` and `reads_out` must already be positionally aligned (same order,
  # same length, straight from filterAndTrim's own output). Apply one boolean
  # mask to both together — never filter+reindex them separately, or a
  # zero-read sample sitting in the middle of the list will silently
  # misalign every entry after it.
  keep      <- file.exists(files)
  files     <- files[keep]
  reads_out <- reads_out[keep]
  if (length(files) == 0) return(files)
  reads_out[is.na(reads_out)] <- 0

  # Skip extreme outlier samples for error-LEARNING purposes only — they are
  # still processed normally everywhere else in the pipeline, this only
  # affects which files feed the error MODEL. A single sample that dwarfs
  # every other sample (e.g. one huge library sitting next to many small
  # ones) forces learnErrors() to swallow it whole just to reach *any*
  # nbases target, defeating the point of the setting entirely. dada2's
  # error model represents the *sequencing run's* error characteristics, not
  # any one sample's biology — using a representative subset (excluding
  # wildly oversized outliers) is standard DADA2 practice, not a shortcut.
  pos <- reads_out[reads_out > 0]
  if (length(pos) >= 3) {
    med <- stats::median(pos)
    normal <- reads_out > 0 & reads_out <= outlier_mult * med
    if (any(normal) && !all(normal)) {
      files_n     <- files[normal]
      reads_out_n <- reads_out[normal]
      cum_n <- cumsum(pmax(reads_out_n, 0) * max(read_len, 1))
      # Only take this path if the non-outlier samples alone can reach a
      # reasonable fraction of the target — otherwise they're too sparse to
      # be a fair substitute, so fall through to the full pool below.
      if (length(cum_n) > 0 && max(cum_n) >= nbases * 0.1) {
        n_skip <- sum(!normal)
        cat(sprintf("  [nbases pre-select] skipping %d outlier sample(s) (>%.0fx median reads) — using %d representative sample(s) totaling ~%.0f bases instead\n",
                    n_skip, outlier_mult, length(files_n), max(cum_n)))
        n_needed <- which(cum_n >= nbases)[1]
        if (is.na(n_needed)) n_needed <- length(files_n)
        return(files_n[seq_len(n_needed)])
      }
    }
  }

  cum <- cumsum(pmax(reads_out, 0) * max(read_len, 1))
  n_needed <- which(cum >= nbases)[1]
  if (is.na(n_needed)) n_needed <- length(files)  # target never reached — use everything available
  files[seq_len(n_needed)]
}

reads_out_vec <- if (!is.null(out) && ncol(out) >= 2 && nrow(out) == length(filtFs)) {
  out[, 2]
} else {
  rep(NA, length(filtFs))  # shape mismatch — fall back to "use everything" below
}
# Paired-end: filterAndTrim truncates to an exact length, so truncLen is a precise
# per-read base estimate. Single-end/ONT mode isn't truncated to a fixed length, so
# use the midpoint of the configured min/max length filter as a rough proxy instead.
readLen_F <- if (is_single) {
  ont_max_est <- if (is.finite(ont_maxLen)) ont_maxLen else 2000
  (ont_minLen + ont_max_est) / 2
} else opt$truncLen_F

# Debug visibility: confirm reads_out_vec actually carries real per-sample counts
# (not silently falling back to all-NA/0, which would force n_needed to always
# equal "every file that exists" regardless of the nbases target).
cat(sprintf("  [nbases pre-select] read_len_est=%.0f  reads_out: n=%d sum=%.0f min=%.0f median=%.0f max=%.0f  (%d/%d NA)\n",
            readLen_F, length(reads_out_vec), sum(reads_out_vec, na.rm=TRUE),
            suppressWarnings(min(reads_out_vec, na.rm=TRUE)),
            stats::median(reads_out_vec, na.rm=TRUE),
            suppressWarnings(max(reads_out_vec, na.rm=TRUE)),
            sum(is.na(reads_out_vec)), length(reads_out_vec)))

filtFs_for_err <- select_files_for_nbases(filtFs, reads_out_vec, readLen_F, opt$nbases)
cat(sprintf("  Using %d/%d sample file(s) to reach the ~%.0f bases target for error learning\n",
            length(filtFs_for_err), length(filtFs), opt$nbases))
errF <- learnErrors(filtFs_for_err, nbases=opt$nbases, multithread=THREADS)
if (!is_single) {
  filtRs_for_err <- select_files_for_nbases(filtRs, reads_out_vec, opt$truncLen_R, opt$nbases)
  errR <- learnErrors(filtRs_for_err, nbases=opt$nbases, multithread=THREADS)
}
cat("  Done.\n\n")

# ── Save Error Model plots → output_dir/Dada2/ ────────────────
tryCatch({
  dada2_dir <- file.path(opt$output, "Dada2")
  dir.create(dada2_dir, showWarnings=FALSE, recursive=TRUE)
  if (has_ggplot2) {
    p_errF <- dada2::plotErrors(errF, nominalQ=TRUE)
    ggplot2::ggsave(file.path(dada2_dir, "ErrorModel_R1.pdf"),
                    p_errF, width=10, height=8, device="pdf")
    if (!is_single) {
      p_errR <- dada2::plotErrors(errR, nominalQ=TRUE)
      ggplot2::ggsave(file.path(dada2_dir, "ErrorModel_R2.pdf"),
                      p_errR, width=10, height=8, device="pdf")
      cat("  ErrorModel_R1.pdf + ErrorModel_R2.pdf saved to Dada2/\n")
    } else {
      cat("  ErrorModel_R1.pdf saved to Dada2/\n")
    }
  }
}, error=function(e) cat("  [skip] Error model plots:", e$message, "\n"))

# ── Step 4: Denoise ───────────────────────────────────────────
prog(48, "Step 3/8 — Denoising & ASV inference")
cat("Step 3/5: Denoising (DADA2 Inference)...\n")

# WORKAROUND: unlike error-LEARNING above, every sample here still needs its
# own ASVs — we can't just skip a big one. But derepFastq()/dada() have no
# built-in cap on reads-per-sample, and both scale with read count AND
# unique-sequence count. A single wildly oversized sample (e.g. one library
# sequenced 20-50x deeper than the rest of the batch — we've seen a sample
# in this pipeline with >1M reads / ~790K unique sequences, i.e. >75%
# uniqueness, vs. a typical sample's few thousand uniques) can stall the
# whole run for many hours on that one file alone while every other sample
# finishes in seconds (matches the low, single-core-pegged CPU% observed
# live: 35 samples done, 1 giant one still grinding through derep/dada()).
# Subsampling any such sample down to a generous cap before denoising keeps
# its ASV composition/proportions essentially unchanged (same statistical
# logic as rarefaction — a large-enough random subset represents the same
# underlying sequence distribution) while bounding worst-case runtime. The
# full, uncapped reads are still what's saved to Trim_seq/ for the user —
# only what feeds derepFastq()/dada() below is capped.
# NOTE on the cap value: 1M raw reads is NOT, by itself, a lot of data — the
# real cost driver is dada()'s own partitioning algorithm, whose runtime
# scales worse than linearly with the number of UNIQUE sequences (not raw
# read count). A 75%+ uniqueness ratio (791K uniques from 1.05M reads, vs. a
# few thousand uniques for a normal sample) is itself unusual for real
# amplicon data and is what actually explains a 10+ hour stall — subsampling
# reads reduces the unique-sequence count roughly proportionally (a sample
# this noisy stays noisy at any depth), so the cap is set low enough (100K,
# not just "a bit less than 1M") to reliably keep the resulting unique-count
# in a tractable range rather than merely delaying the same problem.
MAX_READS_PER_SAMPLE_DENOISE <- 100000
cap_sample_for_denoise <- function(pathF, pathR, max_reads) {
  if (!file.exists(pathF) || !requireNamespace("ShortRead", quietly=TRUE))
    return(list(F=pathF, R=pathR, capped=FALSE))
  fqF <- ShortRead::readFastq(pathF)
  n <- length(fqF)
  if (n <= max_reads) return(list(F=pathF, R=pathR, capped=FALSE))
  set.seed(42)
  keep <- sort(sample.int(n, max_reads))
  outF <- sub("\\.fastq\\.gz$", "_capped.fastq.gz", pathF)
  ShortRead::writeFastq(fqF[keep], outF, compress=TRUE, mode="w")
  result <- list(F=outF, R=pathR, capped=TRUE, n_before=n, n_after=max_reads)
  if (!is.null(pathR) && length(pathR) > 0 && !is.na(pathR) && file.exists(pathR)) {
    fqR <- ShortRead::readFastq(pathR)
    if (length(fqR) == n) {
      outR <- sub("\\.fastq\\.gz$", "_capped.fastq.gz", pathR)
      ShortRead::writeFastq(fqR[keep], outR, compress=TRUE, mode="w")
      result$R <- outR
    }
  }
  result
}
if (exists("reads_out_vec") && length(reads_out_vec) == length(filtFs)) {
  for (.i in seq_along(filtFs)) {
    if (!is.na(reads_out_vec[.i]) && reads_out_vec[.i] > MAX_READS_PER_SAMPLE_DENOISE) {
      .pathR <- if (!is_single && length(filtRs) >= .i) filtRs[.i] else NULL
      .capped <- tryCatch(
        cap_sample_for_denoise(filtFs[.i], .pathR, MAX_READS_PER_SAMPLE_DENOISE),
        error=function(e) { cat("  [denoise cap] skipped (error):", e$message, "\n"); list(capped=FALSE) }
      )
      if (isTRUE(.capped$capped)) {
        cat(sprintf("  [denoise cap] %s: %d -> %d reads (capped to bound denoising runtime)\n",
                    sample_names[.i], .capped$n_before, .capped$n_after))
        filtFs[.i] <- .capped$F
        if (!is_single && !is.null(.capped$R)) filtRs[.i] <- .capped$R
      }
    }
  }
}

pool_val <- if (opt$pool == "TRUE") TRUE else if (opt$pool == "FALSE") FALSE else "pseudo"
derepFs  <- derepFastq(filtFs)
if (!is.list(derepFs)) derepFs <- setNames(list(derepFs), sample_names)
names(derepFs) <- sample_names
if (is_single) {
  # ONT reads: BAND_SIZE widens the alignment band and
  # HOMOPOLYMER_GAP_PENALTY=-1 relaxes the indel penalty in homopolymer
  # runs — both needed because ONT errors are predominantly indels,
  # unlike Illumina's substitution-dominated error profile.
  cat("  [ONT mode] dada() using BAND_SIZE=32, HOMOPOLYMER_GAP_PENALTY=-1 (indel-tolerant)\n")
  dadaFs <- dada(derepFs, err=errF, pool=pool_val, multithread=THREADS,
                 BAND_SIZE=32, HOMOPOLYMER_GAP_PENALTY=-1)
} else {
  dadaFs <- dada(derepFs, err=errF, pool=pool_val, multithread=THREADS)
}
if (inherits(dadaFs, "dada")) dadaFs <- setNames(list(dadaFs), sample_names)

if (!is_single) {
  derepRs <- derepFastq(filtRs)
  if (!is.list(derepRs)) derepRs <- setNames(list(derepRs), sample_names)
  names(derepRs) <- sample_names
  dadaRs <- dada(derepRs, err=errR, pool=pool_val, multithread=THREADS)
  if (inherits(dadaRs, "dada")) dadaRs <- setNames(list(dadaRs), sample_names)
}
cat("  Done.\n\n")

# ── Step 5: Sequence Table & Chimera ─────────────────────────
prog(65, "Step 4/8 — Building sequence table & removing chimeras")
cat("Step 4/5: Sequence Table & Chimera Removal...\n")

if (is_single) {
  # Single-end: build sequence table directly from dadaFs (no merge)
  seqtab <- makeSequenceTable(dadaFs)
} else {
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
  if (is.data.frame(mergers)) mergers <- setNames(list(mergers), sample_names)
  seqtab <- makeSequenceTable(mergers)
}
seqtab_nochim <- removeBimeraDenovo(seqtab, method=opt$chimeraMethod, multithread=THREADS)
cat("  ASVs after chimera removal:", ncol(seqtab_nochim), "\n\n")

# ── Read tracking (early) for checkpoint evaluation ───────────
if (is_single) {
  track_early <- cbind(
    out,
    sapply(dadaFs, function(d) sum(d$denoised)),
    rowSums(seqtab_nochim)
  )
  colnames(track_early) <- c("input","filtered","denoised","nonchim")
} else {
  track_early <- cbind(
    out,
    sapply(dadaFs, function(d) sum(d$denoised)),
    sapply(dadaRs, function(d) sum(d$denoised)),
    sapply(mergers, function(m) sum(m$abundance[m$accept])),
    rowSums(seqtab_nochim)
  )
  colnames(track_early) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
}
rownames(track_early) <- sample_names

# ── CHECKPOINT: warn if reads are very low ────────────────────
# Single-end (ONT) mode: skip checkpoint entirely.
# nonchim/input ratio is always low for ONT because input contains
# many variable-quality reads that are expected to be filtered out.
# The relevant check for ONT is nonchim/filtered, not nonchim/input.
if (is_single) {
  nonchim_of_filtered <- mean(track_early[,"nonchim"] /
                              pmax(track_early[,"filtered"], 1) * 100)
  cat(sprintf("  [ONT single-end] nonchim/filtered: %.1f%%\n", nonchim_of_filtered))
  cat(sprintf("  [ONT single-end] nonchim/input:    %.1f%% (expected low — skipping checkpoint)\n",
      mean(track_early[,"nonchim"] / pmax(track_early[,"input"], 1) * 100)))
  merged_pct  <- 100
  nonchim_pct <- 100  # skip checkpoint for single-end ONT
} else {
  merged_pct  <- mean(track_early[,"merged"]  / pmax(track_early[,"input"], 1) * 100)
  nonchim_pct <- mean(track_early[,"nonchim"] / pmax(track_early[,"input"], 1) * 100)
}

if (merged_pct < 10 || nonchim_pct < 10) {
  track_list <- lapply(seq_len(nrow(track_early)), function(i) {
    as.list(track_early[i,])
  })
  names(track_list) <- rownames(track_early)
  checkpoint_data <- list(
    type        = "low_merge_rate",
    merged_pct  = round(merged_pct,  2),
    nonchim_pct = round(nonchim_pct, 2),
    n_samples   = nrow(track_early),
    track       = track_list
  )
  cat(sprintf("CHECKPOINT:low_merge|%s\n",
      toJSON(checkpoint_data, auto_unbox=TRUE)))
  flush.console()

  # Write checkpoint file for reference
  write(toJSON(checkpoint_data, auto_unbox=TRUE),
        file.path(opt$output, "checkpoint.json"))

  # Wait for signal file (up to 2 hours)
  signal_file <- file.path(opt$output, "checkpoint_signal")
  cat("  [checkpoint] Waiting for user decision (merged: ",
      round(merged_pct,1), "%, nonchim: ", round(nonchim_pct,1), "%)...\n", sep="")
  flush.console()
  waited <- 0
  while (!file.exists(signal_file) && waited < 7200) {
    Sys.sleep(2); waited <- waited + 2
  }
  if (!file.exists(signal_file)) {
    stop("Checkpoint timed out — no user response after 2 hours")
  }
  signal <- trimws(readLines(signal_file, warn=FALSE)[1])
  if (signal == "abort") {
    stop("User chose to abort and adjust settings")
  }
  cat("  [checkpoint] User chose to continue. Proceeding...\n")
  flush.console()
}

# ── Step 6: Taxonomy ──────────────────────────────────────────
prog(78, sprintf("Step 5/8 — Taxonomic assignment (%d ASVs)", ncol(seqtab_nochim)))
cat("Step 5/5: Taxonomic Assignment...\n")
tax <- NULL

db_path <- opt$dbPath
cat("  dbPath received:", if (nchar(db_path) > 0) db_path else "(empty)", "\n")

if (db_path == "" || !file.exists(db_path)) {
  cat("  Searching for database in script directory...\n")
  args      <- commandArgs(trailingOnly=FALSE)
  file_arg  <- args[grepl("^--file=", args)]
  script_dir <- if (length(file_arg) > 0)
    dirname(sub("^--file=", "", file_arg[1]))
  else
    dirname(normalizePath(sys.frame(1)$ofile, mustWork=FALSE))
  db_dir <- file.path(dirname(script_dir), "databases")
  cat("  Looking in:", db_dir, "\n")

  db_files <- list.files(db_dir, pattern="\\.fa(\\.gz)?$|\\.fasta(\\.gz)?$",
                         full.names=TRUE, ignore.case=TRUE)
  cat("  Found", length(db_files), "database file(s)\n")

  togenus_files <- db_files[grepl("togenus|toGenus|train_set|trainset", basename(db_files), ignore.case=TRUE) &
                             !grepl("toSpecies|tospecies|assignSpecies", basename(db_files), ignore.case=TRUE)]
  train_files   <- db_files[grepl("train_set|trainset|nr99", basename(db_files), ignore.case=TRUE) &
                             !grepl("toSpecies|tospecies|assignSpecies", basename(db_files), ignore.case=TRUE)]
  other_files   <- db_files[!grepl("species|Species", basename(db_files), ignore.case=TRUE)]

  if (length(togenus_files) > 0) {
    db_path <- togenus_files[1]
    cat("  Selected toGenus trainset:", basename(db_path), "\n")
  } else if (length(train_files) > 0) {
    db_path <- train_files[1]
    cat("  Selected trainset:", basename(db_path), "\n")
  } else if (length(other_files) > 0) {
    db_path <- other_files[1]
    cat("  Selected fallback:", basename(db_path), "\n")
  } else if (length(db_files) > 0) {
    db_path <- db_files[1]
    cat("  Warning: using first available file:", basename(db_path), "\n")
  }
}

cat("  Using database:", basename(db_path), "\n")
is_species_file <- grepl("assignSpecies|toSpecies", basename(db_path), ignore.case=TRUE)

if (!is.null(db_path) && db_path != "" && file.exists(db_path)) {
  if (is_species_file) {
    cat("  Note: '", basename(db_path), "' is an assignSpecies file.\n", sep="")
    cat("  Skipping genus-level taxonomy (need silva_nr99_v138.2_train_set.fa.gz).\n")
    cat("  ASV table will be saved without taxonomy.\n\n")
  } else {
    gc(verbose=FALSE)
    cat("  Memory freed — starting taxonomy assignment...\n")

    mem_free_gb <- tryCatch({
      as.numeric(system("awk '/MemAvailable/ {printf \"%.1f\", $2/1048576}' /proc/meminfo",
                        intern=TRUE))
    }, error=function(e) NA)
    use_threads <- if (!is.na(mem_free_gb) && mem_free_gb < 8) {
      cat(sprintf("  Available RAM: %.1f GB — using single thread to reduce peak memory.\n",
                  mem_free_gb))
      FALSE
    } else {
      cat(sprintf("  Available RAM: %.1f GB — using multi-thread.\n",
                  if (is.na(mem_free_gb)) 99 else mem_free_gb))
      TRUE
    }

    tryCatch({
      tax <- assignTaxonomy(seqtab_nochim, db_path, minBoot=opt$minBoot,
                            multithread=use_threads, verbose=FALSE)
      cat("  Taxonomy assigned successfully.\n\n")
      db_dir2       <- dirname(db_path)
      sp_candidates <- list.files(db_dir2,
                                  pattern="(assignspecies|species_assignment).*\\.fa(\\.gz)?$",
                                  full.names=TRUE, ignore.case=TRUE)
      if (length(sp_candidates) > 0) {
        tryCatch({
          gc(verbose=FALSE)
          tax <- addSpecies(tax, sp_candidates[1])
          cat("  Species added from:", basename(sp_candidates[1]), "\n\n")
        }, error=function(e) cat("  Species assignment skipped:", e$message, "\n\n"))
      }
    }, error=function(e) {
      msg <- e$message
      if (grepl("cannot allocate|out of memory|memory exhausted", msg, ignore.case=TRUE)) {
        cat("  [MEMORY ERROR] Not enough RAM to load the taxonomy database.\n")
        cat("  Try: sudo fallocate -l 16G /swapfile && sudo mkswap /swapfile && sudo swapon /swapfile\n")
        cat("  Error detail:", msg, "\n\n")
      } else {
        cat("  Taxonomy warning:", msg, "\n\n")
      }
    })
  }
} else if (opt$marker %in% c("18S","18S-nema","nematode")) {
  # ── Nematode 18S: auto-detect NemaBase / PR2 database ────────
  cat("  Marker is 18S/nematode — searching NemaBase / PR2...\n")
  nema_db <- opt$nema_db
  if (nchar(nema_db) == 0 || !file.exists(nema_db)) {
    json_path <- opt$db_paths_json
    if (nchar(json_path) == 0)
      json_path <- path.expand("~/r16s-app/backend/databases/db_paths.json")
    if (file.exists(json_path)) {
      db_json <- tryCatch(jsonlite::fromJSON(json_path), error=function(e) list())
      if (!is.null(db_json$NemaBase_18S) && file.exists(db_json$NemaBase_18S))
        nema_db <- db_json$NemaBase_18S
      else if (!is.null(db_json$PR2_18S) && file.exists(db_json$PR2_18S))
        nema_db <- db_json$PR2_18S
    }
  }
  if (nchar(nema_db) > 0 && file.exists(nema_db)) {
    cat("  Using 18S database:", nema_db, "\n")
    tryCatch({
      tax <- assignTaxonomy(seqtab_nochim, nema_db, minBoot=opt$minBoot,
                            multithread=THREADS, verbose=FALSE)
      cat("  18S taxonomy assigned.\n\n")
    }, error=function(e) cat("  18S taxonomy error:", e$message, "\n"))
  } else {
    cat("  No NemaBase/PR2 database found. Run download_databases.sh\n\n")
  }
} else {
  cat("  No database file found — skipping taxonomy assignment.\n\n")
}

# ── Save Results ──────────────────────────────────────────────
prog(88, "Step 6/8 — Saving results & CSV tables...")
cat("Saving results...\n")

# ASV table
asv_df          <- as.data.frame(t(seqtab_nochim))
asv_df$sequence <- rownames(asv_df)
write.csv(asv_df, file.path(opt$output, "asv_table.csv"), row.names=FALSE)

# Read tracking
if (is_single) {
  track <- cbind(
    out,
    sapply(dadaFs, function(d) sum(d$denoised)),
    rowSums(seqtab_nochim)
  )
  colnames(track) <- c("input","filtered","denoised","nonchim")
} else {
  track <- cbind(
    out,
    sapply(dadaFs, function(d) sum(d$denoised)),
    sapply(dadaRs, function(d) sum(d$denoised)),
    sapply(mergers, function(m) sum(m$abundance[m$accept])),
    rowSums(seqtab_nochim)
  )
  colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
}
rownames(track) <- sample_names
write.csv(as.data.frame(track), file.path(opt$output, "read_tracking.csv"))

# Taxonomy table + summary JSON
tax_summary_list <- list()
if (!is.null(tax)) {
  write.csv(as.data.frame(tax), file.path(opt$output, "taxonomy_table.csv"))
  # colSums(seqtab_nochim) = reads per ASV — matches tax rows (ASVs)
  asv_counts <- colSums(seqtab_nochim)
  for (lvl in c("Phylum","Class","Order","Family","Genus")) {
    col_idx <- which(colnames(tax) == lvl)
    if (length(col_idx) > 0) {
      taxon_vec <- tax[, col_idx]
      not_na    <- !is.na(taxon_vec) & taxon_vec != ""
      if (sum(not_na) == 0) next
      tbl   <- sort(tapply(asv_counts[not_na], taxon_vec[not_na], sum), decreasing=TRUE)
      total <- sum(tbl)
      tax_summary_list[[lvl]] <- lapply(names(tbl), function(n) {
        list(name=n, abundance=round(as.numeric(tbl[n])/total*100, 2))
      })
    }
  }
  write(toJSON(tax_summary_list, auto_unbox=TRUE),
        file.path(opt$output, "taxonomy_summary.json"))
}

# Summary JSON
summary_data <- list(
  marker       = opt$marker,
  database     = opt$taxDatabase,
  n_samples    = nrow(seqtab_nochim),
  n_asvs       = ncol(seqtab_nochim),
  total_reads  = sum(seqtab_nochim),
  has_taxonomy = !is.null(tax),
  output_files = list.files(opt$output, full.names=FALSE),
  status       = "completed"
)
write(toJSON(summary_data, auto_unbox=TRUE), file.path(opt$output, "summary.json"))

prog(92, "Step 7/8 — Generating taxonomy & report plots...")
cat("Generating plots...\n")
# (has_ggplot2/has_reshape2/has_vegan/has_ape/has_pheatmap/has_scales are
# declared near the top of the script now — see the note there.)

# ── Palette helpers ────────────────────────────────────────────
palette20 <- c("#e74c3c","#e67e22","#f1c40f","#2ecc71","#1abc9c",
               "#3498db","#9b59b6","#e91e63","#00bcd4","#8bc34a",
               "#ff5722","#607d8b","#795548","#ffc107","#03a9f4",
               "#4caf50","#673ab7","#ff9800","#009688","#9e9e9e")

# Extended palette for large taxon sets (top 50)
make_pal <- function(n) {
  base <- c(palette20,
            "#a52a2a","#5f9ea0","#d2691e","#6495ed","#dc143c",
            "#00ced1","#ff8c00","#9400d3","#32cd32","#ff1493",
            "#1e90ff","#ffd700","#adff2f","#ff6347","#40e0d0",
            "#ee82ee","#f5deb3","#00ff7f","#87ceeb","#dda0dd",
            "#b0c4de","#ffb6c1","#7b68ee","#20b2aa","#f08080",
            "#e0ffff","#fafad2","#d3d3d3","#90ee90","#ffb347")
  pal <- rep(base, length.out=n)
  pal[seq_len(n)]
}

n_samp    <- nrow(seqtab_nochim)
samp_cols <- setNames(make_pal(n_samp), sample_names)

# ── Load metadata (if provided) ────────────────────────────────
meta_df   <- NULL
group_vec <- NULL   # named char vector: sample → primary group value
group_pal <- NULL   # named colour vector: primary group → colour
grp_col   <- NULL   # column name of primary group

if (!is.null(opt$metadata) && nchar(opt$metadata) > 0 && file.exists(opt$metadata)) {
  tryCatch({
    meta_df <- read.csv(opt$metadata, stringsAsFactors=FALSE)
    sid_col <- names(meta_df)[tolower(names(meta_df)) %in% c("sampleid","sample_id","sample","#sampleid")][1]
    if (!is.na(sid_col)) {
      rownames(meta_df) <- trimws(as.character(meta_df[[sid_col]]))
      meta_df <- meta_df[sample_names, , drop=FALSE]
    } else {
      rownames(meta_df) <- sample_names
    }
    grp_col <- names(meta_df)[tolower(names(meta_df)) %in% c("group","treatment","grp")][1]
    if (!is.na(grp_col)) {
      grp_vals <- trimws(as.character(meta_df[[grp_col]]))
      grp_vals[is.na(grp_vals)] <- ""
      if (any(nchar(grp_vals) > 0)) {
        group_vec <- setNames(grp_vals, sample_names)
        groups    <- unique(group_vec[nchar(group_vec) > 0])
        group_pal <- setNames(palette20[seq_along(groups)], groups)
        cat("  Metadata loaded:", nrow(meta_df), "samples,",
            length(groups), "group(s):", paste(groups, collapse=", "), "\n")
      } else {
        cat("  Metadata loaded but no Group values — using sample colours.\n")
      }
    } else {
      cat("  Metadata loaded but no Group column — using sample colours.\n")
    }
  }, error=function(e) cat("  [metadata] Could not load:", e$message, "\n"))
} else {
  cat("  No metadata provided — all samples treated as one group.\n")
}

# ── Detect ALL categorical metadata columns ────────────────────
# meta_cols: named list col_name → named char vector (sample → value)
# meta_pals: named list col_name → named colour vector (group → colour)
meta_cols <- list()
meta_pals <- list()

if (!is.null(meta_df)) {
  skip_pat <- c("sampleid","sample_id","sample","#sampleid")
  for (i in seq_along(names(meta_df))) {
    col <- names(meta_df)[i]
    if (tolower(col) %in% skip_pat) next
    vals <- trimws(as.character(meta_df[[col]]))
    vals[vals == "NA" | is.na(vals)] <- ""
    n_uniq <- length(unique(vals[nchar(vals) > 0]))
    # Keep as grouping if 2–20 unique values and not all unique (i.e., not a row-ID)
    if (n_uniq >= 2 && n_uniq <= 20 && n_uniq < n_samp) {
      gv  <- setNames(vals, sample_names)
      lvs <- unique(gv[nchar(gv) > 0])
      # Offset colours so each column gets a visually distinct palette
      offset <- (length(meta_cols) * 7) %% length(palette20)
      gp <- setNames(palette20[(seq_along(lvs) + offset - 1) %% length(palette20) + 1], lvs)
      meta_cols[[col]] <- gv
      meta_pals[[col]] <- gp
      cat("  Grouping column detected:", col, "(", n_uniq, "levels:", paste(lvs, collapse=", "), ")\n")
    }
  }
  cat("  Total grouping columns:", length(meta_cols), "\n")
}

# Helper: colour per sample (by primary group if available, else by sample)
point_col <- if (!is.null(group_vec) && !is.null(group_pal)) {
  gv_clean <- group_vec
  gv_clean[!gv_clean %in% names(group_pal)] <- names(group_pal)[1]
  setNames(group_pal[gv_clean], sample_names)
} else {
  samp_cols
}

# ═══════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════

# Build sample × taxon relative-abundance matrix (%) — top_n taxa
make_tax_mat <- function(level, top_n=50) {
  if (is.null(tax)) return(NULL)
  col_idx <- which(colnames(tax) == level)
  if (length(col_idx) == 0) return(NULL)
  taxa_vec <- tax[, col_idx]
  taxa_vec[is.na(taxa_vec) | taxa_vec == ""] <- "Unclassified"
  taxa_vec <- sub("^[a-z]__", "", taxa_vec)
  taxa_uniq <- unique(taxa_vec)
  mat <- sapply(taxa_uniq, function(t)
    rowSums(seqtab_nochim[, taxa_vec == t, drop=FALSE]))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow=1, dimnames=list(sample_names, taxa_uniq))
  mat <- mat[, order(colSums(mat), decreasing=TRUE), drop=FALSE]
  if (ncol(mat) > top_n) {
    other <- rowSums(mat[, (top_n+1):ncol(mat), drop=FALSE])
    mat   <- cbind(mat[, 1:top_n, drop=FALSE], Other=other)
  }
  pct <- sweep(mat, 1, rowSums(mat), "/") * 100
  pct[is.nan(pct)] <- 0
  pct
}

# Stacked bar plot (ggplot2 or base-R fallback)
tax_stacked_bar <- function(pct_mat, title_str, outfile,
                             x_labels=NULL, x_title="Sample") {
  if (is.null(pct_mat)) return(invisible(NULL))
  n_taxa <- ncol(pct_mat)
  n_x    <- nrow(pct_mat)
  cols_t <- make_pal(n_taxa)
  if (is.null(x_labels)) x_labels <- rownames(pct_mat)

  if (has_ggplot2 && has_reshape2) {
    df <- as.data.frame(pct_mat)
    df[[x_title]] <- factor(x_labels, levels=x_labels)
    long <- reshape2::melt(df, id.vars=x_title, variable.name="Taxon", value.name="Pct")
    long$Taxon <- factor(long$Taxon, levels=colnames(pct_mat))
    p <- ggplot2::ggplot(long, ggplot2::aes(x=.data[[x_title]], y=Pct, fill=Taxon)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::scale_fill_manual(values=setNames(cols_t, colnames(pct_mat))) +
      ggplot2::labs(title=title_str, x=x_title, y="Relative Abundance (%)") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                     legend.text=ggplot2::element_text(size=7),
                     legend.key.size=ggplot2::unit(0.35,"cm"),
                     legend.title=ggplot2::element_text(size=9))
    w <- max(7, n_x * 1.1 + 4)
    h <- max(6, ceiling(n_taxa / 4) * 0.4 + 5)
    ggplot2::ggsave(outfile, p, width=w, height=h, device="pdf", limitsize=FALSE)
  } else {
    w <- max(9, n_x * 1.1 + 4)
    pdf(outfile, width=w, height=6)
    par(mar=c(8, 5, 3, max(8, max(nchar(colnames(pct_mat))) * 0.5)), xpd=TRUE)
    barplot(t(pct_mat), beside=FALSE, col=cols_t, border=NA,
            main=title_str, ylab="Relative Abundance (%)",
            names.arg=x_labels, las=2, cex.names=0.8, ylim=c(0,100))
    legend(par("usr")[2]*1.02, par("usr")[4], legend=colnames(pct_mat),
           fill=cols_t, bty="n", cex=0.65, xpd=TRUE)
    dev.off()
  }
  cat("  ", basename(outfile), "\n", sep="")
}

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 1 — READ TRACKING LINE PLOT
# ═══════════════════════════════════════════════════════════════
tryCatch({
  track_df   <- as.data.frame(track)
  step_names <- colnames(track_df)

  if (has_ggplot2 && has_reshape2) {
    df <- track_df; df$Sample <- rownames(df)
    long <- reshape2::melt(df, id.vars="Sample", variable.name="Step", value.name="Reads")
    long$Step <- factor(long$Step, levels=step_names)
    p <- ggplot2::ggplot(long, ggplot2::aes(x=Step, y=Reads, group=Sample, colour=Sample)) +
      ggplot2::geom_line(linewidth=1) +
      ggplot2::geom_point(size=3) +
      ggplot2::scale_colour_manual(values=samp_cols) +
      ggplot2::labs(title="Read Counts Through Pipeline Steps",
                    x="Pipeline Step", y="Number of Reads") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1))
    ggplot2::ggsave(file.path(opt$output, "read_tracking_plot.pdf"),
                    p, width=9, height=5, device="pdf")
  } else {
    pdf(file.path(opt$output, "read_tracking_plot.pdf"), width=9, height=5)
    par(mar=c(7,5,3,2))
    matplot(t(track_df), type="b", lty=1, pch=16, lwd=2, col=samp_cols, xaxt="n",
            main="Read Counts Through Pipeline Steps",
            ylab="Number of Reads", xlab="")
    axis(1, at=seq_along(step_names), labels=step_names, las=2, cex.axis=0.85)
    if (n_samp > 1) legend("topright", legend=sample_names,
                            col=samp_cols, lty=1, lwd=2, pch=16, bty="n", cex=0.8)
    dev.off()
  }
  cat("  read_tracking_plot.pdf\n")
}, error=function(e) cat("  [skip] read_tracking plot:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 2 — QC READ-COUNT BOXPLOT
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (n_samp >= 3) {
    track_df2 <- as.data.frame(track)
    if (has_ggplot2 && has_reshape2) {
      df <- track_df2; df$Sample <- rownames(df)
      long <- reshape2::melt(df, id.vars="Sample", variable.name="Step", value.name="Reads")
      long$Step <- factor(long$Step, levels=colnames(track_df2))
      long$SampleCol <- samp_cols[long$Sample]
      p <- ggplot2::ggplot(long, ggplot2::aes(x=Step, y=Reads)) +
        ggplot2::geom_boxplot(fill="#3b82f620", outlier.shape=NA) +
        ggplot2::geom_jitter(ggplot2::aes(colour=Sample), width=0.15, size=2.5) +
        ggplot2::scale_colour_manual(values=samp_cols) +
        ggplot2::scale_y_continuous(
          labels = if (has_scales) scales::label_comma() else function(x) format(x, big.mark=",", scientific=FALSE)
        ) +
        ggplot2::labs(title="Read Count Distribution Per Pipeline Step",
                      x="Pipeline Step", y="Read Count") +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1))
      ggplot2::ggsave(file.path(opt$output, "qc_readcount_boxplot.pdf"),
                      p, width=9, height=5, device="pdf")
      cat("  qc_readcount_boxplot.pdf\n")
    }
  }
}, error=function(e) cat("  [skip] QC boxplot:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 3 — ASV LENGTH DISTRIBUTION
# ═══════════════════════════════════════════════════════════════
tryCatch({
  seq_lens <- nchar(colnames(seqtab_nochim))
  med_len  <- median(seq_lens)

  if (has_ggplot2) {
    df <- data.frame(Length=seq_lens)
    p  <- ggplot2::ggplot(df, ggplot2::aes(x=Length)) +
      ggplot2::geom_histogram(bins=50, fill="#3b82f6", colour="white", linewidth=0.3) +
      ggplot2::geom_vline(xintercept=med_len, colour="#ef4444", linetype="dashed", linewidth=1) +
      ggplot2::annotate("text", x=med_len+2, y=Inf,
                        label=paste("Median:", med_len, "bp"),
                        colour="#ef4444", vjust=2, hjust=0, size=3.5) +
      ggplot2::labs(title="ASV Length Distribution",
                    x="Sequence Length (bp)", y="Number of ASVs") +
      ggplot2::theme_bw(base_size=11)
    ggplot2::ggsave(file.path(opt$output, "asv_length_distribution.pdf"),
                    p, width=7, height=4, device="pdf")
  } else {
    pdf(file.path(opt$output, "asv_length_distribution.pdf"), width=7, height=4)
    par(mar=c(4,4,3,1))
    hist(seq_lens, breaks=50, main="ASV Length Distribution",
         xlab="Sequence Length (bp)", ylab="Number of ASVs",
         col="#3b82f6", border="white", las=1)
    abline(v=med_len, col="#ef4444", lwd=2, lty=2)
    legend("topright", legend=paste("Median:", med_len, "bp"),
           col="#ef4444", lty=2, lwd=2, bty="n")
    dev.off()
  }
  cat("  asv_length_distribution.pdf\n")
}, error=function(e) cat("  [skip] ASV length:", e$message, "\n"))

prog(95, "Step 8/8 — Diversity & statistical analyses...")

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 4 — ALPHA DIVERSITY
# ═══════════════════════════════════════════════════════════════
tryCatch({
  chao1 <- apply(seqtab_nochim, 1, function(x) {
    s <- sum(x > 0); f1 <- sum(x == 1); f2 <- sum(x == 2)
    if (f2 > 0) s + f1^2/(2*f2) else s + f1*(f1-1)/2
  })
  shannon <- apply(seqtab_nochim, 1, function(x) {
    p <- x[x > 0] / sum(x); -sum(p * log(p))
  })
  simpson <- 1 - apply(seqtab_nochim, 1, function(x) {
    p <- x / sum(x); sum(p^2)
  })
  alpha_df <- data.frame(
    Sample     = sample_names,
    TotalReads = as.integer(rowSums(seqtab_nochim)),
    Observed   = as.integer(rowSums(seqtab_nochim > 0)),
    Chao1      = round(chao1, 2),
    Shannon    = round(shannon, 4),
    Simpson    = round(simpson, 4),
    stringsAsFactors = FALSE
  )
  write.csv(alpha_df, file.path(opt$output, "alpha_diversity.csv"), row.names=FALSE)

  if (has_ggplot2 && has_reshape2) {
    long <- reshape2::melt(alpha_df[, c("Sample","Observed","Chao1","Shannon","Simpson")],
                           id.vars="Sample", variable.name="Metric", value.name="Value")
    long$Sample <- factor(long$Sample, levels=sample_names)

    # Use boxplot if metadata groups exist, else per-sample bar chart
    if (!is.null(group_vec) && !is.null(group_pal)) {
      long$Group <- factor(group_vec[as.character(long$Sample)], levels=names(group_pal))
      p <- ggplot2::ggplot(long, ggplot2::aes(x=Group, y=Value, fill=Group)) +
        ggplot2::geom_boxplot(alpha=0.7, outlier.shape=NA) +
        ggplot2::geom_jitter(ggplot2::aes(colour=Group), width=0.18, size=2.5, alpha=0.85) +
        ggplot2::scale_fill_manual(values=group_pal) +
        ggplot2::scale_colour_manual(values=group_pal) +
        ggplot2::facet_wrap(~Metric, scales="free_y", ncol=2) +
        ggplot2::labs(title=sprintf("Alpha Diversity by %s", grp_col),
                      x=grp_col, y="Value") +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1),
                       legend.position="none")
    } else {
      # No groups — single boxplot per metric showing distribution across all samples
      # with coloured jitter dots so individual samples are still identifiable
      p <- ggplot2::ggplot(long, ggplot2::aes(x=Metric, y=Value)) +
        ggplot2::geom_boxplot(fill="#93c5fd", colour="#1e40af",
                              alpha=0.6, outlier.shape=NA, width=0.45) +
        ggplot2::geom_jitter(ggplot2::aes(colour=Sample),
                             width=0.15, size=3, alpha=0.9) +
        ggplot2::scale_colour_manual(values=samp_cols) +
        ggplot2::facet_wrap(~Metric, scales="free_y", ncol=2) +
        ggplot2::labs(title="Alpha Diversity Metrics",
                      x="", y="Value", colour="Sample") +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       legend.position="right")
    }
    ggplot2::ggsave(file.path(opt$output, "alpha_diversity.pdf"),
                    p, width=9, height=7, device="pdf")
  } else {
    pdf(file.path(opt$output, "alpha_diversity.pdf"), width=11, height=6)
    par(mfrow=c(1,4), mar=c(7,4,3,1))
    for (metric in c("Observed","Chao1","Shannon","Simpson")) {
      barplot(alpha_df[[metric]], names.arg=alpha_df$Sample,
              col=samp_cols, border=NA, main=metric, ylab=metric, las=2, cex.names=0.8)
    }
    dev.off()
  }
  cat("  alpha_diversity.pdf + alpha_diversity.csv\n")
}, error=function(e) cat("  [skip] alpha diversity:", e$message, "\n"))

# ── ALPHA STATS: Pairwise Wilcoxon + BH correction ─────────────
tryCatch({
  gv_alpha <- if (!is.null(group_vec)) group_vec else NULL
  if (!is.null(gv_alpha) && length(unique(gv_alpha[sample_names])) >= 2) {
    grps_alpha <- factor(gv_alpha[sample_names])
    wilcox_list <- lapply(c("Observed","Chao1","Shannon","Simpson"), function(metric) {
      if (!metric %in% colnames(alpha_df)) return(NULL)
      vals <- setNames(alpha_df[[metric]], alpha_df$Sample)
      tryCatch({
        wt <- pairwise.wilcox.test(vals, grps_alpha, p.adjust.method="BH")
        p_mat <- wt$p.value
        rows <- rownames(p_mat); cols <- colnames(p_mat)
        do.call(rbind, lapply(seq_along(rows), function(i)
          do.call(rbind, lapply(seq_along(cols), function(j) {
            v <- p_mat[i, j]
            if (is.na(v)) return(NULL)
            data.frame(metric=metric, group1=rows[i], group2=cols[j],
                       p.adj=round(v,6), stringsAsFactors=FALSE)
          }))
        ))
      }, error=function(e) NULL)
    })
    wilcox_all <- do.call(rbind, Filter(Negate(is.null), wilcox_list))
    if (!is.null(wilcox_all) && nrow(wilcox_all) > 0) {
      write.csv(wilcox_all,
                file.path(opt$output, "alpha_diversity_stats.csv"), row.names=FALSE)
      cat("  alpha_diversity_stats.csv\n")
    }
  }
}, error=function(e) cat("  [skip] Wilcoxon alpha:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 5 — RAREFACTION CURVES (overall, coloured by sample)
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (has_vegan && n_samp >= 2) {
    min_reads <- min(rowSums(seqtab_nochim))
    step_size <- max(1, floor(min_reads / 100))

    pdf(file.path(opt$output, "rarefaction_curves.pdf"), width=9, height=6)
    vegan::rarecurve(seqtab_nochim, step=step_size, col=point_col,
                     lwd=2, label=(n_samp <= 12),
                     main="Rarefaction Curves", xlab="Reads Sampled", ylab="Observed ASVs")
    abline(v=min_reads, lty=2, col="gray50")
    legend("bottomright", paste("Min depth:", min_reads),
           lty=2, col="gray50", bty="n", cex=0.9)
    dev.off()
    cat("  rarefaction_curves.pdf\n")
  }
}, error=function(e) cat("  [skip] rarefaction:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 6 — TAXONOMY STACKED BARS (top N per-sample + per-group)
# ═══════════════════════════════════════════════════════════════
top_n_taxa <- as.integer(opt$topN)   # 30 / 50 / 100 from UI
if (!is.null(tax)) {
  for (lvl in c("Phylum","Class","Order","Family","Genus","Species")) {
    tryCatch({
      m <- make_tax_mat(lvl, top_n=top_n_taxa)
      if (is.null(m)) next
      out_f <- file.path(opt$output, paste0("taxonomy_", tolower(lvl), ".pdf"))
      tax_stacked_bar(m,
                      title_str = sprintf("Relative Abundance — Top %d %s (per sample)", top_n_taxa, lvl),
                      outfile   = out_f)
      write.csv(as.data.frame(m),
                file.path(opt$output, paste0("taxonomy_", tolower(lvl), ".csv")))
    }, error=function(e) cat("  [skip] taxonomy", lvl, ":", e$message, "\n"))
  }
}

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 7 — BETA DIVERSITY (PCoA / UPGMA / Heatmap / PERMANOVA)
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (n_samp >= 3 && has_vegan) {
    rel_ab  <- seqtab_nochim / rowSums(seqtab_nochim)
    bc_dist <- vegan::vegdist(rel_ab, method="bray")

    # ── 7a. PCoA — coloured by each metadata column ─────────────
    if (has_ape) {
      pcoa_res <- ape::pcoa(bc_dist)
      ax1    <- pcoa_res$values$Relative_eig[1] * 100
      ax2    <- pcoa_res$values$Relative_eig[2] * 100
      scores <- pcoa_res$vectors[, 1:min(2, ncol(pcoa_res$vectors)), drop=FALSE]

      # Default PCoA (by sample or primary group)
      pcoa_colour <- if (!is.null(group_vec)) group_vec else sample_names
      pcoa_pal    <- if (!is.null(group_pal)) group_pal else samp_cols

      if (has_ggplot2) {
        df_pc <- data.frame(PC1=scores[,1], PC2=scores[,2],
                            Sample=sample_names,
                            Group=factor(pcoa_colour))
        p <- ggplot2::ggplot(df_pc, ggplot2::aes(x=PC1, y=PC2, colour=Group, label=Sample)) +
          ggplot2::geom_point(size=4) +
          ggplot2::geom_text(vjust=-0.8, size=3) +
          ggplot2::scale_colour_manual(values=pcoa_pal) +
          ggplot2::labs(title="PCoA — Bray-Curtis Distance",
                        x=sprintf("PC1 (%.1f%%)", ax1),
                        y=sprintf("PC2 (%.1f%%)", ax2),
                        colour=if (!is.null(group_vec)) "Group" else "Sample") +
          ggplot2::theme_bw(base_size=11)
        ggplot2::ggsave(file.path(opt$output, "beta_pcoa.pdf"),
                        p, width=7, height=6, device="pdf")
      } else {
        pdf(file.path(opt$output, "beta_pcoa.pdf"), width=7, height=6)
        par(mar=c(5,5,3,7))
        plot(scores[,1], scores[,2], col=pcoa_pal[pcoa_colour], pch=16, cex=2,
             main="PCoA — Bray-Curtis Distance",
             xlab=sprintf("PC1 (%.1f%%)", ax1), ylab=sprintf("PC2 (%.1f%%)", ax2))
        text(scores[,1], scores[,2], labels=sample_names, pos=3, cex=0.8)
        if (!is.null(group_pal)) legend("topright", legend=names(pcoa_pal), fill=pcoa_pal, bty="n")
        dev.off()
      }
      cat("  beta_pcoa.pdf\n")

      # Additional PCoA per extra metadata column
      if (has_ggplot2 && length(meta_cols) > 0) {
        for (col_name in names(meta_cols)) {
          tryCatch({
            gv  <- meta_cols[[col_name]]
            gp  <- meta_pals[[col_name]]
            grps <- unique(gv[nchar(gv) > 0])
            grp_factor <- factor(gv[sample_names], levels=grps)
            df_pc2 <- data.frame(PC1=scores[,1], PC2=scores[,2],
                                 Sample=sample_names, Group=grp_factor)
            p2 <- ggplot2::ggplot(df_pc2, ggplot2::aes(x=PC1, y=PC2, colour=Group, label=Sample)) +
              ggplot2::geom_point(size=4) +
              ggplot2::geom_text(vjust=-0.8, size=3) +
              ggplot2::scale_colour_manual(values=gp) +
              ggplot2::labs(title=sprintf("PCoA — Bray-Curtis (coloured by %s)", col_name),
                            x=sprintf("PC1 (%.1f%%)", ax1),
                            y=sprintf("PC2 (%.1f%%)", ax2), colour=col_name) +
              ggplot2::theme_bw(base_size=11)
            safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)
            ggplot2::ggsave(file.path(opt$output, sprintf("beta_pcoa_%s.pdf", safe_col)),
                            p2, width=7, height=6, device="pdf")
            cat("  beta_pcoa_", safe_col, ".pdf\n", sep="")
          }, error=function(e) cat("  [skip] PCoA", col_name, ":", e$message, "\n"))
        }
      }
    }

    # ── 7b. UPGMA Dendrogram (colored by each metadata group) ────
    hc <- hclust(bc_dist, method="average")

    # Export distance matrix CSV
    dist_csv_path <- file.path(opt$output, "bray_curtis_distance_matrix.csv")
    write.csv(as.matrix(bc_dist), dist_csv_path)
    cat("  bray_curtis_distance_matrix.csv\n")

    # Default UPGMA (colored by primary group or all-black)
    upgma_dir <- file.path(opt$output, "UPGMA")
    dir.create(upgma_dir, showWarnings=FALSE, recursive=TRUE)

    if (has_ggplot2 && has_ape) {
      make_upgma_tree <- function(hc_obj, tip_colors, group_vec_col, col_name, out_prefix) {
        tryCatch({
          phylo_obj <- ape::as.phylo(hc_obj)
          tip_order <- phylo_obj$tip.label
          tip_col_vec <- tip_colors[tip_order]
          tip_col_vec[is.na(tip_col_vec)] <- "#444444"

          grp_vals <- group_vec_col[tip_order]
          uniq_grps <- unique(grp_vals[nchar(grp_vals) > 0])
          grp_pal_local <- tip_colors[uniq_grps]
          grp_pal_local <- grp_pal_local[!is.na(grp_pal_local)]

          w <- max(7, length(tip_order) * 0.22 + 4)
          h <- max(6, length(tip_order) * 0.28 + 2)

          # PDF
          pdf(file.path(upgma_dir, paste0(out_prefix, ".pdf")), width=w, height=h)
          par(mar=c(4, 2, 4, 8))
          plot(phylo_obj, type="phylogram", tip.color=tip_col_vec,
               cex=0.8, label.offset=0.002,
               main=sprintf("UPGMA tree colored by %s", col_name))
          if (length(grp_pal_local) > 0) {
            legend("bottomright", legend=names(grp_pal_local),
                   col=grp_pal_local, pch=19, bty="n", cex=0.8,
                   title=col_name, xpd=TRUE)
          }
          dev.off()
          cat("  UPGMA/", out_prefix, ".pdf\n", sep="")

          # Save tip colors CSV
          tc_df <- data.frame(Sample=tip_order, Group=grp_vals[tip_order],
                              Color=tip_col_vec, stringsAsFactors=FALSE)
          write.csv(tc_df, file.path(upgma_dir, paste0("tip_colors_", out_prefix, ".csv")),
                    row.names=FALSE)
        }, error=function(e) cat("  [skip] UPGMA", col_name, ":", e$message, "\n"))
      }

      if (length(meta_cols) > 0) {
        for (col_name in names(meta_cols)) {
          gv <- meta_cols[[col_name]]
          gp <- meta_pals[[col_name]]
          tip_cols_mc <- setNames(sapply(sample_names, function(s) {
            g <- gv[s]
            if (!is.na(g) && nchar(g) > 0 && g %in% names(gp)) gp[[g]] else "#888888"
          }), sample_names)
          safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)
          make_upgma_tree(hc, tip_cols_mc, gv, col_name, safe_col)
        }
      } else {
        # No metadata — plain dendrogram
        tryCatch({
          phylo_obj <- ape::as.phylo(hc)
          pdf(file.path(upgma_dir, "UPGMA_all_samples.pdf"),
              width=max(7, n_samp*0.22+4), height=max(6, n_samp*0.28+2))
          par(mar=c(4,2,4,2))
          plot(phylo_obj, type="phylogram", cex=0.8, label.offset=0.002,
               main="UPGMA Dendrogram — Bray-Curtis")
          dev.off()
          cat("  UPGMA/UPGMA_all_samples.pdf\n")
        }, error=function(e) cat("  [skip] UPGMA plain:", e$message, "\n"))
      }
    } else {
      # Fallback: base R dendrogram
      pdf(file.path(upgma_dir, "UPGMA_dendrogram.pdf"),
          width=max(7, n_samp*0.8+2), height=5)
      par(mar=c(5, 4, 3, 1))
      plot(as.dendrogram(hc), main="UPGMA Dendrogram — Bray-Curtis",
           ylab="Distance", hang=-1, las=2)
      dev.off()
      cat("  UPGMA/UPGMA_dendrogram.pdf\n")
    }
    # Also keep a top-level file for backward compat
    file.copy(if (length(meta_cols) > 0)
      file.path(upgma_dir, paste0(gsub("[^a-zA-Z0-9]","_",names(meta_cols)[1]), ".pdf"))
      else file.path(upgma_dir, "UPGMA_all_samples.pdf"),
      file.path(opt$output, "beta_upgma.pdf"), overwrite=TRUE)

    # ── 7c. Beta-diversity distance heatmap ──────────────────────
    dist_mat <- as.matrix(bc_dist)
    if (has_pheatmap && length(meta_cols) > 0) {
      ann_col <- do.call(data.frame, c(
        lapply(meta_cols, function(gv) factor(gv[sample_names])),
        list(stringsAsFactors=FALSE, row.names=sample_names)
      ))
      names(ann_col) <- names(meta_cols)
      ann_colors <- setNames(lapply(names(meta_cols), function(cn) {
        gp <- meta_pals[[cn]]
        gp[names(gp) %in% unique(ann_col[[cn]])]
      }), names(meta_cols))
      pheatmap::pheatmap(
        dist_mat,
        annotation_col    = ann_col,
        annotation_row    = ann_col,
        annotation_colors = ann_colors,
        color             = colorRampPalette(c("#ffffff","#3b82f6","#1e3a8a"))(50),
        main              = "Bray-Curtis Distance Heatmap",
        fontsize          = 8,
        filename          = file.path(opt$output, "beta_heatmap.pdf"),
        width             = max(5, n_samp*0.6+4),
        height            = max(5, n_samp*0.6+4)
      )
    } else if (has_pheatmap) {
      pheatmap::pheatmap(
        dist_mat,
        color    = colorRampPalette(c("#ffffff","#3b82f6","#1e3a8a"))(50),
        main     = "Bray-Curtis Distance Heatmap",
        fontsize = 8,
        filename = file.path(opt$output, "beta_heatmap.pdf"),
        width    = max(5, n_samp*0.6+2),
        height   = max(5, n_samp*0.6+2)
      )
    } else {
      pdf(file.path(opt$output, "beta_heatmap.pdf"),
          width=max(5, n_samp*0.6+2), height=max(5, n_samp*0.6+2))
      image(dist_mat, main="Bray-Curtis Distance Heatmap",
            col=colorRampPalette(c("white","#3b82f6","#1e3a8a"))(50), axes=FALSE)
      axis(1, at=seq(0,1,len=n_samp), labels=rownames(dist_mat), las=2, cex.axis=0.8)
      axis(2, at=seq(0,1,len=n_samp), labels=rownames(dist_mat), las=2, cex.axis=0.8)
      dev.off()
    }
    cat("  beta_heatmap.pdf\n")

    # ── 7d. PERMANOVA + beta dispersion ──────────────────────────
    tryCatch({
      sink(file.path(opt$output, "beta_stats.txt"))
      cat("=== Beta Diversity Statistics ===\n\n")
      cat("Bray-Curtis distance matrix:\n")
      print(round(dist_mat, 4))
      if (length(meta_cols) > 0) {
        for (col_name in names(meta_cols)) {
          gv <- meta_cols[[col_name]]
          grps <- unique(gv[nchar(gv) > 0])
          if (length(grps) >= 2) {
            grp_factor <- factor(gv)
            cat(sprintf("\n\n=== PERMANOVA — %s effect ===\n", col_name))
            perm_res <- vegan::adonis2(bc_dist ~ grp_factor, permutations=999)
            print(perm_res)
            cat(sprintf("\n=== Beta Dispersion — %s ===\n", col_name))
            disp_res <- vegan::betadisper(bc_dist, grp_factor)
            print(vegan::permutest(disp_res, permutations=999))
          }
        }
      } else if (!is.null(group_vec) && length(unique(group_vec)) >= 2) {
        grp_factor <- factor(group_vec)
        cat("\n=== PERMANOVA (vegan::adonis2) — Group effect ===\n")
        perm_res <- vegan::adonis2(bc_dist ~ grp_factor, permutations=999)
        print(perm_res)
        cat("\n=== Beta Dispersion (vegan::betadisper) ===\n")
        disp_res <- vegan::betadisper(bc_dist, grp_factor)
        print(vegan::permutest(disp_res, permutations=999))
      } else {
        cat("\n(Skipping PERMANOVA — need ≥2 groups with metadata)\n")
      }
      sink()
      cat("  beta_stats.txt\n")
    }, error=function(e) { if (sink.number() > 0) sink(); cat("  [skip] PERMANOVA:", e$message, "\n") })
  }
}, error=function(e) cat("  [skip] beta diversity:", e$message, "\n"))

# ── BETA EXTRA: Jaccard NMDS ────────────────────────────────────
tryCatch({
  if (n_samp >= 3 && has_vegan) {
    rel_ab_jacc <- seqtab_nochim / rowSums(seqtab_nochim)
    jacc_dist   <- vegan::vegdist(rel_ab_jacc, method="jaccard", binary=TRUE)
    set.seed(42)
    nmds_jacc   <- vegan::metaMDS(as.matrix(jacc_dist), k=2, trymax=50, trace=FALSE)
    nmds_df     <- as.data.frame(nmds_jacc$points)
    colnames(nmds_df) <- c("NMDS1","NMDS2")
    nmds_df$Sample <- rownames(nmds_df)
    if (has_ggplot2) {
      nmds_df$Group <- factor(if (!is.null(group_vec)) group_vec[nmds_df$Sample]
                               else rep("All", nrow(nmds_df)))
      pal_jacc <- if (!is.null(group_pal)) group_pal else samp_cols
      p_jacc <- ggplot2::ggplot(nmds_df,
                  ggplot2::aes(x=NMDS1, y=NMDS2, colour=Group, label=Sample)) +
        ggplot2::geom_point(size=4) +
        ggplot2::geom_text(vjust=-0.8, size=3) +
        ggplot2::scale_colour_manual(values=pal_jacc) +
        ggplot2::annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
                          label=sprintf("stress = %.4f", nmds_jacc$stress), size=3.5) +
        ggplot2::labs(title="NMDS — Binary Jaccard Distance") +
        ggplot2::theme_bw(base_size=11)
      ggplot2::ggsave(file.path(opt$output, "beta_nmds_jaccard.pdf"),
                      p_jacc, width=7, height=6, device="pdf")
    } else {
      pdf(file.path(opt$output, "beta_nmds_jaccard.pdf"), width=7, height=6)
      plot(nmds_jacc$points, pch=16, col=samp_cols, cex=2,
           main=sprintf("NMDS — Binary Jaccard (stress=%.4f)", nmds_jacc$stress))
      text(nmds_jacc$points, labels=sample_names, pos=3, cex=0.8)
      dev.off()
    }
    cat("  beta_nmds_jaccard.pdf\n")
  }
}, error=function(e) cat("  [skip] Jaccard NMDS:", e$message, "\n"))

# ── BETA EXTRA: Pairwise PERMANOVA ─────────────────────────────
tryCatch({
  if (n_samp >= 4 && has_vegan) {
    gv_beta <- if (length(meta_cols) > 0) meta_cols[[1]] else
               if (!is.null(group_vec)) setNames(group_vec, sample_names) else NULL
    if (!is.null(gv_beta)) {
      gv_vals <- gv_beta[sample_names]
      grps    <- unique(gv_vals[!is.na(gv_vals) & nchar(gv_vals) > 0])
      if (length(grps) >= 2) {
        rel_ab_pw <- seqtab_nochim / rowSums(seqtab_nochim)
        bc_dist_pw <- vegan::vegdist(rel_ab_pw, method="bray")
        pairs <- combn(grps, 2, simplify=FALSE)
        pw_rows <- do.call(rbind, lapply(pairs, function(pair) {
          subs <- names(gv_vals)[gv_vals %in% pair]
          if (length(subs) < 4) return(NULL)
          sub_d <- as.dist(as.matrix(bc_dist_pw)[subs, subs])
          sub_g <- factor(gv_vals[subs])
          tryCatch({
            pm <- vegan::adonis2(sub_d ~ sub_g, permutations=999)
            data.frame(group1=pair[1], group2=pair[2],
                       R2=round(pm$R2[1],4), F.stat=round(pm$F[1],4),
                       p=pm$`Pr(>F)`[1], p.adj=NA_real_,
                       stringsAsFactors=FALSE)
          }, error=function(e) NULL)
        }))
        if (!is.null(pw_rows) && nrow(pw_rows) > 0) {
          pw_rows$p.adj <- p.adjust(pw_rows$p, method="BH")
          write.csv(pw_rows, file.path(opt$output, "beta_pairwise_permanova.csv"),
                    row.names=FALSE)
          cat("  beta_pairwise_permanova.csv\n")
        }
      }
    }
  }
}, error=function(e) cat("  [skip] Pairwise PERMANOVA:", e$message, "\n"))

# ── BETA EXTRA: ANOSIM (global + pairwise) ─────────────────────
tryCatch({
  if (n_samp >= 3 && has_vegan) {
    gv_anosim <- if (length(meta_cols) > 0) meta_cols[[1]] else
                 if (!is.null(group_vec)) setNames(group_vec, sample_names) else NULL
    if (!is.null(gv_anosim)) {
      gv_vals_a <- gv_anosim[sample_names]
      grp_fac_a <- factor(gv_vals_a)
      if (nlevels(grp_fac_a) >= 2) {
        rel_ab_an  <- seqtab_nochim / rowSums(seqtab_nochim)
        bc_dist_an <- vegan::vegdist(rel_ab_an, method="bray")
        sink_path  <- file.path(opt$output, "beta_anosim.txt")
        sink(sink_path)
        cat("=== ANOSIM (Bray-Curtis) ===\n\n")
        anosi_res <- vegan::anosim(bc_dist_an, grp_fac_a, permutations=999)
        print(anosi_res)
        cat("\n=== Pairwise ANOSIM ===\n")
        grps_a <- unique(gv_vals_a[!is.na(gv_vals_a) & nchar(gv_vals_a) > 0])
        sink()   # close before loop to avoid nested sink issues
        pw_anosim <- do.call(rbind, lapply(combn(grps_a, 2, simplify=FALSE), function(pair) {
          subs <- names(gv_vals_a)[gv_vals_a %in% pair]
          if (length(subs) < 4) return(NULL)
          sub_d <- as.dist(as.matrix(bc_dist_an)[subs, subs])
          sub_g <- factor(gv_vals_a[subs])
          tryCatch({
            an <- vegan::anosim(sub_d, sub_g, permutations=999)
            data.frame(group1=pair[1], group2=pair[2],
                       R=round(an$statistic,4), p=an$signif,
                       stringsAsFactors=FALSE)
          }, error=function(e) NULL)
        }))
        if (!is.null(pw_anosim) && nrow(pw_anosim) > 0) {
          pw_anosim$p.adj <- p.adjust(pw_anosim$p, method="BH")
          # Append pairwise results to text file
          write("\n=== Pairwise ANOSIM (BH-corrected) ===", sink_path, append=TRUE)
          write(capture.output(print(pw_anosim)), sink_path, append=TRUE)
          write.csv(pw_anosim, file.path(opt$output, "beta_anosim_pairwise.csv"),
                    row.names=FALSE)
          cat("  beta_anosim_pairwise.csv\n")
        }
        cat("  beta_anosim.txt\n")
      }
    }
  }
}, error=function(e) {
  if (sink.number() > 0) sink()
  cat("  [skip] ANOSIM:", e$message, "\n")
})

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 7b — DESeq2 Differential Abundance
# ═══════════════════════════════════════════════════════════════
tryCatch({
  has_phyloseq_da <- requireNamespace("phyloseq", quietly=TRUE)
  has_deseq2      <- requireNamespace("DESeq2",   quietly=TRUE)
  gv_da <- if (!is.null(group_vec)) group_vec else NULL
  if (has_phyloseq_da && has_deseq2 && !is.null(gv_da) &&
      length(unique(gv_da[sample_names])) >= 2 && n_samp >= 4) {

    suppressPackageStartupMessages({
      library(phyloseq); library(DESeq2)
    })
    cat("[DESeq2] Building phyloseq...\n")

    # Build phyloseq from seqtab_nochim
    OTU_ps  <- phyloseq::otu_table(t(seqtab_nochim), taxa_are_rows=TRUE)
    meta_ps <- data.frame(Group=factor(gv_da[sample_names]), row.names=sample_names)
    SAMP_ps <- phyloseq::sample_data(meta_ps)
    if (!is.null(tax)) {
      TAX_ps <- phyloseq::tax_table(as.matrix(tax))
      ps_da  <- phyloseq::phyloseq(OTU_ps, TAX_ps, SAMP_ps)
    } else {
      ps_da  <- phyloseq::phyloseq(OTU_ps, SAMP_ps)
    }

    # Collapse at genus if rank available
    if (!is.null(tax) && "Genus" %in% colnames(tax)) {
      ps_da <- phyloseq::tax_glom(ps_da, taxrank="Genus", NArm=FALSE)
    }

    # Filter: remove features with >90% zeros
    ps_filt <- phyloseq::prune_taxa(
      rowSums(phyloseq::otu_table(ps_da) == 0) <
        ncol(phyloseq::otu_table(ps_da)) * 0.9,
      ps_da
    )
    if (phyloseq::ntaxa(ps_filt) < 2) stop("Too few taxa after filtering for DESeq2")

    cat(sprintf("[DESeq2] %d taxa × %d samples\n",
                phyloseq::ntaxa(ps_filt), phyloseq::nsamples(ps_filt)))

    ps_ds <- phyloseq::phyloseq_to_deseq2(ps_filt, ~ Group)
    ds    <- DESeq2::estimateSizeFactors(ps_ds, type="poscounts")
    ds    <- DESeq2::DESeq(ds, test="Wald", fitType="parametric", quiet=TRUE)
    res   <- DESeq2::results(ds, alpha=0.05)
    res_df <- as.data.frame(res)
    res_df$taxon <- rownames(res_df)

    # Annotate labels
    if (!is.null(tax) && "Genus" %in% colnames(tax)) {
      lab <- tax[res_df$taxon, "Genus"]
      lab[is.na(lab) | lab == ""] <- res_df$taxon[is.na(lab) | lab == ""]
    } else {
      lab <- res_df$taxon
    }
    res_df$label <- sub("^[a-z]__", "", lab)
    res_df       <- res_df[order(res_df$padj, na.last=NA), ]

    write.csv(res_df, file.path(opt$output, "deseq2_results.csv"), row.names=FALSE)
    cat("  deseq2_results.csv\n")

    if (has_ggplot2) {
      # ── Volcano ──
      vdf <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
      vdf$Significant   <- vdf$padj < 0.05
      vdf$neg_log10_p   <- -log10(vdf$padj + 1e-300)
      top_lab_df        <- head(vdf[vdf$Significant, ], 15)

      p_volc <- ggplot2::ggplot(vdf,
                  ggplot2::aes(x=log2FoldChange, y=neg_log10_p, colour=Significant)) +
        ggplot2::geom_point(alpha=0.7, size=2) +
        ggplot2::scale_colour_manual(values=c("grey60","#ef4444")) +
        ggplot2::geom_hline(yintercept=-log10(0.05), linetype=2, colour="#6b7280") +
        ggplot2::geom_vline(xintercept=0,            linetype=2, colour="#6b7280") +
        ggplot2::labs(title="DESeq2 — Differential Abundance",
                      x="Log2 Fold Change", y="-log10(padj)") +
        ggplot2::theme_bw(base_size=11)
      if (requireNamespace("ggrepel", quietly=TRUE) && nrow(top_lab_df) > 0)
        p_volc <- p_volc + ggrepel::geom_text_repel(
          data=top_lab_df,
          ggplot2::aes(label=label), size=2.8, show.legend=FALSE)
      ggplot2::ggsave(file.path(opt$output, "deseq2_volcano.pdf"),
                      p_volc, width=8, height=6, device="pdf")
      cat("  deseq2_volcano.pdf\n")

      # ── Heatmap of top 20 ──
      sig_taxa <- head(res_df$taxon, 20)
      if (length(sig_taxa) >= 2 && has_pheatmap) {
        ps_rel  <- phyloseq::transform_sample_counts(ps_da,
                     function(x) x / sum(x) * 100)
        ps_sig  <- phyloseq::prune_taxa(
                     intersect(sig_taxa, phyloseq::taxa_names(ps_rel)), ps_rel)
        mat     <- as.matrix(phyloseq::otu_table(ps_sig))
        if (!phyloseq::taxa_are_rows(ps_sig)) mat <- t(mat)

        if (!is.null(tax) && "Genus" %in% colnames(tax)) {
          rl <- sub("^[a-z]__", "", tax[rownames(mat), "Genus"])
          rl[is.na(rl) | rl == ""] <- rownames(mat)[is.na(rl) | rl == ""]
          rownames(mat) <- rl
        }

        ann_col_da <- data.frame(
          Group = meta_ps[colnames(mat), "Group"],
          row.names = colnames(mat)
        )
        n_lvl    <- nlevels(ann_col_da$Group)
        grp_pal_da <- if (!is.null(group_pal)) group_pal[seq_len(n_lvl)] else
                        setNames(scales::hue_pal()(n_lvl), levels(ann_col_da$Group))

        # Phylum annotation if available
        if (!is.null(tax) && "Phylum" %in% colnames(tax)) {
          matched <- intersect(sig_taxa, rownames(tax))
          phy_vec <- sub("^[a-z]__", "", tax[matched, "Phylum"])
          phy_vec[is.na(phy_vec)] <- "Unknown"
          phy_uniq <- unique(phy_vec)
          n_phy    <- length(phy_uniq)
          phy_cols <- if (requireNamespace("RColorBrewer", quietly=TRUE))
            setNames(RColorBrewer::brewer.pal(max(3, min(n_phy,12)),"Paired")[seq_len(n_phy)], phy_uniq)
          else setNames(scales::hue_pal()(n_phy), phy_uniq)
          ann_row_da <- data.frame(Phylum=phy_vec, row.names=rownames(mat)[seq_along(phy_vec)])
          ann_colors_da <- list(Group=grp_pal_da, Phylum=phy_cols)
        } else {
          ann_row_da    <- NULL
          ann_colors_da <- list(Group=grp_pal_da)
        }

        pheatmap::pheatmap(
          mat,
          scale             = "row",
          annotation_col    = ann_col_da,
          annotation_row    = ann_row_da,
          annotation_colors = ann_colors_da,
          color             = colorRampPalette(c("#3b82f6","white","#ef4444"))(100),
          main              = "DESeq2 — Top 20 Differential Taxa (row-scaled %)",
          fontsize          = 8,
          filename          = file.path(opt$output, "deseq2_heatmap.pdf"),
          width             = max(7, n_samp * 0.5 + 4),
          height            = max(6, length(sig_taxa) * 0.4 + 3)
        )
        cat("  deseq2_heatmap.pdf\n")
      }
    }

  } else {
    cat("  [DESeq2] Skipped (needs DESeq2 + phyloseq + >=2 groups + >=4 samples)\n")
  }
}, error=function(e) cat("  [skip] DESeq2:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 8 — PREVALENCE vs MEAN ABUNDANCE (Genus, all samples)
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (!is.null(tax) && n_samp >= 2) {
    col_idx <- which(colnames(tax) == "Genus")
    if (length(col_idx) > 0) {
      taxa_vec  <- sub("^[a-z]__", "", tax[, col_idx])
      taxa_vec[is.na(taxa_vec) | taxa_vec == ""] <- "Unclassified"
      rel_ab    <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/")
      taxa_uniq <- unique(taxa_vec)
      mat_genus <- sapply(taxa_uniq, function(t)
        rowSums(rel_ab[, taxa_vec == t, drop=FALSE]))
      if (is.null(dim(mat_genus)))
        mat_genus <- matrix(mat_genus, nrow=1, dimnames=list(sample_names, taxa_uniq))

      prev    <- colMeans(mat_genus > 0) * 100
      mean_ab <- colMeans(mat_genus) * 100
      pv_df   <- data.frame(Genus=taxa_uniq, Prevalence=prev,
                             MeanRelAbundance=mean_ab, stringsAsFactors=FALSE)

      if (has_ggplot2) {
        top_label <- head(taxa_uniq[order(mean_ab, decreasing=TRUE)], 10)
        pv_df$Label <- ifelse(pv_df$Genus %in% top_label, pv_df$Genus, "")
        p <- ggplot2::ggplot(pv_df[pv_df$MeanRelAbundance > 0, ],
                             ggplot2::aes(x=MeanRelAbundance, y=Prevalence, label=Label)) +
          ggplot2::geom_point(alpha=0.6, colour="#3b82f6", size=2.5) +
          ggplot2::geom_text(size=2.8, vjust=-0.6, colour="#1e293b") +
          ggplot2::scale_x_log10() +
          ggplot2::labs(title="Prevalence vs Mean Relative Abundance (Genus)",
                        x="Mean Relative Abundance (%, log scale)",
                        y="Prevalence (% samples)") +
          ggplot2::theme_bw(base_size=11)
        ggplot2::ggsave(file.path(opt$output, "prevalence_abundance.pdf"),
                        p, width=8, height=6, device="pdf")
      } else {
        pdf(file.path(opt$output, "prevalence_abundance.pdf"), width=8, height=6)
        par(mar=c(5,5,3,2))
        plot(mean_ab, prev, log="x", pch=16, col="#3b82f680", cex=1.5,
             main="Prevalence vs Mean Relative Abundance (Genus)",
             xlab="Mean Relative Abundance (%, log scale)",
             ylab="Prevalence (% samples)")
        top10 <- order(mean_ab, decreasing=TRUE)[1:min(10, length(mean_ab))]
        text(mean_ab[top10], prev[top10], labels=taxa_uniq[top10], pos=3, cex=0.7)
        dev.off()
      }
      cat("  prevalence_abundance.pdf\n")
    }
  }
}, error=function(e) cat("  [skip] prevalence plot:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 9 — OBSERVED ASVs PER SAMPLE
# ═══════════════════════════════════════════════════════════════
tryCatch({
  obs <- rowSums(seqtab_nochim > 0)
  if (has_ggplot2) {
    df <- data.frame(Sample=factor(sample_names, levels=sample_names),
                     ASVs=as.integer(obs))
    p  <- ggplot2::ggplot(df, ggplot2::aes(x=Sample, y=ASVs, fill=Sample)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::scale_fill_manual(values=point_col) +
      ggplot2::labs(title="Observed ASVs per Sample", x="Sample", y="Observed ASVs") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                     legend.position="none")
    ggplot2::ggsave(file.path(opt$output, "observed_asvs.pdf"),
                    p, width=max(5, n_samp*0.8+2), height=5, device="pdf")
  } else {
    pdf(file.path(opt$output, "observed_asvs.pdf"), width=max(5, n_samp*0.8+2), height=5)
    par(mar=c(7,5,3,2))
    barplot(obs, names.arg=sample_names, col=point_col, border=NA,
            main="Observed ASVs per Sample", ylab="Observed ASVs", las=2, cex.names=0.85)
    dev.off()
  }
  cat("  observed_asvs.pdf\n")
}, error=function(e) cat("  [skip] observed ASVs:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 10 — GROUP MEAN STACKED BARS (per metadata col × level)
# ═══════════════════════════════════════════════════════════════
if (length(meta_cols) > 0 && !is.null(tax)) {
  cat("Generating group mean stacked bars...\n")
  for (col_name in names(meta_cols)) {
    gv   <- meta_cols[[col_name]]
    grps <- unique(gv[nchar(gv) > 0])
    if (length(grps) < 2) next
    safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)

    for (lvl in c("Phylum","Class","Order","Family","Genus")) {
      tryCatch({
        col_idx <- which(colnames(tax) == lvl)
        if (length(col_idx) == 0) next

        taxa_vec  <- sub("^[a-z]__", "", tax[, col_idx])
        taxa_vec[is.na(taxa_vec) | taxa_vec == ""] <- "Unclassified"
        rel_ab    <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/") * 100
        taxa_uniq <- unique(taxa_vec)

        mat <- sapply(taxa_uniq, function(t)
          rowSums(rel_ab[, taxa_vec == t, drop=FALSE]))
        if (is.null(dim(mat)))
          mat <- matrix(mat, nrow=1, dimnames=list(sample_names, taxa_uniq))

        # Mean abundance per group
        grp_mat <- do.call(rbind, lapply(grps, function(g) {
          idx <- which(gv == g)
          if (length(idx) == 0) return(rep(0, ncol(mat)))
          if (length(idx) == 1) return(mat[idx, ])
          colMeans(mat[idx, , drop=FALSE])
        }))
        rownames(grp_mat) <- grps

        # Top N by overall mean (uses opt$topN from UI)
        top_n  <- min(top_n_taxa, ncol(grp_mat))
        ord    <- order(colMeans(grp_mat), decreasing=TRUE)
        if (ncol(grp_mat) > top_n) {
          other  <- rowSums(grp_mat[, ord[(top_n+1):ncol(grp_mat)], drop=FALSE])
          grp_mat <- cbind(grp_mat[, ord[1:top_n], drop=FALSE], Other=other)
        } else {
          grp_mat <- grp_mat[, ord, drop=FALSE]
        }

        out_f <- file.path(opt$output,
                           sprintf("group_mean_%s_%s.pdf", safe_col, tolower(lvl)))
        tax_stacked_bar(
          grp_mat,
          title_str = sprintf("Mean Relative Abundance by %s — Top %d %s",
                               col_name, top_n_taxa, lvl),
          outfile   = out_f,
          x_title   = col_name
        )
      }, error=function(e) cat("  [skip] group mean", col_name, lvl, ":", e$message, "\n"))
    }
  }
}

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 11 — TAXONOMY HEATMAPS WITH METADATA ANNOTATION
# ═══════════════════════════════════════════════════════════════
if (!is.null(tax) && has_pheatmap) {
  cat("Generating taxonomy heatmaps...\n")

  # Build annotation dataframe from all metadata columns
  ann_col_hm    <- NULL
  ann_colors_hm <- list()
  if (length(meta_cols) > 0) {
    ann_col_hm <- data.frame(
      lapply(meta_cols, function(gv) factor(gv[sample_names])),
      row.names = sample_names, stringsAsFactors = FALSE
    )
    names(ann_col_hm) <- names(meta_cols)
    for (cn in names(meta_cols)) {
      gp <- meta_pals[[cn]]
      lvs <- levels(ann_col_hm[[cn]])
      ann_colors_hm[[cn]] <- gp[lvs[lvs %in% names(gp)]]
    }
  }

  for (hm_lvl in c("Genus","Family","Phylum")) {
    tryCatch({
      col_idx_hm <- which(colnames(tax) == hm_lvl)
      if (length(col_idx_hm) == 0) next

      tv <- sub("^[a-z]__", "", tax[, col_idx_hm])
      tv[is.na(tv) | tv == ""] <- "Unclassified"
      rel_hm <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/") * 100
      tu <- unique(tv)
      mat_hm <- sapply(tu, function(t) rowSums(rel_hm[, tv == t, drop=FALSE]))
      if (is.null(dim(mat_hm)))
        mat_hm <- matrix(mat_hm, nrow=1, dimnames=list(sample_names, tu))

      top_n_hm <- if (hm_lvl == "Phylum") min(top_n_taxa, ncol(mat_hm)) else min(top_n_taxa, ncol(mat_hm))
      top_taxa_hm <- names(sort(colMeans(mat_hm), decreasing=TRUE))[1:top_n_hm]
      heat_mat_hm <- t(mat_hm[, top_taxa_hm, drop=FALSE])   # taxa × samples

      hm_colors <- switch(hm_lvl,
        Genus   = colorRampPalette(c("#f0f4ff","#3b82f6","#1e1b4b"))(100),
        Family  = colorRampPalette(c("#fff7ed","#f97316","#431407"))(100),
        Phylum  = colorRampPalette(c("#f0fff4","#22c55e","#14532d"))(100)
      )
      hm_title <- sprintf("Top %d %s — Relative Abundance Heatmap", top_n_hm, hm_lvl)
      hm_file  <- file.path(opt$output, sprintf("taxonomy_heatmap_%s.pdf", tolower(hm_lvl)))
      hm_w     <- max(8, n_samp * 0.6 + 4 + length(meta_cols) * 0.8)
      hm_h     <- max(8, top_n_hm * 0.35 + 3)

      pheatmap::pheatmap(
        heat_mat_hm,
        annotation_col    = ann_col_hm,
        annotation_colors = if (length(ann_colors_hm) > 0) ann_colors_hm else NULL,
        color             = hm_colors,
        scale             = "row",
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        main              = hm_title,
        fontsize_row      = max(5, min(9, 200/top_n_hm)),
        fontsize_col      = 8,
        filename          = hm_file,
        width             = hm_w,
        height            = hm_h
      )
      cat("  ", basename(hm_file), "\n", sep="")
    }, error=function(e) cat("  [skip] heatmap", hm_lvl, ":", e$message, "\n"))
  }
}

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 12 — GROUPED ALPHA DIVERSITY BOXPLOTS
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (length(meta_cols) > 0 && has_ggplot2 && has_reshape2 &&
      file.exists(file.path(opt$output, "alpha_diversity.csv"))) {
    cat("Generating grouped alpha boxplots...\n")
    alpha_df2 <- read.csv(file.path(opt$output, "alpha_diversity.csv"),
                          stringsAsFactors=FALSE)

    for (col_name in names(meta_cols)) {
      gv <- meta_cols[[col_name]]
      gp <- meta_pals[[col_name]]
      grps <- unique(gv[nchar(gv) > 0])
      if (length(grps) < 2) next

      alpha_df2[["Group_"]] <- gv[alpha_df2$Sample]
      alpha_df2[["Group_"]][is.na(alpha_df2[["Group_"]])] <- ""
      sub_df <- alpha_df2[nchar(alpha_df2[["Group_"]]) > 0, ]
      sub_df[["Group_"]] <- factor(sub_df[["Group_"]], levels=grps)

      long <- reshape2::melt(
        sub_df[, c("Sample","Group_","Observed","Chao1","Shannon","Simpson")],
        id.vars = c("Sample","Group_"),
        variable.name = "Metric", value.name = "Value"
      )

      safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)
      out_f    <- file.path(opt$output, sprintf("alpha_boxplot_%s.pdf", safe_col))

      p <- ggplot2::ggplot(long, ggplot2::aes(x=Group_, y=Value, fill=Group_)) +
        ggplot2::geom_boxplot(alpha=0.7, outlier.shape=NA) +
        ggplot2::geom_jitter(ggplot2::aes(colour=Group_), width=0.18, size=2.5, alpha=0.85) +
        ggplot2::scale_fill_manual(values=gp) +
        ggplot2::scale_colour_manual(values=gp) +
        ggplot2::facet_wrap(~Metric, scales="free_y", ncol=2) +
        ggplot2::labs(title=sprintf("Alpha Diversity by %s", col_name),
                      x=col_name, y="Value") +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1),
                       legend.position="none")
      ggplot2::ggsave(out_f, p,
                      width=max(8, length(grps)*1.8+3), height=8, device="pdf")
      cat("  ", basename(out_f), "\n", sep="")
    }
  }
}, error=function(e) cat("  [skip] grouped alpha boxplots:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 13 — GROUPED RAREFACTION CURVES (coloured by group)
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (has_vegan && length(meta_cols) > 0 && n_samp >= 2) {
    cat("Generating grouped rarefaction curves...\n")
    min_reads <- min(rowSums(seqtab_nochim))
    step_size <- max(1, floor(min_reads / 100))

    for (col_name in names(meta_cols)) {
      gv   <- meta_cols[[col_name]]
      gp   <- meta_pals[[col_name]]
      grps <- unique(gv[nchar(gv) > 0])
      if (length(grps) < 2) next

      cols_samp <- setNames(sapply(sample_names, function(s) {
        g <- gv[s]
        if (nchar(g) > 0 && g %in% names(gp)) gp[[g]] else "#aaaaaa"
      }), sample_names)

      safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)
      out_f    <- file.path(opt$output, sprintf("rarefaction_%s.pdf", safe_col))

      pdf(out_f, width=9, height=6)
      vegan::rarecurve(seqtab_nochim, step=step_size, col=cols_samp, lwd=2,
                       label=(n_samp <= 12),
                       main=sprintf("Rarefaction Curves — coloured by %s", col_name),
                       xlab="Reads Sampled", ylab="Observed ASVs")
      abline(v=min_reads, lty=2, col="gray50")
      legend("bottomright", legend=grps, col=gp[grps], lty=1, lwd=2,
             bty="n", cex=0.9, title=col_name)
      dev.off()
      cat("  ", basename(out_f), "\n", sep="")
    }
  }
}, error=function(e) cat("  [skip] grouped rarefaction:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 14 — PREVALENCE BY GROUP (faceted per metadata col)
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (!is.null(tax) && n_samp >= 2 && has_ggplot2 && length(meta_cols) > 0) {
    cat("Generating prevalence by group plots...\n")
    col_idx_pv <- which(colnames(tax) == "Genus")
    if (length(col_idx_pv) > 0) {
      tv_pv  <- sub("^[a-z]__", "", tax[, col_idx_pv])
      tv_pv[is.na(tv_pv) | tv_pv == ""] <- "Unclassified"
      rel_pv <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/")
      tu_pv  <- unique(tv_pv)
      mat_pv <- sapply(tu_pv, function(t) rowSums(rel_pv[, tv_pv == t, drop=FALSE]))
      if (is.null(dim(mat_pv)))
        mat_pv <- matrix(mat_pv, nrow=1, dimnames=list(sample_names, tu_pv))

      for (col_name in names(meta_cols)) {
        tryCatch({
          gv   <- meta_cols[[col_name]]
          gp   <- meta_pals[[col_name]]
          grps <- unique(gv[nchar(gv) > 0])
          if (length(grps) < 2) next

          prev_list <- lapply(grps, function(g) {
            idx <- which(gv == g)
            if (length(idx) == 0) return(NULL)
            sub_mat <- mat_pv[idx, , drop=FALSE]
            if (is.null(dim(sub_mat))) sub_mat <- matrix(sub_mat, nrow=1)
            data.frame(
              Genus            = tu_pv,
              Prevalence       = colMeans(sub_mat > 0) * 100,
              MeanRelAbundance = colMeans(sub_mat) * 100,
              Group            = g, stringsAsFactors = FALSE
            )
          })
          prev_df <- do.call(rbind, Filter(Negate(is.null), prev_list))
          prev_df$Group <- factor(prev_df$Group, levels=grps)

          # Label top 8 per group by mean abundance
          prev_df$Label <- ""
          for (g in grps) {
            idx2 <- which(prev_df$Group == g & prev_df$MeanRelAbundance > 0)
            if (length(idx2) == 0) next
            top8 <- idx2[order(prev_df$MeanRelAbundance[idx2], decreasing=TRUE)[1:min(8, length(idx2))]]
            prev_df$Label[top8] <- prev_df$Genus[top8]
          }

          safe_col <- gsub("[^a-zA-Z0-9]", "_", col_name)
          out_f    <- file.path(opt$output, sprintf("prevalence_by_%s.pdf", safe_col))
          n_cols_facet <- min(3, length(grps))

          p <- ggplot2::ggplot(prev_df[prev_df$MeanRelAbundance > 0, ],
                               ggplot2::aes(x=MeanRelAbundance, y=Prevalence,
                                            colour=Group, label=Label)) +
            ggplot2::geom_point(alpha=0.7, size=2) +
            ggplot2::geom_text(size=2.4, vjust=-0.6, show.legend=FALSE) +
            ggplot2::scale_x_log10() +
            ggplot2::scale_colour_manual(values=gp) +
            ggplot2::facet_wrap(~Group, ncol=n_cols_facet) +
            ggplot2::labs(
              title = sprintf("Prevalence vs Mean Rel. Abundance by %s (Genus)", col_name),
              x     = "Mean Relative Abundance (%, log scale)",
              y     = "Prevalence (% samples)"
            ) +
            ggplot2::theme_bw(base_size=10) +
            ggplot2::theme(legend.position="none")
          ggplot2::ggsave(out_f, p,
                          width  = max(8, n_cols_facet * 3.5 + 1),
                          height = max(5, ceiling(length(grps)/n_cols_facet) * 3.5),
                          device = "pdf", limitsize=FALSE)
          cat("  ", basename(out_f), "\n", sep="")
        }, error=function(e) cat("  [skip] prevalence by", col_name, ":", e$message, "\n"))
      }
    }
  }
}, error=function(e) cat("  [skip] prevalence by group:", e$message, "\n"))

# ═══════════════════════════════════════════════════════════════
#  PLOT SECTION 15 — ASV RICHNESS BOXPLOT BY GROUP
# ═══════════════════════════════════════════════════════════════
tryCatch({
  if (length(meta_cols) > 0 && has_ggplot2) {
    cat("Generating ASV richness by group...\n")
    obs_df <- data.frame(Sample=sample_names,
                         ObservedASVs=as.integer(rowSums(seqtab_nochim > 0)),
                         stringsAsFactors=FALSE)
    # Merge metadata
    if (nrow(metadata_df) > 0) {
      obs_df <- merge(obs_df, metadata_df, by.x="Sample", by.y="SampleID", all.x=TRUE)
    }
    for (col_name in meta_cols) {
      tryCatch({
        out_f <- file.path(plots_dir, paste0("ASV_richness_by_", col_name, ".pdf"))
        obs_df[[col_name]] <- factor(obs_df[[col_name]])
        p <- ggplot2::ggplot(obs_df, ggplot2::aes_string(x=col_name, y="ObservedASVs", fill=col_name)) +
          ggplot2::geom_boxplot(outlier.size=1.5, alpha=0.85) +
          ggplot2::geom_jitter(width=0.15, size=1.2, alpha=0.6) +
          ggplot2::labs(title=paste("ASV Richness by", col_name),
                        x=col_name, y="Observed ASVs") +
          ggplot2::theme_bw(base_size=13) +
          ggplot2::theme(legend.position="none")
        ggplot2::ggsave(out_f, p, width=7, height=5, device="pdf")
        cat("  ", basename(out_f), "\n", sep="")
      }, error=function(e) cat("  [skip] ASV richness by", col_name, ":", e$message, "\n"))
    }
  }
}, error=function(e) cat("  [skip] ASV richness by group:", e$message, "\n"))

# =============================================================
#  SECTION -- TAX4FUN2 FUNCTIONAL PREDICTION (optional)
# =============================================================
if (isTRUE(opt$tax4fun)) {
  cat("\n--- Tax4Fun2 Functional Prediction ---\n")
  tryCatch({
    if (!requireNamespace("Tax4Fun2", quietly=TRUE)) {
      cat("  [skip] Tax4Fun2 package not installed.\n")
      cat("  Install: Rscript -e \"install.packages('Tax4Fun2',\n")
      cat("    repos=c('https://bwemheu.r-universe.dev','https://cloud.r-project.org'))\"\n")
    } else {
      tf2_dir <- file.path(opt$output_dir, "Tax4Fun2")
      dir.create(tf2_dir, recursive=TRUE, showWarnings=FALSE)

      # Export ASV sequences as FASTA
      asv_seqs <- colnames(seqtab_nochim)
      asv_ids  <- paste0("ASV", seq_along(asv_seqs))
      fasta_lines <- c(rbind(paste0(">", asv_ids), asv_seqs))
      asv_fasta <- file.path(tf2_dir, "ASVs.fasta")
      writeLines(fasta_lines, asv_fasta)
      cat("  ASV FASTA exported\n")

      # Export abundance table
      abund_df <- data.frame(SampleID=rownames(seqtab_nochim),
                             as.data.frame(seqtab_nochim),
                             stringsAsFactors=FALSE)
      colnames(abund_df)[-1] <- asv_ids
      abund_csv <- file.path(tf2_dir, "ASV_table.csv")
      write.csv(abund_df, abund_csv, row.names=FALSE)
      cat("  ASV table exported\n")

      # Reference data path
      ref_dir <- if (nchar(opt$tax4fun_ref) > 0) opt$tax4fun_ref else
        file.path(path.expand("~"), "r16s-app", "databases",
                  "Tax4Fun2_ReferenceData_v2")

      if (!dir.exists(ref_dir)) {
        cat("  [skip] Tax4Fun2 reference data not found at:", ref_dir, "\n")
        cat("  Download: https://zenodo.org/record/6327578\n")
      } else {
        cat("  Running Tax4Fun2...\n")
        Tax4Fun2::Tax4Fun2(
          file_path_otu_table = abund_csv,
          file_path_ref_data  = ref_dir,
          path_to_working_dir = tf2_dir,
          use_parallel        = FALSE
        )
        cat("  Tax4Fun2 complete. Results in:", tf2_dir, "\n")
      }
    }
  }, error=function(e) cat("  [skip] Tax4Fun2:", e$message, "\n"))
}

# =============================================================
#  EXTRA VIZ — interactive chart CSVs for Edit Charts
#  PCoA, NMDS, Jaccard, Rarefaction, ASV lengths
# =============================================================
tryCatch({
  extra_script <- file.path(dirname(sub("--file=", "",
    grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])),
    "dada2_extra_viz.R")
  if (!is.na(extra_script) && file.exists(extra_script)) {
    source(extra_script)
    dada2_extra_viz(opt$output_dir)
  }
}, error=function(e) cat("  [skip] Extra viz:", e$message, "\n"))

# =============================================================
#  PHYLOGENETIC TREE — export ASV FASTA + build NJ/ML tree
#  Produces phylo_tree.nwk for circular tree visualization
# =============================================================
tryCatch({
  cat("\n── Phylogenetic Tree ────────────────────────────────────────────\n")
  asv_seqs <- colnames(seqtab_nochim)
  n_asvs   <- length(asv_seqs)
  cat("  ASVs to align:", n_asvs, "\n")

  # Limit to top 500 ASVs by total reads (performance)
  if (n_asvs > 500) {
    asv_totals <- colSums(seqtab_nochim)
    top_idx    <- order(asv_totals, decreasing=TRUE)[1:500]
    asv_seqs   <- asv_seqs[top_idx]
    cat("  Limiting to top 500 ASVs by abundance\n")
  }

  # Write FASTA file
  fasta_out <- file.path(opt$output_dir, "asvs.fasta")
  asv_ids   <- paste0("ASV", seq_along(asv_seqs))
  fasta_lines <- character(length(asv_seqs) * 2)
  for (i in seq_along(asv_seqs)) {
    fasta_lines[2*i - 1] <- paste0(">", asv_ids[i])
    fasta_lines[2*i]     <- asv_seqs[i]
  }
  writeLines(fasta_lines, fasta_out)
  cat("  ✓ asvs.fasta written (", length(asv_seqs), "sequences)\n")

  # ── Try MAFFT + FastTree (if installed) ──────────────────────────────────
  mafft_bin    <- Sys.which("mafft")
  fasttree_bin <- Sys.which("FastTree")
  if (nchar(fasttree_bin) == 0) fasttree_bin <- Sys.which("fasttree")

  tree_nwk <- file.path(opt$output_dir, "phylo_tree.nwk")

  if (nchar(mafft_bin) > 0 && nchar(fasttree_bin) > 0) {
    aln_out <- file.path(opt$output_dir, "asvs_aligned.fasta")
    cat("  Running MAFFT alignment...\n")
    ret_mafft <- system2(mafft_bin,
                         args = c("--auto", "--thread", "-1", "--quiet", fasta_out),
                         stdout = aln_out, stderr = FALSE)
    if (ret_mafft == 0 && file.exists(aln_out)) {
      cat("  ✓ MAFFT alignment done\n")
      cat("  Running FastTree...\n")
      ret_ft <- system2(fasttree_bin,
                        args = c("-nt", "-gtr", "-quiet", aln_out),
                        stdout = tree_nwk, stderr = FALSE)
      if (ret_ft == 0 && file.exists(tree_nwk)) {
        cat("  ✓ FastTree phylogenetic tree built:", tree_nwk, "\n")
      } else {
        cat("  [warn] FastTree failed — falling back to NJ\n")
        file.remove(tree_nwk)
      }
    } else {
      cat("  [warn] MAFFT failed — falling back to NJ\n")
    }
  }

  # ── Fallback: NJ tree from k-mer distances (no external tools needed) ─────
  if (!file.exists(tree_nwk) && has_ape) {
    cat("  Building NJ tree from k-mer distances (ape)...\n")
    # Use ape::dist.dna requires DNAbin — compute simple edit distance instead
    # Build character matrix from ASV sequences
    seqs_char <- strsplit(asv_seqs, "")
    maxlen    <- max(sapply(seqs_char, length))
    # Pad shorter sequences
    seqs_pad  <- lapply(seqs_char, function(s) c(s, rep("-", maxlen - length(s))))
    seq_mat   <- do.call(rbind, seqs_pad)
    rownames(seq_mat) <- asv_ids
    # Convert to DNAbin
    dna_bin   <- as.DNAbin(seq_mat)
    # Compute distance
    d_mat     <- tryCatch(dist.dna(dna_bin, model="K80", pairwise.deletion=TRUE),
                          error=function(e) dist.dna(dna_bin, model="raw",
                                                     pairwise.deletion=TRUE))
    d_mat[!is.finite(d_mat)] <- 0.5
    nj_tree   <- nj(d_mat)
    write.tree(nj_tree, file=tree_nwk)
    cat("  ✓ NJ phylogenetic tree built (k-mer distance):", tree_nwk, "\n")
  }

  if (!file.exists(tree_nwk))
    cat("  [skip] Could not build phylogenetic tree\n")

}, error=function(e) cat("  [skip] Tree building failed:", e$message, "\n"))

# =============================================================
#  DONE
# =============================================================
cat("\n========================================\n")
cat("  dada2_pipeline.R completed successfully\n")
cat("========================================\n")
cat("Output directory:", opt$output_dir, "\n")
