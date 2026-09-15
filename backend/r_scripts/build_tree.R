#!/usr/bin/env Rscript
# =============================================================================
#  build_tree.R — Build phylo_tree.nwk from asvs.fasta
#
#  Usage: Rscript build_tree.R <output_dir> [threads]
#
#  WHY THIS IS A SEPARATE SCRIPT
#  ----------------------------
#  Tree building is the most memory-hungry optional step in the pipeline, and
#  the libraries it uses (DECIPHER's aligner, phangorn's ML optimiser, ape's
#  distance code) allocate in C. When they run out of memory they do not throw
#  an R condition that tryCatch() can catch — they abort the whole R process
#  with "An irrecoverable exception occurred. R is aborting now". Run inline,
#  that killed dada2_pipeline.R outright at the very last step, after all the
#  real results had already been written, so a finished run was reported as a
#  failure.
#
#  Running it as its own process means the worst case is a missing tree file
#  and a warning line, never a lost run.
#
#  Methods, best first:
#    1. MAFFT + FastTree      (external binaries — fastest and best)
#    2. DECIPHER alignment    (+ optional phangorn ML refinement)
#    3. NJ on gap-padded seqs (approximate; no real alignment)
#
#  Writes: <output_dir>/phylo_tree.nwk
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()
THREADS <- if (length(args) >= 2) suppressWarnings(as.integer(args[2])) else 1L
if (is.na(THREADS) || THREADS < 1) THREADS <- 1L

fasta_in <- file.path(out_dir, "asvs.fasta")
tree_nwk <- file.path(out_dir, "phylo_tree.nwk")

cat("── Phylogenetic Tree ────────────────────────────────────────────\n")

if (!file.exists(fasta_in)) {
  cat("  [skip] asvs.fasta not found in", out_dir, "\n")
  quit(status = 0)
}

has_ape <- requireNamespace("ape", quietly = TRUE)
if (!has_ape) {
  cat("  [skip] 'ape' not installed — cannot build a tree\n")
  quit(status = 0)
}
suppressPackageStartupMessages(library(ape))

# ── Read the ASV FASTA ───────────────────────────────────────────────────────
lines <- readLines(fasta_in, warn = FALSE)
hdr   <- grepl("^>", lines)
asv_ids  <- sub("^>", "", lines[hdr])
asv_seqs <- lines[!hdr]
if (length(asv_ids) != length(asv_seqs) || length(asv_ids) < 3) {
  cat("  [skip] need at least 3 sequences to build a tree (found",
      length(asv_ids), ")\n")
  quit(status = 0)
}
n_tips <- length(asv_seqs)
cat("  Sequences:", n_tips, "\n")

# ── How much memory is actually free? ────────────────────────────────────────
# The heavy paths below are exactly what ran the machine out of memory before,
# so gate them on real headroom rather than hoping.
mem_avail_gb <- function() {
  tryCatch({
    v <- suppressWarnings(as.numeric(system(
      "awk '/MemAvailable/ {printf \"%.2f\", $2/1048576}' /proc/meminfo",
      intern = TRUE, ignore.stderr = TRUE)))
    if (length(v) == 0 || is.na(v[1])) NA_real_ else v[1]
  }, error = function(e) NA_real_)
}
mem_gb <- mem_avail_gb()
cat("  Available RAM:", if (is.na(mem_gb)) "unknown" else sprintf("%.1f GB", mem_gb), "\n")

# ── 1. MAFFT + FastTree ──────────────────────────────────────────────────────
mafft_bin    <- Sys.which("mafft")
fasttree_bin <- Sys.which("FastTree")
if (nchar(fasttree_bin) == 0) fasttree_bin <- Sys.which("fasttree")

if (nchar(mafft_bin) > 0 && nchar(fasttree_bin) > 0) {
  aln_out <- file.path(out_dir, "asvs_aligned.fasta")
  cat("  Running MAFFT alignment...\n")
  ret_mafft <- suppressWarnings(system2(mafft_bin,
      args = c("--auto", "--thread", "-1", "--quiet", fasta_in),
      stdout = aln_out, stderr = FALSE))
  if (ret_mafft == 0 && file.exists(aln_out) && file.size(aln_out) > 0) {
    cat("  ✓ MAFFT alignment done\n  Running FastTree...\n")
    ret_ft <- suppressWarnings(system2(fasttree_bin,
        args = c("-nt", "-gtr", "-quiet", aln_out),
        stdout = tree_nwk, stderr = FALSE))
    if (ret_ft == 0 && file.exists(tree_nwk) && file.size(tree_nwk) > 0) {
      cat("  ✓ FastTree phylogenetic tree built:", tree_nwk, "\n")
      quit(status = 0)
    }
    cat("  [warn] FastTree failed — falling back\n")
    if (file.exists(tree_nwk)) file.remove(tree_nwk)
  } else {
    cat("  [warn] MAFFT failed — falling back\n")
  }
}

# ── 2. DECIPHER alignment (+ optional ML) ────────────────────────────────────
# Needs real headroom: AlignSeqs holds the full profile in memory and
# optim.pml holds per-site likelihoods for every tip.
MIN_GB_ALIGN <- 1.5
MIN_GB_ML    <- 3.0
MAX_TIPS_ML  <- 150

can_align <- requireNamespace("DECIPHER", quietly = TRUE) &&
             requireNamespace("Biostrings", quietly = TRUE) &&
             (is.na(mem_gb) || mem_gb >= MIN_GB_ALIGN)

if (!file.exists(tree_nwk) && can_align) {
  # processors=1 when memory is tight: each worker holds its own buffers.
  procs <- if (!is.na(mem_gb) && mem_gb < 4) 1L else THREADS
  cat(sprintf("  Aligning with DECIPHER (%d processor%s)...\n",
              procs, if (procs == 1) "" else "s"))
  ok <- tryCatch({
    dna_set <- Biostrings::DNAStringSet(setNames(asv_seqs, asv_ids))
    aligned <- DECIPHER::AlignSeqs(dna_set, verbose = FALSE, processors = procs)
    aln_mat <- as.DNAbin(as.matrix(as.character(aligned)))
    rownames(aln_mat) <- asv_ids
    rm(dna_set, aligned); gc(verbose = FALSE)

    d_mat <- tryCatch(dist.dna(aln_mat, model = "K80", pairwise.deletion = TRUE),
                      error = function(e) dist.dna(aln_mat, model = "raw",
                                                   pairwise.deletion = TRUE))
    finite_max <- suppressWarnings(max(d_mat[is.finite(d_mat)], na.rm = TRUE))
    if (!is.finite(finite_max)) finite_max <- 0.5
    d_mat[!is.finite(d_mat)] <- max(finite_max, 0.5)
    tr <- nj(d_mat)

    mem_now <- mem_avail_gb()
    if (requireNamespace("phangorn", quietly = TRUE) &&
        n_tips <= MAX_TIPS_ML &&
        (!is.na(mem_now) && mem_now >= MIN_GB_ML)) {
      tryCatch({
        cat("  Refining with phangorn (ML, GTR+G)...\n")
        phy_dat <- phangorn::as.phyDat(aln_mat)
        fit <- phangorn::pml(phangorn::midpoint(tr), data = phy_dat)
        fit <- phangorn::optim.pml(fit, model = "GTR", optGamma = TRUE,
                                   rearrangement = "NNI",
                                   control = phangorn::pml.control(trace = 0))
        tr <- fit$tree
        cat("  ✓ ML tree optimised\n")
      }, error = function(e) cat("  [warn] ML refinement skipped:", e$message, "\n"))
    } else if (n_tips > MAX_TIPS_ML) {
      cat(sprintf("  ML refinement skipped (%d tips > %d)\n", n_tips, MAX_TIPS_ML))
    } else if (!is.na(mem_now) && mem_now < MIN_GB_ML) {
      cat(sprintf("  ML refinement skipped (%.1f GB free < %.1f GB needed)\n",
                  mem_now, MIN_GB_ML))
    }

    write.tree(tr, file = tree_nwk)
    cat("  ✓ Phylogenetic tree built (DECIPHER alignment):", tree_nwk, "\n")
    TRUE
  }, error = function(e) {
    cat("  [warn] DECIPHER tree failed:", e$message, "\n")
    FALSE
  })
  if (isTRUE(ok)) quit(status = 0)
} else if (!file.exists(tree_nwk) && !is.na(mem_gb) && mem_gb < MIN_GB_ALIGN) {
  cat(sprintf("  DECIPHER skipped (%.1f GB free < %.1f GB needed)\n",
              mem_gb, MIN_GB_ALIGN))
}

# ── 3. NJ from gap-padded sequences (last resort, cheap) ─────────────────────
if (!file.exists(tree_nwk)) {
  cat("  Building NJ tree from padded sequences (ape)...\n")
  cat("  NOTE: sequences are not properly aligned — this tree is approximate.\n")
  cat("  Install MAFFT + FastTree, or the DECIPHER R package, for a real alignment.\n")
  tryCatch({
    seqs_char <- strsplit(asv_seqs, "")
    maxlen    <- max(vapply(seqs_char, length, integer(1)))
    seqs_pad  <- lapply(seqs_char, function(s) c(s, rep("-", maxlen - length(s))))
    seq_mat   <- do.call(rbind, seqs_pad)
    rownames(seq_mat) <- asv_ids
    dna_bin   <- as.DNAbin(seq_mat)
    d_mat     <- tryCatch(dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE),
                          error = function(e) dist.dna(dna_bin, model = "raw",
                                                       pairwise.deletion = TRUE))
    d_mat[!is.finite(d_mat)] <- 0.5
    write.tree(nj(d_mat), file = tree_nwk)
    cat("  ✓ NJ phylogenetic tree built (approximate, unaligned):", tree_nwk, "\n")
  }, error = function(e) cat("  [skip] Tree building failed:", e$message, "\n"))
}

if (!file.exists(tree_nwk)) cat("  [skip] Could not build phylogenetic tree\n")
quit(status = 0)
