# =============================================================================
#  qc_helpers.R — read-length profiling, run BEFORE filterAndTrim
#
#  filterAndTrim() discards outright any read shorter than its truncLen. That is
#  correct behaviour and it reports nothing: the reads are simply gone, and the
#  only symptom is a "filtered" count that came out lower than expected, several
#  minutes later, with no indication of why.
#
#  The maximum usable truncLen is NOT 250 minus the primer length. It is a
#  property of the trimmed reads, and the arithmetic can be off:
#
#    On a real 2x250 V3-V4 library, 250 minus the 21 bp reverse primer says 229.
#    The reads said otherwise — after cutadapt, 91.5% of reverse reads were
#    exactly 228 bp, because cutadapt removed 22 bases, not 21 (the primer as
#    entered was missing its leading G, so the match began one base in). At
#    truncLen_R = 229 only 7.2% of reverse reads survived, and the run fell from
#    22,607 filtered pairs to 1,831 — a wasted hour, over one base.
#
#  So measure instead of deriving. The distribution is bimodal: a tight cluster
#  of properly trimmed full-length reads, and a tail of short fragments. The
#  recommendation is the largest truncLen that keeps nearly all of the first
#  group, found by ignoring the short tail and then walking down until the
#  retained fraction crosses KEEP_TARGET.
# =============================================================================

QC_LONG_FRACTION <- 0.60   # reads below this fraction of the longest are "short tail"
QC_KEEP_TARGET   <- 0.90   # keep at least this share of the long reads


# Length histogram of a gzipped or plain FASTQ. Reads in chunks so a large file
# does not have to be held in memory; `max_reads` caps the work on big runs.
qc_length_hist <- function(path, max_reads = 400000L) {
  if (!file.exists(path)) return(NULL)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  tab <- integer(0); seen <- 0L; line_no <- 0L
  repeat {
    chunk <- readLines(con, n = 40000L, warn = FALSE)
    if (length(chunk) == 0) break
    # sequence lines are every 4th starting at line 2 of the file
    idx <- which(((line_no + seq_along(chunk)) %% 4L) == 2L)
    if (length(idx) > 0) {
      L <- nchar(chunk[idx], type = "bytes")
      t2 <- table(L)
      for (k in names(t2)) {
        key <- as.character(k)
        tab[key] <- (if (is.na(tab[key])) 0L else tab[key]) + as.integer(t2[[k]])
      }
      seen <- seen + length(idx)
    }
    line_no <- line_no + length(chunk)
    if (seen >= max_reads) break
  }
  if (seen == 0L) return(NULL)
  lens <- as.integer(names(tab)); cnt <- as.integer(tab)
  ord  <- order(lens)
  list(len = lens[ord], count = cnt[ord], n = seen, capped = seen >= max_reads)
}


# Largest truncLen that keeps >= QC_KEEP_TARGET of the full-length reads.
qc_recommend_trunclen <- function(h) {
  if (is.null(h)) return(NULL)
  max_len  <- max(h$len)
  long_min <- floor(max_len * QC_LONG_FRACTION)
  long_ok  <- h$len >= long_min
  n_long   <- sum(h$count[long_ok])
  if (n_long == 0L) return(NULL)
  # keep(L) = long reads at least L bases long; monotonically decreasing in L
  keep <- rev(cumsum(rev(h$count)))            # keep[i] = reads with len >= len[i]
  ok   <- long_ok & (keep / n_long) >= QC_KEEP_TARGET
  if (!any(ok)) return(NULL)
  rec <- max(h$len[ok])
  list(
    recommended = rec,
    max_len     = max_len,
    n_reads     = h$n,
    n_long      = n_long,
    keep_at_rec = keep[which(h$len == rec)[1]],
    capped      = isTRUE(h$capped),
    # the few lengths around the recommendation, which is where the cliff is
    profile     = lapply(which(h$len >= rec - 3L & h$len <= rec + 4L), function(i)
      list(len = h$len[i], reads = h$count[i],
           keep = keep[i], keep_pct = round(100 * keep[i] / n_long, 2)))
  )
}


# How many reads a given truncLen would keep, as a share of the long reads.
qc_keep_at <- function(h, L) {
  if (is.null(h) || is.na(L) || L <= 0) return(NA_real_)
  long_min <- floor(max(h$len) * QC_LONG_FRACTION)
  n_long   <- sum(h$count[h$len >= long_min])
  if (n_long == 0L) return(NA_real_)
  round(100 * sum(h$count[h$len >= L]) / n_long, 2)
}
