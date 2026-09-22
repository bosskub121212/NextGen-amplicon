#!/usr/bin/env Rscript
# =============================================================================
#  add_species.R — addSpecies() in a process of its own
#
#  assignTaxonomy() and addSpecies() never need their references at the same
#  time, but run inline they end up holding them at the same time anyway: R's
#  allocator does not hand freed memory back to the OS, so the trainset's
#  footprint is still resident when the species file starts loading. Peak
#  becomes the SUM of the two instead of the larger of them.
#
#  With SILVA 144 that sum is what breaks a small machine:
#
#      trainset       v138.2 133.5 MB   ->  v144 191.8 MB
#      assignSpecies  v138.2  65.8 MB   ->  v144 141.3 MB
#
#  On an 8 GB laptop (about 2.3 GB actually free inside WSL) the v144 pair does
#  not fit, and R spends the difference in garbage collection: measured over
#  four real runs of the same sample, v144 took 49 and 54 minutes where v138.2
#  took 17.5 and 18 — three times slower for a reference only 1.4x larger,
#  because the slowdown is thrashing, not work.
#
#  A child process starts with an empty heap, loads only the species file, and
#  gives everything back on exit. Peak becomes max(trainset, species).
#
#  The taxonomy table is small — a few hundred rows of sequence plus rank
#  strings — so passing it through an RDS file costs nothing.
#
#  Usage:
#    Rscript add_species.R <tax.rds> <species_ref.fa.gz> <out.rds> [helpers.R]
#
#  Exit codes: 0 written, 1 failed (the caller keeps the genus-level table).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("usage: add_species.R <tax.rds> <ref.fa.gz> <out.rds> [helpers.R]\n")
  quit(status = 1)
}
tax_rds <- args[1]; ref_fa <- args[2]; out_rds <- args[3]
helpers <- if (length(args) >= 4) args[4] else ""

ok <- tryCatch({
  suppressPackageStartupMessages(library(dada2))
  if (nzchar(helpers) && file.exists(helpers)) source(helpers)

  tax <- readRDS(tax_rds)
  if (is.null(tax) || nrow(tax) == 0) stop("empty taxonomy table")

  # tax_add_species() also strips SILVA 144's "--other" genus suffix, which
  # DADA2's matchGenera() would otherwise treat as a genus mismatch and use to
  # discard species that matched at 100%. Fall back to a plain addSpecies()
  # when the helpers are not on disk, so this still works standalone.
  res <- if (exists("tax_add_species", mode = "function")) {
    tax_add_species(tax, ref_fa)
  } else {
    addSpecies(tax, ref_fa)
  }

  saveRDS(res, out_rds)
  # gc() returns six columns: used, (Mb), gc trigger, (Mb), max used, (Mb) — and
  # "max used" is a CELL COUNT, not bytes. Column 6 is the same figure already
  # converted to Mb for both Ncells and Vcells, which is what we want; treating
  # the counts as bytes reports numbers roughly 600x too large.
  peak_mb <- tryCatch(sum(gc()[, 6]), error = function(e) NA_real_)
  if (!is.na(peak_mb))
    cat(sprintf("  [species] child process peak memory: %.0f MB\n", peak_mb))
  TRUE
}, error = function(e) {
  cat("  [species] failed in child process:", conditionMessage(e), "\n")
  FALSE
})

quit(status = if (isTRUE(ok)) 0 else 1)
