# =============================================================================
#  maaslin3_run.R — multivariable differential abundance with MaAsLin 3
#
#  Runs on any finished amplicon result folder (DADA2, ITS, COX1, PacBio, Emu,
#  VSEARCH): counts are summed to one taxonomic rank (Genus by default), and every
#  feature is tested for an association with each metadata variable — in
#  ABUNDANCE (how much, where present) and in PREVALENCE (present or absent),
#  adjusted for all the other variables at once. That is what the single-group
#  tests elsewhere in the app (DESeq2, ANCOM-BC2, LEfSe) cannot do.
#
#    Rscript maaslin3_run.R --output_dir RESULT [--metadata metadata.csv]
#            [--fixed group,age] [--random subject] [--reference "group,Control"]
#            [--rank Genus] [--min_prevalence 0.1] [--threads 4]
#
#  Writes RESULT/maaslin3/ (MaAsLin's own output and figures) and
#  RESULT/r_tables/maaslin3_results.csv (all tests, read by the report).
#
#  Function tables (shotgun + HUMAnN): --feature_table func/pathways_unstratified.tsv
#  --tag pathways tests that table (features x samples, CPM; HUMAnN's UNMAPPED /
#  UNINTEGRATED / UNGROUPED rows left out) instead of the taxonomy, into
#  RESULT/maaslin3/func_<tag>/ and RESULT/r_tables/maaslin3_func_<tag>.csv
#  (skip reason: RESULT/maaslin3/SKIPPED_<tag>.txt).
#  Exits 0 with a message (never fails a run) when the design cannot be tested:
#  no metadata, no variable with replicated levels, or MaAsLin 3 not installed.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--output_dir", type = "character"),
  make_option("--metadata",   type = "character", default = ""),
  make_option("--fixed",      type = "character", default = ""),
  make_option("--random",     type = "character", default = ""),
  make_option("--reference",  type = "character", default = ""),
  make_option("--rank",       type = "character", default = "Genus"),
  make_option("--min_prevalence", type = "double", default = 0.1),
  make_option("--max_significance", type = "double", default = 0.1),
  make_option("--exclude",    type = "character", default = ""),   # control / mock sample IDs
  make_option("--threads",    type = "integer", default = 1),
  make_option("--feature_table", type = "character", default = ""),  # function table instead of taxa
  make_option("--tag",        type = "character", default = "")
)))

`%||%` <- function(a, b) if (is.null(a)) b else a
say  <- function(...) cat("[maaslin3]", ..., "\n")
skip <- function(why) { say("skipped:", why)
  dir.create(file.path(opt$output_dir, "maaslin3"), showWarnings = FALSE)
  writeLines(why, file.path(opt$output_dir, "maaslin3",
                            if (nzchar(opt$tag)) sprintf("SKIPPED_%s.txt", opt$tag) else "SKIPPED.txt"))
  quit(status = 0) }

if (!requireNamespace("maaslin3", quietly = TRUE))
  skip("MaAsLin 3 is not installed (BiocManager::install('maaslin3'); needs R >= 4.4)")

od <- opt$output_dir
md_file <- if (nzchar(opt$metadata)) opt$metadata else file.path(od, "metadata.csv")
if (!file.exists(md_file)) skip("no metadata.csv — nothing to test against")

# ── counts: features x samples, whatever the pipeline's layout ───────────────
if (nzchar(opt$feature_table)) {
  ft_file <- opt$feature_table
  if (!file.exists(ft_file)) ft_file <- file.path(od, opt$feature_table)
  if (!file.exists(ft_file)) skip(paste("no table", opt$feature_table))
  ft <- if (grepl("\\.csv$", ft_file)) read.csv(ft_file, check.names = FALSE, stringsAsFactors = FALSE) else
        read.delim(ft_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  feats <- as.character(ft[[1]])
  fm <- as.matrix(ft[, -1, drop = FALSE]); storage.mode(fm) <- "numeric"; fm[is.na(fm)] <- 0
  colnames(fm) <- sub("(_Abundance(-RPKs|-CPM)?|_Coverage|-RPKs|-CPM)$", "", colnames(fm))
  keepf <- !(sub(":.*$", "", feats) %in% c("UNMAPPED", "UNINTEGRATED", "UNGROUPED")) & rowSums(fm) > 0
  if (sum(keepf) < 2) skip("fewer than two features in the function table")
  rownames(fm) <- make.unique(feats)
  counts <- as.data.frame(t(fm[keepf, , drop = FALSE]), check.names = FALSE)
  rank <- opt$tag
  say(sprintf("%d samples, %d function features (%s)", nrow(counts), ncol(counts), opt$tag))
} else {
asv <- read.csv(file.path(od, "asv_table.csv"), check.names = FALSE, stringsAsFactors = FALSE)
tax_file <- if (file.exists(file.path(od, "taxonomy_table.csv"))) file.path(od, "taxonomy_table.csv") else
            file.path(od, "taxonomy.csv")
tax <- read.csv(tax_file, check.names = FALSE, stringsAsFactors = FALSE)
first <- tolower(names(asv)[1])
if (first %in% c("sample_id", "sampleid", "sample")) {            # samples as rows (old COX1/PacBio)
  m <- t(as.matrix(asv[, -1, drop = FALSE])); colnames(m) <- asv[[1]]
  ids <- rownames(m)
} else if ("sequence" %in% names(asv)) {                          # DADA2 family
  ids <- asv$sequence
  m <- as.matrix(asv[, setdiff(names(asv), "sequence"), drop = FALSE])
} else {                                                          # Emu / VSEARCH: first column = id
  ids <- as.character(asv[[1]])
  m <- as.matrix(asv[, -1, drop = FALSE])
}
storage.mode(m) <- "numeric"; m[is.na(m)] <- 0
rownames(m) <- ids
key <- as.character(tax[[1]])
rank <- opt$rank
if (!rank %in% names(tax)) skip(paste("the taxonomy has no", rank, "column"))
lab <- tax[[rank]][match(ids, key)]
lab <- sub("^[a-z]__", "", lab)
lab[is.na(lab) | !nzchar(trimws(lab)) | lab %in% c("NA", "Unclassified", "unidentified")] <- NA
keep <- !is.na(lab)
if (sum(keep) < 2) skip(paste("fewer than two features named at", rank))
agg <- rowsum(m[keep, , drop = FALSE], lab[keep])                 # rank x samples
counts <- as.data.frame(t(agg), check.names = FALSE)              # samples x features
say(sprintf("%d samples, %d %s-level features (%.1f%% of reads named at %s)",
            nrow(counts), ncol(counts), tolower(rank), 100 * sum(agg) / max(1, sum(m)), tolower(rank)))
}

# ── metadata ──────────────────────────────────────────────────────────────────
md <- read.csv(md_file, check.names = FALSE, stringsAsFactors = FALSE)
idc <- which(tolower(names(md)) %in% c("sampleid", "sample_id", "sample", "#sampleid"))[1]
if (is.na(idc)) idc <- 1
rownames(md) <- as.character(md[[idc]])
md <- md[, -idc, drop = FALSE]
common <- intersect(rownames(counts), rownames(md))
if (length(common) < 4) skip(sprintf("only %d samples have metadata (need at least 4)", length(common)))
counts <- counts[common, , drop = FALSE]; md <- md[common, , drop = FALSE]
# control / mock samples are not part of the biological comparison
ctl_col <- which(tolower(names(md)) %in% c("control", "sample_type", "type"))[1]
if (!is.na(ctl_col)) {
  ctl <- tolower(as.character(md[[ctl_col]])) %in% c("negative", "blank", "ntc", "neg", "mock", "positive")
  if (any(ctl)) { say(sprintf("left out %d control/mock sample(s)", sum(ctl)))
    counts <- counts[!ctl, , drop = FALSE]; md <- md[!ctl, , drop = FALSE] }
}

if (nzchar(opt$exclude)) {
  ex <- intersect(trimws(strsplit(opt$exclude, ",")[[1]]), rownames(md))
  if (length(ex)) { say(sprintf("left out %d control/mock sample(s) named in the run", length(ex)))
    keep_s <- !rownames(md) %in% ex
    counts <- counts[keep_s, , drop = FALSE]; md <- md[keep_s, , drop = FALSE] }
}
if (nrow(md) < 4) skip(sprintf("only %d biological samples left (need at least 4)", nrow(md)))

usable <- function(v) {                     # a variable worth testing
  x <- md[[v]]; x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) < 4) return(FALSE)
  if (is.numeric(x)) return(length(unique(x)) > 2)
  tb <- table(x); length(tb) >= 2 && sum(tb >= 2) >= 2 && length(tb) < length(x)
}
skipcols <- c("description", "desc", "notes", "control", "sample_type", "type")
fixed <- if (nzchar(opt$fixed)) trimws(strsplit(opt$fixed, ",")[[1]]) else
         setdiff(names(md), skipcols)
fixed <- fixed[fixed %in% names(md)]
random <- if (nzchar(opt$random)) trimws(strsplit(opt$random, ",")[[1]]) else character(0)
random <- random[random %in% names(md)]
dropped <- fixed[!vapply(fixed, usable, logical(1))]
fixed <- setdiff(fixed[vapply(fixed, usable, logical(1))], random)
if (length(dropped)) say("not tested (no replicated levels):", paste(dropped, collapse = ", "))
if (length(fixed) == 0) skip("no metadata variable has at least two levels with two or more samples each")
for (v in fixed) if (!is.numeric(md[[v]])) md[[v]] <- factor(md[[v]])

ref <- NULL
if (nzchar(opt$reference)) ref <- trimws(strsplit(opt$reference, ";")[[1]])
# multi-level factors need a reference; default: the first level alphabetically
for (v in fixed) if (is.factor(md[[v]]) && nlevels(md[[v]]) > 2 &&
                     !any(startsWith(ref %||% character(0), paste0(v, ",")))) {
  ref <- c(ref, paste0(v, ",", levels(md[[v]])[1]))
}

out <- file.path(od, "maaslin3")
dir.create(out, showWarnings = FALSE)
if (nzchar(opt$tag)) { out <- file.path(out, paste0("func_", opt$tag)); dir.create(out, showWarnings = FALSE) }
say("fixed effects:", paste(fixed, collapse = ", "),
    if (length(random)) paste("| random:", paste(random, collapse = ", ")) else "")
fit <- tryCatch(maaslin3::maaslin3(
  input_data = counts, input_metadata = md, output = out,
  fixed_effects = fixed, random_effects = if (length(random)) random else NULL,
  reference = ref, normalization = "TSS", transform = "LOG",
  min_prevalence = opt$min_prevalence, max_significance = opt$max_significance,
  augment = TRUE, standardize = TRUE, median_comparison_abundance = TRUE,
  plot_summary_plot = TRUE, plot_associations = FALSE, cores = max(1, opt$threads),
  verbosity = "WARN"), error = function(e) { say("MaAsLin 3 error:", conditionMessage(e)); NULL })
if (is.null(fit)) skip("MaAsLin 3 stopped with an error (see the log)")

res <- file.path(out, "all_results.tsv")
if (file.exists(res)) {
  r <- read.delim(res, check.names = FALSE, stringsAsFactors = FALSE)
  r$rank <- rank
  dir.create(file.path(od, "r_tables"), showWarnings = FALSE)
  write.csv(r, file.path(od, "r_tables", if (nzchar(opt$tag)) sprintf("maaslin3_func_%s.csv", opt$tag) else
                                         "maaslin3_results.csv"), row.names = FALSE)
  qcol <- intersect(c("qval_individual", "qval_joint"), names(r))[1]
  nsig <- if (!is.na(qcol)) sum(r[[qcol]] < opt$max_significance, na.rm = TRUE) else NA
  say(sprintf("%d tests, %s with q < %.2f", nrow(r), nsig, opt$max_significance))
  writeLines(jsonlite::toJSON(list(fixed = fixed, random = random, reference = ref, rank = rank,
                                   n_samples = nrow(counts), n_features = ncol(counts),
                                   not_tested = dropped, max_significance = opt$max_significance,
                                   version = as.character(utils::packageVersion("maaslin3"))),
                              auto_unbox = TRUE, pretty = TRUE),
             file.path(out, "design.json"))
}
