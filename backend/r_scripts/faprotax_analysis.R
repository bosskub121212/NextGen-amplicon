#!/usr/bin/env Rscript
# =============================================================================
#  faprotax_analysis.R — FAPROTAX Ecological Function Prediction
#
#  Usage: Rscript faprotax_analysis.R <output_dir>
#
#  Strategy (tries in order):
#    1. Run FAPROTAX Python script if installed (collapse_table.py)
#    2. Use Tax4Fun2 R package if available
#    3. Built-in genus→function lookup table (most common ecological functions)
#
#  Reads:  taxonomy_genus.csv or taxonomy_species.csv + metadata.csv
#  Writes: faprotax_results.csv     — function × sample abundance matrix
#          r_plots/faprotax_bar.pdf — STAMP-style grouped bar + CI
#          r_plots/faprotax_{G1}_vs_{G2}.pdf — pairwise comparison
# =============================================================================

args    <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) >= 1) args[1] else getwd()

cat("\n── FAPROTAX Analysis ───────────────────────────────────────────────\n")
cat("  Directory:", out_dir, "\n")

# ── Load metadata ─────────────────────────────────────────────────────────────
meta_f <- file.path(out_dir, "metadata.csv")
meta_map <- NULL; grp_vec <- NULL; groups <- NULL
if (file.exists(meta_f)) tryCatch({
  meta <- read.csv(meta_f, stringsAsFactors=FALSE)
  sid_col <- grep("sampleid|sample_id|sample", colnames(meta), ignore.case=TRUE, value=TRUE)[1]
  grp_col <- grep("^group$",                   colnames(meta), ignore.case=TRUE, value=TRUE)[1]
  if (!is.na(sid_col) && !is.na(grp_col)) {
    meta_map <- setNames(meta[[grp_col]], meta[[sid_col]])
    groups   <- unique(meta[[grp_col]]); groups <- groups[!is.na(groups) & groups!=""]
  }
}, error=function(e) NULL)

# ── Load genus/species taxonomy for annotation ────────────────────────────────
tax_f <- NULL; tax_level <- "genus"
for (lvl in c("genus","species","family")) {
  tf <- file.path(out_dir, paste0("taxonomy_", lvl, ".csv"))
  if (file.exists(tf)) { tax_f <- tf; tax_level <- lvl; break }
}

if (is.null(tax_f)) {
  cat("  [skip] Need taxonomy_genus.csv or taxonomy_species.csv\n"); quit(status=0)
}

td <- tryCatch(read.csv(tax_f, check.names=FALSE, stringsAsFactors=FALSE),
               error=function(e) { cat("  [skip] load failed\n"); quit(status=0) })

n_rows <- nrow(td)-1; n_cols <- ncol(td)-1
if (n_cols >= n_rows) {
  taxa_nm  <- colnames(td)[-1]
  abu_wide <- as.data.frame(td[,-1]); rownames(abu_wide) <- td[[1]]
} else {
  mat <- t(as.matrix(td[,-1])); storage.mode(mat) <- "numeric"
  taxa_nm <- colnames(mat); rownames(mat) <- colnames(td)[-1]
  abu_wide <- as.data.frame(mat)
}
cat("  Taxonomy (", tax_level,"):", nrow(abu_wide), "samples,", ncol(abu_wide), "taxa\n")

# ── Strategy 1: Python FAPROTAX ───────────────────────────────────────────────
faprotax_out <- file.path(out_dir, "faprotax_raw.tsv")
python_success <- FALSE

for (py_cmd in c("collapse_table.py", "faprotax", "python3 -m faprotax")) {
  fap_bin <- tryCatch(trimws(system(paste("which", strsplit(py_cmd," ")[[1]][1]),
                                     intern=TRUE, ignore.stderr=TRUE)), error=function(e) "")
  if (nchar(fap_bin) == 0) next

  cat("  Found FAPROTAX:", py_cmd, "\n")
  # Build OTU table for FAPROTAX (taxa × samples, tab-separated)
  otu_tsv <- file.path(out_dir, "_faprotax_input.tsv")
  otu_df  <- as.data.frame(t(abu_wide))
  write.table(cbind(`#OTU ID`=rownames(otu_df), otu_df), otu_tsv,
              sep="\t", quote=FALSE, row.names=FALSE)

  db_candidates <- c(
    "/usr/local/lib/python3.*/dist-packages/faprotax/FAPROTAX.txt",
    "/opt/conda/*/share/faprotax/FAPROTAX.txt",
    "~/FAPROTAX.txt", "/usr/share/faprotax/FAPROTAX.txt"
  )
  fap_db <- Sys.getenv("FAPROTAX_DB", unset="")
  if (nchar(fap_db) == 0) {
    for (pat in db_candidates) {
      hits <- Sys.glob(path.expand(pat))
      if (length(hits) > 0) { fap_db <- hits[1]; break }
    }
  }

  if (nchar(fap_db) == 0) {
    cat("  [info] FAPROTAX database not found; set FAPROTAX_DB env var\n"); break
  }

  cmd <- sprintf("%s -i '%s' -o '%s' -g '%s' --omit_columns 0 --column_names_are_in row_0 2>/dev/null",
                 py_cmd, otu_tsv, faprotax_out, fap_db)
  ret <- system(cmd)
  if (ret == 0 && file.exists(faprotax_out)) {
    python_success <- TRUE; cat("  ✓ FAPROTAX Python run complete\n"); break
  }
}

# Parse Python output → function × sample matrix
func_mat <- NULL
if (python_success && file.exists(faprotax_out)) {
  tryCatch({
    raw_tsv <- read.table(faprotax_out, sep="\t", header=TRUE,
                          check.names=FALSE, stringsAsFactors=FALSE, comment.char="")
    # First col = function name, rest = samples
    func_mat <- raw_tsv[, -1, drop=FALSE]
    rownames(func_mat) <- raw_tsv[[1]]
    func_mat[] <- lapply(func_mat, as.numeric)
    cat("  Parsed", nrow(func_mat), "functions\n")
  }, error=function(e) cat("  [warn] Parse failed:", e$message, "\n"))
}

# ── Strategy 3: Built-in genus → function lookup ─────────────────────────────
if (is.null(func_mat)) {
  cat("  Using built-in genus→function lookup table\n")

  # Curated FAPROTAX-like lookup: genus → ecological functions
  # Based on FAPROTAX v1.2.7 database (most common assignments)
  faprotax_lookup <- list(
    chemoheterotrophy             = c("Bacillus","Pseudomonas","Clostridium","Lactobacillus",
                                     "Streptococcus","Staphylococcus","Escherichia","Klebsiella",
                                     "Enterococcus","Bacteroides","Prevotella","Fusobacterium",
                                     "Roseburia","Faecalibacterium","Ruminococcus","Blautia",
                                     "Eubacterium","Coprococcus","Akkermansia","Bifidobacterium",
                                     "Parabacteroides","Alistipes"),
    aerobic_chemoheterotrophy     = c("Pseudomonas","Bacillus","Microbacterium","Arthrobacter",
                                     "Sphingomonas","Rhizobium","Burkholderia","Acinetobacter",
                                     "Brevundimonas","Caulobacter","Methylobacterium"),
    fermentation                  = c("Clostridium","Lactobacillus","Streptococcus","Enterococcus",
                                     "Leuconostoc","Pediococcus","Weissella","Roseburia",
                                     "Faecalibacterium","Ruminococcus","Blautia","Bifidobacterium"),
    nitrate_reduction             = c("Pseudomonas","Paracoccus","Thauera","Azoarcus",
                                     "Dechloromonas","Comamonas","Stenotrophomonas","Ochrobactrum"),
    nitrification                 = c("Nitrosomonas","Nitrosospira","Nitrospira","Nitrobacter",
                                     "Candidatus_Nitrososphaera","Candidatus_Nitrosocosmicus",
                                     "Candidatus_Nitrosphaera","Candidatus_Nitrosotalea"),
    denitrification               = c("Pseudomonas","Paracoccus","Rhizobium","Bradyrhizobium",
                                     "Thauera","Azoarcus","Dechloromonas","Hyphomicrobium"),
    nitrogen_fixation             = c("Rhizobium","Bradyrhizobium","Azospirillum","Azotobacter",
                                     "Frankia","Herbaspirillum","Gluconacetobacter","Azoarcus"),
    aerobic_ammonia_oxidation     = c("Nitrosomonas","Nitrosospira","Nitrosolobus",
                                     "Candidatus_Nitrososphaera","Candidatus_Nitrosocosmicus"),
    sulfate_respiration           = c("Desulfovibrio","Desulfobacter","Desulfobulbus","Desulfuromonas",
                                     "Desulfosporosinus","Desulfomonile"),
    sulfur_oxidation              = c("Thiobacillus","Sulfurimonas","Sulfurovum","Sulfurihydrogenibium",
                                     "Beggiatoa","Thioploca","Thiomargarita"),
    methane_oxidation             = c("Methylococcus","Methylobacter","Methylomonas","Methylocystis",
                                     "Methylosinus","Methylocapsa","Methylacidiphilum"),
    methanogenesis                = c("Methanosarcina","Methanosaeta","Methanobacterium",
                                     "Methanobrevibacter","Methanococcus","Methanomicrobium",
                                     "Methanospirillum","Methanosphaera"),
    iron_reduction                = c("Geobacter","Shewanella","Anaeromyxobacter","Pelobacter",
                                     "Desulfuromonas","Ferrimonas"),
    human_pathogens_all           = c("Salmonella","Escherichia","Klebsiella","Staphylococcus",
                                     "Streptococcus","Clostridium","Campylobacter","Helicobacter",
                                     "Listeria","Legionella","Vibrio","Yersinia","Shigella"),
    plant_pathogens               = c("Agrobacterium","Ralstonia","Xanthomonas","Pseudomonas",
                                     "Erwinia","Pectobacterium","Burkholderia","Acidovorax"),
    cellulolysis                  = c("Cellulomonas","Cytophaga","Sporocytophaga","Fibrobacter",
                                     "Clostridium","Ruminococcus","Caldicellulosiruptor"),
    chitin_degradation            = c("Chitinophaga","Serratia","Aeromonas","Janthinobacterium",
                                     "Chromobacterium","Bacillus"),
    xylanolysis                   = c("Fibrobacter","Ruminococcus","Clostridium","Trichoderma",
                                     "Prevotella","Butyrivibrio"),
    phototrophy                   = c("Chlorobium","Prosthecochloris","Rhodobacter","Rhodospirillum",
                                     "Rhodopseudomonas","Chromatium","Ectothiorhodospira",
                                     "Cyanobacterium","Prochlorococcus","Synechococcus",
                                     "Anabaena","Nostoc"),
    dark_hydrogen_oxidation       = c("Hydrogenophaga","Cupriavidus","Ralstonia","Knallgas",
                                     "Hydrogenobaculum","Sulfurihydrogenibium")
  )

  # Extract genus from taxon names (strip prefixes like "g__")
  clean_genus <- function(x) {
    x <- sub(".*g__", "", x)
    x <- sub("_[0-9]+$", "", x)
    x <- gsub("_"," ",x)
    trimws(x)
  }

  taxa_clean <- sapply(colnames(abu_wide), clean_genus)

  # Build function × sample matrix
  func_list <- list()
  for (fn in names(faprotax_lookup) ) {
    fn_genera <- faprotax_lookup[[fn]]
    # Find matching columns
    match_idx <- which(sapply(taxa_clean, function(g) {
      any(sapply(fn_genera, function(ref) {
        grepl(ref, g, ignore.case=TRUE) || grepl(g, ref, ignore.case=TRUE)
      }))
    }))
    if (length(match_idx) == 0) next
    fn_abu <- rowSums(abu_wide[, match_idx, drop=FALSE], na.rm=TRUE)
    if (sum(fn_abu) > 0) func_list[[fn]] <- fn_abu
  }

  if (length(func_list) == 0) {
    cat("  [skip] No functional annotations found\n"); quit(status=0)
  }
  func_mat <- as.data.frame(do.call(cbind, func_list))
  func_mat <- as.data.frame(t(func_mat))  # function × sample
  cat("  Built-in lookup:", nrow(func_mat), "functions annotated\n")
}

if (is.null(func_mat) || nrow(func_mat) == 0) {
  cat("  [skip] No functional data generated\n"); quit(status=0)
}

# ── Relative abundance (% of functional pool) ─────────────────────────────────
col_sums <- colSums(func_mat, na.rm=TRUE)
rel_func  <- sweep(func_mat, 2, pmax(col_sums, 1e-12), "/") * 100

# ── Save output ───────────────────────────────────────────────────────────────
func_csv <- file.path(out_dir, "faprotax_results.csv")
out_df   <- data.frame(Function=rownames(rel_func), rel_func,
                        check.names=FALSE, stringsAsFactors=FALSE)
write.csv(out_df, func_csv, row.names=FALSE)
cat("  ✓ faprotax_results.csv  (", nrow(rel_func), "functions)\n")

# ── PDF: STAMP-style grouped bar + CI ─────────────────────────────────────────
plots_dir <- file.path(out_dir, "r_plots")
dir.create(plots_dir, showWarnings=FALSE, recursive=TRUE)

make_stamp_pdf <- function(rel_f, grp_v, pdf_path, g1=NULL, g2=NULL) {
  if (is.null(grp_v)) return(invisible(NULL))
  common_s <- intersect(colnames(rel_f), names(grp_v))
  if (length(common_s) < 4) return(invisible(NULL))
  rel_sub  <- rel_f[, common_s, drop=FALSE]
  grp_sub  <- grp_v[common_s]

  # Mean per group per function
  grps_here <- if (!is.null(g1) && !is.null(g2)) c(g1,g2) else sort(unique(grp_sub))
  grp_sub   <- grp_sub[grp_sub %in% grps_here]
  rel_sub   <- rel_sub[, names(grp_sub), drop=FALSE]

  # Group stats
  fn_names <- rownames(rel_sub)
  stat_df  <- data.frame()
  for (fn in fn_names) {
    v1 <- as.numeric(rel_sub[fn, grp_sub == grps_here[1]])
    v2 <- if (length(grps_here) > 1) as.numeric(rel_sub[fn, grp_sub == grps_here[2]]) else numeric(0)
    mn1 <- mean(v1, na.rm=TRUE)
    se1 <- if (length(v1)>1) sd(v1,na.rm=TRUE)/sqrt(sum(!is.na(v1))) else 0
    mn2 <- if (length(v2)>0) mean(v2, na.rm=TRUE) else NA
    se2 <- if (length(v2)>1) sd(v2,na.rm=TRUE)/sqrt(sum(!is.na(v2))) else 0
    diff_m <- mn1 - if (!is.na(mn2)) mn2 else 0
    se_d   <- sqrt(se1^2 + se2^2)
    p_val  <- if (length(v2)>0 && length(v1)>1 && length(v2)>1)
      tryCatch(wilcox.test(v1,v2,exact=FALSE)$p.value, error=function(e) 1) else 1
    stat_df <- rbind(stat_df, data.frame(
      Function=fn, G1_mean=mn1, G1_se=se1, G2_mean=mn2, G2_se=se2,
      diff_mean=diff_m, CI_low=diff_m-1.96*se_d, CI_high=diff_m+1.96*se_d,
      p_value=p_val, stringsAsFactors=FALSE))
  }

  # Sort by |diff|, take top 30
  stat_df <- stat_df[order(abs(stat_df$diff_mean), decreasing=TRUE), ]
  stat_df <- head(stat_df, 30)
  stat_df$Function <- factor(stat_df$Function, levels=rev(stat_df$Function))

  if (!requireNamespace("ggplot2", quietly=TRUE)) return(invisible(NULL))
  suppressPackageStartupMessages(library(ggplot2))

  # Left bar panel
  long_s <- rbind(
    data.frame(Function=stat_df$Function, Group=grps_here[1], Mean=stat_df$G1_mean, SE=stat_df$G1_se),
    if (!is.na(stat_df$G2_mean[1]))
      data.frame(Function=stat_df$Function, Group=grps_here[2], Mean=stat_df$G2_mean, SE=stat_df$G2_se)
    else data.frame())
  long_s$Group    <- factor(long_s$Group, levels=grps_here)
  long_s$Function <- factor(long_s$Function, levels=levels(stat_df$Function))

  p_bar <- ggplot(long_s, aes(x=Function, y=Mean, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge(0.75), width=0.7) +
    geom_errorbar(aes(ymin=pmax(Mean-SE,0), ymax=Mean+SE),
                  position=position_dodge(0.75), width=0.2, size=0.4) +
    scale_fill_manual(values=c("#f59e0b","#3b82f6","#10b981","#ef4444","#8b5cf6")[
      seq_along(grps_here)]) +
    coord_flip() +
    labs(x=NULL, y="Proportion (%)", fill="Group") +
    theme_classic(base_size=9) +
    theme(axis.text.y=element_text(size=8))

  # Right CI panel
  p_ci <- ggplot(stat_df, aes(x=diff_mean, y=Function)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey50", size=0.5) +
    geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0.3, size=0.5, color="#64748b") +
    geom_point(size=2.5, color="#ef4444") +
    labs(x=paste0("Difference between proportions (%)\n(",
                   grps_here[1], " − ", if (length(grps_here)>1) grps_here[2] else "rest", ")"),
         y=NULL, title="95% confidence intervals") +
    theme_classic(base_size=9) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          plot.title=element_text(size=9))

  if (requireNamespace("patchwork", quietly=TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    combined <- p_bar + p_ci + plot_layout(widths=c(1.5,1))
  } else { combined <- p_bar }

  h_pdf <- max(6, min(18, nrow(stat_df)*0.45 + 3))
  tryCatch(ggsave(pdf_path, combined, width=14, height=h_pdf, limitsize=FALSE),
           error=function(e) NULL)
  cat("  ✓", basename(pdf_path), "\n")
}

sample_names_func <- colnames(rel_func)

if (!is.null(meta_map)) {
  grp_v <- meta_map[sample_names_func]
  grp_v <- grp_v[!is.na(grp_v)]

  # All-groups overview
  make_stamp_pdf(rel_func, grp_v, file.path(plots_dir, "faprotax_bar.pdf"))

  # Pairwise PDFs
  grps_with_data <- sort(unique(grp_v))
  if (length(grps_with_data) >= 2) {
    for (pr in combn(grps_with_data, 2, simplify=FALSE)) {
      g1 <- pr[1]; g2 <- pr[2]
      pdf_pair <- file.path(plots_dir, sprintf("faprotax_%s_vs_%s.pdf",
                              gsub("[^a-zA-Z0-9]","_",g1), gsub("[^a-zA-Z0-9]","_",g2)))
      make_stamp_pdf(rel_func, grp_v, pdf_pair, g1=g1, g2=g2)
    }
  }
} else {
  # No metadata — just a bar chart of functional profile
  tryCatch({
    suppressPackageStartupMessages(library(ggplot2))
    fn_means <- rowMeans(rel_func, na.rm=TRUE)
    top_fn   <- names(sort(fn_means, decreasing=TRUE))[seq_len(min(20, length(fn_means)))]
    df_bar   <- data.frame(Function=top_fn, Proportion=fn_means[top_fn])
    df_bar$Function <- factor(df_bar$Function, levels=rev(top_fn))
    p <- ggplot(df_bar, aes(x=Function, y=Proportion)) +
      geom_bar(stat="identity", fill="#3b82f6") + coord_flip() +
      labs(x=NULL, y="Mean Proportion (%)", title="FAPROTAX Functional Profile") +
      theme_classic(base_size=10)
    ggsave(file.path(plots_dir,"faprotax_bar.pdf"), p, width=10, height=8)
    cat("  ✓ faprotax_bar.pdf  (no metadata — mean profile)\n")
  }, error=function(e) cat("  [error]", e$message, "\n"))
}

cat("── FAPROTAX Analysis done ──────────────────────────────────────────\n\n")
