# =============================================================================
#  custom_report_plots.R — the figures of the per-sample customer report
#
#    Rscript custom_report_plots.R --jobs jobs.tsv
#  jobs.tsv: work_dir <TAB> Images_and_Outfiles dir <TAB> sample name [<TAB> unit], one line per sample.
#  unit: ASV (DADA2 runs), taxon (ONT 16S / Emu, estimated reads) or ppm (Emu without counts:
#  no rarefaction). Quality profiles (dada2) are drawn for ASV runs only.
#  work_dir holds what custom_report.py prepared: asv.tsv (SV, label, ranks, count),
#  track.tsv (step, count), krona.tsv, reads.json (raw / filtered FASTQ paths).
#
#  Writes, with the file names of the lab's report:
#    Box_plot.png, Taxonomic profiles/plot_bar_<rank>.png (phyloseq + ggplot2),
#    Heatmap/Heat_map_label_ASV_number_top50ASVs.png (microbiomeutilities, else pheatmap),
#    Krona plot/krona_static.png (still of the Krona rings), Rarefaction curves/
#    Rarefaction_curves.png (vegan), plotQualityProfile[_filtered].pdf (dada2),
#  and in work_dir: alpha.tsv (phyloseq::estimate_richness), plots_info.json (versions).
#  A figure that fails is skipped with a message; the others are still made.
# =============================================================================
suppressPackageStartupMessages({ library(optparse) })
opt <- parse_args(OptionParser(option_list = list(make_option("--jobs", type = "character"))))
say <- function(...) cat("[custom_report]", ..., "\n")
has <- function(p) requireNamespace(p, quietly = TRUE)
ver <- function(p) if (has(p)) as.character(utils::packageVersion(p)) else ""
suppressPackageStartupMessages({ library(ggplot2); library(phyloseq) })
RANKS <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

palette_n <- function(n) {
  base <- c("#7E57C2", "#26A69A", "#EF5350", "#42A5F5", "#FFA726", "#8D6E63", "#EC407A", "#66BB6A",
            "#5C6BC0", "#D4E157", "#26C6DA", "#AB47BC", "#FF7043", "#9CCC65", "#78909C", "#FFCA28",
            "#29B6F6", "#EF9A9A", "#80CBC4", "#B39DDB")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dynamic")
}

jobs <- read.delim(opt$jobs, header = FALSE, stringsAsFactors = FALSE, quote = "")
heat_tool <- "pheatmap"
for (k in seq_len(nrow(jobs))) {
  wd <- jobs[k, 1]; img <- jobs[k, 2]; s <- jobs[k, 3]
  unit <- if (ncol(jobs) >= 4 && !is.na(jobs[k, 4]) && nzchar(jobs[k, 4])) jobs[k, 4] else "ASV"
  say("sample", s)
  asv <- read.delim(file.path(wd, "asv.tsv"), stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  trk <- read.delim(file.path(wd, "track.tsv"), stringsAsFactors = FALSE)
  reads <- tryCatch(jsonlite::fromJSON(file.path(wd, "reads.json")), error = function(e) list())

  # ── Box plot of the read counts per step ───────────────────────────────────
  tryCatch({
    trk$step <- factor(trk$step, levels = trk$step)
    g <- ggplot(trk, aes(step, count)) + geom_boxplot() + geom_point(size = 2) +
      scale_y_continuous(labels = scales::comma) + labs(x = NULL, y = "Number") + theme_gray(base_size = 11)
    ggsave(file.path(img, "Box_plot.png"), g, width = 8, height = 5, dpi = 150)
  }, error = function(e) say("box plot failed:", conditionMessage(e)))

  have_reads <- nrow(asv) > 0 && sum(asv$count) > 0
  # rows that are not a taxon (Emu's "unassigned"): in the profiles, not in richness / rarefaction
  ina <- if ("in_alpha" %in% names(asv)) asv$in_alpha == 1 else rep(TRUE, nrow(asv))
  asv_a <- asv[ina, , drop = FALSE]
  ps <- NULL
  if (have_reads) {
    otu <- matrix(asv$count, ncol = 1, dimnames = list(asv$SV, s))
    tt <- as.matrix(asv[, RANKS, drop = FALSE]); rownames(tt) <- asv$SV
    sd <- data.frame(Sample_name = s, row.names = s)
    ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tt), sample_data(sd))
  }

  # ── Bar charts per rank (phyloseq tax_glom + psmelt, ggplot2) ─────────────
  for (r in RANKS) {
    tryCatch({
      f <- file.path(img, "Taxonomic profiles", sprintf("plot_bar_%s.png", tolower(r)))
      if (is.null(ps)) { png(f, 600, 400); plot.new(); text(.5, .5, "No reads"); dev.off(); next }
      gl <- tax_glom(ps, taxrank = r, NArm = FALSE)
      rel <- transform_sample_counts(gl, function(x) x / sum(x))
      m <- psmelt(rel)
      m[[r]] <- factor(m[[r]], levels = unique(m[[r]][order(-m$Abundance)]))
      n <- length(unique(m[[r]]))
      g <- ggplot(m, aes(x = Sample, y = Abundance, fill = .data[[r]])) +
        geom_col(width = 0.35, colour = NA) +
        scale_fill_manual(values = palette_n(n)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.01)) +
        labs(title = sprintf("%d Taxa (%s)", n, r), x = NULL, y = "Relative Abundance", fill = r) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.key.size = grid::unit(9, "pt"),
              legend.text = element_text(size = if (n > 30) 6 else 8, face = if (r %in% c("Genus", "Species")) "italic" else "plain"),
              plot.title = element_text(size = 10)) +
        guides(fill = guide_legend(ncol = if (n > 30) 2 else 1))
      ggsave(f, g, width = if (n > 30) 7 else 5.5, height = 6, dpi = 150)
    }, error = function(e) say("bar", r, "failed:", conditionMessage(e)))
  }

  # ── Top 50 ASV heatmap (microbiomeutilities, else pheatmap) ───────────────
  tryCatch({
    f <- file.path(img, "Heatmap", "Heat_map_label_ASV_number_top50ASVs.png")
    if (!is.null(ps)) {
      top <- head(asv[order(-asv$count), ], 50)
      done <- FALSE
      if (has("microbiomeutilities") && nrow(top) >= 2) {
        done <- tryCatch({
          pt <- prune_taxa(top$SV, ps)
          taxa_names(pt) <- top$label[match(taxa_names(pt), top$SV)]
          png(f, width = 1600, height = 300 + 42 * nrow(top), res = 150)
          microbiomeutilities::plot_taxa_heatmap(pt, subset.top = 50, VariableA = "Sample_name",
                                                 transformation = "log10", cluster_cols = FALSE,
                                                 heatcolors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100))
          dev.off(); heat_tool <- paste("microbiomeutilities", ver("microbiomeutilities")); TRUE
        }, error = function(e) { try(dev.off(), silent = TRUE); say("microbiomeutilities heatmap failed:", conditionMessage(e)); FALSE })
      }
      if (!done) {
        mat <- matrix(log10(top$count + 1), ncol = 1, dimnames = list(top$label, s))
        ann <- data.frame(Sample_name = s, row.names = s)
        png(f, width = 1600, height = 300 + 42 * nrow(top), res = 150)
        pheatmap::pheatmap(mat, cluster_rows = nrow(mat) > 2, cluster_cols = FALSE, annotation_col = ann,
                           color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100),
                           fontsize_row = 11, cellwidth = 40, cellheight = 20, border_color = "grey60",
                           main = "log10 reads")
        dev.off(); heat_tool <- paste("pheatmap", ver("pheatmap"))
      }
    } else { png(f, 600, 300); plot.new(); text(.5, .5, "No reads"); dev.off() }
  }, error = function(e) { try(dev.off(), silent = TRUE); say("heatmap failed:", conditionMessage(e)) })

  # ── Krona still: taxonomy rings, Phylum inside to Species outside ─────────
  tryCatch({
    f <- file.path(img, "Krona plot", "krona_static.png")
    if (have_reads) {
      a <- asv[order(asv$Phylum, asv$Class, asv$Order, asv$Family, asv$Genus, asv$Species), ]
      tot <- sum(a$count); rings <- list()
      for (i in seq_along(RANKS)) {
        key <- do.call(paste, c(a[, RANKS[1:i], drop = FALSE], sep = "|"))
        agg <- aggregate(a$count, list(key = factor(key, levels = unique(key))), sum)
        agg$ymax <- cumsum(agg$x) / tot; agg$ymin <- c(0, head(agg$ymax, -1))
        agg$name <- sub(".*\\|", "", as.character(agg$key)); agg$ring <- i
        rings[[i]] <- agg
      }
      df <- do.call(rbind, rings)
      df$lab <- ifelse(df$ymax - df$ymin > 0.04, df$name, "")
      g <- ggplot(df) + geom_rect(aes(xmin = ring, xmax = ring + 0.95, ymin = ymin, ymax = ymax, fill = name),
                                  colour = "white", linewidth = 0.2, show.legend = FALSE) +
        geom_text(aes(x = ring + 0.47, y = (ymin + ymax) / 2, label = lab), size = 2.1) +
        coord_polar(theta = "y") + xlim(0, length(RANKS) + 1) + theme_void() +
        scale_fill_manual(values = palette_n(length(unique(df$name)))) +
        labs(title = s) + theme(plot.title = element_text(hjust = 0.5))
      ggsave(f, g, width = 7, height = 7, dpi = 150)
    } else { png(f, 600, 300); plot.new(); text(.5, .5, "No reads"); dev.off() }
  }, error = function(e) say("Krona still failed:", conditionMessage(e)))

  # ── Rarefaction (vegan) ───────────────────────────────────────────────────
  tryCatch({
    f <- file.path(img, "Rarefaction curves", "Rarefaction_curves.png")
    png(f, width = 1400, height = 900, res = 150)
    if (unit == "ppm") { plot.new(); text(.5, .5, "Counts are ppm, not reads: no rarefaction") } else
    if (nrow(asv_a) > 0 && sum(asv_a$count) > 1) {
      m <- matrix(asv_a$count, nrow = 1, dimnames = list(s, asv_a$SV))
      vegan::rarecurve(m, step = max(1, floor(sum(asv_a$count) / 200)), label = TRUE, col = "#1E6E8C",
                       xlab = "Sample size (reads)", ylab = if (unit == "ASV") "ASVs" else "Taxa")
    } else { plot.new(); text(.5, .5, "No reads") }
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); say("rarefaction failed:", conditionMessage(e)) })

  # ── Alpha diversity (phyloseq::estimate_richness) ─────────────────────────
  tryCatch({
    if (!is.null(ps) && nrow(asv_a) > 0 && sum(asv_a$count) > 0) {
      ps_a <- prune_taxa(asv_a$SV, ps)
      er <- estimate_richness(ps_a, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
      out <- data.frame(Sample_ID = s, Sample_name = s, Observed = er$Observed, Chao1 = er$Chao1,
                        se.chao1 = er$se.chao1, Shannon = er$Shannon, Simpson = er$Simpson, check.names = FALSE)
      write.table(out, file.path(wd, "alpha.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }, error = function(e) say("alpha failed:", conditionMessage(e)))

  # ── Quality profiles (dada2) ──────────────────────────────────────────────
  if (has("dada2") && unit == "ASV") {
    for (pair in list(list("raw", "plotQualityProfile.pdf"), list("filtered", "plotQualityProfile_filtered.pdf"))) {
      fs <- reads[[pair[[1]]]]
      fs <- fs[file.exists(fs)]
      if (length(fs)) tryCatch(ggsave(file.path(img, pair[[2]]), dada2::plotQualityProfile(fs), width = 10, height = 5),
                               error = function(e) say("quality profile failed:", conditionMessage(e)))
    }
  }
  jsonlite::write_json(list(phyloseq = paste("phyloseq", ver("phyloseq")), ggplot2 = ver("ggplot2"),
                            vegan = ver("vegan"), dada2 = ver("dada2"), heatmap_tool = heat_tool,
                            venneuler = ver("venneuler"), eulerr = ver("eulerr")),
                       file.path(wd, "plots_info.json"), auto_unbox = TRUE)
}
