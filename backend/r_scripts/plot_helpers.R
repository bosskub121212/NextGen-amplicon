# =============================================================================
#  plot_helpers.R — one taxonomy-bar spec shared by every renderer
#
#  The same stacked bar is drawn three times in this app: by dada2_pipeline.R
#  (taxonomy_<rank>.pdf), by viz_pipeline.R (04_taxonomy_<rank>.pdf), and by the
#  browser with Plotly (PreviewPage.tsx). They were drawing it three different
#  ways, and the differences all showed up at once on the same sample:
#
#    * viz_pipeline.R loaded the taxonomy with tax_mat[is.na(tax_mat)] <- "",
#      then replaced only NA with "Unclassified" — so the 45.7% unassigned group
#      kept an EMPTY label. It rendered as a blank legend row in the PDF and as
#      "Unknown" in the browser, which is the JS fallback for an empty name.
#    * ggplot2 puts the FIRST factor level at the TOP of a stack; Plotly puts the
#      first trace at the BOTTOM. Same data, same order, mirrored bars.
#    * viz_pipeline used scales::hue_pal(), dada2_pipeline used its own 50-colour
#      list, and the browser used a third one — so a genus was a different colour
#      in each of the three views of the same run.
#
#  Anything that draws a taxonomy bar sources this file and uses these three
#  functions, so the PDF a customer receives and the chart on screen agree.
#
#  THE SPEC (the browser already did this; the PDFs now match it)
#    label   unassigned is always "Unclassified", never "", NA or "Unknown"
#    order   descending total abundance, with "Other" pinned last
#    stack   first item at the BOTTOM, so the largest group is the base
#    legend  reversed, so reading the legend top-to-bottom follows the bar
#            top-to-bottom
#    colour  by position in that order, from one shared list; "Unclassified"
#            and "Other" are always grey and never consume a colour slot
# =============================================================================

TAXA_UNASSIGNED <- "Unclassified"
TAXA_OTHER      <- "Other"

# Kept character-for-character in sync with DEFAULT_COLORS in
# frontend/src/pages/PreviewPage.tsx. Editing one without the other puts the
# colours back out of step, which is the bug this file exists to prevent.
TAXA_PALETTE <- c(
  "#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6",
  "#06b6d4","#f97316","#84cc16","#ec4899","#14b8a6",
  "#6366f1","#a855f7","#22d3ee","#fb923c","#a3e635",
  "#f472b6","#2dd4bf","#818cf8","#fb7185","#fbbf24",
  "#0ea5e9","#d946ef","#65a30d","#dc2626","#7c3aed",
  "#0891b2","#be123c","#4d7c0f","#c2410c","#1d4ed8",
  "#9333ea","#047857","#b45309","#0369a1","#a21caf",
  "#15803d","#b91c1c","#4338ca","#0f766e","#c026d3"
)
TAXA_COL_UNASSIGNED <- "#9ca3af"
TAXA_COL_OTHER      <- "#d1d5db"


# Every spelling of "we could not assign this" collapses to one label. Covers
# NA, empty and whitespace-only strings, the literal "NA" that survives a CSV
# round-trip, the JS-side "Unknown", and QIIME-style "g__" / "g__unassigned"
# prefixes that otherwise appear as their own taxon.
tax_clean_labels <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(sub("^[a-z]__", "", x))
  blank <- x == "" |
    grepl("^(na|nan|null|unknown|unassigned|unclassified|uncultured|undetermined)$",
          x, ignore.case = TRUE)
  x[blank] <- TAXA_UNASSIGNED
  x
}


# Factor levels in draw order: biggest first, "Other" last. `abund` is any
# numeric vector parallel to `labels` (counts or percentages both work).
tax_order_levels <- function(labels, abund) {
  tot <- tapply(as.numeric(abund), labels, sum, na.rm = TRUE)
  tot <- tot[order(tot, decreasing = TRUE)]
  lv  <- names(tot)
  if (TAXA_OTHER %in% lv) lv <- c(setdiff(lv, TAXA_OTHER), TAXA_OTHER)
  lv
}


# Named colour vector for those levels. The two grey buckets are pinned and
# skipped in the rotation, so adding or removing an "Other" bucket does not
# shift every other taxon's colour.
tax_colors <- function(levels) {
  out <- character(length(levels))
  names(out) <- levels
  i <- 0L
  for (k in seq_along(levels)) {
    lv <- levels[k]
    if (identical(lv, TAXA_UNASSIGNED))      out[k] <- TAXA_COL_UNASSIGNED
    else if (identical(lv, TAXA_OTHER))      out[k] <- TAXA_COL_OTHER
    else { i <- i + 1L; out[k] <- TAXA_PALETTE[((i - 1L) %% length(TAXA_PALETTE)) + 1L] }
  }
  out
}
