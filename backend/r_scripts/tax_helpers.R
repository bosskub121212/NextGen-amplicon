# =============================================================================
#  tax_helpers.R — shared taxonomy-assignment helpers
#
#  source()d by dada2_pipeline.R and pacbio_pipeline.R. These three problems are
#  properties of the REFERENCE FILE, not of any one pipeline, so both callers hit
#  them identically — and when the fixes lived only in dada2_pipeline.R, the
#  PacBio pipeline silently kept all three. Keeping them here means a fix lands
#  everywhere at once.
#
#    1. Rank depth   — SILVA 144 inserts a Kingdom rank between Domain and
#                      Phylum (LPSN scheme). With DADA2's default 7 rank names
#                      every rank below Domain is labelled one level too high:
#                      families end up in the "Genus" column and the real genus
#                      falls off the end. Nothing errors — the table just means
#                      something other than what its headers say.
#    2. Species file — picking the first match in the folder pairs whichever
#                      file the directory listing returns first with the
#                      trainset, silently mixing reference versions.
#    3. Duplicate    — addSpecies() APPENDS a Species column. When the trainset
#       Species        already carried one, two columns end up named "Species",
#                      and every later which(colnames(tax)=="Species") returns
#                      two indices, which breaks downstream indexing with
#                      "(subscript) logical subscript too long".
# =============================================================================

# Maximum lineage depth over many headers — NOT the first one. Truncated
# lineages are common (SILVA 144's own leading records stop at family), so a
# single header under-reports the reference's real depth.
tax_detect_db_depth <- function(p, n_headers = 3000) {
  tryCatch({
    con <- if (grepl("\\.gz$", p)) gzfile(p, "rt") else file(p, "rt")
    on.exit(close(con), add = TRUE)
    seen <- 0L; best <- 0L
    while (seen < n_headers) {
      chunk <- readLines(con, n = 2000, warn = FALSE)
      if (length(chunk) == 0) break
      hdrs <- chunk[startsWith(chunk, ">")]
      if (length(hdrs) > 0) {
        d <- vapply(hdrs, function(h) {
          parts <- trimws(strsplit(sub("^>", "", h), ";")[[1]])
          sum(nzchar(parts))
        }, integer(1))
        best <- max(best, max(d)); seen <- seen + length(hdrs)
      }
    }
    if (best > 0L) best else NA_integer_
  }, error = function(e) NA_integer_)
}


# Rank names matching the reference's actual depth. Returns NULL to mean
# "use the DADA2 default", so callers can pass it straight through.
tax_pick_levels <- function(db_path, verbose = TRUE) {
  depth <- tax_detect_db_depth(db_path)
  if (is.na(depth)) {
    if (verbose) cat("  [warn] Could not read rank depth from the reference — using DADA2 defaults\n")
    return(NULL)
  }
  is_to_species <- grepl("tospecies", basename(db_path), ignore.case = TRUE)
  # Four shapes in circulation:
  #   6  Kingdom..Genus              SILVA 138.x toGenus
  #   7  Kingdom..Species            SILVA 138.x toSpecies
  #   7  Domain..Genus               SILVA 144  toGenus    (extra Kingdom rank)
  #   8  Domain..Species             SILVA 144  toSpecies
  # Depth alone cannot separate the two 7s, so the filename breaks the tie.
  has_domain <- depth >= 8 || (depth == 7 && !is_to_species)
  lv <- if (depth >= 8) {
    c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:min(depth, 8)]
  } else if (depth == 7 && !is_to_species) {
    c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")
  } else if (depth == 6) {
    c("Kingdom","Phylum","Class","Order","Family","Genus")
  } else {
    c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:min(depth, 7)]
  }
  if (verbose) {
    cat(sprintf("  Reference has %d rank level(s) -> %s\n", depth, paste(lv, collapse = ", ")))
    if (has_domain)
      cat("  (extra Kingdom rank detected — SILVA 144 / LPSN naming)\n")
  }
  lv
}


# The assignSpecies file whose version matches the trainset. Returns NULL when
# the folder has none.
tax_pick_species_file <- function(db_path, verbose = TRUE) {
  cands <- list.files(dirname(db_path),
                      pattern = "(assignspecies|species_assignment|_species)\\.?.*\\.fa(sta)?(\\.gz)?$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(cands) == 0) return(NULL)

  ver <- function(f) {
    n <- tolower(basename(f))
    m <- regmatches(n, regexpr("v[0-9]+(\\.[0-9]+)*", n))        # SILVA: v138.2
    if (length(m) > 0) return(sub("^v", "", m[1]))
    m <- regmatches(n, regexpr("_r[0-9]+", n))                    # GTDB: _r220
    if (length(m) > 0) return(sub("^_r", "r", m[1]))
    NA_character_
  }
  train_ver <- ver(db_path)
  cand_vers <- vapply(cands, ver, character(1))

  if (!is.na(train_ver)) {
    exact <- cands[!is.na(cand_vers) & cand_vers == train_ver]
    if (length(exact) > 0) {
      if (verbose && length(cands) > 1)
        cat(sprintf("  Version-matched species file (trainset v%s) out of %d candidates\n",
                    train_ver, length(cands)))
      return(exact[1])
    }
  }
  if (verbose && !is.na(train_ver) && length(cands) > 1) {
    cat(sprintf("  [WARN] No species file matching trainset version '%s'.\n", train_ver))
    cat(sprintf("         Using '%s' — species names come from a different reference build.\n",
                basename(cands[1])))
    cat("         Available:", paste(basename(cands), collapse = ", "), "\n")
  }
  cands[1]
}


# Collapse duplicate "Species" columns, keeping the one with real assignments.
tax_merge_dup_species <- function(tax, verbose = TRUE) {
  if (is.null(tax)) return(tax)
  dup <- which(colnames(tax) == "Species")
  if (length(dup) <= 1) return(tax)
  filled <- vapply(dup, function(i) sum(!is.na(tax[, i])), integer(1))
  keep   <- dup[which.max(filled)]
  drop   <- setdiff(dup, keep)
  if (verbose)
    cat(sprintf("  Merged %d duplicate 'Species' column(s); kept the one with %d assignment(s)\n",
                length(drop), max(filled)))
  tax[, -drop, drop = FALSE]
}


# ---------------------------------------------------------------------------
#  tax_add_species() — addSpecies() that survives SILVA 144's genus suffixes
#
#  addSpecies() does NOT assign on sequence identity alone. After the exact
#  100% match it calls DADA2's internal matchGenera() to check the binomial's
#  genus against the genus already in the table, and drops the species when
#  they disagree:
#
#      matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
#        if (is.na(gen.tax) || is.na(gen.binom)) return(FALSE)
#        if ((gen.tax == gen.binom) ||
#            grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
#            grepl(paste0("[ _", split.glyph, "]", gen.binom, "$"), gen.tax))
#          TRUE else FALSE
#      }
#
#  Only space, underscore and "/" count as separators. SILVA 144 marks an
#  uncertain genus by appending "--other" (Chryseobacterium--other,
#  Staphylococcus--other, Coprococcus--other), so:
#
#      matchGenera("Chryseobacterium--other", "Chryseobacterium")
#        "==" ? no.  "^Chryseobacterium[ _/]" ? no — next char is "-".
#        -> FALSE, and Chryseobacterium hominis is discarded despite a
#           100% sequence match.
#
#  On a real V3-V4 run this silently removed every species sitting on one of
#  20 such ASVs (8.3% of reads). Nothing errors; the species column is simply
#  emptier than it should be, which reads as "this sample has no species" —
#  a database problem that is not one.
#
#  So: strip the suffix, match, then PUT THE ORIGINAL NAME BACK. The "--other"
#  is real information — it is SILVA saying it could not resolve the genus —
#  and dropping it permanently would turn a hedged call into a confident one.
# ---------------------------------------------------------------------------
tax_add_species <- function(tax, sp_file, verbose = TRUE) {
  if (is.null(tax) || is.null(sp_file)) return(tax)

  # Suffixes seen in the wild. Anchored, so a genus that legitimately contains
  # a hyphen (e.g. "Candidatus Hepatoplasma") is left alone.
  suffix_re <- "--+(other|unclassified|uncultured|unknown)$"

  gcol <- which(colnames(tax) == "Genus")
  genus_orig <- NULL
  if (length(gcol) == 1) {
    g <- tax[, gcol]
    hit <- !is.na(g) & grepl(suffix_re, g, ignore.case = TRUE)
    if (any(hit)) {
      genus_orig   <- g
      tax[, gcol]  <- sub(suffix_re, "", g, ignore.case = TRUE)
      if (verbose)
        cat(sprintf("  Normalized %d genus name(s) carrying a '--other'-style suffix so addSpecies() can match them\n",
                    sum(hit)))
    }
  } else if (length(gcol) > 1 && verbose) {
    cat("  [warn] more than one 'Genus' column — skipping genus normalization\n")
  }

  gc(verbose = FALSE)
  tax <- addSpecies(tax, sp_file)

  # Restore the hedge. addSpecies() appends a column rather than reordering,
  # so the Genus column is still at gcol and still row-aligned.
  if (!is.null(genus_orig) && length(gcol) == 1 &&
      identical(colnames(tax)[gcol], "Genus") && nrow(tax) == length(genus_orig)) {
    tax[, gcol] <- genus_orig
  }

  tax <- tax_merge_dup_species(tax, verbose = verbose)

  # How many of the assignments landed on a normalized genus? Those are the
  # ones plain addSpecies() would have dropped. Counted directly rather than
  # inferred from a before/after delta, which would also sweep in species that
  # were never at risk.
  if (verbose) {
    scol <- which(colnames(tax) == "Species")
    if (length(scol) == 1) {
      assigned <- !is.na(tax[, scol])
      rescued  <- if (!is.null(genus_orig))
        sum(assigned & grepl(suffix_re, genus_orig, ignore.case = TRUE)) else 0L
      cat(sprintf("  Species assigned to %d ASV(s)%s\n", sum(assigned),
                  if (rescued > 0)
                    sprintf("; %d would have been dropped by the genus-name check", rescued) else ""))
    }
  }
  tax
}
