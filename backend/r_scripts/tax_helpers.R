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


# What KIND of reference file this is, and how deep it goes — read from the
# file's own headers, never from its name.
#
# Two things used to be decided by filename and both got them wrong:
#
#   * "is this a species-assignment file?"  matched assignSpecies|toSpecies,
#     so silva_nr99_v138.2_toSpecies_trainset.fa.gz — a perfectly good 7-rank
#     trainset — was treated as a species file and skipped, while a real
#     assignSpecies file handed in under any other name sailed straight into
#     assignTaxonomy() and produced a run with no taxonomy at all.
#   * "does depth 7 mean Domain..Genus or Kingdom..Species?"  matched
#     "tospecies" in the name, so a renamed or third-party file shifted every
#     rank by one with no error.
#
# Both are answered by the headers. A trainset header is a ';'-separated
# lineage; an assignSpecies header is "ACCESSION Genus species" with no ';'.
# And a trainset that carries a Species rank writes a binomial whose first word
# repeats the genus above it ("...;Bacillus;Bacillus subtilis;"), which the
# majority test below detects without trusting any one record.
tax_inspect_reference <- function(p, n_headers = 2000) {
  default <- list(depth = NA_integer_, has_species = NA, kind = "unknown", n = 0L)
  tryCatch({
    con <- if (grepl("\\.gz$", p)) gzfile(p, "rt") else file(p, "rt")
    on.exit(close(con), add = TRUE)
    hdrs <- character(0)
    while (length(hdrs) < n_headers) {
      chunk <- readLines(con, n = 2000, warn = FALSE)
      if (length(chunk) == 0) break
      hdrs <- c(hdrs, chunk[startsWith(chunk, ">")])
    }
    if (length(hdrs) == 0) return(default)
    parts <- lapply(hdrs, function(h) {
      x <- trimws(strsplit(sub("^>", "", h), ";")[[1]]); x[nzchar(x)]
    })
    d <- lengths(parts)
    depth <- max(d)
    med   <- stats::median(d)

    if (med <= 2) {
      first  <- vapply(parts, function(x) if (length(x)) x[1] else "", character(1))
      binom  <- mean(grepl("^\\S+\\s+\\S+\\s+\\S+", first))   # "ACC Genus species"
      return(list(depth = depth, has_species = NA,
                  kind = if (binom > 0.5) "species_assignment" else "unknown",
                  n = length(hdrs)))
    }

    full <- parts[d == depth]
    has_sp <- NA
    if (length(full) >= 20) {
      # a species rank is written as "Genus species", "Genus_species" or the bare
      # epithet ("coli") depending on who built the file — any of the three counts
      hit <- vapply(full, function(x) {
        g <- x[length(x) - 1]; sp <- x[length(x)]
        isTRUE(nzchar(g) && (startsWith(sp, paste0(g, " ")) || startsWith(sp, paste0(g, "_")) ||
                             grepl("^[a-z][a-z-]+$", sp)))
      }, logical(1))
      has_sp <- mean(hit) > 0.5
    }
    # SILVA 144 (LPSN) puts a kingdom between domain and phylum, and LPSN
    # kingdom names end in "-ati" (Bacillati, Pseudomonadati, Methanobacteriati).
    # A phylum never does, so the second level settles Domain..Genus vs
    # Kingdom..Species for a 7-level file even when the species test cannot.
    second <- vapply(parts[d >= 3], function(x) x[2], character(1))
    has_kingdom <- if (length(second) >= 20) mean(grepl("ati$", second)) > 0.5 else NA
    # MIDORI2 writes the NCBI taxid after every name ("Arthropoda_6656;...;
    # Musca domestica_7370;") and starts at phylum: 6 levels = Phylum..Species.
    # Read with DADA2's default names every rank lands one column too high
    # (phylum as "Kingdom", species as "Genus").
    last <- vapply(full, function(x) x[length(x)], character(1))
    midori <- length(last) >= 20 && mean(grepl("_[0-9]+$", last)) > 0.8
    list(depth = depth, has_species = has_sp, has_kingdom = has_kingdom,
         midori = midori, kind = "trainset", n = length(hdrs))
  }, error = function(e) default)
}


# Rank names matching the reference's actual depth. Returns NULL to mean
# "use the DADA2 default", so callers can pass it straight through.
tax_pick_levels <- function(db_path, verbose = TRUE) {
  info  <- tax_inspect_reference(db_path)
  depth <- info$depth
  if (is.na(depth)) {
    if (verbose) cat("  [warn] Could not read rank depth from the reference — using DADA2 defaults\n")
    return(NULL)
  }
  # Measured from the headers; the filename is only a last resort when the
  # binomial test had too few full-depth records to decide.
  hk <- if (is.null(info$has_kingdom)) NA else info$has_kingdom
  is_to_species <- if (depth == 7 && !is.na(hk)) !isTRUE(hk)
                   else if (!is.na(info$has_species)) isTRUE(info$has_species)
                   else grepl("tospecies|wspecies", basename(db_path), ignore.case = TRUE)
  # Four shapes in circulation:
  #   6  Kingdom..Genus              SILVA 138.x toGenus
  #   7  Kingdom..Species            SILVA 138.x toSpecies
  #   7  Domain..Genus               SILVA 144  toGenus    (extra Kingdom rank)
  #   8  Domain..Species             SILVA 144  toSpecies
  # Depth alone cannot separate the two 7s; the binomial test above breaks the tie.
  has_domain <- depth >= 8 || (depth == 7 && !is_to_species)
  # PR2 v5 (18S): domain > supergroup > division > subdivision > class > order >
  # family > genus > species. Division is the phylum-level rank, so it is named
  # Phylum here — everything downstream (plots, report) reads Phylum..Species.
  is_pr2 <- depth >= 9 || grepl("pr2", basename(db_path), ignore.case = TRUE)
  is_midori <- isTRUE(info$midori) || grepl("midori", basename(db_path), ignore.case = TRUE)
  lv <- if (is_midori && depth <= 7) {
    # MIDORI2: phylum;class;order;family;genus;species (6), or with kingdom (7)
    c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[(8 - depth):7]
  } else if (is_pr2 && depth >= 8) {
    c("Kingdom","Supergroup","Phylum","Subdivision","Class","Order","Family","Genus","Species")[1:min(depth, 9)]
  } else if (depth >= 8) {
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
    if (is_midori && depth <= 7)
      cat("  (MIDORI2 ranks: phylum to species; NCBI taxids are stripped from the names)\n")
    else if (is_pr2 && depth >= 8)
      cat("  (PR2 ranks: Division reported as Phylum)\n")
    else if (has_domain)
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


# MIDORI2 names carry the NCBI taxid ("Arthropoda_6656") and fill a missing rank
# with a placeholder built from the rank below ("order_Sciaenidae_30870"). The
# taxid is stripped; a placeholder means the rank is unknown, so it becomes NA.
tax_clean_midori <- function(tax) {
  if (is.null(tax)) return(tax)
  v <- tax
  v[] <- sub("_[0-9]+$", "", v)
  v[grepl("^(kingdom|phylum|class|order|family|genus|species)_", v)] <- NA
  v
}


# ── When the reference is too big for assignTaxonomy() ───────────────────────
# DADA2's classifier keeps one 8-mer profile (65,536 floats) per distinct lineage
# in the reference. Memory is therefore (lineages x 256 KB): SILVA at genus level
# ~1 GB, UNITE fungi with species ~9 GB, and MIDORI2 COI with species ~64 GB —
# which no analysis machine has, so a COX1 run against MIDORI2 died out of memory.
# These helpers measure that up front and, when it will not fit, classify with
# the best VSEARCH hit instead, keeping each rank only where the identity
# supports it.

tax_lineage_of <- function(h) {                      # vectorised
  h <- sub("^>", "", h)
  bar <- grepl("|", h, fixed = TRUE)                  # UNITE: name|acc|SH|ref|lineage
  h[bar] <- sub("^.*\\|", "", h[bar])
  h
}

# number of distinct lineages, read from the headers (UNITE, MIDORI2, SILVA, PR2)
tax_count_lineages <- function(db) {
  tryCatch({
    con <- if (grepl("\\.gz$", db)) gzfile(db, "rt") else file(db, "rt")
    on.exit(close(con), add = TRUE)
    seen <- new.env(hash = TRUE, size = 200000L)
    repeat {
      ch <- readLines(con, n = 200000L, warn = FALSE)
      if (length(ch) == 0) break
      for (l in unique(tax_lineage_of(ch[startsWith(ch, ">")])))
        assign(l, TRUE, envir = seen)
    }
    length(ls(seen, all.names = TRUE))
  }, error = function(e) NA_integer_)
}

tax_dada2_mem_gb <- function(db) {
  n <- tax_count_lineages(db)
  if (is.na(n)) NA_real_ else n * 65536 * 4 / 1e9
}

tax_avail_mem_gb <- function() {
  tryCatch({
    m <- readLines("/proc/meminfo", warn = FALSE)
    kb <- as.numeric(sub("\\D+(\\d+).*", "\\1", m[grepl("^MemAvailable:", m)]))
    kb / 1e6
  }, error = function(e) NA_real_)
}

tax_find_vsearch <- function() {
  v <- Sys.which("vsearch")
  if (nzchar(v)) return(unname(v))
  home <- Sys.getenv("HOME")
  for (b in c("miniconda3", "anaconda3", "miniforge3", "mambaforge"))
    for (e in c("emu", "ngamp-wgs", "ngamp-phylo", "wgs")) {
      p <- file.path(home, b, "envs", e, "bin", "vsearch")
      if (file.exists(p)) return(p)
    }
  ""
}

# Best-hit classification. `levels` names the ranks of the reference's lineage
# (as tax_pick_levels returns them); `cutoffs` is a named vector of the minimum
# identity (%) for keeping each rank. Returns a matrix like assignTaxonomy's:
# one row per sequence (rownames = sequences), one column per rank.
# A reference FASTA as pieces of at most n records (the file itself when smaller).
tax_split_fasta <- function(db, n = 250000L) {
  con <- if (grepl("\\.gz$", db)) gzfile(db, "rt") else file(db, "rt")
  on.exit(close(con), add = TRUE)
  files <- character(0); k <- 0L
  repeat {
    ch <- readLines(con, n = 1000000L, warn = FALSE)
    if (length(ch) == 0) break
    rec <- k + cumsum(startsWith(ch, ">"))            # record number of every line
    k <- rec[length(rec)]
    part <- (pmax(rec, 1L) - 1L) %/% n + 1L
    for (pi in unique(part)) {
      while (length(files) < pi) files <- c(files, tempfile(fileext = ".fasta"))
      cat(ch[part == pi], file = files[pi], sep = "\n", append = TRUE)
    }
  }
  if (length(files) <= 1L) { unlink(files); return(db) }
  files
}

tax_vsearch_tophit <- function(seqs, db, levels, cutoffs, threads = 1, verbose = TRUE,
                               near = 0.5, agree = 0.8) {
  # Best-hit naming with a consensus over near-best hits: every hit within
  # `near` % identity of the best is considered, uninformative names
  # ("unidentified", "Incertae_sedis", "Genus_sp", MIDORI2 placeholders) are
  # ignored, and a rank is named only when at least `agree` of the informative
  # names there agree — otherwise it and every rank below stay unassigned. A
  # single top hit alone let a reference record labelled "Ascomycota_sp" erase
  # the genus of a sequence that matched named Aspergillus records just as well.
  vs <- tax_find_vsearch()
  if (!nzchar(vs)) stop("vsearch not found — it is needed to classify against ", basename(db),
                        " (too large for DADA2's classifier on this machine)")
  fa  <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  writeLines(as.vector(rbind(paste0(">q", seq_along(seqs)), seqs)), fa)
  minid <- max(0.5, min(cutoffs[cutoffs > 0]) / 100)
  # vsearch indexes the whole reference in memory (~4 bytes per base: MIDORI2
  # COI is ~5 GB). A large reference is searched in pieces of 250,000 records
  # and the hits pooled, so the peak is a fraction of that.
  parts <- tax_split_fasta(db, 250000L)
  on.exit(unlink(setdiff(parts, db)), add = TRUE)
  if (verbose) cat(sprintf("  vsearch --usearch_global against %s (%d piece%s) ...\n",
                           basename(db), length(parts), if (length(parts) == 1) "" else "s"))
  hits <- NULL
  for (part in parts) {
    out <- tempfile(fileext = ".tsv")
    args <- c("--usearch_global", fa, "--db", part, "--id", sprintf("%.3f", minid),
              "--maxaccepts", "50", "--maxrejects", "200", "--maxhits", "50",
              "--strand", "both", "--userout", out, "--userfields", "query+target+id",
              "--threads", as.character(threads), "--quiet")
    rc <- system2(vs, args, stdout = FALSE, stderr = FALSE)
    if (!identical(rc, 0L) || !file.exists(out))
      stop("vsearch --usearch_global failed (exit ", rc, if (identical(rc, 137L)) ", out of memory" else "", ")")
    h <- tryCatch(utils::read.delim(out, header = FALSE, quote = "", comment.char = "",
                                    stringsAsFactors = FALSE),
                  error = function(e) NULL)
    unlink(out)
    if (!is.null(h) && nrow(h)) hits <- rbind(hits, h)
  }
  if (is.null(hits)) hits <- data.frame(V1 = character(0), V2 = character(0), V3 = numeric(0))
  clean <- function(lin) {
    x <- trimws(strsplit(tax_lineage_of(lin), ";")[[1]])
    x <- x[nzchar(x)]
    x <- sub("^[a-z]__", "", x)
    x <- sub("_[0-9]+$", "", x)                                      # MIDORI2 taxid
    x <- gsub("_", " ", x)
    x[grepl("^(kingdom|phylum|class|order|family|genus|species) ", x) |
      grepl("unidentified|incertae sedis", x, ignore.case = TRUE) |
      grepl("( sp\\.?$)|( sp [0-9])", x)] <- NA
    x[seq_len(min(length(x), length(levels)))]
  }
  tax <- matrix(NA_character_, nrow = length(seqs), ncol = length(levels),
                dimnames = list(seqs, levels))
  n_hit <- 0L
  for (q in split(hits, hits$V1)) {
    qi  <- as.integer(sub("^q", "", q$V1[1]))
    top <- max(q$V3)
    near_hits <- q[q$V3 >= top - near, , drop = FALSE]
    lins <- lapply(near_hits$V2, clean)
    n_hit <- n_hit + 1L
    for (j in seq_along(levels)) {
      cut <- if (!is.na(cutoffs[levels[j]])) cutoffs[[levels[j]]] else 0
      if (top < cut) break
      v <- vapply(lins, function(x) if (length(x) >= j) x[j] else NA_character_, character(1))
      v <- v[!is.na(v)]
      if (length(v) == 0) next                 # no informative name here; deeper ranks may have one
      tab <- sort(table(v), decreasing = TRUE)
      if (tab[1] / length(v) < agree) break    # informative names disagree: stop (LCA)
      tax[qi, j] <- names(tab)[1]
    }
  }
  if (verbose) cat(sprintf("  %d of %d sequences had a hit at >= %.0f%% identity\n",
                           n_hit, length(seqs), 100 * minid))
  tax
}

# Rank cut-offs (% identity) for best-hit naming.
#   COI: Ratnasingham & Hebert BIN / common metabarcoding practice — species 97,
#        genus 95, family 90, order 85, class 80, phylum 75.
#   ITS: Vu et al. 2019 (Stud. Mycol. 92) — species 98.5, genus 94.3, family 88.5,
#        order 81.2, class 80.9, phylum 75.
TAX_CUTOFFS_COI <- c(Kingdom = 0, Phylum = 75, Class = 80, Order = 85, Family = 90, Genus = 95, Species = 97)
TAX_CUTOFFS_ITS <- c(Kingdom = 0, Phylum = 75, Class = 80.9, Order = 81.2, Family = 88.5, Genus = 94.3, Species = 98.5)
