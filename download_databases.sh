#!/usr/bin/env bash
# ============================================================
#  NextGen-Amplicon — Reference Database Downloader
#  Downloads:
#    1. SILVA 138.2   (16S / 12S rRNA)
#    2. UNITE v10     (ITS fungi)
#    3. PR2 v5        (18S eukaryotes / nematodes)
#    4. MIDORI2 COX1  (animal COX1 metabarcoding)
#    5. Tax4Fun2      (16S functional prediction)
#
#  Usage:  bash ~/r16s-app/download_databases.sh
# ============================================================
set -uo pipefail   # -u = undefined vars, -o pipefail; NO -e so downloads don't abort whole script

DB_DIR="$HOME/r16s-app/backend/databases"
mkdir -p "$DB_DIR"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "  ${GREEN}✔ $*${NC}"; }
warn() { echo -e "  ${YELLOW}⚠ $*${NC}"; }
fail() { echo -e "  ${RED}✘ $*${NC}"; }

echo ""
echo "============================================================"
echo "  NextGen-Amplicon  —  Reference Database Downloader"
echo "  Target: $DB_DIR"
echo "============================================================"
echo ""

# ── Robust download helper ─────────────────────────────────────
# try_download <label> <dest_file> <url> [<url2> ...]
# Returns 0 if file downloaded OK, 1 if all URLs failed
try_download() {
  local label="$1"; local dest="$2"; shift 2
  if [[ -f "$dest" && -s "$dest" ]]; then
    ok "Already exists: $(basename "$dest")"; return 0
  fi
  mkdir -p "$(dirname "$dest")"
  for url in "$@"; do
    echo "  Trying: $url"
    if wget -q --show-progress --timeout=60 --tries=2 \
            --user-agent="Mozilla/5.0" \
            -O "$dest" "$url" 2>/dev/null && [[ -s "$dest" ]]; then
      ok "$label"
      return 0
    fi
    rm -f "$dest"
  done
  fail "$label — all URLs failed"
  return 1
}

# ── 1. SILVA 138.2 — 16S / 12S rRNA ──────────────────────────
echo "[1/5] SILVA 138.2  (16S / 12S — ~130 MB)"
SILVA_DIR="$DB_DIR/SILVA"
mkdir -p "$SILVA_DIR"

try_download "SILVA 138.1 genus-level (16S train set)" \
  "$SILVA_DIR/silva_nr99_v138.1_train_set.fa.gz" \
  "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" \
  "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza" || \
warn "Download silva_nr99_v138.1_train_set.fa.gz from https://zenodo.org/record/4587955"

try_download "SILVA 138.1 species assignment" \
  "$SILVA_DIR/silva_species_assignment_v138.1.fa.gz" \
  "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz" || \
warn "Download silva_species_assignment_v138.1.fa.gz from https://zenodo.org/record/4587955"

# ── 2. UNITE v10 — ITS Fungi ──────────────────────────────────
echo ""
echo "[2/5] UNITE v10  (ITS fungi — ~300 MB)"
UNITE_DIR="$DB_DIR/UNITE"
mkdir -p "$UNITE_DIR"

UNITE_FASTA="$UNITE_DIR/sh_general_release_dynamic_all.fasta"
if [[ -f "$UNITE_FASTA" && -s "$UNITE_FASTA" ]]; then
  ok "Already exists: UNITE v10"
else
  UNITE_TGZ="$UNITE_DIR/unite_v10.tgz"
  DOWNLOADED=0
  for url in \
    "https://files.plutof.ut.ee/doi/10.15156/BIO/2938079" \
    "https://files.plutof.ut.ee/doi/10.15156/BIO/2959332" \
    "https://unite.ut.ee/sh_files/UNITE_public_all_04.04.2024.tgz"
  do
    echo "  Trying: $url"
    if wget -q --show-progress --timeout=120 --tries=2 \
            --user-agent="Mozilla/5.0" \
            -O "$UNITE_TGZ" "$url" 2>/dev/null && [[ -s "$UNITE_TGZ" ]]; then
      DOWNLOADED=1; break
    fi
    rm -f "$UNITE_TGZ"
  done

  if [[ $DOWNLOADED -eq 1 ]]; then
    echo "  Extracting UNITE..."
    tar -xzf "$UNITE_TGZ" -C "$UNITE_DIR" --strip-components=1 2>/dev/null || \
    tar -xzf "$UNITE_TGZ" -C "$UNITE_DIR" 2>/dev/null || true
    rm -f "$UNITE_TGZ"
    # Rename the extracted fasta to a known name if needed
    FOUND=$(find "$UNITE_DIR" -name "sh_general_release*.fasta" | head -1)
    if [[ -n "$FOUND" && "$FOUND" != "$UNITE_FASTA" ]]; then
      cp "$FOUND" "$UNITE_FASTA"
    fi
    if [[ -f "$UNITE_FASTA" ]]; then
      ok "UNITE v10"
    else
      warn "Extracted but FASTA not found at expected path. Check: $UNITE_DIR"
    fi
  else
    fail "UNITE v10 — all URLs failed"
    warn "Manual download: https://unite.ut.ee/repository.php"
    warn "Extract and place .fasta at: $UNITE_FASTA"
  fi
fi

# ── 3. PR2 v5 — 18S eukaryotes / nematodes ───────────────────
echo ""
echo "[3/5] PR2 v5  (18S eukaryotes / nematodes — ~50 MB)"
NEMA_DIR="$DB_DIR/NemaBase"
mkdir -p "$NEMA_DIR"

try_download "PR2 v5 SSU DADA2" \
  "$NEMA_DIR/pr2_v5_dada2.fasta.gz" \
  "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz" \
  "https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_dada2.fasta.gz" || \
warn "Download from https://github.com/pr2database/pr2database/releases"

# 18S NemaBase (Zenodo) — optional, PR2 is the reliable fallback
try_download "18S-NemaBase (nematode-specific)" \
  "$NEMA_DIR/18S-NemaBase_train.fasta.gz" \
  "https://zenodo.org/record/7660142/files/18S-NemaBase_train.fasta.gz" \
  "https://zenodo.org/record/7660142/files/NemaBase_18S_DADA2.fasta.gz" || \
warn "NemaBase not available — PR2 v5 will be used as fallback for 18S"

# ── 4. MIDORI2 COX1 — Animal metabarcoding ────────────────────
echo ""
echo "[4/5] MIDORI2 COX1  (animal COX1 — ~2 GB, may take a while)"
COX1_DIR="$DB_DIR/MIDORI2"
mkdir -p "$COX1_DIR"

MIDORI_FILE="MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz"
MIDORI_OUT="$COX1_DIR/$MIDORI_FILE"
MIDORI_BASE="https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank264/DADA2"

try_download "MIDORI2 COX1 genus" \
  "$MIDORI_OUT" \
  "${MIDORI_BASE}/${MIDORI_FILE}" \
  "https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank264/DADA2_sp/${MIDORI_FILE}" || {
  warn "MIDORI2 auto-download failed (their server blocks bots)"
  warn "Manual download:"
  warn "  1. Go to: https://www.reference-midori.info/download.php"
  warn "  2. Download: $MIDORI_FILE"
  warn "  3. Place at: $MIDORI_OUT"
}

# Species-level file (optional, for addSpecies step)
MIDORI_SP="MIDORI2_UNIQ_NUC_GB264_CO1_DADA2_sp.fasta.gz"
MIDORI_SP_OUT="$COX1_DIR/$MIDORI_SP"
if [[ ! -f "$MIDORI_SP_OUT" ]]; then
  wget -q --show-progress --timeout=120 --tries=1 \
       --user-agent="Mozilla/5.0" \
       -O "$MIDORI_SP_OUT" \
       "${MIDORI_BASE}/${MIDORI_SP}" 2>/dev/null && ok "MIDORI2 COX1 species" || \
  { rm -f "$MIDORI_SP_OUT"; warn "MIDORI2 species-level file skipped (optional)"; }
fi

# ── 5. Tax4Fun2 — 16S functional prediction ───────────────────
echo ""
echo "[5/5] Tax4Fun2 Ref99NR  (16S functional prediction — ~1 GB)"
TF2_DIR="$DB_DIR/Tax4Fun2"
mkdir -p "$TF2_DIR"

if [[ -d "$TF2_DIR/Tax4Fun2_ReferenceData_v2" ]]; then
  ok "Already exists: Tax4Fun2 reference data"
else
  TF2_ZIP="$TF2_DIR/Tax4Fun2_ReferenceData_v2.zip"
  try_download "Tax4Fun2 reference data" \
    "$TF2_ZIP" \
    "https://zenodo.org/record/5794357/files/Tax4Fun2_ReferenceData_v2.zip" \
    "https://cloudstor.aarnet.edu.au/plus/s/DnBdEDxRzE6FPrR/download" && {
    echo "  Extracting Tax4Fun2..."
    unzip -q "$TF2_ZIP" -d "$TF2_DIR" 2>/dev/null && rm -f "$TF2_ZIP" && ok "Tax4Fun2 extracted"
  } || {
    warn "Tax4Fun2 auto-download failed"
    warn "Alternative: run in R:"
    warn "  Tax4Fun2::downloadTax4Fun2Ref(path='$TF2_DIR')"
  }
fi

# ── Write db_paths.json ───────────────────────────────────────
echo ""
echo "Writing database paths config..."

# Find SILVA fasta (genus-level)
SILVA_FASTA=$(find "$SILVA_DIR" -name "silva_nr99*train*.fa*" | head -1)
SILVA_SP=$(find "$SILVA_DIR" -name "silva_species*.fa*" | head -1)

# Find UNITE fasta
UNITE_FASTA_FOUND=$(find "$UNITE_DIR" -name "*.fasta" | head -1)

# Find PR2
PR2_FASTA=$(find "$NEMA_DIR" -name "pr2_v5*.fasta.gz" | head -1)
NEMA_FASTA=$(find "$NEMA_DIR" -name "18S*.fasta.gz" | head -1)
[[ -z "$NEMA_FASTA" ]] && NEMA_FASTA="$PR2_FASTA"  # fallback to PR2 if NemaBase missing

cat > "$DB_DIR/db_paths.json" << EOF
{
  "SILVA_16S":       "${SILVA_FASTA:-$SILVA_DIR/silva_nr99_v138.1_train_set.fa.gz}",
  "SILVA_16S_sp":    "${SILVA_SP:-$SILVA_DIR/silva_species_assignment_v138.1.fa.gz}",
  "UNITE_ITS1":      "${UNITE_FASTA_FOUND:-$UNITE_DIR/sh_general_release_dynamic_all.fasta}",
  "UNITE_ITS2":      "${UNITE_FASTA_FOUND:-$UNITE_DIR/sh_general_release_dynamic_all.fasta}",
  "NemaBase_18S":    "${NEMA_FASTA:-$NEMA_DIR/18S-NemaBase_train.fasta.gz}",
  "PR2_18S":         "${PR2_FASTA:-$NEMA_DIR/pr2_v5_dada2.fasta.gz}",
  "MIDORI2_COX1":    "$COX1_DIR/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz",
  "MIDORI2_COX1_sp": "$COX1_DIR/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2_sp.fasta.gz",
  "Tax4Fun2_ref":    "$TF2_DIR/Tax4Fun2_ReferenceData_v2"
}
EOF
ok "db_paths.json saved to $DB_DIR/db_paths.json"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  Download summary"
echo ""
check_db() {
  local label="$1"; local path="$2"
  if [[ -f "$path" && -s "$path" ]] || [[ -d "$path" ]]; then
    echo -e "  ${GREEN}✔${NC} $label"
  else
    echo -e "  ${RED}✘${NC} $label  (${path})"
  fi
}

check_db "SILVA 138.2  (16S/12S)" \
  "${SILVA_FASTA:-$SILVA_DIR/silva_nr99_v138.1_train_set.fa.gz}"
check_db "UNITE v10    (ITS)" \
  "${UNITE_FASTA_FOUND:-$UNITE_DIR/sh_general_release_dynamic_all.fasta}"
check_db "PR2 v5       (18S)" \
  "${PR2_FASTA:-$NEMA_DIR/pr2_v5_dada2.fasta.gz}"
check_db "MIDORI2 COX1 (COX1)" \
  "$COX1_DIR/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz"
check_db "Tax4Fun2     (16S functional)" \
  "$TF2_DIR/Tax4Fun2_ReferenceData_v2"

echo ""
echo "  Config: $DB_DIR/db_paths.json"
echo "============================================================"
echo ""
