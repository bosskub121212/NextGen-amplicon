#!/usr/bin/env bash
# =============================================================
#  setup_emu_db.sh — Install Emu + build SILVA database
#  for NextGen-Amplicon ONT 16S pipeline
#
#  Run from any machine that already has ~/r16s-app installed:
#    bash /mnt/c/Claude/setup_emu_db.sh
#  OR after cloning the repo:
#    bash ~/r16s-app/setup_emu_db.sh
# =============================================================
set -euo pipefail

GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[1;33m'
RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YELLOW}  ⚠${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }

APP_DIR="$HOME/r16s-app"
DB_DIR="$APP_DIR/backend/databases/emu_silva"
PREP_DIR="$DB_DIR/prep"
SILVA_DB_DIR="$PREP_DIR/silva_db"
DB_PATHS_JSON="$APP_DIR/backend/databases/db_paths.json"

# Detect SILVA DADA2 trainset (v138.1 or v138.2)
SILVA_FASTA=""
for candidate in \
    "$APP_DIR/backend/databases/SILVA/silva_nr99_v138.2_toGenus_trainset.fa.gz" \
    "$APP_DIR/backend/databases/SILVA/silva_nr99_v138.2_wSpecies_train_set.fa.gz" \
    "$APP_DIR/backend/databases/SILVA/silva_nr99_v138.1_wSpecies_train_set.fa.gz" \
    "$APP_DIR/backend/databases/SILVA/silva_nr99_v138_wSpecies_train_set.fa.gz"; do
    if [[ -f "$candidate" ]]; then
        SILVA_FASTA="$candidate"
        break
    fi
done

echo ""
echo -e "${BOLD}${CYAN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║   NextGen-Amplicon — Emu ONT 16S Database Setup         ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ── Step 1: conda check ───────────────────────────────────────
step "Checking conda"
if ! command -v conda &>/dev/null; then
    die "conda not found. Install Miniconda3 first:
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh"
fi
ok "conda found: $(conda --version)"

# ── Step 2: Install Emu conda env ────────────────────────────
step "Emu conda environment"
if conda run -n emu emu --version &>/dev/null 2>&1; then
    EMU_VER=$(conda run -n emu emu --version 2>&1 | head -1)
    ok "Emu already installed ($EMU_VER)"
else
    info "Creating conda env 'emu' (bioconda + conda-forge)…"
    _CONDA_CMD="conda"
    command -v mamba &>/dev/null && _CONDA_CMD="mamba"
    $_CONDA_CMD create -n emu -c bioconda -c conda-forge emu -y
    if conda run -n emu emu --version &>/dev/null 2>&1; then
        EMU_VER=$(conda run -n emu emu --version 2>&1 | head -1)
        ok "Emu installed ($EMU_VER)"
    else
        die "Emu install failed. Try: conda create -n emu -c bioconda -c conda-forge emu -y"
    fi
fi

# ── Step 3: Check if DB already built ────────────────────────
step "Checking existing Emu database"
if [[ -f "$SILVA_DB_DIR/species_taxid.fasta" && -f "$SILVA_DB_DIR/taxonomy.tsv" ]]; then
    N=$(grep -c "^>" "$SILVA_DB_DIR/species_taxid.fasta" 2>/dev/null || echo "?")
    ok "Database already exists at $SILVA_DB_DIR"
    ok "  species_taxid.fasta — $N sequences"
    info "Skipping build. Delete $SILVA_DB_DIR to rebuild."
    _SKIP_BUILD=true
else
    _SKIP_BUILD=false
fi

# ── Step 4: Locate SILVA FASTA ───────────────────────────────
if [[ "$_SKIP_BUILD" == "false" ]]; then
    step "Locating SILVA DADA2 trainset"

    if [[ -z "$SILVA_FASTA" ]]; then
        # Try to download SILVA v138.2 toGenus from Zenodo
        SILVA_DEST="$APP_DIR/backend/databases/SILVA/silva_nr99_v138.2_toGenus_trainset.fa.gz"
        mkdir -p "$(dirname "$SILVA_DEST")"
        echo ""
        info "SILVA DADA2 trainset not found locally."
        info "Attempting download from Zenodo (~130 MB)…"
        ZENODO_URL="https://zenodo.org/records/8392695/files/silva_nr99_v138.2_toGenus_train_set.fa.gz"
        if wget -q --show-progress -O "$SILVA_DEST" "$ZENODO_URL"; then
            ok "Downloaded to $SILVA_DEST"
            SILVA_FASTA="$SILVA_DEST"
        else
            rm -f "$SILVA_DEST"
            die "Download failed.
    Please download manually:
      wget -O \"$SILVA_DEST\" \\
        \"https://zenodo.org/records/8392695/files/silva_nr99_v138.2_toGenus_train_set.fa.gz\"
    Then re-run this script."
        fi
    fi

    ok "Using: $SILVA_FASTA"

    # ── Step 5: Build intermediate files ─────────────────────
    step "Building intermediate files (sequences, seq2taxid, taxonomy_list)"
    mkdir -p "$PREP_DIR"
    BUILD_SCRIPT="$APP_DIR/backend/python_scripts/build_emu_db.py"
    if [[ ! -f "$BUILD_SCRIPT" ]]; then
        die "build_emu_db.py not found at $BUILD_SCRIPT — is ~/r16s-app up to date?"
    fi
    python3 "$BUILD_SCRIPT" "$SILVA_FASTA" "$PREP_DIR"
    ok "Intermediate files written to $PREP_DIR"

    # ── Step 6: Run emu build-database ───────────────────────
    step "Running emu build-database (~2-5 min)"
    conda run --no-capture-output -n emu \
        emu build-database \
        --sequences     "$PREP_DIR/sequences.fasta" \
        --seq2tax       "$PREP_DIR/seq2taxid.tsv" \
        --taxonomy-list "$PREP_DIR/taxonomy_list.tsv" \
        --db-name       silva_db \
        --output-dir    "$PREP_DIR"

    if [[ -f "$SILVA_DB_DIR/species_taxid.fasta" && -f "$SILVA_DB_DIR/taxonomy.tsv" ]]; then
        N=$(grep -c "^>" "$SILVA_DB_DIR/species_taxid.fasta")
        ok "Database built: $SILVA_DB_DIR  ($N sequences)"
    else
        die "build-database completed but output files not found in $SILVA_DB_DIR"
    fi
fi

# ── Step 7: Update db_paths.json ─────────────────────────────
step "Updating db_paths.json"
python3 - <<PYEOF
import json, pathlib
p = pathlib.Path("$DB_PATHS_JSON")
if p.exists():
    d = json.loads(p.read_text())
else:
    d = {}
d["emu_silva"] = "$SILVA_DB_DIR"
p.write_text(json.dumps(d, indent=2))
print(f"  emu_silva → $SILVA_DB_DIR")
PYEOF
ok "db_paths.json updated"

# ── Done ─────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║   ✅  Emu database ready!                                ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"
echo "  Database : $SILVA_DB_DIR"
echo "  db_paths : $DB_PATHS_JSON"
echo ""
echo "  Quick test:"
echo "    conda run --no-capture-output -n emu \\"
echo "      emu abundance --help"
echo ""
echo "  Then restart the backend:"
echo "    pkill -f 'uvicorn main:app'; sleep 2; bash ~/r16s-app/start_backend.sh &"
echo ""
