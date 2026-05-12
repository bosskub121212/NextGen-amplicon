#!/usr/bin/env bash
# Download SILVA 138.1 taxonomy databases for DADA2
# Run from WSL: bash /mnt/c/Claude/download_silva.sh

SILVA_DIR="$HOME/r16s-app/backend/databases/SILVA"
mkdir -p "$SILVA_DIR"
cd "$SILVA_DIR"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Downloading SILVA 138.1 databases for DADA2           ║"
echo "║   Destination: $SILVA_DIR"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# ── File 1: Genus-level taxonomy training set (137 MB) ──────────
FILE1="silva_nr99_v138.1_train_set.fa.gz"
if [[ -f "$FILE1" ]]; then
  echo "  ✓ $FILE1 already exists — skipping"
else
  echo "  ↓ Downloading $FILE1 (137 MB)..."
  wget -q --show-progress \
    "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1" \
    -O "$FILE1"
  echo "  ✓ $FILE1 done"
fi

echo ""

# ── File 2: Species-level assignment (79 MB) ────────────────────
FILE2="silva_species_assignment_v138.1.fa.gz"
if [[ -f "$FILE2" ]]; then
  echo "  ✓ $FILE2 already exists — skipping"
else
  echo "  ↓ Downloading $FILE2 (79 MB)..."
  wget -q --show-progress \
    "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1" \
    -O "$FILE2"
  echo "  ✓ $FILE2 done"
fi

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Done! Files in:                                       ║"
echo "║   ~/r16s-app/backend/databases/SILVA/                   ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "  Next: In the app → Settings → Database Paths"
echo "        Set SILVA path to: ~/r16s-app/backend/databases/SILVA/silva_nr99_v138.1_train_set.fa.gz"
echo ""
ls -lh "$SILVA_DIR"/*.fa.gz 2>/dev/null
