#!/usr/bin/env bash
# =============================================================================
#  NextGen-Amplicon — repair the emu_silva database (ONT 16S)
#
#    bash ~/r16s-app/fix_emu_silva.sh
#
#  An Emu database is a folder with species_taxid.fasta and a taxonomy.tsv
#  whose header starts with "tax_id". On some machines databases/emu_silva holds
#  only build_emu_db.py's INTERMEDIATE files (sequences.fasta, seq2taxid.tsv,
#  taxonomy_list.tsv, taxonomy.tsv -> taxonomy_list.tsv) — `emu build-database`
#  ran into emu_silva/prep/silva_db instead, or never ran. Emu then stops with
#      KeyError: "None of ['tax_id'] are in the columns"
#
#  This script makes emu_silva a working Emu database:
#    1. uses emu_silva/prep/silva_db when it is a complete build, else
#    2. builds one with `emu build-database` from the intermediate files,
#  then moves the built files into emu_silva/ and everything else into
#  ~/r16s-db-backup/emu_silva_src/ (moved, not deleted). Safe to re-run.
# =============================================================================
set -uo pipefail
GREEN='\033[0;32m'; YEL='\033[1;33m'; CYAN='\033[0;36m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YEL}  !${NC}  $1"; }
fail() { echo -e "${RED}  ✗${NC}  $1"; exit 1; }

D="$HOME/r16s-app/backend/databases"
S="$D/emu_silva"
BK="$HOME/r16s-db-backup"          # outside the app: the database list never shows it
SRC="$BK/emu_silva_src"
[[ -d "$S" ]] || fail "$S not found"

is_built() {   # is_built <dir>
  [[ -s "$1/species_taxid.fasta" ]] && head -1 "$1/taxonomy.tsv" 2>/dev/null | grep -q "^tax_id"
}

if is_built "$S" && [[ ! -L "$S/taxonomy.tsv" ]]; then
  ok "emu_silva is already a built Emu database — nothing to do"
  exit 0
fi

BUILT=""
if is_built "$S/prep/silva_db"; then
  BUILT="$S/prep/silva_db"
  ok "found a complete build in emu_silva/prep/silva_db"
else
  info "no complete build found — running emu build-database (10–30 min)"
  EMU=""
  for b in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3"; do
    [[ -x "$b/envs/emu/bin/emu" ]] && EMU="$b/envs/emu/bin" && break
  done
  [[ -n "$EMU" ]] || fail "Emu conda env not found (expected ~/miniconda3/envs/emu)"
  # the intermediate files: emu_silva/ first, emu_silva/prep/ second
  for d in "$S" "$S/prep"; do
    if [[ -s "$d/sequences.fasta" && -s "$d/seq2taxid.tsv" && -s "$d/taxonomy_list.tsv" ]]; then
      IN="$d"; break
    fi
  done
  [[ -n "${IN:-}" ]] || fail "no sequences.fasta + seq2taxid.tsv + taxonomy_list.tsv to build from"
  mkdir -p "$BK"
  ( cd "$BK" && rm -rf emu_silva_built && \
    PATH="$EMU:$PATH" "$EMU/python3" "$EMU/emu" build-database emu_silva_built \
      --sequences "$IN/sequences.fasta" --seq2tax "$IN/seq2taxid.tsv" \
      --taxonomy-list "$IN/taxonomy_list.tsv" ) || fail "emu build-database failed"
  is_built "$BK/emu_silva_built" || fail "build finished but the output is incomplete"
  BUILT="$BK/emu_silva_built"
  ok "built: $BUILT"
fi

# move everything that is not the database out of emu_silva/
mkdir -p "$SRC"
shopt -s dotglob nullglob
for f in "$S"/*; do
  [[ "$(basename "$f")" == "prep" && "$BUILT" == "$S/prep/silva_db" ]] && continue
  if [[ -L "$f" ]]; then rm -f "$f"; continue; fi        # the taxonomy.tsv symlink
  mv "$f" "$SRC/" && info "moved $(basename "$f") → ~/r16s-db-backup/emu_silva_src/"
done
# the built database into emu_silva/
for f in "$BUILT"/*; do
  mv "$f" "$S/" && ok "emu_silva/$(basename "$f")"
done
# what is left of prep (its own intermediate files) goes to the source folder
if [[ -d "$S/prep" ]]; then
  rmdir "$S/prep/silva_db" 2>/dev/null || true
  mv "$S/prep" "$SRC/prep" && info "moved prep/ → ~/r16s-db-backup/emu_silva_src/prep/"
fi

if is_built "$S"; then
  ok "emu_silva is a working Emu database:"
  ls -la "$S"
  echo "      taxonomy.tsv: $(head -1 "$S/taxonomy.tsv" | cut -c1-90)"
  echo "      $(grep -c '>' "$S/species_taxid.fasta") reference sequences"
  info "Intermediate files kept in $SRC — delete that folder once a run with emu_silva has worked."
  info "No backend restart needed. In Beta pick: emu_silva → species_taxid.fasta"
else
  fail "emu_silva is still incomplete — send this output"
fi
