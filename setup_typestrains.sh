#!/usr/bin/env bash
# =============================================================================
#  setup_typestrains.sh — NCBI type-strain genomes for OFFLINE ONT-WGS runs
#
#  The type-strain comparison (ANIb / ANIm / TETRA, genome tree, reference
#  mapping) normally downloads the genomes of the isolate's genus from NCBI
#  during the run. This downloads them ahead of time, for every genus or for the
#  genera you name, into the cache the pipeline already reads (~/.ngamp_wgs_refs).
#
#    bash ~/r16s-app/setup_typestrains.sh --dry-run          # how many, how big
#    bash ~/r16s-app/setup_typestrains.sh                    # everything, in the background
#    bash ~/r16s-app/setup_typestrains.sh --genera Bacillus,Lactiplantibacillus
#    bash ~/r16s-app/setup_typestrains.sh --status
#    tail -f ~/.ngamp_wgs_refs/build.log                      # progress of a background build
#
#  Stop it any time (or lose the connection): running it again downloads only
#  what is missing. Re-run every few months to pick up newly described species.
#  Analyses can keep running while it downloads; no backend restart needed.
# =============================================================================
set -uo pipefail
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="$APP/venv/bin/python"; [[ -x "$PY" ]] || PY="$(command -v python3)"
TOOL="$APP/backend/wgs_typestrain_db.py"
CACHE="$HOME/.ngamp_wgs_refs"          # where ont_wgs_pipeline.py looks (--ref_cache)
LOG="$CACHE/build.log"
mkdir -p "$CACHE"

case "${1:-}" in
  --status)  exec "$PY" "$TOOL" status --cache "$CACHE" ;;
  --dry-run) exec "$PY" "$TOOL" build --dry-run --cache "$CACHE" ;;
  --genera)
    [[ -n "${2:-}" ]] || { echo "usage: $0 --genera Genus1,Genus2"; exit 2; }
    exec "$PY" "$TOOL" build --genera "$2" --cache "$CACHE" ;;
  --foreground)
    exec "$PY" "$TOOL" build --cache "$CACHE" ;;
  "")
    if pgrep -f "wgs_typestrain_db.py build" >/dev/null; then
      echo "A build is already running — progress: tail -f $LOG"; exit 0
    fi
    avail=$(df -BG --output=avail "$CACHE" | tail -1 | tr -dc 0-9)
    if [[ "${avail:-0}" -lt 60 ]]; then
      echo "Only ${avail} GB free where the cache lives ($CACHE); about 40 GB is needed."
      echo "Free some space, or download only the genera you need:  bash $0 --genera A,B"
      exit 1
    fi
    nohup "$PY" "$TOOL" build --cache "$CACHE" >>"$LOG" 2>&1 &
    echo "Downloading every NCBI type-strain genome in the background (several hours)."
    echo "  progress:  tail -f $LOG"
    echo "  status:    bash $0 --status"
    echo "Safe to close this window. Run the same command again to resume after a stop."
    ;;
  -h|--help) sed -n 2,20p "$0" ;;
  *) echo "unknown option $1 — see: bash $0 --help"; exit 2 ;;
esac
