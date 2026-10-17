#!/usr/bin/env bash
# Kept for old instructions that still call it. The reference downloads now live
# in one place — download_databases.sh — which fetches SILVA 144 (and 138.2 with
# --silva138) with checksums. This script used to say "138.1" in some places and
# "138.2" in others for the same 138.1 files.
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$APP_DIR/download_databases.sh" ]] || APP_DIR="$HOME/r16s-app"
exec bash "$APP_DIR/download_databases.sh" --only silva "$@"
