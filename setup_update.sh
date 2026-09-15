#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# NextGen-Amplicon — Initial Setup Script (run once on a new machine)
# ══════════════════════════════════════════════════════════════════════════════
# This script:
#   1. Checks for Git, Python 3, Node / npm
#   2. Clones (or verifies) the private GitHub repo
#   3. Saves your GitHub Personal Access Token for auto-updates
#   4. Installs Python dependencies (venv)
#   5. Installs Node dependencies (npm install)
#
# Usage:
#   bash setup_update.sh
#   bash setup_update.sh --token ghp_xxxx --skip-clone
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

REPO_URL="https://github.com/bosskub121212/NextGen-amplicon.git"
INSTALL_DIR="$HOME/r16s-app"
TOKEN_FILE="$HOME/.config/amplicon/github_token"
TOKEN=""
SKIP_CLONE=false

# ── Parse flags ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --token)      TOKEN="$2"; shift 2 ;;
    --skip-clone) SKIP_CLONE=true;  shift ;;
    --dir)        INSTALL_DIR="$2"; shift 2 ;;
    *) echo "Unknown flag: $1"; exit 1 ;;
  esac
done

echo ""
echo "═══════════════════════════════════════════════════════"
echo "  NextGen-Amplicon  —  Setup"
echo "═══════════════════════════════════════════════════════"

# ── 1. Dependency checks ──────────────────────────────────────────────────────
echo ""
echo "[1/5] Checking dependencies..."

check_cmd() {
  if command -v "$1" &>/dev/null; then
    echo "  ✅ $1 found ($(command -v "$1"))"
  else
    echo "  ❌ $1 NOT found — please install it first"
    exit 1
  fi
}

check_cmd git
check_cmd python3
check_cmd node || check_cmd nodejs
# npm may need NVM
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && source "$NVM_DIR/nvm.sh"
check_cmd npm

# ── 2. GitHub token ───────────────────────────────────────────────────────────
echo ""
echo "[2/5] GitHub token setup..."

if [[ -z "$TOKEN" && -f "$TOKEN_FILE" ]]; then
  echo "  ℹ️  Token already saved at $TOKEN_FILE — skipping prompt."
  echo "     (delete that file to re-enter a new token)"
  TOKEN=$(cat "$TOKEN_FILE")
elif [[ -z "$TOKEN" ]]; then
  echo "  Enter your GitHub Personal Access Token (PAT)."
  echo "  The token needs read access to the private repo."
  echo "  (Generate one at: https://github.com/settings/tokens)"
  echo ""
  read -rsp "  Token (input hidden): " TOKEN
  echo ""
fi

if [[ -z "$TOKEN" ]]; then
  echo "  ❌ No token provided. Cannot access private repo."
  exit 1
fi

# Save token
mkdir -p "$(dirname "$TOKEN_FILE")"
echo "$TOKEN" > "$TOKEN_FILE"
chmod 600 "$TOKEN_FILE"
echo "  ✅ Token saved to $TOKEN_FILE"

# ── 3. Clone / update repo ───────────────────────────────────────────────────
echo ""
echo "[3/5] Repository setup..."

if [[ "$SKIP_CLONE" == true ]]; then
  echo "  ⏭  --skip-clone set — skipping git clone."
elif [[ -d "$INSTALL_DIR/.git" ]]; then
  echo "  ℹ️  Repo already exists at $INSTALL_DIR."
  echo "  Running git pull to ensure it's up to date..."
  git -C "$INSTALL_DIR" remote set-url origin "https://${TOKEN}@github.com/bosskub121212/NextGen-amplicon.git"
  git -C "$INSTALL_DIR" pull origin main
else
  echo "  Cloning into $INSTALL_DIR ..."
  git clone "https://${TOKEN}@github.com/bosskub121212/NextGen-amplicon.git" "$INSTALL_DIR"
fi

# Set remote URL with token so future git pulls work
if [[ -d "$INSTALL_DIR/.git" ]]; then
  git -C "$INSTALL_DIR" remote set-url origin "https://${TOKEN}@github.com/bosskub121212/NextGen-amplicon.git"
fi

# ── 4. Python venv ────────────────────────────────────────────────────────────
echo ""
echo "[4/5] Python virtual environment..."

VENV="$INSTALL_DIR/venv"
if [[ -d "$VENV" ]]; then
  echo "  ℹ️  venv already exists."
else
  python3 -m venv "$VENV"
  echo "  ✅ Created venv at $VENV"
fi

source "$VENV/bin/activate"
pip install --upgrade pip -q
pip install -r "$INSTALL_DIR/backend/requirements.txt" -q
echo "  ✅ Python dependencies installed"

# ── 5. Node dependencies ──────────────────────────────────────────────────────
echo ""
echo "[5/5] Node dependencies..."

cd "$INSTALL_DIR/frontend"
npm install --legacy-peer-deps
echo "  ✅ Node dependencies installed"

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════"
echo "  ✅  Setup complete!"
echo ""
echo "  To launch the app:"
echo "    • Double-click launch_app.bat       (from Windows, for other machines)"
echo "    • Double-click launch_app_local.bat (from Windows, for this machine)"
echo ""
echo "  Or manually in WSL:"
echo "    cd ~/r16s-app && source venv/bin/activate"
echo "    uvicorn backend.main:app --reload --port 8000"
echo "    (in another terminal)"
echo "    cd ~/r16s-app/frontend && npm run dev"
echo "═══════════════════════════════════════════════════════"
echo ""
