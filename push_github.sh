#!/usr/bin/env bash
# =============================================================================
#  push_github.sh — Commit and push the app repo to GitHub
#
#  Run from WSL:
#    bash ~/r16s-app/push_github.sh "v2.7.1: fix tree crash"
#    bash ~/r16s-app/push_github.sh                    # message from version.json
#    bash ~/r16s-app/push_github.sh -m "msg" --yes     # no confirmation prompt
#    bash ~/r16s-app/push_github.sh --dry-run          # show what would go, push nothing
#
#  Deliberately separate from deploy_dev.sh: updating the running app and
#  publishing source are different decisions with different risks, and you
#  often want the first without the second.
#
#  TOKEN HANDLING
#  --------------
#  The token is read from ~/.config/amplicon/github_token and never printed.
#  The remote is written as https://oauth2:<token>@github.com/... — the
#  oauth2: username matters: with the bare https://<token>@github.com/ form,
#  git treats the token as a username, prompts for a password, and ECHOES THE
#  WHOLE URL (token included) to the terminal, where it lands in scrollback,
#  screenshots and pasted logs. Every line this script prints masks it.
# =============================================================================
set -euo pipefail

GREEN='\033[0;32m'; CYAN='\033[0;36m'; YEL='\033[1;33m'; RED='\033[0;31m'
BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YEL}  !${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }
mask() { sed -E 's#://[^/@]*@#://***@#g'; }   # hide credentials in any URL

APP_DIR="$HOME/r16s-app"
TOKEN_FILE="$HOME/.config/amplicon/github_token"
REPO_PATH="bosskub121212/NextGen-amplicon"

MSG=""
ASSUME_YES=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--message) MSG="${2:-}"; shift 2 ;;
    -y|--yes)     ASSUME_YES=1; shift ;;
    --dry-run)    DRY_RUN=1; shift ;;
    -h|--help)    sed -n '2,18p' "$0"; exit 0 ;;
    -*)           die "Unknown option: $1" ;;
    *)            MSG="$1"; shift ;;
  esac
done

cd "$APP_DIR" || die "Repo not found: $APP_DIR"
git rev-parse --git-dir >/dev/null 2>&1 || die "$APP_DIR is not a git repository"

# Default message from version.json
if [[ -z "$MSG" ]]; then
  V=$(grep -oP '"version"\s*:\s*"\K[^"]+' version.json 2>/dev/null || echo "")
  C=$(grep -oP '"changelog"\s*:\s*"\K[^"]+' version.json 2>/dev/null || echo "")
  if [[ -n "$V" ]]; then
    MSG="v$V${C:+: $C}"
    info "Message from version.json"
  else
    die "No commit message given and version.json has no version"
  fi
fi

# ── 1. Show exactly what would be committed ───────────────────────────────
# A lockfile once slipped into a release commit because the old script ran
# `git add -A` and nobody looked. Always show the list before writing history.
step "Changes to commit"
git add -A
if git diff --cached --quiet; then
  info "Nothing to commit — working tree matches HEAD"
  git status -sb | head -1 | mask
  # With no upstream configured, `git rev-list @{u}..HEAD` errors out. Treating
  # that as "0 unpushed" would report the branch as already in sync and skip the
  # push entirely — exactly wrong on a fresh clone, where nothing is published yet.
  if git rev-parse --abbrev-ref '@{u}' >/dev/null 2>&1; then
    UNPUSHED=$(git rev-list --count '@{u}..HEAD' 2>/dev/null || echo 0)
    if [[ "$UNPUSHED" -gt 0 ]]; then
      info "$UNPUSHED commit(s) committed but not pushed — continuing to push"
    else
      ok "Already in sync with the remote"
      exit 0
    fi
  else
    info "No upstream branch set — will publish and set it on push"
  fi
else
  git diff --cached --stat | sed 's/^/    /'
  echo ""
  UNTRACKED_NEW=$(git diff --cached --name-only --diff-filter=A | wc -l)
  [[ "$UNTRACKED_NEW" -gt 0 ]] && info "$UNTRACKED_NEW new file(s) will be added"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo ""
  info "--dry-run: unstaging and exiting without committing"
  git reset >/dev/null
  exit 0
fi

if [[ "$ASSUME_YES" -eq 0 ]]; then
  echo ""
  echo -e "  ${BOLD}Commit message:${NC} $MSG"
  read -r -p "  Commit and push this? [y/N] " REPLY
  if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
    git reset >/dev/null
    info "Cancelled — nothing committed, changes left in the working tree"
    exit 0
  fi
fi

# ── 2. Commit ─────────────────────────────────────────────────────────────
step "Committing"
if git diff --cached --quiet; then
  info "No staged changes — pushing existing commits only"
else
  git commit -q -m "$MSG"
  ok "$(git log -1 --oneline)"
fi

# ── 3. Configure the remote from the stored token ─────────────────────────
step "Preparing remote"
if [[ ! -r "$TOKEN_FILE" ]]; then
  warn "Token file not found: $TOKEN_FILE"
  echo ""
  echo "  Create a Personal Access Token (scope: repo) at:"
  echo "    https://github.com/settings/tokens"
  echo "  then save it — the file, never this terminal:"
  echo "    mkdir -p ~/.config/amplicon"
  echo "    printf '%s' 'PASTE_TOKEN_HERE' > $TOKEN_FILE"
  echo "    chmod 600 $TOKEN_FILE"
  echo ""
  die "Cannot push without a token. The commit above is saved locally."
fi

chmod 600 "$TOKEN_FILE" 2>/dev/null || true
TOKEN="$(tr -d '[:space:]' < "$TOKEN_FILE")"
[[ -n "$TOKEN" ]] || die "Token file is empty: $TOKEN_FILE"

git remote set-url origin "https://oauth2:${TOKEN}@github.com/${REPO_PATH}.git"
unset TOKEN
ok "Remote set: $(git remote get-url origin | mask)"

# ── 4. Push ───────────────────────────────────────────────────────────────
step "Pushing"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
set +e
PUSH_OUT="$(git push -u origin "$BRANCH" 2>&1)"
PUSH_RC=$?
set -e
echo "$PUSH_OUT" | mask | sed 's/^/    /'

if [[ "$PUSH_RC" -ne 0 ]]; then
  echo ""
  if echo "$PUSH_OUT" | grep -qi "authentication failed\|could not read\|invalid username"; then
    warn "Authentication failed — the token is probably expired or revoked."
    warn "Generate a new one and overwrite: $TOKEN_FILE"
  fi
  die "Push failed. Your commit is safe locally; fix the cause and re-run."
fi

echo ""
echo -e "${BOLD}${GREEN}  ✅  Pushed to $REPO_PATH ($BRANCH)${NC}"
echo ""
echo "  Other machines can now pull it:"
echo "    bash ~/r16s-app/update.sh"
echo ""
