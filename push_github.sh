#!/usr/bin/env bash
# =============================================================================
#  push_github.sh — Commit and push the app repo to GitHub over SSH
#
#  Run from WSL:
#    bash ~/r16s-app/push_github.sh "v2.8.2: report builder"
#    bash ~/r16s-app/push_github.sh                    # message from version.json
#    bash ~/r16s-app/push_github.sh -m "msg" --yes     # no confirmation prompt
#    bash ~/r16s-app/push_github.sh --dry-run          # show what would go, push nothing
#
#  Deliberately separate from deploy_dev.sh: updating the running app and
#  publishing source are different decisions with different risks, and you
#  often want the first without the second.
#
#  WHY SSH AND NOT A TOKEN
#  -----------------------
#  The previous version read a Personal Access Token and wrote it into the
#  remote URL. That is how the old token leaked: when authentication fails or
#  git prompts, git prints the WHOLE remote URL — token included — to the
#  terminal, where it lands in scrollback, screenshots and pasted logs. Masking
#  the script's own output does not help, because the line came from git.
#
#  With SSH the secret never appears in any URL, any argument, or any message:
#  the private key stays in ~/.ssh and ssh does the authentication. There is
#  also no expiry date to manage.
#
#  If a token-bearing remote is still configured, this script rewrites it to
#  SSH and says so — leaving it in place would keep the old secret sitting in
#  .git/config in plain text.
# =============================================================================
set -euo pipefail

GREEN='\033[0;32m'; CYAN='\033[0;36m'; YEL='\033[1;33m'; RED='\033[0;31m'
BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓${NC}  $1"; }
info() { echo -e "${CYAN}  ℹ${NC}  $1"; }
warn() { echo -e "${YEL}  !${NC}  $1"; }
die()  { echo -e "${RED}  ✗${NC}  $1"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $1 ══${NC}"; }
# Kept as defence in depth: if any path ever reintroduces a credentialed URL,
# nothing this script prints should carry it.
mask() { sed -E 's#://[^/@]*@#://***@#g'; }

APP_DIR="$HOME/r16s-app"
REPO_PATH="bosskub121212/NextGen-amplicon"
SSH_URL="git@github.com:${REPO_PATH}.git"
OLD_TOKEN_FILE="$HOME/.config/amplicon/github_token"

MSG=""
ASSUME_YES=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--message) MSG="${2:-}"; shift 2 ;;
    -y|--yes)     ASSUME_YES=1; shift ;;
    --dry-run)    DRY_RUN=1; shift ;;
    -h|--help)    sed -n '2,10p' "$0"; exit 0 ;;
    -*)           die "Unknown option: $1" ;;
    *)            MSG="$1"; shift ;;
  esac
done

cd "$APP_DIR" || die "Repo not found: $APP_DIR"
git rev-parse --git-dir >/dev/null 2>&1 || die "$APP_DIR is not a git repository"

# Default message from version.json
#
# git expects a short subject line, a blank line, then the detail. The earlier
# version pasted the whole changelog into the subject, which made
# `git log --oneline` a wall of text and left GitHub's commit list truncating
# mid-sentence. Take the first sentence as the subject, keep the rest as the
# body.
build_message() {
  local V="$1" C="$2" SUBJ FIRST TRUNCATED=0
  if [[ -z "$C" ]]; then printf '%s' "v$V"; return; fi
  # Split on ". " followed by a capital, so "e.g." and version numbers do not
  # count as the end of a sentence.
  FIRST="$(printf '%s' "$C" | sed -E 's/^(.*?[^A-Z][.])[[:space:]]+[A-Z].*$/\1/')"
  SUBJ="v$V: $FIRST"
  if [[ ${#SUBJ} -gt 72 ]]; then          # git convention: subject <= 72 chars
    SUBJ="$(printf '%s' "${SUBJ:0:69}" | sed -E 's/[[:space:],;:.]+[^[:space:]]*$//')…"
    TRUNCATED=1
  fi
  # Repeating the changelog as a body adds nothing when the subject already
  # carries all of it.
  if [[ "$TRUNCATED" -eq 0 && "$FIRST" == "$C" ]]; then
    printf '%s' "$SUBJ"
  else
    printf '%s\n\n%s' "$SUBJ" "$(printf '%s' "$C" | fold -s -w 72)"
  fi
}

if [[ -z "$MSG" ]]; then
  V=$(grep -oP '"version"\s*:\s*"\K[^"]+' version.json 2>/dev/null || echo "")
  C=$(grep -oP '"changelog"\s*:\s*"\K[^"]+' version.json 2>/dev/null || echo "")
  if [[ -n "$V" ]]; then
    MSG="$(build_message "$V" "$C")"
    info "Message from version.json"
  else
    die "No commit message given and version.json has no version"
  fi
fi

# ── 0. Remote and key ─────────────────────────────────────────────────────
step "Checking SSH access"

CUR_URL="$(git remote get-url origin 2>/dev/null || echo "")"
if [[ -z "$CUR_URL" ]]; then
  git remote add origin "$SSH_URL"
  ok "Remote added: $SSH_URL"
elif [[ "$CUR_URL" == *"@github.com:"* ]]; then
  ok "Remote already SSH: $CUR_URL"
else
  # An https:// remote, with or without an embedded token. Either way it is
  # replaced — and if it carried a token, that token has been sitting in
  # .git/config in plain text and must be treated as compromised.
  if [[ "$CUR_URL" == *"@github.com"* ]]; then
    warn "The remote had a credential embedded in its URL."
    warn "That token was stored in .git/config in plain text — REVOKE IT at"
    warn "  https://github.com/settings/tokens"
  fi
  git remote set-url origin "$SSH_URL"
  ok "Remote switched to SSH: $SSH_URL"
fi

if [[ -e "$OLD_TOKEN_FILE" ]]; then
  warn "An old token file is still on disk: $OLD_TOKEN_FILE"
  warn "It is no longer used. Revoke the token, then:  rm $OLD_TOKEN_FILE"
fi

# ssh -T github.com exits 1 even on success ("You've successfully
# authenticated, but GitHub does not provide shell access"), so test the
# message, not the exit code.
# -n matters: without it ssh inherits the terminal as its stdin and swallows
# it, so the "Commit and push this?" prompt further down reads EOF instead of
# the keystroke and the script reports "Cancelled" no matter what you type.
SSH_OUT="$(ssh -n -o BatchMode=yes -o StrictHostKeyChecking=accept-new \
              -T git@github.com 2>&1 </dev/null || true)"
if echo "$SSH_OUT" | grep -q "successfully authenticated"; then
  ok "$(echo "$SSH_OUT" | head -1)"
else
  echo ""
  warn "GitHub did not accept an SSH key. It said:"
  echo "$SSH_OUT" | sed 's/^/      /'
  echo ""
  echo "  Set one up once, then re-run this script:"
  echo ""
  echo "    ssh-keygen -t ed25519 -C \"nextgen-amplicon\$(date +%Y%m%d)\" -f ~/.ssh/id_ed25519"
  echo "    cat ~/.ssh/id_ed25519.pub"
  echo ""
  echo "  Paste that public key at  https://github.com/settings/keys"
  echo "  (New SSH key → Authentication key). The .pub file is safe to paste;"
  echo "  the file WITHOUT .pub is the private key and never leaves this machine."
  echo ""
  echo "  If your network blocks port 22, use GitHub's SSH over 443:"
  echo "    printf 'Host github.com\\n  Hostname ssh.github.com\\n  Port 443\\n  User git\\n' >> ~/.ssh/config"
  echo ""
  die "No SSH access yet. Nothing was committed."
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
  echo -e "  ${BOLD}Commit message:${NC} $(printf '%s' "$MSG" | head -1)"
  if [[ "$(printf '%s' "$MSG" | wc -l)" -gt 0 ]]; then
    printf '%s' "$MSG" | tail -n +3 | sed 's/^/      /'
  fi
  # Read the answer from the terminal itself, not from whatever stdin happens
  # to be — and accept "yes" as well as "y".
  ANSWER=""
  if [[ -r /dev/tty ]]; then
    read -r -p "  Commit and push this? [y/N] " ANSWER </dev/tty || ANSWER=""
  else
    read -r -p "  Commit and push this? [y/N] " ANSWER || ANSWER=""
  fi
  if [[ ! "$ANSWER" =~ ^([Yy]|[Yy][Ee][Ss])$ ]]; then
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

# ── 3. Push ───────────────────────────────────────────────────────────────
step "Pushing"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
set +e
PUSH_OUT="$(git push -u origin "$BRANCH" 2>&1)"
PUSH_RC=$?
set -e
echo "$PUSH_OUT" | mask | sed 's/^/    /'

if [[ "$PUSH_RC" -ne 0 ]]; then
  echo ""
  if echo "$PUSH_OUT" | grep -qi "permission denied\|publickey"; then
    warn "The key reached GitHub but was refused for this repository."
    warn "Check the key is on the account that owns $REPO_PATH, or add it as a"
    warn "deploy key with write access on the repo itself."
  elif echo "$PUSH_OUT" | grep -qi "rejected\|non-fast-forward"; then
    warn "The remote has commits you do not have locally. Pull first:"
    warn "  git -C $APP_DIR pull --rebase"
  fi
  die "Push failed. Your commit is safe locally; fix the cause and re-run."
fi

echo ""
echo -e "${BOLD}${GREEN}  ✅  Pushed to $REPO_PATH ($BRANCH)${NC}"
echo ""
echo "  Other machines can now pull it:"
echo "    bash ~/r16s-app/update.sh"
echo ""
