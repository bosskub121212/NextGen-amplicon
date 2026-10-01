#!/bin/bash
# =============================================================================
#  NextGen-Amplicon — Update Script
#  Run after each new release to pull code + rebuild frontend:
#
#    bash ~/r16s-app/update.sh
#
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   NextGen-Amplicon — Updating                           ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Pull latest code
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 1: Pulling latest code ─────────────────────────────"

# jobs_history.json is local-only data (job queue/history per machine).
# If git is still tracking it (was committed before .gitignore was added),
# we preserve the content, reset the tracked copy, pull, then restore.
JOBS_FILE="backend/jobs_history.json"
JOBS_BACKUP=""
if git ls-files --error-unmatch "$JOBS_FILE" &>/dev/null; then
  echo "  ℹ  Preserving local job history across update..."
  JOBS_BACKUP=$(cat "$JOBS_FILE" 2>/dev/null || echo "")
  git checkout -- "$JOBS_FILE" 2>/dev/null || true
  # Permanently stop tracking this file so this never happens again
  git rm --cached "$JOBS_FILE" 2>/dev/null || true
fi

# Files that are MACHINE OUTPUT but are tracked in git, because they were
# committed once before anyone noticed what they were. Every time this machine
# regenerates one, the next "git pull" aborts on it — an update that breaks the
# update. Discarding is always right for these (nobody edits build output), and
# untracking them is righter still, so do both: reset the copy, then drop it
# from the index permanently. The removal travels with the next push, and until
# then this sweep keeps every other machine unblocked.
#
#   frontend/package-lock.json — rewritten by "npm install". Step 2 now uses
#     "npm ci", which never writes it, so this should not recur; the sweep
#     stays for machines coming from an older version.
#   backend/Rplots.pdf — R's default graphics device. R opens it and writes to
#     it whenever a base-graphics call is made with no device already open, so
#     it is rewritten by EVERY pipeline run. On a machine whose job is running
#     pipelines, that means update.sh could never pull again. Nothing reads
#     this file; it is the place stray plots go to be lost.
REGENERATED=("frontend/package-lock.json" "backend/Rplots.pdf")

for f in "${REGENERATED[@]}"; do
  # REPAIR. An earlier version of this script ran `git rm --cached` here, which
  # leaves a STAGED DELETION in the index for good. On the next run the file is
  # no longer in the index, so the sweep below skips it — while the staged
  # deletion still shows up as a change and stops the update. That version
  # worked once per machine and then wedged it permanently, which is worse than
  # the problem it was written for. Unstage it and the machine recovers itself.
  if git diff --cached --name-only 2>/dev/null | grep -qxF "$f"; then
    git reset -q HEAD -- "$f" 2>/dev/null || true
    echo "  ℹ  Unstaged a leftover index change for $f"
  fi

  # Discard the working-tree copy if it is tracked and has been rewritten.
  # DISCARD, not untrack: whether a file belongs in the repo is one decision,
  # made once, in a commit from the development machine — not something each
  # runtime machine does to its own index behind the other machines' backs.
  if git ls-files --error-unmatch "$f" &>/dev/null && ! git diff --quiet -- "$f" 2>/dev/null; then
    echo "  ℹ  Discarding local changes to $f (machine output, not source)"
    git checkout -- "$f" 2>/dev/null || true
  fi
done

# Anything else modified is a real local edit. git pull would abort on it with a
# wall of text that does not say what to do; say it here, before the pull.
DIRTY=$(git diff --name-only; git diff --cached --name-only)
DIRTY=$(printf '%s\n' "$DIRTY" | sed '/^$/d' | sort -u)
# Nothing in REGENERATED is ever a local edit worth keeping — that is what being
# on that list means — so none of it can be a reason to stop. Excluding the whole
# list, rather than only what this run happened to touch, is what makes the check
# the same on the first run and the hundredth.
for f in "${REGENERATED[@]}"; do
  DIRTY=$(printf '%s\n' "$DIRTY" | grep -vxF "$f" || true)
done
if [ -n "$DIRTY" ]; then
  echo ""
  echo "  ✗  Update stopped — these tracked files have local changes:"
  printf '       %s\n' $DIRTY
  echo ""
  # Print the commands with the ACTUAL file list already in them. Telling
  # someone to substitute <file> for four paths they now have to retype, while
  # they are blocked, is how a stash gets run on three of the four.
  DIRTY_ARGS=$(printf '%s\n' $DIRTY | tr '\n' ' ')
  echo "     Keep them:"
  echo "       git stash push -- $DIRTY_ARGS"
  echo "     Discard them:"
  echo "       git checkout -- $DIRTY_ARGS"
  echo "     Then re-run:  bash update.sh"
  echo ""
  echo "     If a file above is listed with no visible difference, it is a"
  echo "     permission change, not an edit — check with:  git diff --summary"
  echo ""
  echo "     On the DEVELOPMENT machine this is expected — deploy_dev.sh rsyncs"
  echo "     source into this folder. Use deploy_dev.sh + push_github.sh there,"
  echo "     not update.sh."
  echo ""
  exit 1
fi

git pull origin main
echo "  ✓ Code updated"

# Restore job history if we had one
if [ -n "$JOBS_BACKUP" ]; then
  echo "$JOBS_BACKUP" > "$JOBS_FILE"
  echo "  ✓ Job history restored"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Rebuild frontend
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 2: Rebuilding frontend ─────────────────────────────"
cd frontend
# Load nvm if available (needed on machines that install Node via nvm)
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && source "$NVM_DIR/nvm.sh"
nvm use 20 >/dev/null 2>&1 || true
# npm ci, not npm install: "install" rewrites package-lock.json, which is
# tracked, so the next update.sh run finds a dirty tree and git pull aborts.
# "ci" installs exactly what the lockfile says and never writes to it. Fall back
# to install only when ci refuses (lockfile genuinely out of sync with
# package.json), and put the lockfile back afterwards so the tree stays clean.
if ! npm ci --silent 2>/dev/null; then
  echo "  ℹ  npm ci failed (lockfile out of sync) — falling back to npm install"
  npm install --silent
  git checkout -- package-lock.json 2>/dev/null || true
fi
npm run build
cd ..
echo "  ✓ Frontend rebuilt"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Ensure start_backend.sh exists and is executable
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 3: Checking helper scripts ────────────────────────"
if [[ ! -f "$SCRIPT_DIR/start_backend.sh" ]]; then
  cat > "$SCRIPT_DIR/start_backend.sh" <<'STARTSH'
#!/usr/bin/env bash
cd "$(dirname "$0")"
source venv/bin/activate
cd backend
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
STARTSH
  echo "  ✓ start_backend.sh created"
fi
chmod +x "$SCRIPT_DIR/start_backend.sh"
echo "  ✓ start_backend.sh ready"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Check if new R packages are needed
# ─────────────────────────────────────────────────────────────────────────────
echo "── Step 4: Checking R packages ─────────────────────────────"
Rscript backend/r_scripts/install_packages.R
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Update complete!                                       ║"
echo "║                                                          ║"
echo "║   Restart the backend to apply changes:                  ║"
echo "║     bash ~/r16s-app/start_backend.sh                    ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
