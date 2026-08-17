"""
NextGen-Amplicon — Auto Update System
Checks GitHub (public repo) for new versions and applies updates.
No GitHub token required — repo is public.
"""
import json, subprocess, os
from pathlib import Path
from urllib import request, error
from urllib.request import Request

REPO       = "bosskub121212/NextGen-amplicon"
REPO_URL   = f"https://github.com/{REPO}.git"
APP_DIR    = Path(__file__).parent.parent
VER_FILE   = APP_DIR / "version.json"
TOKEN_FILE = Path.home() / ".config" / "amplicon" / "github_token"

# ── Token helpers (optional — only needed for private repos) ──
def get_token() -> str:
    if TOKEN_FILE.exists():
        return TOKEN_FILE.read_text().strip()
    return os.getenv("GITHUB_TOKEN", "")

def save_token(token: str):
    TOKEN_FILE.parent.mkdir(parents=True, exist_ok=True)
    TOKEN_FILE.write_text(token.strip())
    TOKEN_FILE.chmod(0o600)

# ── Version helpers ───────────────────────────────────────────
def get_local_version() -> dict:
    try:
        return json.loads(VER_FILE.read_text())
    except Exception:
        return {"version": "0.0.0", "release_date": "", "changelog": ""}

def _fetch_remote_version() -> dict:
    """Fetch version.json from public GitHub repo — no token needed."""
    # Use raw.githubusercontent.com for public repos (no auth required)
    url = f"https://raw.githubusercontent.com/{REPO}/main/version.json"
    req = Request(url, headers={
        "User-Agent": "NextGen-Amplicon-Updater/1.0",
        "Cache-Control": "no-cache",
    })
    with request.urlopen(req, timeout=8) as resp:
        return json.loads(resp.read().decode())

# ── Public API ────────────────────────────────────────────────
def check_update() -> dict:
    """Compare local version.json with remote. Returns update info."""
    try:
        remote    = _fetch_remote_version()
        local     = get_local_version()
        available = remote.get("version", "0") != local.get("version", "0")
        return {
            "available":       available,
            "current_version": local.get("version"),
            "latest_version":  remote.get("version"),
            "release_date":    remote.get("release_date", ""),
            "changelog":       remote.get("changelog", ""),
        }
    except error.HTTPError as e:
        # 429 (rate-limited) and 5xx (transient GitHub/CDN hiccup) are not
        # something the user can act on and not a real problem with the app —
        # skip silently instead of showing a scary red error banner. Only
        # surface genuinely actionable errors (404 = repo/branch moved, etc).
        if e.code == 429 or e.code >= 500:
            return {"available": False}
        return {"available": False, "error": f"http_{e.code}", "message": f"GitHub error {e.code}"}
    except Exception as e:
        # Network hiccups (DNS, timeout, offline) are equally not actionable —
        # skip silently too rather than nagging on every transient blip.
        return {"available": False}


def apply_update() -> dict:
    """Run git pull + npm install + npm build. Returns step-by-step results."""
    steps = []

    # ── Set public HTTPS remote (no token needed) ─────────────
    subprocess.run(
        ["git", "remote", "set-url", "origin", REPO_URL],
        cwd=APP_DIR, capture_output=True
    )

    # ── git pull ──────────────────────────────────────────────
    r = subprocess.run(
        ["git", "pull", "origin", "main"],
        cwd=APP_DIR, capture_output=True, text=True, timeout=60
    )
    steps.append({
        "step":    "git pull",
        "success": r.returncode == 0,
        "output":  (r.stdout + r.stderr).strip()
    })
    if r.returncode != 0:
        return {"success": False, "steps": steps}

    # ── npm install (pick up new packages) ────────────────────
    npm_cmd = "source ~/.nvm/nvm.sh 2>/dev/null || true; npm install --legacy-peer-deps"
    r2 = subprocess.run(
        ["bash", "-lc", npm_cmd],
        cwd=APP_DIR / "frontend", capture_output=True, text=True, timeout=120
    )
    steps.append({
        "step":    "npm install",
        "success": r2.returncode == 0,
        "output":  (r2.stdout + r2.stderr).strip()[-300:]
    })

    # ── npm build ─────────────────────────────────────────────
    r3 = subprocess.run(
        ["bash", "-lc", "source ~/.nvm/nvm.sh 2>/dev/null || true; npm run build"],
        cwd=APP_DIR / "frontend", capture_output=True, text=True, timeout=180
    )
    steps.append({
        "step":    "npm build",
        "success": r3.returncode == 0,
        "output":  (r3.stdout + r3.stderr).strip()[-300:]
    })

    return {
        "success": True,
        "steps":   steps,
        "message": "Update complete — please restart the app to apply changes."
    }
