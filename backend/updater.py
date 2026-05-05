"""
NextGen-Amplicon — Auto Update System
Checks GitHub (private repo) for new versions and applies updates.
"""
import json, subprocess, os
from pathlib import Path
from urllib import request, error
from urllib.request import Request

REPO      = "bosskub121212/NextGen-amplicon"
APP_DIR   = Path(__file__).parent.parent
VER_FILE  = APP_DIR / "version.json"
TOKEN_FILE = Path.home() / ".config" / "amplicon" / "github_token"

# ── Token helpers ─────────────────────────────────────────────
def get_token() -> str:
    """Read GitHub PAT from ~/.config/amplicon/github_token"""
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

def _fetch_remote_version(token: str) -> dict:
    url = f"https://api.github.com/repos/{REPO}/contents/version.json"
    req = Request(url, headers={
        "Authorization": f"token {token}",
        "Accept":        "application/vnd.github.v3.raw",
        "User-Agent":    "NextGen-Amplicon-Updater/1.0",
    })
    with request.urlopen(req, timeout=8) as resp:
        return json.loads(resp.read().decode())

# ── Public API ────────────────────────────────────────────────
def check_update() -> dict:
    """Compare local version.json with remote. Returns update info."""
    token = get_token()
    if not token:
        return {
            "available": False,
            "error": "no_token",
            "message": "GitHub token not configured — run setup_update.sh"
        }
    try:
        remote  = _fetch_remote_version(token)
        local   = get_local_version()
        available = remote.get("version", "0") != local.get("version", "0")
        return {
            "available":       available,
            "current_version": local.get("version"),
            "latest_version":  remote.get("version"),
            "release_date":    remote.get("release_date", ""),
            "changelog":       remote.get("changelog", ""),
        }
    except error.HTTPError as e:
        code = e.code
        msg  = "Invalid token" if code == 401 else f"GitHub API error {code}"
        return {"available": False, "error": f"http_{code}", "message": msg}
    except Exception as e:
        return {"available": False, "error": "network", "message": str(e)}


def apply_update() -> dict:
    """Run git pull + npm install. Returns step-by-step results."""
    token = get_token()
    steps = []

    # ── Configure git credentials (HTTPS with token) ──────────
    if token:
        creds_url = f"https://{token}@github.com"
        subprocess.run(
            ["git", "config", "credential.helper",
             f"!echo password={token}; echo username=x-token"],
            cwd=APP_DIR, capture_output=True
        )
        # Simpler: embed token in remote URL
        subprocess.run(
            ["git", "remote", "set-url", "origin",
             f"https://{token}@github.com/{REPO}.git"],
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

    # ── npm build (only if dist/ already exists → production mode) ─
    if (APP_DIR / "frontend" / "dist").exists():
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
