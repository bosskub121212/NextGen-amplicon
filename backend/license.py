"""
NextGen-Amplicon — Offline License System
==========================================
All validation is done locally using HMAC-SHA256.
No internet connection required after key activation.

Key format:  NGAMP-{BASE32_PAYLOAD}-{HMAC8}

Payload bytes:
  [0]      — pipeline bitmask (1 byte)
  [1..3]   — days since 2000-01-01 (uint24 big-endian)
  [4..15]  — machine ID bytes (12 bytes, optional)

Pipeline bit positions:
  bit 0 = 16S      bit 1 = 12S      bit 2 = ITS1     bit 3 = ITS2
  bit 4 = COX1     bit 5 = 18S-nema bit 6 = PacBio

DEV BYPASS:  Place ".dev_bypass" in the app root — always fully licensed.
"""

import base64
import hashlib
import hmac
import json
import os
import platform
import struct
from datetime import date, datetime, timedelta
from pathlib import Path

# ── Optional encryption (graceful fallback if not installed) ──────────────────
try:
    from cryptography.fernet import Fernet, InvalidToken
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    _CRYPTO_OK = True
except ImportError:
    _CRYPTO_OK = False

# ── Master switch ─────────────────────────────────────────────────────────────
LICENSE_ENABLED = True

# ── HMAC secret — MUST match tools/keygen.py ─────────────────────────────────
HMAC_SECRET = "NGAMP-OFFLINE-KEY-2026-CHANGE-BEFORE-PROD"

# ── Dev-bypass signing secret (internal only — never expose) ──────────────────
_DEV_BYPASS_SECRET = "NGAMP-DEV-INTERNAL-7f3a9c2e-BYPASS-SIGN"

# ── Cache encryption salt (fixed, non-secret) ─────────────────────────────────
_CACHE_SALT = b"ngamp_cache_v2_salt_2026"

# ── Pipeline registry ─────────────────────────────────────────────────────────
PIPELINES: list = [
    ("16S",       0),
    ("12S",       1),
    ("ITS1",      2),
    ("ITS2",      3),
    ("COX1",      4),
    ("18S-nema",  5),
    ("PacBio",    6),
]
ALL_MASK = (1 << len(PIPELINES)) - 1

PIPELINE_NAMES = [name for name, _ in PIPELINES]
EPOCH = date(2000, 1, 1)

# ── File paths ────────────────────────────────────────────────────────────────
APP_DIR    = Path(__file__).parent.parent
DEV_BYPASS = APP_DIR / ".dev_bypass"
CACHE_FILE = APP_DIR / ".license_cache.json"
EXPIRY_WARN_DAYS = 14


# ── Machine ID ────────────────────────────────────────────────────────────────

def get_machine_id() -> str:
    raw = f"{platform.node()}:{platform.machine()}:{platform.system()}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


# ── Low-level key codec ───────────────────────────────────────────────────────

def _payload_to_b32(b: bytes) -> str:
    return base64.b32encode(b).decode().rstrip("=")


def _b32_to_payload(s: str) -> bytes:
    pad = (8 - len(s) % 8) % 8
    return base64.b32decode(s.upper() + "=" * pad)


def _sign(payload_b32: str) -> str:
    return hmac.new(
        HMAC_SECRET.encode(), payload_b32.encode(), hashlib.sha256
    ).hexdigest()[:8].upper()


def _make_payload(pipe_mask: int, expiry: date, machine_id: str = "") -> bytes:
    days = (expiry - EPOCH).days
    b = bytes([pipe_mask & 0xFF]) + days.to_bytes(3, "big")
    if machine_id.strip():
        mid_hex = machine_id.strip().replace("-", "").replace(":", "").lower()
        mid_hex = (mid_hex + "0" * 24)[:24]
        try:
            b += bytes.fromhex(mid_hex)
        except ValueError:
            pass
    return b


def generate_key(pipe_mask: int, expiry: date, machine_id: str = "") -> str:
    b   = _make_payload(pipe_mask, expiry, machine_id)
    p   = _payload_to_b32(b)
    sig = _sign(p)
    return f"NGAMP-{p}-{sig}"


# ── Key validation ────────────────────────────────────────────────────────────

def _validate_and_decode(key: str) -> dict:
    try:
        parts = key.strip().upper().split("-")
        if len(parts) < 3 or parts[0] != "NGAMP":
            return {"ok": False, "error": "Invalid key format (bad prefix)"}

        sig_given   = parts[-1]
        payload_b32 = "".join(parts[1:-1])

        sig_expected = _sign(payload_b32)
        if not hmac.compare_digest(sig_given, sig_expected):
            return {"ok": False, "error": "Invalid key — signature mismatch"}

        raw = _b32_to_payload(payload_b32)
        if len(raw) < 4:
            return {"ok": False, "error": "Key payload too short"}

        pipe_mask = raw[0]
        days      = int.from_bytes(raw[1:4], "big")
        expiry    = EPOCH + timedelta(days=days)
        mid       = raw[4:16].hex() if len(raw) >= 16 else ""

        return {
            "ok":         True,
            "pipe_mask":  pipe_mask,
            "expiry":     expiry,
            "machine_id": mid,
            "error":      None,
        }
    except Exception as exc:
        return {"ok": False, "error": f"Key parse error: {exc}"}


def _mask_to_pipelines(mask: int) -> list:
    return [name for name, bit in PIPELINES if mask & (1 << bit)]


# ── Cache encryption helpers ──────────────────────────────────────────────────

def _get_cache_cipher():
    """Derive a Fernet cipher key from the machine ID. Requires cryptography lib."""
    if not _CRYPTO_OK:
        return None
    mid = get_machine_id().encode()
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=_CACHE_SALT,
        iterations=200_000,
    )
    key = base64.urlsafe_b64encode(kdf.derive(mid))
    return Fernet(key)


def _load_cache() -> dict:
    if not CACHE_FILE.exists():
        return {}
    raw = CACHE_FILE.read_bytes()
    # Try encrypted first
    if _CRYPTO_OK:
        try:
            cipher = _get_cache_cipher()
            decrypted = cipher.decrypt(raw)
            return json.loads(decrypted)
        except Exception:
            pass
    # Fallback: try plain JSON (old format — invalidate by returning empty)
    try:
        data = json.loads(raw.decode())
        # Old plain-text cache found — re-encrypt and save securely
        if data:
            _save_cache(data)
        return data
    except Exception:
        return {}


def _save_cache(data: dict):
    try:
        payload = json.dumps(data).encode()
        if _CRYPTO_OK:
            cipher = _get_cache_cipher()
            CACHE_FILE.write_bytes(cipher.encrypt(payload))
        else:
            # No crypto lib — store plain (degraded security)
            CACHE_FILE.write_bytes(payload)
    except Exception as exc:
        print(f"[license] Cache write failed: {exc}")


# ── Dev-bypass authentication ─────────────────────────────────────────────────

def _verify_dev_bypass() -> bool:
    """Check .dev_bypass contains a valid machine-bound HMAC signature."""
    if not DEV_BYPASS.exists():
        return False
    try:
        content = DEV_BYPASS.read_text().strip()
        # Format: NGAMP-DEV-{machine_id}-{sig8}
        parts = content.split("-")
        if len(parts) < 4 or parts[0] != "NGAMP" or parts[1] != "DEV":
            return False
        given_sig = parts[-1].upper()
        mid_and_prefix = "-".join(parts[:-1])
        expected_sig = hmac.new(
            _DEV_BYPASS_SECRET.encode(),
            mid_and_prefix.encode(),
            hashlib.sha256,
        ).hexdigest()[:8].upper()
        if not hmac.compare_digest(given_sig, expected_sig):
            print("[license] WARNING: .dev_bypass present but signature invalid — ignored")
            return False
        # Also verify machine ID matches this machine
        claimed_mid = parts[2] if len(parts) >= 3 else ""
        if claimed_mid and claimed_mid.lower() != get_machine_id().lower():
            print("[license] WARNING: .dev_bypass machine ID mismatch — ignored")
            return False
        return True
    except Exception:
        return False


def generate_dev_bypass_token() -> str:
    """Generate a signed .dev_bypass token for this machine (run by developer)."""
    mid = get_machine_id()
    prefix = f"NGAMP-DEV-{mid}"
    sig = hmac.new(
        _DEV_BYPASS_SECRET.encode(),
        prefix.encode(),
        hashlib.sha256,
    ).hexdigest()[:8].upper()
    return f"{prefix}-{sig}"


# ── Public API ────────────────────────────────────────────────────────────────

def check_license() -> dict:
    if _verify_dev_bypass():
        return {
            "status":          "dev",
            "message":         "Developer machine — license bypass active.",
            "pipelines":       PIPELINE_NAMES,
            "days_remaining":  None,
            "expiry_warning":  False,
            "machine_bound":   False,
        }

    if not LICENSE_ENABLED:
        return {
            "status":          "disabled",
            "message":         "License system is not yet active.",
            "pipelines":       PIPELINE_NAMES,
            "days_remaining":  None,
            "expiry_warning":  False,
            "machine_bound":   False,
        }

    cache = _load_cache()
    if not cache:
        return {
            "status":          "no_license",
            "message":         "No license found. Please enter your license key.",
            "pipelines":       [],
            "days_remaining":  None,
            "expiry_warning":  False,
            "machine_bound":   False,
        }

    key = cache.get("license_key", "")
    decoded = _validate_and_decode(key)
    if not decoded["ok"]:
        return {
            "status":          "invalid",
            "message":         f"Cached key is invalid: {decoded['error']}",
            "pipelines":       [],
            "days_remaining":  None,
            "expiry_warning":  False,
            "machine_bound":   False,
        }

    bound_mid = decoded["machine_id"]
    if bound_mid:
        this_mid = get_machine_id()
        if this_mid != bound_mid:
            return {
                "status":          "invalid",
                "message":         "This key is bound to a different machine.",
                "pipelines":       [],
                "days_remaining":  None,
                "expiry_warning":  False,
                "machine_bound":   True,
            }

    expiry    = decoded["expiry"]
    today     = date.today()
    if today > expiry:
        return {
            "status":          "expired",
            "message":         f"License expired on {expiry.strftime('%d %b %Y')}.",
            "pipelines":       [],
            "days_remaining":  0,
            "expiry_warning":  False,
            "machine_bound":   bool(bound_mid),
        }

    days_left = (expiry - today).days
    warn      = days_left <= EXPIRY_WARN_DAYS
    pipes     = _mask_to_pipelines(decoded["pipe_mask"])

    return {
        "status":          "active",
        "message":         f"License active — {days_left} day(s) remaining.",
        "pipelines":       pipes,
        "days_remaining":  days_left,
        "expiry_warning":  warn,
        "expiry_date":     expiry.isoformat(),
        "machine_bound":   bool(bound_mid),
    }


def get_license_pipelines() -> list:
    if _verify_dev_bypass() or not LICENSE_ENABLED:
        return PIPELINE_NAMES
    status = check_license()
    return status.get("pipelines", [])


def activate_license(license_key: str) -> dict:
    if not LICENSE_ENABLED:
        return {
            "success": False,
            "message": "License system is not active.",
            "status":  check_license(),
        }

    if not license_key or not license_key.strip():
        return {"success": False, "message": "Key cannot be empty.", "status": check_license()}

    decoded = _validate_and_decode(license_key.strip())
    if not decoded["ok"]:
        return {"success": False, "message": decoded["error"], "status": check_license()}

    if date.today() > decoded["expiry"]:
        return {
            "success": False,
            "message": f"This key expired on {decoded['expiry'].strftime('%d %b %Y')}.",
            "status":  check_license(),
        }

    bound_mid = decoded["machine_id"]
    if bound_mid:
        this_mid = get_machine_id()
        if this_mid != bound_mid:
            return {
                "success": False,
                "message": (
                    "This key is bound to a different machine.\n"
                    f"Key machine ID:  {bound_mid}\n"
                    f"This machine ID: {this_mid}"
                ),
                "status": check_license(),
            }

    pipes = _mask_to_pipelines(decoded["pipe_mask"])
    cache = {
        "license_key":  license_key.strip().upper(),
        "machine_id":   get_machine_id(),
        "pipelines":    pipes,
        "expiry_date":  decoded["expiry"].isoformat(),
        "activated_at": datetime.utcnow().isoformat(),
    }
    _save_cache(cache)

    status = check_license()
    return {
        "success": True,
        "message": (
            f"License activated!  "
            f"Expires {decoded['expiry'].strftime('%d %b %Y')}  |  "
            f"Pipelines: {', '.join(pipes)}"
        ),
        "status": status,
    }


def deactivate_license() -> dict:
    if CACHE_FILE.exists():
        CACHE_FILE.unlink()
    return {"success": True, "message": "License removed from this machine."}


def is_pipeline_allowed(marker: str) -> bool:
    allowed = get_license_pipelines()
    marker_upper = marker.upper()
    for name in allowed:
        if name.upper() == marker_upper:
            return True
        if marker_upper == "ITS" and name.upper().startswith("ITS"):
            return True
    return False
