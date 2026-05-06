"""
NextGen-Amplicon — License Key Generator
=========================================
Desktop GUI tool (Python + tkinter) for generating offline HMAC-SHA256
license keys.

Key format:  NGAMP-{BASE32_PAYLOAD}-{HMAC8}

Payload bytes (4 or 16 bytes):
  [0]      — pipeline bitmask (1 byte)
  [1..3]   — days since 2000-01-01 (uint24 big-endian)
  [4..15]  — machine ID bytes (12 bytes, only when binding to a machine)

Pipeline bits:
  bit 0 = 16S
  bit 1 = 12S
  bit 2 = ITS1
  bit 3 = ITS2
  bit 4 = COX1
  bit 5 = 18S-nema
  bit 6 = PacBio

Usage:
  python keygen.py
  (or double-click on Windows)
"""

import tkinter as tk
from tkinter import ttk, messagebox
import base64
import hashlib
import hmac
import platform
from datetime import date, timedelta

# ── HMAC secret — MUST match license.py on the backend ───────────────────────
DEFAULT_SECRET = "NGAMP-OFFLINE-KEY-2026-CHANGE-BEFORE-PROD"

# ── Pipeline definitions ──────────────────────────────────────────────────────
PIPELINES = [
    ("16S",       0),
    ("12S",       1),
    ("ITS1",      2),
    ("ITS2",      3),
    ("COX1",      4),
    ("18S-nema",  5),
    ("PacBio",    6),
]
ALL_MASK = (1 << len(PIPELINES)) - 1   # 0b1111111 = 127

EPOCH = date(2000, 1, 1)


# ── Key generation helpers ────────────────────────────────────────────────────

def _make_payload_bytes(pipe_mask: int, expiry: date, machine_id: str = "") -> bytes:
    """Encode pipeline mask + expiry (+ optional machine binding) into bytes."""
    days = (expiry - EPOCH).days
    b = bytes([pipe_mask & 0xFF]) + days.to_bytes(3, "big")
    if machine_id.strip():
        # Normalise: hex string, pad/truncate to 24 chars (= 12 bytes)
        mid_hex = machine_id.strip().replace("-", "").replace(":", "").lower()
        mid_hex = (mid_hex + "0" * 24)[:24]
        try:
            b += bytes.fromhex(mid_hex)
        except ValueError:
            pass   # ignore bad machine ID
    return b


def _payload_to_b32(b: bytes) -> str:
    """Base32-encode without padding '=' characters."""
    return base64.b32encode(b).decode().rstrip("=")


def _sign(payload_b32: str, secret: str) -> str:
    """HMAC-SHA256 of the payload, returning first 8 hex chars (uppercase)."""
    return hmac.new(
        secret.encode(), payload_b32.encode(), hashlib.sha256
    ).hexdigest()[:8].upper()


def generate_key(pipe_mask: int, expiry: date,
                 machine_id: str = "", secret: str = DEFAULT_SECRET) -> str:
    """Return a formatted NGAMP-… key string."""
    b   = _make_payload_bytes(pipe_mask, expiry, machine_id)
    p   = _payload_to_b32(b)
    sig = _sign(p, secret)
    return f"NGAMP-{p}-{sig}"


def get_machine_id_local() -> str:
    """Return the machine ID of the current computer (same algo as license.py)."""
    raw = f"{platform.node()}:{platform.machine()}:{platform.system()}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


# ── GUI ───────────────────────────────────────────────────────────────────────

class KeygenApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("🔑 NextGen-Amplicon — License Key Generator")
        self.resizable(False, False)
        self._build_ui()
        self._on_expiry_change()

    # ── UI construction ────────────────────────────────────────────────────────
    def _build_ui(self):
        pad = {"padx": 12, "pady": 6}

        # ── Header ────────────────────────────────────────────────────────────
        hdr = tk.Frame(self, bg="#1e293b")
        hdr.pack(fill="x")
        tk.Label(
            hdr, text="🔑  NextGen-Amplicon  –  License Key Generator",
            font=("Segoe UI", 14, "bold"), fg="#f8fafc", bg="#1e293b",
            pady=14, padx=16
        ).pack(side="left")

        # ── Main container ────────────────────────────────────────────────────
        main = ttk.Frame(self, padding=16)
        main.pack(fill="both")

        # Section helper
        def section(parent, title):
            f = ttk.LabelFrame(parent, text=title, padding=10)
            f.pack(fill="x", pady=(0, 10))
            return f

        # ── HMAC Secret ──────────────────────────────────────────────────────
        sec_frm = section(main, "🔐  HMAC Secret  (must match license.py on backend)")
        self._secret_var = tk.StringVar(value=DEFAULT_SECRET)
        self._show_secret = tk.BooleanVar(value=False)
        sec_row = ttk.Frame(sec_frm)
        sec_row.pack(fill="x")
        self._secret_entry = ttk.Entry(
            sec_row, textvariable=self._secret_var,
            show="•", width=52, font=("Consolas", 10)
        )
        self._secret_entry.pack(side="left", fill="x", expand=True)
        ttk.Button(
            sec_row, text="👁", width=3,
            command=self._toggle_secret
        ).pack(side="left", padx=(6, 0))

        # ── Pipelines ─────────────────────────────────────────────────────────
        pip_frm = section(main, "🧬  Pipelines enabled in this key")
        self._pipe_vars = {}
        grid = ttk.Frame(pip_frm)
        grid.pack(fill="x")
        for i, (name, bit) in enumerate(PIPELINES):
            v = tk.BooleanVar(value=True)
            self._pipe_vars[bit] = v
            ttk.Checkbutton(grid, text=name, variable=v, width=10).grid(
                row=i // 4, column=i % 4, sticky="w", padx=4, pady=2
            )
        btn_row = ttk.Frame(pip_frm)
        btn_row.pack(fill="x", pady=(6, 0))
        ttk.Button(btn_row, text="Select All",
                   command=lambda: [v.set(True)  for v in self._pipe_vars.values()]
                   ).pack(side="left", padx=(0, 6))
        ttk.Button(btn_row, text="Clear All",
                   command=lambda: [v.set(False) for v in self._pipe_vars.values()]
                   ).pack(side="left")

        # ── Expiry ────────────────────────────────────────────────────────────
        exp_frm = section(main, "📅  Expiry")
        self._expiry_mode = tk.StringVar(value="trial")
        modes = [
            ("7-day Trial",  "trial"),
            ("1 Year",       "1y"),
            ("2 Years",      "2y"),
            ("Custom date:", "custom"),
        ]
        for text, val in modes:
            rb = ttk.Radiobutton(
                exp_frm, text=text, variable=self._expiry_mode,
                value=val, command=self._on_expiry_change
            )
            rb.pack(anchor="w")

        # Custom date entry + preview
        cust_row = ttk.Frame(exp_frm)
        cust_row.pack(fill="x", pady=(4, 0))
        ttk.Label(cust_row, text="  Date (YYYY-MM-DD):").pack(side="left")
        self._custom_date_var = tk.StringVar(
            value=(date.today() + timedelta(days=365)).strftime("%Y-%m-%d")
        )
        self._custom_entry = ttk.Entry(
            cust_row, textvariable=self._custom_date_var, width=14,
            font=("Consolas", 10), state="disabled"
        )
        self._custom_entry.pack(side="left", padx=(6, 0))

        self._expiry_preview = ttk.Label(
            exp_frm, text="", font=("Segoe UI", 9), foreground="#64748b"
        )
        self._expiry_preview.pack(anchor="w", pady=(4, 0))

        # ── Machine binding ───────────────────────────────────────────────────
        mid_frm = section(main, "🖥️  Machine binding  (optional — leave blank for floating key)")
        mid_row = ttk.Frame(mid_frm)
        mid_row.pack(fill="x")
        ttk.Label(mid_row, text="Machine ID:").pack(side="left")
        self._machine_var = tk.StringVar()
        ttk.Entry(
            mid_row, textvariable=self._machine_var, width=36,
            font=("Consolas", 10)
        ).pack(side="left", padx=(8, 8), fill="x", expand=True)
        ttk.Button(
            mid_row, text="Use This Machine",
            command=self._fill_this_machine
        ).pack(side="left")
        ttk.Label(
            mid_frm,
            text="When set, the key only works on the machine with that ID.",
            font=("Segoe UI", 8), foreground="#94a3b8"
        ).pack(anchor="w", pady=(4, 0))

        # ── Generate button ───────────────────────────────────────────────────
        ttk.Button(
            main, text="🔑   Generate Key",
            command=self._generate, style="Accent.TButton"
        ).pack(fill="x", ipady=8, pady=(4, 10))

        # ── Output ────────────────────────────────────────────────────────────
        out_frm = section(main, "📋  Generated Key")
        self._key_var = tk.StringVar()
        key_out = ttk.Entry(
            out_frm, textvariable=self._key_var, state="readonly",
            font=("Consolas", 12), foreground="#1e40af"
        )
        key_out.pack(fill="x", ipady=6)

        out_btns = ttk.Frame(out_frm)
        out_btns.pack(fill="x", pady=(8, 0))
        ttk.Button(
            out_btns, text="📋  Copy to Clipboard",
            command=self._copy
        ).pack(side="left", padx=(0, 8))
        self._info_label = ttk.Label(
            out_btns, text="", font=("Segoe UI", 9), foreground="#16a34a"
        )
        self._info_label.pack(side="left")

        # ── Footer ────────────────────────────────────────────────────────────
        ttk.Label(
            main,
            text="⚠  Keep the HMAC secret private — it controls key validity.",
            font=("Segoe UI", 8), foreground="#ef4444"
        ).pack()

    # ── Event handlers ─────────────────────────────────────────────────────────

    def _toggle_secret(self):
        self._show_secret.set(not self._show_secret.get())
        self._secret_entry.config(show="" if self._show_secret.get() else "•")

    def _on_expiry_change(self):
        mode = self._expiry_mode.get()
        self._custom_entry.config(
            state="normal" if mode == "custom" else "disabled"
        )
        try:
            exp = self._calc_expiry()
            self._expiry_preview.config(
                text=f"Expires: {exp.strftime('%d %B %Y')}  "
                     f"({(exp - date.today()).days} days from today)"
            )
        except Exception:
            self._expiry_preview.config(text="")

    def _fill_this_machine(self):
        mid = get_machine_id_local()
        self._machine_var.set(mid)

    def _calc_expiry(self) -> date:
        mode = self._expiry_mode.get()
        today = date.today()
        if mode == "trial":
            return today + timedelta(days=7)
        elif mode == "1y":
            return today.replace(year=today.year + 1)
        elif mode == "2y":
            return today.replace(year=today.year + 2)
        else:
            return date.fromisoformat(self._custom_date_var.get().strip())

    def _get_pipe_mask(self) -> int:
        mask = 0
        for bit, var in self._pipe_vars.items():
            if var.get():
                mask |= (1 << bit)
        return mask

    def _generate(self):
        self._info_label.config(text="")
        # Validate pipelines
        mask = self._get_pipe_mask()
        if mask == 0:
            messagebox.showerror("No Pipelines", "Select at least one pipeline.")
            return

        # Validate / compute expiry
        try:
            expiry = self._calc_expiry()
        except ValueError:
            messagebox.showerror("Invalid Date",
                                 "Custom date must be YYYY-MM-DD format.")
            return
        if expiry <= date.today():
            messagebox.showerror("Invalid Expiry",
                                 "Expiry date must be in the future.")
            return

        # Validate secret
        secret = self._secret_var.get().strip()
        if not secret:
            messagebox.showerror("No Secret", "HMAC secret cannot be empty.")
            return

        key = generate_key(mask, expiry, self._machine_var.get(), secret)
        self._key_var.set(key)

        # Build summary label
        pipe_names = [n for n, b in PIPELINES if mask & (1 << b)]
        days = (expiry - date.today()).days
        summary = f"Pipelines: {', '.join(pipe_names)}   |   Expires in {days} days"
        self._info_label.config(text=summary, foreground="#64748b")

    def _copy(self):
        key = self._key_var.get()
        if not key:
            messagebox.showinfo("Empty", "Generate a key first.")
            return
        self.clipboard_clear()
        self.clipboard_append(key)
        self._info_label.config(text="✅ Copied to clipboard!", foreground="#16a34a")
        self.after(2500, lambda: self._info_label.config(text=""))


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    try:
        # Create root first so we can configure styles BEFORE building the UI
        root = tk.Tk()
        root.withdraw()   # hide briefly while we set up styles

        # Apply a modern ttk theme
        style = ttk.Style(root)
        for theme in ("vista", "winnative", "clam", "alt", "default"):
            if theme in style.theme_names():
                style.theme_use(theme)
                break

        # Define Accent.TButton BEFORE KeygenApp builds the UI
        style.configure(
            "Accent.TButton",
            font=("Segoe UI", 11, "bold"),
            foreground="#ffffff",
            background="#4f46e5",
        )
        style.map(
            "Accent.TButton",
            background=[("active", "#4338ca"), ("!disabled", "#4f46e5")],
            foreground=[("!disabled", "#ffffff")],
        )

        root.destroy()   # discard the temp root

        app = KeygenApp()
        app.mainloop()

    except Exception:
        import traceback
        # Show error in a messagebox so it doesn't just vanish
        try:
            tk.Tk().withdraw()
            messagebox.showerror(
                "KeyGen Startup Error",
                traceback.format_exc()
            )
        except Exception:
            print(traceback.format_exc())
            input("Press Enter to exit...")
