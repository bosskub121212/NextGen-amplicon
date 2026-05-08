/**
 * LicenseModal — Offline HMAC license activation UI
 *
 * Props:
 *   required  — when true the modal cannot be dismissed (first-run / expired).
 *               The user MUST enter a valid key before using the app.
 *   onClose   — called when the user closes the modal (only possible when !required)
 *   onActivated — called with the new status after a successful activation
 */

import { useEffect, useState } from "react";

const API = ""; // relative URL — works on any host/port

// ── Types ─────────────────────────────────────────────────────────────────────
export interface LicenseStatus {
  status: "active" | "expired" | "invalid" | "no_license" | "disabled" | "dev";
  message?: string;
  pipelines?: string[];
  days_remaining?: number | null;
  expiry_warning?: boolean;
  expiry_date?: string;
  machine_bound?: boolean;
}

interface Props {
  required?:     boolean;              // force-show — no close button
  onClose?:      () => void;
  onActivated?:  (status: LicenseStatus) => void;
}

// Pipeline display order
const ALL_PIPELINES = ["16S", "12S", "ITS1", "ITS2", "COX1", "18S-nema", "PacBio"];

// ── Component ─────────────────────────────────────────────────────────────────
export default function LicenseModal({ required = false, onClose, onActivated }: Props) {
  const [status,       setStatus]      = useState<LicenseStatus | null>(null);
  const [keyInput,     setKeyInput]    = useState("");
  const [loading,      setLoading]     = useState(false);
  const [deactivating, setDeactivating] = useState(false);
  const [msg,          setMsg]         = useState<{ text: string; ok: boolean } | null>(null);
  const [tab,          setTab]         = useState<"activate" | "info">("activate");

  // ── Fetch status on mount ──────────────────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/license/status`)
      .then(r => r.json())
      .then((s: LicenseStatus) => {
        setStatus(s);
        if (s.status === "active" || s.status === "dev") setTab("info");
      })
      .catch(() =>
        setStatus({ status: "no_license", message: "Could not reach backend." })
      );
  }, []);

  // ── Activate ───────────────────────────────────────────────────────────────
  const activate = async () => {
    if (!keyInput.trim()) return;
    setLoading(true);
    setMsg(null);
    try {
      const res = await fetch(`${API}/license/activate`, {
        method:  "POST",
        headers: { "Content-Type": "application/json" },
        body:    JSON.stringify({ license_key: keyInput.trim() }),
      }).then(r => r.json());

      if (res.success) {
        setMsg({ text: res.message, ok: true });
        const newStatus: LicenseStatus = res.status;
        setStatus(newStatus);
        setTab("info");
        onActivated?.(newStatus);
      } else {
        setMsg({ text: res.message, ok: false });
      }
    } catch {
      setMsg({ text: "❌ Could not reach backend.", ok: false });
    }
    setLoading(false);
  };

  // ── Deactivate ─────────────────────────────────────────────────────────────
  const deactivate = async () => {
    if (!window.confirm(
      "This will remove the license from this machine.\n" +
      "You can then activate it on another machine.\n\nContinue?"
    )) return;
    setDeactivating(true);
    setMsg(null);
    try {
      await fetch(`${API}/license/deactivate`, { method: "POST" });
      const newStatus: LicenseStatus = { status: "no_license" };
      setStatus(newStatus);
      setTab("activate");
      setMsg({ text: "License removed from this machine.", ok: true });
    } catch {
      setMsg({ text: "❌ Could not reach backend.", ok: false });
    }
    setDeactivating(false);
  };

  // ── Helpers ────────────────────────────────────────────────────────────────
  const isActive  = status?.status === "active" || status?.status === "dev";
  const isDev     = status?.status === "dev";
  const isExpired = status?.status === "expired";
  const canClose  = !required && (isActive || isDev);

  const statusIcon =
    isDev      ? "🛠️" :
    isActive   ? "✅" :
    isExpired  ? "⏰" :
                 "🔒";

  const statusColor: Record<string, string> = {
    active:     "#16a34a",
    dev:        "#0891b2",
    expired:    "#dc2626",
    invalid:    "#dc2626",
    no_license: "#6b7280",
    disabled:   "#6b7280",
  };

  return (
    <div className="popup-overlay" style={{ zIndex: 9999 }}>
      <div className="license-modal" style={{ maxWidth: 520, width: "95vw" }}>

        {/* ── Header ─────────────────────────────────────────────────────── */}
        <div className="license-modal__header">
          <h2>🔐 License</h2>
          {canClose && onClose && (
            <button className="license-modal__close" onClick={onClose}>✕</button>
          )}
        </div>

        {/* ── Status block ────────────────────────────────────────────────── */}
        <div
          className={`license-status-block license-status-block--${status?.status ?? "loading"}`}
          style={{ borderLeft: `4px solid ${statusColor[status?.status ?? "no_license"] ?? "#6b7280"}` }}
        >
          {!status ? (
            <span style={{ color: "#6b7280" }}>Checking license…</span>
          ) : (
            <>
              <div className="license-status-row">
                <span className="license-status-icon">{statusIcon}</span>
                <div style={{ flex: 1 }}>
                  <div className="license-status-label" style={{ color: statusColor[status.status] }}>
                    {isDev      ? "Developer Mode"         :
                     isActive   ? "License Active"         :
                     isExpired  ? "License Expired"        :
                     status.status === "no_license" ? "No License" :
                     status.status === "disabled"   ? "License Disabled" :
                                   "Invalid License"}
                  </div>
                  {status.message && (
                    <div className="license-status-msg">{status.message}</div>
                  )}
                </div>
                {isActive && !isDev && status.days_remaining != null && (
                  <div className={`license-days ${status.expiry_warning ? "license-days--warn" : ""}`}>
                    {status.expiry_warning && "⚠️ "}
                    {status.days_remaining}d left
                  </div>
                )}
              </div>

              {/* ── Pipeline badges ─────────────────────────────────────── */}
              {isActive && status.pipelines && status.pipelines.length > 0 && (
                <div className="license-pipelines" style={{ marginTop: 10 }}>
                  <span style={{ fontSize: 12, color: "#6b7280", marginRight: 6 }}>
                    Pipelines:
                  </span>
                  {ALL_PIPELINES.map(p => {
                    const enabled = status.pipelines!.includes(p);
                    return (
                      <span
                        key={p}
                        className={`license-pipe-badge ${enabled ? "license-pipe-badge--on" : "license-pipe-badge--off"}`}
                      >
                        {enabled ? "✓" : "✗"} {p}
                      </span>
                    );
                  })}
                </div>
              )}

              {/* ── Expiry date ─────────────────────────────────────────── */}
              {isActive && !isDev && status.expiry_date && (
                <div style={{ fontSize: 12, color: "#6b7280", marginTop: 6 }}>
                  Expires: {new Date(status.expiry_date).toLocaleDateString("en-GB", {
                    day: "numeric", month: "long", year: "numeric"
                  })}
                  {status.machine_bound && "  ·  Machine-bound key"}
                </div>
              )}
            </>
          )}
        </div>

        {/* ── Tabs (only when both make sense) ────────────────────────────── */}
        {status && (
          <div className="license-tabs">
            <button
              className={`license-tab ${tab === "activate" ? "license-tab--active" : ""}`}
              onClick={() => setTab("activate")}
            >
              🔑 Enter Key
            </button>
            <button
              className={`license-tab ${tab === "info" ? "license-tab--active" : ""}`}
              onClick={() => setTab("info")}
            >
              ℹ️ Info
            </button>
          </div>
        )}

        {/* ── Tab: Activate ─────────────────────────────────────────────── */}
        {tab === "activate" && (
          <div className="license-activate-form">
            {required && !isActive && (
              <div className="license-required-notice">
                ⚠️ A valid license key is required to use this application.
                Please enter your key below.
              </div>
            )}

            <label className="license-label">License Key</label>
            <input
              className="license-key-input"
              type="text"
              placeholder="NGAMP-XXXXXXXX-XXXXXXXX"
              value={keyInput}
              onChange={e => setKeyInput(e.target.value.toUpperCase())}
              onKeyDown={e => e.key === "Enter" && activate()}
              spellCheck={false}
              autoComplete="off"
            />
            <button
              className="ub-btn ub-btn--primary"
              onClick={activate}
              disabled={loading || !keyInput.trim()}
              style={{ marginTop: 8 }}
            >
              {loading ? "Activating…" : "Activate License"}
            </button>

            {msg && (
              <div
                className="license-activate-msg"
                style={{ color: msg.ok ? "#16a34a" : "#dc2626" }}
              >
                {msg.text}
              </div>
            )}
          </div>
        )}

        {/* ── Tab: Info / Management ────────────────────────────────────── */}
        {tab === "info" && (
          <div className="license-info-tab">
            {msg && (
              <div
                className="license-activate-msg"
                style={{ color: msg.ok ? "#16a34a" : "#dc2626", marginBottom: 8 }}
              >
                {msg.text}
              </div>
            )}

            <div className="license-machine-id">
              <span style={{ fontSize: 12, color: "#6b7280" }}>This machine's ID:</span>
              <MachineIdDisplay api={API} />
            </div>

            {(isActive && !isDev) && (
              <button
                className="ub-btn ub-btn--ghost"
                onClick={deactivate}
                disabled={deactivating}
                style={{ marginTop: 12 }}
              >
                {deactivating ? "Removing…" : "🔄 Transfer to another machine"}
              </button>
            )}

            {!isActive && (
              <p style={{ fontSize: 12, color: "#6b7280", marginTop: 10 }}>
                Switch to the "Enter Key" tab to activate a license.
              </p>
            )}
          </div>
        )}

      </div>
    </div>
  );
}

// ── Small helper: async machine ID fetch ───────────────────────────────────
function MachineIdDisplay({ api }: { api: string }) {
  const [mid, setMid] = useState<string>("…");
  useEffect(() => {
    fetch(`${api}/license/machine-id`)
      .then(r => r.json())
      .then(d => setMid(d.machine_id ?? "unknown"))
      .catch(() => setMid("(could not reach backend)"));
  }, [api]);
  return (
    <code style={{
      fontFamily: "Consolas, monospace", fontSize: 12,
      background: "#f1f5f9", padding: "2px 6px", borderRadius: 4,
      marginLeft: 6, userSelect: "all",
    }}>
      {mid}
    </code>
  );
}
