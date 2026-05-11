import { useState, useEffect } from "react";

export type ThemeId = "default" | "orange" | "purple";

interface ThemeOption {
  id: ThemeId;
  label: string;
  sub: string;
  gradient: string;
  btn: string;
  accent: string;
}

const THEMES: ThemeOption[] = [
  {
    id: "default",
    label: "Indigo / Blue",
    sub: "Classic (default)",
    gradient: "linear-gradient(135deg, #1e1b4b 0%, #4338ca 100%)",
    btn: "#4f46e5",
    accent: "#818cf8",
  },
  {
    id: "orange",
    label: "Orange / Gray",
    sub: "Warm & bold",
    gradient: "linear-gradient(135deg, #1c1917 0%, #ea580c 100%)",
    btn: "#ea580c",
    accent: "#fb923c",
  },
  {
    id: "purple",
    label: "Black / Purple",
    sub: "Dark & rich",
    gradient: "linear-gradient(135deg, #0f0a1e 0%, #6d28d9 100%)",
    btn: "#7c3aed",
    accent: "#a78bfa",
  },
];

interface UpdateInfo {
  current: string;
  latest: string | null;
  has_update: boolean;
  checking: boolean;
  error: string | null;
  updating: boolean;
  update_done: boolean;
}

interface Props {
  onClose: () => void;
  theme: ThemeId;
  onTheme: (t: ThemeId) => void;
}

export default function SettingsPanel({ onClose, theme, onTheme }: Props) {
  const [update, setUpdate] = useState<UpdateInfo>({
    current: "—", latest: null,
    has_update: false, checking: false,
    error: null, updating: false, update_done: false,
  });

  // Load current version on mount
  useEffect(() => {
    fetch("/update/check")
      .then(r => r.ok ? r.json() : null)
      .then(d => {
        if (!d) return;
        setUpdate(prev => ({
          ...prev,
          current: d.current_version ?? "—",
          latest: d.latest_version ?? null,
          has_update: !!d.update_available,
        }));
      })
      .catch(() => {});
  }, []);

  const checkUpdate = () => {
    setUpdate(prev => ({ ...prev, checking: true, error: null }));
    fetch("/update/check")
      .then(r => r.ok ? r.json() : null)
      .then(d => {
        if (!d) { setUpdate(prev => ({ ...prev, checking: false, error: "Could not reach server" })); return; }
        setUpdate(prev => ({
          ...prev,
          checking: false,
          current: d.current_version ?? prev.current,
          latest: d.latest_version ?? null,
          has_update: !!d.update_available,
        }));
      })
      .catch(() => setUpdate(prev => ({ ...prev, checking: false, error: "Network error" })));
  };

  const doUpdate = () => {
    setUpdate(prev => ({ ...prev, updating: true, error: null }));
    fetch("/update/apply", { method: "POST" })
      .then(r => r.ok ? r.json() : Promise.reject())
      .then(() => setUpdate(prev => ({ ...prev, updating: false, update_done: true })))
      .catch(() => setUpdate(prev => ({ ...prev, updating: false, error: "Update failed — check logs" })));
  };

  return (
    <div className="sp-overlay" onClick={e => { if (e.target === e.currentTarget) onClose(); }}>
      <div className="sp-modal">

        {/* Header */}
        <div className="sp-header">
          <div className="sp-header-left">
            <span className="sp-header-icon">⚙️</span>
            <div>
              <div className="sp-header-title">Settings</div>
              <div className="sp-header-sub">NextGen-Amplicon preferences</div>
            </div>
          </div>
          <button className="sp-close" onClick={onClose}>✕</button>
        </div>

        <div className="sp-body">

          {/* ── Section 1: Theme ───────────────────────────────────── */}
          <div className="sp-section">
            <div className="sp-section-title">🎨 Appearance</div>
            <p className="sp-section-hint">Choose a colour theme for the interface</p>
            <div className="sp-theme-grid">
              {THEMES.map(t => (
                <button
                  key={t.id}
                  className={`sp-theme-card ${theme === t.id ? "sp-theme-card--active" : ""}`}
                  onClick={() => onTheme(t.id)}
                >
                  {/* Mini gradient preview */}
                  <div className="sp-theme-preview" style={{ background: t.gradient }}>
                    <div className="sp-theme-preview-btn" style={{ background: t.btn }} />
                    <div className="sp-theme-preview-chip" style={{ background: t.accent }} />
                    <div className="sp-theme-preview-bar" />
                  </div>
                  <div className="sp-theme-label">{t.label}</div>
                  <div className="sp-theme-sub">{t.sub}</div>
                  {theme === t.id && <div className="sp-theme-check">✓</div>}
                </button>
              ))}
            </div>
          </div>

          {/* ── Section 2: Updates ─────────────────────────────────── */}
          <div className="sp-section">
            <div className="sp-section-title">🔄 Software Updates</div>
            <div className="sp-update-box">
              <div className="sp-update-row">
                <span className="sp-update-label">Current version</span>
                <span className="sp-ver-badge">{update.current}</span>
              </div>
              {update.latest && (
                <div className="sp-update-row">
                  <span className="sp-update-label">Latest available</span>
                  <span className={`sp-ver-badge ${update.has_update ? "sp-ver-badge--new" : ""}`}>
                    {update.latest}
                  </span>
                </div>
              )}
              {update.error && (
                <div className="sp-update-error">⚠ {update.error}</div>
              )}
              {update.update_done && (
                <div className="sp-update-ok">
                  ✅ Update applied — please restart the app
                </div>
              )}
            </div>

            <div className="sp-update-actions">
              <button
                className="sp-btn sp-btn--ghost"
                onClick={checkUpdate}
                disabled={update.checking || update.updating}
              >
                {update.checking ? "⏳ Checking…" : "🔍 Check for Update"}
              </button>
              {update.has_update && !update.update_done && (
                <button
                  className="sp-btn sp-btn--primary"
                  onClick={doUpdate}
                  disabled={update.updating}
                >
                  {update.updating ? "⏳ Updating…" : "⬆️ Update Now"}
                </button>
              )}
            </div>
          </div>

          {/* ── Section 3: About ───────────────────────────────────── */}
          <div className="sp-section sp-section--about">
            <div className="sp-section-title">ℹ️ About</div>
            <div className="sp-about-logo">🧬</div>
            <div className="sp-about-name">NextGen-Amplicon</div>
            <div className="sp-about-ver">Version {update.current}</div>
            <div className="sp-about-company">
              NEXTGEN NETWORK CORPORATION COMPANY LIMITED
            </div>
            <div className="sp-about-copy">© {new Date().getFullYear()} All rights reserved</div>
          </div>

        </div>
      </div>
    </div>
  );
}
