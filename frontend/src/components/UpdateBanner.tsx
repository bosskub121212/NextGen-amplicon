/**
 * UpdateBanner — polls /update/check on mount and shows a dismissible
 * notification when a new version is available.
 *
 * Public repo — no GitHub token required for checking or applying updates.
 */

import { useEffect, useState } from "react";

const API = ""; // relative URL — works on any host/port
const CHECK_INTERVAL_MS = 60 * 60 * 1000; // re-check every hour

interface UpdateInfo {
  available: boolean;
  current_version?: string;
  latest_version?: string;
  release_date?: string;
  changelog?: string;
  error?: string;
  message?: string;
}

interface ApplyResult {
  success: boolean;
  message?: string;
  steps?: { step: string; success: boolean; output: string }[];
}

export default function UpdateBanner() {
  const [info, setInfo]               = useState<UpdateInfo | null>(null);
  const [applying, setApplying]       = useState(false);
  const [applyResult, setApplyResult] = useState<ApplyResult | null>(null);
  const [dismissed, setDismissed]     = useState(false);
  const [checking, setChecking]       = useState(false);

  // ── Check version (no token needed — public repo) ─────────────
  const doCheck = async () => {
    setChecking(true);
    try {
      const data: UpdateInfo = await fetch(`${API}/update/check`).then(r => r.json());
      setInfo(data);
    } catch {
      // silently ignore — backend may not be up yet
    }
    setChecking(false);
  };

  useEffect(() => {
    doCheck();
    const iv = setInterval(doCheck, CHECK_INTERVAL_MS);
    return () => clearInterval(iv);
  }, []);

  // ── Apply update ──────────────────────────────────────────────
  const applyUpdate = async () => {
    if (!window.confirm(
      "This will run git pull + npm install + npm build.\n\nMake sure no analysis is running!\n\nContinue?"
    )) return;

    setApplying(true);
    setApplyResult(null);
    try {
      const res: ApplyResult = await fetch(`${API}/update/apply`, {
        method: "POST"
      }).then(r => r.json());
      setApplyResult(res);
    } catch {
      setApplyResult({ success: false, message: "Could not reach backend." });
    }
    setApplying(false);
  };

  // ── Nothing to show ───────────────────────────────────────────
  if (dismissed) return null;
  if (!info) return null;

  // ── Error banner ──────────────────────────────────────────────
  if (info.error) {
    return (
      <div className="update-banner update-banner--error">
        <span className="update-banner__icon">⚠️</span>
        <span className="update-banner__sub">
          Update check failed: {info.message || info.error}
        </span>
        <button className="ub-btn ub-btn--ghost"
          onClick={() => setDismissed(true)}>✕</button>
      </div>
    );
  }

  // ── Up to date — show nothing ─────────────────────────────────
  if (!info.available) return null;

  // ── Update-available banner ───────────────────────────────────
  return (
    <div className="update-banner update-banner--available">
      <div className="update-banner__left">
        <span className="update-banner__icon">🚀</span>
        <div>
          <span className="update-banner__title">
            Update available — v{info.latest_version}
          </span>
          <span className="update-banner__sub">
            Current: v{info.current_version}
            {info.release_date ? `  ·  Released ${info.release_date}` : ""}
          </span>
          {info.changelog && (
            <span className="update-banner__changelog">{info.changelog}</span>
          )}
        </div>
      </div>

      <div className="update-banner__right">
        {!applyResult ? (
          <>
            <button className="ub-btn ub-btn--primary"
              onClick={applyUpdate} disabled={applying}>
              {applying ? "⏳ Updating…" : "⬆️ Update Now"}
            </button>
            <button className="ub-btn ub-btn--ghost"
              onClick={() => setDismissed(true)}>Later</button>
          </>
        ) : (
          <div className="ub-apply-result">
            {applyResult.success ? (
              <span className="ub-apply-ok">
                ✅ Update complete — please restart the app
              </span>
            ) : (
              <span className="ub-apply-err">❌ {applyResult.message}</span>
            )}
            {applyResult.steps?.map((s, i) => (
              <div key={i} className={`ub-step ${s.success ? "ok" : "fail"}`}>
                <strong>{s.step}:</strong> {s.success ? "✓" : "✗"}
                {!s.success && s.output && (
                  <pre className="ub-step-out">{s.output.slice(-400)}</pre>
                )}
              </div>
            ))}
            <button className="ub-btn ub-btn--ghost"
              onClick={() => setDismissed(true)}>Dismiss</button>
          </div>
        )}
      </div>
    </div>
  );
}
