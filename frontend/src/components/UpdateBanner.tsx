/**
 * UpdateBanner — polls /update/check on mount and shows a dismissible
 * notification when a new version is available.
 *
 * Also shows a token-setup prompt if no GitHub token is configured.
 */

import { useEffect, useState, useRef } from "react";

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
  const [info, setInfo]             = useState<UpdateInfo | null>(null);
  const [tokenNeeded, setTokenNeeded] = useState(false);
  const [tokenInput, setTokenInput] = useState("");
  const [tokenSaving, setTokenSaving] = useState(false);
  const [tokenMsg, setTokenMsg]     = useState("");
  const [applying, setApplying]     = useState(false);
  const [applyResult, setApplyResult] = useState<ApplyResult | null>(null);
  const [dismissed, setDismissed]   = useState(false);
  const [showTokenForm, setShowTokenForm] = useState(false);

  // ── Check token then version ──────────────────────────────────
  const doCheck = async () => {
    try {
      const ts = await fetch(`${API}/update/token-status`).then(r => r.json());
      if (!ts.configured) {
        setTokenNeeded(true);
        return;
      }
      setTokenNeeded(false);
      const data: UpdateInfo = await fetch(`${API}/update/check`).then(r => r.json());
      setInfo(data);
    } catch {
      // silently ignore network errors — backend might not be up yet
    }
  };

  useEffect(() => {
    doCheck();
    const iv = setInterval(doCheck, CHECK_INTERVAL_MS);
    return () => clearInterval(iv);
  }, []);

  // ── Save token ────────────────────────────────────────────────
  const saveToken = async () => {
    if (!tokenInput.trim()) return;
    setTokenSaving(true);
    setTokenMsg("");
    try {
      const res = await fetch(`${API}/update/save-token`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ token: tokenInput.trim() })
      }).then(r => r.json());
      if (res.ok || res.success) {
        setTokenMsg("✅ Token saved! Checking for updates...");
        setTokenInput("");
        setShowTokenForm(false);
        setTimeout(() => { setTokenMsg(""); doCheck(); }, 1500);
      } else {
        setTokenMsg(`❌ ${res.message}`);
      }
    } catch {
      setTokenMsg("❌ Could not reach backend.");
    }
    setTokenSaving(false);
  };

  // ── Apply update ──────────────────────────────────────────────
  const applyUpdate = async () => {
    if (!window.confirm(
      "This will run git pull and npm install.\n\nMake sure no analysis is running!\n\nContinue?"
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

  // ── Nothing to show ──────────────────────────────────────────
  if (dismissed) return null;
  if (!tokenNeeded && (!info || (!info.available && !info.error))) return null;

  // ── Token-needed banner ───────────────────────────────────────
  if (tokenNeeded) {
    return (
      <div className="update-banner update-banner--token">
        <div className="update-banner__left">
          <span className="update-banner__icon">🔑</span>
          <div>
            <span className="update-banner__title">Auto-update setup needed</span>
            <span className="update-banner__sub">
              Add a GitHub token to enable automatic updates
            </span>
          </div>
        </div>
        <div className="update-banner__right">
          {!showTokenForm ? (
            <>
              <button className="ub-btn ub-btn--primary"
                onClick={() => setShowTokenForm(true)}>
                Setup Token
              </button>
              <button className="ub-btn ub-btn--ghost"
                onClick={() => setDismissed(true)}>✕</button>
            </>
          ) : (
            <div className="ub-token-form">
              <input
                className="ub-token-input"
                type="password"
                placeholder="ghp_xxxxxxxxxxxx"
                value={tokenInput}
                onChange={e => setTokenInput(e.target.value)}
                onKeyDown={e => e.key === "Enter" && saveToken()}
              />
              <button className="ub-btn ub-btn--primary"
                onClick={saveToken} disabled={tokenSaving}>
                {tokenSaving ? "Saving…" : "Save"}
              </button>
              <button className="ub-btn ub-btn--ghost"
                onClick={() => setShowTokenForm(false)}>Cancel</button>
              {tokenMsg && <span className="ub-token-msg">{tokenMsg}</span>}
            </div>
          )}
        </div>
      </div>
    );
  }

  // ── Error banner ──────────────────────────────────────────────
  if (info?.error) {
    if (info.error === "no_token") return null; // already handled above
    return (
      <div className="update-banner update-banner--error">
        <span className="update-banner__icon">⚠️</span>
        <span className="update-banner__sub">{info.message}</span>
        <button className="ub-btn ub-btn--ghost" onClick={() => setDismissed(true)}>✕</button>
      </div>
    );
  }

  // ── Update-available banner ───────────────────────────────────
  if (!info?.available) return null;

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
              {applying ? "Updating…" : "Update Now"}
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
                  <pre className="ub-step-out">{s.output}</pre>
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
