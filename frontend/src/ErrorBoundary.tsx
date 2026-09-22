import React from "react";

/**
 * Catches a render error anywhere below it and shows what went wrong instead of
 * an empty page.
 *
 * React 18 unmounts the whole tree when a render throws, so before this existed
 * one undefined field took down the entire app: a pipeline emitted a new kind of
 * checkpoint whose payload had no `merged_pct`, the warning modal called
 * `.toFixed()` on it, and every screen — job list, history, settings, all of it —
 * went blank. Nothing on screen said why, and the app looked dead rather than
 * like one broken panel.
 *
 * A boundary cannot prevent that bug, but it turns a blank page into a message
 * with the error and a way back, which is the difference between "the app is
 * broken" and "this panel is broken".
 */
interface Props  { children: React.ReactNode; label?: string }
interface State  { error: Error | null }

export default class ErrorBoundary extends React.Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State { return { error }; }

  componentDidCatch(error: Error, info: React.ErrorInfo) {
    // Keep the detail in the console for whoever opens devtools next.
    console.error(`[${this.props.label ?? "app"}] render error:`, error, info.componentStack);
  }

  render() {
    const { error } = this.state;
    if (!error) return this.props.children;
    return (
      <div className="eb-wrap">
        <div className="eb-card">
          <div className="eb-icon">⚠️</div>
          <h2 className="eb-title">
            {this.props.label ? `${this.props.label} failed to render` : "Something went wrong"}
          </h2>
          <p className="eb-msg">{error.message || String(error)}</p>
          <pre className="eb-stack">{(error.stack || "").split("\n").slice(0, 8).join("\n")}</pre>
          <div className="eb-actions">
            <button className="eb-btn eb-btn-primary" onClick={() => this.setState({ error: null })}>
              ↺ Try again
            </button>
            <button className="eb-btn" onClick={() => window.location.reload()}>
              ⟳ Reload the app
            </button>
          </div>
          <p className="eb-hint">
            The rest of your data is safe — nothing on the server changed.
            Running jobs keep going.
          </p>
        </div>
      </div>
    );
  }
}
