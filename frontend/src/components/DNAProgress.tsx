import { useState, useEffect, useRef } from "react";

type ThemeId = "default" | "orange" | "purple";

interface Props {
  percent:      number;
  logs:         string[];
  currentStep?: string;
  theme?:       ThemeId;
}

// Strand colors per theme
const THEME_COLORS: Record<ThemeId, { a: string; b: string; glow1: string; glow2: string }> = {
  default: { a: "#10b981", b: "#3b82f6", glow1: "#10b98133", glow2: "#3b82f633" },
  orange:  { a: "#ea580c", b: "#9ca3af", glow1: "#ea580c33", glow2: "#9ca3af33" },
  purple:  { a: "#7c3aed", b: "#1e1040", glow1: "#7c3aed44", glow2: "#7c3aed22" },
};

export default function DNAProgress({ percent, logs, currentStep, theme = "default" }: Props) {
  const [phase, setPhase] = useState(0);
  const animRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const logBodyRef = useRef<HTMLDivElement>(null);
  const stickToBottom = useRef(true);

  useEffect(() => {
    animRef.current = setInterval(() => setPhase(p => (p + 3) % 360), 40);
    return () => { if (animRef.current) clearInterval(animRef.current); };
  }, []);

  // Auto-scroll to the newest log line, but only if the user hasn't
  // scrolled up to read earlier output — don't yank them back down.
  useEffect(() => {
    const el = logBodyRef.current;
    if (el && stickToBottom.current) el.scrollTop = el.scrollHeight;
  }, [logs]);

  const handleLogScroll = () => {
    const el = logBodyRef.current;
    if (!el) return;
    stickToBottom.current = el.scrollHeight - el.scrollTop - el.clientHeight < 24;
  };

  const { a: COL_A, b: COL_B, glow1, glow2 } = THEME_COLORS[theme];
  const GRAY = theme === "purple" ? "#312e81" : "#cbd5e1";

  const pctGrad = theme === "orange"
    ? "linear-gradient(90deg,#ea580c,#f97316)"
    : theme === "purple"
    ? "linear-gradient(90deg,#7c3aed,#a78bfa)"
    : "linear-gradient(90deg,#10b981,#3b82f6)";

  // ── Helix geometry ──────────────────────────────────────────
  const W = 380, H = 82, cy = H / 2, A = 28, wl = 58;
  const rad = (phase * Math.PI) / 180;

  const strandPath = (offset: number) => {
    const pts: string[] = [];
    for (let i = 0; i <= 240; i++) {
      const x = (W * i) / 240;
      const y = cy + A * Math.sin((2 * Math.PI * x) / wl + offset);
      pts.push(`${i === 0 ? "M" : "L"} ${x.toFixed(1)} ${y.toFixed(1)}`);
    }
    return pts.join(" ");
  };

  const rungs: { x: number; y1: number; y2: number }[] = [];
  for (let x = wl / 4; x < W; x += wl / 2) {
    rungs.push({
      x,
      y1: cy + A * Math.sin((2 * Math.PI * x) / wl + rad),
      y2: cy + A * Math.sin((2 * Math.PI * x) / wl + rad + Math.PI),
    });
  }

  const fillX = W * Math.min(percent, 100) / 100;

  return (
    <div className="dna-wrap">

      {/* ── DNA Helix ── */}
      <div className="dna-helix-row">
        <svg width="100%" viewBox={`0 0 ${W} ${H}`}
          style={{ display: "block", filter: `drop-shadow(0 0 8px ${glow1}) drop-shadow(0 0 8px ${glow2})` }}>
          <defs>
            <clipPath id="dna-h-clip">
              <rect x={0} y={0} width={fillX} height={H} />
            </clipPath>
          </defs>

          {/* Gray unfilled layer */}
          <g opacity={0.2}>
            <path d={strandPath(rad)}           stroke={GRAY} strokeWidth={3.5} fill="none" strokeLinecap="round" />
            <path d={strandPath(rad + Math.PI)} stroke={GRAY} strokeWidth={3.5} fill="none" strokeLinecap="round" />
            {rungs.map((r, i) => (
              <g key={i}>
                <line x1={r.x} y1={r.y1} x2={r.x} y2={r.y2} stroke={GRAY} strokeWidth={2.5} />
                <circle cx={r.x} cy={r.y1} r={3.5} fill={GRAY} />
                <circle cx={r.x} cy={r.y2} r={3.5} fill={GRAY} />
              </g>
            ))}
          </g>

          {/* Colored filled layer */}
          <g clipPath="url(#dna-h-clip)">
            <path d={strandPath(rad)}           stroke={COL_A} strokeWidth={3.5} fill="none" strokeLinecap="round" />
            <path d={strandPath(rad + Math.PI)} stroke={COL_B} strokeWidth={3.5} fill="none" strokeLinecap="round" />
            {rungs.map((r, i) => {
              const cTop = i % 2 === 0 ? COL_A : COL_B;
              const cBot = i % 2 === 0 ? COL_B : COL_A;
              const mid  = (r.y1 + r.y2) / 2;
              return (
                <g key={i}>
                  <line x1={r.x} y1={r.y1} x2={r.x} y2={mid}  stroke={cTop} strokeWidth={2.5} />
                  <line x1={r.x} y1={mid}  x2={r.x} y2={r.y2} stroke={cBot} strokeWidth={2.5} />
                  <circle cx={r.x} cy={r.y1} r={4} fill={cTop} />
                  <circle cx={r.x} cy={r.y2} r={4} fill={cBot} />
                </g>
              );
            })}
          </g>

          {/* Progress frontier line */}
          {percent > 0 && percent < 100 && (
            <line x1={fillX} y1={4} x2={fillX} y2={H - 4}
              stroke="white" strokeWidth={1.5} strokeDasharray="3,3" opacity={0.6} />
          )}
        </svg>
      </div>

      {/* ── % + step label ── */}
      <div className="dna-info-row">
        <div className="dna-pct"
          style={{ background: pctGrad, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" }}>
          {Math.round(percent)}%
        </div>
        {currentStep && (
          <div className="dna-step-label"
            style={{ background: pctGrad, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" }}>
            {currentStep}
          </div>
        )}
      </div>

      {/* ── Terminal log ── */}
      <div className="dna-log-box">
        <div className="dna-log-header">
          <span className="dna-log-dot red" />
          <span className="dna-log-dot yellow" />
          <span className="dna-log-dot green" />
          <span className="dna-log-title">pipeline log{logs.length > 0 ? ` (${logs.length} lines)` : ""}</span>
        </div>
        <div className="dna-log-body" ref={logBodyRef} onScroll={handleLogScroll}>
          {logs.length === 0 ? (
            <div className="dna-log-line dim">Waiting for output...</div>
          ) : (
            logs.map((line, i, arr) => (
              <div key={i} className={`dna-log-line ${i === arr.length - 1 ? "latest" : "dim"}`}>
                <span className="dna-log-prompt">$</span> {line}
              </div>
            ))
          )}
          <span className="dna-log-cursor">▋</span>
        </div>
      </div>
    </div>
  );
}
