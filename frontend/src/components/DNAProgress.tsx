import { useState, useEffect, useRef } from "react";

interface Props {
  percent: number;
  logs: string[];
  currentStep?: string;
}

const GREEN = "#10b981";
const BLUE  = "#3b82f6";
const GRAY  = "#cbd5e1";

export default function DNAProgress({ percent, logs, currentStep }: Props) {
  const [phase, setPhase] = useState(0);
  const animRef = useRef<ReturnType<typeof setInterval> | null>(null);

  useEffect(() => {
    animRef.current = setInterval(() => {
      setPhase((p) => (p + 3) % 360);
    }, 40);
    return () => { if (animRef.current) clearInterval(animRef.current); };
  }, []);

  // ── Horizontal helix dimensions ──────────────────────────────
  const W  = 380;
  const H  = 82;
  const cy = H / 2;
  const A  = 28;          // amplitude (half-height of wave)
  const wl = 58;          // wavelength in pixels
  const rad = (phase * Math.PI) / 180;

  // Strand path — runs left→right along X axis
  const strandPath = (offset: number) => {
    const pts: string[] = [];
    for (let i = 0; i <= 240; i++) {
      const x = (W * i) / 240;
      const y = cy + A * Math.sin((2 * Math.PI * x) / wl + offset);
      pts.push(`${i === 0 ? "M" : "L"} ${x.toFixed(1)} ${y.toFixed(1)}`);
    }
    return pts.join(" ");
  };

  // Rungs — vertical crossbars between the two strands
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

      {/* ── Horizontal DNA Helix SVG ── */}
      <div className="dna-helix-row">
        <svg
          width="100%" viewBox={`0 0 ${W} ${H}`}
          style={{
            display: "block",
            filter: `drop-shadow(0 0 8px #10b98133) drop-shadow(0 0 8px #3b82f633)`,
          }}
        >
          <defs>
            {/* Fill clip: left portion up to fillX */}
            <clipPath id="dna-h-clip">
              <rect x={0} y={0} width={fillX} height={H} />
            </clipPath>
          </defs>

          {/* ── Gray (unfilled) layer ── */}
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

          {/* ── Colored (filled) layer — clipped left→right ── */}
          <g clipPath="url(#dna-h-clip)">
            {/* Strand 1 — GREEN */}
            <path d={strandPath(rad)}
              stroke={GREEN} strokeWidth={3.5} fill="none" strokeLinecap="round" />
            {/* Strand 2 — BLUE */}
            <path d={strandPath(rad + Math.PI)}
              stroke={BLUE} strokeWidth={3.5} fill="none" strokeLinecap="round" />

            {/* Rungs — split color at center */}
            {rungs.map((r, i) => {
              const cTop = i % 2 === 0 ? GREEN : BLUE;
              const cBot = i % 2 === 0 ? BLUE  : GREEN;
              const mid  = (r.y1 + r.y2) / 2;
              return (
                <g key={i}>
                  <line x1={r.x} y1={r.y1} x2={r.x} y2={mid} stroke={cTop} strokeWidth={2.5} />
                  <line x1={r.x} y1={mid}  x2={r.x} y2={r.y2} stroke={cBot} strokeWidth={2.5} />
                  <circle cx={r.x} cy={r.y1} r={4} fill={cTop} />
                  <circle cx={r.x} cy={r.y2} r={4} fill={cBot} />
                </g>
              );
            })}
          </g>

          {/* Fill frontier indicator */}
          {percent > 0 && percent < 100 && (
            <line
              x1={fillX} y1={4} x2={fillX} y2={H - 4}
              stroke="#ffffff" strokeWidth={1.5} strokeDasharray="3,3" opacity={0.6}
            />
          )}
        </svg>
      </div>

      {/* ── % + Step label ── */}
      <div className="dna-info-row">
        <div
          className="dna-pct"
          style={{
            background: `linear-gradient(90deg, ${GREEN}, ${BLUE})`,
            WebkitBackgroundClip: "text",
            WebkitTextFillColor: "transparent",
          }}
        >
          {Math.round(percent)}%
        </div>
        {currentStep && (
          <div
            className="dna-step-label"
            style={{
              background: `linear-gradient(90deg, ${GREEN}, ${BLUE})`,
              WebkitBackgroundClip: "text",
              WebkitTextFillColor: "transparent",
            }}
          >
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
          <span className="dna-log-title">pipeline log</span>
        </div>
        <div className="dna-log-body">
          {logs.length === 0 ? (
            <div className="dna-log-line dim">Waiting for output...</div>
          ) : (
            logs.slice(-3).map((line, i, arr) => (
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
