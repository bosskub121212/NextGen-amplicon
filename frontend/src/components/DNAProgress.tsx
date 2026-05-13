import { useState, useEffect, useRef, useMemo } from "react";

type ThemeId = "default" | "orange" | "purple";

interface Props {
  percent:      number;
  logs:         string[];
  currentStep?: string;
  theme?:       ThemeId;
}

interface Asteroid {
  id:       number;
  y:        number;
  size:     number;
  fromLeft: boolean;
  duration: number;
  born:     number;
  rotDir:   number;
}

/* ═══════════════════════════════════════════════════════════════
   SPACE PROGRESS  (orange = rocket → moon | purple = UFO → Earth)
═══════════════════════════════════════════════════════════════ */
function SpaceProgress({ percent, isOrange }: { percent: number; isOrange: boolean }) {
  const W = 380, H = 90;
  const DEST_X   = W - 28;
  const ORIGIN_X = 28;
  const shipX    = ORIGIN_X + (Math.min(percent, 100) / 100) * (DEST_X - ORIGIN_X - 44);
  const arrived  = percent >= 100;

  const [tick,      setTick]      = useState(0);
  const [asteroids, setAsteroids] = useState<Asteroid[]>([]);
  const aidRef = useRef(0);

  // Main animation tick — 80 ms
  useEffect(() => {
    const id = setInterval(() => setTick(t => t + 1), 80);
    return () => clearInterval(id);
  }, []);

  // Asteroid spawner — first after 3-8 s, then every 20-30 s
  useEffect(() => {
    let timer: ReturnType<typeof setTimeout>;
    const spawnAsteroid = () => {
      const now = Date.now();
      setAsteroids(prev => {
        const active = prev.filter(a => now - a.born < a.duration + 300);
        return [...active, {
          id:       aidRef.current++,
          y:        8 + Math.random() * (H - 20),
          size:     5 + Math.random() * 7,
          fromLeft: Math.random() > 0.5,
          duration: 2800 + Math.random() * 2400,
          born:     now,
          rotDir:   Math.random() > 0.5 ? 1 : -1,
        }];
      });
    };
    const scheduleNext = () => {
      timer = setTimeout(() => { spawnAsteroid(); scheduleNext(); },
        20000 + Math.random() * 10000);
    };
    timer = setTimeout(() => { spawnAsteroid(); scheduleNext(); },
      3000 + Math.random() * 5000);
    return () => clearTimeout(timer);
  }, []);

  // Fixed star positions (deterministic so they don't jump on re-render)
  const stars = useMemo(() => Array.from({ length: 26 }, (_, i) => ({
    x:  (i * 14.7  + 11) % W,
    y:  (i * 9.3   +  7) % H,
    r:  0.6 + (i % 3) * 0.35,
    op: 0.3 + (i % 5) * 0.13,
  })), []);

  const now          = Date.now();
  const flameFlicker = 0.72 + 0.28 * Math.sin(tick * 0.88);
  const beamFlicker  = 0.42 + 0.58 * Math.abs(Math.sin(tick * 0.36));
  const g = isOrange ? "o" : "p";   // gradient ID prefix

  return (
    <svg width="100%" viewBox={`0 0 ${W} ${H}`} style={{ display: "block" }}>
      <defs>
        {/* Moon shade */}
        <radialGradient id={`${g}-ms`}>
          <stop offset="0%"   stopColor="transparent" />
          <stop offset="100%" stopColor="rgba(0,0,0,.42)" />
        </radialGradient>
        {/* Earth blue */}
        <radialGradient id={`${g}-eb`} cx="44%" cy="44%">
          <stop offset="0%"   stopColor="#1976d2" />
          <stop offset="100%" stopColor="#0d47a1" />
        </radialGradient>
        {/* Earth shade */}
        <radialGradient id={`${g}-es`}>
          <stop offset="0%"   stopColor="transparent" />
          <stop offset="100%" stopColor="rgba(0,0,0,.3)" />
        </radialGradient>
        {/* Rocket flame */}
        <linearGradient id={`${g}-fl`} x1="0%" y1="50%" x2="100%" y2="50%">
          <stop offset="0%"   stopColor="#fbbf24" stopOpacity="0" />
          <stop offset="60%"  stopColor="#f97316" stopOpacity=".95" />
          <stop offset="100%" stopColor="#fef3c7" />
        </linearGradient>
        {/* Ship body gloss */}
        <radialGradient id={`${g}-sb`} cx="64%" cy="28%">
          <stop offset="0%"   stopColor="rgba(255,255,255,.28)" />
          <stop offset="100%" stopColor="rgba(0,0,0,.22)" />
        </radialGradient>
        {/* Saucer gloss */}
        <radialGradient id={`${g}-sc`} cx="50%" cy="20%">
          <stop offset="0%"   stopColor="rgba(255,255,255,.32)" />
          <stop offset="100%" stopColor="rgba(0,0,0,.28)" />
        </radialGradient>
        {/* Tractor beam */}
        <linearGradient id={`${g}-bm`} x1="50%" y1="0%" x2="50%" y2="100%">
          <stop offset="0%"   stopColor="#c4b5fd" stopOpacity=".95" />
          <stop offset="100%" stopColor="#c4b5fd" stopOpacity="0" />
        </linearGradient>
      </defs>

      {/* ── Space background ── */}
      <rect width={W} height={H} rx={10}
        fill={isOrange ? "#0d0a07" : "#060310"} />

      {/* ── Stars ── */}
      {stars.map((s, i) => (
        <circle key={i} cx={s.x} cy={s.y} r={s.r} fill="white" opacity={s.op} />
      ))}

      {/* ── Track ── */}
      <line x1={ORIGIN_X} y1={H / 2} x2={DEST_X - 20} y2={H / 2}
        stroke={isOrange ? "#ea580c" : "#7c3aed"}
        strokeWidth={1.2} strokeDasharray="3,8" opacity={0.28} />

      {/* ── Progress fill ── */}
      <line x1={ORIGIN_X} y1={H / 2}
        x2={Math.min(shipX + 5, DEST_X - 42)} y2={H / 2}
        stroke={isOrange ? "#f97316" : "#8b5cf6"}
        strokeWidth={2.5} strokeLinecap="round" opacity={0.72} />

      {/* ══ Destination ══ */}
      {isOrange ? (
        /* ─ Moon ─ */
        <g transform={`translate(${DEST_X},${H / 2})`}>
          <circle r={15} fill="#c0c0c0" />
          <circle cx={-4}  cy={-5} r={3.5}  fill="#999" opacity={0.55} />
          <circle cx={5}   cy={3}  r={2.5}  fill="#999" opacity={0.45} />
          <circle cx={-5}  cy={6}  r={1.8}  fill="#999" opacity={0.4} />
          <circle r={15} fill={`url(#${g}-ms)`} />
          {arrived && [0,45,90,135,180,225,270,315].map(a => (
            <line key={a}
              x1={Math.cos(a*Math.PI/180)*17} y1={Math.sin(a*Math.PI/180)*17}
              x2={Math.cos(a*Math.PI/180)*(22+3*Math.abs(Math.sin(tick*.22)))}
              y2={Math.sin(a*Math.PI/180)*(22+3*Math.abs(Math.sin(tick*.22)))}
              stroke="#fbbf24" strokeWidth={1.8} opacity={0.82} />
          ))}
        </g>
      ) : (
        /* ─ Earth ─ */
        <g transform={`translate(${DEST_X},${H / 2})`}>
          <circle r={16} fill={`url(#${g}-eb)`} />
          <ellipse cx={-4} cy={-5} rx={6}   ry={4}   fill="#2e7d32" opacity={0.95} />
          <ellipse cx={5}  cy={2}  rx={5}   ry={3}   fill="#388e3c" opacity={0.9} />
          <ellipse cx={-6} cy={5}  rx={3}   ry={2}   fill="#1b5e20" opacity={0.85} />
          <ellipse cx={3}  cy={-8} rx={2.5} ry={1.8} fill="#43a047" opacity={0.7} />
          <circle r={16} fill={`url(#${g}-es)`} />
          {arrived && [0,60,120,180,240,300].map(a => (
            <line key={a}
              x1={Math.cos(a*Math.PI/180)*18} y1={Math.sin(a*Math.PI/180)*18}
              x2={Math.cos(a*Math.PI/180)*(25+4*Math.abs(Math.sin(tick*.25)))}
              y2={Math.sin(a*Math.PI/180)*(25+4*Math.abs(Math.sin(tick*.25)))}
              stroke="#a78bfa" strokeWidth={2} opacity={0.85} />
          ))}
        </g>
      )}

      {/* ══ Asteroids ══ */}
      {asteroids.map(a => {
        const t = (now - a.born) / a.duration;
        if (t > 1.06) return null;
        const ax = a.fromLeft
          ? -a.size*3 + t*(W + a.size*6)
          :  W+a.size*3 - t*(W + a.size*6);
        const op  = Math.min(t*5,1) * Math.min((1-t)*5,1);
        const rot = t * 190 * a.rotDir;
        const c1 = isOrange ? "#92400e" : "#3b0764";
        const c2 = isOrange ? "#b45309" : "#5b21b6";
        return (
          <g key={a.id} transform={`translate(${ax},${a.y})`} opacity={op}>
            <g transform={`rotate(${rot})`}>
              <ellipse rx={a.size} ry={a.size*.78} fill={c1} />
              <ellipse cx={a.size*.3}  cy={-a.size*.25}
                rx={a.size*.65} ry={a.size*.5}  fill={c2} opacity={0.62} />
              <ellipse cx={-a.size*.35} cy={a.size*.3}
                rx={a.size*.38} ry={a.size*.28} fill="rgba(255,255,255,.06)" />
            </g>
          </g>
        );
      })}

      {/* ══ Rocket (orange) ══ */}
      {isOrange && (
        <g transform={`translate(${shipX},${H / 2})`}>
          {/* Flame */}
          {!arrived && <>
            <ellipse cx={-18} cy={0}
              rx={10*flameFlicker} ry={5.5*flameFlicker}
              fill={`url(#${g}-fl)`} />
            <ellipse cx={-14} cy={0}
              rx={4.8*flameFlicker} ry={2.2}
              fill="#fef9c3" opacity={0.78} />
          </>}
          {/* Body */}
          <ellipse rx={17} ry={9} fill="#ea580c" />
          <ellipse rx={17} ry={9} fill={`url(#${g}-sb)`} />
          {/* Nose cone */}
          <path d="M17,0 L10,-5 L10,5Z" fill="#c2410c" />
          <ellipse cx={13} cy={0} rx={5.5} ry={4.5} fill="#c2410c" />
          {/* Cockpit window */}
          <circle cx={3} cy={-.5} r={5} fill="#bfdbfe" opacity={0.88} />
          <ellipse cx={1.5} cy={-1.5} rx={2} ry={1.5} fill="white" opacity={0.5} />
          {/* Top fin */}
          <path d="M-2,-9 L-10,-18 L-16,-9Z" fill="#9a3412" />
          {/* Bottom fin */}
          <path d="M-2,9 L-10,18 L-16,9Z" fill="#9a3412" />
          {/* Highlight stripe */}
          <ellipse cx={3} cy={-4} rx={10} ry={2.5} fill="white" opacity={0.13} />
        </g>
      )}

      {/* ══ UFO (purple) ══ */}
      {!isOrange && (
        <g transform={`translate(${shipX},${H / 2})`}>
          {/* Tractor beam */}
          {!arrived && (
            <polygon points="0,8 -19,34 19,34"
              fill={`url(#${g}-bm)`} opacity={beamFlicker} />
          )}
          {/* Saucer rim */}
          <ellipse rx={22} ry={8} fill="#6d28d9" />
          <ellipse rx={22} ry={8} fill={`url(#${g}-sc)`} />
          {/* Dome */}
          <ellipse cy={-7.5} rx={11} ry={9} fill="#7c3aed" />
          <ellipse cy={-7.5} rx={11} ry={9} fill={`url(#${g}-sc)`} />
          {/* Rim lights */}
          {[-12,-6,0,6,12].map((lx, i) => (
            <circle key={i} cx={lx} cy={2} r={2.3}
              fill={i%2===0 ? "#a78bfa" : "#34d399"}
              opacity={.78 + .22*Math.sin(tick*.42 + i*1.2)} />
          ))}
          {/* Dome window */}
          <ellipse cy={-8.5} rx={5.5} ry={4.5} fill="#ede9fe" opacity={0.65} />
          <ellipse cx={-1}   cy={-10} rx={2}   ry={1.5}  fill="white" opacity={0.45} />
          {/* Highlight */}
          <ellipse cx={0} cy={-.5} rx={17} ry={3} fill="white" opacity={0.1} />
        </g>
      )}
    </svg>
  );
}

/* ═══════════════════════════════════════════════════════════════
   DNA HELIX  (default theme)
═══════════════════════════════════════════════════════════════ */
const GREEN = "#10b981";
const BLUE  = "#3b82f6";
const GRAY  = "#cbd5e1";

export default function DNAProgress({ percent, logs, currentStep, theme = "default" }: Props) {
  const [phase, setPhase] = useState(0);
  const animRef = useRef<ReturnType<typeof setInterval> | null>(null);

  useEffect(() => {
    animRef.current = setInterval(() => setPhase(p => (p + 3) % 360), 40);
    return () => { if (animRef.current) clearInterval(animRef.current); };
  }, []);

  const isSpace = theme === "orange" || theme === "purple";

  const pctGrad = theme === "orange"
    ? "linear-gradient(90deg,#ea580c,#f97316)"
    : theme === "purple"
    ? "linear-gradient(90deg,#7c3aed,#a78bfa)"
    : `linear-gradient(90deg,${GREEN},${BLUE})`;

  // DNA geometry
  const W = 380, H = 82, cy = H/2, A = 28, wl = 58;
  const rad = (phase * Math.PI) / 180;

  const strandPath = (offset: number) => {
    const pts: string[] = [];
    for (let i = 0; i <= 240; i++) {
      const x = (W * i) / 240;
      const y = cy + A * Math.sin((2*Math.PI*x)/wl + offset);
      pts.push(`${i===0?"M":"L"} ${x.toFixed(1)} ${y.toFixed(1)}`);
    }
    return pts.join(" ");
  };

  const rungs: {x:number;y1:number;y2:number}[] = [];
  for (let x = wl/4; x < W; x += wl/2) {
    rungs.push({
      x,
      y1: cy + A*Math.sin((2*Math.PI*x)/wl + rad),
      y2: cy + A*Math.sin((2*Math.PI*x)/wl + rad + Math.PI),
    });
  }
  const fillX = W * Math.min(percent, 100) / 100;

  return (
    <div className="dna-wrap">

      {/* ── Helix / Space animation ── */}
      <div className="dna-helix-row">
        {isSpace ? (
          <SpaceProgress percent={percent} isOrange={theme === "orange"} />
        ) : (
          <svg width="100%" viewBox={`0 0 ${W} ${H}`}
            style={{ display:"block", filter:"drop-shadow(0 0 8px #10b98133) drop-shadow(0 0 8px #3b82f633)" }}>
            <defs>
              <clipPath id="dna-h-clip">
                <rect x={0} y={0} width={fillX} height={H} />
              </clipPath>
            </defs>
            {/* Gray unfilled layer */}
            <g opacity={0.2}>
              <path d={strandPath(rad)}           stroke={GRAY} strokeWidth={3.5} fill="none" strokeLinecap="round"/>
              <path d={strandPath(rad+Math.PI)}   stroke={GRAY} strokeWidth={3.5} fill="none" strokeLinecap="round"/>
              {rungs.map((r,i)=>(
                <g key={i}>
                  <line x1={r.x} y1={r.y1} x2={r.x} y2={r.y2} stroke={GRAY} strokeWidth={2.5}/>
                  <circle cx={r.x} cy={r.y1} r={3.5} fill={GRAY}/>
                  <circle cx={r.x} cy={r.y2} r={3.5} fill={GRAY}/>
                </g>
              ))}
            </g>
            {/* Colored filled layer */}
            <g clipPath="url(#dna-h-clip)">
              <path d={strandPath(rad)}         stroke={GREEN} strokeWidth={3.5} fill="none" strokeLinecap="round"/>
              <path d={strandPath(rad+Math.PI)} stroke={BLUE}  strokeWidth={3.5} fill="none" strokeLinecap="round"/>
              {rungs.map((r,i)=>{
                const cTop = i%2===0 ? GREEN : BLUE;
                const cBot = i%2===0 ? BLUE  : GREEN;
                const mid  = (r.y1+r.y2)/2;
                return (
                  <g key={i}>
                    <line x1={r.x} y1={r.y1} x2={r.x} y2={mid}  stroke={cTop} strokeWidth={2.5}/>
                    <line x1={r.x} y1={mid}  x2={r.x} y2={r.y2} stroke={cBot} strokeWidth={2.5}/>
                    <circle cx={r.x} cy={r.y1} r={4} fill={cTop}/>
                    <circle cx={r.x} cy={r.y2} r={4} fill={cBot}/>
                  </g>
                );
              })}
            </g>
            {percent>0 && percent<100 && (
              <line x1={fillX} y1={4} x2={fillX} y2={H-4}
                stroke="white" strokeWidth={1.5} strokeDasharray="3,3" opacity={0.6}/>
            )}
          </svg>
        )}
      </div>

      {/* ── % + step label ── */}
      <div className="dna-info-row">
        <div className="dna-pct"
          style={{ background:pctGrad, WebkitBackgroundClip:"text", WebkitTextFillColor:"transparent" }}>
          {Math.round(percent)}%
        </div>
        {currentStep && (
          <div className="dna-step-label"
            style={{ background:pctGrad, WebkitBackgroundClip:"text", WebkitTextFillColor:"transparent" }}>
            {currentStep}
          </div>
        )}
      </div>

      {/* ── Terminal log ── */}
      <div className="dna-log-box">
        <div className="dna-log-header">
          <span className="dna-log-dot red"/>
          <span className="dna-log-dot yellow"/>
          <span className="dna-log-dot green"/>
          <span className="dna-log-title">pipeline log</span>
        </div>
        <div className="dna-log-body">
          {logs.length===0 ? (
            <div className="dna-log-line dim">Waiting for output...</div>
          ) : (
            logs.slice(-3).map((line,i,arr) => (
              <div key={i} className={`dna-log-line ${i===arr.length-1?"latest":"dim"}`}>
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
