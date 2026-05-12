import { useState, useEffect, useCallback } from "react";

export interface TaxaColor {
  name: string;
  color: string;
  abundance: number;
}

type LevelColors  = Record<string, string>;
type AllColors    = Record<string, LevelColors>;

// ── Per-taxon swatch palette ─────────────────────────────────
const PALETTE = [
  "#ef4444","#f97316","#f59e0b","#eab308","#84cc16",
  "#22c55e","#10b981","#06b6d4","#3b82f6","#6366f1",
  "#8b5cf6","#ec4899","#14b8a6","#0ea5e9","#a855f7",
  "#64748b","#d1d5db","#1d4ed8","#b45309","#065f46",
];

function randomColor(seed: string): string {
  let hash = 0;
  for (let i = 0; i < seed.length; i++) hash = seed.charCodeAt(i) + ((hash << 5) - hash);
  return PALETTE[Math.abs(hash) % PALETTE.length];
}

// ── Heatmap palettes ─────────────────────────────────────────
export interface HeatmapPaletteOption {
  id:    string;
  label: string;
  stops: [string, string, string];
}

export const HEATMAP_PALETTE_OPTIONS: HeatmapPaletteOption[] = [
  { id: "blue",      label: "Blue",      stops: ["#f0f4ff","#3b82f6","#1e1b4b"] },
  { id: "orange",    label: "Orange",    stops: ["#fff7ed","#f97316","#431407"] },
  { id: "green",     label: "Green",     stops: ["#f0fff4","#22c55e","#14532d"] },
  { id: "purple",    label: "Purple",    stops: ["#faf5ff","#a855f7","#3b0764"] },
  { id: "red",       label: "Red",       stops: ["#fff1f2","#ef4444","#7f1d1d"] },
  { id: "teal",      label: "Teal",      stops: ["#f0fdfa","#14b8a6","#042f2e"] },
  { id: "pink",      label: "Pink",      stops: ["#fdf2f8","#ec4899","#500724"] },
  { id: "viridis",   label: "Viridis",   stops: ["#440154","#21918c","#fde725"] },
  { id: "plasma",    label: "Plasma",    stops: ["#0d0887","#cc4778","#f0f921"] },
  { id: "rdbu",      label: "Red↔Blue",  stops: ["#1e40af","#ffffff","#b91c1c"] },
  { id: "greyscale", label: "Greyscale", stops: ["#f8fafc","#64748b","#0f172a"] },
  { id: "brown",     label: "Brown",     stops: ["#fdf8f0","#d97706","#451a03"] },
];

interface HeatmapTypeConfig { key: string; label: string; default: string; }
const HEATMAP_TYPES: HeatmapTypeConfig[] = [
  { key: "taxonomy_genus",  label: "Taxonomy — Genus",   default: "blue"   },
  { key: "taxonomy_family", label: "Taxonomy — Family",  default: "orange" },
  { key: "taxonomy_phylum", label: "Taxonomy — Phylum",  default: "green"  },
  { key: "beta_bray",       label: "Beta — Bray-Curtis", default: "blue"   },
  { key: "beta_jaccard",    label: "Beta — Jaccard",     default: "blue"   },
];

type TaxaLevel = "Phylum"|"Class"|"Order"|"Genus"|"Species";
type ActiveTab = TaxaLevel | "Heatmaps";

interface Props { jobId: string; apiBase: string; }

// ── Gradient preview swatch ───────────────────────────────────
function GradientSwatch({ stops, selected, onClick, label }:
  { stops: [string,string,string]; selected: boolean; onClick:()=>void; label:string }) {
  const grad = `linear-gradient(to right, ${stops[0]}, ${stops[1]}, ${stops[2]})`;
  return (
    <div className="tc-pal-swatch-wrap" onClick={onClick}>
      <div style={{
        background:   grad,
        height:       22,
        borderRadius: 5,
        border:       selected ? "2.5px solid #6366f1" : "2px solid #334155",
        boxShadow:    selected ? "0 0 0 2px #a5b4fc" : "none",
        cursor:       "pointer",
        transition:   "box-shadow 0.15s, border 0.1s",
      }} />
      <div className="tc-pal-swatch-label" style={{ color: selected ? "#818cf8" : "#94a3b8" }}>
        {label}
      </div>
    </div>
  );
}

export default function TaxonomyColorPicker({ jobId, apiBase }: Props) {
  const [taxa,    setTaxa]    = useState<TaxaColor[]>([]);
  const [loading, setLoading] = useState(true);
  const [level,   setLevel]   = useState<TaxaLevel>("Phylum");
  const [tab,     setTab]     = useState<ActiveTab>("Phylum");

  const [allColors,       setAllColors]       = useState<AllColors>({});
  const [heatmapPalettes, setHeatmapPalettes] = useState<Record<string,string>>(
    () => Object.fromEntries(HEATMAP_TYPES.map(t => [t.key, t.default]))
  );

  const [saving,  setSaving]  = useState(false);
  const [saveMsg, setSaveMsg] = useState<{type:"ok"|"err"; text:string}|null>(null);

  // Load saved colours from taxonomy_summary (custom_colors embedded by replot.R)
  useEffect(() => {
    fetch(`${apiBase}/taxonomy/${jobId}?level=Phylum`)
      .then(r => r.json())
      .then(data => {
        if (data?.custom_colors) {
          const cc = data.custom_colors as AllColors & { heatmaps?: Record<string,string> };
          const { heatmaps, ...taxaColors } = cc as any;
          if (Object.keys(taxaColors).length) setAllColors(taxaColors);
          if (heatmaps) setHeatmapPalettes(prev => ({ ...prev, ...heatmaps }));
        }
      })
      .catch(() => {});
  }, [jobId, apiBase]);

  const loadTaxa = useCallback(async (lvl: string) => {
    setLoading(true);
    try {
      const res  = await fetch(`${apiBase}/taxonomy/${jobId}?level=${lvl}`);
      const data = res.ok ? await res.json() : null;
      const rawTaxa: Array<{name:string; abundance:number}> = data?.taxa ?? getMock(lvl);
      const saved = allColors[lvl] ?? {};
      setTaxa(rawTaxa.map((t: any) => ({
        name: t.name, abundance: t.abundance,
        color: saved[t.name] ?? randomColor(t.name),
      })));
    } catch {
      setTaxa(getMock(lvl).map(t => ({
        ...t, color: (allColors[lvl] ?? {})[t.name] ?? randomColor(t.name),
      })));
    }
    setLoading(false);
  }, [jobId, apiBase, allColors]);

  const switchTab = (newTab: ActiveTab) => {
    // Persist current taxa colours before switching away
    if (tab !== "Heatmaps") {
      setAllColors(prev => ({
        ...prev,
        [level]: Object.fromEntries(taxa.map(t => [t.name, t.color])),
      }));
    }
    setTab(newTab);
    if (newTab !== "Heatmaps") setLevel(newTab as TaxaLevel);
  };

  useEffect(() => {
    if (tab !== "Heatmaps") loadTaxa(level);
  }, [level, tab]);

  const updateColor = (name: string, color: string) =>
    setTaxa(p => p.map(t => t.name === name ? { ...t, color } : t));

  const saveColors = async () => {
    setSaving(true);
    setSaveMsg(null);
    const taxaSnapshot = tab !== "Heatmaps"
      ? { ...allColors, [level]: Object.fromEntries(taxa.map(t => [t.name, t.color])) }
      : allColors;
    const merged: AllColors = { ...taxaSnapshot, heatmaps: heatmapPalettes as any };
    try {
      const res = await fetch(`${apiBase}/replot/${jobId}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ colors: merged }),
      });
      if (res.ok) {
        setAllColors(merged);
        setSaveMsg({ type:"ok", text:"✅ Plots regenerated! Re-download results to get the new PDFs." });
      } else {
        const err = await res.json().catch(() => ({}));
        setSaveMsg({ type:"err", text:`❌ ${err.error ?? "Replot failed"}` });
      }
    } catch (e: any) {
      setSaveMsg({ type:"err", text:`❌ Network error: ${e.message}` });
    }
    setSaving(false);
  };

  const total = taxa.reduce((s, t) => s + t.abundance, 0) || 1;

  return (
    <div className="tc-wrap">
      {/* ── Tab bar ── */}
      <div className="tc-header">
        <h3>🎨 Color Settings</h3>
        <div className="tc-level-group">
          {(["Phylum","Class","Order","Genus","Species"] as TaxaLevel[]).map(l => (
            <button key={l}
              className={`tc-level-btn ${tab === l ? "active" : ""}`}
              onClick={() => switchTab(l)}>{l}
            </button>
          ))}
          <button
            className={`tc-level-btn ${tab === "Heatmaps" ? "active" : ""}`}
            style={tab === "Heatmaps" ? {} : { background:"#1e293b", border:"1px solid #334155" }}
            onClick={() => switchTab("Heatmaps")}>
            🌡 Heatmaps
          </button>
        </div>
      </div>

      {/* ══════════════════════════════════
          HEATMAP PALETTE PICKER
      ══════════════════════════════════ */}
      {tab === "Heatmaps" && (
        <div className="tc-heatmap-section">
          <p className="tc-heatmap-desc">
            เลือกโทนสีสำหรับแต่ละ Heatmap แล้วกด Save — PDF จะถูก regenerate อัตโนมัติ
          </p>

          {HEATMAP_TYPES.map(ht => (
            <div key={ht.key} className="tc-heatmap-row">
              <div className="tc-heatmap-label">{ht.label}</div>
              <div className="tc-heatmap-palettes">
                {HEATMAP_PALETTE_OPTIONS.map(p => (
                  <GradientSwatch
                    key={p.id}
                    stops={p.stops}
                    label={p.label}
                    selected={heatmapPalettes[ht.key] === p.id}
                    onClick={() => setHeatmapPalettes(prev => ({ ...prev, [ht.key]: p.id }))}
                  />
                ))}
              </div>
            </div>
          ))}
        </div>
      )}

      {/* ══════════════════════════════════
          TAXONOMY COLOR EDITOR
      ══════════════════════════════════ */}
      {tab !== "Heatmaps" && (
        <>
          {loading ? (
            <div className="tc-loading">Loading taxonomy data…</div>
          ) : (
            <>
              <div className="tc-bar-label">Stacked Bar Chart Preview</div>
              <div className="tc-stacked-bar">
                {taxa.map(t => (
                  <div key={t.name} className="tc-bar-segment"
                    style={{ width:`${(t.abundance/total)*100}%`, background:t.color }}
                    title={`${t.name}: ${t.abundance.toFixed(1)}%`} />
                ))}
              </div>

              <div className="tc-list">
                {taxa.map(t => (
                  <div key={t.name} className="tc-item">
                    <div className="tc-item-left">
                      <div className="tc-swatch" style={{ background:t.color }} />
                      <div>
                        <div className="tc-name">{t.name}</div>
                        <div className="tc-pct">{t.abundance.toFixed(1)}%</div>
                      </div>
                    </div>
                    <div className="tc-item-right">
                      <div className="tc-palette">
                        {PALETTE.map(c => (
                          <div key={c}
                            className={`tc-dot ${t.color === c ? "selected" : ""}`}
                            style={{ background:c }}
                            onClick={() => updateColor(t.name, c)} />
                        ))}
                      </div>
                      <input type="color" className="tc-color-input"
                        value={t.color} title="Custom color"
                        onChange={e => updateColor(t.name, e.target.value)} />
                    </div>
                  </div>
                ))}
              </div>

              <div className="tc-legend">
                <div className="tc-legend-title">Legend Preview</div>
                <div className="tc-legend-items">
                  {taxa.map(t => (
                    <div key={t.name} className="tc-legend-item">
                      <div className="tc-legend-dot" style={{ background:t.color }} />
                      <span>{t.name} ({t.abundance.toFixed(1)}%)</span>
                    </div>
                  ))}
                </div>
              </div>
            </>
          )}
        </>
      )}

      {/* ── Save button ── */}
      <button className="tc-apply-btn"
        onClick={saveColors} disabled={saving}
        style={{ opacity: saving ? 0.7 : 1, cursor: saving ? "wait" : "pointer" }}>
        {saving ? "⏳ Regenerating plots…" : "💾 Save Colors for Report"}
      </button>

      {saveMsg && (
        <div className={`tc-save-msg ${saveMsg.type==="ok" ? "tc-save-ok" : "tc-save-err"}`}>
          {saveMsg.text}
        </div>
      )}

      <div className="tc-hint">
        {tab === "Heatmaps"
          ? "โทนสีที่เลือกใช้กับ Heatmap PDFs ทั้งหมด — กด Save เพื่อ regenerate"
          : "สีที่เลือกใช้กับ stacked bar charts ทุก PDF — เปลี่ยน level ได้ด้านบน"}
      </div>
    </div>
  );
}

// ── Mock data ─────────────────────────────────────────────────
function getMock(level: string): Array<{name:string; abundance:number}> {
  const mocks: Record<string, Array<{name:string; abundance:number}>> = {
    Phylum: [
      {name:"Bacillota",abundance:38.5},{name:"Bacteroidota",abundance:28.1},
      {name:"Actinomycetota",abundance:18.4},{name:"Pseudomonadota",abundance:10.5},
      {name:"Other",abundance:4.5},
    ],
    Class: [
      {name:"Clostridia",abundance:32.1},{name:"Bacteroidia",abundance:25.3},
      {name:"Actinomycetia",abundance:15.8},{name:"Bacilli",abundance:12.4},
      {name:"Other",abundance:14.4},
    ],
    Order: [
      {name:"Lachnospirales",abundance:22.0},{name:"Bacteroidales",abundance:25.3},
      {name:"Oscillospirales",abundance:14.0},{name:"Actinomycetales",abundance:12.2},
      {name:"Other",abundance:26.5},
    ],
    Genus: [
      {name:"Blautia",abundance:14.2},{name:"Bacteroides",abundance:18.7},
      {name:"Faecalibacterium",abundance:10.1},{name:"Lachnospiraceae",abundance:9.3},
      {name:"Ruminococcus",abundance:8.5},{name:"Other",abundance:39.2},
    ],
    Species: [
      {name:"Blautia obeum",abundance:9.1},{name:"Bacteroides uniformis",abundance:11.3},
      {name:"Faecalibacterium prausnitzii",abundance:8.8},
      {name:"Ruminococcus gnavus",abundance:6.4},{name:"Other",abundance:64.4},
    ],
  };
  return mocks[level] ?? mocks["Phylum"];
}
