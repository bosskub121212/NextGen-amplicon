import { useState, useEffect, useCallback } from "react";

export interface TaxaColor {
  name: string;
  color: string;
  abundance: number;
}

// Colors per level, keyed by taxon name
type LevelColors = Record<string, string>;
type AllColors   = Record<string, LevelColors>; // { Phylum: {...}, Class: {...}, ... }

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

interface Props { jobId: string; apiBase: string; }

export default function TaxonomyColorPicker({ jobId, apiBase }: Props) {
  const [taxa,    setTaxa]    = useState<TaxaColor[]>([]);
  const [loading, setLoading] = useState(true);
  const [level,   setLevel]   = useState<"Phylum"|"Class"|"Order"|"Genus"|"Species">("Phylum");

  // Store colours for ALL levels simultaneously, so we can send all at once
  const [allColors, setAllColors] = useState<AllColors>({});

  const [saving,   setSaving]  = useState(false);
  const [saveMsg,  setSaveMsg] = useState<{type:"ok"|"err"; text:string} | null>(null);

  // Load taxa from API and merge with already-saved colours for this level
  const loadTaxa = useCallback(async (lvl: string) => {
    setLoading(true);
    try {
      const res = await fetch(`${apiBase}/taxonomy/${jobId}?level=${lvl}`);
      const data = res.ok ? await res.json() : null;
      const rawTaxa: Array<{name:string; abundance:number}> =
        data?.taxa ?? getMock(lvl);

      // Restore previously saved colours for this level, else random
      const savedForLevel = allColors[lvl] ?? {};
      setTaxa(rawTaxa.map((t: any) => ({
        name:      t.name,
        abundance: t.abundance,
        color:     savedForLevel[t.name] ?? randomColor(t.name),
      })));
    } catch {
      setTaxa(getMock(lvl).map((t) => ({
        ...t,
        color: (allColors[lvl] ?? {})[t.name] ?? randomColor(t.name),
      })));
    }
    setLoading(false);
  }, [jobId, apiBase, allColors]);

  // When level changes: first persist current level colours, then load new level
  const switchLevel = (newLvl: "Phylum"|"Class"|"Order"|"Genus"|"Species") => {
    // Save current taxa colours into allColors before switching
    setAllColors((prev) => ({
      ...prev,
      [level]: Object.fromEntries(taxa.map((t) => [t.name, t.color])),
    }));
    setLevel(newLvl);
  };

  useEffect(() => { loadTaxa(level); }, [level]);   // re-load when level changes

  const updateColor = (name: string, color: string) =>
    setTaxa((p) => p.map((t) => t.name === name ? { ...t, color } : t));

  // "Save Colors for Report" — send ALL levels' colours to backend → replot
  const saveColors = async () => {
    setSaving(true);
    setSaveMsg(null);

    // Merge current view into allColors before sending
    const merged: AllColors = {
      ...allColors,
      [level]: Object.fromEntries(taxa.map((t) => [t.name, t.color])),
    };

    try {
      const res = await fetch(`${apiBase}/replot/${jobId}`, {
        method:  "POST",
        headers: { "Content-Type": "application/json" },
        body:    JSON.stringify({ colors: merged }),
      });
      if (res.ok) {
        setAllColors(merged);   // persist to state
        setSaveMsg({ type:"ok", text:"✅ Plots regenerated! Re-download results to get the new PDFs." });
      } else {
        const err = await res.json().catch(() => ({}));
        setSaveMsg({ type:"err", text: `❌ ${err.error ?? "Replot failed"}` });
      }
    } catch (e: any) {
      setSaveMsg({ type:"err", text:`❌ Network error: ${e.message}` });
    }
    setSaving(false);
  };

  const total = taxa.reduce((s, t) => s + t.abundance, 0) || 1;

  if (loading) return <div className="tc-loading">Loading taxonomy data…</div>;

  return (
    <div className="tc-wrap">
      <div className="tc-header">
        <h3>🎨 Taxonomy Color Settings</h3>
        <div className="tc-level-group">
          {(["Phylum","Class","Order","Genus","Species"] as const).map((l) => (
            <button key={l}
              className={`tc-level-btn ${level === l ? "active" : ""}`}
              onClick={() => switchLevel(l)}>{l}</button>
          ))}
        </div>
      </div>

      {/* Stacked bar preview */}
      <div className="tc-bar-label">Stacked Bar Chart Preview</div>
      <div className="tc-stacked-bar">
        {taxa.map((t) => (
          <div key={t.name} className="tc-bar-segment"
            style={{ width:`${(t.abundance/total)*100}%`, background:t.color }}
            title={`${t.name}: ${t.abundance.toFixed(1)}%`} />
        ))}
      </div>

      {/* Color editor */}
      <div className="tc-list">
        {taxa.map((t) => (
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
                {PALETTE.map((c) => (
                  <div key={c}
                    className={`tc-dot ${t.color === c ? "selected" : ""}`}
                    style={{ background:c }}
                    onClick={() => updateColor(t.name, c)} />
                ))}
              </div>
              <input type="color" className="tc-color-input"
                value={t.color} title="Custom color"
                onChange={(e) => updateColor(t.name, e.target.value)} />
            </div>
          </div>
        ))}
      </div>

      {/* Legend */}
      <div className="tc-legend">
        <div className="tc-legend-title">Legend Preview</div>
        <div className="tc-legend-items">
          {taxa.map((t) => (
            <div key={t.name} className="tc-legend-item">
              <div className="tc-legend-dot" style={{ background:t.color }} />
              <span>{t.name} ({t.abundance.toFixed(1)}%)</span>
            </div>
          ))}
        </div>
      </div>

      {/* Save button */}
      <button className="tc-apply-btn"
        onClick={saveColors}
        disabled={saving}
        style={{ opacity: saving ? 0.7 : 1, cursor: saving ? "wait" : "pointer" }}>
        {saving ? "⏳ Regenerating plots…" : "💾 Save Colors for Report"}
      </button>

      {saveMsg && (
        <div className={`tc-save-msg ${saveMsg.type === "ok" ? "tc-save-ok" : "tc-save-err"}`}>
          {saveMsg.text}
        </div>
      )}

      <div className="tc-hint">
        Colours apply to all taxonomy stacked bar PDFs.<br/>
        Switch levels above to customise Phylum, Class, Order, Genus, and Species separately.
      </div>
    </div>
  );
}

// ── Mock data for offline preview ─────────────────────────────
function getMock(level: string): Array<{name:string; abundance:number}> {
  const mocks: Record<string, Array<{name:string; abundance:number}>> = {
    Phylum: [
      {name:"Bacillota",       abundance:38.5},
      {name:"Bacteroidota",    abundance:28.1},
      {name:"Actinomycetota",  abundance:18.4},
      {name:"Pseudomonadota",  abundance:10.5},
      {name:"Other",           abundance:4.5},
    ],
    Class: [
      {name:"Clostridia",      abundance:32.1},
      {name:"Bacteroidia",     abundance:25.3},
      {name:"Actinomycetia",   abundance:15.8},
      {name:"Bacilli",         abundance:12.4},
      {name:"Other",           abundance:14.4},
    ],
    Order: [
      {name:"Lachnospirales",  abundance:22.0},
      {name:"Bacteroidales",   abundance:25.3},
      {name:"Oscillospirales", abundance:14.0},
      {name:"Actinomycetales", abundance:12.2},
      {name:"Other",           abundance:26.5},
    ],
    Genus: [
      {name:"Blautia",         abundance:14.2},
      {name:"Bacteroides",     abundance:18.7},
      {name:"Faecalibacterium",abundance:10.1},
      {name:"Lachnospiraceae", abundance:9.3},
      {name:"Ruminococcus",    abundance:8.5},
      {name:"Other",           abundance:39.2},
    ],
    Species: [
      {name:"Blautia obeum",          abundance:9.1},
      {name:"Bacteroides uniformis",  abundance:11.3},
      {name:"Faecalibacterium prausnitzii", abundance:8.8},
      {name:"Ruminococcus gnavus",    abundance:6.4},
      {name:"Other",                  abundance:64.4},
    ],
  };
  return mocks[level] ?? mocks["Phylum"];
}
