import { useState } from "react";

// ── Full parameter set for all pipeline types ─────────────────────────────
export interface PipelineParams {
  // 16S / 12S / 18S-nema
  truncLen_F:      number;
  truncLen_R:      number;
  maxEE_F:         number;
  maxEE_R:         number;
  trimLeft_F:      number;
  trimLeft_R:      number;
  nbases:          number;
  pool:            "FALSE" | "TRUE" | "pseudo";
  chimeraMethod:   "consensus" | "per-sample" | "pooled";
  taxDatabase:     string;
  dbPath:          string;
  minBoot:         number;
  topN:            30 | 50 | 100;
  // ITS
  its_region:      "ITS1" | "ITS2";
  primer_f:        string;
  primer_r:        string;
  // COX1
  truncLen_cox1_f: number;
  truncLen_cox1_r: number;
  codon_table:     5 | 2 | 1;     // 5=invertebrate mt, 2=vertebrate mt, 1=standard
  cox1_min_len:    number;
  cox1_max_len:    number;
  run_lulu:        boolean;
  // PacBio
  pb_min_len:      number;
  pb_max_len:      number;
  pb_maxEE:        number;
  pb_region:       string;
  // Functional
  run_tax4fun:     boolean;
  run_picrust2:    boolean;
}

export const defaultParams: PipelineParams = {
  truncLen_F: 240, truncLen_R: 200,
  maxEE_F: 2, maxEE_R: 2,
  trimLeft_F: 0, trimLeft_R: 0,
  nbases: 1e8,
  pool: "FALSE",
  chimeraMethod: "consensus",
  taxDatabase: "SILVA_138.2_toGenus",
  dbPath: "/home/boss/r16s-app/backend/databases/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  minBoot: 50,
  topN: 30,
  // ITS defaults
  its_region: "ITS1",
  primer_f: "",
  primer_r: "",
  // COX1 defaults
  truncLen_cox1_f: 230,
  truncLen_cox1_r: 200,
  codon_table: 5,
  cox1_min_len: 300,
  cox1_max_len: 330,
  run_lulu: true,
  // PacBio defaults
  pb_min_len: 1000,
  pb_max_len: 1600,
  pb_maxEE: 3.0,
  pb_region: "V1-V9",
  // Functional
  run_tax4fun: false,
  run_picrust2: false,
};

// ── All marker types ──────────────────────────────────────────────────────
export type MarkerType = "16S" | "12S" | "ITS1" | "ITS2" | "COX1" | "18S-nema" | "PacBio";

const MARKER_OPTIONS: { value: MarkerType; label: string; description: string; icon: string }[] = [
  { value: "16S",     icon: "🦠", label: "16S rRNA",         description: "Bacteria & Archaea — SILVA 138 database" },
  { value: "12S",     icon: "🐟", label: "12S rRNA",         description: "Vertebrates / Fish — PR2 or SILVA" },
  { value: "ITS1",    icon: "🍄", label: "ITS1 Fungi",       description: "Fungal ITS1 — UNITE v10 database" },
  { value: "ITS2",    icon: "🍄", label: "ITS2 Fungi",       description: "Fungal ITS2 — UNITE v10 database" },
  { value: "COX1",    icon: "🦑", label: "COX1 / CO1",       description: "Animal metabarcoding — MIDORI2 database" },
  { value: "18S-nema",icon: "🪱", label: "18S Nematode",     description: "Nematode 18S — 18S-NemaBase / PR2" },
  { value: "PacBio",  icon: "🧬", label: "PacBio CCS 16S",   description: "Full-length 16S V1–V9 long reads" },
];

const DB_OPTIONS = [
  {
    value: "SILVA_138.2_toGenus",
    label: "SILVA 138.2 — Genus level (recommended)",
    path: "/home/boss/r16s-app/backend/databases/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    marker: ["16S", "18S-nema", "PacBio"],
  },
  {
    value: "SILVA_138.2_toSpecies",
    label: "SILVA 138.2 — Species level (slower)",
    path: "/home/boss/r16s-app/backend/databases/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
    marker: ["16S", "18S-nema", "PacBio"],
  },
  {
    value: "PR2",
    label: "PR2 v5 SSU (12S / Eukaryotes)",
    path: "/home/boss/r16s-app/backend/databases/pr2_version_5.0.0_SSU_dada2.fasta.gz",
    marker: ["12S", "18S-nema"],
  },
  { value: "custom", label: "Custom path...", path: "", marker: ["16S","12S","18S-nema","PacBio"] },
];

// Default primers per marker
const DEFAULT_PRIMERS: Record<string, { f: string; r: string; label: string }[]> = {
  ITS1: [
    { f: "CTTGGTCATTTAGAGGAAGTAA", r: "GCTGCGTTCTTCATCGATGC", label: "ITS1F / ITS2 (standard)" },
    { f: "GTGYCAGCMGCCGCGGTAA",    r: "GGACTACNVGGGTWTCTAAT", label: "ITS1F / ITS4 (long)" },
  ],
  ITS2: [
    { f: "GTGAATCATCGAATCTTTGAA", r: "TCCTCCGCTTATTGATATGC", label: "ITS3 / ITS4" },
    { f: "GTGYCAGCMGCCGCGGTAA",   r: "GGACTACNVGGGTWTCTAAT", label: "515F / 806R (fallback)" },
  ],
  COX1: [
    { f: "GGWACWGGWTGAACWGTWTAYCCYCC", r: "TANACYTCNGGRTGNCCRAARAAYCA", label: "mlCOIintF / jgHCO2198 (invertebrate, 313 bp)" },
    { f: "CGAATTTGAGGAGGAGGAGGA",      r: "GCTCGTGTGTCTACGAAATCGG",     label: "BF2 / BR2 (broader)" },
  ],
  "18S-nema": [
    { f: "CGCGAATRGCTCATTACAACAGC", r: "GGGCGGTGTGTACAAAGGGCAGGG", label: "NF1 / 18Sr2b (nematode-specific)" },
    { f: "TTGTACACACCGCCC",          r: "CCTTCYGCAGGTTCACCTAC",     label: "TAReuk454FWD1 / TAReukREV3 (eukaryote)" },
  ],
  PacBio: [
    { f: "AGRGTTYGATYMTGGCTCAG", r: "RGYTACCTTGTTACGACTT", label: "27F / 1492R (V1–V9 full-length)" },
    { f: "AGAGTTTGATCMTGGCTCAG", r: "GYTACCTTGTTACGACTT",  label: "8F / 1492R (alternative)" },
  ],
};

interface Props {
  params:    PipelineParams;
  onChange:  (p: PipelineParams) => void;
  marker:    MarkerType;
  onMarker?: (m: MarkerType) => void;
}

export default function PipelineSettings({ params, onChange, marker, onMarker }: Props) {
  const [openStep, setOpenStep] = useState<number | null>(1);
  const set = (key: keyof PipelineParams, val: any) =>
    onChange({ ...params, [key]: val });

  const isITS     = marker === "ITS1" || marker === "ITS2";
  const isCOX1    = marker === "COX1";
  const isPacBio  = marker === "PacBio";
  const isNema    = marker === "18S-nema";
  const is16S     = marker === "16S";
  const isStandard = !isITS && !isCOX1 && !isPacBio;  // 16S/12S/18S-nema use dada2_pipeline.R

  const availableDBs = DB_OPTIONS.filter(db =>
    (db.marker as string[]).includes(marker)
  );

  const defaultPrimers = DEFAULT_PRIMERS[marker] || [];

  return (
    <div className="settings-wrap">

      {/* ── Marker selector ──────────────────────────────────────────────── */}
      {onMarker && (
        <div className="marker-selector-wrap">
          <p className="settings-hint">Select the amplicon marker / pipeline type:</p>
          <div className="marker-grid">
            {MARKER_OPTIONS.map(m => (
              <button
                key={m.value}
                className={`marker-btn ${marker === m.value ? "active" : ""}`}
                onClick={() => onMarker(m.value)}
                title={m.description}
              >
                <span className="marker-icon">{m.icon}</span>
                <span className="marker-label">{m.label}</span>
              </button>
            ))}
          </div>
          {/* Description */}
          <p className="marker-desc">
            {MARKER_OPTIONS.find(m => m.value === marker)?.description}
          </p>
        </div>
      )}

      {/* ── ITS settings ─────────────────────────────────────────────────── */}
      {isITS && (
        <div className="ext-section">
          <h4 className="ext-section-title">ITS Fungal Settings</h4>

          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">ITS Region</label>
              <span className="param-hint">Choose the amplified region</span>
              <div className="topn-toggle">
                {(["ITS1","ITS2"] as const).map(r => (
                  <button key={r}
                    className={`topn-btn ${params.its_region === r ? "active" : ""}`}
                    onClick={() => set("its_region", r)}>
                    {r}
                  </button>
                ))}
              </div>
            </div>
            <ParamNumber label="Max EE Forward" hint="Maximum expected errors (forward reads)"
              value={params.maxEE_F} min={1} max={10} step={0.5}
              onChange={v => set("maxEE_F", v)} />
            <ParamNumber label="Max EE Reverse" hint="Maximum expected errors (reverse reads)"
              value={params.maxEE_R} min={1} max={10} step={0.5}
              onChange={v => set("maxEE_R", v)} />
          </div>

          <div className="primer-section">
            <label className="param-label">Primers</label>
            <span className="param-hint">
              Leave blank to use defaults, or select / enter custom primers
            </span>
            <div className="primer-presets">
              {defaultPrimers.map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => { set("primer_f", pr.f); set("primer_r", pr.r); }}>
                  {pr.label}
                </button>
              ))}
              <button className="primer-preset-btn"
                onClick={() => { set("primer_f", ""); set("primer_r", ""); }}>
                Use defaults
              </button>
            </div>
            {params.primer_f && (
              <div className="primer-inputs">
                <div className="param-item">
                  <label className="param-label">Forward primer (5'→3')</label>
                  <input type="text" className="param-input"
                    value={params.primer_f} onChange={e => set("primer_f", e.target.value)} />
                </div>
                <div className="param-item">
                  <label className="param-label">Reverse primer (5'→3')</label>
                  <input type="text" className="param-input"
                    value={params.primer_r} onChange={e => set("primer_r", e.target.value)} />
                </div>
              </div>
            )}
          </div>

          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">Top Taxa to Show</label>
              <div className="topn-toggle">
                {([30,50,100] as const).map(n => (
                  <button key={n} className={`topn-btn ${params.topN === n ? "active" : ""}`}
                    onClick={() => set("topN", n)}>Top {n}</button>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* ── COX1 settings ────────────────────────────────────────────────── */}
      {isCOX1 && (
        <div className="ext-section">
          <h4 className="ext-section-title">COX1 Metabarcoding Settings</h4>
          <div className="param-grid">
            <ParamNumber label="Truncate Forward (bp)" hint="Truncate forward reads"
              value={params.truncLen_cox1_f} min={0} max={350}
              onChange={v => set("truncLen_cox1_f", v)} />
            <ParamNumber label="Truncate Reverse (bp)" hint="Truncate reverse reads (0 = no truncation)"
              value={params.truncLen_cox1_r} min={0} max={350}
              onChange={v => set("truncLen_cox1_r", v)} />
            <ParamNumber label="Max EE Forward" hint="Maximum expected errors (forward)"
              value={params.maxEE_F} min={1} max={10} step={0.5}
              onChange={v => set("maxEE_F", v)} />
            <ParamNumber label="Max EE Reverse" hint="Maximum expected errors (reverse; COX1 R-reads are often lower quality)"
              value={params.maxEE_R} min={1} max={10} step={0.5}
              onChange={v => set("maxEE_R", v)} />
            <ParamNumber label="Min ASV Length (bp)" hint="Remove ASVs shorter than this (typical mlCOIintF/jgHCO2198 = 300 bp)"
              value={params.cox1_min_len} min={100} max={500}
              onChange={v => set("cox1_min_len", v)} />
            <ParamNumber label="Max ASV Length (bp)" hint="Remove ASVs longer than this"
              value={params.cox1_max_len} min={100} max={500}
              onChange={v => set("cox1_max_len", v)} />
          </div>

          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">Genetic Code (codon table)</label>
              <span className="param-hint">
                Used to filter NUMTs (pseudogenes) via stop-codon translation check
              </span>
              <select className="param-input" value={params.codon_table}
                onChange={e => set("codon_table", Number(e.target.value))}>
                <option value={5}>5 — Invertebrate Mitochondrial (default)</option>
                <option value={2}>2 — Vertebrate Mitochondrial</option>
                <option value={1}>1 — Standard (nuclear)</option>
              </select>
            </div>

            <div className="param-item">
              <label className="param-label">LULU Post-curation</label>
              <span className="param-hint">
                BLAST-based removal of artefact / daughter ASVs (requires BLAST+)
              </span>
              <div className="topn-toggle">
                <button className={`topn-btn ${params.run_lulu ? "active" : ""}`}
                  onClick={() => set("run_lulu", true)}>Enabled</button>
                <button className={`topn-btn ${!params.run_lulu ? "active" : ""}`}
                  onClick={() => set("run_lulu", false)}>Disabled</button>
              </div>
            </div>
          </div>

          {/* Primer presets */}
          <div className="primer-section">
            <label className="param-label">Primer Pair</label>
            <span className="param-hint">Select a preset or enter custom primers</span>
            <div className="primer-presets">
              {defaultPrimers.map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => { set("primer_f", pr.f); set("primer_r", pr.r); }}>
                  {pr.label}
                </button>
              ))}
            </div>
            {params.primer_f && (
              <div className="primer-inputs">
                <div className="param-item">
                  <label className="param-label">Forward primer</label>
                  <input type="text" className="param-input"
                    value={params.primer_f} onChange={e => set("primer_f", e.target.value)} />
                </div>
                <div className="param-item">
                  <label className="param-label">Reverse primer</label>
                  <input type="text" className="param-input"
                    value={params.primer_r} onChange={e => set("primer_r", e.target.value)} />
                </div>
              </div>
            )}
          </div>

          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">Top Taxa to Show</label>
              <div className="topn-toggle">
                {([30,50,100] as const).map(n => (
                  <button key={n} className={`topn-btn ${params.topN === n ? "active" : ""}`}
                    onClick={() => set("topN", n)}>Top {n}</button>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* ── PacBio settings ───────────────────────────────────────────────── */}
      {isPacBio && (
        <div className="ext-section">
          <h4 className="ext-section-title">PacBio CCS Long-Read 16S Settings</h4>
          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">16S Region Amplified</label>
              <span className="param-hint">For record-keeping and metadata output</span>
              <select className="param-input" value={params.pb_region}
                onChange={e => set("pb_region", e.target.value)}>
                <option value="V1-V9">V1–V9 (full-length, ~1500 bp)</option>
                <option value="V1-V3">V1–V3 (~450 bp)</option>
                <option value="V3-V5">V3–V5 (~400 bp)</option>
                <option value="V4">V4 (~250 bp)</option>
                <option value="custom">Custom</option>
              </select>
            </div>
            <ParamNumber label="Min Read Length (bp)" hint="Minimum length for filtering"
              value={params.pb_min_len} min={100} max={2000}
              onChange={v => set("pb_min_len", v)} />
            <ParamNumber label="Max Read Length (bp)" hint="Maximum length for filtering"
              value={params.pb_max_len} min={100} max={2500}
              onChange={v => set("pb_max_len", v)} />
            <ParamNumber label="Max EE" hint="Maximum expected errors (CCS reads are high accuracy; 2-5 typical)"
              value={params.pb_maxEE} min={1} max={20} step={0.5}
              onChange={v => set("pb_maxEE", v)} />
          </div>

          <div className="primer-section">
            <label className="param-label">Primer Pair</label>
            <span className="param-hint">Select a preset for your full-length 16S primers</span>
            <div className="primer-presets">
              {defaultPrimers.map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => { set("primer_f", pr.f); set("primer_r", pr.r); }}>
                  {pr.label}
                </button>
              ))}
            </div>
          </div>

          {/* Taxonomy database */}
          <div className="param-item" style={{ marginTop: "12px" }}>
            <label className="param-label">Taxonomy Database</label>
            <span className="param-hint">SILVA recommended for full-length 16S</span>
            <div className="db-options">
              {availableDBs.map(db => (
                <div key={db.value}
                  className={`db-option ${params.taxDatabase === db.value ? "active" : ""}`}
                  onClick={() => { set("taxDatabase", db.value); if (db.value !== "custom") set("dbPath", db.path); }}>
                  <span className="db-radio">{params.taxDatabase === db.value ? "🔵" : "⚪"}</span>
                  <span className="db-label">{db.label}</span>
                </div>
              ))}
            </div>
          </div>

          <div className="param-grid">
            <ParamNumber label="Min Bootstrap (%)" hint="Minimum confidence for taxonomy"
              value={params.minBoot} min={0} max={100} step={5}
              onChange={v => set("minBoot", v)} />
            <div className="param-item">
              <label className="param-label">Top Taxa to Show</label>
              <div className="topn-toggle">
                {([30,50,100] as const).map(n => (
                  <button key={n} className={`topn-btn ${params.topN === n ? "active" : ""}`}
                    onClick={() => set("topN", n)}>Top {n}</button>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* ── Standard settings (16S / 12S / 18S-nema) ─────────────────────── */}
      {isStandard && (
        <>
          <p className="settings-hint">
            Configure DADA2 parameters for {marker} analysis.
            Default values are suitable for most datasets.
          </p>

          {[
            { id: 1, icon: "🧹", title: "Preprocessing & Filtering",    color: "#3b82f6" },
            { id: 2, icon: "📉", title: "Learn Error Rates",             color: "#8b5cf6" },
            { id: 3, icon: "🔬", title: "Denoising / ASV Inference",     color: "#ec4899" },
            { id: 4, icon: "📊", title: "Sequence Table & Chimera",      color: "#f59e0b" },
            { id: 5, icon: "🏷️", title: "Taxonomic Assignment",         color: "#10b981" },
          ].map(s => (
            <div key={s.id} className="step-accordion">
              <div className="step-acc-header"
                style={{ borderLeft: `4px solid ${s.color}` }}
                onClick={() => setOpenStep(openStep === s.id ? null : s.id)}>
                <span className="step-acc-icon">{s.icon}</span>
                <span className="step-acc-title">Step {s.id} — {s.title}</span>
                <span className="step-acc-chev">{openStep === s.id ? "▲" : "▼"}</span>
              </div>

              {openStep === s.id && (
                <div className="step-acc-body">

                  {s.id === 1 && (
                    <div className="param-grid">
                      <ParamNumber label="Truncate Forward (bp)" hint="Truncate forward reads at this length"
                        value={params.truncLen_F} min={50} max={300}
                        onChange={v => set("truncLen_F", v)} />
                      <ParamNumber label="Truncate Reverse (bp)" hint="Truncate reverse reads at this length"
                        value={params.truncLen_R} min={50} max={300}
                        onChange={v => set("truncLen_R", v)} />
                      <ParamNumber label="Max EE Forward" hint="Maximum expected errors (forward)"
                        value={params.maxEE_F} min={1} max={10} step={0.5}
                        onChange={v => set("maxEE_F", v)} />
                      <ParamNumber label="Max EE Reverse" hint="Maximum expected errors (reverse)"
                        value={params.maxEE_R} min={1} max={10} step={0.5}
                        onChange={v => set("maxEE_R", v)} />
                      <ParamNumber label="Trim Left Forward (bp)" hint="Trim primers from left — forward"
                        value={params.trimLeft_F} min={0} max={50}
                        onChange={v => set("trimLeft_F", v)} />
                      <ParamNumber label="Trim Left Reverse (bp)" hint="Trim primers from left — reverse"
                        value={params.trimLeft_R} min={0} max={50}
                        onChange={v => set("trimLeft_R", v)} />
                      {isNema && (
                        <div className="primer-section" style={{ gridColumn: "1/-1" }}>
                          <label className="param-label">Nematode 18S Primers</label>
                          <span className="param-hint">Select standard primers or enter custom</span>
                          <div className="primer-presets">
                            {(DEFAULT_PRIMERS["18S-nema"] || []).map((pr, i) => (
                              <button key={i}
                                className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                                onClick={() => { set("primer_f", pr.f); set("primer_r", pr.r); }}>
                                {pr.label}
                              </button>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>
                  )}

                  {s.id === 2 && (
                    <div className="param-grid">
                      <ParamSelect label="Bases to learn from"
                        hint="More bases = higher accuracy but slower"
                        value={String(params.nbases)}
                        options={[
                          { value: "1e7", label: "10M bases (fast)" },
                          { value: "1e8", label: "100M bases (default)" },
                          { value: "5e8", label: "500M bases (accurate)" },
                        ]}
                        onChange={v => set("nbases", Number(v))} />
                    </div>
                  )}

                  {s.id === 3 && (
                    <div className="param-grid">
                      <ParamSelect label="Pooling Mode"
                        hint="How to pool samples for inference"
                        value={params.pool}
                        options={[
                          { value: "FALSE",  label: "FALSE — fast (independent)" },
                          { value: "pseudo", label: "pseudo — balanced (recommended)" },
                          { value: "TRUE",   label: "TRUE — most sensitive, slowest" },
                        ]}
                        onChange={v => set("pool", v)} />
                    </div>
                  )}

                  {s.id === 4 && (
                    <div className="param-grid">
                      <ParamSelect label="Chimera Detection Method"
                        hint="Method for identifying chimeric sequences"
                        value={params.chimeraMethod}
                        options={[
                          { value: "consensus",  label: "consensus (recommended)" },
                          { value: "per-sample", label: "per-sample" },
                          { value: "pooled",     label: "pooled (slowest)" },
                        ]}
                        onChange={v => set("chimeraMethod", v)} />
                    </div>
                  )}

                  {s.id === 5 && (
                    <div className="param-grid-full">
                      <div className="param-item">
                        <label className="param-label">Taxonomy Database</label>
                        <span className="param-hint">
                          {marker === "16S" ? "For Bacteria / Archaea"
                            : marker === "12S" ? "For Vertebrates / Fish"
                            : "For Nematodes / Eukaryotes"}
                        </span>
                        <div className="db-options">
                          {availableDBs.map(db => (
                            <div key={db.value}
                              className={`db-option ${params.taxDatabase === db.value ? "active" : ""}`}
                              onClick={() => {
                                set("taxDatabase", db.value);
                                if (db.value !== "custom") set("dbPath", db.path);
                              }}>
                              <div className="db-option-top">
                                <span className="db-radio">{params.taxDatabase === db.value ? "🔵" : "⚪"}</span>
                                <span className="db-label">{db.label}</span>
                              </div>
                              {db.path && db.value !== "custom" && (
                                <code className="db-path">{db.path.split("/").pop()}</code>
                              )}
                            </div>
                          ))}
                        </div>
                      </div>

                      {params.taxDatabase === "custom" && (
                        <div className="param-item">
                          <label className="param-label">Custom Database Path</label>
                          <input type="text" className="param-input"
                            placeholder="/home/boss/r16s-app/backend/databases/mydb.fa.gz"
                            value={params.dbPath}
                            onChange={e => set("dbPath", e.target.value)} />
                        </div>
                      )}

                      <div className="param-item">
                        <ParamNumber label="Min Bootstrap (%)"
                          hint="Minimum confidence for taxonomy assignment"
                          value={params.minBoot} min={0} max={100} step={5}
                          onChange={v => set("minBoot", v)} />
                      </div>

                      <div className="param-item">
                        <label className="param-label">Top Taxa to Show</label>
                        <div className="topn-toggle">
                          {([30,50,100] as const).map(n => (
                            <button key={n}
                              className={`topn-btn ${params.topN === n ? "active" : ""}`}
                              onClick={() => set("topN", n)}>Top {n}</button>
                          ))}
                        </div>
                      </div>
                    </div>
                  )}
                </div>
              )}
            </div>
          ))}
        </>
      )}

      {/* ── Functional prediction section (16S / PacBio / 18S-nema) ─────── */}
      {(is16S || isPacBio) && (
        <div className="ext-section" style={{ marginTop: "16px" }}>
          <h4 className="ext-section-title">Functional Prediction (optional)</h4>
          <p className="param-hint" style={{ marginBottom: "10px" }}>
            Predict functional gene content from 16S composition.
            Requires additional setup — see install_extensions.sh
          </p>
          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">Tax4Fun2</label>
              <span className="param-hint">
                KEGG KO prediction from 16S (R package; fast, no extra conda needed)
              </span>
              <div className="topn-toggle">
                <button className={`topn-btn ${params.run_tax4fun ? "active" : ""}`}
                  onClick={() => set("run_tax4fun", true)}>Enable</button>
                <button className={`topn-btn ${!params.run_tax4fun ? "active" : ""}`}
                  onClick={() => set("run_tax4fun", false)}>Skip</button>
              </div>
            </div>

            <div className="param-item">
              <label className="param-label">PICRUSt2</label>
              <span className="param-hint">
                Phylogenetic placement + KEGG/MetaCyc pathway prediction (conda env required)
              </span>
              <div className="topn-toggle">
                <button className={`topn-btn ${params.run_picrust2 ? "active" : ""}`}
                  onClick={() => set("run_picrust2", true)}>Enable</button>
                <button className={`topn-btn ${!params.run_picrust2 ? "active" : ""}`}
                  onClick={() => set("run_picrust2", false)}>Skip</button>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

// ── Sub-components ────────────────────────────────────────────────────────

function ParamNumber({ label, hint, value, min, max, step = 1, onChange }: {
  label: string; hint: string; value: number;
  min: number; max: number; step?: number;
  onChange: (v: number) => void;
}) {
  return (
    <div className="param-item">
      <label className="param-label">{label}</label>
      <span className="param-hint">{hint}</span>
      <input type="number" className="param-input"
        value={value} min={min} max={max} step={step}
        onChange={e => onChange(Number(e.target.value))} />
    </div>
  );
}

function ParamSelect({ label, hint, value, options, onChange }: {
  label: string; hint: string; value: string;
  options: { value: string; label: string }[];
  onChange: (v: string) => void;
}) {
  return (
    <div className="param-item">
      <label className="param-label">{label}</label>
      <span className="param-hint">{hint}</span>
      <select className="param-input" value={value}
        onChange={e => onChange(e.target.value)}>
        {options.map(o => (
          <option key={o.value} value={o.value}>{o.label}</option>
        ))}
      </select>
    </div>
  );
}
