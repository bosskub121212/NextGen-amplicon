import { useState, useEffect } from "react";
import DbFileBrowser from "./DbFileBrowser";
import MetadataUpload, { type MetadataInfo } from "./MetadataUpload";

// ── Full parameter set ────────────────────────────────────────────────────────
export interface PipelineParams {
  // DADA2 core
  truncLen_F:        number;
  truncLen_R:        number;
  maxEE_F:           number;
  maxEE_R:           number;
  trimLeft_F:        number;
  trimLeft_R:        number;
  nbases:            number;
  pool:              "FALSE" | "TRUE" | "pseudo";
  chimeraMethod:     "consensus" | "per-sample" | "pooled";
  nThreads:          number;
  // Taxonomy
  taxDatabase:       string;
  dbPath:            string;
  minBoot:           number;
  topN:              30 | 50 | 100;
  taxMethod:         string;
  confidence:        number;
  collapseLevels:    number[];
  // ITS
  its_region:        "ITS1" | "ITS2";
  primer_f:          string;
  primer_r:          string;
  // COX1
  truncLen_cox1_f:   number;
  truncLen_cox1_r:   number;
  codon_table:       5 | 2 | 1;
  cox1_min_len:      number;
  cox1_max_len:      number;
  run_lulu:          boolean;
  // PacBio
  pb_min_len:        number;
  pb_max_len:        number;
  pb_maxEE:          number;
  pb_region:         string;
  // ONT
  ont_region:        string;
  ont_min_abundance: number;
  ont_db_path:       string;
  // Cutadapt
  cutadaptErrorRate: number;
  cutadaptOverlap:   number;
  discardUntrimmed:  boolean;
  // Phylogeny
  runPhylogeny:      boolean;
  phyloThreads:      number;
  // Diversity
  runDiversity:               boolean;
  samplingDepth:              number;
  groupCol:                   string;
  additionalGroupCols:        string;
  diversityPhyloMetrics:      string[];
  diversityNonPhyloMetrics:   string[];
  // Differential Abundance
  runDiffAbund:      boolean;
  diffAbundLevel:    string;
  diffAbundMinFreq:  number;
  diffAbundPval:     number;
  // Functional
  run_tax4fun:       boolean;
  run_picrust2:      boolean;
  picrust2PlaceTool: string;
  picrust2MaxNSTI:   number;
  picrust2Threads:   number;
  picrust2Databases: string[];
  // Metadata
  metadataPath:      string;
  // Custom Classifier
  customClassifierMode: "default" | "train" | "upload";
  customClassifierPath: string;
  trainAmpliconMinLen:  number;
  trainAmpliconMaxLen:  number;
  // Sequencer type (for DADA2 pipelines)
  sequencerType: "illumina" | "ont";
  // ONT long-read length filter (used instead of truncLen_F when sequencerType === "ont")
  ontMinLen: number;
  ontMaxLen: number;
}

export const defaultParams: PipelineParams = {
  truncLen_F: 240, truncLen_R: 200,
  maxEE_F: 2,      maxEE_R: 2,
  trimLeft_F: 0,   trimLeft_R: 0,
  nbases: 1e8,
  pool: "FALSE",
  chimeraMethod: "consensus",
  nThreads: 4,
  taxDatabase: "SILVA_16S",
  dbPath: "",
  minBoot: 50,
  topN: 30,
  taxMethod: "sklearn",
  confidence: 0.7,
  collapseLevels: [5, 6, 7],
  its_region: "ITS1",
  primer_f: "",
  primer_r: "",
  truncLen_cox1_f: 230, truncLen_cox1_r: 200,
  codon_table: 5,
  cox1_min_len: 300, cox1_max_len: 330,
  run_lulu: true,
  pb_min_len: 1000, pb_max_len: 1600, pb_maxEE: 3.0, pb_region: "V1-V9",
  ont_region: "V1-V9", ont_min_abundance: 0.0001, ont_db_path: "",
  cutadaptErrorRate: 0.1, cutadaptOverlap: 3, discardUntrimmed: false,
  runPhylogeny: true, phyloThreads: 4,
  runDiversity: true,
  samplingDepth: 10000, groupCol: "treatment", additionalGroupCols: "",
  diversityPhyloMetrics: ["weighted_unifrac", "unweighted_unifrac", "faith_pd"],
  diversityNonPhyloMetrics: ["shannon", "observed_features", "evenness", "bray_curtis", "jaccard"],
  runDiffAbund: true,
  diffAbundLevel: "L6", diffAbundMinFreq: 10, diffAbundPval: 0.05,
  run_tax4fun: false, run_picrust2: false,
  picrust2PlaceTool: "epa-ng", picrust2MaxNSTI: 2, picrust2Threads: 8,
  picrust2Databases: ["metacyc", "ec", "ko"],
  metadataPath: "",
  customClassifierMode: "default",
  customClassifierPath: "",
  trainAmpliconMinLen: 200,
  trainAmpliconMaxLen: 600,
  sequencerType: "illumina",
  ontMinLen: 300, ontMaxLen: 600,
};

export type MarkerType = "16S" | "12S" | "ITS1" | "ITS2" | "COX1" | "18S-nema" | "PacBio" | "ONT-16S";

// ── Primer presets (Cutadapt step) ────────────────────────────────────────────
interface PrimerPreset {
  label:      string;
  f:          string;
  r:          string;
  markers:    string[];
  region?:    string;   // named region — used for smart hints
  truncF?:    number;   // suggested DADA2 truncLen_F when this preset is chosen
  truncR?:    number;   // suggested DADA2 truncLen_R
  ampMinLen?: number;   // suggested classifier-training amplicon min length
  ampMaxLen?: number;   // suggested classifier-training amplicon max length
}

const CUTADAPT_PRESETS: PrimerPreset[] = [
  { label: "16S V3–V4",  f: "CCTACGGGNGGCWGCAG",          r: "GACTACHVGGGTATCTAATCC",
    markers: ["16S"],    region: "V3-V4", truncF: 240, truncR: 200, ampMinLen: 350, ampMaxLen: 500 },
  { label: "16S V4",     f: "GTGYCAGCMGCCGCGGTAA",         r: "GGACTACNVGGGTWTCTAAT",
    markers: ["16S"],    region: "V4",    truncF: 240, truncR: 160, ampMinLen: 250, ampMaxLen: 300 },
  { label: "16S V7–V8",  f: "AACMGGATTAGATACCCKG",         r: "ACGTCATCCCCACCTTCC",
    markers: ["16S"],    region: "V7-V8", truncF: 250, truncR: 200, ampMinLen: 300, ampMaxLen: 400 },
  { label: "ITS1",       f: "TCCGTAGGTGAACCTGCGG",         r: "GCTGCGTTCTTCATCGATGC",     markers: ["ITS1"] },
  { label: "ITS2",       f: "GTGAATCATCGAATCTTTGAA",       r: "TCCTCCGCTTATTGATATGC",     markers: ["ITS2"] },
  { label: "18S (V4)",   f: "CCAGCASCYGCGGTAATTCC",        r: "ACTTTCGTTCTTGATYRA",       markers: ["18S-nema"] },
  { label: "12S (fish)", f: "GTCGGTAAAACTCGTGCCAGC",       r: "CATAGTGGGGTATCTAATCCCAGTTTG", markers: ["12S"] },
  { label: "COX1 mlCOI", f: "GGWACWGGWTGAACWGTWTAYCCYCC",  r: "TANACYTCNGGRTGNCCRAARAAYCA",  markers: ["COX1"] },
];

const MARKER_OPTIONS: { value: MarkerType; label: string; description: string; icon: string }[] = [
  { value: "16S",      icon: "🦠", label: "16S rRNA",       description: "Bacteria & Archaea — SILVA 138 database" },
  { value: "12S",      icon: "🐟", label: "12S rRNA",       description: "Vertebrates / Fish — PR2 or SILVA" },
  { value: "ITS1",     icon: "🍄", label: "ITS1 Fungi",     description: "Fungal ITS1 — UNITE v10 database" },
  { value: "ITS2",     icon: "🍄", label: "ITS2 Fungi",     description: "Fungal ITS2 — UNITE v10 database" },
  { value: "COX1",     icon: "🦑", label: "COX1 / CO1",     description: "Animal metabarcoding — MIDORI2 database" },
  { value: "18S-nema", icon: "🐛", label: "18S Nematode",   description: "Nematode 18S — NemaBase / PR2" },
  { value: "PacBio",   icon: "🧬", label: "PacBio CCS 16S", description: "Full-length 16S V1–V9 long reads" },
  { value: "ONT-16S",  icon: "🧫", label: "ONT 16S",        description: "Oxford Nanopore V7-V8 / V1-V9 — Emu pipeline" },
];

const DB_OPTIONS = [
  { value: "SILVA_16S",    label: "SILVA 138.1 — Genus level",       pathKey: "SILVA_16S",    marker: ["16S","18S-nema","PacBio"] },
  { value: "SILVA_16S_sp", label: "SILVA 138.1 — Species level",     pathKey: "SILVA_16S_sp", marker: ["16S","18S-nema","PacBio"] },
  { value: "PR2_18S",      label: "PR2 v5 SSU (12S / Eukaryotes)",   pathKey: "PR2_18S",      marker: ["12S","18S-nema"] },
  { value: "NemaBase_18S", label: "18S-NemaBase (nematode-specific)", pathKey: "NemaBase_18S", marker: ["18S-nema"] },
  { value: "UNITE_ITS1",   label: "UNITE v10 (ITS Fungi)",           pathKey: "UNITE_ITS1",   marker: ["ITS1","ITS2"] },
  { value: "MIDORI2_COX1", label: "MIDORI2 COX1 (Animal barcoding)", pathKey: "MIDORI2_COX1", marker: ["COX1"] },
  { value: "EMU_SILVA",    label: "Emu — SILVA database",             pathKey: "emu_silva",    marker: ["ONT-16S"] },
  { value: "custom",       label: "Custom path...",                   pathKey: "",             marker: ["16S","12S","18S-nema","PacBio","ITS1","ITS2","COX1","ONT-16S"] },
];

const PACBIO_PRIMERS: Record<string, { f: string; r: string; label: string }[]> = {
  "V1-V9": [{ f: "AGRGTTYGATYMTGGCTCAG", r: "RGYTACCTTGTTACGACTT", label: "27F / 1492R (V1–V9)" }],
  "V1-V3": [{ f: "AGAGTTTGATCMTGGCTCAG", r: "ATTACCGCGGCTGCTGG",   label: "8F / 534R (V1–V3)" }],
  "V4":    [{ f: "GTGYCAGCMGCCGCGGTAA",  r: "GGACTACNVGGGTWTCTAAT",label: "515F / 806R (V4)" }],
  "custom":[{ f: "AGRGTTYGATYMTGGCTCAG", r: "RGYTACCTTGTTACGACTT", label: "27F / 1492R (default)" }],
};

const ITS_PRIMERS: Record<string, { f: string; r: string; label: string }[]> = {
  ITS1: [
    { f: "CTTGGTCATTTAGAGGAAGTAA", r: "GCTGCGTTCTTCATCGATGC", label: "ITS1F / ITS2 (standard)" },
    { f: "TCCGTAGGTGAACCTGCGG",   r: "GCTGCGTTCTTCATCGATGC", label: "ITS1 / ITS2 (alt)" },
  ],
  ITS2: [
    { f: "GTGAATCATCGAATCTTTGAA", r: "TCCTCCGCTTATTGATATGC", label: "ITS3 / ITS4" },
  ],
};

const COX1_PRIMERS = [
  { f: "GGWACWGGWTGAACWGTWTAYCCYCC", r: "TANACYTCNGGRTGNCCRAARAAYCA", label: "mlCOIintF / jgHCO2198 (313 bp)" },
  { f: "CGAATTTGAGGAGGAGGAGGA",      r: "GCTCGTGTGTCTACGAAATCGG",     label: "BF2 / BR2 (broader)" },
];

const NEMA_PRIMERS = [
  { f: "CGCGAATRGCTCATTACAACAGC", r: "GGGCGGTGTGTACAAAGGGCAGGG", label: "NF1 / 18Sr2b (nematode)" },
];

// ── Quality plot (simulated) ──────────────────────────────────────────────────
function QualityPlot({ truncF, truncR }: { truncF: number; truncR: number }) {
  const W = 560, H = 190, PL = 38, PR = 12, PT = 12, PB = 28;
  const maxPos = 300, maxQ = 42;
  const fwdQ = (p: number) => Math.max(5, 38 - Math.max(0, (p - 160) * 0.09));
  const revQ = (p: number) => Math.max(5, 35 - p * 0.07 - (p > 180 ? (p - 180) * 0.06 : 0));
  const toX = (p: number) => PL + (p / maxPos) * (W - PL - PR);
  const toY = (q: number) => PT + (1 - q / maxQ) * (H - PT - PB);
  const pts = Array.from({ length: 31 }, (_, i) => i * 10);
  const fwdPath = pts.map((p, i) => `${i ? "L" : "M"}${toX(p).toFixed(1)},${toY(fwdQ(p)).toFixed(1)}`).join(" ");
  const revPath = pts.map((p, i) => `${i ? "L" : "M"}${toX(p).toFixed(1)},${toY(revQ(p)).toFixed(1)}`).join(" ");
  return (
    <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", background: "#0f172a", borderRadius: 8, border: "1px solid #1e293b", display: "block" }}>
      {/* Grid lines */}
      {[0,10,20,30,40].map(q => (
        <g key={q}>
          <line x1={PL} x2={W - PR} y1={toY(q)} y2={toY(q)}
            stroke={q === 25 ? "#334155" : "#1e293b"} strokeWidth={q === 25 ? 1 : 0.5}
            strokeDasharray={q === 25 ? "5,4" : "2,3"} />
          <text x={PL - 4} y={toY(q) + 3} fill="#64748b" fontSize={8} textAnchor="end">{q}</text>
        </g>
      ))}
      <text x={PL - 4} y={toY(25) + 3} fill="#94a3b8" fontSize={8} textAnchor="end">Q25</text>
      {/* X axis */}
      {[0,50,100,150,200,250,300].map(p => (
        <text key={p} x={toX(p)} y={H - 4} fill="#475569" fontSize={8} textAnchor="middle">{p}bp</text>
      ))}
      {/* Curves */}
      <path d={fwdPath} fill="none" stroke="#22d3ee" strokeWidth={1.8} />
      <path d={revPath} fill="none" stroke="#818cf8" strokeWidth={1.8} />
      {pts.map(p => (
        <g key={p}>
          <circle cx={toX(p)} cy={toY(fwdQ(p))} r={2.2} fill="#22d3ee" />
          <circle cx={toX(p)} cy={toY(revQ(p))} r={2.2} fill="#818cf8" />
        </g>
      ))}
      {/* TruncLen markers */}
      {truncF > 0 && truncF <= maxPos && (
        <g>
          <line x1={toX(truncF)} x2={toX(truncF)} y1={PT} y2={H - PB}
            stroke="#22d3ee" strokeWidth={1.5} strokeDasharray="5,3" opacity={0.9} />
          <text x={toX(truncF) + 3} y={PT + 10} fill="#22d3ee" fontSize={9}>F:{truncF}</text>
        </g>
      )}
      {truncR > 0 && truncR <= maxPos && (
        <g>
          <line x1={toX(truncR)} x2={toX(truncR)} y1={PT} y2={H - PB}
            stroke="#818cf8" strokeWidth={1.5} strokeDasharray="5,3" opacity={0.9} />
          <text x={toX(truncR) + 3} y={PT + 22} fill="#818cf8" fontSize={9}>R:{truncR}</text>
        </g>
      )}
      {/* Legend */}
      <circle cx={W - 110} cy={PT + 8} r={3} fill="#22d3ee" />
      <text x={W - 104} y={PT + 12} fill="#94a3b8" fontSize={9}>Forward</text>
      <circle cx={W - 60} cy={PT + 8} r={3} fill="#818cf8" />
      <text x={W - 54} y={PT + 12} fill="#94a3b8" fontSize={9}>Reverse</text>
      {/* Title */}
      <text x={W / 2} y={PT + 8} fill="#64748b" fontSize={9} textAnchor="middle">
        Simulated Read Quality — Phred Score by Position
      </text>
    </svg>
  );
}

// ── Component ─────────────────────────────────────────────────────────────────
interface Props {
  params:     PipelineParams;
  onChange:   (p: PipelineParams) => void;
  marker:     MarkerType;
  onMarker?:  (m: MarkerType) => void;
  onMetadata?: (path: string, groupCol: string) => void;
}

export default function PipelineSettings({ params, onChange, marker, onMarker, onMetadata }: Props) {
  const [openStep,       setOpenStep]       = useState<number | null>(1);
  const [showBrowser,    setShowBrowser]    = useState(false);
  const [dbPaths,        setDbPaths]        = useState<Record<string, string>>({});
  const [dbDir,          setDbDir]          = useState<string>("");
  // Local primer selection — gives instant visual feedback without waiting for parent re-render
  const [selectedPrimerF, setSelectedPrimerF] = useState<string>(params.primer_f);

  // Sync local state if params.primer_f is changed from outside (e.g. reset)
  useEffect(() => {
    setSelectedPrimerF(params.primer_f);
  }, [params.primer_f]);

  const set = (key: keyof PipelineParams, val: any) => onChange({ ...params, [key]: val });
  const toggle = <T,>(arr: T[], val: T): T[] =>
    arr.includes(val) ? arr.filter(x => x !== val) : [...arr, val];

  const selectPreset = (
    f: string,
    r: string,
    extras?: Partial<Pick<PipelineParams,
      "truncLen_F" | "truncLen_R" | "trainAmpliconMinLen" | "trainAmpliconMaxLen">>
  ) => {
    setSelectedPrimerF(f);
    onChange({ ...params, primer_f: f, primer_r: r, ...(extras ?? {}) });
  };

  useEffect(() => {
    fetch("/databases").then(r => r.ok ? r.json() : null).then(data => {
      if (data?.databases) setDbPaths(data.databases);
      if (data?.db_dir)    setDbDir(data.db_dir);
    }).catch(() => {});
  }, []);

  const isITS     = marker === "ITS1" || marker === "ITS2";
  const isCOX1    = marker === "COX1";
  const isPacBio  = marker === "PacBio";
  const isONT     = marker === "ONT-16S";
  const isNema    = marker === "18S-nema";
  const is16S     = marker === "16S";
  const isStandard = !isITS && !isCOX1 && !isPacBio && !isONT;
  const hasPhylo  = is16S || marker === "12S" || isPacBio;

  const availableDBs = DB_OPTIONS.filter(db => (db.marker as string[]).includes(marker));
  const resolveDbPath = (opt: typeof DB_OPTIONS[0]) =>
    opt.value === "custom" ? params.dbPath : (dbPaths[opt.pathKey] || "");

  const pacBioPrimers = PACBIO_PRIMERS[params.pb_region] ?? PACBIO_PRIMERS["V1-V9"];
  const setPbRegion = (region: string) => {
    const pr = (PACBIO_PRIMERS[region] ?? PACBIO_PRIMERS["V1-V9"])[0];
    onChange({ ...params, pb_region: region, primer_f: pr.f, primer_r: pr.r });
  };

  // Filtered cutadapt presets for current marker
  const cutadaptPresets = CUTADAPT_PRESETS.filter(p => p.markers.includes(marker));
  const allPresets = cutadaptPresets.length > 0 ? cutadaptPresets : CUTADAPT_PRESETS;

  // Derive selected region from current primer_f (no extra state needed)
  const activePreset = CUTADAPT_PRESETS.find(p => p.f === params.primer_f && p.markers.includes(marker));
  const selectedRegion = activePreset?.region ?? "";

  const renderDbOptions = (hint: string) => (
    <div className="param-item" style={{ marginTop: 12 }}>
      <label className="param-label">Taxonomy Database</label>
      <span className="param-hint">{hint}</span>
      <div className="db-options">
        {availableDBs.map(db => {
          const realPath = resolveDbPath(db);
          const isSelected = params.taxDatabase === db.value;
          return (
            <label key={db.value} className={`db-option ${isSelected ? "active" : ""}`} style={{ cursor: "pointer", display: "block" }}>
              <input type="radio" name="taxDatabase" value={db.value} checked={isSelected}
                onChange={() => { set("taxDatabase", db.value); if (db.value !== "custom" && realPath) set("dbPath", realPath); }}
                style={{ marginRight: 8 }} />
              <span className="db-label">{db.label}</span>
              {isSelected && db.value !== "custom" && (
                <div className="db-path-info" style={{ color: realPath ? "#10b981" : "#ef4444", fontSize: 12 }}>
                  {realPath ? <code className="db-path">{realPath}</code> : "⚠ ไม่พบไฟล์ — รัน download_databases.sh ก่อน"}
                </div>
              )}
            </label>
          );
        })}
      </div>
      {params.taxDatabase === "custom" && (
        <div className="param-item" style={{ marginTop: 8 }}>
          <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
            <input type="text" className="param-input" style={{ flex: 1 }}
              placeholder={`${dbDir || "~/r16s-app/backend/databases"}/mydb.fa.gz`}
              value={params.dbPath} onChange={e => set("dbPath", e.target.value)} />
            <button className="browse-btn" onClick={() => setShowBrowser(true)}>📂 Browse</button>
          </div>
        </div>
      )}
      {showBrowser && (
        <DbFileBrowser onSelect={path => { set("dbPath", path); }} onClose={() => setShowBrowser(false)} />
      )}
    </div>
  );

  // ── Standard pipeline steps ─────────────────────────────────────────────────
  const STEPS_STANDARD = [
    { id: 1, icon: "✂️",  title: "Cutadapt — Primer Trimming",    color: "#06b6d4" },
    { id: 2, icon: "📊",  title: "Demux QC — Quality Review",      color: "#8b5cf6" },
    { id: 3, icon: "🔬",  title: "DADA2 — ASV Denoising",          color: "#ec4899" },
    { id: 4, icon: "🌳",  title: "Phylogenetic Tree",               color: "#10b981" },
    { id: 5, icon: "🏷️", title: "Taxonomy Classification",         color: "#f59e0b" },
    { id: 6, icon: "🌐",  title: "Diversity Analysis",              color: "#3b82f6" },
    { id: 7, icon: "📈",  title: "Differential Abundance",          color: "#ef4444" },
    { id: 8, icon: "🔗",  title: "Functional Analysis",             color: "#a855f7" },
  ];

  return (
    <div className="settings-wrap">

      {/* ── Marker selector ────────────────────────────────────────────── */}
      {onMarker && (
        <div className="marker-selector-wrap">
          <p className="settings-hint">Select the amplicon marker / pipeline type:</p>
          <div className="marker-grid">
            {MARKER_OPTIONS.map(m => (
              <button key={m.value}
                className={`marker-btn ${marker === m.value ? "active" : ""}`}
                onClick={() => onMarker(m.value)} title={m.description}>
                <span className="marker-icon">{m.icon}</span>
                <span className="marker-label">{m.label}</span>
              </button>
            ))}
          </div>
          <p className="marker-desc">{MARKER_OPTIONS.find(m => m.value === marker)?.description}</p>
        </div>
      )}

      {/* ── PacBio specific ────────────────────────────────────────────── */}
      {isPacBio && (
        <div className="ext-section">
          <h4 className="ext-section-title">PacBio CCS Long-Read Settings</h4>
          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">16S Region</label>
              <select className="param-input" value={params.pb_region} onChange={e => setPbRegion(e.target.value)}>
                <option value="V1-V9">V1–V9 (full-length, ~1500 bp)</option>
                <option value="V1-V3">V1–V3 (~450 bp)</option>
                <option value="V4">V4 (~250 bp)</option>
                <option value="custom">Custom</option>
              </select>
            </div>
            <ParamNumber label="Min Read Length (bp)" hint="Minimum for filtering"
              value={params.pb_min_len} min={100} max={2000} onChange={v => set("pb_min_len", v)} />
            <ParamNumber label="Max Read Length (bp)" hint="Maximum for filtering"
              value={params.pb_max_len} min={100} max={2500} onChange={v => set("pb_max_len", v)} />
            <ParamNumber label="Max EE" hint="CCS reads are high accuracy; 2–5 typical"
              value={params.pb_maxEE} min={1} max={20} step={0.5} onChange={v => set("pb_maxEE", v)} />
          </div>
          <div className="primer-section">
            <label className="param-label">Primer Pair</label>
            <div className="primer-presets">
              {pacBioPrimers.map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => onChange({ ...params, primer_f: pr.f, primer_r: pr.r })}>
                  {pr.label}
                </button>
              ))}
            </div>
          </div>
          {renderDbOptions("SILVA recommended for full-length 16S")}
          <div className="param-grid" style={{ marginTop: 12 }}>
            <ParamNumber label="Min Bootstrap (%)" hint="Minimum confidence for taxonomy"
              value={params.minBoot} min={0} max={100} step={5} onChange={v => set("minBoot", v)} />
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

      {/* ── ONT-16S specific ──────────────────────────────────────────── */}
      {isONT && (
        <div className="ext-section">
          <h4 className="ext-section-title">🧫 ONT 16S Settings — Emu Pipeline</h4>

          {/* Info box */}
          <div className="ps-info-box ps-info-box--ok" style={{ marginBottom: 16 }}>
            ✅ ONT reads (single FASTQ per sample) → Emu abundance → same visualizations as Illumina
          </div>

          <div className="param-grid">
            {/* Region selector */}
            <div className="param-item">
              <label className="param-label">16S Region</label>
              <span className="param-hint">V7-V8 = ~337 bp (genus-level) · V1-V9 = ~1500 bp (species-level)</span>
              <select className="param-input" value={params.ont_region}
                onChange={e => {
                  const region = e.target.value;
                  // Auto-fill primer presets
                  if (region === "V7-V8") {
                    onChange({ ...params, ont_region: region, primer_f: "AACMGGATTAGATACCCKG", primer_r: "ACGTCATCCCCACCTTCC" });
                  } else if (region === "V1-V9") {
                    onChange({ ...params, ont_region: region, primer_f: "AGRGTTYGATYMTGGCTCAG", primer_r: "RGYTACCTTGTTACGACTT" });
                  } else {
                    onChange({ ...params, ont_region: region });
                  }
                }}>
                <option value="V1-V9">V1–V9 (full-length, ~1500 bp, species-level)</option>
                <option value="V7-V8">V7–V8 (~337 bp, genus-level)</option>
                <option value="custom">Custom</option>
              </select>
            </div>

            {/* Min abundance */}
            <div className="param-item">
              <label className="param-label">Min Abundance Filter</label>
              <span className="param-hint">Discard taxa below this relative abundance (0.0001 = 0.01%)</span>
              <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
                <input type="range" min={0.00001} max={0.01} step={0.00001}
                  value={params.ont_min_abundance}
                  onChange={e => set("ont_min_abundance", Number(e.target.value))}
                  style={{ flex: 1, accentColor: "#22d3ee" }} />
                <span style={{ color: "#22d3ee", fontWeight: 700, width: 60, textAlign: "right", fontSize: 12 }}>
                  {params.ont_min_abundance.toExponential(1)}
                </span>
              </div>
            </div>
          </div>

          {/* Primer presets */}
          <div className="primer-section" style={{ marginTop: 12 }}>
            <label className="param-label">Primer Pair</label>
            <div className="primer-presets">
              <button className={`primer-preset-btn ${params.primer_f === "AGRGTTYGATYMTGGCTCAG" ? "active" : ""}`}
                onClick={() => onChange({ ...params, primer_f: "AGRGTTYGATYMTGGCTCAG", primer_r: "RGYTACCTTGTTACGACTT" })}>
                27F / 1492R (V1–V9)
              </button>
              <button className={`primer-preset-btn ${params.primer_f === "AACMGGATTAGATACCCKG" ? "active" : ""}`}
                onClick={() => onChange({ ...params, primer_f: "AACMGGATTAGATACCCKG", primer_r: "ACGTCATCCCCACCTTCC" })}>
                1055F / 1392R (V7–V8, ~337 bp)
              </button>
              <button className={`primer-preset-btn ${!params.primer_f ? "active" : ""}`}
                onClick={() => onChange({ ...params, primer_f: "", primer_r: "" })}>
                No primers (skip cutadapt)
              </button>
            </div>
          </div>

          {/* Primer inputs */}
          <div className="param-grid" style={{ marginTop: 8 }}>
            <div className="param-item">
              <label className="param-label">Forward Primer (5'→3')</label>
              <input type="text" className="param-input"
                style={{ fontFamily: "monospace", letterSpacing: "0.5px" }}
                placeholder="e.g. AGRGTTYGATYMTGGCTCAG"
                value={params.primer_f}
                onChange={e => set("primer_f", e.target.value.toUpperCase())} />
            </div>
            <div className="param-item">
              <label className="param-label">Reverse Primer (5'→3')</label>
              <input type="text" className="param-input"
                style={{ fontFamily: "monospace", letterSpacing: "0.5px" }}
                placeholder="e.g. RGYTACCTTGTTACGACTT"
                value={params.primer_r}
                onChange={e => set("primer_r", e.target.value.toUpperCase())} />
            </div>
          </div>

          {/* Emu database */}
          <div className="param-item" style={{ marginTop: 12 }}>
            <label className="param-label">Emu Database</label>
            <span className="param-hint">เลือก Emu database directory (ต้องมี species_taxid.fasta + taxonomy.tsv)</span>
            {(() => {
              // Collect all emu_* entries from db_paths.json
              const emuDbs: { key: string; path: string; label: string }[] = Object.entries(dbPaths)
                .filter(([k, v]) => k.startsWith("emu_") && v)
                .map(([k, v]) => ({
                  key: k,
                  path: v as string,
                  label: k === "emu_silva" ? "Emu SILVA (full-length)"
                       : k === "emu_db_mar2026" ? "Emu DB Mar 2026"
                       : k.replace("emu_", "Emu "),
                }));

              const currentVal = params.ont_db_path || (emuDbs[0]?.path ?? "");
              const isKnown = emuDbs.some(d => d.path === currentVal);

              if (emuDbs.length === 0) {
                return (
                  <div style={{ fontSize: 12, color: "#f59e0b", marginBottom: 6 }}>
                    ⚠ ไม่พบ Emu database — ดูวิธีสร้างด้านล่าง
                  </div>
                );
              }

              return (
                <>
                  <select
                    className="param-input"
                    style={{ width: "100%", marginBottom: 6 }}
                    value={isKnown ? currentVal : "__custom__"}
                    onChange={e => {
                      set("ont_db_path", e.target.value === "__custom__" ? "" : e.target.value);
                    }}
                  >
                    {emuDbs.map(d => (
                      <option key={d.key} value={d.path}>{d.label} — {d.path.split("/").slice(-2).join("/")}</option>
                    ))}
                    <option value="__custom__">Custom path...</option>
                  </select>
                  {!isKnown && (
                    <input type="text" className="param-input" style={{ width: "100%", marginBottom: 4 }}
                      placeholder="Path to Emu database directory"
                      autoFocus
                      value={currentVal}
                      onChange={e => set("ont_db_path", e.target.value)} />
                  )}
                  {currentVal && (
                    <div style={{ fontSize: 11, color: "#10b981" }}>
                      ✅ {currentVal}
                    </div>
                  )}
                </>
              );
            })()}
          </div>

          {/* Region-specific Emu DB hint */}
          {params.ont_region === "V7-V8" && (
            <div style={{ marginTop: 14, border: "1px solid #854d0e", borderRadius: 8, padding: 14, background: "#1c1408" }}>
              <div style={{ fontWeight: 700, fontSize: 12, color: "#fbbf24", marginBottom: 6 }}>
                💡 V7-V8 Region-Specific Database (Recommended)
              </div>
              <p className="param-hint" style={{ marginBottom: 10 }}>
                Emu default database ใช้ full-length SILVA (~1500 bp) ซึ่งทำให้ alignment score ต่ำเมื่อ read ยาวแค่ 300–450 bp
                ควร build V7-V8-specific database เพื่อ accuracy สูงขึ้น
              </p>
              <div className="ps-cmd-block">
                <div className="ps-cmd-header"><span>build V7-V8 Emu database</span></div>
                <pre className="ps-cmd-body">{`# ขั้นตอนที่ 1: สร้าง intermediate files (ตัด V7-V8 region ออกจาก SILVA)
python3 ~/r16s-app/backend/python_scripts/build_emu_db.py \\
  --silva-db ~/r16s-app/backend/databases/SILVA/silva_nr99_v138.1_train_set.fa.gz \\
  --output-dir ~/r16s-app/backend/databases/emu_v7v8 \\
  --region V7-V8

# ขั้นตอนที่ 2: build Emu database จาก intermediate files (ต้อง cd เข้าไปก่อน)
cd ~/r16s-app/backend/databases/emu_v7v8 && \\
conda run -n emu emu build-database \\
  --sequences  sequences.fasta \\
  --seq2tax    seq2taxid.tsv \\
  --taxonomy-list taxonomy.tsv \\
  silva_v7v8

# ขั้นตอนที่ 3: ใส่ path ด้านบน → ~/r16s-app/backend/databases/emu_v7v8/silva_v7v8`}</pre>
              </div>
              <div style={{ fontSize: 11, color: "#78716c", marginTop: 8 }}>
                หรือใช้ custom primer อื่นได้: เพิ่ม --primer-f XXXXXXXX --primer-r XXXXXXXX --min-len 250 --max-len 450
              </div>
            </div>
          )}

          {/* TopN + CPU Threads */}
          <div className="param-grid" style={{ marginTop: 12 }}>
            <div className="param-item">
              <label className="param-label">Top Taxa to Show in Plots</label>
              <div className="topn-toggle">
                {([30,50,100] as const).map(n => (
                  <button key={n} className={`topn-btn ${params.topN === n ? "active" : ""}`}
                    onClick={() => set("topN", n)}>Top {n}</button>
                ))}
              </div>
            </div>
            <ParamNumber label="CPU Threads" hint="Cores for cutadapt + Emu minimap2 (default 4)"
              value={params.nThreads} min={1} max={128} onChange={v => set("nThreads", v)} />
          </div>

          {/* Pipeline summary */}
          <div className="ps-cmd-block" style={{ marginTop: 12 }}>
            <div className="ps-cmd-header"><span>emu_pipeline.py — steps</span></div>
            <pre className="ps-cmd-body">{`1. Scan input FASTQ files (single file per sample)
2. Trim primers with cutadapt${params.primer_f ? ` (${params.primer_f.slice(0,10)}… / ${params.primer_r.slice(0,10)}…)` : " — SKIPPED (no primers set)"}
3. Run Emu abundance on each sample
4. Combine per-sample TSVs → asv_table.csv + taxonomy.csv
5. Write read_tracking.csv + summary.json
6. Generate visualizations (viz_pipeline.R)`}</pre>
          </div>

          {/* Metadata */}
          <div className="ext-section" style={{ marginTop: 16, borderTop: "1px solid #1e293b", paddingTop: 16 }}>
            <h4 className="ext-section-title">
              📊 Metadata{" "}
              <span style={{ fontWeight: 400, fontSize: 13, color: "#6b7280" }}>
                (optional — for diversity grouping)
              </span>
            </h4>
            <MetadataUpload
              inline
              onMetadataReady={(info: MetadataInfo, grp: string) => {
                onChange({ ...params, metadataPath: info.path, groupCol: grp });
                onMetadata?.(info.path, grp);
              }}
              onGroupColChange={(col: string) => onChange({ ...params, groupCol: col })}
            />
            {params.metadataPath && (
              <p style={{ fontSize: 12, color: "#10b981", marginTop: 8 }}>
                ✅ Metadata loaded · Group column: <strong>{params.groupCol}</strong>
              </p>
            )}
          </div>
        </div>
      )}

      {/* ── ITS specific ───────────────────────────────────────────────── */}
      {isITS && (
        <div className="ext-section">
          <h4 className="ext-section-title">ITS Fungal Settings</h4>
          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">ITS Region</label>
              <div className="topn-toggle">
                {(["ITS1","ITS2"] as const).map(r => (
                  <button key={r} className={`topn-btn ${params.its_region === r ? "active" : ""}`}
                    onClick={() => set("its_region", r)}>{r}</button>
                ))}
              </div>
            </div>
            <ParamNumber label="Max EE Forward" hint="Maximum expected errors (forward)"
              value={params.maxEE_F} min={1} max={10} step={0.5} onChange={v => set("maxEE_F", v)} />
            <ParamNumber label="Max EE Reverse" hint="Maximum expected errors (reverse)"
              value={params.maxEE_R} min={1} max={10} step={0.5} onChange={v => set("maxEE_R", v)} />
          </div>
          <div className="primer-section">
            <label className="param-label">Primers</label>
            <div className="primer-presets">
              {(ITS_PRIMERS[params.its_region] || []).map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => onChange({ ...params, primer_f: pr.f, primer_r: pr.r })}>
                  {pr.label}
                </button>
              ))}
              <button className="primer-preset-btn" onClick={() => onChange({ ...params, primer_f: "", primer_r: "" })}>
                Use defaults
              </button>
            </div>
            {params.primer_f && (
              <div className="primer-inputs">
                <div className="param-item">
                  <label className="param-label">Forward primer (5'→3')</label>
                  <input type="text" className="param-input" value={params.primer_f}
                    onChange={e => set("primer_f", e.target.value)} />
                </div>
                <div className="param-item">
                  <label className="param-label">Reverse primer (5'→3')</label>
                  <input type="text" className="param-input" value={params.primer_r}
                    onChange={e => set("primer_r", e.target.value)} />
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ── COX1 specific ──────────────────────────────────────────────── */}
      {isCOX1 && (
        <div className="ext-section">
          <h4 className="ext-section-title">COX1 Metabarcoding Settings</h4>
          <div className="param-grid">
            <ParamNumber label="Truncate Forward (bp)" hint="Truncate forward reads"
              value={params.truncLen_cox1_f} min={0} max={350} onChange={v => set("truncLen_cox1_f", v)} />
            <ParamNumber label="Truncate Reverse (bp)" hint="0 = no truncation"
              value={params.truncLen_cox1_r} min={0} max={350} onChange={v => set("truncLen_cox1_r", v)} />
            <ParamNumber label="Max EE Forward" hint="Maximum expected errors (forward)"
              value={params.maxEE_F} min={1} max={10} step={0.5} onChange={v => set("maxEE_F", v)} />
            <ParamNumber label="Max EE Reverse" hint="Maximum expected errors (reverse)"
              value={params.maxEE_R} min={1} max={10} step={0.5} onChange={v => set("maxEE_R", v)} />
            <ParamNumber label="Min ASV Length (bp)" hint="Typical mlCOI: 300 bp"
              value={params.cox1_min_len} min={100} max={500} onChange={v => set("cox1_min_len", v)} />
            <ParamNumber label="Max ASV Length (bp)" hint="Remove ASVs longer than this"
              value={params.cox1_max_len} min={100} max={500} onChange={v => set("cox1_max_len", v)} />
          </div>
          <div className="param-grid">
            <div className="param-item">
              <label className="param-label">Genetic Code</label>
              <select className="param-input" value={params.codon_table}
                onChange={e => set("codon_table", Number(e.target.value))}>
                <option value={5}>5 — Invertebrate Mitochondrial (default)</option>
                <option value={2}>2 — Vertebrate Mitochondrial</option>
                <option value={1}>1 — Standard (nuclear)</option>
              </select>
            </div>
            <div className="param-item">
              <label className="param-label">LULU Post-curation</label>
              <span className="param-hint">BLAST-based removal of artefact ASVs</span>
              <div className="topn-toggle">
                <button className={`topn-btn ${params.run_lulu ? "active" : ""}`}
                  onClick={() => set("run_lulu", true)}>Enabled</button>
                <button className={`topn-btn ${!params.run_lulu ? "active" : ""}`}
                  onClick={() => set("run_lulu", false)}>Disabled</button>
              </div>
            </div>
          </div>
          <div className="primer-section">
            <label className="param-label">Primer Pair</label>
            <div className="primer-presets">
              {COX1_PRIMERS.map((pr, i) => (
                <button key={i}
                  className={`primer-preset-btn ${params.primer_f === pr.f ? "active" : ""}`}
                  onClick={() => { set("primer_f", pr.f); set("primer_r", pr.r); }}>
                  {pr.label}
                </button>
              ))}
            </div>
          </div>
        </div>
      )}

      {/* ══════════════════════════════════════════════════════════════════
          STANDARD PIPELINE STEPS (16S / 12S / 18S-nema)
         ══════════════════════════════════════════════════════════════════ */}
      {isStandard && (
        <>
          {STEPS_STANDARD.map(s => (
            <div key={s.id} className="step-accordion">
              <div className="step-acc-header"
                style={{ borderLeft: `4px solid ${s.color}` }}
                onClick={() => setOpenStep(openStep === s.id ? null : s.id)}>
                <span className="step-acc-icon">{s.icon}</span>
                <span className="step-acc-title">Step {s.id} — {s.title}</span>
                {/* Skip badges */}
                {s.id === 4 && !params.runPhylogeny && (
                  <span className="step-skip-badge">SKIP</span>
                )}
                {s.id === 6 && !params.runDiversity && (
                  <span className="step-skip-badge">SKIP</span>
                )}
                {s.id === 7 && !params.runDiffAbund && (
                  <span className="step-skip-badge">SKIP</span>
                )}
                <span className="step-acc-chev">{openStep === s.id ? "▲" : "▼"}</span>
              </div>

              {openStep === s.id && (
                <div className="step-acc-body">

                  {/* ── Step 1: Cutadapt ─────────────────────────────── */}
                  {s.id === 1 && (
                    <div>
                      <p className="param-hint" style={{ marginBottom: 12 }}>
                        Remove primers from reads before DADA2. Leave primers blank to skip cutadapt.
                      </p>

                      {/* Sequencer type toggle */}
                      <div className="ps-section-label">SEQUENCER TYPE</div>
                      <div style={{ display: "flex", gap: 8, marginBottom: 10 }}>
                        <button type="button"
                          style={{
                            padding: "6px 14px", borderRadius: 6, fontSize: 13, cursor: "pointer",
                            border: params.sequencerType !== "ont" ? "2px solid #6366f1" : "1px solid #334155",
                            background: params.sequencerType !== "ont" ? "#1e1b4b" : "#1e293b",
                            color: params.sequencerType !== "ont" ? "#a5b4fc" : "#94a3b8",
                            fontWeight: params.sequencerType !== "ont" ? 600 : 400,
                          }}
                          onClick={() => set("sequencerType", "illumina")}>
                          🔬 Illumina (paired-end)
                        </button>
                        <button type="button"
                          style={{
                            padding: "6px 14px", borderRadius: 6, fontSize: 13, cursor: "pointer",
                            border: params.sequencerType === "ont" ? "2px solid #f59e0b" : "1px solid #334155",
                            background: params.sequencerType === "ont" ? "#1c1200" : "#1e293b",
                            color: params.sequencerType === "ont" ? "#fcd34d" : "#94a3b8",
                            fontWeight: params.sequencerType === "ont" ? 600 : 400,
                          }}
                          onClick={() => onChange({
                            ...params, sequencerType: "ont",
                            maxEE_F: params.maxEE_F < 5 ? 5 : params.maxEE_F,
                          })}>
                          🧬 ONT R10.4+ (single-end)
                        </button>
                      </div>
                      {params.sequencerType === "ont" && (
                        <div style={{
                          background: "#1c1200", border: "1px solid #92400e",
                          borderRadius: 6, padding: "8px 12px", marginBottom: 12, fontSize: 12,
                          color: "#fcd34d", lineHeight: 1.6,
                        }}>
                          ⚠️ <strong>ONT mode</strong> — reads processed as single-end (no R2, no merge step),
                          filtered by <strong>length range</strong> instead of truncLen (ONT reads vary in length),
                          and denoised with indel-tolerant settings (BAND_SIZE=32, homopolymer-lenient).
                          Requires <strong>R10.4.1 flowcell</strong> + <strong>Q20+ basecalling</strong>.
                          Recommended for short amplicons ≤500 bp (e.g. V7–V8 ~337 bp). DADA2 still expects
                          much lower error rates than raw ONT reads — for best accuracy prefer the
                          <strong> ONT 16S (Emu)</strong> marker instead when available.
                        </div>
                      )}

                      {/* Primer presets */}
                      <div className="ps-section-label">PRIMER PRESETS</div>
                      <div className="ps-primer-grid">
                        {allPresets.map((pr, i) => (
                          <button key={i} type="button"
                            className={`ps-primer-card ${selectedPrimerF === pr.f ? "ps-primer-card--active" : ""}`}
                            onClick={() => selectPreset(pr.f, pr.r, {
                              ...(pr.truncF    !== undefined ? { truncLen_F: pr.truncF, truncLen_R: pr.truncR ?? params.truncLen_R } : {}),
                              ...(pr.ampMinLen !== undefined ? { trainAmpliconMinLen: pr.ampMinLen, trainAmpliconMaxLen: pr.ampMaxLen ?? params.trainAmpliconMaxLen } : {}),
                              ...(pr.ampMinLen !== undefined ? { ontMinLen: pr.ampMinLen, ontMaxLen: pr.ampMaxLen ?? params.ontMaxLen } : {}),
                            })}>
                            <div className="ps-primer-name">{pr.label}</div>
                            <div className="ps-primer-seq">F: {pr.f.length > 18 ? pr.f.slice(0, 18) + "…" : pr.f}</div>
                            <div className="ps-primer-seq">R: {pr.r.length > 18 ? pr.r.slice(0, 18) + "…" : pr.r}</div>
                          </button>
                        ))}
                        <button type="button"
                          className={`ps-primer-card ${!selectedPrimerF ? "ps-primer-card--active" : ""}`}
                          onClick={() => selectPreset("", "")}>
                          <div className="ps-primer-name">Skip</div>
                          <div className="ps-primer-seq">No primer removal</div>
                        </button>
                        <button type="button"
                          className={`ps-primer-card ps-primer-card--custom ${
                            selectedPrimerF === "__custom__" || (selectedPrimerF && !allPresets.some(p => p.f === selectedPrimerF))
                              ? "ps-primer-card--active" : ""
                          }`}
                          onClick={() => { setSelectedPrimerF("__custom__"); }}>
                          <div className="ps-primer-name">✏️ Custom</div>
                          <div className="ps-primer-seq">Define your own primers</div>
                        </button>
                      </div>

                      {/* Primer inputs */}
                      <div className="ps-section-label" style={{ marginTop: 16 }}>PRIMER SETTINGS</div>
                      <div className="param-grid">
                        <div className="param-item">
                          <label className="param-label">Forward Primer (5'→3')</label>
                          <input type="text" className="param-input"
                            style={{ fontFamily: "monospace", letterSpacing: "0.5px" }}
                            placeholder="e.g. CCTACGGGNGGCWGCAG"
                            value={params.primer_f}
                            onChange={e => set("primer_f", e.target.value.toUpperCase())} />
                        </div>
                        <div className="param-item">
                          <label className="param-label">Reverse Primer (5'→3')</label>
                          <input type="text" className="param-input"
                            style={{ fontFamily: "monospace", letterSpacing: "0.5px" }}
                            placeholder="e.g. GACTACHVGGGTATCTAATCC"
                            value={params.primer_r}
                            onChange={e => set("primer_r", e.target.value.toUpperCase())} />
                        </div>
                      </div>

                      <div className="param-grid" style={{ marginTop: 4 }}>
                        <div className="param-item">
                          <label className="param-label">Error Rate</label>
                          <span className="param-hint">Fraction of mismatches allowed in primer match</span>
                          <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
                            <input type="range" min={0.05} max={0.30} step={0.01}
                              value={params.cutadaptErrorRate}
                              onChange={e => set("cutadaptErrorRate", Number(e.target.value))}
                              style={{ flex: 1, accentColor: "#22d3ee" }} />
                            <span style={{ color: "#22d3ee", fontWeight: 700, width: 36, textAlign: "right" }}>
                              {params.cutadaptErrorRate.toFixed(2)}
                            </span>
                          </div>
                        </div>
                        <ParamNumber label="Min Overlap (bp)" hint="Minimum bases that must match primer"
                          value={params.cutadaptOverlap} min={1} max={30}
                          onChange={v => set("cutadaptOverlap", v)} />
                      </div>

                      <label className="ps-checkbox-row" style={{ marginTop: 8 }}>
                        <input type="checkbox" checked={params.discardUntrimmed}
                          onChange={e => set("discardUntrimmed", e.target.checked)} />
                        <span>Discard untrimmed reads (reads where primer was not found)</span>
                      </label>

                      {/* Command preview */}
                      {params.primer_f && (
                        <div className="ps-cmd-block">
                          <div className="ps-cmd-header">
                            <span>cutadapt trim-paired</span>
                          </div>
                          <pre className="ps-cmd-body">{`cutadapt \\
  --front ${params.primer_f} \\
  --adapter ${params.primer_r} \\
  --error-rate ${params.cutadaptErrorRate} \\
  --overlap ${params.cutadaptOverlap}${params.discardUntrimmed ? " \\\n  --discard-untrimmed" : ""} \\
  -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \\
  input_R1.fastq.gz input_R2.fastq.gz`}</pre>
                        </div>
                      )}
                    </div>
                  )}

                  {/* ── Step 2: Demux QC ─────────────────────────────── */}
                  {s.id === 2 && (
                    <div>
                      <p className="param-hint" style={{ marginBottom: 12 }}>
                        Inspect read quality before setting DADA2 truncation lengths.
                        The plot below is simulated — use your actual .qzv file at{" "}
                        <a href="https://view.qiime2.org" target="_blank" rel="noreferrer"
                          style={{ color: "#22d3ee" }}>view.qiime2.org</a>{" "}
                        for real quality scores. Set truncLen values in Step 3 based on where quality drops below Q25.
                      </p>

                      <div className="ps-section-label">SIMULATED QUALITY PLOT</div>
                      <div style={{ marginBottom: 12 }}>
                        <QualityPlot truncF={params.truncLen_F} truncR={params.truncLen_R} />
                      </div>

                      <div className="ps-demux-tips">
                        {params.sequencerType === "ont" ? (
                          <>
                            <div className="ps-tip">🧬 <strong>ONT single-end</strong> — only forward read quality shown, no truncLen (length-range filter instead)</div>
                            <div className="ps-tip">💡 Set Min/Max Length to the expected amplicon size (e.g. 300–400 bp for V7–V8)</div>
                            <div className="ps-tip">💡 Recommend <strong>Max EE ≥ 5</strong> for ONT (more lenient than Illumina — raw ONT error is much higher)</div>
                            <div className="ps-tip" style={{ color: "#fcd34d" }}>
                              Current: length filter = {params.ontMinLen}–{params.ontMaxLen} bp (single-end, no merge needed)
                            </div>
                          </>
                        ) : (
                          <>
                            <div className="ps-tip">💡 <strong>Forward reads</strong> usually stay high quality to ~250 bp</div>
                            <div className="ps-tip">💡 <strong>Reverse reads</strong> often drop at ~150–200 bp — truncate here</div>
                            <div className="ps-tip">💡 <strong>Overlap</strong>: truncLen_F + truncLen_R must exceed amplicon length by ≥20 bp to merge</div>
                            <div className="ps-tip" style={{ color: "#22d3ee" }}>
                              Current: F={params.truncLen_F} + R={params.truncLen_R} = {params.truncLen_F + params.truncLen_R} bp
                              {selectedRegion === "V3-V4" && " (16S V3–V4 amplicon ~460 bp → need ≥480 combined)"}
                              {selectedRegion === "V4"    && " (16S V4 amplicon ~253 bp → need ≥273 combined)"}
                              {selectedRegion === "V7-V8" && " (16S V7–V8 amplicon ~337 bp → need ≥357 combined)"}
                              {!selectedRegion && marker === "16S" && " (16S V3–V4 amplicon ~460 bp → need ≥480 combined)"}
                            </div>
                          </>
                        )}
                      </div>
                    </div>
                  )}

                  {/* ── Step 3: DADA2 ────────────────────────────────── */}
                  {s.id === 3 && (
                    <div className="param-grid">
                      {params.sequencerType === "ont" && (
                        <div style={{
                          gridColumn: "1/-1", background: "#1c1200",
                          border: "1px solid #92400e", borderRadius: 6,
                          padding: "8px 12px", marginBottom: 4, fontSize: 12, color: "#fcd34d",
                        }}>
                          🧬 <strong>ONT single-end mode</strong> — R2 settings are disabled. Length-range filter replaces truncLen (ONT reads vary in length).
                        </div>
                      )}
                      {params.sequencerType === "ont" ? (
                        <>
                          <ParamNumber label="Min Length (bp)" hint="Discard reads shorter than this (expected amplicon size)"
                            value={params.ontMinLen} min={50} max={2000} onChange={v => set("ontMinLen", v)} />
                          <ParamNumber label="Max Length (bp)" hint="Discard reads longer than this (expected amplicon size)"
                            value={params.ontMaxLen} min={50} max={2000} onChange={v => set("ontMaxLen", v)} />
                        </>
                      ) : (
                        <>
                          <ParamNumber label="Truncate Forward (bp)" hint="Truncate at this position (set 0 to disable)"
                            value={params.truncLen_F} min={0} max={500} onChange={v => set("truncLen_F", v)} />
                          <ParamNumber label="Truncate Reverse (bp)" hint="Truncate at this position (set 0 to disable)"
                            value={params.truncLen_R} min={0} max={350} onChange={v => set("truncLen_R", v)} />
                        </>
                      )}
                      <ParamNumber label="Max EE Forward" hint="Maximum expected errors — forward reads"
                        value={params.maxEE_F} min={1} max={params.sequencerType === "ont" ? 30 : 10} step={0.5} onChange={v => set("maxEE_F", v)} />
                      {params.sequencerType !== "ont" && (
                        <ParamNumber label="Max EE Reverse" hint="Maximum expected errors — reverse reads"
                          value={params.maxEE_R} min={1} max={10} step={0.5} onChange={v => set("maxEE_R", v)} />
                      )}
                      <ParamNumber label="Trim Left Forward (bp)" hint="Trim from 5' end (if cutadapt not used)"
                        value={params.trimLeft_F} min={0} max={50} onChange={v => set("trimLeft_F", v)} />
                      {params.sequencerType !== "ont" && (
                        <ParamNumber label="Trim Left Reverse (bp)" hint="Trim from 5' end (if cutadapt not used)"
                          value={params.trimLeft_R} min={0} max={50} onChange={v => set("trimLeft_R", v)} />
                      )}
                      <ParamSelect label="Pooling Mode" hint="How to pool samples for ASV inference"
                        value={params.pool}
                        options={[
                          { value: "FALSE",  label: "FALSE — fast (independent)" },
                          { value: "pseudo", label: "pseudo — balanced (recommended)" },
                          { value: "TRUE",   label: "TRUE — most sensitive, slowest" },
                        ]}
                        onChange={v => set("pool", v)} />
                      <ParamSelect label="Chimera Method" hint="Method for detecting chimeric sequences"
                        value={params.chimeraMethod}
                        options={[
                          { value: "consensus",  label: "consensus (recommended)" },
                          { value: "per-sample", label: "per-sample" },
                          { value: "pooled",     label: "pooled (slowest)" },
                        ]}
                        onChange={v => set("chimeraMethod", v)} />
                      <ParamSelect label="Bases to Learn From" hint="More = higher accuracy but slower"
                        value={String(params.nbases)}
                        options={[
                          { value: "1e7", label: "10M bases (fast)" },
                          { value: "1e8", label: "100M bases (default)" },
                          { value: "5e8", label: "500M bases (accurate)" },
                        ]}
                        onChange={v => set("nbases", Number(v))} />
                      <ParamNumber label="CPU Threads" hint="Parallel workers for DADA2"
                        value={params.nThreads} min={1} max={32} onChange={v => set("nThreads", v)} />
                      {isNema && (
                        <div className="param-item" style={{ gridColumn: "1/-1" }}>
                          <label className="param-label">Nematode 18S Primers</label>
                          <div className="primer-presets">
                            {NEMA_PRIMERS.map((pr, i) => (
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

                  {/* ── Step 4: Phylogeny ─────────────────────────────── */}
                  {s.id === 4 && (
                    <div>
                      {hasPhylo ? (
                        <>
                          <div className="ps-toggle-row">
                            <span className="param-label">Build Phylogenetic Tree</span>
                            <div className="topn-toggle">
                              <button className={`topn-btn ${params.runPhylogeny ? "active" : ""}`}
                                onClick={() => set("runPhylogeny", true)}>Enable</button>
                              <button className={`topn-btn ${!params.runPhylogeny ? "active" : ""}`}
                                onClick={() => set("runPhylogeny", false)}>Skip</button>
                            </div>
                          </div>
                          {params.runPhylogeny && (
                            <>
                              <div className="param-grid" style={{ marginTop: 12 }}>
                                <ParamNumber label="Threads (MAFFT)" hint="Parallel threads for alignment"
                                  value={params.phyloThreads} min={1} max={32}
                                  onChange={v => set("phyloThreads", v)} />
                              </div>
                              <div className="ps-section-label" style={{ marginTop: 16 }}>OUTPUTS GENERATED</div>
                              <div className="ps-cmd-block">
                                <pre className="ps-cmd-body">{`qza/
├── aligned-rep-seqs.qza
├── masked-aligned-rep-seqs.qza
├── unrooted-tree.qza
└── rooted-tree.qza  ← used for UniFrac, Faith PD`}</pre>
                              </div>
                              <div className="ps-cmd-block" style={{ marginTop: 8 }}>
                                <div className="ps-cmd-header"><span>align-to-tree-mafft-fasttree</span></div>
                                <pre className="ps-cmd-body">{`mafft --thread ${params.phyloThreads} --auto rep-seqs.fasta > aligned.fasta
FastTree -gtr -nt aligned.fasta > unrooted-tree.nwk
# Then mid-point rooting → rooted-tree.nwk`}</pre>
                              </div>
                            </>
                          )}
                        </>
                      ) : (
                        <div className="ps-info-box ps-info-box--warn">
                          ⚠️ Phylogenetic tree is not built for <strong>{marker}</strong> (variable-length amplicon).
                          UniFrac metrics will be skipped in diversity analysis.
                        </div>
                      )}
                    </div>
                  )}

                  {/* ── Step 5: Taxonomy ──────────────────────────────── */}
                  {s.id === 5 && (
                    <div>
                      {/* Classifier */}
                      <div className="ps-section-label">CLASSIFIER</div>
                      <div className="ps-clf-chips">
                        {availableDBs.filter(db => db.value !== "custom").map(db => (
                          <button key={db.value}
                            className={`ps-clf-chip ${params.taxDatabase === db.value ? "ps-clf-chip--active" : ""}`}
                            onClick={() => {
                              const realPath = resolveDbPath(db);
                              set("taxDatabase", db.value);
                              if (realPath) set("dbPath", realPath);
                            }}>
                            {db.label.split(" — ")[0]}
                          </button>
                        ))}
                        <button
                          className={`ps-clf-chip ${params.taxDatabase === "custom" ? "ps-clf-chip--active" : ""}`}
                          onClick={() => set("taxDatabase", "custom")}>
                          Custom
                        </button>
                      </div>

                      {params.taxDatabase === "custom" && (
                        <div style={{ marginTop: 8, display: "flex", gap: 8 }}>
                          <input type="text" className="param-input" style={{ flex: 1 }}
                            placeholder="Path to database (.fa.gz / .fasta.gz)"
                            value={params.dbPath} onChange={e => set("dbPath", e.target.value)} />
                          <button className="browse-btn" onClick={() => setShowBrowser(true)}>📂</button>
                        </div>
                      )}
                      {params.taxDatabase !== "custom" && (() => {
                        const db = availableDBs.find(d => d.value === params.taxDatabase);
                        const realPath = db ? resolveDbPath(db) : "";
                        return realPath ? (
                          <div style={{ fontSize: 12, color: "#10b981", marginTop: 6 }}>
                            ✅ {realPath.split("/").pop()}
                          </div>
                        ) : (
                          <div style={{ fontSize: 12, color: "#ef4444", marginTop: 6 }}>
                            ⚠ Database not found — run download_databases.sh
                          </div>
                        );
                      })()}

                      {showBrowser && (
                        <DbFileBrowser onSelect={p => { set("dbPath", p); }} onClose={() => setShowBrowser(false)} />
                      )}

                      {/* Custom Classifier Training */}
                      {(marker === "16S" || marker === "12S" || marker === "18S-nema") && (
                        <div style={{ marginTop: 18, border: "1px solid #334155", borderRadius: 8, padding: 14, background: "#0f172a" }}>
                          <div className="ps-section-label" style={{ marginBottom: 10 }}>
                            🧬 CUSTOM CLASSIFIER (Region-specific training)
                          </div>
                          <p className="param-hint" style={{ marginBottom: 10 }}>
                            Train a Naive Bayes classifier for your exact primer region (e.g. 300–500 bp amplicons).
                            Outperforms full-length classifiers for partial-region data.
                          </p>
                          <div style={{ display: "flex", gap: 6, marginBottom: 12 }}>
                            {(["default", "train", "upload"] as const).map(mode => (
                              <button key={mode} type="button"
                                style={{
                                  padding: "5px 14px", borderRadius: 6, fontSize: 12, fontWeight: 600, cursor: "pointer",
                                  border: params.customClassifierMode === mode ? "1px solid #6366f1" : "1px solid #334155",
                                  background: params.customClassifierMode === mode ? "#312e81" : "#1e293b",
                                  color: params.customClassifierMode === mode ? "#c7d2fe" : "#94a3b8",
                                }}
                                onClick={() => set("customClassifierMode", mode)}>
                                {mode === "default" ? "🔵 Auto-detect" : mode === "train" ? "🟣 Train from primers" : "📁 Upload .qza"}
                              </button>
                            ))}
                          </div>

                          {params.customClassifierMode === "default" && (
                            <div>
                              <div className="ps-info-box" style={{ fontSize: 12 }}>
                                Auto-detect classifier from <code>~/r16s-app/backend/classifiers/</code> based on marker type.
                                Place your <code>.qza</code> file there to use it automatically.
                              </div>
                              {selectedRegion === "V7-V8" && (
                                <div className="ps-info-box ps-info-box--warn" style={{ fontSize: 12, marginTop: 8 }}>
                                  ⚠️ <strong>V7–V8 region ต้องการ classifier เฉพาะ region</strong> — classifier มาตรฐาน (V3–V4/V4)
                                  จะให้ผลไม่ดีกับ amplicon ที่ตัดจาก 1055F/1392R
                                  แนะนำเปลี่ยนเป็น <strong>🟣 Train from primers</strong> เพื่อ train V7–V8 classifier จาก SILVA อัตโนมัติ
                                  (ครั้งแรกใช้เวลา 10–30 นาที, ครั้งต่อไป reuse ได้เลย)
                                </div>
                              )}
                            </div>
                          )}

                          {params.customClassifierMode === "train" && (
                            <div>
                              <div className="ps-info-box ps-info-box--warn" style={{ fontSize: 12, marginBottom: 10 }}>
                                ⚙️ Will run <strong>qiime feature-classifier extract-reads</strong> + <strong>fit-classifier-naive-bayes</strong>
                                on first job. Training may take 10–30 min and requires SILVA database.
                              </div>
                              <div className="param-grid">
                                <div className="param-item">
                                  <label className="param-label">Forward Primer</label>
                                  <span className="param-hint">Same as Step 1 — uses primer_f automatically</span>
                                  <input type="text" className="param-input"
                                    style={{ fontFamily: "monospace", background: "#0f172a", color: "#64748b" }}
                                    value={params.primer_f || "(from Step 1)"} readOnly />
                                </div>
                                <div className="param-item">
                                  <label className="param-label">Reverse Primer</label>
                                  <span className="param-hint">Same as Step 1 — uses primer_r automatically</span>
                                  <input type="text" className="param-input"
                                    style={{ fontFamily: "monospace", background: "#0f172a", color: "#64748b" }}
                                    value={params.primer_r || "(from Step 1)"} readOnly />
                                </div>
                                <div className="param-item">
                                  <label className="param-label">Min Amplicon Length (bp)</label>
                                  <span className="param-hint">Trim extracted reads shorter than this</span>
                                  <input type="number" className="param-input"
                                    min={50} max={2000} value={params.trainAmpliconMinLen}
                                    onChange={e => set("trainAmpliconMinLen", Number(e.target.value))} />
                                </div>
                                <div className="param-item">
                                  <label className="param-label">Max Amplicon Length (bp)</label>
                                  <span className="param-hint">Trim extracted reads longer than this</span>
                                  <input type="number" className="param-input"
                                    min={50} max={2000} value={params.trainAmpliconMaxLen}
                                    onChange={e => set("trainAmpliconMaxLen", Number(e.target.value))} />
                                </div>
                              </div>
                              <div style={{ fontSize: 11, color: "#475569", marginTop: 8 }}>
                                💡 Tip: V3–V4 (341F/806R) → min 350 max 500 &nbsp;|&nbsp;
                                V4 (515F/806R) → min 250 max 300 &nbsp;|&nbsp;
                                V7–V8 (1055F/1392R) → min 300 max 400 &nbsp;|&nbsp;
                                V4–V5 (515F/926R) → min 300 max 450 &nbsp;|&nbsp;
                                V1–V3 (27F/534R) → min 400 max 550
                                {selectedRegion === "V7-V8" && (
                                  <span style={{ color: "#a78bfa", fontWeight: 600 }}>
                                    &nbsp;← ค่า V7–V8 ถูกตั้งให้แล้วอัตโนมัติ
                                  </span>
                                )}
                              </div>
                            </div>
                          )}

                          {params.customClassifierMode === "upload" && (
                            <div>
                              <label className="param-label">Path to pre-trained classifier (.qza)</label>
                              <div style={{ display: "flex", gap: 8, marginTop: 4 }}>
                                <input type="text" className="param-input" style={{ flex: 1 }}
                                  placeholder="/path/to/classifier.qza"
                                  value={params.customClassifierPath}
                                  onChange={e => set("customClassifierPath", e.target.value)} />
                                <button className="browse-btn" onClick={() => setShowBrowser(true)}>📂</button>
                              </div>
                              <div style={{ fontSize: 11, color: "#475569", marginTop: 6 }}>
                                Download pre-trained classifiers from{" "}
                                <span style={{ color: "#6366f1" }}>docs.qiime2.org → Data resources</span>
                              </div>
                            </div>
                          )}
                        </div>
                      )}

                      {/* Method */}
                      <div className="ps-section-label" style={{ marginTop: 16 }}>METHOD</div>
                      <div className="param-grid">
                        <ParamSelect label="Classification Method" hint="Algorithm for assigning taxonomy"
                          value={params.taxMethod}
                          options={[
                            { value: "sklearn",        label: "sklearn Naive Bayes (fast, default)" },
                            { value: "blast_consensus", label: "BLAST consensus" },
                            { value: "vsearch",        label: "VSEARCH global alignment" },
                          ]}
                          onChange={v => set("taxMethod", v)} />

                        <div className="param-item">
                          <label className="param-label">Confidence Threshold</label>
                          <span className="param-hint">Min confidence to retain assignment (0.5–1.0)</span>
                          <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
                            <input type="range" min={0.5} max={1.0} step={0.05}
                              value={params.confidence}
                              onChange={e => set("confidence", Number(e.target.value))}
                              style={{ flex: 1, accentColor: "#22d3ee" }} />
                            <span style={{ color: "#22d3ee", fontWeight: 700, width: 32, textAlign: "right" }}>
                              {params.confidence.toFixed(2)}
                            </span>
                          </div>
                        </div>

                        <ParamNumber label="Min Bootstrap (%)" hint="Minimum confidence for assignment"
                          value={params.minBoot} min={0} max={100} step={5}
                          onChange={v => set("minBoot", v)} />

                        <div className="param-item">
                          <label className="param-label">Top Taxa to Show in Plots</label>
                          <div className="topn-toggle">
                            {([30,50,100] as const).map(n => (
                              <button key={n} className={`topn-btn ${params.topN === n ? "active" : ""}`}
                                onClick={() => set("topN", n)}>Top {n}</button>
                            ))}
                          </div>
                        </div>
                      </div>

                      {/* Collapse levels */}
                      <div className="ps-section-label" style={{ marginTop: 16 }}>COLLAPSE LEVELS</div>
                      <p className="param-hint" style={{ marginBottom: 8 }}>
                        Taxonomic levels to collapse for bar plots and diversity. L2=Phylum, L3=Class, L4=Order, L5=Family, L6=Genus, L7=Species
                      </p>
                      <div className="ps-level-chips">
                        {[2,3,4,5,6,7].map(l => (
                          <button key={l}
                            className={`ps-level-chip ${params.collapseLevels.includes(l) ? "ps-level-chip--active" : ""}`}
                            onClick={() => set("collapseLevels", toggle(params.collapseLevels, l))}>
                            L{l}
                          </button>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* ── Step 6: Diversity ─────────────────────────────── */}
                  {s.id === 6 && (
                    <div>
                      <div className="ps-toggle-row">
                        <span className="param-label">Run Diversity Analysis</span>
                        <div className="topn-toggle">
                          <button className={`topn-btn ${params.runDiversity ? "active" : ""}`}
                            onClick={() => set("runDiversity", true)}>Enable</button>
                          <button className={`topn-btn ${!params.runDiversity ? "active" : ""}`}
                            onClick={() => set("runDiversity", false)}>Skip</button>
                        </div>
                      </div>

                      {params.runDiversity && (
                        <>
                          <div className="ps-section-label" style={{ marginTop: 16 }}>CORE SETTINGS</div>
                          <div className="param-grid">
                            <ParamNumber label="Sampling Depth (Rarefaction)"
                              hint="Samples below this depth will be excluded"
                              value={params.samplingDepth} min={100} max={200000} step={1000}
                              onChange={v => set("samplingDepth", v)} />
                            <div className="param-item">
                              <label className="param-label">Primary Group Column</label>
                              <span className="param-hint">Main metadata column for group comparisons</span>
                              <input type="text" className="param-input"
                                value={params.groupCol} placeholder="e.g. treatment"
                                onChange={e => set("groupCol", e.target.value)} />
                            </div>
                            <div className="param-item" style={{ gridColumn: "1/-1" }}>
                              <label className="param-label">Additional Metadata Columns <span style={{ fontWeight: 400, color: "#6b7280" }}>(comma-separated)</span></label>
                              <input type="text" className="param-input"
                                value={params.additionalGroupCols}
                                placeholder="sex,age_group,site"
                                onChange={e => set("additionalGroupCols", e.target.value)} />
                            </div>
                          </div>

                          {/* Phylogenetic metrics */}
                          {hasPhylo && params.runPhylogeny && (
                            <>
                              <div className="ps-section-label" style={{ marginTop: 16 }}>
                                PHYLOGENETIC METRICS <span style={{ color: "#6b7280", fontSize: 11 }}>({marker} only)</span>
                              </div>
                              <div className="ps-info-box ps-info-box--ok" style={{ marginBottom: 8 }}>
                                ✅ rooted-tree available — core-metrics-phylogenetic will be used
                              </div>
                              <div className="ps-metric-chips">
                                {[
                                  { key: "weighted_unifrac",   label: "Weighted UniFrac" },
                                  { key: "unweighted_unifrac", label: "Unweighted UniFrac" },
                                  { key: "faith_pd",           label: "Faith PD" },
                                ].map(m => (
                                  <button key={m.key}
                                    className={`ps-metric-chip ${params.diversityPhyloMetrics.includes(m.key) ? "ps-metric-chip--active" : ""}`}
                                    onClick={() => set("diversityPhyloMetrics", toggle(params.diversityPhyloMetrics, m.key))}>
                                    {m.label}
                                  </button>
                                ))}
                              </div>
                            </>
                          )}

                          {/* Non-phylogenetic */}
                          <div className="ps-section-label" style={{ marginTop: 16 }}>NON-PHYLOGENETIC METRICS</div>
                          <div className="ps-metric-chips">
                            {[
                              { key: "shannon",           label: "Shannon" },
                              { key: "observed_features", label: "Observed Features" },
                              { key: "evenness",          label: "Evenness" },
                              { key: "bray_curtis",       label: "Bray-Curtis" },
                              { key: "jaccard",           label: "Jaccard" },
                            ].map(m => (
                              <button key={m.key}
                                className={`ps-metric-chip ${params.diversityNonPhyloMetrics.includes(m.key) ? "ps-metric-chip--active" : ""}`}
                                onClick={() => set("diversityNonPhyloMetrics", toggle(params.diversityNonPhyloMetrics, m.key))}>
                                {m.label}
                              </button>
                            ))}
                          </div>

                          {/* Outputs */}
                          <div className="ps-section-label" style={{ marginTop: 16 }}>OUTPUTS GENERATED</div>
                          <div className="ps-cmd-block">
                            <pre className="ps-cmd-body">{`diversity/
├── alpha/   shannon.tsv  observed_features.tsv  evenness.tsv  faith_pd.tsv
├── beta/    bray_curtis_pcoa.pdf  weighted_unifrac_pcoa.pdf
│            permanova_results.tsv
└── rarefaction_curves.pdf`}</pre>
                          </div>
                        </>
                      )}
                    </div>
                  )}

                  {/* ── Step 7: Differential Abundance ───────────────── */}
                  {s.id === 7 && (
                    <div>
                      <div className="ps-toggle-row">
                        <span className="param-label">Run ANCOM-BC (Differential Abundance)</span>
                        <div className="topn-toggle">
                          <button className={`topn-btn ${params.runDiffAbund ? "active" : ""}`}
                            onClick={() => set("runDiffAbund", true)}>Enable</button>
                          <button className={`topn-btn ${!params.runDiffAbund ? "active" : ""}`}
                            onClick={() => set("runDiffAbund", false)}>Skip</button>
                        </div>
                      </div>

                      {params.runDiffAbund && (
                        <>
                          {!params.metadataPath && (
                            <div className="ps-info-box ps-info-box--warn" style={{ marginTop: 12 }}>
                              ⚠️ Metadata required for ANCOM-BC — upload metadata in the Metadata section above.
                            </div>
                          )}
                          <div className="ps-section-label" style={{ marginTop: 16 }}>ANCOM-BC SETTINGS</div>
                          <div className="param-grid">
                            <div className="param-item">
                              <label className="param-label">Formula (metadata columns)</label>
                              <span className="param-hint">Defaults to primary group column</span>
                              <input type="text" className="param-input"
                                value={params.groupCol} placeholder="e.g. treatment"
                                onChange={e => set("groupCol", e.target.value)} />
                            </div>
                            <div className="param-item">
                              <label className="param-label">Taxonomic Level</label>
                              <select className="param-input" value={params.diffAbundLevel}
                                onChange={e => set("diffAbundLevel", e.target.value)}>
                                <option value="L2">Phylum (L2)</option>
                                <option value="L3">Class (L3)</option>
                                <option value="L4">Order (L4)</option>
                                <option value="L5">Family (L5)</option>
                                <option value="L6">Genus (L6)</option>
                                <option value="L7">Species (L7)</option>
                              </select>
                            </div>
                            <ParamNumber label="Min Frequency Filter"
                              hint="Remove features with fewer than this many total counts"
                              value={params.diffAbundMinFreq} min={1} max={1000}
                              onChange={v => set("diffAbundMinFreq", v)} />
                            <ParamNumber label="P-value Threshold"
                              hint="Adjusted p-value cutoff for significance"
                              value={params.diffAbundPval} min={0.001} max={0.1} step={0.005}
                              onChange={v => set("diffAbundPval", v)} />
                          </div>

                          <div className="ps-section-label" style={{ marginTop: 16 }}>OUTPUTS GENERATED</div>
                          <div className="ps-cmd-block">
                            <pre className="ps-cmd-body">{`diffabund/
├── ancombc_results.tsv       (log-fold-change + p-values)
├── ancombc_significant.tsv   (significant taxa only, p < ${params.diffAbundPval})
└── ancombc_barplot.pdf       (lollipop / waterfall plot)`}</pre>
                          </div>
                        </>
                      )}
                    </div>
                  )}

                  {/* ── Step 8: Functional Analysis ───────────────────── */}
                  {s.id === 8 && (
                    <div>
                      <div className="ps-func-grid">

                        {/* PICRUSt2 */}
                        <div className="ps-func-card">
                          <div className="ps-func-title">🔬 PICRUSt2</div>
                          <div className="ps-func-desc">
                            Phylogenetic placement + KEGG/MetaCyc pathway prediction (16S only, requires conda env)
                          </div>
                          <div className="topn-toggle" style={{ marginTop: 8 }}>
                            <button className={`topn-btn ${params.run_picrust2 ? "active" : ""}`}
                              onClick={() => set("run_picrust2", true)}>Enable</button>
                            <button className={`topn-btn ${!params.run_picrust2 ? "active" : ""}`}
                              onClick={() => set("run_picrust2", false)}>Skip</button>
                          </div>

                          {params.run_picrust2 && (
                            <div style={{ marginTop: 12 }}>
                              <div className="ps-section-label">PICRUST2 SETTINGS</div>
                              <div className="param-grid" style={{ gridTemplateColumns: "1fr 1fr" }}>
                                <div className="param-item">
                                  <label className="param-label">Placement Tool</label>
                                  <select className="param-input" value={params.picrust2PlaceTool}
                                    onChange={e => set("picrust2PlaceTool", e.target.value)}>
                                    <option value="epa-ng">epa-ng (fast, default)</option>
                                    <option value="sepp">SEPP (accurate)</option>
                                  </select>
                                </div>
                                <ParamNumber label="Max NSTI" hint="Max Nearest Sequenced Taxon Index"
                                  value={params.picrust2MaxNSTI} min={0.5} max={5} step={0.5}
                                  onChange={v => set("picrust2MaxNSTI", v)} />
                                <ParamNumber label="Threads" hint="Parallel workers for PICRUSt2"
                                  value={params.picrust2Threads} min={1} max={32}
                                  onChange={v => set("picrust2Threads", v)} />
                              </div>
                              <div className="ps-section-label" style={{ marginTop: 12 }}>OUTPUT DATABASES</div>
                              <div className="ps-metric-chips">
                                {[
                                  { key: "metacyc", label: "MetaCyc pathways" },
                                  { key: "ec",      label: "EC numbers" },
                                  { key: "ko",      label: "KO (KEGG)" },
                                ].map(db => (
                                  <button key={db.key}
                                    className={`ps-metric-chip ${params.picrust2Databases.includes(db.key) ? "ps-metric-chip--active" : ""}`}
                                    onClick={() => set("picrust2Databases", toggle(params.picrust2Databases, db.key))}>
                                    {db.label}
                                  </button>
                                ))}
                              </div>
                              <div className="ps-section-label" style={{ marginTop: 12 }}>OUTPUTS GENERATED</div>
                              <div className="ps-cmd-block">
                                <pre className="ps-cmd-body">{`picrust2/
├── pathway_abundance.tsv
├── ec_metagenome.tsv
└── ko_metagenome.tsv`}</pre>
                              </div>
                            </div>
                          )}
                        </div>

                        {/* Tax4Fun2 */}
                        <div className="ps-func-card">
                          <div className="ps-func-title">📊 Tax4Fun2</div>
                          <div className="ps-func-desc">
                            KEGG KO prediction from 16S (R package; fast, no extra conda needed)
                          </div>
                          <div className="topn-toggle" style={{ marginTop: 8 }}>
                            <button className={`topn-btn ${params.run_tax4fun ? "active" : ""}`}
                              onClick={() => set("run_tax4fun", true)}>Enable</button>
                            <button className={`topn-btn ${!params.run_tax4fun ? "active" : ""}`}
                              onClick={() => set("run_tax4fun", false)}>Skip</button>
                          </div>
                          {params.run_tax4fun && (
                            <div className="ps-info-box ps-info-box--ok" style={{ marginTop: 10 }}>
                              ✅ Tax4Fun2 will run using the SILVA reference database
                            </div>
                          )}
                        </div>

                        {/* FUNGuildR (ITS only) */}
                        {isITS && (
                          <div className="ps-func-card">
                            <div className="ps-func-title">🍄 FUNGuildR</div>
                            <div className="ps-func-desc">
                              Ecological guild annotation for fungi (ITS only)
                            </div>
                            <div className="topn-toggle" style={{ marginTop: 8 }}>
                              <button className="topn-btn active">Auto-enabled for ITS</button>
                            </div>
                          </div>
                        )}
                      </div>
                    </div>
                  )}

                </div>
              )}
            </div>
          ))}
        </>
      )}

      {/* ── Metadata upload (all pipeline types) ──────────────────────── */}
      <div className="ext-section">
        <h4 className="ext-section-title">
          📊 Metadata{" "}
          <span style={{ fontWeight: 400, fontSize: 13, color: "#6b7280" }}>
            (optional — for diversity grouping &amp; ANCOM-BC)
          </span>
        </h4>
        <MetadataUpload
          inline
          onMetadataReady={(info: MetadataInfo, grp: string) => {
            onChange({ ...params, metadataPath: info.path, groupCol: grp });
            onMetadata?.(info.path, grp);
          }}
          onGroupColChange={(col: string) => onChange({ ...params, groupCol: col })}
        />
        {params.metadataPath && (
          <p style={{ fontSize: 12, color: "#10b981", marginTop: 8 }}>
            ✅ Metadata loaded · Group column: <strong>{params.groupCol}</strong>
          </p>
        )}
      </div>

    </div>
  );
}

// ── Sub-components ────────────────────────────────────────────────────────────
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
      <select className="param-input" value={value} onChange={e => onChange(e.target.value)}>
        {options.map(o => <option key={o.value} value={o.value}>{o.label}</option>)}
      </select>
    </div>
  );
}
