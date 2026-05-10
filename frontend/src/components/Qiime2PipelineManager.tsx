/**
 * Qiime2PipelineManager
 * ─────────────────────
 * Full QIIME2 pipeline UI with sidebar step navigation:
 *   SETUP  → Environment Check / Project Setup / Input Validation
 *   PIPELINE → 9 steps (Import → Cutadapt → Demux QC → DADA2 →
 *               Phylogeny → Taxonomy → Diversity → Diff. Abundance → Functional)
 *   OUTPUT → Export Results / Output Files
 */

import { useState, useEffect, useRef, useCallback } from "react";
import MetadataUpload, { MetadataInfo } from "./MetadataUpload";

const API = "";

// ── Types ─────────────────────────────────────────────────────────────────────
type StepId =
  | "env-check" | "project-setup" | "input-validation"
  | "import" | "cutadapt" | "demux-qc" | "dada2"
  | "phylogeny" | "taxonomy" | "diversity" | "diff-abundance" | "functional"
  | "export" | "output-files";

type StepStatus = "pending" | "active" | "done" | "error" | "skip";

export type MarkerType = "16S" | "12S" | "ITS1" | "ITS2" | "COX1" | "18S-nema" | "PacBio";

interface SidebarItem {
  id:    StepId;
  label: string;
  tag?:  string;
  icon?: string;
}

const SETUP_STEPS: SidebarItem[] = [
  { id: "env-check",          label: "Environment Check", icon: "🔍" },
  { id: "project-setup",      label: "Project Setup",     icon: "⚙️" },
  { id: "input-validation",   label: "Input Validation",  icon: "✅" },
];

const PIPELINE_STEPS: SidebarItem[] = [
  { id: "import",        label: "1. Import Data",       icon: "📥" },
  { id: "cutadapt",      label: "2. Cutadapt",          icon: "✂️",  tag: "primers" },
  { id: "demux-qc",      label: "3. Demux QC",          icon: "📊" },
  { id: "dada2",         label: "4. DADA2",             icon: "🔬",  tag: "ASV" },
  { id: "phylogeny",     label: "5. Phylogeny",         icon: "🌳",  tag: "16S" },
  { id: "taxonomy",      label: "6. Taxonomy",          icon: "🏷️" },
  { id: "diversity",     label: "7. Diversity",         icon: "🌐" },
  { id: "diff-abundance",label: "8. Diff. Abundance",  icon: "📈" },
  { id: "functional",    label: "9. Functional",        icon: "🔗",  tag: "PICRUSt2" },
];

const OUTPUT_STEPS: SidebarItem[] = [
  { id: "export",       label: "Export Results", icon: "📦" },
  { id: "output-files", label: "Output Files",   icon: "📂" },
];

const ALL_STEPS: StepId[] = [
  ...SETUP_STEPS.map(s => s.id),
  ...PIPELINE_STEPS.map(s => s.id),
  ...OUTPUT_STEPS.map(s => s.id),
];

const MARKER_OPTIONS: { value: MarkerType; label: string; icon: string; hint: string }[] = [
  { value: "16S",      icon: "🦠", label: "16S rRNA",      hint: "Bacteria & Archaea — SILVA 138" },
  { value: "12S",      icon: "🐟", label: "12S rRNA",      hint: "Vertebrates / Fish eDNA" },
  { value: "ITS1",     icon: "🍄", label: "ITS1 Fungi",    hint: "Fungal ITS1 — UNITE v10" },
  { value: "ITS2",     icon: "🍄", label: "ITS2 Fungi",    hint: "Fungal ITS2 — UNITE v10" },
  { value: "COX1",     icon: "🦑", label: "COX1 / CO1",    hint: "Animal metabarcoding — MIDORI2" },
  { value: "18S-nema", icon: "🐛", label: "18S Nematode",  hint: "Nematode 18S — NemaBase / PR2" },
  { value: "PacBio",   icon: "🧬", label: "PacBio CCS 16S",hint: "Full-length V1–V9 long reads" },
];

interface Params {
  jobName:        string;
  marker:         MarkerType;
  // primers
  primerF:        string;
  primerR:        string;
  // DADA2
  truncLenF:      number;
  truncLenR:      number;
  maxEEF:         number;
  maxEER:         number;
  trimLeftF:      number;
  trimLeftR:      number;
  chimeraMethod:  "consensus" | "per-sample" | "pooled";
  nThreads:       number;
  // taxonomy
  classifierPath: string;
  confidence:     number;
  // diversity
  samplingDepth:  number;
  groupCol:       string;
  // diff abundance
  runDiffAbund:   boolean;
  diffFormula:    string;
  // functional
  runPicrust2:    boolean;
  runRViz:        boolean;
  // metadata
  metadataPath:   string;
  // manifest
  manifestPath:   string;
}

const DEFAULT_PRIMERS: Record<string, { f: string; r: string }> = {
  "16S":      { f: "GTGYCAGCMGCCGCGGTAA", r: "GGACTACNVGGGTWTCTAAT" },
  "12S":      { f: "ACTGGGATTAGATACCCC",   r: "TAGAACAGGCTCCTCTAG" },
  "ITS1":     { f: "CTTGGTCATTTAGAGGAAGTAA", r: "GCTGCGTTCTTCATCGATGC" },
  "ITS2":     { f: "GTGAATCATCGAATCTTTGAA", r: "TCCTCCGCTTATTGATATGC" },
  "COX1":     { f: "GGWACWGGWTGAACWGTWTAYCCYCC", r: "TANACYTCNGGRTGNCCRAARAAYCA" },
  "18S-nema": { f: "CGCGAATRGCTCATTACAACAGC", r: "GGGCGGTGTGTACAAAGGGCAGGG" },
  "PacBio":   { f: "AGRGTTYGATYMTGGCTCAG", r: "RGYTACCTTGTTACGACTT" },
};

// ── Component ─────────────────────────────────────────────────────────────────
interface Props {
  onBack?: () => void;
}

export default function Qiime2PipelineManager({ onBack }: Props) {
  const [activeStep, setActiveStep]       = useState<StepId>("env-check");
  const [stepStatus, setStepStatus]       = useState<Record<StepId, StepStatus>>(
    Object.fromEntries(ALL_STEPS.map(id => [id, "pending"])) as Record<StepId, StepStatus>
  );

  const [params, setParams] = useState<Params>({
    jobName: "", marker: "16S",
    primerF: DEFAULT_PRIMERS["16S"].f, primerR: DEFAULT_PRIMERS["16S"].r,
    truncLenF: 240, truncLenR: 200, maxEEF: 2, maxEER: 2,
    trimLeftF: 0, trimLeftR: 0,
    chimeraMethod: "consensus", nThreads: 4,
    classifierPath: "", confidence: 0.7,
    samplingDepth: 10000, groupCol: "treatment",
    runDiffAbund: true, diffFormula: "",
    runPicrust2: false, runRViz: true,
    metadataPath: "", manifestPath: "",
  });

  const [envInfo,      setEnvInfo]      = useState<any>(null);
  const [envLoading,   setEnvLoading]   = useState(false);
  const [manifestInfo, setManifestInfo] = useState<any>(null);
  const [metaInfo,     setMetaInfo]     = useState<MetadataInfo | null>(null);

  // Job tracking
  const [jobId,      setJobId]      = useState<string | null>(null);
  const [jobStatus,  setJobStatus]  = useState<string>("idle");
  const [logs,       setLogs]       = useState<string[]>([]);
  const [progress,   setProgress]   = useState(0);
  const [progLabel,  setProgLabel]  = useState("");
  const [running,    setRunning]    = useState(false);
  const [outputDir,  setOutputDir]  = useState<string | null>(null);
  const [outputFiles, setOutputFiles] = useState<string[]>([]);

  const logRef    = useRef<HTMLDivElement>(null);
  const pollRef   = useRef<any>(null);
  const abortRef  = useRef<AbortController | null>(null);

  // Active pipeline step count (for "Pipeline N/9" counter)
  const pipelineStepStatuses = PIPELINE_STEPS.map(s => stepStatus[s.id]);
  const doneCount = pipelineStepStatuses.filter(s => s === "done").length;

  // ── Sync step status → active when visiting ──────────────────────────────
  useEffect(() => {
    setStepStatus(prev => {
      const next = { ...prev };
      if (next[activeStep] === "pending") next[activeStep] = "active";
      return next;
    });
  }, [activeStep]);

  // ── Env check ─────────────────────────────────────────────────────────────
  const checkEnv = useCallback(async () => {
    setEnvLoading(true);
    try {
      const res = await fetch(`${API}/qiime2/env/check`).then(r => r.json());
      setEnvInfo(res);
      if (res.found) {
        setStepStatus(prev => ({ ...prev, "env-check": "done" }));
      } else {
        setStepStatus(prev => ({ ...prev, "env-check": "error" }));
      }
    } catch {
      setEnvInfo({ found: false, errors: ["Cannot reach backend"] });
      setStepStatus(prev => ({ ...prev, "env-check": "error" }));
    }
    setEnvLoading(false);
  }, []);

  useEffect(() => { checkEnv(); }, []);

  // ── Manifest upload ────────────────────────────────────────────────────────
  const uploadManifest = async (file: File) => {
    const fd = new FormData();
    fd.append("file", file);
    try {
      const res = await fetch(`${API}/qiime2/upload/manifest`, { method: "POST", body: fd }).then(r => r.json());
      setManifestInfo(res);
      if (res.valid) {
        setParams(p => ({ ...p, manifestPath: res.path }));
        setStepStatus(prev => ({ ...prev, "input-validation": "done", "import": "pending" }));
      }
    } catch { alert("Manifest upload failed"); }
  };

  // ── Change marker — reset primers ─────────────────────────────────────────
  const setMarker = (m: MarkerType) => {
    const pr = DEFAULT_PRIMERS[m] || { f: "", r: "" };
    setParams(p => ({ ...p, marker: m, primerF: pr.f, primerR: pr.r }));
    setStepStatus(prev => ({
      ...prev,
      "phylogeny": m === "ITS1" || m === "ITS2" || m === "COX1" || m === "18S-nema" ? "skip" : "pending",
    }));
  };

  // ── Run pipeline ──────────────────────────────────────────────────────────
  const runPipeline = async () => {
    if (!params.manifestPath) { alert("Please upload a manifest file first (Input Validation step)"); return; }
    setRunning(true);
    setLogs([]);
    setProgress(0);
    setJobStatus("running");
    setOutputDir(null);

    // Reset pipeline step statuses
    setStepStatus(prev => {
      const next = { ...prev };
      PIPELINE_STEPS.forEach(s => { next[s.id] = "pending"; });
      next["export"] = "pending"; next["output-files"] = "pending";
      return next;
    });

    try {
      const res = await fetch(`${API}/qiime2/run`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          marker:            params.marker,
          manifest_path:     params.manifestPath,
          metadata_path:     params.metadataPath || null,
          trim_left_f:       params.trimLeftF,
          trim_left_r:       params.trimLeftR,
          trunc_len_f:       params.truncLenF,
          trunc_len_r:       params.truncLenR,
          max_ee_f:          params.maxEEF,
          max_ee_r:          params.maxEER,
          chimera_method:    params.chimeraMethod,
          n_threads:         params.nThreads,
          primer_f:          params.primerF,
          primer_r:          params.primerR,
          classifier_path:   params.classifierPath,
          confidence:        params.confidence,
          sampling_depth:    params.samplingDepth,
          group_col:         params.groupCol,
          run_diffabund:     params.runDiffAbund,
          diffabund_formula: params.diffFormula,
          run_r_viz:         params.runRViz,
        }),
      }).then(r => r.json());

      const jid = res.job_id;
      setJobId(jid);
      setActiveStep("import");

      // Poll status every 3 s
      pollRef.current = setInterval(async () => {
        try {
          const st = await fetch(`${API}/status/${jid}`).then(r => r.json());
          setProgress(st.progress ?? 0);
          setProgLabel(st.step_label ?? "");
          setJobStatus(st.status);

          // Map step_label → sidebar step
          const label = (st.step_label || "").toLowerCase();
          if      (label.includes("import"))        markStep("import");
          else if (label.includes("cutadapt"))      markStep("cutadapt");
          else if (label.includes("demux"))         markStep("demux-qc");
          else if (label.includes("dada2"))         markStep("dada2");
          else if (label.includes("phylogeny"))     markStep("phylogeny");
          else if (label.includes("taxonomy"))      markStep("taxonomy");
          else if (label.includes("diversity"))     markStep("diversity");
          else if (label.includes("ancom"))         markStep("diff-abundance");
          else if (label.includes("r visual") || label.includes("picrust")) markStep("functional");
          else if (label.includes("export"))        markStep("export");

          if (st.logs) setLogs(st.logs);
          if (st.step_log) setLogs(prev => [...prev, ...st.step_log]);

          if (st.status === "done" || st.status === "completed") {
            clearInterval(pollRef.current);
            setRunning(false);
            setOutputDir(st.output_dir ?? null);
            setStepStatus(prev => ({ ...prev, "export": "done", "output-files": "done" }));
            setActiveStep("output-files");
            if (st.output_dir) fetchOutputFiles(st.output_dir);
          }
          if (st.status === "error") {
            clearInterval(pollRef.current);
            setRunning(false);
          }
        } catch { /* ignore */ }
      }, 3000);

    } catch (e: any) {
      alert("Failed to start pipeline: " + e.message);
      setRunning(false);
      setJobStatus("error");
    }
  };

  const markStep = (id: StepId) => {
    setActiveStep(id);
    setStepStatus(prev => {
      const next = { ...prev };
      // Mark prior pipeline steps as done
      const idx = PIPELINE_STEPS.findIndex(s => s.id === id);
      PIPELINE_STEPS.slice(0, idx).forEach(s => {
        if (next[s.id] !== "skip") next[s.id] = "done";
      });
      next[id] = "active";
      return next;
    });
  };

  const fetchOutputFiles = async (dir: string) => {
    try {
      const res = await fetch(`${API}/results/${jobId}`).then(r => r.json());
      setOutputFiles(res.files ?? []);
    } catch { /* noop */ }
  };

  // Auto-scroll logs
  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight;
  }, [logs]);

  // Cleanup
  useEffect(() => () => { if (pollRef.current) clearInterval(pollRef.current); }, []);

  const set = (key: keyof Params, val: any) => setParams(p => ({ ...p, [key]: val }));

  // ── Sidebar item renderer ─────────────────────────────────────────────────
  const SidebarItem = ({ item }: { item: SidebarItem }) => {
    const st = stepStatus[item.id];
    const isActive = activeStep === item.id;
    return (
      <button
        className={`q2-sidebar-item ${isActive ? "q2-sidebar-item--active" : ""} ${st === "done" ? "q2-sidebar-item--done" : ""} ${st === "error" ? "q2-sidebar-item--error" : ""} ${st === "skip" ? "q2-sidebar-item--skip" : ""}`}
        onClick={() => setActiveStep(item.id)}
      >
        <span className="q2-sidebar-dot">
          {st === "done"   ? "✓" :
           st === "error"  ? "✗" :
           st === "skip"   ? "–" :
           st === "active" ? "●" : "○"}
        </span>
        <span className="q2-sidebar-label">{item.label}</span>
        {item.tag && <span className="q2-sidebar-tag">{item.tag}</span>}
      </button>
    );
  };

  // ── Step content panels ───────────────────────────────────────────────────
  const renderStepContent = () => {
    switch (activeStep) {

      // ── SETUP ──────────────────────────────────────────────────────────
      case "env-check":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🔍 Environment Check</h2>
            <p className="q2-step-desc">Verify QIIME2 conda environment and pre-trained classifiers</p>
            {envLoading && <div className="q2-loading">Checking environment…</div>}
            {envInfo && (
              <div className={`q2-env-box ${envInfo.found ? "q2-env-box--ok" : "q2-env-box--err"}`}>
                <div className="q2-env-row">
                  <strong>{envInfo.found ? "✅ QIIME2 Found" : "❌ QIIME2 Not Found"}</strong>
                  {envInfo.version && <span className="q2-env-ver">{envInfo.version}</span>}
                </div>
                {envInfo.errors?.length > 0 && (
                  <div className="q2-env-errors">
                    {envInfo.errors.map((e: string, i: number) => <div key={i}>{e}</div>)}
                    <div className="q2-env-fix">Run: <code>bash ~/r16s-app/setup_qiime2.sh</code></div>
                  </div>
                )}
                {envInfo.classifiers && (
                  <div className="q2-classifiers">
                    <div className="q2-clf-title">Classifiers</div>
                    {Object.entries(envInfo.classifiers).map(([marker, info]: [string, any]) => (
                      <div key={marker} className={`q2-clf-row ${info.found ? "q2-clf-row--ok" : "q2-clf-row--miss"}`}>
                        <span className="q2-clf-marker">{marker}</span>
                        <span>{info.found ? `✓ ${info.path.split("/").pop()}` : "✗ Not found"}</span>
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}
            <button className="q2-btn q2-btn--secondary" onClick={checkEnv} disabled={envLoading}>
              ↻ Re-check
            </button>
            {envInfo?.found && (
              <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("project-setup")}>
                Next →
              </button>
            )}
          </div>
        );

      case "project-setup":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">⚙️ Project Setup</h2>
            <p className="q2-step-desc">Configure job name, marker type, and metadata</p>

            <div className="q2-field">
              <label className="q2-label">Job Name</label>
              <input className="q2-input" placeholder="e.g. Soil_16S_Run1 (optional)"
                value={params.jobName} onChange={e => set("jobName", e.target.value)} />
            </div>

            <div className="q2-field">
              <label className="q2-label">Amplicon Marker</label>
              <div className="q2-marker-grid">
                {MARKER_OPTIONS.map(m => (
                  <button key={m.value}
                    className={`q2-marker-btn ${params.marker === m.value ? "q2-marker-btn--active" : ""}`}
                    onClick={() => setMarker(m.value)}
                    title={m.hint}>
                    <span className="q2-marker-icon">{m.icon}</span>
                    <span className="q2-marker-label">{m.label}</span>
                  </button>
                ))}
              </div>
              <p className="q2-hint">{MARKER_OPTIONS.find(m => m.value === params.marker)?.hint}</p>
            </div>

            <div className="q2-field">
              <label className="q2-label">Metadata <span className="q2-optional">(optional)</span></label>
              <MetadataUpload inline
                onMetadataReady={(info, grp) => {
                  setMetaInfo(info);
                  set("metadataPath", info.path);
                  set("groupCol", grp);
                }}
                onGroupColChange={col => set("groupCol", col)}
              />
            </div>

            <button className="q2-btn q2-btn--primary" onClick={() => {
              setStepStatus(prev => ({ ...prev, "project-setup": "done" }));
              setActiveStep("input-validation");
            }}>Next →</button>
          </div>
        );

      case "input-validation":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">✅ Input Validation</h2>
            <p className="q2-step-desc">
              Upload a QIIME2 manifest TSV. Format: <code>sample-id</code>, <code>forward-absolute-filepath</code>, <code>reverse-absolute-filepath</code>
            </p>
            <pre className="q2-code-block">{`sample-id\tforward-absolute-filepath\treverse-absolute-filepath
S1\t/home/boss/data/S1_R1.fastq.gz\t/home/boss/data/S1_R2.fastq.gz
S2\t/home/boss/data/S2_R1.fastq.gz\t/home/boss/data/S2_R2.fastq.gz`}</pre>

            <div className="q2-upload-zone"
              onClick={() => document.getElementById("manifest-file-input")?.click()}
              onDragOver={e => e.preventDefault()}
              onDrop={e => {
                e.preventDefault();
                const f = e.dataTransfer.files[0];
                if (f) uploadManifest(f);
              }}>
              <input id="manifest-file-input" type="file" accept=".tsv,.txt,.csv"
                style={{ display: "none" }}
                onChange={e => { const f = e.target.files?.[0]; if (f) uploadManifest(f); }} />
              📋 Drop manifest TSV here or click to browse
            </div>

            {manifestInfo && (
              <div className={`q2-manifest-info ${manifestInfo.valid ? "q2-manifest-info--ok" : "q2-manifest-info--err"}`}>
                {manifestInfo.valid
                  ? `✅ ${manifestInfo.filename} — ${manifestInfo.n_samples} samples loaded`
                  : `❌ ${manifestInfo.issues?.join(" | ")}`}
                {manifestInfo.valid && manifestInfo.samples?.length > 0 && (
                  <div className="q2-sample-list">
                    {manifestInfo.samples.slice(0, 5).map((s: any) => (
                      <div key={s.id} className="q2-sample-row">
                        <span className="q2-sample-id">{s.id}</span>
                        <span className="q2-sample-path">{s.forward?.split("/").pop()}</span>
                      </div>
                    ))}
                    {manifestInfo.n_samples > 5 && <div>… and {manifestInfo.n_samples - 5} more</div>}
                  </div>
                )}
              </div>
            )}

            {manifestInfo?.valid && (
              <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("import")}>
                Next →
              </button>
            )}
          </div>
        );

      // ── PIPELINE ───────────────────────────────────────────────────────
      case "import":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">📥 1. Import Data</h2>
            <p className="q2-step-desc">Import FASTQ files into QIIME2 artifact (.qza) using the manifest</p>
            <div className="q2-info-box">
              <div>Manifest: <code>{params.manifestPath || "not set"}</code></div>
              <div>Artifact type: <code>SampleData[{params.marker === "PacBio" ? "SequencesWithQuality" : "PairedEndSequencesWithQuality"}]</code></div>
            </div>
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("cutadapt")}>Next →</button>
          </div>
        );

      case "cutadapt":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">✂️ 2. Cutadapt — Primer Trimming</h2>
            <p className="q2-step-desc">Remove primers from reads. Leave blank to skip.</p>
            <div className="q2-param-grid">
              <div className="q2-field">
                <label className="q2-label">Forward Primer (5'→3')</label>
                <input className="q2-input q2-input--mono" value={params.primerF}
                  onChange={e => set("primerF", e.target.value.toUpperCase())} />
              </div>
              <div className="q2-field">
                <label className="q2-label">Reverse Primer (5'→3')</label>
                <input className="q2-input q2-input--mono" value={params.primerR}
                  onChange={e => set("primerR", e.target.value.toUpperCase())} />
              </div>
            </div>
            <div className="q2-preset-row">
              <span className="q2-label" style={{ alignSelf: "center" }}>Presets:</span>
              {Object.entries(DEFAULT_PRIMERS).filter(([m]) => m === params.marker).map(([m, pr]) => (
                <button key={m} className="q2-preset-btn"
                  onClick={() => { set("primerF", pr.f); set("primerR", pr.r); }}>
                  {m} default
                </button>
              ))}
              <button className="q2-preset-btn" onClick={() => { set("primerF", ""); set("primerR", ""); }}>
                Clear (skip)
              </button>
            </div>
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("demux-qc")}>Next →</button>
          </div>
        );

      case "demux-qc":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">📊 3. Demux QC — Quality Trimming</h2>
            <p className="q2-step-desc">Set truncation lengths and maximum expected errors per read</p>
            {(params.marker === "ITS1" || params.marker === "ITS2") && (
              <div className="q2-info-box q2-info-box--warn">
                ℹ️ ITS amplicons have variable length — truncation is set to 0 (disabled) automatically.
              </div>
            )}
            <div className="q2-param-grid">
              <ParamNum label="Truncate Forward (bp)" hint="0 = no truncation"
                value={params.marker === "ITS1" || params.marker === "ITS2" ? 0 : params.truncLenF}
                disabled={params.marker === "ITS1" || params.marker === "ITS2"}
                onChange={v => set("truncLenF", v)} min={0} max={350} />
              <ParamNum label="Truncate Reverse (bp)" hint="0 = no truncation"
                value={params.marker === "ITS1" || params.marker === "ITS2" ? 0 : params.truncLenR}
                disabled={params.marker === "ITS1" || params.marker === "ITS2"}
                onChange={v => set("truncLenR", v)} min={0} max={350} />
              <ParamNum label="Max EE Forward" hint="Maximum expected errors" value={params.maxEEF}
                onChange={v => set("maxEEF", v)} min={1} max={10} step={0.5} />
              <ParamNum label="Max EE Reverse" hint="Maximum expected errors" value={params.maxEER}
                onChange={v => set("maxEER", v)} min={1} max={10} step={0.5} />
              <ParamNum label="Trim Left Forward (bp)" hint="Trim from 5' end (primer bases)"
                value={params.trimLeftF} onChange={v => set("trimLeftF", v)} min={0} max={50} />
              <ParamNum label="Trim Left Reverse (bp)" hint="Trim from 5' end (primer bases)"
                value={params.trimLeftR} onChange={v => set("trimLeftR", v)} min={0} max={50} />
            </div>
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("dada2")}>Next →</button>
          </div>
        );

      case "dada2":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🔬 4. DADA2 — ASV Inference</h2>
            <p className="q2-step-desc">Denoise reads and infer exact Amplicon Sequence Variants</p>
            <div className="q2-param-grid">
              <div className="q2-field">
                <label className="q2-label">Chimera Method</label>
                <select className="q2-select" value={params.chimeraMethod}
                  onChange={e => set("chimeraMethod", e.target.value)}>
                  <option value="consensus">consensus (recommended)</option>
                  <option value="per-sample">per-sample</option>
                  <option value="pooled">pooled</option>
                </select>
              </div>
              <ParamNum label="CPU Threads" hint="Parallel threads for DADA2"
                value={params.nThreads} onChange={v => set("nThreads", v)} min={1} max={32} />
            </div>
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("phylogeny")}>Next →</button>
          </div>
        );

      case "phylogeny": {
        const hasPhylo = !["ITS1","ITS2","COX1","18S-nema"].includes(params.marker);
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🌳 5. Phylogeny — MAFFT + FastTree</h2>
            {hasPhylo ? (
              <>
                <p className="q2-step-desc">Build multiple sequence alignment and phylogenetic tree for UniFrac diversity</p>
                <div className="q2-info-box">Method: MAFFT alignment → FastTree unrooted tree → mid-point rooting</div>
                <ParamNum label="CPU Threads (phylogeny)" hint="Threads for MAFFT"
                  value={params.nThreads} onChange={v => set("nThreads", v)} min={1} max={32} />
              </>
            ) : (
              <div className="q2-info-box q2-info-box--warn">
                ⚠️ Phylogenetic tree is not generated for <strong>{params.marker}</strong> (variable region / non-16S marker).
                UniFrac diversity metrics will be skipped.
              </div>
            )}
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("taxonomy")}>Next →</button>
          </div>
        );
      }

      case "taxonomy":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🏷️ 6. Taxonomy — classify-sklearn</h2>
            <p className="q2-step-desc">Assign taxonomy using pre-trained Naive Bayes classifier</p>
            <div className="q2-field">
              <label className="q2-label">Classifier Path <span className="q2-optional">(auto-detected)</span></label>
              <input className="q2-input q2-input--mono" placeholder="Leave blank to auto-detect from classifiers/"
                value={params.classifierPath} onChange={e => set("classifierPath", e.target.value)} />
              {envInfo?.classifiers?.[params.marker]?.found && (
                <div className="q2-clf-found">
                  ✅ Auto-detected: {envInfo.classifiers[params.marker].path.split("/").pop()}
                </div>
              )}
            </div>
            <ParamNum label="Confidence Threshold" hint="Minimum confidence to assign (0.7 default)"
              value={params.confidence} onChange={v => set("confidence", v)} min={0.5} max={1} step={0.05} />
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("diversity")}>Next →</button>
          </div>
        );

      case "diversity":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🌐 7. Diversity Analysis</h2>
            <p className="q2-step-desc">Alpha diversity (Shannon, Chao1, Faith PD) and beta diversity (PCoA, PERMANOVA)</p>
            {!params.metadataPath && (
              <div className="q2-info-box q2-info-box--warn">
                ⚠️ No metadata loaded — diversity analysis will run but group comparisons and significance tests will be skipped.
                Add metadata in <button className="q2-link-btn" onClick={() => setActiveStep("project-setup")}>Project Setup</button>.
              </div>
            )}
            <ParamNum label="Rarefaction Depth" hint="Samples below this depth will be excluded"
              value={params.samplingDepth} onChange={v => set("samplingDepth", v)} min={100} max={100000} step={1000} />
            {params.metadataPath && (
              <div className="q2-field">
                <label className="q2-label">Group Column</label>
                <input className="q2-input" value={params.groupCol}
                  onChange={e => set("groupCol", e.target.value)} placeholder="e.g. treatment" />
              </div>
            )}
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("diff-abundance")}>Next →</button>
          </div>
        );

      case "diff-abundance":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">📈 8. Differential Abundance — ANCOM-BC</h2>
            <p className="q2-step-desc">Identify taxa significantly different between groups (requires metadata)</p>
            <div className="q2-toggle-row">
              <label className="q2-label">Run ANCOM-BC</label>
              <div className="q2-toggle-btns">
                <button className={`q2-toggle ${params.runDiffAbund ? "q2-toggle--on" : ""}`}
                  onClick={() => set("runDiffAbund", true)}>Enable</button>
                <button className={`q2-toggle ${!params.runDiffAbund ? "q2-toggle--on" : ""}`}
                  onClick={() => set("runDiffAbund", false)}>Skip</button>
              </div>
            </div>
            {params.runDiffAbund && (
              <>
                {!params.metadataPath && (
                  <div className="q2-info-box q2-info-box--warn">⚠️ Metadata required for ANCOM-BC</div>
                )}
                <div className="q2-field">
                  <label className="q2-label">Formula <span className="q2-optional">(defaults to group column)</span></label>
                  <input className="q2-input" value={params.diffFormula}
                    onChange={e => set("diffFormula", e.target.value)}
                    placeholder={params.groupCol || "treatment"} />
                </div>
              </>
            )}
            <button className="q2-btn q2-btn--primary" onClick={() => setActiveStep("functional")}>Next →</button>
          </div>
        );

      case "functional":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">🔗 9. Functional Prediction & R Visualization</h2>
            <p className="q2-step-desc">Optional: predict functional pathways and generate publication-quality plots</p>
            <div className="q2-func-grid">
              <div className="q2-func-card">
                <div className="q2-func-title">🔬 PICRUSt2</div>
                <div className="q2-func-desc">Phylogenetic placement + KEGG/MetaCyc pathway prediction (16S only, requires conda env)</div>
                <div className="q2-toggle-btns">
                  <button className={`q2-toggle ${params.runPicrust2 ? "q2-toggle--on" : ""}`}
                    onClick={() => set("runPicrust2", true)}>Enable</button>
                  <button className={`q2-toggle ${!params.runPicrust2 ? "q2-toggle--on" : ""}`}
                    onClick={() => set("runPicrust2", false)}>Skip</button>
                </div>
              </div>
              <div className="q2-func-card">
                <div className="q2-func-title">📊 R Visualization</div>
                <div className="q2-func-desc">phyloseq + ggplot2 + ANCOMBC2 + FUNGuildR — publication-ready PDFs</div>
                <div className="q2-toggle-btns">
                  <button className={`q2-toggle ${params.runRViz ? "q2-toggle--on" : ""}`}
                    onClick={() => set("runRViz", true)}>Enable</button>
                  <button className={`q2-toggle ${!params.runRViz ? "q2-toggle--on" : ""}`}
                    onClick={() => set("runRViz", false)}>Skip</button>
                </div>
              </div>
            </div>

            {/* ── RUN BUTTON ── */}
            <div className="q2-run-section">
              <div className="q2-run-summary">
                <div><strong>Marker:</strong> {params.marker}</div>
                <div><strong>Samples:</strong> {manifestInfo?.n_samples ?? "?"}</div>
                {params.metadataPath && <div><strong>Metadata:</strong> ✅ loaded ({params.groupCol})</div>}
                <div><strong>Threads:</strong> {params.nThreads}</div>
              </div>
              <button className="q2-run-btn" onClick={runPipeline} disabled={running || !params.manifestPath}>
                {running ? "⏳ Running…" : "🚀 Run Full Pipeline"}
              </button>
              {!params.manifestPath && (
                <p className="q2-run-warn">Upload a manifest in Input Validation before running</p>
              )}
            </div>
          </div>
        );

      // ── OUTPUT ─────────────────────────────────────────────────────────
      case "export":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">📦 Export Results</h2>
            <p className="q2-step-desc">QIIME2 artifacts exported to TSV/BIOM/Newick format</p>
            {outputDir ? (
              <div className="q2-info-box q2-info-box--ok">
                ✅ Output directory: <code>{outputDir}</code>
              </div>
            ) : running ? (
              <div className="q2-loading">Pipeline running… exports will appear here when done</div>
            ) : (
              <div className="q2-info-box">Run the pipeline first</div>
            )}
          </div>
        );

      case "output-files":
        return (
          <div className="q2-step-panel">
            <h2 className="q2-step-title">📂 Output Files</h2>
            {outputDir && (
              <div className="q2-output-dir">📁 {outputDir}</div>
            )}
            {outputFiles.length > 0 ? (
              <div className="q2-file-list">
                {outputFiles.map(f => (
                  <div key={f} className="q2-file-row">
                    <span className="q2-file-name">{f.split("/").pop()}</span>
                    <a href={`${API}/download/${encodeURIComponent(f)}`}
                      className="q2-file-dl" download>↓ Download</a>
                  </div>
                ))}
              </div>
            ) : running ? (
              <div className="q2-loading">Waiting for pipeline to finish…</div>
            ) : (
              <div className="q2-info-box">No files yet — run the pipeline first</div>
            )}
          </div>
        );

      default:
        return <div className="q2-step-panel"><p>Select a step from the sidebar</p></div>;
    }
  };

  // ── Main layout ───────────────────────────────────────────────────────────
  return (
    <div className="q2-manager">

      {/* ── Sidebar ──────────────────────────────────────────────────── */}
      <aside className="q2-sidebar">
        <div className="q2-sidebar-header">
          <div className="q2-sidebar-logo">Q2</div>
          <div className="q2-sidebar-title">QIIME2 Pipeline Manager</div>
        </div>

        <div className="q2-sidebar-section-label">SETUP</div>
        {SETUP_STEPS.map(item => <SidebarItem key={item.id} item={item} />)}

        <div className="q2-sidebar-section-label">PIPELINE</div>
        {PIPELINE_STEPS.map(item => <SidebarItem key={item.id} item={item} />)}

        <div className="q2-sidebar-section-label">OUTPUT</div>
        {OUTPUT_STEPS.map(item => <SidebarItem key={item.id} item={item} />)}

        {/* ── Footer progress ── */}
        <div className="q2-sidebar-footer">
          <div className="q2-footer-label">Pipeline  {doneCount}/{PIPELINE_STEPS.length}</div>
          <div className="q2-footer-bar">
            <div className="q2-footer-fill"
              style={{ width: `${(doneCount / PIPELINE_STEPS.length) * 100}%` }} />
          </div>
        </div>
      </aside>

      {/* ── Main content ──────────────────────────────────────────────── */}
      <main className="q2-main">

        {/* Top bar */}
        <div className="q2-topbar">
          {onBack && (
            <button className="q2-back-btn" onClick={onBack}>← Back</button>
          )}
          <div className="q2-topbar-status">
            {running && (
              <div className="q2-run-progress">
                <div className="q2-run-bar">
                  <div className="q2-run-fill" style={{ width: `${progress}%` }} />
                </div>
                <span className="q2-run-label">{progLabel}</span>
              </div>
            )}
            {jobStatus === "done" && <span className="q2-badge q2-badge--done">✅ Complete</span>}
            {jobStatus === "error" && <span className="q2-badge q2-badge--err">❌ Error</span>}
          </div>
        </div>

        {/* Step content */}
        <div className="q2-content">
          {renderStepContent()}
        </div>

        {/* Log panel (shown while running) */}
        {(running || logs.length > 0) && (
          <div className="q2-log-panel">
            <div className="q2-log-header">
              <span>📋 Pipeline Log</span>
              <button className="q2-log-clear" onClick={() => setLogs([])}>Clear</button>
            </div>
            <div className="q2-log-body" ref={logRef}>
              {logs.map((line, i) => <div key={i} className="q2-log-line">{line}</div>)}
              {running && <div className="q2-log-line q2-log-line--blink">▋</div>}
            </div>
          </div>
        )}
      </main>
    </div>
  );
}

// ── Helper ────────────────────────────────────────────────────────────────────
function ParamNum({ label, hint, value, onChange, min, max, step = 1, disabled = false }: {
  label: string; hint: string; value: number;
  onChange: (v: number) => void; min: number; max: number; step?: number; disabled?: boolean;
}) {
  return (
    <div className="q2-field">
      <label className="q2-label">{label}</label>
      <span className="q2-hint">{hint}</span>
      <input type="number" className="q2-input" value={value} min={min} max={max} step={step}
        disabled={disabled} onChange={e => onChange(Number(e.target.value))} />
    </div>
  );
}
