import { useState, useRef, useEffect, useCallback } from "react";
import axios from "axios";
import PipelineSettings, { PipelineParams, defaultParams, MarkerType } from "./components/PipelineSettings";
import TaxonomyColorPicker from "./components/TaxonomyColorPicker";
import DNAProgress from "./components/DNAProgress";
import MetadataEditor, { MetaRow } from "./components/MetadataEditor";
import UpdateBanner from "./components/UpdateBanner";
import LicenseModal, { LicenseStatus } from "./components/LicenseModal";
import "./App.css";

const API = "http://localhost:8000";
type Screen = "home" | "new-job" | "history";

const STATUS_COLOR: Record<string, string> = {
  uploaded: "#6b7280", queued: "#f59e0b", running: "#3b82f6",
  completed: "#10b981", error: "#ef4444", cancelled: "#9ca3af",
  waiting_checkpoint: "#f59e0b",
};
const STATUS_ICON: Record<string, string> = {
  uploaded: "📁", queued: "⏳", running: "🔬",
  completed: "✅", error: "❌", cancelled: "🚫",
  waiting_checkpoint: "⚠️",
};

interface CheckpointData {
  type:        string;
  merged_pct:  number;
  nonchim_pct: number;
  n_samples:   number;
  track:       Record<string, Record<string, number>>;
}

interface JobSummary {
  job_id:     string;
  job_name:   string;
  status:     "uploaded"|"queued"|"running"|"completed"|"error"|"cancelled"|"waiting_checkpoint";
  progress:   number;
  marker:     string;
  database:   string;
  files:      string[];
  step_label: string;
  started_at?: number;
}

interface StepRecord {
  label:      string;
  started_at: number;
  ended_at:   number | null;
  status:     "running" | "completed" | "error";
}
interface JobDetail {
  status:       string;
  started_at:   number | null;
  finished_at:  number | null;
  total_secs:   number | null;
  step_history: StepRecord[];
  error:        string;
}

// ── Helpers ──────────────────────────────────────────────────────────
const formatElapsed = (startedAt?: number): string => {
  if (!startedAt) return "—";
  const secs = Math.floor(Date.now() / 1000 - startedAt);
  return fmtSecs(secs);
};

const fmtSecs = (secs: number): string => {
  if (secs < 60) return `${Math.round(secs)}s`;
  const m = Math.floor(secs / 60);
  const s = Math.round(secs % 60);
  if (m < 60) return `${m}m ${s.toString().padStart(2, "0")}s`;
  const h = Math.floor(m / 60);
  return `${h}h ${(m % 60).toString().padStart(2, "0")}m`;
};

export default function App() {
  const [screen, setScreen]             = useState<Screen>("home");
  const [showSubmitPopup, setShowSubmitPopup] = useState(false);

  // Job form data
  const [jobName, setJobName]           = useState("");
  const [marker, setMarker]             = useState<MarkerType>("16S");
  const [params, setParams]             = useState<PipelineParams>(defaultParams);
  const [selectedFiles, setSelectedFiles] = useState<FileList|null>(null);
  const [pendingJobId, setPendingJobId] = useState<string|null>(null);
  const [sampleNames, setSampleNames]   = useState<string[]>([]);
  const [metadata, setMetadata]         = useState<MetaRow[]>([]);
  const [loading, setLoading]           = useState(false);
  const [dragOver, setDragOver]         = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);

  // System / jobs
  const [allJobs, setAllJobs]           = useState<JobSummary[]>([]);
  const [maxWorkers, setMaxWorkers]     = useState(2);
  const [focusJob, setFocusJob]         = useState<string|null>(null);
  const [focusLogs, setFocusLogs]       = useState<string[]>([]);
  const [focusLabel, setFocusLabel]     = useState("");
  const [focusCpu, setFocusCpu]         = useState<number|null>(null);
  const [focusRam, setFocusRam]         = useState<number|null>(null);
  const [focusRamTotal, setFocusRamTotal] = useState<number|null>(null);
  const [cpuLogical, setCpuLogical]     = useState<number|null>(null);
  const [cpuPhysical, setCpuPhysical]   = useState<number|null>(null);
  const [showColorPicker, setShowColorPicker] = useState<string|null>(null);
  const [detailJob, setDetailJob]   = useState<string|null>(null);
  const [detailData, setDetailData] = useState<Record<string, JobDetail>>({});
  // Delete confirmation popup
  const [deleteConfirm, setDeleteConfirm] = useState<{jobId: string; jobName: string} | null>(null);
  const [clearConfirm, setClearConfirm]   = useState(false);
  // Checkpoint warning
  const [checkpointJobId, setCheckpointJobId] = useState<string|null>(null);
  const [checkpointData, setCheckpointData]   = useState<CheckpointData|null>(null);
  const [checkpointLoading, setCheckpointLoading] = useState(false);
  const [checkpointParams, setCheckpointParams]   = useState<PipelineParams>(defaultParams);

  // ── License ───────────────────────────────────────────────────────
  const [licenseStatus, setLicenseStatus]   = useState<LicenseStatus|null>(null);
  const [showLicense,   setShowLicense]     = useState(false);
  const [licenseChecked, setLicenseChecked] = useState(false);

  const fileInputRef = useRef<HTMLInputElement>(null);
  const pollRef      = useRef<ReturnType<typeof setInterval>|null>(null);
  const progRef      = useRef<ReturnType<typeof setInterval>|null>(null);

  // ── Load step detail for a finished job ──────────────────────────
  const loadDetail = async (jobId: string) => {
    try {
      const res = await axios.get(`${API}/detail/${jobId}`);
      setDetailData(prev => ({ ...prev, [jobId]: res.data }));
    } catch {}
  };

  // ── Checkpoint handlers ────────────────────────────────────────────
  const handleCheckpointContinue = async () => {
    if (!checkpointJobId) return;
    setCheckpointLoading(true);
    try {
      await axios.post(`${API}/jobs/${checkpointJobId}/checkpoint/continue`);
      setCheckpointJobId(null);
      setCheckpointData(null);
      await refreshJobs();
    } catch { alert("Failed to continue pipeline."); }
    setCheckpointLoading(false);
  };

  const handleCheckpointAbort = async () => {
    if (!checkpointJobId) return;
    setCheckpointLoading(true);
    try {
      await axios.post(`${API}/jobs/${checkpointJobId}/checkpoint/abort`);
      setCheckpointJobId(null);
      setCheckpointData(null);
      await refreshJobs();
    } catch { alert("Failed to abort pipeline."); }
    setCheckpointLoading(false);
  };

  const handleCheckpointRerun = async () => {
    if (!checkpointJobId) return;
    setCheckpointLoading(true);
    try {
      // 1. Reset job back to "uploaded" state
      await axios.post(`${API}/jobs/${checkpointJobId}/reset`);
      // 2. Re-run with updated params
      const job = allJobs.find(j => j.job_id === checkpointJobId);
      await axios.post(`${API}/run/${checkpointJobId}`, {
        job_name: job?.job_name ?? "",
        marker:   job?.marker ?? marker,
        ...checkpointParams,
      });
      setCheckpointJobId(null);
      setCheckpointData(null);
      await refreshJobs();
    } catch { alert("Failed to re-run pipeline."); }
    setCheckpointLoading(false);
  };

  // ── Poll jobs list ────────────────────────────────────────────────
  const refreshJobs = useCallback(async () => {
    try {
      const res = await axios.get(`${API}/jobs`);
      setAllJobs(res.data.jobs);
      setMaxWorkers(res.data.max_workers);
      // Detect any job that just hit a checkpoint
      const chkJob = res.data.jobs.find((j: JobSummary) => j.status === "waiting_checkpoint");
      if (chkJob && chkJob.job_id !== checkpointJobId) {
        setCheckpointJobId(chkJob.job_id);
        setCheckpointParams(params);   // seed re-run form from current params
        // Fetch checkpoint data
        try {
          const chk = await axios.get(`${API}/jobs/${chkJob.job_id}/checkpoint`);
          setCheckpointData(chk.data.checkpoint_data);
        } catch {}
      }
    } catch {}
  }, [checkpointJobId]);

  useEffect(() => {
    refreshJobs();
    pollRef.current = setInterval(refreshJobs, 5000);
    return () => { if (pollRef.current) clearInterval(pollRef.current); };
  }, [refreshJobs]);

  // ── License check on startup ──────────────────────────────────────
  useEffect(() => {
    fetch(`${API}/license/status`)
      .then(r => r.json())
      .then((s: LicenseStatus) => {
        setLicenseStatus(s);
        setLicenseChecked(true);
        // Force modal if no valid license
        const needsKey = s.status === "no_license" || s.status === "expired" || s.status === "invalid";
        if (needsKey) setShowLicense(true);
      })
      .catch(() => {
        setLicenseChecked(true);  // backend unreachable — don't block
      });
  }, []);

  // ── Poll focused job progress + resource stats ────────────────────
  useEffect(() => {
    if (progRef.current) clearInterval(progRef.current);
    if (!focusJob) return;
    progRef.current = setInterval(async () => {
      try {
        const res = await axios.get(`${API}/progress/${focusJob}`);
        setFocusLogs(res.data.logs ?? []);
        setFocusLabel(res.data.step_label ?? "");
        setFocusCpu(res.data.cpu_pct ?? null);
        setFocusRam(res.data.ram_mb ?? null);
        setFocusRamTotal(res.data.ram_total_mb ?? null);
        if (res.data.logical  != null) setCpuLogical(res.data.logical);
        if (res.data.physical != null) setCpuPhysical(res.data.physical);
      } catch {}
    }, 2000);
    return () => { if (progRef.current) clearInterval(progRef.current); };
  }, [focusJob]);

  // ── Upload ────────────────────────────────────────────────────────
  const handleUpload = async () => {
    if (!selectedFiles?.length) { alert("Please select FASTQ files first."); return; }
    setLoading(true);
    const form = new FormData();
    Array.from(selectedFiles).forEach(f => form.append("files", f));
    try {
      const res = await axios.post(`${API}/upload`, form);
      setPendingJobId(res.data.job_id);
      const r1   = Array.from(selectedFiles).filter(f => /_R1|_1\.(fq|fastq)/i.test(f.name));
      const base = r1.length > 0 ? r1 : Array.from(selectedFiles);
      const names = base.map(f =>
        f.name.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "")
      );
      setSampleNames(names);
      setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
    } catch { alert("Upload failed — is the backend running?"); }
    setLoading(false);
  };

  // ── Run ───────────────────────────────────────────────────────────
  const handleRun = async () => {
    if (!pendingJobId) return;
    setLoading(true);
    try {
      await axios.post(`${API}/run/${pendingJobId}`, {
        job_name: jobName,
        marker,
        ...params,
        metadata,
      });
      resetWizard();
      setShowSubmitPopup(false);
      setScreen("history");
      await refreshJobs();
    } catch { alert("Failed to start analysis."); }
    setLoading(false);
  };

  // ── Reset ─────────────────────────────────────────────────────────
  const resetWizard = () => {
    setJobName(""); setSelectedFiles(null); setPendingJobId(null);
    setSampleNames([]); setMetadata([]);
    setShowAdvanced(false); setShowSubmitPopup(false);
    if (fileInputRef.current) fileInputRef.current.value = "";
  };

  // ── Cancel / delete ───────────────────────────────────────────────
  const handleCancel = async (id: string) => {
    await axios.delete(`${API}/jobs/${id}`); await refreshJobs();
  };
  const handleDeleteHistory = async (id: string) => {
    await axios.delete(`${API}/history/${id}`); await refreshJobs();
  };
  const handleClearHistory = async () => {
    await axios.delete(`${API}/history`); await refreshJobs();
  };

  // Auto-navigate to history when a checkpoint is detected
  useEffect(() => {
    if (checkpointJobId) setScreen("history");
  }, [checkpointJobId]);

  const runningCount = allJobs.filter(j => j.status === "running").length;
  const queuedCount  = allJobs.filter(j => j.status === "queued").length;
  const activeJobs   = [...allJobs].reverse().filter(j =>
    ["running","queued","uploaded","waiting_checkpoint"].includes(j.status));
  const historyJobs  = [...allJobs].reverse().filter(j =>
    ["completed","cancelled","error"].includes(j.status));

  // ── Job Card ──────────────────────────────────────────────────────
  const renderJobCard = (j: JobSummary) => {
    const displayName = j.job_name || "(unnamed)";

    return (
      <div key={j.job_id}
        className={`job-card ${focusJob === j.job_id ? "focused" : ""}`}
        style={{ borderLeft: `4px solid ${STATUS_COLOR[j.status]}` }}>

        <div className="job-card-header">
          <div className="job-card-left">
            <span className="job-status-icon">{STATUS_ICON[j.status]}</span>
            <div>
              <div className="job-id">
                <span className="job-name-text">{displayName}</span>
                <span className="job-id-badge">({j.job_id})</span>
                <span className="job-marker-tag">{j.marker}</span>
              </div>
              <div className="job-db">
                {j.database} · {j.files.length} file(s)
                {j.status === "running" && j.started_at && (
                  <span className="job-elapsed"> · ⏱ {formatElapsed(j.started_at)}</span>
                )}
              </div>
            </div>
          </div>
          <div className="job-card-right">
            <span className="job-status-label" style={{ color: STATUS_COLOR[j.status] }}>
              {j.status.toUpperCase()}
            </span>
            {(j.status === "running" || j.status === "queued" || j.status === "waiting_checkpoint") &&
              <button className="btn-cancel" onClick={() => handleCancel(j.job_id)}>✕</button>}
            {["completed","cancelled","error"].includes(j.status) &&
              <button className="btn-del-hist"
                onClick={() => setDeleteConfirm({ jobId: j.job_id, jobName: j.job_name || j.job_id })}>
                🗑
              </button>}
          </div>
        </div>

        {(j.status === "running" || j.status === "queued") && (
          <div className="job-progress-row">
            <div className="job-progress-track">
              <div className="job-progress-fill"
                style={{ width: `${j.progress}%`, background: STATUS_COLOR[j.status] }} />
            </div>
            <span className="job-progress-pct">{j.progress}%</span>
          </div>
        )}
        {j.status === "waiting_checkpoint" && (
          <div className="job-progress-row">
            <div className="job-progress-track">
              <div className="job-progress-fill chk-pulse"
                style={{ width: "60%", background: "#f59e0b" }} />
            </div>
            <span className="job-progress-pct">⏸</span>
          </div>
        )}
        {j.status === "running" && <div className="job-step">{j.step_label}</div>}
        {j.status === "waiting_checkpoint" && (
          <div className="job-step" style={{ color: "#f59e0b", fontWeight: 600 }}>
            ⚠️ Pipeline paused — low read survival detected
          </div>
        )}

        <div className="job-card-actions">
          {j.status === "waiting_checkpoint" && (
            <button className="btn-checkpoint-review"
              onClick={() => {
                setCheckpointJobId(j.job_id);
                axios.get(`${API}/jobs/${j.job_id}/checkpoint`)
                  .then(res => setCheckpointData(res.data.checkpoint_data))
                  .catch(() => {});
              }}>
              ⚠️ Review Warning &amp; Decide
            </button>
          )}
          {j.status === "running" && (
            <button className="btn-view"
              onClick={() => {
                if (focusJob === j.job_id) {
                  setFocusJob(null);
                  setFocusCpu(null); setFocusRam(null);
                } else {
                  setFocusJob(j.job_id);
                }
              }}>
              {focusJob === j.job_id ? "▲ Hide" : "▼ View Progress"}
            </button>
          )}
          {j.status === "completed" && (
            <>
              <button className="btn-view"
                onClick={() => setShowColorPicker(showColorPicker === j.job_id ? null : j.job_id)}>
                🎨 Taxonomy Colors
              </button>
              <a href={`${API}/download-zip/${j.job_id}`} download className="btn-view">
                📥 Download Results
              </a>
            </>
          )}
          {(j.status === "completed" || j.status === "error") && (
            <button className="btn-view btn-detail"
              onClick={() => {
                if (detailJob === j.job_id) {
                  setDetailJob(null);
                } else {
                  setDetailJob(j.job_id);
                  loadDetail(j.job_id);
                }
              }}>
              {detailJob === j.job_id ? "▲ Hide Detail" : "📊 Detail"}
            </button>
          )}
        </div>

        {/* ── Expanded: DNA helix + Resource stats (running) ── */}
        {focusJob === j.job_id && j.status === "running" && (
          <div className="job-dna-expand">
            <div className="job-expand-row">
              <div className="job-dna-left">
                <DNAProgress percent={j.progress} logs={focusLogs}
                  currentStep={focusLabel || j.step_label} />
              </div>
              <div className="job-stats-panel">
                <div className="stat-item">
                  <span className="stat-icon">⏱️</span>
                  <span className="stat-label">RUNNING TIME</span>
                  <strong className="stat-value">{formatElapsed(j.started_at)}</strong>
                </div>
                <div className="stat-item">
                  <span className="stat-icon">🖥️</span>
                  <span className="stat-label">
                    CPU USAGE
                    {cpuPhysical != null && cpuLogical != null && (
                      <span className="stat-sublabel">
                        {cpuPhysical}c/{cpuLogical}t
                      </span>
                    )}
                  </span>
                  <strong className="stat-value">
                    {focusCpu != null ? `${focusCpu.toFixed(1)}%` : "—"}
                  </strong>
                  {focusCpu != null && cpuLogical != null && (
                    <div className="stat-bar-wrap">
                      <div className="stat-bar"
                        style={{ width: `${Math.min(focusCpu, 100)}%`,
                                 background: focusCpu > 80 ? "#ef4444"
                                           : focusCpu > 50 ? "#f59e0b"
                                           : "#4f46e5" }} />
                    </div>
                  )}
                </div>
                <div className="stat-item">
                  <span className="stat-icon">💾</span>
                  <span className="stat-label">RAM USED</span>
                  <strong className="stat-value">
                    {focusRam != null
                      ? focusRam >= 1024
                        ? `${(focusRam / 1024).toFixed(1)} GB`
                        : `${focusRam} MB`
                      : "—"}
                  </strong>
                  {focusRam != null && focusRamTotal != null && (
                    <>
                      <span className="stat-of">
                        of {focusRamTotal >= 1024
                          ? `${(focusRamTotal / 1024).toFixed(1)} GB`
                          : `${focusRamTotal} MB`}
                      </span>
                      <div className="stat-bar-wrap">
                        <div className="stat-bar"
                          style={{
                            width: `${Math.min((focusRam / focusRamTotal) * 100, 100)}%`,
                            background: (focusRam / focusRamTotal) > 0.85 ? "#ef4444"
                                      : (focusRam / focusRamTotal) > 0.65 ? "#f59e0b"
                                      : "#10b981"
                          }} />
                      </div>
                    </>
                  )}
                </div>
              </div>
            </div>
          </div>
        )}

        {/* ── Expanded: Step Detail (completed / error) ── */}
        {detailJob === j.job_id && (
          <div className="job-detail-panel">
            {!detailData[j.job_id] ? (
              <div className="detail-loading">⏳ Loading...</div>
            ) : (() => {
              const d = detailData[j.job_id];
              return (
                <>
                  <div className="detail-total-row">
                    <span className="detail-total-icon">⏱️</span>
                    <span className="detail-total-label">Total time</span>
                    <strong className="detail-total-value">
                      {d.total_secs != null ? fmtSecs(d.total_secs) : "—"}
                    </strong>
                  </div>

                  {d.step_history.length > 0 ? (
                    <div className="detail-step-list">
                      <div className="detail-step-header">
                        <span>Step</span>
                        <span>Duration</span>
                        <span>Status</span>
                      </div>
                      {d.step_history.map((s, i) => {
                        const dur = s.ended_at && s.started_at
                          ? fmtSecs(s.ended_at - s.started_at) : "—";
                        const icon = s.status === "completed" ? "✅"
                          : s.status === "error" ? "❌" : "🔄";
                        const statusText = s.status === "completed" ? "Completed"
                          : s.status === "error" ? "Error" : "Running";
                        return (
                          <div key={i}
                            className={`detail-step-row detail-step-${s.status}`}>
                            <span className="detail-step-name">
                              <span className="detail-step-icon">{icon}</span>
                              {s.label}
                            </span>
                            <span className="detail-step-dur">{dur}</span>
                            <span className={`detail-step-status detail-status-${s.status}`}>
                              {statusText}
                            </span>
                          </div>
                        );
                      })}
                    </div>
                  ) : (
                    <p className="detail-no-steps">
                      ℹ️ Step breakdown not available — re-run this job to capture step details.
                    </p>
                  )}

                  {d.error && (
                    <div className="detail-error-section">
                      <div className="detail-error-title">❌ Error log</div>
                      <pre className="detail-error-body">{d.error}</pre>
                    </div>
                  )}
                </>
              );
            })()}
          </div>
        )}

        {showColorPicker === j.job_id && (
          <TaxonomyColorPicker jobId={j.job_id} apiBase={API} />
        )}
        {j.status === "error" && detailJob !== j.job_id && (
          <div className="job-error-hint">
            ❌ Job failed — click <strong>📊 Detail</strong> to see the error log
          </div>
        )}
      </div>
    );
  };

  // ══════════════════════════════════════════════════════════════
  //  HOME SCREEN
  // ══════════════════════════════════════════════════════════════
  if (screen === "home") {
    return (
      <div className="home-screen">
        {checkpointJobId && (
          <div className="chk-home-banner" onClick={() => setScreen("history")}>
            ⚠️ A pipeline has paused and needs your decision — <strong>click to review</strong>
          </div>
        )}
        <div className="home-logo-card">
          <span className="home-dna">🧬</span>
          <h1 className="home-title">16S / 12S Amplicon Analysis</h1>
          <p className="home-desc">
            DADA2 Pipeline — Professional Microbiome &amp; Metabarcoding Analysis
          </p>
        </div>
        <div className="home-panels">
          <button className="home-panel new-job-panel"
            onClick={() => { resetWizard(); setScreen("new-job"); }}>
            <span className="panel-icon">➕</span>
            <span className="panel-title">New Job</span>
            <span className="panel-sub">Start analysis / upload FASTQ files</span>
          </button>
          <button className="home-panel history-panel"
            onClick={() => setScreen("history")}>
            <span className="panel-icon">📋</span>
            <span className="panel-title">Job History</span>
            <span className="panel-sub">
              {allJobs.length} total job{allJobs.length !== 1 ? "s" : ""}
              {(runningCount + queuedCount) > 0
                ? ` · ${runningCount + queuedCount} active` : ""}
            </span>
            {(runningCount + queuedCount) > 0 &&
              <span className="panel-badge">{runningCount + queuedCount}</span>}
          </button>
        </div>
      </div>
    );
  }

  // ══════════════════════════════════════════════════════════════
  //  HISTORY SCREEN
  // ══════════════════════════════════════════════════════════════
  if (screen === "history") {
    return (
      <div className="hist-screen">
        <div className="hist-topbar">
          <button className="back-btn" onClick={() => setScreen("home")}>← Home</button>
          <h2 className="hist-title">📋 Job History</h2>
          <div className="hist-actions">
            <span className="worker-info">
              🖥️ {runningCount}/{maxWorkers} running
              {queuedCount > 0 && <span className="queue-badge">⏳ {queuedCount}</span>}
            </span>
            {historyJobs.length > 0 &&
              <button className="btn-clear-hist" onClick={() => setClearConfirm(true)}>🗑 Clear History</button>}
            <button className="btn-refresh" onClick={refreshJobs}>↻ Refresh</button>
          </div>
        </div>
        <div className="hist-body">
          {allJobs.length === 0 ? (
            <div className="card empty-state">
              <p>No jobs yet — go to <strong>New Job</strong> to start an analysis.</p>
            </div>
          ) : (
            <>
              {activeJobs.length > 0 && (
                <>
                  <div className="hist-section-label">🔄 Active ({activeJobs.length})</div>
                  <div className="job-list">{activeJobs.map(renderJobCard)}</div>
                </>
              )}
              {historyJobs.length > 0 && (
                <>
                  <div className="hist-section-label">🕒 Completed ({historyJobs.length})</div>
                  <div className="job-list">{historyJobs.map(renderJobCard)}</div>
                </>
              )}
            </>
          )}
        </div>

        {/* ── Delete single job confirmation ── */}
        {deleteConfirm && (
          <div className="popup-overlay">
            <div className="step-popup del-confirm-popup">
              <div className="popup-step-icon">🗑️</div>
              <h2 className="popup-title">ยืนยันการลบ</h2>
              <p className="del-confirm-msg">
                คุณต้องการลบโปรเจค <strong>"{deleteConfirm.jobName}"</strong> ออกจาก history
                และ<span className="del-warn-text"> ลบไฟล์ผลลัพธ์ทั้งหมดบนดิสก์</span> ด้วยหรือไม่?
              </p>
              <p className="del-confirm-sub">⚠️ การดำเนินการนี้ไม่สามารถย้อนกลับได้</p>
              <div className="popup-footer">
                <button className="popup-back" onClick={() => setDeleteConfirm(null)}>← ยกเลิก</button>
                <button className="btn-del-confirm"
                  onClick={async () => {
                    await handleDeleteHistory(deleteConfirm.jobId);
                    setDeleteConfirm(null);
                  }}>
                  🗑 ลบถาวร
                </button>
              </div>
            </div>
          </div>
        )}

        {/* ── Clear all history confirmation ── */}
        {clearConfirm && (
          <div className="popup-overlay">
            <div className="step-popup del-confirm-popup">
              <div className="popup-step-icon">🗑️</div>
              <h2 className="popup-title">ลบ History ทั้งหมด</h2>
              <p className="del-confirm-msg">
                คุณต้องการลบ <strong>ทุกโปรเจคที่เสร็จแล้ว</strong> ออกจาก history
                และลบไฟล์ผลลัพธ์บนดิสก์ทั้งหมดหรือไม่?
              </p>
              <p className="del-confirm-sub">⚠️ การดำเนินการนี้ไม่สามารถย้อนกลับได้</p>
              <div className="popup-footer">
                <button className="popup-back" onClick={() => setClearConfirm(false)}>← ยกเลิก</button>
                <button className="btn-del-confirm"
                  onClick={async () => {
                    await handleClearHistory();
                    setClearConfirm(false);
                  }}>
                  🗑 ลบทั้งหมด
                </button>
              </div>
            </div>
          </div>
        )}

        {/* ── Checkpoint Warning Modal ── */}
        {checkpointJobId && (
          <div className="popup-overlay">
            <div className="step-popup checkpoint-popup">
              <div className="popup-step-icon">⚠️</div>
              <h2 className="popup-title" style={{ color: "#f59e0b" }}>
                Low Read Survival — Pipeline Paused
              </h2>

              {checkpointData ? (
                <>
                  <div className="checkpoint-warning-box">
                    <p>
                      Only <strong style={{ color: "#ef4444" }}>{checkpointData.merged_pct.toFixed(1)}%</strong> of
                      reads survived merging and <strong style={{ color: "#ef4444" }}>{checkpointData.nonchim_pct.toFixed(1)}%</strong> survived
                      chimera removal. This is far below the expected ≥70%.
                    </p>
                    <p style={{ marginTop: 6 }}>
                      💡 <strong>Common causes:</strong> incorrect <code>truncLen</code> (reads cut too short to overlap),
                      wrong <code>trimLeft</code> (primer not fully removed), or poor sequencing quality.
                    </p>
                  </div>

                  <div className="checkpoint-track-wrap">
                    <table className="checkpoint-track-table">
                      <thead>
                        <tr>
                          <th>Sample</th>
                          <th>Input</th>
                          <th>Filtered</th>
                          <th>Merged</th>
                          <th>Non-chimeric</th>
                          <th>% Retained</th>
                        </tr>
                      </thead>
                      <tbody>
                        {Object.entries(checkpointData.track).map(([sample, row]) => {
                          const pct = row.input > 0 ? (row.nonchim / row.input * 100) : 0;
                          return (
                            <tr key={sample} className={pct < 10 ? "chk-row-warn" : ""}>
                              <td className="chk-sample-cell">{sample}</td>
                              <td>{(row.input ?? 0).toLocaleString()}</td>
                              <td>{(row.filtered ?? 0).toLocaleString()}</td>
                              <td>{(row.merged ?? 0).toLocaleString()}</td>
                              <td>{(row.nonchim ?? 0).toLocaleString()}</td>
                              <td style={{ fontWeight: 700,
                                color: pct < 10 ? "#ef4444" : pct < 30 ? "#f59e0b" : "#10b981" }}>
                                {pct.toFixed(1)}%
                              </td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </>
              ) : (
                <div style={{ textAlign: "center", padding: "24px 0", color: "#9ca3af" }}>
                  ⏳ Loading checkpoint data...
                </div>
              )}

              {/* ── Re-run parameter adjustment form ── */}
              <div className="chk-rerun-section">
                <div className="chk-rerun-title">🔧 ปรับค่าและรันใหม่</div>
                <p className="chk-rerun-hint">
                  แก้ไขค่า trim / truncation ด้านล่าง แล้วกด <strong>Re-run with New Settings</strong><br/>
                  ไฟล์ที่อัปโหลดไว้จะถูกนำมาใช้ใหม่โดยไม่ต้องอัปโหลดซ้ำ
                </p>
                <div className="chk-param-grid">

                  <div className="chk-param-group">
                    <div className="chk-param-group-label">✂️ Trim Left (bp)</div>
                    <div className="chk-param-row">
                      <label>Forward (trimLeft_F)</label>
                      <input type="number" min={0} max={50}
                        value={checkpointParams.trimLeft_F}
                        onChange={e => setCheckpointParams(p => ({ ...p, trimLeft_F: +e.target.value }))} />
                    </div>
                    <div className="chk-param-row">
                      <label>Reverse (trimLeft_R)</label>
                      <input type="number" min={0} max={50}
                        value={checkpointParams.trimLeft_R}
                        onChange={e => setCheckpointParams(p => ({ ...p, trimLeft_R: +e.target.value }))} />
                    </div>
                  </div>

                  <div className="chk-param-group">
                    <div className="chk-param-group-label">✂️ Truncate Length (bp)</div>
                    <div className="chk-param-row">
                      <label>Forward (truncLen_F)</label>
                      <input type="number" min={0} max={300}
                        value={checkpointParams.truncLen_F}
                        onChange={e => setCheckpointParams(p => ({ ...p, truncLen_F: +e.target.value }))} />
                    </div>
                    <div className="chk-param-row">
                      <label>Reverse (truncLen_R)</label>
                      <input type="number" min={0} max={300}
                        value={checkpointParams.truncLen_R}
                        onChange={e => setCheckpointParams(p => ({ ...p, truncLen_R: +e.target.value }))} />
                    </div>
                  </div>

                  <div className="chk-param-group">
                    <div className="chk-param-group-label">📉 Max Expected Errors</div>
                    <div className="chk-param-row">
                      <label>Forward (maxEE_F)</label>
                      <input type="number" min={1} max={10} step={0.5}
                        value={checkpointParams.maxEE_F}
                        onChange={e => setCheckpointParams(p => ({ ...p, maxEE_F: +e.target.value }))} />
                    </div>
                    <div className="chk-param-row">
                      <label>Reverse (maxEE_R)</label>
                      <input type="number" min={1} max={10} step={0.5}
                        value={checkpointParams.maxEE_R}
                        onChange={e => setCheckpointParams(p => ({ ...p, maxEE_R: +e.target.value }))} />
                    </div>
                  </div>

                </div>
              </div>

              <div className="popup-footer chk-footer-3btn" style={{ marginTop: 20 }}>
                <button className="btn-checkpoint-abort"
                  onClick={handleCheckpointAbort}
                  disabled={checkpointLoading}>
                  {checkpointLoading ? "⏳..." : "🗑 Abort"}
                </button>
                <button className="btn-checkpoint-continue"
                  onClick={handleCheckpointContinue}
                  disabled={checkpointLoading}>
                  {checkpointLoading ? "⏳..." : "▶ Continue Anyway"}
                </button>
                <button className="btn-checkpoint-rerun"
                  onClick={handleCheckpointRerun}
                  disabled={checkpointLoading}>
                  {checkpointLoading ? "⏳ Re-running..." : "🔄 Re-run with New Settings"}
                </button>
              </div>
            </div>
          </div>
        )}
      </div>
    );
  }

  // ══════════════════════════════════════════════════════════════
  //  NEW JOB — Single scrollable page
  // ══════════════════════════════════════════════════════════════
  return (
    <div className="wizard-wrap">

      {/* ── Update Banner ── */}
      <UpdateBanner />

      {/* ── License Modal (first-run gate or manual open) ── */}
      {showLicense && (
        <LicenseModal
          required={
            licenseStatus?.status === "no_license" ||
            licenseStatus?.status === "expired"    ||
            licenseStatus?.status === "invalid"
          }
          onClose={() => setShowLicense(false)}
          onActivated={(s) => {
            setLicenseStatus(s);
            if (s.status === "active") setShowLicense(false);
          }}
        />
      )}

      {/* ── CPU / Submit Popup ── */}
      {showSubmitPopup && (
        <div className="popup-overlay">
          <div className="step-popup confirm-popup">
            <div className="popup-step-icon">🚀</div>
            <h2 className="popup-title">Submit Analysis Job</h2>

            <div className="confirm-summary">
              <div className="confirm-row">
                <span>Job name</span>
                <strong>{jobName || <em style={{ color: "#9ca3af" }}>unnamed</em>}</strong>
              </div>
              <div className="confirm-row">
                <span>Marker</span><strong>{marker}</strong>
              </div>
              <div className="confirm-row">
                <span>Files uploaded</span>
                <strong>{selectedFiles?.length ?? 0} file(s)</strong>
              </div>
              <div className="confirm-row">
                <span>Samples</span><strong>{sampleNames.length}</strong>
              </div>
              {metadata.filter(m => m.group).length > 0 && (
                <div className="confirm-row">
                  <span>Groups</span>
                  <strong>
                    {Array.from(new Set(metadata.map(m => m.group).filter(Boolean))).join(", ")}
                  </strong>
                </div>
              )}
            </div>

            <div className="cpu-section">
              <label className="cpu-label">🖥️ CPU Workers</label>
              <div className="worker-controls">
                {[1, 2, 4, 6, 8].map(n => (
                  <button key={n}
                    className={`worker-btn ${maxWorkers === n ? "active" : ""}`}
                    onClick={async () => {
                      setMaxWorkers(n);
                      await axios.post(`${API}/config`, { max_workers: n });
                    }}>{n}</button>
                ))}
              </div>
              <p className="cpu-hint">
                {runningCount}/{maxWorkers} running · {queuedCount} queued
                {runningCount >= maxWorkers && " — this job will be queued"}
              </p>
            </div>

            <div className="popup-footer">
              <button className="popup-back" onClick={() => setShowSubmitPopup(false)}>
                ← Back
              </button>
              <button className="popup-next" onClick={handleRun} disabled={loading}>
                {loading ? "⏳ Submitting..." : "🚀 Submit Job"}
              </button>
            </div>
          </div>
        </div>
      )}

      {/* ── Main new-job page ── */}
      <div className="wizard-page">

        {/* Header */}
        <div className="wizard-header">
          <div className="logo-wrap">
            <span className="logo-dna">🧬</span>
            <div>
              <div className="app-title">NextGen-Amplicon</div>
              <div className="app-subtitle">16S / ITS / COX1 / PacBio Microbiome Pipeline</div>
            </div>
          </div>
          <div className="header-actions">
            <button
              className="btn-history"
              title="License"
              style={{
                background: licenseStatus?.status === "active" || licenseStatus?.status === "dev"
                  ? "rgba(16,185,129,0.15)" : "rgba(239,68,68,0.15)",
              }}
              onClick={() => setShowLicense(true)}
            >
              🔐 {licenseStatus?.status === "active" || licenseStatus?.status === "dev"
                ? `Licensed (${licenseStatus?.days_remaining != null ? licenseStatus.days_remaining + "d" : "∞"})`
                : "Enter License"}
            </button>
            <button className="btn-history" onClick={() => setScreen("history")}>
              🗓 View Jobs
            </button>
          </div>
        </div>

        {/* Job name */}
        <div className="section-card">
          <div className="section-title">💼 Job Name</div>
          <input
            className="job-name-input"
            placeholder="e.g. Soil_16S_Run1 (optional)"
            value={jobName}
            onChange={e => setJobName(e.target.value)}
          />
        </div>

        {/* Marker pill + Pipeline settings */}
        <div className="section-card">
          <div className="section-title">🧬 Pipeline Settings</div>
          {marker && (
            <div className="marker-selected-pill">
              Selected: <strong>{marker}</strong>
            </div>
          )}
          <PipelineSettings
            params={params}
            onChange={setParams}
            marker={marker}
            onMarker={setMarker}
          />
        </div>

        {/* Upload */}
        <div className="section-card">
          <div className="section-title">📂 Upload FASTQ Files</div>
          <div
            className={`drop-zone ${dragOver ? "drag-over" : ""}`}
            onDragOver={e => { e.preventDefault(); setDragOver(true); }}
            onDragLeave={() => setDragOver(false)}
            onDrop={e => {
              e.preventDefault(); setDragOver(false);
              const dt = e.dataTransfer.files;
              if (dt.length) {
                setSelectedFiles(dt);
                const r1 = Array.from(dt).filter(f => /_R1|_1\.(fq|fastq)/i.test(f.name));
                const base = r1.length > 0 ? r1 : Array.from(dt);
                const names = base.map(f =>
                  f.name.replace(/_R1.*|_1\.(fq|fastq).*/i,"").replace(/\.(fastq|fq)(\.gz)?$/i,"")
                );
                setSampleNames(names);
                setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
              }
            }}
            onClick={() => fileInputRef.current?.click()}
          >
            {selectedFiles?.length
              ? <>📁 <strong>{selectedFiles.length} file(s)</strong> selected — click to change</>
              : <>📂 Drop .fastq / .fastq.gz files here, or <u>click to browse</u></>
            }
          </div>
          <input
            ref={fileInputRef} type="file" multiple accept=".fastq,.fastq.gz,.fq,.fq.gz"
            style={{ display: "none" }}
            onChange={e => {
              const files = e.target.files;
              if (files?.length) {
                setSelectedFiles(files);
                const r1 = Array.from(files).filter(f => /_R1|_1\.(fq|fastq)/i.test(f.name));
                const base = r1.length > 0 ? r1 : Array.from(files);
                const names = base.map(f =>
                  f.name.replace(/_R1.*|_1\.(fq|fastq).*/i,"").replace(/\.(fastq|fq)(\.gz)?$/i,"")
                );
                setSampleNames(names);
                setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
              }
            }}
          />

          {sampleNames.length > 0 && (
            <div className="sample-list">
              <div className="sample-list-title">🧪 Detected samples ({sampleNames.length})</div>
              <div className="sample-chips">
                {sampleNames.map(s => <span key={s} className="sample-chip">{s}</span>)}
              </div>
            </div>
          )}
        </div>

        {/* Metadata */}
        {metadata.length > 0 && (
          <div className="section-card">
            <div className="section-title">📊 Sample Metadata (optional)</div>
            <MetadataEditor sampleNames={sampleNames} onChange={setMetadata} />
          </div>
        )}

        {/* Advanced toggle */}
        <div className="section-card">
          <button className="btn-advanced-toggle"
            onClick={() => setShowAdvanced(v => !v)}>
            {showAdvanced ? "▲ Hide Advanced" : "▼ Show Advanced Settings"}
          </button>
          {showAdvanced && (
            <div className="advanced-section">
              <TaxonomyColorPicker
                jobId={focusJob ?? ""}
                apiBase={API}
              />
            </div>
          )}
        </div>

        {/* Upload + Submit buttons */}
        <div className="wizard-actions">
          {!pendingJobId ? (
            <button className="btn-upload" onClick={handleUpload} disabled={loading || !selectedFiles?.length}>
              {loading ? "⏳ Uploading..." : "⬆️ Upload Files"}
            </button>
          ) : (
            <div className="upload-done-row">
              <span className="upload-done-badge">✅ Files uploaded</span>
              <button className="btn-submit" onClick={() => setShowSubmitPopup(true)}>
                🚀 Configure & Submit →
              </button>
            </div>
          )}
        </div>

      </div>
    </div>
  );
}
