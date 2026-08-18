import { useState, useRef, useEffect, useCallback } from "react";
import axios from "axios";
import PipelineSettings, { PipelineParams, defaultParams, MarkerType } from "./components/PipelineSettings";
import TaxonomyColorPicker from "./components/TaxonomyColorPicker";
import DNAProgress from "./components/DNAProgress";
import MetadataEditor, { MetaRow } from "./components/MetadataEditor";
import UpdateBanner from "./components/UpdateBanner";
import LicenseModal, { LicenseStatus } from "./components/LicenseModal";
import SettingsPanel, { ThemeId } from "./components/SettingsPanel";
import PreviewPage from "./pages/PreviewPage";
import "./App.css";

const API = "http://localhost:8000";
type Screen = "home" | "new-job" | "history" | "preview";

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
  params?:      Record<string, any>;
  run_at?:      number | null;
}

// ── Friendly labels + which param keys to surface in Job Detail ────────
// (only shown when the value is non-empty/meaningful for that job's marker)
const PARAM_LABELS: Record<string, string> = {
  marker: "Marker", sequencerType: "Sequencer Type",
  truncLen_F: "Truncate Fwd (bp)", truncLen_R: "Truncate Rev (bp)",
  ontMinLen: "ONT Min Length (bp)", ontMaxLen: "ONT Max Length (bp)",
  maxEE_F: "Max EE Fwd", maxEE_R: "Max EE Rev",
  trimLeft_F: "Trim Left Fwd", trimLeft_R: "Trim Left Rev",
  pool: "Pooling", chimeraMethod: "Chimera Method",
  taxDatabase: "Taxonomy DB", dbPath: "DB Path", minBoot: "Min Bootstrap",
  primer_f: "Forward Primer", primer_r: "Reverse Primer",
  cutadaptErrorRate: "Cutadapt Error Rate", cutadaptOverlap: "Cutadapt Min Overlap",
  discardUntrimmed: "Discard Untrimmed",
  its_region: "ITS Region",
  truncLen_cox1_f: "COX1 Truncate Fwd", truncLen_cox1_r: "COX1 Truncate Rev",
  codon_table: "Codon Table", cox1_min_len: "COX1 Min Len", cox1_max_len: "COX1 Max Len",
  run_lulu: "LULU Curation",
  pb_min_len: "PacBio Min Len", pb_max_len: "PacBio Max Len",
  pb_maxEE: "PacBio Max EE", pb_region: "PacBio Region",
  ont_region: "16S Region (ONT)", ont_min_abundance: "Min Abundance Filter",
  ont_db_path: "Reference Database",
  otuSimilarity: "OTU Similarity",
  run_tax4fun: "Tax4Fun2", run_picrust2: "PICRUSt2",
  customClassifierMode: "Classifier Mode", customClassifierPath: "Classifier Path",
  nThreads: "CPU Threads",
};
const PARAM_KEY_ORDER = Object.keys(PARAM_LABELS);

// Keys that only make sense for certain markers/pipelines — shown in Job Detail
// only when they were actually consulted by the pipeline that ran this job.
// (e.g. truncLen/maxEE only apply to the DADA2 branch — a QIIME2/VSEARCH job
// under the same marker (16S/12S/18S-NEMA) never reads those fields, so
// showing "Truncate Fwd: 120" there is just a leftover value from whatever
// mode was last configured in the shared params object, not something that
// actually affected this run.)
// These keys are only meaningful when sequencerType === "illumina" (the DADA2 branch).
const ILLUMINA_ONLY_KEYS = new Set([
  "truncLen_F", "truncLen_R",
  "maxEE_F", "maxEE_R", "trimLeft_F", "trimLeft_R", "pool", "chimeraMethod",
]);
// ontMinLen/ontMaxLen are read by BOTH the "ont" sequencerType (as read-length
// filters) and the "qiime2_vsearch" sequencerType (as --min_len/--max_len for
// vsearch), so they're relevant for either — just not for "illumina".
const ONT_OR_VSEARCH_LEN_KEYS = new Set(["ontMinLen", "ontMaxLen"]);
// otuSimilarity is only read by the QIIME2/VSEARCH clustering step.
const VSEARCH_ONLY_KEYS = new Set(["otuSimilarity"]);
const ITS_ONLY_KEYS    = new Set(["its_region", "run_lulu"]);
const COX1_ONLY_KEYS   = new Set(["truncLen_cox1_f", "truncLen_cox1_r", "codon_table", "cox1_min_len", "cox1_max_len", "run_lulu"]);
const PACBIO_ONLY_KEYS = new Set(["pb_min_len", "pb_max_len", "pb_maxEE", "pb_region"]);
const ONT16S_ONLY_KEYS = new Set(["ont_region", "ont_min_abundance", "ont_db_path"]);
// Markers that can run either the DADA2 (illumina) or QIIME2/VSEARCH pipeline,
// switched via sequencerType — this is where the truncLen leak bug lived.
const DADA2_OR_VSEARCH_MARKERS = new Set(["16S", "12S", "18S-NEMA"]);

function isParamRelevantForMarker(key: string, marker: string, sequencerType?: string): boolean {
  const m  = (marker || "").toUpperCase();
  const st = sequencerType || "illumina";
  // ont_db_path is shared: used by the ONT-16S (Emu) marker AND by the
  // "QIIME2 (VSEARCH OTU)" sequencer-type option under 16S/12S/18S-nema.
  if (key === "ont_db_path") return m === "ONT-16S" || m === "ONT16S" || m === "ONT" ||
                                     (DADA2_OR_VSEARCH_MARKERS.has(m) && st === "qiime2_vsearch");
  if (key === "sequencerType") return DADA2_OR_VSEARCH_MARKERS.has(m);
  if (ILLUMINA_ONLY_KEYS.has(key))
    return DADA2_OR_VSEARCH_MARKERS.has(m) && st === "illumina";
  if (ONT_OR_VSEARCH_LEN_KEYS.has(key))
    return DADA2_OR_VSEARCH_MARKERS.has(m) && (st === "ont" || st === "qiime2_vsearch");
  if (VSEARCH_ONLY_KEYS.has(key))
    return DADA2_OR_VSEARCH_MARKERS.has(m) && st === "qiime2_vsearch";
  if (ITS_ONLY_KEYS.has(key))    return m === "ITS1" || m === "ITS2" || m === "ITS";
  if (COX1_ONLY_KEYS.has(key))   return m === "COX1";
  if (PACBIO_ONLY_KEYS.has(key)) return m === "PACBIO";
  if (ONT16S_ONLY_KEYS.has(key)) return m === "ONT-16S" || m === "ONT16S" || m === "ONT";
  return true; // shared/common keys (marker, taxDatabase, dbPath, primers, threads, etc.)
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
  const [previewJobId, setPreviewJobId] = useState<string>("");
  const [showSubmitPopup, setShowSubmitPopup] = useState(false);

  // Job form data
  const [jobName, setJobName]           = useState("");
  const [marker, setMarker]             = useState<MarkerType>("16S");
  const [params, setParams]             = useState<PipelineParams>(defaultParams);
  const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
  const [pendingJobId, setPendingJobId] = useState<string|null>(null);
  const [sampleNames, setSampleNames]   = useState<string[]>([]);
  const [metadata, setMetadata]         = useState<MetaRow[]>([]);
  const [loading, setLoading]           = useState(false);
  const [dragOver, setDragOver]         = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);
  // Manual sample↔file pairing (bypasses filename auto-detection)
  const [serverFileList, setServerFileList] = useState<string[]>([]);
  const [useManualPairing, setUseManualPairing] = useState(false);
  const [pairReadMode, setPairReadMode] = useState<"single"|"paired">("paired");
  const [fileMap, setFileMap] = useState<{sample:string; file1:string; file2:string}[]>([]);

  // Keep "Detected samples" / Sample Metadata in sync with manual pairing edits
  useEffect(() => {
    if (!useManualPairing) return;
    const names = fileMap.map(r => r.sample.trim()).filter(Boolean);
    if (names.length > 0) setSampleNames(names);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [useManualPairing, fileMap.map(r => r.sample).join("|")]);

  // Revert "Detected samples" back to filename auto-detect when manual pairing is turned off
  useEffect(() => {
    if (useManualPairing || serverFileList.length === 0) return;
    const isONT = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
    const r1 = isONT ? [] : serverFileList.filter(n => /_R1|_1\.(fq|fastq)/i.test(n));
    const base = r1.length > 0 ? r1 : serverFileList.filter(n => /\.(fastq|fq)(\.gz)?$/i.test(n));
    const names = base.map(n =>
      isONT
        ? n.replace(/\.(fastq|fq)(\.gz)?$/i, "")
        : n.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "")
    );
    if (names.length > 0) setSampleNames(names);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [useManualPairing]);

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
  const [sideFilesExpanded, setSideFilesExpanded] = useState(false);
  // Delete confirmation popup
  const [deleteConfirm, setDeleteConfirm] = useState<{jobId: string; jobName: string} | null>(null);
  const [clearConfirm, setClearConfirm]   = useState(false);
  const [cleanupConfirm, setCleanupConfirm] = useState(false);
  const [cleanupResult, setCleanupResult]   = useState<{removed: number; paths: string[]} | null>(null);
  // Checkpoint warning
  const [checkpointJobId, setCheckpointJobId] = useState<string|null>(null);
  const [checkpointData, setCheckpointData]   = useState<CheckpointData|null>(null);
  const [checkpointLoading, setCheckpointLoading] = useState(false);
  const [checkpointParams, setCheckpointParams]   = useState<PipelineParams>(defaultParams);

  // ── License ───────────────────────────────────────────────────────
  const [licenseStatus, setLicenseStatus]   = useState<LicenseStatus|null>(null);
  const [showLicense,   setShowLicense]     = useState(false);
  const [licenseChecked, setLicenseChecked] = useState(false);

  // ── Theme ─────────────────────────────────────────────────────────
  const [theme, setTheme] = useState<ThemeId>(() =>
    (localStorage.getItem("app-theme") as ThemeId) || "default"
  );
  const [showSettings, setShowSettings] = useState(false);

  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme === "default" ? "" : theme);
    localStorage.setItem("app-theme", theme);
  }, [theme]);

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
        // Seed the re-run form from THIS job's own saved params — not the
        // top-level `params` state, which reflects whatever is currently
        // loaded in the main Settings form and may belong to a totally
        // different job (or still be untouched defaults). Using the wrong
        // source here was silently dropping primers/database/etc on
        // "Re-run with New Settings" for any job that wasn't the one
        // currently open in Settings.
        try {
          const jobDetail = await axios.get(`${API}/detail/${chkJob.job_id}`);
          if (jobDetail.data?.params) {
            setCheckpointParams({ ...defaultParams, ...jobDetail.data.params });
          } else {
            setCheckpointParams(params);
          }
        } catch { setCheckpointParams(params); }
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

  // ── Right-hand detail sidebar: which job is it showing? ────────────
  // Priority: whatever's expanded via "View Progress" (focusJob) or
  // "Detail" (detailJob), else the first active job, else the most
  // recent job overall — so the sidebar is never blank if any job exists.
  const sidebarJobId = focusJob || detailJob ||
    allJobs.find(j => ["running", "queued", "waiting_checkpoint"].includes(j.status))?.job_id ||
    (allJobs.length > 0 ? allJobs[allJobs.length - 1].job_id : null);
  const sidebarJobStatus = allJobs.find(j => j.job_id === sidebarJobId)?.status;

  // ── Keep the sidebar's step/settings data fresh (live while running) ──
  useEffect(() => {
    if (!sidebarJobId) return;
    loadDetail(sidebarJobId);
    if (!sidebarJobStatus || !["running", "queued", "waiting_checkpoint"].includes(sidebarJobStatus)) return;
    const iv = setInterval(() => loadDetail(sidebarJobId), 3000);
    return () => clearInterval(iv);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sidebarJobId, sidebarJobStatus]);

  // Collapse the "show all files" toggle whenever the sidebar switches to a different job
  useEffect(() => { setSideFilesExpanded(false); }, [sidebarJobId]);

  // ── Select a job from the left nav list ─────────────────────────────
  const selectJobFromNav = (j: JobSummary) => {
    // sidebarJobId prefers focusJob over detailJob, so focusJob must always
    // reflect the job actually clicked — otherwise clicking a different job
    // in the nav list while another job is still "focused" (e.g. an earlier
    // running job whose progress view was opened) silently does nothing,
    // because the stale focusJob keeps winning the priority check below.
    if (["running", "queued", "waiting_checkpoint"].includes(j.status)) {
      setFocusJob(j.job_id);
    } else {
      setFocusJob(null); setFocusCpu(null); setFocusRam(null);
    }
    setDetailJob(j.job_id);
    loadDetail(j.job_id);
    setTimeout(() => {
      document.getElementById(`job-card-${j.job_id}`)
        ?.scrollIntoView({ behavior: "smooth", block: "center" });
    }, 50);
  };

  // ── Add files helper (accumulate, deduplicate by name) ───────────
  const addFiles = (incoming: FileList | File[]) => {
    const newArr = Array.from(incoming);
    setSelectedFiles(prev => {
      const merged = [...prev];
      for (const f of newArr) {
        if (!merged.some(m => m.name === f.name)) merged.push(f);
      }
      // Recompute sample names from non-ZIP files
      const isONT = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
      const nonZip = merged.filter(f => !f.name.toLowerCase().endsWith(".zip"));
      const r1 = isONT ? [] : nonZip.filter(f => /_R1|_1\.(fq|fastq)/i.test(f.name));
      const base = r1.length > 0 ? r1 : nonZip;
      const names = base.map(f =>
        isONT
          ? f.name.replace(/\.(fastq|fq)(\.gz)?$/i, "")
          : f.name.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "")
      );
      setSampleNames(names);
      setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
      return merged;
    });
  };

  // ── Upload ────────────────────────────────────────────────────────
  const handleUpload = async () => {
    if (!selectedFiles.length) { alert("Please select FASTQ or ZIP files first."); return; }
    setLoading(true);
    const form = new FormData();
    selectedFiles.forEach(f => form.append("files", f));
    try {
      const res = await axios.post(`${API}/upload`, form);
      setPendingJobId(res.data.job_id);
      // Use server-returned file list (ZIPs were extracted server-side)
      const serverFiles: string[] = res.data.files || [];
      setServerFileList(serverFiles);
      const isONT = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
      const r1 = isONT ? [] : serverFiles.filter(n => /_R1|_1\.(fq|fastq)/i.test(n));
      const base = r1.length > 0 ? r1 : serverFiles.filter(n => /\.(fastq|fq)(\.gz)?$/i.test(n));
      const names = base.map(n =>
        isONT
          ? n.replace(/\.(fastq|fq)(\.gz)?$/i, "")
          : n.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "")
      );
      if (names.length > 0) {
        setSampleNames(names);
        setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
      }
      setPairReadMode(isONT ? "single" : "paired");
      setFileMap(buildFileMapGuess(serverFiles));
    } catch { alert("Upload failed — is the backend running?"); }
    setLoading(false);
  };

  // ── Guess sample↔file pairs from filenames (starting point for manual editing) ──
  const buildFileMapGuess = (files: string[]): {sample:string; file1:string; file2:string}[] => {
    const fastq = files.filter(f => /\.(fastq|fq)(\.gz)?$/i.test(f)).sort();
    const used = new Set<string>();
    const rows: {sample:string; file1:string; file2:string}[] = [];
    // Single-end / long-read modes (ONT-16S marker, or QIIME2-VSEARCH sequencer type inside
    // the 16S marker card): every file is its own independent sample — never try to guess a
    // R1/R2 or _1/_2 mate, since "_1"/"_2" here means replicate number, not read-pair.
    const isSingleEnd = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
    for (const f of fastq) {
      if (used.has(f)) continue;
      if (isSingleEnd) {
        used.add(f);
        rows.push({ sample: f.replace(/\.(fastq|fq)(\.gz)?$/i, ""), file1: f, file2: "" });
        continue;
      }
      // Try to find this file's mate via common R1/R2 or _1/_2 conventions
      let mate = "";
      let sample = f.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "");
      if (/_R1|_1\.(fq|fastq)/i.test(f)) {
        const candidate = f.replace(/_R1/i, "_R2").replace(/_1\.(fq|fastq)/i, "_2.$1");
        if (fastq.includes(candidate)) mate = candidate;
      }
      if (mate) { used.add(f); used.add(mate); rows.push({ sample, file1: f, file2: mate }); }
      else {
        sample = f.replace(/\.(fastq|fq)(\.gz)?$/i, "");
        used.add(f);
        rows.push({ sample, file1: f, file2: "" });
      }
    }
    return rows;
  };

  // Re-derive sample↔file detection whenever the pairing mode changes (Marker or
  // Sequencer Type) for files that are already uploaded — e.g. user uploads files first,
  // then switches to "QIIME2 (VSEARCH OTU)" or ONT-16S afterward. Without this, the
  // "Detected samples" chips / manual-pairing table keep using whichever single-end vs
  // paired-end mode was active at the moment of upload and never update, even though the
  // single-end/paired-end hint banner above the drop-zone re-renders correctly.
  useEffect(() => {
    if (serverFileList.length === 0) return;
    const isONT = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
    setPairReadMode(isONT ? "single" : "paired");
    if (useManualPairing) {
      setFileMap(buildFileMapGuess(serverFileList));
    } else {
      const r1 = isONT ? [] : serverFileList.filter(n => /_R1|_1\.(fq|fastq)/i.test(n));
      const base = r1.length > 0 ? r1 : serverFileList.filter(n => /\.(fastq|fq)(\.gz)?$/i.test(n));
      const names = base.map(n =>
        isONT
          ? n.replace(/\.(fastq|fq)(\.gz)?$/i, "")
          : n.replace(/_R1.*|_1\.(fq|fastq).*/i, "").replace(/\.(fastq|fq)(\.gz)?$/i, "")
      );
      if (names.length > 0) setSampleNames(names);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [marker, params.sequencerType]);

  // ── Run ───────────────────────────────────────────────────────────
  const handleRun = async () => {
    if (!pendingJobId) return;

    // Validate manual pairing table before submit
    let sampleFileMap: {sample:string; file1:string; file2:string}[] = [];
    if (useManualPairing) {
      const rows = fileMap.filter(r => r.sample.trim() || r.file1);
      if (rows.length === 0) { alert("Add at least one sample row, or turn off manual pairing."); return; }
      for (const r of rows) {
        if (!r.sample.trim()) { alert("Every row needs a sample name."); return; }
        if (!r.file1)         { alert(`Sample "${r.sample}" is missing File 1.`); return; }
        if (pairReadMode === "paired" && !r.file2) {
          alert(`Sample "${r.sample}" is missing File 2 (required for Paired-end mode).`); return;
        }
      }
      const dupes = rows.map(r => r.sample).filter((s, i, a) => a.indexOf(s) !== i);
      if (dupes.length > 0) { alert(`Duplicate sample name(s): ${Array.from(new Set(dupes)).join(", ")}`); return; }
      sampleFileMap = rows.map(r => ({ sample: r.sample.trim(), file1: r.file1, file2: pairReadMode === "paired" ? r.file2 : "" }));
    }

    setLoading(true);
    try {
      await axios.post(`${API}/run/${pendingJobId}`, {
        job_name: jobName,
        marker,
        ...params,
        metadata,
        sampleFileMap,
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
    setJobName(""); setSelectedFiles([]); setPendingJobId(null);
    setSampleNames([]); setMetadata([]);
    setShowAdvanced(false); setShowSubmitPopup(false);
    setServerFileList([]); setUseManualPairing(false); setFileMap([]); setPairReadMode("paired");
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
  const handleCleanupOrphans = async () => {
    const res = await axios.delete(`${API}/results/cleanup`);
    setCleanupResult({ removed: res.data.removed, paths: res.data.paths });
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

  // ── Right sidebar: step checklist + run settings for the selected job ──
  const renderSidebarDetail = () => {
    if (!sidebarJobId) {
      return <div className="side-detail-empty">Select a job from the list to see its full<br/>status, step-by-step progress, and settings.</div>;
    }
    const job = allJobs.find(j => j.job_id === sidebarJobId);
    if (!job) return null;
    const d = detailData[sidebarJobId];
    const isActive = ["running", "queued", "waiting_checkpoint"].includes(job.status);
    const elapsed = isActive
      ? formatElapsed(job.started_at)
      : (d?.total_secs != null ? fmtSecs(d.total_secs) : "—");

    return (
      <>
        <div className="side-detail-header">
          <span className="job-status-icon">{STATUS_ICON[job.status]}</span>
          <div>
            <div className="side-detail-name">{job.job_name || "(unnamed)"}</div>
            <div className="side-detail-sub">{job.marker} · {job.job_id}</div>
          </div>
        </div>

        <div className="side-detail-time">
          <span>⏱️ {isActive ? "Running time" : "Total time"}</span>
          <strong>{elapsed}</strong>
        </div>

        <div className="side-detail-section-title">🧬 Pipeline Steps</div>
        {!d ? (
          <div className="detail-loading">⏳ Loading...</div>
        ) : d.step_history.length === 0 ? (
          <p className="detail-no-steps">
            {isActive ? "ℹ️ Starting up — steps will appear here shortly." : "ℹ️ Step breakdown not available — re-run this job to capture step details."}
          </p>
        ) : (
          <div className="side-step-list">
            {d.step_history.map((s, i) => {
              const barPct  = s.status === "running" ? 55 : 100;
              const barColor = s.status === "completed" ? "#10b981"
                             : s.status === "error"     ? "#ef4444" : "#3b82f6";
              const icon = s.status === "completed" ? "✅" : s.status === "error" ? "❌" : "🔄";
              return (
                <div key={i} className="side-step-item">
                  <div className="side-step-row">
                    <span className="side-step-icon">{icon}</span>
                    <span className="side-step-name">{s.label}</span>
                    <span className="side-step-dur">
                      {s.ended_at && s.started_at ? fmtSecs(s.ended_at - s.started_at) : "…"}
                    </span>
                  </div>
                  <div className="side-step-bar-track">
                    <div className={`side-step-bar-fill ${s.status === "running" ? "pulse" : ""}`}
                      style={{ width: `${barPct}%`, background: barColor }} />
                  </div>
                </div>
              );
            })}
          </div>
        )}

        <div className="side-detail-section-title">⚙️ Run Settings</div>
        <div className="side-settings-list">
          <div className="side-settings-row">
            <span>Files</span>
            <span>
              {job.files.length} file(s)
              {job.files.length > 4 && (
                <button className="side-files-toggle"
                  onClick={() => setSideFilesExpanded(v => !v)}>
                  {sideFilesExpanded ? "▲ hide" : "▼ show all"}
                </button>
              )}
            </span>
          </div>
          {job.files.length > 0 && (
            sideFilesExpanded ? (
              <div className="side-files-full-list">
                {job.files.map((f, i) => <div key={i} className="side-files-full-item">{i + 1}. {f}</div>)}
              </div>
            ) : (
              <div className="side-files-preview">
                {job.files.slice(0, 4).join(", ")}{job.files.length > 4 ? "…" : ""}
              </div>
            )
          )}
          <div className="side-settings-row">
            <span>Database</span>
            <span>{job.database || "—"}</span>
          </div>
          {d?.params && PARAM_KEY_ORDER.filter(k => {
            const v = d.params![k];
            if (v === undefined || v === null || v === "" || v === 0) return false;
            if (Array.isArray(v) && v.length === 0) return false;
            return isParamRelevantForMarker(k, d.params!.marker || job.marker, d.params!.sequencerType);
          }).map(k => (
            <div key={k} className="side-settings-row">
              <span>{PARAM_LABELS[k]}</span>
              <span>{typeof d.params![k] === "boolean" ? (d.params![k] ? "✅ Yes" : "✕ No") : String(d.params![k])}</span>
            </div>
          ))}
          {d?.params && Array.isArray(d.params.sampleFileMap) && d.params.sampleFileMap.length > 0 && (
            <div className="side-settings-row">
              <span>Manual Pairing</span>
              <span>{d.params.sampleFileMap.map((r: any) => r.sample).join(", ")}</span>
            </div>
          )}
        </div>

        {d?.error && (
          <div className="detail-error-section">
            <div className="detail-error-title">❌ Error log</div>
            <pre className="detail-error-body">{d.error}</pre>
          </div>
        )}
      </>
    );
  };

  // ── Job Card ──────────────────────────────────────────────────────
  const renderJobCard = (j: JobSummary) => {
    const displayName = j.job_name || "(unnamed)";

    return (
      <div key={j.job_id} id={`job-card-${j.job_id}`}
        className={`job-card ${focusJob === j.job_id ? "focused" : ""} ${sidebarJobId === j.job_id ? "sidebar-selected" : ""}`}
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
            {(j.status === "running" || j.status === "queued" || j.status === "waiting_checkpoint" || j.status === "uploaded") &&
              <button className="btn-cancel" onClick={() => handleCancel(j.job_id)}
                title={j.status === "uploaded" ? "Remove — never run" : "Cancel"}>✕</button>}
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
              <div
                className={`job-progress-fill${j.status === "running" ? " job-progress-fill--brand" : ""}`}
                style={{
                  width: `${j.progress}%`,
                  background: j.status === "running" ? undefined : STATUS_COLOR[j.status],
                }}
              />
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
                // Seed the re-run form with THIS job's own saved params
                // (primers, database, etc.) — without this, checkpointParams
                // was left at whatever it last happened to be (often still
                // the untouched hardcoded defaults), so "Re-run with New
                // Settings" silently dropped primer_f/primer_r/taxDatabase
                // and reset to defaults instead of the job's real settings.
                axios.get(`${API}/detail/${j.job_id}`)
                  .then(res => {
                    if (res.data?.params) {
                      setCheckpointParams({ ...defaultParams, ...res.data.params });
                    }
                  })
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
              <button className="btn-view btn-preview"
                onClick={() => { setPreviewJobId(j.job_id); setScreen("preview"); }}>
                ✏️ Edit Charts
              </button>
              <button className="btn-view"
                onClick={() => setShowColorPicker(showColorPicker === j.job_id ? null : j.job_id)}>
                🎨 Taxonomy Colors
              </button>
              <a href={`${API}/download/${j.job_id}`} download className="btn-view">
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
                  currentStep={focusLabel || j.step_label} theme={theme} />
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

                  {d.params && Object.keys(d.params).length > 0 && (
                    <div className="detail-step-list" style={{ marginTop: 12 }}>
                      <div className="detail-step-header" style={{ gridTemplateColumns: "1fr 1fr" }}>
                        <span>⚙️ Setting</span>
                        <span>Value</span>
                      </div>
                      {PARAM_KEY_ORDER.filter(k => {
                        const v = d.params![k];
                        if (v === undefined || v === null || v === "" || v === 0) return false;
                        if (Array.isArray(v) && v.length === 0) return false;
                        return isParamRelevantForMarker(k, d.params!.marker, d.params!.sequencerType);
                      }).map(k => (
                        <div key={k} className="detail-step-row" style={{ gridTemplateColumns: "1fr 1fr" }}>
                          <span className="detail-step-name">{PARAM_LABELS[k]}</span>
                          <span style={{ wordBreak: "break-all", fontSize: 12 }}>
                            {typeof d.params![k] === "boolean" ? (d.params![k] ? "✅ Yes" : "✕ No") : String(d.params![k])}
                          </span>
                        </div>
                      ))}
                      {Array.isArray(d.params!.sampleFileMap) && d.params!.sampleFileMap.length > 0 && (
                        <div className="detail-step-row" style={{ gridTemplateColumns: "1fr 1fr" }}>
                          <span className="detail-step-name">Manual Sample/File Pairing</span>
                          <span style={{ fontSize: 12 }}>
                            {d.params!.sampleFileMap.map((r: any) => r.sample).join(", ")}
                          </span>
                        </div>
                      )}
                    </div>
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
  // ══════════════════════════════════════════════════════════════
  //  PREVIEW SCREEN
  // ══════════════════════════════════════════════════════════════
  if (screen === "preview") {
    return (
      <PreviewPage
        initialJobId={previewJobId}
        onClose={() => setScreen("history")}
      />
    );
  }

  //  HOME SCREEN
  // ══════════════════════════════════════════════════════════════
  if (screen === "home") {
    return (
      <div className="home-screen">
        {/* Settings gear — top-right corner */}
        <button
          className="btn-settings"
          style={{ position: "absolute", top: 18, right: 22 }}
          title="Settings"
          onClick={() => setShowSettings(true)}
        >⚙️</button>

        {checkpointJobId && (
          <div className="chk-home-banner" onClick={() => setScreen("history")}>
            ⚠️ A pipeline has paused and needs your decision — <strong>click to review</strong>
          </div>
        )}
        <div className="home-logo-card">
          <span className="home-dna">🧬</span>
          <h1 className="home-title">16S / 12S / ONT Amplicon Analysis</h1>
          <p className="home-desc">
            DADA2 · Emu — Professional Microbiome &amp; Metabarcoding Analysis
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
            <button className="btn-cleanup-orphans" onClick={() => setCleanupConfirm(true)} title="ลบไฟล์ที่เหลือค้างบนดิสก์ (orphaned files)">🧹 Clean Disk</button>
            <button className="btn-refresh" onClick={refreshJobs}>↻ Refresh</button>
          </div>
        </div>
        <div className="hist-layout">
          <aside className="job-nav-sidebar">
            <div className="job-nav-title">All Jobs ({allJobs.length})</div>
            <div className="job-nav-list">
              {allJobs.length === 0 && <div className="job-nav-empty">No jobs yet</div>}
              {[...allJobs].reverse().map(j => {
                const dotClass =
                  ["running", "queued", "waiting_checkpoint"].includes(j.status) ? "active"
                  : j.status === "completed" ? "success"
                  : j.status === "error" ? "error" : "neutral";
                return (
                  <button key={j.job_id}
                    className={`job-nav-item ${sidebarJobId === j.job_id ? "selected" : ""}`}
                    onClick={() => selectJobFromNav(j)}
                    title={`${j.job_name || "(unnamed)"} — ${j.status}`}>
                    <span className={`job-nav-dot job-nav-dot-${dotClass}`} />
                    <span className="job-nav-name">{j.job_name || "(unnamed)"}</span>
                    <span className="job-nav-marker">{j.marker}</span>
                  </button>
                );
              })}
            </div>
          </aside>

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

          <aside className="job-side-detail">
            {renderSidebarDetail()}
          </aside>
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

        {/* ── Cleanup orphaned files confirmation ── */}
        {cleanupConfirm && (
          <div className="popup-overlay">
            <div className="step-popup del-confirm-popup">
              <div className="popup-step-icon">🧹</div>
              <h2 className="popup-title">ล้างไฟล์ค้างบนดิสก์</h2>
              <p className="del-confirm-msg">
                ลบโฟลเดอร์ผลลัพธ์บนดิสก์ที่<strong> ไม่มีรายการใน history</strong> แล้ว
                (orphaned files จาก job ที่ถูกลบไปแล้ว)
              </p>
              <p className="del-confirm-sub">⚠️ การดำเนินการนี้ไม่สามารถย้อนกลับได้</p>
              <div className="popup-footer">
                <button className="popup-back" onClick={() => setCleanupConfirm(false)}>← ยกเลิก</button>
                <button className="btn-del-confirm"
                  onClick={async () => {
                    await handleCleanupOrphans();
                    setCleanupConfirm(false);
                  }}>
                  🧹 ล้างไฟล์
                </button>
              </div>
            </div>
          </div>
        )}

        {/* ── Cleanup result ── */}
        {cleanupResult && (
          <div className="popup-overlay">
            <div className="step-popup del-confirm-popup">
              <div className="popup-step-icon">{cleanupResult.removed > 0 ? "✅" : "ℹ️"}</div>
              <h2 className="popup-title">ผลการล้างไฟล์</h2>
              {cleanupResult.removed === 0 ? (
                <p className="del-confirm-msg">ไม่พบไฟล์ค้างบนดิสก์ — ดิสก์สะอาดดีอยู่แล้ว 👍</p>
              ) : (
                <>
                  <p className="del-confirm-msg">
                    ลบไฟล์ค้าง <strong>{cleanupResult.removed} โฟลเดอร์</strong> เรียบร้อยแล้ว
                  </p>
                  <div className="cleanup-path-list">
                    {cleanupResult.paths.map(p => (
                      <div key={p} className="cleanup-path-item">📁 {p}</div>
                    ))}
                  </div>
                </>
              )}
              <div className="popup-footer">
                <button className="btn-del-confirm" onClick={() => setCleanupResult(null)}>ตกลง</button>
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

      {showSettings && (
        <SettingsPanel
          theme={theme}
          onTheme={(t) => setTheme(t)}
          onClose={() => setShowSettings(false)}
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
                <strong>{selectedFiles.length} file(s)</strong>
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
                <input
                  type="number"
                  className={`worker-custom-input ${![1,2,4,6,8].includes(maxWorkers) ? "active" : ""}`}
                  min={1} max={256}
                  value={maxWorkers}
                  title="กำหนดจำนวน CPU เอง"
                  placeholder="เอง"
                  onChange={async (e) => {
                    const v = Math.max(1, Math.min(256, parseInt(e.target.value) || 1));
                    setMaxWorkers(v);
                    await axios.post(`${API}/config`, { max_workers: v });
                  }}
                />
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
            <button className="btn-settings" title="Settings" onClick={() => setShowSettings(true)}>
              ⚙️
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

        {/* Upload */}
        <div className="section-card">
          <div className="section-title">📂 Upload FASTQ Files</div>
          {(marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch")) && (
            <div className="ps-info-box ps-info-box--ok" style={{ marginBottom: 10, fontSize: 13 }}>
              🧫 <strong>Single-end mode:</strong> Upload one FASTQ file per sample (single reads, not paired R1/R2).
              Files like <code>AC_1.fastq.gz</code> / <code>AC_2.fastq.gz</code> are treated as two separate
              samples, not an R1/R2 pair.
            </div>
          )}
          <div
            className={`drop-zone ${dragOver ? "drag-over" : ""}`}
            onDragOver={e => { e.preventDefault(); setDragOver(true); }}
            onDragLeave={() => setDragOver(false)}
            onDrop={e => { e.preventDefault(); setDragOver(false); if (e.dataTransfer.files.length) addFiles(e.dataTransfer.files); }}
            onClick={() => fileInputRef.current?.click()}
          >
            {selectedFiles.length
              ? <>➕ <strong>{selectedFiles.length} file(s)</strong> selected — click/drop to add more</>
              : (marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch"))
                ? <>🧫 Drop single-end .fastq / .fastq.gz / .zip here, or <u>click to browse</u></>
                : <>📂 Drop .fastq / .fastq.gz / .zip here, or <u>click to browse</u></>
            }
          </div>
          <input
            ref={fileInputRef} type="file" multiple
            accept=".fastq,.fastq.gz,.fq,.fq.gz,.zip"
            style={{ display: "none" }}
            onChange={e => { if (e.target.files?.length) { addFiles(e.target.files); e.target.value = ""; } }}
          />
          {selectedFiles.length > 0 && (
            <div className="sample-list" style={{ marginTop: 8 }}>
              <div className="sample-list-title" style={{ display:"flex", justifyContent:"space-between", alignItems:"center" }}>
                <span>📁 Selected files ({selectedFiles.length})</span>
                <button className="prev-btn-reset" style={{ fontSize:11, padding:"2px 8px" }}
                  onClick={() => { setSelectedFiles([]); setSampleNames([]); setMetadata([]); }}>
                  ✕ Clear all
                </button>
              </div>
              <div className="sample-chips" style={{ marginTop: 4 }}>
                {selectedFiles.map(f => (
                  <span key={f.name} className="sample-chip" style={{ display:"inline-flex", alignItems:"center", gap:4 }}>
                    {f.name.toLowerCase().endsWith(".zip") ? "📦 " : ""}{f.name}
                    <button style={{ background:"none", border:"none", cursor:"pointer", color:"#94a3b8", fontSize:12, padding:0, lineHeight:1 }}
                      onClick={ev => { ev.stopPropagation();
                        setSelectedFiles(prev => {
                          const next = prev.filter(x => x.name !== f.name);
                          const isONT = marker === "ONT-16S" || (marker === "16S" && params.sequencerType === "qiime2_vsearch");
                          const nonZip = next.filter(x => !x.name.toLowerCase().endsWith(".zip"));
                          const r1 = isONT ? [] : nonZip.filter(x => /_R1|_1\.(fq|fastq)/i.test(x.name));
                          const base = r1.length > 0 ? r1 : nonZip;
                          const names = base.map(x => isONT ? x.name.replace(/\.(fastq|fq)(\.gz)?$/i,"") : x.name.replace(/_R1.*|_1\.(fq|fastq).*/i,"").replace(/\.(fastq|fq)(\.gz)?$/i,""));
                          setSampleNames(names);
                          setMetadata(names.map(s => ({ sampleId: s, group: "", description: "" })));
                          return next;
                        });
                      }}>✕</button>
                  </span>
                ))}
              </div>
            </div>
          )}

          {sampleNames.length > 0 && (
            <div className="sample-list">
              <div className="sample-list-title">🧪 Detected samples ({sampleNames.length})</div>
              <div className="sample-chips">
                {sampleNames.map(s => <span key={s} className="sample-chip">{s}</span>)}
              </div>
            </div>
          )}

          {/* ── Manual sample↔file pairing (advanced) ─────────────────────── */}
          {pendingJobId && serverFileList.length > 0 && (
            <div style={{ marginTop: 14, borderTop: "1px solid #334155", paddingTop: 12 }}>
              <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer", fontSize: 13 }}>
                <input type="checkbox" checked={useManualPairing}
                  onChange={e => setUseManualPairing(e.target.checked)} />
                ✏️ <strong>Manually assign files per sample</strong>
                <span style={{ color: "#94a3b8", fontWeight: 400 }}>
                  (fixes auto-detect mistakes — e.g. ONT vs paired-end mix-ups)
                </span>
              </label>

              {useManualPairing && (
                <div style={{ marginTop: 10 }}>
                  <div style={{ display: "flex", gap: 8, marginBottom: 10 }}>
                    <button type="button"
                      style={{
                        padding: "5px 14px", borderRadius: 6, fontSize: 12, cursor: "pointer",
                        border: pairReadMode === "single" ? "2px solid #f59e0b" : "1px solid #334155",
                        background: pairReadMode === "single" ? "#1c1200" : "#1e293b",
                        color: pairReadMode === "single" ? "#fcd34d" : "#94a3b8",
                        fontWeight: pairReadMode === "single" ? 600 : 400,
                      }}
                      onClick={() => setPairReadMode("single")}>
                      🧬 Single-end / long-read (1 file/sample)
                    </button>
                    <button type="button"
                      style={{
                        padding: "5px 14px", borderRadius: 6, fontSize: 12, cursor: "pointer",
                        border: pairReadMode === "paired" ? "2px solid #6366f1" : "1px solid #334155",
                        background: pairReadMode === "paired" ? "#1e1b4b" : "#1e293b",
                        color: pairReadMode === "paired" ? "#a5b4fc" : "#94a3b8",
                        fontWeight: pairReadMode === "paired" ? 600 : 400,
                      }}
                      onClick={() => setPairReadMode("paired")}>
                      🔗 Paired-end (merge 2 files/sample)
                    </button>
                    <button type="button" className="prev-btn-reset" style={{ fontSize: 12 }}
                      onClick={() => setFileMap(buildFileMapGuess(serverFileList))}>
                      ↺ Re-guess from filenames
                    </button>
                  </div>

                  <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                    <div style={{ display: "grid",
                      gridTemplateColumns: pairReadMode === "paired" ? "1fr 1.4fr 1.4fr 28px" : "1fr 1.4fr 28px",
                      gap: 6, fontSize: 11, color: "#94a3b8", fontWeight: 600, padding: "0 4px" }}>
                      <span>Sample Name</span>
                      <span>File 1{pairReadMode === "single" ? "" : " (Forward)"}</span>
                      {pairReadMode === "paired" && <span>File 2 (Reverse)</span>}
                      <span></span>
                    </div>
                    {fileMap.map((row, i) => (
                      <div key={i} style={{ display: "grid",
                        gridTemplateColumns: pairReadMode === "paired" ? "1fr 1.4fr 1.4fr 28px" : "1fr 1.4fr 28px",
                        gap: 6, alignItems: "center" }}>
                        <input value={row.sample}
                          onChange={e => setFileMap(prev => prev.map((r, j) => j === i ? { ...r, sample: e.target.value } : r))}
                          style={{ padding: "5px 8px", borderRadius: 5, border: "1px solid #334155",
                                   background: "#0f172a", color: "#e2e8f0", fontSize: 12 }} />
                        <select value={row.file1}
                          onChange={e => setFileMap(prev => prev.map((r, j) => j === i ? { ...r, file1: e.target.value } : r))}
                          style={{ padding: "5px 8px", borderRadius: 5, border: "1px solid #334155",
                                   background: "#0f172a", color: "#e2e8f0", fontSize: 12 }}>
                          <option value="">— select file —</option>
                          {serverFileList.map(f => <option key={f} value={f}>{f}</option>)}
                        </select>
                        {pairReadMode === "paired" && (
                          <select value={row.file2}
                            onChange={e => setFileMap(prev => prev.map((r, j) => j === i ? { ...r, file2: e.target.value } : r))}
                            style={{ padding: "5px 8px", borderRadius: 5, border: "1px solid #334155",
                                     background: "#0f172a", color: "#e2e8f0", fontSize: 12 }}>
                            <option value="">— select file —</option>
                            {serverFileList.map(f => <option key={f} value={f}>{f}</option>)}
                          </select>
                        )}
                        <button type="button" title="Remove row"
                          onClick={() => setFileMap(prev => prev.filter((_, j) => j !== i))}
                          style={{ background: "none", border: "none", color: "#94a3b8", cursor: "pointer", fontSize: 14 }}>
                          ✕
                        </button>
                      </div>
                    ))}
                  </div>

                  <button type="button" className="prev-btn-reset" style={{ marginTop: 8, fontSize: 12 }}
                    onClick={() => setFileMap(prev => [...prev, { sample: "", file1: "", file2: "" }])}>
                    + Add sample row
                  </button>

                  <div style={{ marginTop: 8, fontSize: 11.5, color: "#94a3b8", lineHeight: 1.5 }}>
                    💡 This overrides filename auto-detection. {pairReadMode === "paired"
                      ? "Every row needs both File 1 and File 2 — they'll be merged with DADA2's paired-end workflow."
                      : "Only File 1 is needed per row — reads are analyzed as single/long reads (no merge step)."}
                  </div>
                </div>
              )}
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
