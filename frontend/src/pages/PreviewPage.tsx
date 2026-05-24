// ============================================================
//  Result Preview Page — Interactive Chart Viewer
//  Left 20%: job list  |  Right 80%: chart + customize
// ============================================================
import { useState, useEffect, useRef, useCallback } from "react";
import "./PreviewPage.css";

const API = "http://localhost:8000";
declare const Plotly: any;

// ── Default color palette ─────────────────────────────────────
const DEFAULT_COLORS = [
  "#3b82f6","#ef4444","#10b981","#f59e0b","#8b5cf6",
  "#06b6d4","#f97316","#84cc16","#ec4899","#14b8a6",
  "#6366f1","#a855f7","#22d3ee","#fb923c","#a3e635",
  "#f472b6","#2dd4bf","#818cf8","#fb7185","#fbbf24",
  "#0ea5e9","#d946ef","#65a30d","#dc2626","#7c3aed",
];

// ── Chart tab definitions ─────────────────────────────────────
interface TabDef {
  id: string; label: string; file: string;
  type: "taxonomy"|"alpha"|"multialpha"|"bar"|"line"|"scatter"|"box"|"heatmap"|"taxheatmap"|"scree"|"longline"|"specaccum"|"pdf"|"readtrack";
  metric?: string; group?: string;
  pdfFile?: string;   // for type:"pdf" — filename inside r_plots/ or job root
  altFile?: string;   // fallback CSV filename when primary doesn't exist
}
const ALL_TABS: TabDef[] = [
  // ── Taxonomy ───────────────────────────────────────────────
  { id:"phylum",      label:"Phylum",           file:"taxonomy_phylum.csv",    type:"taxonomy", group:"Taxonomy" },
  { id:"class",       label:"Class",            file:"taxonomy_class.csv",     type:"taxonomy", group:"Taxonomy" },
  { id:"order",       label:"Order",            file:"taxonomy_order.csv",     type:"taxonomy", group:"Taxonomy" },
  { id:"family",      label:"Family",           file:"taxonomy_family.csv",    type:"taxonomy", group:"Taxonomy" },
  { id:"genus",       label:"Genus",            file:"taxonomy_genus.csv",     type:"taxonomy", group:"Taxonomy" },
  { id:"species",     label:"Species",          file:"taxonomy_species.csv",   type:"taxonomy", group:"Taxonomy" },
  // Taxonomy heatmaps — interactive (sample × taxon abundance heatmap)
  { id:"tax_hm_phy",  label:"Tax. Heatmap Phylum", file:"taxonomy_phylum.csv", type:"taxheatmap", group:"Taxonomy" },
  { id:"tax_hm_fam",  label:"Tax. Heatmap Family", file:"taxonomy_family.csv", type:"taxheatmap", group:"Taxonomy" },
  { id:"tax_hm_gen",  label:"Tax. Heatmap Genus",  file:"taxonomy_genus.csv",  type:"taxheatmap", group:"Taxonomy" },
  // ── Alpha diversity ────────────────────────────────────────
  { id:"shannon",     label:"Shannon",          file:"alpha_diversity.csv",    type:"alpha", metric:"Shannon",  group:"Alpha" },
  { id:"observed",    label:"Observed",         file:"alpha_diversity.csv",    type:"alpha", metric:"Observed", group:"Alpha" },
  { id:"chao1",       label:"Chao1",            file:"alpha_diversity.csv",    type:"alpha", metric:"Chao1",    group:"Alpha" },
  { id:"simpson",     label:"Simpson",          file:"alpha_diversity.csv",    type:"alpha", metric:"Simpson",  group:"Alpha" },
  { id:"faiths_pd",   label:"Faith's PD",       file:"faiths_pd.csv",         type:"alpha", metric:"PD",       group:"Alpha" },
  { id:"shan_rar",    label:"Shannon Rar.",      file:"shannon_rarefaction.csv",type:"longline",                group:"Alpha" },
  // All alpha metrics together — interactive 2×2 subplot
  { id:"multialpha",  label:"All Alpha Metrics", file:"alpha_diversity.csv", type:"multialpha", group:"Alpha" },
  // PDFs — alpha fallback (DADA2 root or r_plots/)
  { id:"pdf_obs_pdf", label:"Observed (PDF)",    file:"_pdf", type:"pdf", pdfFile:"observed_asvs.pdf", group:"Alpha" },
  // ── Beta diversity ─────────────────────────────────────────
  { id:"pca",         label:"PCoA",             file:"pca_scores.csv",         type:"scatter", group:"Beta" },
  { id:"pca_scree",   label:"PCA Scree",        file:"pca_scree.csv",          type:"scree",   group:"Beta" },
  { id:"nmds_bray",   label:"NMDS Bray",        file:"nmds_bray.csv",          type:"scatter", group:"Beta" },
  { id:"nmds_jac",    label:"NMDS Jaccard",     file:"nmds_jaccard.csv",       type:"scatter", group:"Beta" },
  // Heatmap: viz_pipeline → beta_heatmap.csv (similarity); DADA2 → bray_curtis_distance_matrix.csv (distance)
  { id:"heatmap",     label:"Heatmap (Bray)",   file:"beta_heatmap.csv",       type:"heatmap", altFile:"bray_curtis_distance_matrix.csv", group:"Beta" },
  { id:"jac_heatmap", label:"Heatmap (Jac)",    file:"jaccard_heatmap.csv",    type:"heatmap", group:"Beta" },
  // PDFs — beta (viz_pipeline r_plots/ names)
  { id:"upgma_bray",  label:"UPGMA Bray",       file:"_pdf", type:"pdf", pdfFile:"05_beta_UPGMA.pdf",           group:"Beta" },
  { id:"upgma_jac",   label:"UPGMA Jaccard",    file:"_pdf", type:"pdf", pdfFile:"12_upgma_jaccard.pdf",        group:"Beta" },
  // PDFs — beta (DADA2 root names)
  { id:"pdf_pcoa",    label:"PCoA (PDF)",        file:"_pdf", type:"pdf", pdfFile:"beta_pcoa.pdf",              group:"Beta" },
  { id:"pdf_upgma_d", label:"UPGMA (PDF)",       file:"_pdf", type:"pdf", pdfFile:"beta_upgma.pdf",            group:"Beta" },
  { id:"pdf_hm_bray", label:"Heatmap Bray (PDF)",file:"_pdf", type:"pdf", pdfFile:"beta_heatmap.pdf",          group:"Beta" },
  { id:"pdf_hm_jac",  label:"Heatmap Jac (PDF)", file:"_pdf", type:"pdf", pdfFile:"beta_heatmap_jaccard.pdf",  group:"Beta" },
  { id:"pdf_nmds_jac",label:"NMDS Jac (PDF)",    file:"_pdf", type:"pdf", pdfFile:"beta_nmds_jaccard.pdf",     group:"Beta" },
  // ── Other ──────────────────────────────────────────────────
  { id:"reads",       label:"Read Counts",      file:"asv_summary.csv",        type:"bar",       group:"Other" },
  { id:"read_track",  label:"Read Tracking",    file:"read_tracking.csv",      type:"readtrack", group:"Other" },
  { id:"rarefaction", label:"Rarefaction",      file:"rarefaction.csv",        type:"line",      group:"Other" },
  { id:"specaccum",   label:"Spec. Accum.",     file:"specaccum.csv",          type:"specaccum", group:"Other" },
  { id:"rank_abund",  label:"Rank Abund.",      file:"rank_abundance.csv",     type:"longline",  group:"Other" },
  { id:"otu",         label:"OTU Dist",         file:"otu_distribution.csv",   type:"box",       group:"Other" },
  // PDFs — other (DADA2 root)
  { id:"pdf_rar",     label:"Rarefaction (PDF)",file:"_pdf", type:"pdf", pdfFile:"rarefaction_curves.pdf",      group:"Other" },
  { id:"pdf_asv_len", label:"ASV Length (PDF)", file:"_pdf", type:"pdf", pdfFile:"asv_length_distribution.pdf", group:"Other" },
  { id:"pdf_prev",    label:"Prevalence (PDF)", file:"_pdf", type:"pdf", pdfFile:"prevalence_abundance.pdf",    group:"Other" },
  { id:"pdf_qc",      label:"QC Read Count (PDF)",file:"_pdf",type:"pdf", pdfFile:"qc_readcount_boxplot.pdf",  group:"Other" },
  { id:"pdf_rtrack",  label:"Read Tracking (PDF)",file:"_pdf",type:"pdf", pdfFile:"read_tracking_plot.pdf",    group:"Other" },
];

// ── Font config ───────────────────────────────────────────────
interface FontConfig {
  titleSize: number; axisSize: number; legendSize: number;
  titleText: string; titleBold: boolean; titleItalic: boolean;
  chartType: "bar"|"bar100"|"line"; showGrid: boolean;
  xItalic: boolean; xBold: boolean; xSize: number;
  colorscale: string;  // for heatmap/taxheatmap
}
const defaultFont = (): FontConfig => ({
  titleSize: 16, axisSize: 12, legendSize: 11,
  titleText: "", titleBold: false, titleItalic: false,
  chartType: "bar", showGrid: true,
  xItalic: false, xBold: false, xSize: 12,
  colorscale: "Viridis",
});

// ── Colorscale options for heatmap tabs ──────────────────────
const COLORSCALES = [
  { name:"Viridis", gradient:"linear-gradient(to right,#440154,#31688e,#35b779,#fde725)" },
  { name:"Plasma",  gradient:"linear-gradient(to right,#0d0887,#7e03a8,#f89540,#f0f921)" },
  { name:"Blues",   gradient:"linear-gradient(to right,#f7fbff,#6baed6,#08306b)" },
  { name:"Reds",    gradient:"linear-gradient(to right,#fff5f0,#fc8d59,#67000d)" },
  { name:"YlOrRd",  gradient:"linear-gradient(to right,#ffffcc,#fd8d3c,#800026)" },
  { name:"RdYlGn",  gradient:"linear-gradient(to right,#d73027,#ffffbf,#1a9850)" },
  { name:"Hot",     gradient:"linear-gradient(to right,#000,#f00,#ff0,#fff)" },
  { name:"Greens",  gradient:"linear-gradient(to right,#f7fcf5,#74c476,#00441b)" },
];

// ── CSV parser ────────────────────────────────────────────────
function parseCSV(text: string): string[][] {
  const rows: string[][] = [];
  for (const line of text.trim().split("\n")) {
    const row: string[] = [];
    let cur = "", inQ = false;
    for (const ch of line) {
      if (ch === '"') { inQ = !inQ; }
      else if (ch === "," && !inQ) { row.push(cur.trim()); cur = ""; }
      else cur += ch;
    }
    row.push(cur.trim());
    rows.push(row);
  }
  return rows;
}
const num = (v: string) => { const n = parseFloat(v); return isNaN(n) ? 0 : n; };
const shortName = (s: string) => s.replace(/FBE\d+_pass_/i,"").replace(/_\w{8}_\w{8}_\d+/,"").slice(0,20);

// ── Get raw sample names for x-axis (used in customize panel) ─
function getUniqueSamples(tab: TabDef, rows: string[][]): string[] {
  if (!rows.length) return [];
  if (tab.type === "taxonomy" || tab.type === "taxheatmap") return rows.slice(1).map(r => r[0]);
  if (tab.type === "alpha" || tab.type === "multialpha") {
    const sIdx = rows[0].indexOf("Sample");
    return sIdx >= 0 ? rows.slice(1).map(r => r[sIdx]) : rows.slice(1).map(r => r[r.length - 1]);
  }
  if (tab.type === "bar" || tab.type === "box") return rows.slice(1).map(r => r[0]);
  return [];
}

// ── Interfaces ────────────────────────────────────────────────
interface Job { job_id: string; job_name: string; status: string; marker: string; }
interface TableInfo { tables: string[]; plots: string[]; summary: any; }

// ── Props ─────────────────────────────────────────────────────
interface PreviewPageProps { initialJobId?: string; onClose: () => void; }

export default function PreviewPage({ initialJobId, onClose }: PreviewPageProps) {
  const [jobs, setJobs]           = useState<Job[]>([]);
  const [search, setSearch]       = useState("");
  const [selJob, setSelJob]       = useState<string>(initialJobId || "");
  const [tableInfo, setTableInfo] = useState<TableInfo | null>(null);
  const [availTabs, setAvailTabs] = useState<TabDef[]>([]);
  const [activeTab, setActiveTab] = useState<TabDef | null>(null);
  const [csvData, setCsvData]     = useState<string[][] | null>(null);
  const [colors, setColors]       = useState<Record<string, string>>({});
  const [font, setFont]           = useState<FontConfig>(defaultFont());
  const [sampleAliases, setSampleAliases] = useState<Record<string, string>>({});
  const [knownSamples, setKnownSamples]   = useState<string[]>([]);  // persists across PDF tabs
  const [showCustomize, setShowCustomize] = useState(true);
  const [legendPicker, setLegendPicker]   = useState<{name:string; x:number; y:number; color:string} | null>(null);
  const lastMousePos = useRef({ x: 100, y: 100 });
  const [loading, setLoading]     = useState(false);
  const [pdfPanel, setPdfPanel]   = useState(false);
  const [saved, setSaved]         = useState(false);
  const [rerunning, setRerunning] = useState(false);
  const [exporting, setExporting] = useState(false);
  const [exportStep, setExportStep] = useState("");
  const chartRef = useRef<HTMLDivElement>(null);
  const plotted  = useRef(false);

  // ── Load jobs ─────────────────────────────────────────────
  useEffect(() => {
    axios_get<{jobs: Job[]}>(`${API}/jobs`)
      .then(d => setJobs(d.jobs.filter(j => j.status === "completed")));
  }, []);

  // ── Load tables when job changes ──────────────────────────
  useEffect(() => {
    if (!selJob) return;
    setTableInfo(null); setAvailTabs([]); setActiveTab(null); setCsvData(null);

    // Load settings + tables in parallel; apply settings before activating first tab
    const loadSettings = async (): Promise<boolean> => {
      // Try server first — check edit_charts/settings then legacy preview_settings
      try {
        const s = await fetch(`${API}/results/${selJob}/edit_charts/settings`).then(r => r.json());
        if (s && Object.keys(s).length > 0) {
          if (s.font)          setFont({ ...defaultFont(), ...s.font });
          if (s.colors)        setColors(s.colors);
          if (s.sampleAliases) setSampleAliases(s.sampleAliases);
          return true;
        }
      } catch { /* fall through */ }
      // Fall back to localStorage
      try {
        const local = localStorage.getItem(`prev_settings_${selJob}`);
        if (local) {
          const s = JSON.parse(local);
          if (s.font)          setFont({ ...defaultFont(), ...s.font });
          if (s.colors)        setColors(s.colors);
          if (s.sampleAliases) setSampleAliases(s.sampleAliases);
          return true;
        }
      } catch { /* ignore */ }
      // No saved settings
      setFont(defaultFont()); setColors({}); setSampleAliases({});
      return false;
    };

    const loadTables = async () => {
      const info: TableInfo = await fetch(`${API}/results/${selJob}/tables`).then(r => r.json());
      setTableInfo(info);
      const avail = ALL_TABS.filter(t => {
        if (t.type === "pdf") return (info.plots || []).includes(t.pdfFile || "");
        return info.tables.includes(t.file) || (!!t.altFile && info.tables.includes(t.altFile));
      });
      setAvailTabs(avail);
      return { info, avail };
    };

    // Run both in parallel, then activate the first tab
    Promise.all([loadSettings(), loadTables()])
      .then(([hadSettings, { info, avail }]) => {
        if (avail.length > 0) setActiveTab(avail[0]);
        if (info.summary?.job_name && !hadSettings) {
          setFont(f => ({ ...f, titleText: info.summary.job_name }));
        }
      })
      .catch(() => {
        fetch(`${API}/results/${selJob}/tables`).then(r => r.json()).then((info: TableInfo) => {
          setTableInfo(info);
          const avail = ALL_TABS.filter(t => info.tables.includes(t.file) || (!!t.altFile && info.tables.includes(t.altFile)));
          setAvailTabs(avail);
          if (avail.length > 0) setActiveTab(avail[0]);
        }).catch(() => {});
      });
  }, [selJob]);

  // ── Load CSV when tab changes ─────────────────────────────
  useEffect(() => {
    if (!activeTab || !selJob) return;
    if (activeTab.type === "pdf") { setCsvData(null); setLoading(false); return; }
    setLoading(true);
    // Use altFile when primary doesn't exist in tableInfo
    const fileToLoad = (tableInfo?.tables || []).includes(activeTab.file)
      ? activeTab.file
      : (activeTab.altFile || activeTab.file);
    fetch(`${API}/results/${selJob}/table/${fileToLoad}`)
      .then(r => r.text())
      .then(text => {
        const rows = parseCSV(text);
        setCsvData(rows);
        // auto-assign colors to new series
        setColors(prev => {
          const next = { ...prev };
          const series = getSeries(activeTab, rows);
          series.forEach((name, i) => {
            if (!next[name]) next[name] = DEFAULT_COLORS[i % DEFAULT_COLORS.length];
          });
          return next;
        });
        // Persist sample names for use when viewing PDF tabs
        const samplesFromRows = getUniqueSamples(activeTab, rows);
        if (samplesFromRows.length > 0) setKnownSamples(samplesFromRows);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, [activeTab, selJob]);

  // ── Re-render chart when data / colors / font / aliases changes ────
  // (renderChart is declared below but useCallback deps already cover it)
  const triggerRender = useCallback(() => {
    if (!csvData || !chartRef.current || !activeTab) return;
    if (activeTab.type === "pdf") return; // PDF tabs render inline, not via Plotly
    if (!window.Plotly) {
      // Plotly CDN not ready yet — retry in 500 ms
      const el = chartRef.current;
      setTimeout(() => {
        if (window.Plotly && el) {
          const p2 = buildPlotly();
          if (p2) { Plotly.newPlot(el, p2.data, p2.layout, { responsive: true }); plotted.current = true; }
        }
      }, 500);
      return;
    }
    const p = buildPlotly();
    if (!p) return;
    if (plotted.current) {
      Plotly.react(chartRef.current, p.data, p.layout, { responsive: true });
    } else {
      Plotly.newPlot(chartRef.current, p.data, p.layout, { responsive: true });
      plotted.current = true;
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [csvData, colors, font, activeTab, sampleAliases]);

  useEffect(() => { triggerRender(); }, [triggerRender]);

  // ── Close legend picker when clicking outside (delayed to avoid self-close) ─
  useEffect(() => {
    if (!legendPicker) return;
    let cleanup: (() => void) | undefined;
    // 200ms delay: legendclick fires a window-click too — don't catch it
    const t = setTimeout(() => {
      const close = (e: MouseEvent) => {
        const target = e.target as HTMLElement;
        if (target.closest?.(".prev-legend-picker")) return; // click inside → keep open
        setLegendPicker(null);
      };
      window.addEventListener("click", close);
      cleanup = () => window.removeEventListener("click", close);
    }, 200);
    return () => { clearTimeout(t); cleanup?.(); };
  }, [legendPicker]);

  // ── Track mouse position for legend picker placement ─────────
  useEffect(() => {
    const track = (e: MouseEvent) => { lastMousePos.current = { x: e.clientX, y: e.clientY }; };
    window.addEventListener("mousemove", track);
    return () => window.removeEventListener("mousemove", track);
  }, []);

  // ── Plotly legend click → floating color picker ────────────
  useEffect(() => {
    const el = chartRef.current as any;
    if (!el || !activeTab || activeTab.type === "pdf") return;

    const handler = (data: any) => {
      const traceName = data.data?.[data.curveNumber]?.name;
      if (!traceName) return;
      // Try data.event first; fall back to tracked mouse position
      const ex = (data.event as any)?.clientX ?? lastMousePos.current.x;
      const ey = (data.event as any)?.clientY ?? lastMousePos.current.y;
      const fallback = DEFAULT_COLORS[data.curveNumber % DEFAULT_COLORS.length];
      setLegendPicker({ name: traceName, x: ex, y: ey, color: colors[traceName] || fallback });
      return false; // prevents Plotly from toggling trace visibility
    };
    el.on?.("plotly_legendclick", handler);
    return () => { try { el.removeAllListeners?.("plotly_legendclick"); } catch {} };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [csvData, activeTab, colors]);

  // ── Get series names for a tab ────────────────────────────
  function getSeries(tab: TabDef, rows: string[][]): string[] {
    if (!rows.length) return [];
    if (tab.type === "taxonomy")    return rows[0].slice(1).slice(0, 30);  // taxon names → color assignment
    if (tab.type === "taxheatmap") return rows.slice(1).map(r => r[0]);   // sample names → for color (unused but consistent)
    if (tab.type === "readtrack")  return rows.slice(1).map(r => r[0]);   // sample names → color per line
    if (tab.type === "alpha" || tab.type === "multialpha") return rows.length > 1 ? rows.slice(1).map(r => {
      const idx = rows[0].indexOf("Sample"); return idx >= 0 ? r[idx] : r[r.length-1];
    }) : [];
    if (tab.type === "line")      return rows[0]?.slice(1) || [];
    if (tab.type === "scatter")   return rows.slice(1).map(r => r[0]) || [];
    if (tab.type === "bar")       return rows.slice(1).map(r => shortName(r[0]));
    if (tab.type === "longline")  {
      const idx = rows[0].indexOf("Sample");
      return idx >= 0 ? Array.from(new Set(rows.slice(1).map(r => r[idx]))) : [];
    }
    return [];
  }

  // ── Build Plotly traces & layout ──────────────────────────
  // Accepts optional overrides so exportResults can render any tab without switching state
  function buildPlotly(overrideData?: string[][], overrideTab?: TabDef): { data: any[], layout: any } | null {
    const rows = overrideData ?? csvData;
    const activeTab2 = overrideTab ?? activeTab;
    if (!rows || !activeTab2) return null;
    const title = font.titleText || activeTab2.label;
    const titleStr = (font.titleBold ? "<b>" : "") +
                     (font.titleItalic ? "<i>" : "") +
                     title +
                     (font.titleItalic ? "</i>" : "") +
                     (font.titleBold ? "</b>" : "");
    // alias: check raw name first, then shortName version
    const alias = (rawName: string) =>
      sampleAliases[rawName] || sampleAliases[shortName(rawName)] || shortName(rawName);
    // Format x-axis tick label with bold/italic HTML (Plotly renders HTML in ticktext)
    const xFmt = (name: string) => {
      let t = name;
      if (font.xBold)   t = `<b>${t}</b>`;
      if (font.xItalic) t = `<i>${t}</i>`;
      return t;
    };
    const base = {
      title: { text: titleStr, font: { size: font.titleSize, family: "Arial" } },
      xaxis: {
        tickfont: { size: font.xSize },
        showgrid: font.showGrid, tickangle: -35,
      },
      yaxis: { tickfont: { size: font.axisSize }, showgrid: font.showGrid },
      legend: { font: { size: font.legendSize }, orientation: "v" as const },
      paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff",
      margin: { l: 60, r: 20, t: 60, b: 120 },
    };

    // TAXONOMY — stacked bar: X = samples, series = taxa (phylum/genus/etc.)
    if (activeTab2.type === "taxonomy" && rows.length > 1) {
      const phylumNames = rows[0].slice(1);          // column headers = taxa names
      const dataRows    = rows.slice(1);             // each row = one sample
      const sampleNames = dataRows.map(r => alias(r[0]));
      const is100 = font.chartType === "bar100";
      // normalize per row (sample) so each sample sums to 100%
      const rowSums = dataRows.map(r => r.slice(1).reduce((s, v) => s + num(v), 0));
      const data = phylumNames.map((taxon, ci) => ({
        name: taxon || "Unknown",
        x: sampleNames,
        y: dataRows.map((row, ri) => {
          const v = num(row[ci + 1]);
          return is100 && rowSums[ri] > 0 ? (v / rowSums[ri]) * 100 : v;
        }),
        type: "bar",
        marker: { color: colors[taxon] || DEFAULT_COLORS[ci % DEFAULT_COLORS.length] },
      }));
      return { data, layout: {
        ...base, barmode: "stack",
        xaxis: { ...base.xaxis, tickmode: "array", tickvals: sampleNames, ticktext: sampleNames.map(xFmt) },
        yaxis: { ...base.yaxis, title: { text: is100 ? "Relative Abundance (%)" : "Abundance (%)" } },
      }};
    }

    // ALPHA DIVERSITY — bar per sample
    if (activeTab2.type === "alpha" && rows.length > 1) {
      const metric = activeTab2.metric || "Shannon";
      const mIdx   = rows[0].indexOf(metric);
      const sIdx   = rows[0].indexOf("Sample");
      if (mIdx < 0) return null;
      const samples = rows.slice(1).map(r => alias(sIdx >= 0 ? r[sIdx] : r[r.length-1]));
      const values  = rows.slice(1).map(r => num(r[mIdx]));
      const data = [{ x: samples, y: values, type: "bar",
        marker: { color: samples.map((s, i) => colors[s] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]) },
        text: values.map(v => v.toFixed(3)), textposition: "outside",
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, tickmode: "array", tickvals: samples, ticktext: samples.map(xFmt) },
        yaxis: { ...base.yaxis, title: { text: metric } },
      }};
    }

    // BAR (read counts / asv summary)
    if (activeTab2.type === "bar" && rows.length > 1) {
      const samples = rows.slice(1).map(r => alias(r[0]));
      const reads   = rows.slice(1).map(r => num(r[1]));
      const asvs    = rows.slice(1).map(r => num(r[2]));
      const data = [
        { name: "Reads", x: samples, y: reads, type: "bar",
          marker: { color: samples.map((s,i) => colors[s+"_reads"] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]) } },
        { name: "ASVs",  x: samples, y: asvs,  type: "bar", yaxis: "y2",
          marker: { color: samples.map((_,i) => colors["asvs_"+i] || DEFAULT_COLORS[(i+10) % DEFAULT_COLORS.length]) } },
      ];
      return { data, layout: { ...base, barmode: "group",
        xaxis:  { ...base.xaxis, tickmode: "array", tickvals: samples, ticktext: samples.map(xFmt) },
        yaxis:  { ...base.yaxis, title: { text: "Reads" } },
        yaxis2: { title: { text: "ASVs" }, overlaying: "y", side: "right",
                  tickfont: { size: font.axisSize } },
      }};
    }

    // LINE (rarefaction)
    if (activeTab2.type === "line" && rows.length > 1) {
      const depths  = rows.slice(1).map(r => num(r[0]));
      const samples = rows[0].slice(1);
      const data = samples.map((name, ci) => ({
        name: shortName(name),
        x: depths,
        y: rows.slice(1).map(r => num(r[ci + 1])),
        type: "scatter", mode: "lines",
        line: { color: colors[name] || DEFAULT_COLORS[ci % DEFAULT_COLORS.length], width: 2 },
      }));
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: "Sequencing Depth" } },
        yaxis: { ...base.yaxis, title: { text: "Observed OTUs" } },
      }};
    }

    // SCATTER — PCoA (pca_scores.csv) or NMDS (nmds_bray.csv / nmds_jaccard.csv)
    // Columns: Sample, PC1/NMDS1, PC2/NMDS2[, PC3][, Group, PC1_var, PC2_var, Stress]
    if (activeTab2.type === "scatter" && rows.length > 1) {
      const hdr = rows[0];
      const axis1Name = hdr[1] || "PC1";   // e.g. "PC1" or "NMDS1"
      const axis2Name = hdr[2] || "PC2";
      const varIdx1 = hdr.indexOf(axis1Name + "_var");
      const varIdx2 = hdr.indexOf(axis2Name + "_var");
      const grpIdx   = hdr.indexOf("Group");
      const stressIdx = hdr.indexOf("Stress");
      const var1 = varIdx1 >= 0 ? num(rows[1][varIdx1]) : 0;
      const var2 = varIdx2 >= 0 ? num(rows[1][varIdx2]) : 0;
      const stress = stressIdx >= 0 ? num(rows[1][stressIdx]) : 0;
      const samples = rows.slice(1).map(r => alias(r[0]));
      const xVals = rows.slice(1).map(r => num(r[1]));
      const yVals = rows.slice(1).map(r => num(r[2]));
      const groups = grpIdx >= 0 ? rows.slice(1).map(r => r[grpIdx] || "no group") : null;
      const data = groups
        ? Array.from(new Set(groups)).map(g => {
            const idx = groups.map((gg,i) => gg === g ? i : -1).filter(i => i >= 0);
            return {
              name: g, x: idx.map(i => xVals[i]), y: idx.map(i => yVals[i]),
              mode: "markers+text", type: "scatter",
              text: idx.map(i => samples[i]), textposition: "top center",
              marker: { size: 12, color: colors[g] || DEFAULT_COLORS[Array.from(new Set(groups)).indexOf(g) % DEFAULT_COLORS.length] },
            };
          })
        : [{ x: xVals, y: yVals, mode: "markers+text", type: "scatter",
             text: samples, textposition: "top center",
             marker: { size: 12, color: samples.map((s,i) => colors[s] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]) },
           }];
      const xLabel = var1 > 0 ? `${axis1Name} [${var1}%]` : axis1Name;
      const yLabel = var2 > 0 ? `${axis2Name} [${var2}%]` : axis2Name;
      const annotations = stress > 0
        ? [{ text: `stress = ${stress.toFixed(4)}`, showarrow: false,
             xref: "paper", yref: "paper", x: 0.02, y: 0.97,
             font: { size: 11, color: "#64748b" }, align: "left" as const }]
        : [];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: xLabel }, zeroline: true },
        yaxis: { ...base.yaxis, title: { text: yLabel }, zeroline: true },
        annotations,
      }};
    }

    // SCREE — pca_scree.csv: PC (integer), Variance (%)
    if (activeTab2.type === "scree" && rows.length > 1) {
      const pcs = rows.slice(1).map(r => `PC${r[0]}`);
      const variances = rows.slice(1).map(r => num(r[1]));
      const data = [{
        x: pcs, y: variances, type: "bar",
        marker: { color: "#3b82f6" },
        text: variances.map(v => v.toFixed(1) + "%"),
        textposition: "outside",
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: "Principal Component" } },
        yaxis: { ...base.yaxis, title: { text: "Variance Explained (%)" }, rangemode: "tozero" },
        margin: { l: 60, r: 20, t: 60, b: 80 },
      }};
    }

    // LONGLINE — long-format line chart (Sample, x_col, y_col)
    // Handles: shannon_rarefaction.csv (Depth, Shannon) and rank_abundance.csv (Rank, RelAbundance)
    if (activeTab2.type === "longline" && rows.length > 1) {
      const hdr = rows[0];
      const sampleIdx = hdr.indexOf("Sample");
      if (sampleIdx < 0) return null;
      const isShannon = activeTab2.id === "shan_rar";
      const isRank    = activeTab2.id === "rank_abund";
      const xColName  = isShannon ? "Depth" : "Rank";
      const yColName  = isShannon ? "Shannon" : "RelAbundance";
      const xIdx = hdr.indexOf(xColName);
      const yIdx = hdr.indexOf(yColName);
      if (xIdx < 0 || yIdx < 0) return null;
      // Group rows by sample
      const sampleMap: Record<string, {x: number[], y: number[]}> = {};
      for (const r of rows.slice(1)) {
        const rawName = r[sampleIdx];
        const sName = alias(rawName);
        if (!sampleMap[sName]) sampleMap[sName] = {x:[], y:[]};
        sampleMap[sName].x.push(num(r[xIdx]));
        sampleMap[sName].y.push(num(r[yIdx]));
      }
      const entries = Object.entries(sampleMap);
      const data = entries.map(([name, {x, y}], ci) => ({
        name, x, y,
        type: "scatter", mode: "lines",
        line: { color: colors[name] || DEFAULT_COLORS[ci % DEFAULT_COLORS.length], width: 2 },
      }));
      const xTitle = isShannon ? "Sequencing Depth" : "ASV Rank";
      const yTitle = isShannon ? "Shannon Index" : "Relative Abundance (%)";
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: xTitle } },
        yaxis: { ...base.yaxis, title: { text: yTitle }, ...(isRank ? {type: "log" as const} : {}) },
      }};
    }

    // SPECACCUM — species accumulation curve: Sites, Richness, SD
    if (activeTab2.type === "specaccum" && rows.length > 1) {
      const sites    = rows.slice(1).map(r => num(r[0]));
      const richness = rows.slice(1).map(r => num(r[1]));
      const sd       = rows.slice(1).map(r => num(r[2]));
      const upper    = richness.map((r,i) => r + sd[i]);
      const lower    = richness.map((r,i) => Math.max(0, r - sd[i]));
      const data = [
        {
          x: [...sites, ...sites.slice().reverse()],
          y: [...upper, ...lower.slice().reverse()],
          fill: "toself" as const, fillcolor: "rgba(59,130,246,0.15)",
          line: { color: "transparent" }, type: "scatter" as const, mode: "lines" as const,
          name: "± SD", showlegend: true,
        },
        {
          x: sites, y: richness, type: "scatter" as const, mode: "lines+markers" as const,
          line: { color: colors["Richness"] || "#3b82f6", width: 2 },
          marker: { size: 6, color: colors["Richness"] || "#3b82f6" },
          name: "Richness",
        },
      ];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: "Number of Samples" } },
        yaxis: { ...base.yaxis, title: { text: "Cumulative ASV Richness" } },
      }};
    }

    // MULTIALPHA — alpha_diversity.csv: all metrics in 2×N subplot grid, one scatter trace per sample
    // Each sample gets a colored dot in each subplot; legend shows sample names → click to change color
    if (activeTab2.type === "multialpha" && rows.length > 1) {
      const ALL_METRICS = ["Observed", "Chao1", "Shannon", "Simpson", "PD"];
      const presentMetrics = ALL_METRICS.filter(m => rows[0].indexOf(m) >= 0);
      const sIdx    = rows[0].indexOf("Sample");
      const rawNames = rows.slice(1).map(r => sIdx >= 0 ? r[sIdx] : r[0]);

      const cols  = Math.min(presentMetrics.length, 2);
      const rowsN = Math.ceil(presentMetrics.length / cols);
      const gap   = 0.08;
      const cellW = (1 - gap * (cols - 1)) / cols;
      const cellH = (1 - gap * (rowsN - 1) - 0.04) / rowsN;

      const data: any[] = [];
      const axisLayout: any = {};
      const annotations: any[] = [];

      presentMetrics.forEach((metric, mi) => {
        const mIdx   = rows[0].indexOf(metric);
        const col    = mi % cols;
        const row    = Math.floor(mi / cols);
        const xStart = col * (cellW + gap);
        const yStart = 1 - (row + 1) * cellH - row * gap;
        const axN    = mi === 0 ? "" : String(mi + 1);

        // One trace per sample for this subplot — enables per-sample legend + color picker
        rawNames.forEach((rawName, si) => {
          const sampleLabel = alias(rawName);
          const color = colors[rawName] || DEFAULT_COLORS[si % DEFAULT_COLORS.length];
          data.push({
            x: [sampleLabel],
            y: [num(rows[si + 1][mIdx])],
            type: "scatter", mode: "markers",
            name: rawName,              // raw name so legendclick can look up colors[rawName]
            legendgroup: rawName,       // same group across all subplots → appears once in legend
            showlegend: mi === 0,       // show legend entry only for first subplot
            marker: { size: 11, color, opacity: 0.9, line: { width: 1, color: "rgba(255,255,255,0.4)" } },
            xaxis: `x${axN}`, yaxis: `y${axN}`,
          });
        });

        axisLayout[`xaxis${axN}`] = {
          domain: [xStart, xStart + cellW], anchor: `y${axN}`,
          showticklabels: false, showgrid: false, zeroline: false,
        };
        axisLayout[`yaxis${axN}`] = {
          domain: [yStart, yStart + cellH], anchor: `x${axN}`,
          title: { text: metric, font: { size: font.axisSize } },
          showgrid: font.showGrid, gridcolor: "rgba(148,163,184,0.15)",
          zeroline: false, tickfont: { size: font.axisSize - 1 },
          rangemode: "tozero" as const,
        };
        annotations.push({
          text: `<b>${metric}</b>`,
          xref: "paper", yref: "paper",
          x: xStart + cellW / 2, y: yStart + cellH + 0.01,
          xanchor: "center", yanchor: "bottom", showarrow: false,
          font: { size: font.axisSize + 1, color: "#475569" },
        });
      });

      return { data, layout: {
        ...base, ...axisLayout,
        annotations, showlegend: true,
        legend: { font: { size: font.legendSize }, orientation: "v" as const,
                  x: 1.01, y: 1, xanchor: "left" },
        margin: { l: 60, r: 140, t: 50, b: 20 },
      }};
    }

    // HEATMAP (beta_heatmap.csv: similarity; bray_curtis_distance_matrix.csv: distance)
    if (activeTab2.type === "heatmap" && rows.length > 1) {
      const sampleNames = rows.slice(1).map(r => alias(r[0]));
      const zValues = rows.slice(1).map(r => r.slice(1).map(v => num(v)));
      // Auto-detect: diagonal ≈ 0 → distance matrix; diagonal ≈ 1 → similarity
      const diag0 = zValues[0] ? zValues[0][0] : 1;
      const isDist = diag0 < 0.1;
      // Use user-selected colorscale if set, otherwise auto-select by type
      const cs = font.colorscale && font.colorscale !== "Viridis"
        ? font.colorscale
        : (isDist ? "Reds" : "Blues");
      const data = [{ z: zValues, x: sampleNames, y: sampleNames,
        type: "heatmap",
        colorscale: cs,
        zmin: 0, zmax: 1,
        hoverongaps: false,
        colorbar: { title: isDist ? "Distance" : "Similarity" },
        reversescale: cs === "Blues",  // Blues: dark=high similarity
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, tickangle: -45, tickmode: "array", tickvals: sampleNames, ticktext: sampleNames.map(xFmt) },
        yaxis: { ...base.yaxis, tickmode: "array", tickvals: sampleNames, ticktext: sampleNames.map(xFmt), autorange: "reversed" },
        margin: { l: 100, r: 20, t: 60, b: 120 },
      }};
    }

    // TAXHEATMAP — taxonomy_phylum/family/genus.csv: rows=samples, cols=taxa
    // Renders as sample × taxon abundance heatmap (interactive, editable)
    if (activeTab2.type === "taxheatmap" && rows.length > 1) {
      const taxonNames  = rows[0].slice(1);
      const sampleNames = rows.slice(1).map(r => alias(r[0]));
      const zValues     = rows.slice(1).map(r => r.slice(1).map(v => num(v)));
      // Sort taxa by total abundance (descending) so most dominant are first
      const totals = taxonNames.map((_, ci) =>
        zValues.reduce((s, row) => s + (row[ci] || 0), 0));
      const order = totals.map((_, i) => i).sort((a, b) => totals[b] - totals[a]);
      const sortedTaxa   = order.map(i => taxonNames[i]);
      const sortedZ      = zValues.map(row => order.map(i => row[i]));
      const level = activeTab2.label.replace("Tax. Heatmap ", "");
      const data = [{
        z: sortedZ, x: sortedTaxa, y: sampleNames,
        type: "heatmap" as const,
        colorscale: font.colorscale || "Viridis",
        hoverongaps: false,
        colorbar: { title: "Abundance (%)", thickness: 16 },
        hovertemplate: `<b>%{y}</b><br>${level}: %{x}<br>Abundance: %{z:.2f}%<extra></extra>`,
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, tickangle: -50, title: { text: level } },
        yaxis: { ...base.yaxis, title: { text: "Sample" }, autorange: "reversed" as const },
        margin: { l: 110, r: 60, t: 60, b: 160 },
      }};
    }

    // READ TRACKING — read_tracking.csv: "","input","filtered","denoised[FR]","merged","nonchim"
    // Shows read retention through each pipeline step as lines per sample
    if (activeTab2.type === "readtrack" && rows.length > 1) {
      const steps   = rows[0].slice(1);     // pipeline step names (variable per pipeline)
      const data = rows.slice(1).map((row, ri) => ({
        name: alias(row[0]),
        x: steps,
        y: steps.map((_, si) => num(row[si + 1])),
        type: "scatter" as const,
        mode: "lines+markers" as const,
        line: { color: colors[row[0]] || DEFAULT_COLORS[ri % DEFAULT_COLORS.length], width: 2 },
        marker: { size: 7, color: colors[row[0]] || DEFAULT_COLORS[ri % DEFAULT_COLORS.length] },
      }));
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: "Pipeline Step" } },
        yaxis: { ...base.yaxis, title: { text: "Read Count" }, rangemode: "tozero" as const },
      }};
    }

    // BOX (OTU distribution)
    if (activeTab2.type === "box" && rows.length > 1) {
      // otu_distribution.csv: Sample, mean, median, ...
      const samples = rows.slice(1).map(r => alias(r[0]));
      const means   = rows.slice(1).map(r => num(r[1]));
      const data = [{ x: samples, y: means, type: "bar",
        marker: { color: samples.map((s,i) => colors[s] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]) },
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, tickmode: "array", tickvals: samples, ticktext: samples.map(xFmt) },
        yaxis: { ...base.yaxis, title: { text: "Mean reads/OTU" } },
      }};
    }

    return null;
  }

  useEffect(() => { plotted.current = false; }, [activeTab, selJob]);

  // ── Export PNG (download + save to edit_charts/ on server) ──
  const exportPNG = async () => {
    if (!chartRef.current || !activeTab) return;
    // Download to browser
    Plotly.downloadImage(chartRef.current, {
      format: "png", filename: `${selJob}_${activeTab.id}`,
      width: 1400, height: 900, scale: 2,
    });
    // Also upload to server → edit_charts/{tab_id}.png (included in ZIP)
    try {
      const imgUrl: string = await Plotly.toImage(chartRef.current,
        { format: "png", width: 1400, height: 900, scale: 2 });
      const b64 = imgUrl.split(",")[1];
      await fetch(`${API}/results/${selJob}/edit_charts/png`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ tab_id: activeTab.id, b64 }),
      });
    } catch { /* non-fatal — local download still works */ }
  };

  // ── Save settings → edit_charts/settings.json on server ──
  const saveSettings = async () => {
    if (!selJob) return;
    const payload = { font, colors, sampleAliases };
    localStorage.setItem(`prev_settings_${selJob}`, JSON.stringify(payload));
    // Persist to server inside edit_charts/ folder (included in ZIP export)
    try {
      await fetch(`${API}/results/${selJob}/edit_charts/settings`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
    } catch { /* non-fatal — localStorage backup still works */ }
    setSaved(true);
    setTimeout(() => setSaved(false), 2000);
  };

  // ── Re-run R visualization for old/missing charts ─────────
  const rerunViz = async () => {
    if (!selJob || rerunning) return;
    setRerunning(true);
    try {
      const r = await fetch(`${API}/results/${selJob}/rerun-viz`, { method: "POST" });
      if (r.ok) {
        // Reload table info after re-run
        const info: TableInfo = await fetch(`${API}/results/${selJob}/tables`).then(x => x.json());
        setTableInfo(info);
        const avail = ALL_TABS.filter(t => info.tables.includes(t.file));
        setAvailTabs(avail);
        if (avail.length > 0) setActiveTab(avail[0]);
      }
    } catch { /* ignore */ }
    setRerunning(false);
  };

  // ── Export full results as ZIP (with custom-named PNG charts) ──
  // Flow: render PNGs → upload to server (saved in preview_charts/) → server ZIPs everything
  const exportResults = async () => {
    if (!selJob || exporting) return;
    setExporting(true);

    // Hidden off-screen div for rendering Plotly charts
    const tmpDiv = document.createElement("div");
    tmpDiv.style.cssText = "position:fixed;top:-9999px;left:-9999px;width:1200px;height:700px;background:#fff;";
    document.body.appendChild(tmpDiv);

    try {
      // ── 1. Render each available tab as PNG ──────────────────
      const charts: { name: string; b64: string }[] = [];
      for (let i = 0; i < availTabs.length; i++) {
        const tab = availTabs[i];
        if (tab.type === "pdf") continue; // PDF tabs are already in r_plots/ — no Plotly render needed
        setExportStep(`Rendering chart ${i + 1}/${availTabs.length}: ${tab.label}…`);
        try {
          const text = await fetch(`${API}/results/${selJob}/table/${tab.file}`).then(r => r.text());
          const tabData = parseCSV(text);
          const plotData = buildPlotly(tabData, tab);
          if (!plotData) continue;
          await Plotly.newPlot(tmpDiv, plotData.data, plotData.layout,
            { responsive: false, displayModeBar: false });
          const imgUrl: string = await Plotly.toImage(tmpDiv,
            { format: "png", width: 1400, height: 900, scale: 2 });
          charts.push({ name: `${i + 1}_${tab.id}.png`, b64: imgUrl.split(",")[1] });
        } catch { /* skip this tab */ }
      }
      Plotly.purge(tmpDiv);

      // ── 2. Upload PNGs to server (saved into preview_charts/) ─
      setExportStep("Uploading charts to server…");
      await fetch(`${API}/results/${selJob}/preview-charts`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ charts }),
      });

      // ── 3. Trigger server-side ZIP download (all files) ───────
      setExportStep("Building ZIP…");
      const resp = await fetch(`${API}/results/${selJob}/download`);
      if (!resp.ok) throw new Error("Download failed");
      const blob = await resp.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${job?.job_name || selJob}_results.zip`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(url), 30000);

    } catch (e) {
      alert("Export failed: " + String(e));
    } finally {
      document.body.removeChild(tmpDiv);
      setExporting(false);
      setExportStep("");
    }
  };

  // ── Helpers ───────────────────────────────────────────────
  const job      = jobs.find(j => j.job_id === selJob);
  const filtered = jobs.filter(j =>
    j.job_name.toLowerCase().includes(search.toLowerCase()) ||
    j.job_id.toLowerCase().includes(search.toLowerCase()) ||
    j.marker.toLowerCase().includes(search.toLowerCase())
  );

  // group tabs by "group"
  const tabGroups: Record<string, TabDef[]> = {};
  availTabs.forEach(t => {
    const g = (t as any).group || "Other";
    if (!tabGroups[g]) tabGroups[g] = [];
    tabGroups[g].push(t);
  });
  const GROUP_ORDER = ["Taxonomy","Alpha","Beta","Other"];
  // Beta only shows when charts exist (needs ≥2 samples)

  // Series for customization panel
  // For PDF tabs csvData is null — fall back to knownSamples for color assignment
  const series = csvData && activeTab
    ? getSeries(activeTab, csvData)
    : (activeTab?.type === "pdf" ? knownSamples : []);
  // Unique sample names (raw) for alias editing panel
  const uniqueSamples = csvData && activeTab
    ? getUniqueSamples(activeTab, csvData)
    : (activeTab?.type === "pdf" ? knownSamples : []);

  return (
    <div className="prev-root">

      {/* ── Left sidebar ──────────────────────────────────── */}
      <aside className="prev-sidebar">
        <div className="prev-sidebar-header">
          <button className="prev-back" onClick={onClose}>← Back</button>
          <span className="prev-sidebar-title">Edit Charts</span>
        </div>

        <div className="prev-search-wrap">
          <input className="prev-search" placeholder="🔍 Search jobs…"
            value={search} onChange={e => setSearch(e.target.value)} />
        </div>

        <div className="prev-job-list">
          {filtered.length === 0 && (
            <div className="prev-no-jobs">No completed jobs</div>
          )}
          {filtered.map(j => (
            <div key={j.job_id}
              className={`prev-job-item ${j.job_id === selJob ? "active" : ""}`}
              onClick={() => setSelJob(j.job_id)}>
              <div className="prev-job-name">{j.job_name || j.job_id.slice(0,8)}</div>
              <div className="prev-job-meta">
                <span className="prev-marker-badge">{j.marker}</span>
                <span className="prev-job-id">{j.job_id.slice(0,8)}</span>
              </div>
            </div>
          ))}
        </div>
      </aside>

      {/* ── Right panel ───────────────────────────────────── */}
      <main className="prev-main">
        {!selJob ? (
          <div className="prev-empty">
            <div className="prev-empty-icon">📊</div>
            <div>Select a completed job from the left to preview results</div>
          </div>
        ) : (
          <>
            {/* Top bar */}
            <div className="prev-topbar">
              <div className="prev-topbar-left">
                <span className="prev-job-title">{job?.job_name || selJob.slice(0,8)}</span>
                {tableInfo?.summary && (
                  <span className="prev-summary-chips">
                    <span className="prev-chip">{tableInfo.summary.marker}</span>
                    <span className="prev-chip">{tableInfo.summary.n_samples} samples</span>
                    <span className="prev-chip">{tableInfo.summary.n_taxa?.toLocaleString()} taxa</span>
                  </span>
                )}
              </div>
              <div className="prev-topbar-right">
                <button className="prev-btn prev-btn-pdf"
                  onClick={() => setPdfPanel(p => !p)}>
                  📄 Original PDFs
                </button>
                <button className={`prev-btn prev-btn-rerun ${rerunning ? "running" : ""}`}
                  onClick={rerunViz} disabled={rerunning} title="Re-run R visualization to regenerate charts">
                  {rerunning ? "⏳ Re-running…" : "🔄 Re-run Viz"}
                </button>
                <button className="prev-btn prev-btn-png" onClick={exportPNG}>
                  🖼 Save PNG
                </button>
                <button
                  className={`prev-btn prev-btn-save ${saved ? "saved" : ""}`}
                  onClick={saveSettings}>
                  {saved ? "✓ Saved!" : "💾 Save"}
                </button>
                <button
                  className={`prev-btn prev-btn-export ${exporting ? "exporting" : ""}`}
                  onClick={exportResults} disabled={exporting}
                  title={exporting ? exportStep : "Export results + custom-named PNG charts as ZIP"}>
                  {exporting
                    ? (exportStep.startsWith("Rendering")
                        ? `⏳ Chart ${exportStep.match(/(\d+\/\d+)/)?.[1] ?? "…"}`
                        : "⏳ Packing…")
                    : "📦 Export Result"}
                </button>
              </div>
            </div>

            {/* Original PDFs panel */}
            {pdfPanel && tableInfo && tableInfo.plots.length > 0 && (
              <div className="prev-pdf-panel">
                <div className="prev-pdf-panel-title">📄 Pre-generated PDF files (R)</div>
                <div className="prev-pdf-list">
                  {tableInfo.plots.map(f => (
                    <a key={f} className="prev-pdf-link"
                      href={`${API}/results/${selJob}/plot/${f}`}
                      target="_blank" rel="noopener noreferrer">
                      {f.replace(/^\d+_/, "").replace(".pdf","").replace(/_/g," ")}
                    </a>
                  ))}
                </div>
              </div>
            )}

            {/* Tab bar */}
            {/* Banner when no charts at all */}
            {tableInfo && availTabs.length === 0 && (
              <div className="prev-rerun-banner">
                <span>⚠️ No chart tables found — click <b>🔄 Re-run Viz</b> (top-right) to regenerate from existing results.</span>
              </div>
            )}

            {availTabs.length > 0 ? (
              <div className="prev-tabbar">
                {GROUP_ORDER.filter(g => tabGroups[g]).map(g => (
                  <div key={g} className="prev-tab-group">
                    <span className="prev-tab-group-label">{g}</span>
                    {tabGroups[g].map(t => (
                      <button key={t.id}
                        className={`prev-tab ${activeTab?.id === t.id ? "active" : ""}`}
                        onClick={() => { setActiveTab(t); setCsvData(null); }}>
                        {t.label}
                      </button>
                    ))}
                  </div>
                ))}
              </div>
            ) : (
              <div className="prev-no-tables">
                No chart tables found for this job yet — run the pipeline to generate results.
              </div>
            )}

            {/* Chart area + Customize */}
            {activeTab && (
              <div className="prev-content">
                {/* Chart — Plotly for interactive types, iframe for PDF type */}
                <div className="prev-chart-wrap">
                  {loading && <div className="prev-loading">⏳ Loading data…</div>}
                  {activeTab.type === "pdf" ? (
                    <iframe
                      key={activeTab.pdfFile}
                      src={`${API}/results/${selJob}/plot/${activeTab.pdfFile}`}
                      className="prev-pdf-iframe"
                      title={activeTab.label}
                    />
                  ) : (
                    <div ref={chartRef} className="prev-chart" />
                  )}
                </div>

                {/* Customize panel — available for all tab types */}
                <div className="prev-customize">
                  <div className="prev-cust-header"
                    onClick={() => setShowCustomize(s => !s)}>
                    🎨 Customize
                    <span className="prev-cust-toggle">{showCustomize ? "▲" : "▼"}</span>
                  </div>

                  {showCustomize && (
                    <div className="prev-cust-body">

                      {/* Chart style */}
                      <div className="prev-cust-section">
                        <div className="prev-cust-section-title">Chart Style</div>
                        {activeTab.type === "pdf" && (
                          <div className="prev-cust-pdf-note">
                            ℹ️ Style settings apply to interactive chart tabs (not PDF). Sample name aliases and colors are saved globally for this job.
                          </div>
                        )}
                        <div className="prev-cust-row">
                          <label>Title</label>
                          <input className="prev-cust-input"
                            value={font.titleText}
                            onChange={e => setFont(f => ({ ...f, titleText: e.target.value }))}
                            placeholder="Chart title…" />
                        </div>
                        <div className="prev-cust-row">
                          <label>Title size</label>
                          <input type="range" min={10} max={28} value={font.titleSize}
                            onChange={e => setFont(f => ({ ...f, titleSize: +e.target.value }))} />
                          <span className="prev-cust-val">{font.titleSize}px</span>
                        </div>
                        <div className="prev-cust-row">
                          <label>Axis size</label>
                          <input type="range" min={8} max={20} value={font.axisSize}
                            onChange={e => setFont(f => ({ ...f, axisSize: +e.target.value }))} />
                          <span className="prev-cust-val">{font.axisSize}px</span>
                        </div>
                        <div className="prev-cust-row">
                          <label>Legend size</label>
                          <input type="range" min={8} max={18} value={font.legendSize}
                            onChange={e => setFont(f => ({ ...f, legendSize: +e.target.value }))} />
                          <span className="prev-cust-val">{font.legendSize}px</span>
                        </div>
                        <div className="prev-cust-row">
                          <label>Style</label>
                          <div className="prev-toggle-group">
                            <button className={`prev-toggle ${font.titleBold ? "on":""}`}
                              onClick={() => setFont(f=>({...f, titleBold: !f.titleBold}))}>
                              <b>B</b>
                            </button>
                            <button className={`prev-toggle ${font.titleItalic ? "on":""}`}
                              onClick={() => setFont(f=>({...f, titleItalic: !f.titleItalic}))}>
                              <i>I</i>
                            </button>
                            <button className={`prev-toggle ${font.showGrid ? "on":""}`}
                              onClick={() => setFont(f=>({...f, showGrid: !f.showGrid}))}>
                              Grid
                            </button>
                          </div>
                        </div>
                        {(activeTab.type === "taxonomy") && (
                          <div className="prev-cust-row">
                            <label>Bar type</label>
                            <div className="prev-toggle-group">
                              {(["bar","bar100"] as const).map(ct => (
                                <button key={ct}
                                  className={`prev-toggle ${font.chartType === ct ? "on":""}`}
                                  onClick={() => setFont(f=>({...f, chartType: ct}))}>
                                  {ct === "bar" ? "Stacked" : "100%"}
                                </button>
                              ))}
                            </div>
                          </div>
                        )}
                        {/* X-axis label controls */}
                        <div className="prev-cust-row">
                          <label>X label size</label>
                          <input type="range" min={7} max={20} value={font.xSize}
                            onChange={e => setFont(f => ({ ...f, xSize: +e.target.value }))} />
                          <span className="prev-cust-val">{font.xSize}px</span>
                        </div>
                        <div className="prev-cust-row">
                          <label>X label style</label>
                          <div className="prev-toggle-group">
                            <button className={`prev-toggle ${font.xBold ? "on":""}`}
                              onClick={() => setFont(f => ({ ...f, xBold: !f.xBold }))}>
                              <b>B</b>
                            </button>
                            <button className={`prev-toggle ${font.xItalic ? "on":""}`}
                              onClick={() => setFont(f => ({ ...f, xItalic: !f.xItalic }))}>
                              <i>I</i>
                            </button>
                          </div>
                        </div>
                      </div>

                      {/* Colorscale picker — for heatmap and taxheatmap tabs */}
                      {(activeTab.type === "taxheatmap" || activeTab.type === "heatmap") && (
                        <div className="prev-cust-section">
                          <div className="prev-cust-section-title">Color Scale</div>
                          <div className="prev-cs-grid">
                            {COLORSCALES.map(cs => (
                              <button key={cs.name}
                                className={`prev-cs-btn${(font.colorscale || "Viridis") === cs.name ? " active" : ""}`}
                                onClick={() => setFont(f => ({ ...f, colorscale: cs.name }))}>
                                <div className="prev-cs-swatch" style={{ background: cs.gradient }} />
                                <span>{cs.name}</span>
                              </button>
                            ))}
                          </div>
                        </div>
                      )}

                      {/* Sample name aliases */}
                      {uniqueSamples.length > 0 && (
                        <div className="prev-cust-section">
                          <div className="prev-cust-section-title">Sample Names (X-axis)</div>
                          <div className="prev-alias-list">
                            {uniqueSamples.map(raw => (
                              <div key={raw} className="prev-alias-row">
                                <span className="prev-alias-orig" title={raw}>
                                  {shortName(raw)}
                                </span>
                                <span className="prev-alias-arrow">→</span>
                                <input
                                  className="prev-alias-input"
                                  value={sampleAliases[raw] ?? ""}
                                  placeholder={shortName(raw)}
                                  onChange={e => setSampleAliases(a => ({
                                    ...a,
                                    [raw]: e.target.value,
                                  }))}
                                />
                              </div>
                            ))}
                          </div>
                          <button className="prev-btn-reset"
                            onClick={() => setSampleAliases({})}>
                            ↺ Reset names
                          </button>
                        </div>
                      )}

                      {/* Colors — click legend items on the chart to change colors */}
                      {activeTab.type !== "pdf" && (
                        <div className="prev-cust-section prev-cust-color-hint">
                          <div className="prev-cust-section-title">Colors</div>
                          <div className="prev-color-hint-msg">
                            🎨 Click any item in the chart legend to change its color
                          </div>
                          <button className="prev-btn-reset" onClick={() => setColors({})}>
                            ↺ Reset all colors to defaults
                          </button>
                        </div>
                      )}
                    </div>
                  )}
                </div>
              </div>
            )}
          </>
        )}
      </main>

      {/* ── Floating legend color picker ────────────────────── */}
      {legendPicker && (
        <div className="prev-legend-picker"
          style={{ left: legendPicker.x + 12, top: legendPicker.y - 8 }}
          onMouseDown={e => e.stopPropagation()}>
          <input
            type="color"
            className="prev-lp-swatch"
            autoFocus
            value={legendPicker.color}
            onChange={e => {
              const c = e.target.value;
              setLegendPicker(s => s ? { ...s, color: c } : null);
              setColors(prev => ({ ...prev, [legendPicker.name]: c }));
            }}
          />
          <span className="prev-lp-name" title={legendPicker.name}>
            {legendPicker.name.length > 26 ? legendPicker.name.slice(0, 26) + "…" : legendPicker.name}
          </span>
          <button className="prev-lp-close" onClick={() => setLegendPicker(null)}>✕</button>
        </div>
      )}
    </div>
  );
}

// ── Tiny fetch helper (avoid adding axios dep just for one call) ──
async function axios_get<T>(url: string): Promise<T> {
  const r = await fetch(url);
  return r.json();
}
// Make window.Plotly available
declare global { interface Window { Plotly: any; } }
