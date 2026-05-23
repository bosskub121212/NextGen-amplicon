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
  type: "taxonomy"|"alpha"|"bar"|"line"|"scatter"|"box"|"heatmap";
  metric?: string; group?: string;
}
const ALL_TABS: TabDef[] = [
  { id:"phylum",      label:"Phylum",       file:"taxonomy_phylum.csv",   type:"taxonomy", group:"Taxonomy" },
  { id:"class",       label:"Class",        file:"taxonomy_class.csv",    type:"taxonomy", group:"Taxonomy" },
  { id:"order",       label:"Order",        file:"taxonomy_order.csv",    type:"taxonomy", group:"Taxonomy" },
  { id:"family",      label:"Family",       file:"taxonomy_family.csv",   type:"taxonomy", group:"Taxonomy" },
  { id:"genus",       label:"Genus",        file:"taxonomy_genus.csv",    type:"taxonomy", group:"Taxonomy" },
  { id:"shannon",     label:"Shannon",      file:"alpha_diversity.csv",   type:"alpha", metric:"Shannon",  group:"Alpha" },
  { id:"observed",    label:"Observed",     file:"alpha_diversity.csv",   type:"alpha", metric:"Observed", group:"Alpha" },
  { id:"chao1",       label:"Chao1",        file:"alpha_diversity.csv",   type:"alpha", metric:"Chao1",    group:"Alpha" },
  { id:"simpson",     label:"Simpson",      file:"alpha_diversity.csv",   type:"alpha", metric:"Simpson",  group:"Alpha" },
  { id:"reads",       label:"Read Counts",  file:"asv_summary.csv",       type:"bar",                      group:"Other" },
  { id:"rarefaction", label:"Rarefaction",  file:"rarefaction.csv",       type:"line",                     group:"Other" },
  { id:"pca",         label:"PCoA",         file:"pca_scores.csv",        type:"scatter",                  group:"Beta" },
  { id:"heatmap",     label:"Heatmap",      file:"beta_heatmap.csv",      type:"heatmap",                  group:"Beta" },
  { id:"otu",         label:"OTU Dist",     file:"otu_distribution.csv",  type:"box",                      group:"Other" },
];

// ── Font config ───────────────────────────────────────────────
interface FontConfig {
  titleSize: number; axisSize: number; legendSize: number;
  titleText: string; titleBold: boolean; titleItalic: boolean;
  chartType: "bar"|"bar100"|"line"; showGrid: boolean;
  xItalic: boolean; xBold: boolean; xSize: number;
}
const defaultFont = (): FontConfig => ({
  titleSize: 16, axisSize: 12, legendSize: 11,
  titleText: "", titleBold: false, titleItalic: false,
  chartType: "bar", showGrid: true,
  xItalic: false, xBold: false, xSize: 12,
});

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
  if (tab.type === "taxonomy") return rows[0].slice(1);
  if (tab.type === "alpha") {
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
  const [showCustomize, setShowCustomize] = useState(true);
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
      // Try server first
      try {
        const s = await fetch(`${API}/results/${selJob}/settings`).then(r => r.json());
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
      const avail = ALL_TABS.filter(t => info.tables.includes(t.file));
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
          const avail = ALL_TABS.filter(t => info.tables.includes(t.file));
          setAvailTabs(avail);
          if (avail.length > 0) setActiveTab(avail[0]);
        }).catch(() => {});
      });
  }, [selJob]);

  // ── Load CSV when tab changes ─────────────────────────────
  useEffect(() => {
    if (!activeTab || !selJob) return;
    setLoading(true);
    fetch(`${API}/results/${selJob}/table/${activeTab.file}`)
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
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, [activeTab, selJob]);

  // ── Re-render chart when data / colors / font / aliases changes ────
  // (renderChart is declared below but useCallback deps already cover it)
  const triggerRender = useCallback(() => {
    if (!csvData || !chartRef.current || !activeTab) return;
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

  // ── Get series names for a tab ────────────────────────────
  function getSeries(tab: TabDef, rows: string[][]): string[] {
    if (!rows.length) return [];
    if (tab.type === "taxonomy")  return rows.slice(1).map(r => r[0] || "Unknown").slice(0, 30);
    if (tab.type === "alpha")     return rows.length > 1 ? rows.slice(1).map(r => {
      const idx = rows[0].indexOf("Sample"); return idx >= 0 ? r[idx] : r[r.length-1];
    }) : [];
    if (tab.type === "line")      return rows[0]?.slice(1) || [];
    if (tab.type === "scatter")   return rows.slice(1).map(r => r[0]) || [];
    if (tab.type === "bar")       return rows.slice(1).map(r => shortName(r[0]));
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

    // TAXONOMY — stacked bar (absolute %) or 100% normalized
    if (activeTab2.type === "taxonomy" && rows.length > 1) {
      const headers = rows[0];
      const samples = headers.slice(1).map(alias);
      const taxaRows = rows.slice(1).slice(0, 30); // top 30
      const is100 = font.chartType === "bar100";
      // compute column sums for normalization
      const colSums = samples.map((_, ci) =>
        taxaRows.reduce((s, r) => s + num(r[ci + 1]), 0)
      );
      const data = taxaRows.map(row => ({
        name: row[0] || "Unknown",
        x: samples,
        y: samples.map((_, ci) => {
          const v = num(row[ci + 1]);
          return is100 && colSums[ci] > 0 ? (v / colSums[ci]) * 100 : v;
        }),
        type: "bar",
        marker: { color: colors[row[0]] || DEFAULT_COLORS[0] },
      }));
      return { data, layout: {
        ...base, barmode: "stack",
        xaxis: { ...base.xaxis, tickmode: "array", tickvals: samples, ticktext: samples.map(xFmt) },
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

    // SCATTER (PCoA / PCA) — pca_scores.csv: Sample, PC1, PC2[, Group, PC1_var, PC2_var]
    if (activeTab2.type === "scatter" && rows.length > 1) {
      const hdr = rows[0];
      const pc1varIdx = hdr.indexOf("PC1_var");
      const pc2varIdx = hdr.indexOf("PC2_var");
      const grpIdx    = hdr.indexOf("Group");
      const pc1Var = pc1varIdx >= 0 ? num(rows[1][pc1varIdx]) : 0;
      const pc2Var = pc2varIdx >= 0 ? num(rows[1][pc2varIdx]) : 0;
      const samples = rows.slice(1).map(r => alias(r[0]));
      const pc1 = rows.slice(1).map(r => num(r[1]));
      const pc2 = rows.slice(1).map(r => num(r[2]));
      const groups = grpIdx >= 0 ? rows.slice(1).map(r => r[grpIdx] || "no group") : null;
      const data = groups
        ? Array.from(new Set(groups)).map(g => {
            const idx = groups.map((gg,i) => gg === g ? i : -1).filter(i => i >= 0);
            return {
              name: g, x: idx.map(i => pc1[i]), y: idx.map(i => pc2[i]),
              mode: "markers+text", type: "scatter",
              text: idx.map(i => samples[i]), textposition: "top center",
              marker: { size: 12, color: colors[g] || DEFAULT_COLORS[Array.from(new Set(groups)).indexOf(g) % DEFAULT_COLORS.length] },
            };
          })
        : [{ x: pc1, y: pc2, mode: "markers+text", type: "scatter",
             text: samples, textposition: "top center",
             marker: { size: 12, color: samples.map((s,i) => colors[s] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]) },
           }];
      const pc1Label = pc1Var > 0 ? `PC1 [${pc1Var}%]` : "PC1";
      const pc2Label = pc2Var > 0 ? `PC2 [${pc2Var}%]` : "PC2";
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, title: { text: pc1Label }, zeroline: true },
        yaxis: { ...base.yaxis, title: { text: pc2Label }, zeroline: true },
      }};
    }

    // HEATMAP (beta_heatmap.csv: Sample, col1, col2, ...)
    if (activeTab2.type === "heatmap" && rows.length > 1) {
      const sampleNames = rows.slice(1).map(r => alias(r[0]));
      const zValues = rows.slice(1).map(r => r.slice(1).map(v => num(v)));
      const data = [{ z: zValues, x: sampleNames, y: sampleNames,
        type: "heatmap",
        colorscale: "Blues",
        zmin: 0, zmax: 1,
        hoverongaps: false,
        colorbar: { title: "Similarity" },
      }];
      return { data, layout: { ...base,
        xaxis: { ...base.xaxis, tickangle: -45, tickmode: "array", tickvals: sampleNames, ticktext: sampleNames.map(xFmt) },
        yaxis: { ...base.yaxis, tickmode: "array", tickvals: sampleNames, ticktext: sampleNames.map(xFmt), autorange: "reversed" },
        margin: { l: 100, r: 20, t: 60, b: 120 },
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

  // ── Export PNG ────────────────────────────────────────────
  const exportPNG = () => {
    if (!chartRef.current) return;
    Plotly.downloadImage(chartRef.current, {
      format: "png", filename: `${selJob}_${activeTab?.id || "chart"}`,
      width: 1400, height: 900, scale: 2,
    });
  };

  // ── Save settings to localStorage + server ───────────────
  const saveSettings = async () => {
    if (!selJob) return;
    const payload = { font, colors, sampleAliases };
    localStorage.setItem(`prev_settings_${selJob}`, JSON.stringify(payload));
    // Also persist to server so the ZIP export includes the settings
    try {
      await fetch(`${API}/results/${selJob}/settings`, {
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
  const series = csvData && activeTab ? getSeries(activeTab, csvData) : [];
  // Unique sample names (raw) for alias editing panel
  const uniqueSamples = csvData && activeTab ? getUniqueSamples(activeTab, csvData) : [];

  return (
    <div className="prev-root">

      {/* ── Left sidebar ──────────────────────────────────── */}
      <aside className="prev-sidebar">
        <div className="prev-sidebar-header">
          <button className="prev-back" onClick={onClose}>← Back</button>
          <span className="prev-sidebar-title">Results</span>
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
                {/* Chart */}
                <div className="prev-chart-wrap">
                  {loading && <div className="prev-loading">⏳ Loading data…</div>}
                  <div ref={chartRef} className="prev-chart" />
                </div>

                {/* Customize panel */}
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

                      {/* Color pickers */}
                      {series.length > 0 && (
                        <div className="prev-cust-section">
                          <div className="prev-cust-section-title">Colors</div>
                          <div className="prev-color-grid">
                            {series.slice(0,30).map((name, i) => (
                              <div key={name} className="prev-color-row">
                                <input type="color"
                                  value={colors[name] || DEFAULT_COLORS[i % DEFAULT_COLORS.length]}
                                  onChange={e => setColors(c => ({ ...c, [name]: e.target.value }))} />
                                <span className="prev-color-label" title={name}>
                                  {name.length > 22 ? name.slice(0,22)+"…" : name}
                                </span>
                              </div>
                            ))}
                          </div>
                          <button className="prev-btn-reset"
                            onClick={() => setColors({})}>
                            ↺ Reset to defaults
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
