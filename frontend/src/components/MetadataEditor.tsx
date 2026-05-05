import { useState, useRef, useCallback } from "react";

export interface MetaRow {
  sampleId: string;
  group: string;
  description: string;
  [key: string]: string;   // custom columns
}

interface Props {
  sampleNames: string[];
  onChange: (rows: MetaRow[]) => void;
}

export default function MetadataEditor({ sampleNames, onChange }: Props) {
  const [rows, setRows]             = useState<MetaRow[]>(() =>
    sampleNames.map((s) => ({ sampleId: s, group: "", description: "" }))
  );
  const [customCols, setCustomCols] = useState<string[]>([]);
  // Display labels for editable headers (key → display label)
  const [colLabels, setColLabels]   = useState<Record<string, string>>({
    group: "Group", description: "Description",
  });
  const [editingHeader, setEditingHeader] = useState<string | null>(null);
  const [editingHeaderVal, setEditingHeaderVal] = useState("");
  const [newColName, setNewColName] = useState("");
  const [showAddCol, setShowAddCol] = useState(false);
  const [csvStatus, setCsvStatus]   = useState<string | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);

  // ── Update a single cell ────────────────────────────────────
  const updateCell = useCallback((idx: number, key: string, val: string) => {
    setRows((prev) => {
      const next = prev.map((r, i) => i === idx ? { ...r, [key]: val } : r);
      onChange(next);
      return next;
    });
  }, [onChange]);

  // ── Add custom column ────────────────────────────────────────
  const addColumn = () => {
    const col = newColName.trim();
    if (!col || customCols.includes(col)) return;
    const colKey = col.toLowerCase().replace(/\s+/g, "_");
    setCustomCols((c) => [...c, colKey]);
    setColLabels((prev) => ({ ...prev, [colKey]: col }));
    setRows((prev) => {
      const next = prev.map((r) => ({ ...r, [colKey]: "" }));
      onChange(next);
      return next;
    });
    setNewColName("");
    setShowAddCol(false);
  };

  const removeColumn = (colKey: string) => {
    setCustomCols((c) => c.filter((x) => x !== colKey));
    setColLabels((prev) => { const n = { ...prev }; delete n[colKey]; return n; });
    setRows((prev) => {
      const next = prev.map((r) => {
        const { [colKey]: _, ...rest } = r;
        return rest as MetaRow;
      });
      onChange(next);
      return next;
    });
  };

  // ── Rename column header ──────────────────────────────────────
  const commitHeaderRename = (colKey: string) => {
    const newLabel = editingHeaderVal.trim();
    if (newLabel && newLabel !== colLabels[colKey]) {
      setColLabels((prev) => ({ ...prev, [colKey]: newLabel }));
    }
    setEditingHeader(null);
    setEditingHeaderVal("");
  };

  // ── Download CSV template ────────────────────────────────────
  const downloadTemplate = () => {
    const allCols = ["SampleID",
      colLabels["group"] || "Group",
      colLabels["description"] || "Description",
      ...customCols.map((c) => colLabels[c] || c)];
    const header  = allCols.join(",");
    const body    = rows.map((r) =>
      [r.sampleId, r.group, r.description, ...customCols.map((c) => r[c] ?? "")]
        .map((v) => `"${v}"`)
        .join(",")
    ).join("\n");
    const blob = new Blob([header + "\n" + body], { type: "text/csv" });
    const a    = document.createElement("a");
    a.href     = URL.createObjectURL(blob);
    a.download = "metadata_template.csv";
    a.click();
  };

  // ── Upload CSV and parse ─────────────────────────────────────
  const uploadCSV = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      const text  = ev.target?.result as string;
      const lines = text.trim().split(/\r?\n/);
      if (lines.length < 1) return;

      // Parse header
      const parseRow = (line: string) =>
        line.split(",").map((c) => c.trim().replace(/^"|"$/g, ""));

      const rawHeaders = parseRow(lines[0]);
      const headers    = rawHeaders.map((h) => h.toLowerCase());
      const sidIdx  = headers.findIndex((h) => h === "sampleid" || h === "sample_id" || h === "sample");
      const grpIdx  = headers.findIndex((h) => h === "group" || h === "treatment");
      const dscIdx  = headers.findIndex((h) => h === "description" || h === "desc");

      // Extra columns (beyond sampleId, group, description)
      const extraIdxs = rawHeaders.map((_, i) => i).filter(
        (i) => i !== sidIdx && i !== grpIdx && i !== dscIdx
      );
      const extraKeys   = extraIdxs.map((i) =>
        headers[i].replace(/\s+/g, "_") || `col_${i}`
      );
      const extraLabels = extraIdxs.map((i) => rawHeaders[i]);

      // Detect new custom columns from CSV
      const newCustom = extraKeys.filter((k) => !customCols.includes(k));
      if (newCustom.length > 0) {
        setCustomCols((c) => [...c, ...newCustom]);
        const newLabels: Record<string, string> = {};
        newCustom.forEach((k, ni) => {
          newLabels[k] = extraLabels[extraKeys.indexOf(k)] || k;
        });
        setColLabels((prev) => ({ ...prev, ...newLabels }));
      }

      // Update column label from CSV header names
      if (grpIdx >= 0 && rawHeaders[grpIdx])
        setColLabels((prev) => ({ ...prev, group: rawHeaders[grpIdx] }));
      if (dscIdx >= 0 && rawHeaders[dscIdx])
        setColLabels((prev) => ({ ...prev, description: rawHeaders[dscIdx] }));

      // Build a lookup by sampleId
      const csvMap: Record<string, MetaRow> = {};
      for (let i = 1; i < lines.length; i++) {
        const cells  = parseRow(lines[i]);
        const sid    = sidIdx >= 0 ? cells[sidIdx] : "";
        if (!sid) continue;
        const row: MetaRow = {
          sampleId:    sid,
          group:       grpIdx >= 0 ? (cells[grpIdx] ?? "") : "",
          description: dscIdx >= 0 ? (cells[dscIdx] ?? "") : "",
        };
        extraKeys.forEach((k, ki) => {
          const origIdx = extraIdxs[ki];
          row[k] = cells[origIdx] ?? "";
        });
        csvMap[sid] = row;
      }

      // Merge into existing rows (keep sampleId order from FASTQ names)
      let matched = 0;
      setRows((prev) => {
        const next = prev.map((r) => {
          const match = csvMap[r.sampleId];
          if (!match) return r;
          matched++;
          return { ...r, ...match, sampleId: r.sampleId };
        });
        onChange(next);
        return next;
      });

      const totalCols = 2 + newCustom.length + extraKeys.filter((k) => customCols.includes(k)).length;
      setCsvStatus(
        `✅ Metadata loaded from "${file.name}": ` +
        `${Object.keys(csvMap).length} sample(s) in CSV, ` +
        `${matched} matched · ${rawHeaders.length - 1} data column(s)`
      );
    };
    reader.readAsText(file);
    // Reset file input so same file can be re-uploaded
    e.target.value = "";
  };

  const allColKeys = ["group", "description", ...customCols];

  return (
    <div className="meta-editor">
      {/* ── CSV Upload confirmation banner ── */}
      {csvStatus && (
        <div className="meta-csv-status">
          {csvStatus}
          <button className="meta-csv-dismiss" onClick={() => setCsvStatus(null)}>✕</button>
        </div>
      )}

      {/* ── Toolbar ── */}
      <div className="meta-toolbar">
        <button className="btn-meta-action" onClick={() => fileRef.current?.click()}>
          📂 Upload CSV
        </button>
        <input ref={fileRef} type="file" accept=".csv" style={{ display:"none" }} onChange={uploadCSV} />

        <button className="btn-meta-action" onClick={downloadTemplate}>
          ⬇️ Download Template
        </button>

        <button className="btn-meta-action" onClick={() => setShowAddCol((v) => !v)}>
          ➕ Add Column
        </button>

        {showAddCol && (
          <div className="meta-addcol">
            <input
              className="meta-addcol-input"
              placeholder="Column name…"
              value={newColName}
              onChange={(e) => setNewColName(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && addColumn()}
              autoFocus
            />
            <button className="btn-meta-ok" onClick={addColumn}>OK</button>
          </div>
        )}
      </div>

      {/* ── Table ── */}
      <div className="meta-table-wrap">
        <table className="meta-table">
          <thead>
            <tr>
              <th className="meta-th meta-th-sid">Sample ID</th>
              {allColKeys.map((colKey) => (
                <th key={colKey} className="meta-th">
                  {editingHeader === colKey ? (
                    <input
                      className="meta-header-input"
                      value={editingHeaderVal}
                      autoFocus
                      onChange={(e) => setEditingHeaderVal(e.target.value)}
                      onBlur={() => commitHeaderRename(colKey)}
                      onKeyDown={(e) => {
                        if (e.key === "Enter") commitHeaderRename(colKey);
                        if (e.key === "Escape") { setEditingHeader(null); setEditingHeaderVal(""); }
                      }}
                    />
                  ) : (
                    <span
                      className="meta-th-label"
                      title="Click to rename this column"
                      onClick={() => {
                        setEditingHeader(colKey);
                        setEditingHeaderVal(colLabels[colKey] || colKey);
                      }}
                    >
                      {colLabels[colKey] || colKey}
                      <span className="meta-th-edit-hint">✏️</span>
                    </span>
                  )}
                  {customCols.includes(colKey) && editingHeader !== colKey && (
                    <button className="meta-del-col" onClick={() => removeColumn(colKey)} title="Remove column">✕</button>
                  )}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((row, idx) => (
              <tr key={row.sampleId} className={idx % 2 === 0 ? "meta-row-even" : ""}>
                <td className="meta-td meta-sid">{row.sampleId}</td>
                <td className="meta-td">
                  <input
                    className="meta-input"
                    placeholder="e.g. Control"
                    value={row.group}
                    onChange={(e) => updateCell(idx, "group", e.target.value)}
                  />
                </td>
                <td className="meta-td">
                  <input
                    className="meta-input"
                    placeholder="Optional description"
                    value={row.description}
                    onChange={(e) => updateCell(idx, "description", e.target.value)}
                  />
                </td>
                {customCols.map((colKey) => (
                  <td key={colKey} className="meta-td">
                    <input
                      className="meta-input"
                      value={row[colKey] ?? ""}
                      onChange={(e) => updateCell(idx, colKey, e.target.value)}
                    />
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <p className="meta-hint">
        💡 <strong>Group</strong> ใช้สำหรับเปรียบเทียบ beta diversity (PERMANOVA) และสี PCoA
        — ต้องมีอย่างน้อย 2 กลุ่มและ 3 samples ขึ้นไป
        · <em>คลิกชื่อคอลัมน์เพื่อเปลี่ยนชื่อ</em>
      </p>
    </div>
  );
}
