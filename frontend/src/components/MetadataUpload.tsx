/**
 * MetadataUpload — QIIME2 metadata TSV upload & group-column picker
 *
 * Expected file format (QIIME2 metadata):
 *   sample-id  treatment  timepoint  ...
 *   #q2:types  categorical  numeric  ...
 *   S1         control      0        ...
 *   S2         treated      1        ...
 *
 * Props:
 *   onMetadataReady(path, columns, groupCol)
 *     — called after successful upload with server path and available columns
 *   onGroupColChange(col)
 *     — called when the user picks a different group column
 */

import { useRef, useState } from "react";

const API = "";   // relative — proxied by Vite / served by FastAPI

export interface MetadataInfo {
  path:       string;
  filename:   string;
  columns:    string[];
  n_samples:  number;
  sample_ids: string[];
}

interface Props {
  onMetadataReady?: (info: MetadataInfo, groupCol: string) => void;
  onGroupColChange?: (col: string) => void;
  /** If true, the upload widget is rendered inline (no card wrapper) */
  inline?: boolean;
}

export default function MetadataUpload({ onMetadataReady, onGroupColChange, inline }: Props) {
  const fileRef   = useRef<HTMLInputElement>(null);
  const [info,     setInfo]     = useState<MetadataInfo | null>(null);
  const [groupCol, setGroupCol] = useState<string>("");
  const [loading,  setLoading]  = useState(false);
  const [error,    setError]    = useState<string | null>(null);
  const [warnings, setWarnings] = useState<string[]>([]);
  const [dragging, setDragging] = useState(false);

  // ── Upload handler ──────────────────────────────────────────────────────────
  const uploadFile = async (file: File) => {
    if (!file.name.match(/\.(tsv|txt|csv)$/i)) {
      setError("Please upload a .tsv or .txt metadata file");
      return;
    }
    setLoading(true);
    setError(null);
    setWarnings([]);
    setInfo(null);

    const fd = new FormData();
    fd.append("file", file);

    try {
      const res = await fetch(`${API}/qiime2/upload/metadata`, {
        method: "POST",
        body: fd,
      }).then(r => r.json());

      if (!res.valid && res.issues?.length) {
        setError(res.issues.join(" | "));
        setLoading(false);
        return;
      }

      if (res.warnings?.length) setWarnings(res.warnings);

      const metaInfo: MetadataInfo = {
        path:       res.path,
        filename:   res.filename,
        columns:    (res.columns || []).filter((c: string) => !c.startsWith("#")),
        n_samples:  res.n_samples || 0,
        sample_ids: res.sample_ids || [],
      };
      setInfo(metaInfo);

      // Auto-select first non-sample-id column as group
      const dataCols = metaInfo.columns.filter(
        c => !["sample-id","#SampleID","id"].includes(c.toLowerCase())
      );
      const defaultGroup = dataCols[0] || "";
      setGroupCol(defaultGroup);
      onMetadataReady?.(metaInfo, defaultGroup);

    } catch (e: any) {
      setError("Upload failed — could not reach backend");
    }
    setLoading(false);
  };

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) uploadFile(file);
  };

  // ── Drag & Drop ─────────────────────────────────────────────────────────────
  const onDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setDragging(false);
    const file = e.dataTransfer.files?.[0];
    if (file) uploadFile(file);
  };

  // ── Group column change ─────────────────────────────────────────────────────
  const handleGroupChange = (col: string) => {
    setGroupCol(col);
    onGroupColChange?.(col);
    if (info) onMetadataReady?.(info, col);
  };

  // ── UI helpers ──────────────────────────────────────────────────────────────
  const dataCols = (info?.columns || []).filter(
    c => !["sample-id","#SampleID","id"].includes(c.toLowerCase())
  );

  const content = (
    <div className="metadata-upload">

      {/* Drop zone */}
      {!info && (
        <div
          className={`metadata-dropzone ${dragging ? "metadata-dropzone--active" : ""}`}
          onDragEnter={e => { e.preventDefault(); setDragging(true); }}
          onDragLeave={e => { e.preventDefault(); setDragging(false); }}
          onDragOver={e => e.preventDefault()}
          onDrop={onDrop}
          onClick={() => fileRef.current?.click()}
        >
          <input
            ref={fileRef}
            type="file"
            accept=".tsv,.txt,.csv"
            style={{ display: "none" }}
            onChange={handleFileChange}
          />
          {loading ? (
            <div className="metadata-loading">
              <span className="metadata-spinner" />
              Uploading & validating…
            </div>
          ) : (
            <>
              <div className="metadata-dropzone__icon">📋</div>
              <div className="metadata-dropzone__text">
                Drop metadata TSV here or <span className="link">click to browse</span>
              </div>
              <div className="metadata-dropzone__hint">
                QIIME2 metadata format: sample-id + group columns
              </div>
            </>
          )}
        </div>
      )}

      {/* Error */}
      {error && (
        <div className="metadata-error">
          ❌ {error}
          <button
            className="metadata-clear-btn"
            onClick={() => { setError(null); setInfo(null); }}
          >
            Try again
          </button>
        </div>
      )}

      {/* Warnings */}
      {warnings.length > 0 && !error && (
        <div className="metadata-warnings">
          {warnings.map((w, i) => <div key={i}>⚠️ {w}</div>)}
        </div>
      )}

      {/* Success — summary + group picker */}
      {info && (
        <div className="metadata-success">
          <div className="metadata-file-row">
            <span className="metadata-file-icon">📋</span>
            <div className="metadata-file-info">
              <div className="metadata-file-name">{info.filename}</div>
              <div className="metadata-file-meta">
                {info.n_samples} samples · {dataCols.length} metadata columns
              </div>
            </div>
            <button
              className="metadata-remove-btn"
              title="Remove metadata"
              onClick={() => { setInfo(null); setGroupCol(""); setWarnings([]); setError(null); }}
            >
              ✕
            </button>
          </div>

          {/* Group column selector */}
          {dataCols.length > 0 && (
            <div className="metadata-group-row">
              <label className="metadata-group-label">
                Group column
                <span className="metadata-group-hint">
                  (used for diversity & ANCOMBC2)
                </span>
              </label>
              <select
                className="metadata-group-select"
                value={groupCol}
                onChange={e => handleGroupChange(e.target.value)}
              >
                {dataCols.map(c => (
                  <option key={c} value={c}>{c}</option>
                ))}
              </select>
            </div>
          )}

          {/* Sample preview */}
          {info.sample_ids.length > 0 && (
            <details className="metadata-samples-details">
              <summary className="metadata-samples-summary">
                Sample IDs ({info.n_samples} total)
              </summary>
              <div className="metadata-samples-list">
                {info.sample_ids.slice(0, 20).join(", ")}
                {info.n_samples > 20 && ` … and ${info.n_samples - 20} more`}
              </div>
            </details>
          )}

          {/* Format hint */}
          <div className="metadata-format-hint">
            ✅ Metadata ready — <strong>{groupCol || dataCols[0] || "—"}</strong> column will be used for grouping
          </div>
        </div>
      )}

      {/* Format guide */}
      {!info && !loading && !error && (
        <details className="metadata-guide">
          <summary>Metadata format guide</summary>
          <pre className="metadata-format-example">{`sample-id\ttreatment\ttimepoint
#q2:types\tcategorical\tnumeric
S1\tcontrol\t0
S2\tcontrol\t0
S3\ttreated\t1
S4\ttreated\t1`}</pre>
          <p style={{ fontSize: 12, color: "#6b7280", marginTop: 4 }}>
            First column must be <code>sample-id</code> matching your FASTQ filenames.
            The <code>#q2:types</code> row is recommended but optional.
          </p>
        </details>
      )}
    </div>
  );

  if (inline) return content;

  return (
    <div className="metadata-upload-card">
      <div className="metadata-upload-card__header">
        <span className="metadata-upload-card__icon">📊</span>
        <div>
          <div className="metadata-upload-card__title">Metadata</div>
          <div className="metadata-upload-card__subtitle">
            Optional — enables diversity grouping & differential abundance
          </div>
        </div>
      </div>
      {content}
    </div>
  );
}
