import { useState, useEffect } from "react";

interface DbFile {
  name:     string;
  path:     string;
  rel_path: string;
  size_mb:  number;
}

interface Props {
  onSelect: (path: string) => void;
  onClose:  () => void;
}

export default function DbFileBrowser({ onSelect, onClose }: Props) {
  const [files,    setFiles]    = useState<DbFile[]>([]);
  const [dbDir,    setDbDir]    = useState("");
  const [loading,  setLoading]  = useState(true);
  const [error,    setError]    = useState("");
  const [search,   setSearch]   = useState("");
  const [selected, setSelected] = useState("");

  useEffect(() => {
    fetch("/databases/browse")
      .then(r => r.ok ? r.json() : Promise.reject(r.statusText))
      .then(data => {
        setFiles(data.files ?? []);
        setDbDir(data.db_dir ?? "");
        setLoading(false);
      })
      .catch(e => {
        setError("ไม่สามารถโหลดรายการไฟล์ได้: " + e);
        setLoading(false);
      });
  }, []);

  const filtered = files.filter(f =>
    f.name.toLowerCase().includes(search.toLowerCase()) ||
    f.rel_path.toLowerCase().includes(search.toLowerCase())
  );

  return (
    <div className="db-browser-overlay" onClick={onClose}>
      <div className="db-browser-modal" onClick={e => e.stopPropagation()}>

        {/* Header */}
        <div className="db-browser-header">
          <div>
            <h3 className="db-browser-title">📂 เลือก Database File</h3>
            {dbDir && (
              <div className="db-browser-path">
                <span>📁 </span><code>{dbDir}</code>
              </div>
            )}
          </div>
          <button className="db-browser-close" onClick={onClose}>✕</button>
        </div>

        {/* Hint */}
        <div className="db-browser-hint">
          วางไฟล์ <code>.fa.gz</code> / <code>.fasta.gz</code> ของคุณไว้ใน folder ด้านบน แล้วกด Refresh
          <button className="db-browser-refresh"
            onClick={() => { setLoading(true); setError("");
              fetch("/databases/browse").then(r => r.json()).then(d => {
                setFiles(d.files ?? []); setLoading(false);
              }).catch(() => setLoading(false));
            }}>
            🔄 Refresh
          </button>
        </div>

        {/* Search */}
        <div className="db-browser-search-wrap">
          <input
            className="db-browser-search"
            placeholder="🔍 ค้นหาชื่อไฟล์..."
            value={search}
            onChange={e => setSearch(e.target.value)}
            autoFocus
          />
        </div>

        {/* File list */}
        <div className="db-browser-list">
          {loading && <div className="db-browser-status">⏳ กำลังโหลด...</div>}
          {error   && <div className="db-browser-status db-browser-error">{error}</div>}
          {!loading && !error && filtered.length === 0 && (
            <div className="db-browser-status">
              {files.length === 0
                ? <>ไม่พบไฟล์ database — วางไฟล์ไว้ใน <code>{dbDir}</code> แล้วกด Refresh</>
                : "ไม่พบไฟล์ที่ตรงกับการค้นหา"}
            </div>
          )}
          {filtered.map(f => (
            <div key={f.path}
              className={`db-browser-item ${selected === f.path ? "selected" : ""}`}
              onClick={() => setSelected(f.path)}>
              <div className="db-browser-item-icon">🗄️</div>
              <div className="db-browser-item-info">
                <div className="db-browser-item-name">{f.name}</div>
                <div className="db-browser-item-meta">
                  <span className="db-browser-item-rel">{f.rel_path}</span>
                  <span className="db-browser-item-size">{f.size_mb} MB</span>
                </div>
              </div>
              {selected === f.path && <div className="db-browser-item-check">✓</div>}
            </div>
          ))}
        </div>

        {/* Footer */}
        <div className="db-browser-footer">
          <button className="db-browser-cancel" onClick={onClose}>ยกเลิก</button>
          <button
            className="db-browser-select"
            disabled={!selected}
            onClick={() => { if (selected) { onSelect(selected); onClose(); } }}>
            เลือกไฟล์นี้
          </button>
        </div>
      </div>
    </div>
  );
}
