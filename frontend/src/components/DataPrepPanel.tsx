import { useState } from "react";
import axios from "axios";

/**
 * Data preparation — read orientation check & repair.
 *
 * Paired-end amplicon libraries sometimes carry a large share of read pairs in
 * the opposite orientation (R1 starting at the reverse primer instead of the
 * forward one). DADA2 cannot merge those — R1 and revcomp(R2) are biologically
 * misaligned — and cutadapt's paired-end mode only searches one fixed
 * orientation, so it leaves them untrimmed rather than fixing them. The only
 * symptom is a merge rate far below expectation with no error anywhere, which
 * is close to impossible to diagnose from the UI. This panel makes the problem
 * visible before committing to a full run.
 */

interface SampleReport {
  sample?:             string;
  total_pairs?:        number;
  already_correct?:    number;
  flipped?:            number;
  ambiguous?:          number;
  written_pairs?:      number;
  pct_already_correct?: number;
  pct_flipped?:        number;
  pct_ambiguous?:      number;
  error?:              string;
}

interface CheckResult {
  paired:              boolean;
  samples:             SampleReport[];
  max_pct_flipped?:    number;
  needs_reorientation?: boolean;
  recommendation?:     string;
  message?:            string;
}

interface FixResult {
  samples:             SampleReport[];
  replaced_originals:  boolean;
  originals_kept_in:   string | null;
  total_pairs:         number;
  flipped_pairs:       number;
  written_pairs:       number;
  summary:             string;
}

interface Props {
  apiBase:   string;
  jobId:     string | null;
  primerF:   string;
  primerR:   string;
  isPaired:  boolean;
  onFilesChanged?: () => void;
}

export default function DataPrepPanel({
  apiBase, jobId, primerF, primerR, isPaired, onFilesChanged,
}: Props) {
  const [checking,  setChecking]  = useState(false);
  const [fixing,    setFixing]    = useState(false);
  const [check,     setCheck]     = useState<CheckResult | null>(null);
  const [fixed,     setFixed]     = useState<FixResult | null>(null);
  const [err,       setErr]       = useState<string>("");
  const [ambiguous, setAmbiguous] = useState<"keep" | "drop">("keep");

  const ready = Boolean(jobId && primerF && primerR && isPaired);

  const runCheck = async () => {
    if (!jobId) return;
    setChecking(true); setErr(""); setFixed(null);
    try {
      const res = await axios.post(`${apiBase}/prep/${jobId}/orientation-check`, {
        primer_f: primerF, primer_r: primerR,
      });
      setCheck(res.data);
    } catch (e: any) {
      setErr(e?.response?.data?.error ?? "Orientation check failed.");
    }
    setChecking(false);
  };

  const runFix = async () => {
    if (!jobId) return;
    setFixing(true); setErr("");
    try {
      const res = await axios.post(`${apiBase}/prep/${jobId}/reorient`, {
        primer_f: primerF, primer_r: primerR,
        ambiguous, replace: true,
      });
      setFixed(res.data);
      setCheck(null);
      onFilesChanged?.();
    } catch (e: any) {
      setErr(e?.response?.data?.error ?? "Reorientation failed.");
    }
    setFixing(false);
  };

  const pct = check?.max_pct_flipped ?? 0;
  const bad = Boolean(check?.needs_reorientation);

  return (
    <div style={{ marginTop: 4 }}>
      <p className="param-hint" style={{ marginBottom: 10 }}>
        Checks whether your R1/R2 files are in a consistent orientation. Pairs
        sequenced the other way round cannot be merged by DADA2 and are lost
        at the merge step without any error message — this finds them first.
      </p>

      {!isPaired && (
        <div style={{ fontSize: 13, color: "#64748b" }}>
          Orientation repair applies to paired-end data only.
        </div>
      )}

      {isPaired && !jobId && (
        <div style={{ fontSize: 13, color: "#f59e0b" }}>
          ⬆️ Upload your FASTQ files first.
        </div>
      )}

      {isPaired && jobId && (!primerF || !primerR) && (
        <div style={{ fontSize: 13, color: "#f59e0b" }}>
          ⚠ Set both primers in Pipeline Settings → Step 1 first — the check
          works by looking for them at the start of each read.
        </div>
      )}

      {ready && (
        <>
          <div style={{ display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
            <button type="button" className="btn-advanced-toggle"
              style={{ width: "auto", padding: "8px 16px" }}
              onClick={runCheck} disabled={checking || fixing}>
              {checking ? "⏳ Checking…" : "🔍 Check read orientation"}
            </button>
            <span style={{ fontSize: 12, color: "#64748b", fontFamily: "monospace" }}>
              {primerF.slice(0, 12)}… / {primerR.slice(0, 12)}…
            </span>
          </div>

          {err && (
            <div style={{ marginTop: 10, fontSize: 13, color: "#ef4444" }}>⚠ {err}</div>
          )}

          {check && !check.paired && (
            <div style={{ marginTop: 10, fontSize: 13, color: "#94a3b8" }}>
              {check.message}
            </div>
          )}

          {check && check.paired && (
            <div style={{
              marginTop: 12, padding: 12, borderRadius: 8,
              border: `1px solid ${bad ? "#f59e0b" : "#22c55e"}`,
              background: bad ? "#2a1f06" : "#052e16",
            }}>
              <div style={{
                fontWeight: 700, marginBottom: 6,
                color: bad ? "#fbbf24" : "#4ade80",
              }}>
                {bad ? `⚠ Mixed orientation — ${pct.toFixed(1)}% of pairs are flipped`
                     : "✅ Orientation is consistent"}
              </div>
              <div style={{ fontSize: 13, color: "#cbd5e1", marginBottom: 8 }}>
                {check.recommendation}
              </div>

              <table style={{ width: "100%", fontSize: 12, borderCollapse: "collapse" }}>
                <thead>
                  <tr style={{ color: "#94a3b8", textAlign: "right" }}>
                    <th style={{ textAlign: "left", padding: "4px 6px" }}>Sample</th>
                    <th style={{ padding: "4px 6px" }}>Pairs scanned</th>
                    <th style={{ padding: "4px 6px" }}>Forward</th>
                    <th style={{ padding: "4px 6px" }}>Flipped</th>
                    <th style={{ padding: "4px 6px" }}>Ambiguous</th>
                  </tr>
                </thead>
                <tbody>
                  {check.samples.map((s, i) => (
                    <tr key={i} style={{ borderTop: "1px solid #334155", textAlign: "right" }}>
                      <td style={{ textAlign: "left", padding: "4px 6px", color: "#e2e8f0" }}>
                        {s.sample}
                      </td>
                      {s.error ? (
                        <td colSpan={4} style={{ padding: "4px 6px", color: "#ef4444" }}>
                          {s.error}
                        </td>
                      ) : (
                        <>
                          <td style={{ padding: "4px 6px", color: "#94a3b8" }}>
                            {s.total_pairs?.toLocaleString()}
                          </td>
                          <td style={{ padding: "4px 6px", color: "#4ade80" }}>
                            {s.pct_already_correct?.toFixed(1)}%
                          </td>
                          <td style={{
                            padding: "4px 6px", fontWeight: 700,
                            color: (s.pct_flipped ?? 0) >= 5 ? "#fbbf24" : "#94a3b8",
                          }}>
                            {s.pct_flipped?.toFixed(1)}%
                          </td>
                          <td style={{ padding: "4px 6px", color: "#94a3b8" }}>
                            {s.pct_ambiguous?.toFixed(1)}%
                          </td>
                        </>
                      )}
                    </tr>
                  ))}
                </tbody>
              </table>

              {bad && (
                <div style={{ marginTop: 12, borderTop: "1px solid #334155", paddingTop: 10 }}>
                  <div style={{ display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap" }}>
                    <label style={{ fontSize: 12, color: "#94a3b8" }}>
                      Pairs matching neither primer:{" "}
                      <select value={ambiguous}
                        onChange={e => setAmbiguous(e.target.value as "keep" | "drop")}
                        style={{
                          background: "#1e293b", color: "#e2e8f0",
                          border: "1px solid #334155", borderRadius: 4, padding: "3px 6px",
                        }}>
                        <option value="keep">Keep</option>
                        <option value="drop">Drop</option>
                      </select>
                    </label>
                    <button type="button" className="btn-submit"
                      style={{ padding: "8px 16px" }}
                      onClick={runFix} disabled={fixing}>
                      {fixing ? "⏳ Reorienting…" : "🔧 Fix orientation now"}
                    </button>
                  </div>
                  <div style={{ fontSize: 11, color: "#64748b", marginTop: 6 }}>
                    Your original files are kept in a <code>_original</code> folder
                    inside the job; the repaired files take over the same filenames
                    so the run picks them up automatically.
                  </div>
                </div>
              )}
            </div>
          )}

          {fixed && (
            <div style={{
              marginTop: 12, padding: 12, borderRadius: 8,
              border: "1px solid #22c55e", background: "#052e16",
            }}>
              <div style={{ fontWeight: 700, color: "#4ade80", marginBottom: 6 }}>
                ✅ Reads reoriented
              </div>
              <div style={{ fontSize: 13, color: "#cbd5e1" }}>{fixed.summary}</div>
              {fixed.replaced_originals && (
                <div style={{ fontSize: 12, color: "#94a3b8", marginTop: 6 }}>
                  Originals kept in <code>_original/</code>. You can run the
                  pipeline now — it will use the repaired files.
                </div>
              )}
              <button type="button" className="btn-advanced-toggle"
                style={{ width: "auto", padding: "6px 12px", marginTop: 10 }}
                onClick={runCheck} disabled={checking}>
                {checking ? "⏳ Verifying…" : "🔁 Verify the result"}
              </button>
            </div>
          )}
        </>
      )}
    </div>
  );
}
