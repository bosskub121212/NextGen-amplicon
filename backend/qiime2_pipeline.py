"""
NextGen-Amplicon — QIIME2 Pipeline Backend
==========================================
Handles full QIIME2 pipeline for all marker types:
  16S / 12S / ITS1 / ITS2 / COX1 / 18S-nema / PacBio

Pipeline steps per marker:
  Step 1  — Import (manifest TSV → .qza)
  Step 2  — Demux summarize (quality check)
  Step 3  — Cutadapt (primer trim)
  Step 4  — DADA2 denoise (paired / single-end for PacBio)
  Step 5  — Taxonomy (classify-sklearn)
  Step 6  — Phylogeny (mafft-fasttree) [16S / 12S / PacBio only]
  Step 7  — Diversity (core-metrics + alpha/beta significance)
  Step 8  — Differential abundance (ANCOM-BC)
  Step 9  — Export (BIOM, TSV, FASTA, tree)
  Step 10 — R Visualization (phyloseq, ggplot2, vegan, ANCOMBC2, FUNGuildR, LULU)
"""

import asyncio
import json
import os
import shutil
import subprocess
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, BackgroundTasks, UploadFile, File, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel

router = APIRouter(prefix="/qiime2", tags=["qiime2"])

# ─── Config ──────────────────────────────────────────────────────────────────
QIIME2_ENV    = os.getenv("QIIME2_ENV", "qiime2-amplicon-2024.5")
BASE_DIR      = Path(__file__).parent
RESULTS_DIR   = BASE_DIR / "results"
UPLOAD_DIR    = BASE_DIR / "uploads"
CLASSIFIER_DIR = BASE_DIR / "classifiers"
R_SCRIPTS_DIR = BASE_DIR / "r_scripts"

RESULTS_DIR.mkdir(exist_ok=True)
UPLOAD_DIR.mkdir(exist_ok=True)
CLASSIFIER_DIR.mkdir(exist_ok=True)

# ─── Classifier auto-detect map ──────────────────────────────────────────────
CLASSIFIER_PATTERNS = {
    "16S":       ["silva_16S_classifier.qza", "silva-138*515-806*.qza", "silva*16S*.qza"],
    "12S":       ["silva_12S_classifier.qza", "pr2*12S*.qza", "silva*12S*.qza"],
    "ITS1":      ["unite_ITS_classifier.qza", "unite*ITS1*.qza", "unite*.qza"],
    "ITS2":      ["unite_ITS_classifier.qza", "unite*ITS2*.qza", "unite*.qza"],
    "COX1":      ["midori2_COX1_classifier.qza", "midori*.qza", "cox1*.qza"],
    "18S-nema":  ["nemabase_18S_classifier.qza", "pr2*18S*.qza", "silva*18S*.qza"],
    "PacBio":    ["silva_16S_full_classifier.qza", "silva*full*.qza", "silva*16S*.qza"],
}

# ─── Marker settings ─────────────────────────────────────────────────────────
MARKER_CONFIG = {
    "16S": {
        "seq_type": "paired", "use_phylo": True,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
    },
    "12S": {
        "seq_type": "paired", "use_phylo": True,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
    },
    "ITS1": {
        "seq_type": "paired", "use_phylo": False,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
        "no_trunc": True,
    },
    "ITS2": {
        "seq_type": "paired", "use_phylo": False,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
        "no_trunc": True,
    },
    "COX1": {
        "seq_type": "paired", "use_phylo": False,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
    },
    "18S-nema": {
        "seq_type": "paired", "use_phylo": False,
        "artifact_type": "SampleData[PairedEndSequencesWithQuality]",
        "input_format":  "PairedEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-paired",
    },
    "PacBio": {
        "seq_type": "single", "use_phylo": True,
        "artifact_type": "SampleData[SequencesWithQuality]",
        "input_format":  "SingleEndFastqManifestPhred33V2",
        "dada2_method":  "denoise-single",
    },
}

# ─── Job store (shared with main.py) ─────────────────────────────────────────
# Imported from main.py at runtime
_jobs: dict = {}
_log_queues: dict = {}

def set_shared_stores(jobs: dict, log_queues: dict):
    global _jobs, _log_queues
    _jobs = jobs
    _log_queues = log_queues


# ══════════════════════════════════════════════════════════════════════════════
# MODELS
# ══════════════════════════════════════════════════════════════════════════════

class Q2RunParams(BaseModel):
    job_id:          str
    marker:          str = "16S"
    manifest_path:   str
    metadata_path:   Optional[str] = None
    output_dir:      str
    # DADA2
    trim_left_f:     int   = 0
    trim_left_r:     int   = 0
    trunc_len_f:     int   = 240
    trunc_len_r:     int   = 200
    max_ee_f:        float = 2.0
    max_ee_r:        float = 2.0
    chimera_method:  str   = "consensus"
    n_threads:       int   = 4
    # Primers
    primer_f:        str   = ""
    primer_r:        str   = ""
    # Taxonomy
    classifier_path: str   = ""
    confidence:      float = 0.7
    # Custom classifier training
    custom_classifier_mode: str = "default"  # default | train | upload
    custom_classifier_path: str = ""          # path to .qza (upload mode)
    train_amplicon_min_len: int = 200
    train_amplicon_max_len: int = 600
    # Diversity
    sampling_depth:  int   = 10000
    group_col:       str   = "treatment"
    # ANCOM-BC
    run_diffabund:   bool  = True
    diffabund_formula: str = ""
    # Phylogeny
    n_threads_phylo: int   = 4
    # R visualization
    run_r_viz:       bool  = True


class Q2StartParams(BaseModel):
    marker:          str   = "16S"
    manifest_path:   str
    metadata_path:   Optional[str] = None
    output_dir:      Optional[str] = None
    trim_left_f:     int   = 0
    trim_left_r:     int   = 0
    trunc_len_f:     int   = 240
    trunc_len_r:     int   = 200
    max_ee_f:        float = 2.0
    max_ee_r:        float = 2.0
    chimera_method:  str   = "consensus"
    n_threads:       int   = 4
    primer_f:        str   = ""
    primer_r:        str   = ""
    classifier_path: str   = ""
    confidence:      float = 0.7
    custom_classifier_mode: str = "default"
    custom_classifier_path: str = ""
    train_amplicon_min_len: int = 200
    train_amplicon_max_len: int = 600
    sampling_depth:  int   = 10000
    group_col:       str   = "treatment"
    run_diffabund:   bool  = True
    diffabund_formula: str = ""
    run_r_viz:       bool  = True


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def q2(cmd: str) -> str:
    """Prefix command with conda run for QIIME2 env."""
    return f"conda run -n {QIIME2_ENV} --no-capture-output {cmd}"


def find_classifier(marker: str, override: str = "") -> str:
    if override and Path(override).exists():
        return override
    patterns = CLASSIFIER_PATTERNS.get(marker, [])
    import glob
    for pat in patterns:
        matches = glob.glob(str(CLASSIFIER_DIR / pat))
        if matches:
            return matches[0]
    return ""


async def train_custom_classifier(
    primer_f: str,
    primer_r: str,
    marker: str,
    min_len: int,
    max_len: int,
    job_id: str,
    silva_ref_qza: str = "",
    silva_tax_qza: str = "",
) -> str:
    """
    Train a Naive Bayes classifier for a specific primer region.
    Steps:
      1. qiime feature-classifier extract-reads  (primer trimming on SILVA)
      2. qiime feature-classifier fit-classifier-naive-bayes
    Returns path to trained classifier .qza, or "" on failure.
    """
    import glob as _glob
    CLASSIFIER_DIR.mkdir(parents=True, exist_ok=True)

    # Auto-locate SILVA reference sequences and taxonomy if not provided
    if not silva_ref_qza:
        candidates = _glob.glob(str(CLASSIFIER_DIR.parent / "databases" / "silva*seqs*.qza"))
        candidates += _glob.glob(str(CLASSIFIER_DIR.parent / "databases" / "*silva*seqs*.qza"))
        if candidates:
            silva_ref_qza = candidates[0]
    if not silva_tax_qza:
        candidates = _glob.glob(str(CLASSIFIER_DIR.parent / "databases" / "silva*tax*.qza"))
        candidates += _glob.glob(str(CLASSIFIER_DIR.parent / "databases" / "*silva*tax*.qza"))
        if candidates:
            silva_tax_qza = candidates[0]

    if not silva_ref_qza or not Path(silva_ref_qza).exists():
        q = _log_queues.get(job_id)
        if q:
            await q.put({"type": "log", "msg": "⚠️  train_custom_classifier: SILVA seqs .qza not found — skipping training"})
        return ""
    if not silva_tax_qza or not Path(silva_tax_qza).exists():
        q = _log_queues.get(job_id)
        if q:
            await q.put({"type": "log", "msg": "⚠️  train_custom_classifier: SILVA taxonomy .qza not found — skipping training"})
        return ""

    # Build a unique name for this primer combination
    import hashlib
    pkey = hashlib.md5(f"{primer_f}_{primer_r}_{min_len}_{max_len}".encode()).hexdigest()[:8]
    out_reads_qza  = str(CLASSIFIER_DIR / f"ref_reads_{marker}_{pkey}.qza")
    out_clf_qza    = str(CLASSIFIER_DIR / f"classifier_{marker}_{pkey}.qza")

    if Path(out_clf_qza).exists():
        q = _log_queues.get(job_id)
        if q:
            await q.put({"type": "log", "msg": f"✅ Custom classifier already exists: {out_clf_qza}"})
        return out_clf_qza

    q = _log_queues.get(job_id)
    if q:
        await q.put({"type": "log", "msg": f"🧬 Training custom classifier for {marker} region (primers: {primer_f[:12]}… / {primer_r[:12]}…)"})
        await q.put({"type": "log", "msg": "   Step 1/2: extract-reads (may take 10–20 min)…"})

    rc = await stream_cmd(
        q2(f"qiime feature-classifier extract-reads"
           f" --i-sequences {silva_ref_qza}"
           f" --p-f-primer {primer_f}"
           f" --p-r-primer {primer_r}"
           f" --p-min-length {min_len}"
           f" --p-max-length {max_len}"
           f" --p-n-jobs -1"
           f" --o-reads {out_reads_qza}"
           f" --quiet"),
        job_id,
    )
    if rc != 0:
        if q:
            await q.put({"type": "log", "msg": f"❌ extract-reads failed (rc={rc}) — falling back to default classifier"})
        return ""

    if q:
        await q.put({"type": "log", "msg": "   Step 2/2: fit-classifier-naive-bayes (may take 10–20 min)…"})

    rc = await stream_cmd(
        q2(f"qiime feature-classifier fit-classifier-naive-bayes"
           f" --i-reference-reads {out_reads_qza}"
           f" --i-reference-taxonomy {silva_tax_qza}"
           f" --o-classifier {out_clf_qza}"
           f" --quiet"),
        job_id,
    )
    if rc != 0:
        if q:
            await q.put({"type": "log", "msg": f"❌ fit-classifier-naive-bayes failed (rc={rc}) — falling back to default classifier"})
        return ""

    if q:
        await q.put({"type": "log", "msg": f"✅ Custom classifier trained: {out_clf_qza}"})
    return out_clf_qza


async def stream_cmd(cmd: str, job_id: str) -> int:
    q = _log_queues.get(job_id)
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        env=os.environ.copy(),
    )
    async for line in proc.stdout:
        text = line.decode(errors="replace").rstrip()
        if text and q:
            await q.put({"type": "log", "msg": text})
    await proc.wait()
    return proc.returncode


async def run_step(cmd: str, job_id: str, step_name: str,
                   step_num: int, total: int, must_succeed: bool = True) -> bool:
    q = _log_queues.get(job_id)
    pct = int((step_num - 1) / total * 100)
    if q:
        await q.put({"type": "progress", "pct": pct,
                     "label": f"Step {step_num}/{total} — {step_name}"})
        await q.put({"type": "log", "msg": f"\n{'='*60}"})
        await q.put({"type": "log", "msg": f"[STEP {step_num}/{total}] {step_name}"})
        await q.put({"type": "log", "msg": f"[CMD] {cmd}"})
        await q.put({"type": "log", "msg": f"{'='*60}"})

    rc = await stream_cmd(cmd, job_id)
    success = rc == 0

    if q:
        status = "✓ Done" if success else f"✗ Failed (rc={rc})"
        await q.put({"type": "log", "msg": f"[{status}] {step_name}"})

    if must_succeed and not success:
        if _jobs.get(job_id):
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error_step"] = step_name
        if q:
            await q.put({"type": "pipeline_done", "success": False,
                         "error": f"Step failed: {step_name}"})
        raise RuntimeError(f"QIIME2 step failed: {step_name}")

    return success


# ══════════════════════════════════════════════════════════════════════════════
# METADATA UPLOAD
# ══════════════════════════════════════════════════════════════════════════════

@router.post("/upload/metadata")
async def upload_metadata(file: UploadFile = File(...)):
    """Upload QIIME2 metadata TSV and validate format."""
    tmp = UPLOAD_DIR / f"metadata_{uuid.uuid4().hex[:8]}.tsv"
    content = await file.read()
    tmp.write_bytes(content)

    # Basic validation
    issues, warnings, columns, sample_ids = [], [], [], []
    try:
        lines = tmp.read_text(errors="replace").splitlines()
        lines = [l for l in lines if l.strip()]
        if not lines:
            issues.append("File is empty")
        else:
            # Auto-detect delimiter (TSV vs CSV)
            delim = "\t" if "\t" in lines[0] else ","
            header = [c.strip().strip('"').strip("'") for c in lines[0].split(delim)]
            columns = header
            first_col = header[0].lower()
            if first_col not in ("sample-id", "#sampleid", "sampleid", "id", "#sample-id"):
                warnings.append(f"First column should be 'sample-id', got '{header[0]}'")

            has_types = len(lines) > 1 and lines[1].startswith("#q2:types")
            if not has_types:
                warnings.append("No '#q2:types' row — recommended for QIIME2")

            data_start = 2 if has_types else 1
            for line in lines[data_start:]:
                if line.startswith("#"):
                    continue
                cols = line.split(delim)
                sid = cols[0].strip().strip('"').strip("'")
                if sid:
                    sample_ids.append(sid)
    except Exception as e:
        issues.append(str(e))

    return {
        "path": str(tmp),
        "filename": file.filename,
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
        "columns": columns,
        "n_samples": len(sample_ids),
        "sample_ids": sample_ids[:30],
    }


@router.post("/upload/manifest")
async def upload_manifest(file: UploadFile = File(...)):
    """Upload manifest TSV (sample-id, forward, reverse paths)."""
    tmp = UPLOAD_DIR / f"manifest_{uuid.uuid4().hex[:8]}.tsv"
    tmp.write_bytes(await file.read())

    issues, samples = [], []
    try:
        lines = tmp.read_text().splitlines()
        lines = [l for l in lines if l.strip() and not l.startswith("#")]
        if not lines:
            issues.append("Empty file")
        else:
            header = lines[0].split("\t")
            if "sample-id" not in header:
                issues.append("Missing 'sample-id' column")
            for line in lines[1:]:
                cols = line.split("\t")
                if len(cols) >= 2:
                    samples.append({
                        "id": cols[0],
                        "forward": cols[1] if len(cols) > 1 else "",
                        "reverse": cols[2] if len(cols) > 2 else "",
                    })
    except Exception as e:
        issues.append(str(e))

    return {
        "path": str(tmp),
        "filename": file.filename,
        "valid": len(issues) == 0,
        "issues": issues,
        "n_samples": len(samples),
        "samples": samples[:20],
    }


# ══════════════════════════════════════════════════════════════════════════════
# QIIME2 ENVIRONMENT CHECK
# ══════════════════════════════════════════════════════════════════════════════

@router.get("/env/check")
async def check_qiime2_env():
    """Check QIIME2 environment availability."""
    result = {
        "env_name": QIIME2_ENV,
        "found": False,
        "version": None,
        "plugins": {},
        "classifiers": {},
        "errors": [],
    }
    try:
        r = subprocess.run(
            f"conda run -n {QIIME2_ENV} --no-capture-output qiime --version",
            shell=True, capture_output=True, text=True, timeout=30
        )
        if r.returncode == 0:
            result["found"] = True
            result["version"] = (r.stdout + r.stderr).strip()
        else:
            result["errors"].append(f"QIIME2 not found in env '{QIIME2_ENV}'")
            result["errors"].append("Run: bash ~/r16s-app/setup_qiime2.sh")
    except Exception as e:
        result["errors"].append(str(e))

    # Check classifiers
    import glob
    for marker in CLASSIFIER_PATTERNS:
        clf = find_classifier(marker)
        result["classifiers"][marker] = {"found": bool(clf), "path": clf}

    return result


# ══════════════════════════════════════════════════════════════════════════════
# FULL PIPELINE — main entry point
# ══════════════════════════════════════════════════════════════════════════════

@router.post("/run")
async def start_qiime2_pipeline(params: Q2StartParams,
                                 background_tasks: BackgroundTasks):
    """Start full QIIME2 pipeline (all steps) in background."""
    job_id = str(uuid.uuid4())[:8]

    out_dir = params.output_dir or str(RESULTS_DIR / job_id)
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    _jobs[job_id] = {
        "id": job_id,
        "type": "qiime2",
        "marker": params.marker,
        "status": "running",
        "started": datetime.now().isoformat(),
        "output_dir": out_dir,
        "steps": {},
    }
    _log_queues[job_id] = asyncio.Queue()

    run_params = Q2RunParams(
        job_id=job_id,
        marker=params.marker,
        manifest_path=params.manifest_path,
        metadata_path=params.metadata_path,
        output_dir=out_dir,
        trim_left_f=params.trim_left_f,
        trim_left_r=params.trim_left_r,
        trunc_len_f=params.trunc_len_f,
        trunc_len_r=params.trunc_len_r,
        max_ee_f=params.max_ee_f,
        max_ee_r=params.max_ee_r,
        chimera_method=params.chimera_method,
        n_threads=params.n_threads,
        primer_f=params.primer_f,
        primer_r=params.primer_r,
        classifier_path=params.classifier_path,
        confidence=params.confidence,
        custom_classifier_mode=params.custom_classifier_mode,
        custom_classifier_path=params.custom_classifier_path,
        train_amplicon_min_len=params.train_amplicon_min_len,
        train_amplicon_max_len=params.train_amplicon_max_len,
        sampling_depth=params.sampling_depth,
        group_col=params.group_col,
        run_diffabund=params.run_diffabund,
        diffabund_formula=params.diffabund_formula,
        run_r_viz=params.run_r_viz,
    )

    background_tasks.add_task(_run_full_pipeline, run_params)
    return {"job_id": job_id, "output_dir": out_dir}


async def _run_full_pipeline(p: Q2RunParams):
    """Execute all QIIME2 pipeline steps sequentially."""
    q = _log_queues.get(p.job_id)
    cfg = MARKER_CONFIG.get(p.marker, MARKER_CONFIG["16S"])
    is_ITS = p.marker in ("ITS1", "ITS2")
    use_phylo = cfg["use_phylo"]

    # Count total steps
    TOTAL = 9
    if not use_phylo:
        TOTAL -= 1  # no phylogeny
    if not p.run_diffabund or not p.metadata_path:
        TOTAL -= 1  # no ANCOM-BC
    if p.run_r_viz:
        TOTAL += 1  # R viz

    step_n = [0]  # mutable counter

    def nxt():
        step_n[0] += 1
        return step_n[0]

    qza  = Path(p.output_dir) / "qza"
    qzv  = Path(p.output_dir) / "qzv"
    exp  = Path(p.output_dir) / "exported"
    for d in [qza, qzv, exp]:
        d.mkdir(parents=True, exist_ok=True)

    if q:
        await q.put({"type": "log",
                     "msg": f"=== NextGen-Amplicon QIIME2 Pipeline ===\n"
                            f"Marker   : {p.marker}\n"
                            f"Output   : {p.output_dir}\n"
                            f"Steps    : {TOTAL} total\n"})

    try:
        # ── Step 1: Import ────────────────────────────────────────────
        demux_qza = str(qza / "demux.qza")
        await run_step(
            q2(f"qiime tools import "
               f"--type '{cfg['artifact_type']}' "
               f"--input-path '{p.manifest_path}' "
               f"--input-format {cfg['input_format']} "
               f"--output-path '{demux_qza}'"),
            p.job_id, "Import data", nxt(), TOTAL
        )

        # ── Step 2: Demux summarize ───────────────────────────────────
        demux_qzv = str(qzv / "demux-summary.qzv")
        await run_step(
            q2(f"qiime demux summarize "
               f"--i-data '{demux_qza}' "
               f"--o-visualization '{demux_qzv}'"),
            p.job_id, "Demux summarize (quality check)", nxt(), TOTAL,
            must_succeed=False
        )

        # ── Step 3: Cutadapt ─────────────────────────────────────────
        trimmed_qza = demux_qza  # default: use raw if no primers
        if p.primer_f or p.primer_r:
            trimmed_qza = str(qza / "demux-trimmed.qza")
            primer_args = ""
            if p.primer_f:
                primer_args += f"--p-front-f '{p.primer_f}' "
            if p.primer_r:
                primer_args += f"--p-front-r '{p.primer_r}' "

            trim_cmd = (
                f"qiime cutadapt trim-paired "
                if cfg["seq_type"] == "paired"
                else "qiime cutadapt trim-single "
            )
            await run_step(
                q2(f"{trim_cmd} "
                   f"--i-demultiplexed-sequences '{demux_qza}' "
                   f"{primer_args}"
                   f"--p-error-rate 0.1 "
                   f"--p-overlap 10 "
                   f"--p-discard-untrimmed "
                   f"--p-cores {p.n_threads} "
                   f"--o-trimmed-sequences '{trimmed_qza}'"),
                p.job_id, "Cutadapt primer trimming", nxt(), TOTAL
            )

        # ── Step 4: DADA2 ────────────────────────────────────────────
        table_qza   = str(qza / "table.qza")
        repseqs_qza = str(qza / "rep-seqs.qza")
        stats_qza   = str(qza / "denoising-stats.qza")

        dada2_method = cfg["dada2_method"]

        # ITS: no truncation (variable length)
        trunc_f = 0 if is_ITS else p.trunc_len_f
        trunc_r = 0 if is_ITS else p.trunc_len_r

        rev_args = ""
        if cfg["seq_type"] == "paired":
            rev_args = (
                f"--p-trim-left-r {p.trim_left_r} "
                f"--p-trunc-len-r {trunc_r} "
                f"--p-max-ee-r {p.max_ee_r} "
            )

        await run_step(
            q2(f"qiime dada2 {dada2_method} "
               f"--i-demultiplexed-seqs '{trimmed_qza}' "
               f"--p-trim-left-f {p.trim_left_f} "
               f"--p-trunc-len-f {trunc_f} "
               f"--p-max-ee-f {p.max_ee_f} "
               f"{rev_args}"
               f"--p-chimera-method {p.chimera_method} "
               f"--p-n-threads {p.n_threads} "
               f"--o-table '{table_qza}' "
               f"--o-representative-sequences '{repseqs_qza}' "
               f"--o-denoising-stats '{stats_qza}'"),
            p.job_id, "DADA2 denoising & ASV inference", nxt(), TOTAL
        )

        # Visualize denoising stats
        await run_step(
            q2(f"qiime metadata tabulate "
               f"--m-input-file '{stats_qza}' "
               f"--o-visualization '{qzv}/denoising-stats.qzv'"),
            p.job_id, "DADA2 stats visualization", nxt(), TOTAL,
            must_succeed=False
        )

        # ── Step 5: Taxonomy ─────────────────────────────────────────
        taxonomy_qza = str(qza / "taxonomy.qza")

        # Resolve classifier: upload > train > default auto-detect
        if p.custom_classifier_mode == "upload" and p.custom_classifier_path and Path(p.custom_classifier_path).exists():
            clf = p.custom_classifier_path
        elif p.custom_classifier_mode == "train" and p.primer_f and p.primer_r:
            trained = await train_custom_classifier(
                primer_f=p.primer_f,
                primer_r=p.primer_r,
                marker=p.marker,
                min_len=p.train_amplicon_min_len,
                max_len=p.train_amplicon_max_len,
                job_id=p.job_id,
            )
            clf = trained if trained else find_classifier(p.marker, p.classifier_path)
        else:
            clf = find_classifier(p.marker, p.classifier_path)

        if clf:
            await run_step(
                q2(f"qiime feature-classifier classify-sklearn "
                   f"--i-classifier '{clf}' "
                   f"--i-reads '{repseqs_qza}' "
                   f"--p-confidence {p.confidence} "
                   f"--p-n-jobs {p.n_threads} "
                   f"--o-classification '{taxonomy_qza}'"),
                p.job_id, "Taxonomy classification", nxt(), TOTAL
            )

            # Taxa bar plot (needs metadata)
            if p.metadata_path and Path(p.metadata_path).exists():
                await run_step(
                    q2(f"qiime taxa barplot "
                       f"--i-table '{table_qza}' "
                       f"--i-taxonomy '{taxonomy_qza}' "
                       f"--m-metadata-file '{p.metadata_path}' "
                       f"--o-visualization '{qzv}/taxa-bar-plots.qzv'"),
                    p.job_id, "Taxonomy bar plot", nxt(), TOTAL,
                    must_succeed=False
                )
        else:
            if q:
                await q.put({"type": "log",
                             "msg": "[WARN] No classifier found — skipping taxonomy.\n"
                                    f"       Download classifier to: {CLASSIFIER_DIR}"})

        # ── Step 6: Phylogeny (16S / 12S / PacBio only) ──────────────
        rooted_tree_qza = None
        if use_phylo:
            rooted_tree_qza = str(qza / "rooted-tree.qza")
            await run_step(
                q2(f"qiime phylogeny align-to-tree-mafft-fasttree "
                   f"--i-sequences '{repseqs_qza}' "
                   f"--o-alignment '{qza}/aligned-rep-seqs.qza' "
                   f"--o-masked-alignment '{qza}/masked-aligned-rep-seqs.qza' "
                   f"--o-tree '{qza}/unrooted-tree.qza' "
                   f"--o-rooted-tree '{rooted_tree_qza}' "
                   f"--p-n-threads {p.n_threads_phylo}"),
                p.job_id, "Phylogenetic tree (MAFFT + FastTree)", nxt(), TOTAL
            )

        # ── Step 7: Diversity ─────────────────────────────────────────
        div_dir = str(exp / "diversity")
        if p.metadata_path and Path(p.metadata_path).exists():
            if use_phylo and rooted_tree_qza and Path(rooted_tree_qza).exists():
                div_cmd = (
                    f"qiime diversity core-metrics-phylogenetic "
                    f"--i-phylogeny '{rooted_tree_qza}' "
                    f"--i-table '{table_qza}' "
                    f"--p-sampling-depth {p.sampling_depth} "
                    f"--m-metadata-file '{p.metadata_path}' "
                    f"--output-dir '{qzv}/core-metrics/'"
                )
            else:
                div_cmd = (
                    f"qiime diversity core-metrics "
                    f"--i-table '{table_qza}' "
                    f"--p-sampling-depth {p.sampling_depth} "
                    f"--m-metadata-file '{p.metadata_path}' "
                    f"--output-dir '{qzv}/core-metrics/'"
                )
            await run_step(
                q2(div_cmd),
                p.job_id, "Diversity analysis (alpha + beta)", nxt(), TOTAL
            )

            # Alpha significance
            core_dir = f"{qzv}/core-metrics"
            for metric in ["shannon_vector", "observed_features_vector",
                           "evenness_vector"] + (["faith_pd_vector"] if use_phylo else []):
                metric_qza = f"{core_dir}/{metric}.qza"
                if Path(metric_qza).exists():
                    await run_step(
                        q2(f"qiime diversity alpha-group-significance "
                           f"--i-alpha-diversity '{metric_qza}' "
                           f"--m-metadata-file '{p.metadata_path}' "
                           f"--o-visualization '{qzv}/{metric}-significance.qzv'"),
                        p.job_id, f"Alpha significance — {metric}",
                        nxt(), TOTAL, must_succeed=False
                    )

            # Beta significance (PERMANOVA)
            for dist in ["bray_curtis_distance_matrix", "jaccard_distance_matrix"] + \
                        (["weighted_unifrac_distance_matrix"] if use_phylo else []):
                dist_qza = f"{core_dir}/{dist}.qza"
                if Path(dist_qza).exists():
                    await run_step(
                        q2(f"qiime diversity beta-group-significance "
                           f"--i-distance-matrix '{dist_qza}' "
                           f"--m-metadata-file '{p.metadata_path}' "
                           f"--m-metadata-column '{p.group_col}' "
                           f"--p-method permanova "
                           f"--o-visualization '{qzv}/{dist}-permanova.qzv'"),
                        p.job_id, f"Beta significance PERMANOVA — {dist}",
                        nxt(), TOTAL, must_succeed=False
                    )

        # ── Step 8: ANCOM-BC differential abundance ───────────────────
        if p.run_diffabund and p.metadata_path and Path(p.metadata_path).exists():
            formula = p.diffabund_formula or p.group_col
            collapsed_qza = str(qza / "collapsed-table-l6.qza")
            differentials_qza = str(qza / "ancombc-differentials.qza")

            if Path(taxonomy_qza).exists():
                await run_step(
                    q2(f"qiime taxa collapse "
                       f"--i-table '{table_qza}' "
                       f"--i-taxonomy '{taxonomy_qza}' "
                       f"--p-level 6 "
                       f"--o-collapsed-table '{collapsed_qza}'"),
                    p.job_id, "Collapse taxonomy (genus level)", nxt(), TOTAL,
                    must_succeed=False
                )

                if Path(collapsed_qza).exists():
                    await run_step(
                        q2(f"qiime composition ancombc "
                           f"--i-table '{collapsed_qza}' "
                           f"--m-metadata-file '{p.metadata_path}' "
                           f"--p-formula '{formula}' "
                           f"--o-differentials '{differentials_qza}'"),
                        p.job_id, "ANCOM-BC differential abundance", nxt(), TOTAL,
                        must_succeed=False
                    )

                    if Path(differentials_qza).exists():
                        await run_step(
                            q2(f"qiime composition da-barplot "
                               f"--i-data '{differentials_qza}' "
                               f"--p-significance-threshold 0.05 "
                               f"--o-visualization '{qzv}/ancombc-barplot.qzv'"),
                            p.job_id, "ANCOM-BC visualization", nxt(), TOTAL,
                            must_succeed=False
                        )

        # ── Step 9: Export ────────────────────────────────────────────
        exports = [
            ("feature-table", str(qza / "table.qza"),       str(exp / "feature-table")),
            ("rep-seqs",      str(qza / "rep-seqs.qza"),     str(exp / "rep-seqs")),
            ("taxonomy",      str(qza / "taxonomy.qza"),     str(exp / "taxonomy")),
        ]
        if rooted_tree_qza and Path(rooted_tree_qza).exists():
            exports.append(("rooted-tree", rooted_tree_qza, str(exp / "tree")))

        for name, src, dst in exports:
            if Path(src).exists():
                Path(dst).mkdir(parents=True, exist_ok=True)
                await run_step(
                    q2(f"qiime tools export --input-path '{src}' --output-path '{dst}/'"),
                    p.job_id, f"Export {name}", nxt(), TOTAL,
                    must_succeed=False
                )

        # Convert BIOM → TSV
        biom_path = exp / "feature-table" / "feature-table.biom"
        if biom_path.exists():
            tsv_path = exp / "feature-table" / "feature-table.tsv"
            await run_step(
                q2(f"biom convert -i '{biom_path}' -o '{tsv_path}' --to-tsv"),
                p.job_id, "Export BIOM → TSV", nxt(), TOTAL,
                must_succeed=False
            )

        # ── Step 10: R Visualization ──────────────────────────────────
        if p.run_r_viz:
            r_script = R_SCRIPTS_DIR / "viz_pipeline.R"
            if r_script.exists():
                meta_arg = f"--metadata '{p.metadata_path}'" if p.metadata_path else ""
                group_arg = f"--group_col '{p.group_col}'"
                marker_arg = f"--marker '{p.marker}'"
                await run_step(
                    f"Rscript '{r_script}' "
                    f"--output_dir '{p.output_dir}' "
                    f"{meta_arg} {group_arg} {marker_arg}",
                    p.job_id, "R visualization (phyloseq + ggplot2 + ANCOMBC2)",
                    nxt(), TOTAL, must_succeed=False
                )
            else:
                if q:
                    await q.put({"type": "log",
                                 "msg": "[WARN] viz_pipeline.R not found — skipping R viz"})

        # ── Done ──────────────────────────────────────────────────────
        _jobs[p.job_id]["status"] = "done"
        _jobs[p.job_id]["output_dir"] = p.output_dir
        if q:
            await q.put({"type": "progress", "pct": 100, "label": "Pipeline complete!"})
            await q.put({"type": "pipeline_done", "success": True,
                         "output_dir": p.output_dir})

    except RuntimeError as e:
        _jobs[p.job_id]["status"] = "error"
        if q:
            await q.put({"type": "pipeline_done", "success": False, "error": str(e)})
    except Exception as e:
        _jobs[p.job_id]["status"] = "error"
        if q:
            await q.put({"type": "pipeline_done", "success": False,
                         "error": f"Unexpected error: {e}"})
