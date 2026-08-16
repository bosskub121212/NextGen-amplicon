from fastapi import FastAPI, UploadFile, File, BackgroundTasks, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel
from concurrent.futures import ThreadPoolExecutor
import subprocess, os, sys, uuid, json, threading, time, shutil

# Ensure backend/ directory is on sys.path so local modules (license, updater, etc.) can be imported
sys.path.insert(0, os.path.dirname(__file__))

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False
from pathlib import Path
import hashlib as _hashlib

# ── Integrity check: hash license.py at startup ───────────────────────────────
_BACKEND_DIR = Path(__file__).parent
_LICENSE_FILE = _BACKEND_DIR / "license.py"
_LICENSE_HASH_AT_STARTUP: str = ""

def _compute_file_hash(path: Path) -> str:
    try:
        return _hashlib.sha256(path.read_bytes()).hexdigest()
    except Exception:
        return ""

if _LICENSE_FILE.exists():
    _LICENSE_HASH_AT_STARTUP = _compute_file_hash(_LICENSE_FILE)

def _integrity_ok() -> bool:
    """Return False if license.py was modified after the server started."""
    if not _LICENSE_HASH_AT_STARTUP:
        return True  # Can't check → allow (first boot)
    current = _compute_file_hash(_LICENSE_FILE)
    return current == _LICENSE_HASH_AT_STARTUP

app = FastAPI(title="16S/12S Analysis API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:5173",
        "http://localhost:4173",
        "http://127.0.0.1:3000",
        "http://127.0.0.1:5173",
    ],
    allow_origin_regex=r"http://localhost:\d+",   # allow any localhost port
    allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

BASE_DIR      = Path(__file__).parent
UPLOAD_DIR    = BASE_DIR / "uploads"
RESULTS_DIR   = BASE_DIR / "results"
R_SCRIPTS_DIR = BASE_DIR / "r_scripts"
JOBS_FILE     = BASE_DIR / "jobs_history.json"   # ← persistent storage
UPLOAD_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# ── Job store + thread pool ───────────────────────────────────────────────────
jobs: dict = {}
jobs_lock = threading.Lock()
# log_queues for async SSE streaming (used by qiime2_pipeline)
log_queues: dict = {}
MAX_WORKERS = int(os.getenv("MAX_WORKERS", "2"))
executor    = ThreadPoolExecutor(max_workers=MAX_WORKERS)

# ── psutil Process cache (keeps objects alive between polls so cpu_percent works) ──
_proc_cache: dict[int, "psutil.Process"] = {}  # pid → Process

# Warm up system-wide cpu_percent baseline so the first real reading is accurate
if HAS_PSUTIL:
    try:
        import psutil as _ps_init
        _ps_init.cpu_percent(interval=None)
    except Exception:
        pass

# ── Persist helpers ───────────────────────────────────────────────────────────
def load_jobs():
    """Load job history from disk on startup."""
    global jobs
    if JOBS_FILE.exists():
        try:
            with open(JOBS_FILE, "r") as f:
                loaded = json.load(f)
            # Mark any jobs that were running/queued when we crashed as error
            for j in loaded.values():
                if j.get("status") in ("running", "queued"):
                    j["status"]     = "error"
                    j["error"]      = "Server restarted while job was running"
                    j["step_label"] = "Interrupted (server restart)"
            jobs = loaded
            print(f"[history] Loaded {len(jobs)} jobs from disk.")
        except Exception as e:
            print(f"[history] Could not load jobs: {e}")

def save_jobs():
    """Persist current jobs dict to disk."""
    try:
        snapshot = {}
        for jid, j in jobs.items():
            snapshot[jid] = {k: v for k, v in j.items() if k != "log_lines"}
        with open(JOBS_FILE, "w") as f:
            json.dump(snapshot, f, indent=2)
    except Exception as e:
        print(f"[history] Save failed: {e}")

# Load on startup
load_jobs()

# ── Params schema ─────────────────────────────────────────────────────────────
class MetaRow(BaseModel):
    sampleId:    str = ""
    group:       str = ""
    description: str = ""
    model_config = {"extra": "allow"}   # accept custom columns

class RunParams(BaseModel):
    job_name:      str   = ""
    # --- Marker / pipeline type ---
    # 16S, 12S   → dada2_pipeline.R
    # ITS1, ITS2 → its_pipeline.R
    # COX1       → cox1_pipeline.R
    # 18S-nema   → dada2_pipeline.R with nema_db
    # PacBio     → pacbio_pipeline.R
    # ONT-16S    → emu_pipeline.py  (Oxford Nanopore, Emu)
    marker:        str   = "16S"
    # --- 16S / 12S / nema params ---
    truncLen_F:    int   = 240
    truncLen_R:    int   = 200
    maxEE_F:       float = 2.0
    maxEE_R:       float = 2.0
    trimLeft_F:    int   = 0
    trimLeft_R:    int   = 0
    nbases:        float = 1e8
    pool:          str   = "FALSE"
    chimeraMethod: str   = "consensus"
    taxDatabase:   str   = "SILVA_138.2"
    dbPath:        str   = ""
    minBoot:       int   = 50
    topN:          int   = 30
    metadata:      list[MetaRow] = []
    # --- ITS params ---
    its_region:    str   = "ITS1"   # ITS1 or ITS2
    primer_f:           str   = ""    # forward primer (blank = skip cutadapt)
    primer_r:           str   = ""    # reverse primer
    discardUntrimmed:   bool  = False # discard reads where primer not found
    cutadaptErrorRate:  float = 0.1  # -e  fraction mismatches allowed
    cutadaptOverlap:    int   = 3    # -O  minimum overlap bases
    # --- COX1 params ---
    truncLen_cox1_f: int   = 230
    truncLen_cox1_r: int   = 200
    codon_table:     int   = 5      # 5=invertebrate mt, 2=vertebrate mt
    cox1_min_len:    int   = 300
    cox1_max_len:    int   = 330
    run_lulu:        bool  = True
    # --- PacBio params ---
    pb_min_len:    int   = 1000
    pb_max_len:    int   = 1600
    pb_maxEE:      float = 3.0
    pb_region:     str   = "V1-V9"
    # --- ONT params ---
    ont_region:        str   = "V1-V9"   # V1-V9 or V7-V8 etc.
    ont_min_abundance: float = 0.0001    # Emu min abundance threshold
    ont_db_path:       str   = ""        # path to Emu database directory
    # --- Functional prediction ---
    run_tax4fun:   bool  = False
    run_picrust2:  bool  = False
    # --- Sequencer type (for 16S/12S pipelines) ---
    sequencerType: str   = "illumina"   # illumina | ont | qiime2_vsearch
    # --- VSEARCH OTU clustering identity, used when sequencerType == "qiime2_vsearch" ---
    otuSimilarity: float = 0.97
    # --- ONT long-read filtering for DADA2 (used when sequencerType == "ont") ---
    # DADA2 was built for short, high-accuracy Illumina reads with a hard truncLen.
    # ONT reads vary in length and have far higher (mostly indel) error rates, so
    # instead of truncLen we filter by length range; BAND_SIZE/HOMOPOLYMER_GAP_PENALTY
    # in the R script's dada() call compensate for indel-heavy errors at denoising time.
    ontMinLen: int = 300   # minimum read length to keep (bp)
    ontMaxLen: int = 600   # maximum read length to keep (bp)
    # --- Manual sample↔file pairing (overrides filename auto-detection) ---
    # Each entry: {"sample": "AC", "file1": "AC_1.fastq.gz", "file2": "AC_2.fastq.gz"}
    # file2 blank/omitted → that sample is treated as single-end/long-read.
    # If any entry has a non-empty file2, ALL entries must (fully paired job).
    sampleFileMap: list[dict] = []
    # --- Custom classifier (QIIME2) ---
    customClassifierMode: str   = "default"  # default | train | upload
    customClassifierPath: str   = ""          # path to .qza (upload mode)
    trainAmpliconMinLen:  int   = 200         # min len for extract-reads (train mode)
    trainAmpliconMaxLen:  int   = 600         # max len for extract-reads (train mode)
    # --- CPU threads (applies to all pipelines) ---
    nThreads:      int   = 4
    # --- Shared db paths override ---
    db_paths_json: str   = ""

class WorkerConfig(BaseModel):
    max_workers: int = 2

# ── License endpoints ─────────────────────────────────────────────────────────
try:
    from license import check_license, activate_license, deactivate_license, get_machine_id as _get_mid

    @app.get("/license/status")
    def license_status():
        return check_license()

    class LicenseKeyBody(BaseModel):
        license_key: str

    @app.post("/license/activate")
    def license_activate(body: LicenseKeyBody):
        return activate_license(body.license_key)

    @app.post("/license/deactivate")
    def license_deactivate():
        return deactivate_license()

    @app.get("/license/machine-id")
    def license_machine_id():
        return {"machine_id": _get_mid()}

except ImportError:
    pass

# ─────────────────────────────────────────────────────────────────────────────
@app.get("/")
def root():
    # Serve frontend if built, otherwise return API info
    for d in [BASE_DIR.parent / "frontend" / "dist", BASE_DIR.parent / "frontend" / "build"]:
        idx = d / "index.html"
        if idx.exists():
            return FileResponse(str(idx))
    return {"message": "16S/12S Analysis API 🧬", "max_workers": MAX_WORKERS}

# ── 1. Upload ─────────────────────────────────────────────────────────────────
@app.post("/upload")
async def upload_files(files: list[UploadFile] = File(...)):
    import zipfile as _zf, io as _io
    job_id  = str(uuid.uuid4())[:8]
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    saved = []
    FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    for file in files:
        content = await file.read()
        fname = file.filename or ""
        if fname.lower().endswith(".zip"):
            # Extract FASTQ files from ZIP
            try:
                with _zf.ZipFile(_io.BytesIO(content)) as zf:
                    for member in zf.namelist():
                        bname = Path(member).name
                        if not bname:
                            continue
                        if any(bname.lower().endswith(ext) for ext in FASTQ_EXTS):
                            out_path = job_dir / bname
                            with zf.open(member) as src, open(out_path, "wb") as dst:
                                dst.write(src.read())
                            saved.append(bname)
            except Exception as e:
                logging.warning(f"ZIP extract failed for {fname}: {e}")
        else:
            path = job_dir / fname
            with open(path, "wb") as f:
                f.write(content)
            saved.append(fname)
    with jobs_lock:
        jobs[job_id] = {
            "status": "uploaded", "files": saved,
            "progress": 0, "log_lines": [], "step_label": "Waiting...",
            "marker": "—", "database": "—",
        }
        save_jobs()
    return {"job_id": job_id, "files": saved}

# ── 2. Run (submit to thread pool) ────────────────────────────────────────────
@app.post("/run/{job_id}")
async def run_analysis(job_id: str, params: RunParams):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    # ── Integrity gate (detect tampered license.py) ───────────────────────────
    if not _integrity_ok():
        return JSONResponse(
            status_code=403,
            content={"error": "integrity_violation",
                     "message": "Application integrity check failed. Please reinstall."}
        )

    # ── License gate ──────────────────────────────────────────────────────────
    try:
        from license import LICENSE_ENABLED, check_license, is_pipeline_allowed
        if LICENSE_ENABLED:
            lic = check_license()
            if lic["status"] not in ("active", "dev", "disabled"):
                return JSONResponse(
                    status_code=403,
                    content={"error": "license_required", "details": lic}
                )
            if lic["status"] == "active" and not is_pipeline_allowed(params.marker):
                return JSONResponse(
                    status_code=403,
                    content={
                        "error": "pipeline_not_licensed",
                        "message": (
                            f"Pipeline '{params.marker}' is not included in your license. "
                            f"Allowed: {', '.join(lic.get('pipelines', []))}"
                        ),
                    }
                )
    except ImportError:
        pass  # license module not available → allow all

    # Snapshot the exact settings this run was submitted with, so it can be
    # reviewed later from Job Detail — regardless of outcome (success/error).
    try:
        params_snapshot = params.model_dump()   # pydantic v2
    except AttributeError:
        params_snapshot = params.dict()         # pydantic v1 fallback

    running = sum(1 for j in jobs.values() if j["status"] == "running")
    with jobs_lock:
        jobs[job_id]["status"]     = "queued" if running >= MAX_WORKERS else "running"
        jobs[job_id]["step_label"] = "Queued — waiting for a free slot..." if running >= MAX_WORKERS else "Starting..."
        jobs[job_id]["marker"]     = params.marker
        jobs[job_id]["database"]   = params.taxDatabase
        jobs[job_id]["job_name"]   = params.job_name
        jobs[job_id]["params"]     = params_snapshot
        jobs[job_id]["run_at"]     = time.time()
        save_jobs()

    # Also persist to disk in the results folder — survives jobs.json pruning
    try:
        run_params_file = RESULTS_DIR / job_id
        run_params_file.mkdir(parents=True, exist_ok=True)
        (run_params_file / "run_params.json").write_text(
            json.dumps(params_snapshot, indent=2, ensure_ascii=False))
    except Exception as _pe:
        print(f"[params] WARNING: failed to write run_params.json: {_pe}")

    executor.submit(run_r_pipeline, job_id, params)
    return {"job_id": job_id, "status": jobs[job_id]["status"]}

def run_r_pipeline(job_id: str, params: RunParams):
    with jobs_lock:
        jobs[job_id]["status"]      = "running"
        jobs[job_id]["step_label"]  = "Initializing pipeline..."
        jobs[job_id]["started_at"]  = time.time()
        jobs[job_id]["pid"]         = None
        save_jobs()

    input_dir  = str(UPLOAD_DIR / job_id)
    output_dir = str(RESULTS_DIR / job_id)
    os.makedirs(output_dir, exist_ok=True)
    log_file = Path(output_dir) / "pipeline.log"

    # ── Manual sample↔file manifest (from frontend pairing UI) ────────────────
    # Written into input_dir so dada2_pipeline.R picks it up and skips
    # filename-based auto-detection entirely (fixes ONT-vs-paired-end mixups).
    if params.sampleFileMap:
        try:
            manifest_path = Path(input_dir) / "sample_manifest.json"
            manifest_path.write_text(json.dumps(params.sampleFileMap, ensure_ascii=False, indent=2))
            print(f"[manifest] Wrote sample_manifest.json ({len(params.sampleFileMap)} samples) for job {job_id}")
        except Exception as _me:
            print(f"[manifest] WARNING: failed to write sample_manifest.json: {_me}")

    # Persist output_dir so replot endpoint can find it
    with jobs_lock:
        jobs[job_id]["output_dir"] = output_dir

    # Write metadata.csv if provided
    metadata_path = ""
    if params.metadata:
        import csv
        meta_file = Path(output_dir) / "metadata.csv"
        # Collect all column names (including custom ones)
        all_keys = ["sampleId", "group", "description"]
        for row in params.metadata:
            for k in row.model_fields_set | set(row.model_extra or {}):
                if k not in all_keys:
                    all_keys.append(k)
        with open(meta_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_keys, extrasaction="ignore")
            writer.writeheader()
            for row in params.metadata:
                d = {"sampleId": row.sampleId, "group": row.group, "description": row.description}
                if row.model_extra:
                    d.update(row.model_extra)
                writer.writerow(d)
        metadata_path = str(meta_file)
        with jobs_lock:
            jobs[job_id]["metadata_path"] = metadata_path

    # Auto-detect database path from databases/ folder if not specified
    db_path = params.dbPath

    # ── For DADA2 pipelines: reject EMU paths and resolve SILVA from db_paths.json ──
    # EMU databases (species_taxid.fasta / emu_* dirs) cannot be used by DADA2's
    # assignTaxonomy() which expects SILVA-format .fa.gz files.
    _is_dada2_marker = params.marker.upper() not in ("ONT-16S", "ONT16S", "ONT", "PACBIO",
                                                      "ITS1", "ITS2", "ITS", "COX1")
    def _is_bad_dada2_db(p: str) -> bool:
        """Return True if path is unsuitable as primary DADA2 assignTaxonomy database."""
        if not p:
            return False
        n = Path(p).name.lower()
        # EMU-style files/dirs
        if "species_taxid" in p:
            return True
        if any(seg.startswith("emu_") or seg == "emu" for seg in Path(p).parts):
            return True
        # Not a FASTA archive
        if not (p.endswith(".fa.gz") or p.endswith(".fasta.gz")):
            return True
        # assignSpecies/toSpecies files are for addSpecies() only, NOT for assignTaxonomy()
        # (must match the same check dada2_pipeline.R does: grepl("assignSpecies|toSpecies", ...))
        if "assignspecies" in n or "assign_species" in n or "tospecies" in n or "to_species" in n:
            return True
        return False

    if _is_dada2_marker and _is_bad_dada2_db(db_path):
        print(f"[db] dbPath not suitable for DADA2 assignTaxonomy ({Path(db_path).name}) — resolving SILVA from db_paths.json")
        db_path = ""  # reset to trigger SILVA lookup below

    # ── Resolve SILVA from db_paths.json when db_path is unset ─────────────
    if _is_dada2_marker and not db_path:
        _db_paths_file = params.db_paths_json or str(
            Path.home() / "r16s-app" / "backend" / "databases" / "db_paths.json"
        )
        try:
            _dp = json.loads(Path(_db_paths_file).read_text()) if Path(_db_paths_file).exists() else {}
            _marker_up = params.marker.upper()
            _key_prefs = {
                "16S":      ["SILVA_16S_sp", "SILVA_16S"],
                "12S":      ["PR2_18S", "SILVA_16S"],
                "18S-NEMA": ["NemaBase_18S", "PR2_18S"],
            }.get(_marker_up, ["SILVA_16S_sp", "SILVA_16S"])
            for _k in _key_prefs:
                _v = _dp.get(_k, "")
                # Skip if this db_paths entry is itself an assignSpecies-only file
                if _v and Path(_v).exists() and not _is_bad_dada2_db(_v):
                    db_path = _v
                    print(f"[db] Resolved from db_paths.json [{_k}]: {Path(_v).name}")
                    break
            if not db_path:
                print(f"[db] WARNING: No suitable SILVA train_set found in db_paths.json for marker {_marker_up}")
        except Exception as _e:
            print(f"[db] db_paths.json lookup failed: {_e}")

    # ── Legacy fallback: scan root databases/ folder ─────────────────────
    if not db_path or not os.path.exists(db_path):
        db_dir = BASE_DIR / "databases"
        if db_dir.exists():
            db_files = list(db_dir.glob("*.fa.gz")) + list(db_dir.glob("*.fasta.gz"))
            # Priority 1: toGenus trainset (best for full taxonomy plots)
            genus_files = [f for f in db_files if "togenus" in f.name.lower()
                           and any(k in f.name.lower() for k in ("trainset", "train_set"))]
            # Priority 2: any other trainset (toSpecies etc.)
            other_train = [f for f in db_files if any(k in f.name.lower()
                           for k in ("trainset", "train_set", "nr99"))
                           and "togenus" not in f.name.lower()]
            # Priority 3: anything that isn't assignSpecies
            non_assign  = [f for f in db_files if "assignspecies" not in f.name.lower()]

            if genus_files:
                db_path = str(sorted(genus_files)[0])
                print(f"[db] Auto-detected (toGenus): {sorted(genus_files)[0].name}")
            elif other_train:
                db_path = str(sorted(other_train)[0])
                print(f"[db] Auto-detected (trainset): {sorted(other_train)[0].name}")
            elif non_assign:
                db_path = str(sorted(non_assign)[0])
                print(f"[db] Auto-detected: {sorted(non_assign)[0].name}")
            elif db_files:
                db_path = str(sorted(db_files)[0])
                print(f"[db] Fallback db: {sorted(db_files)[0].name}")

    # ── Route to the correct pipeline script ────────────────────────────────
    marker = params.marker.upper()

    db_paths_json = params.db_paths_json or str(
        Path.home() / "r16s-app" / "backend" / "databases" / "db_paths.json"
    )

    if marker in ("ITS1", "ITS2", "ITS"):
        # ── ITS fungal pipeline ─────────────────────────────────────────────
        region = "ITS1" if marker in ("ITS1", "ITS") else "ITS2"
        cmd = [
            "Rscript", str(R_SCRIPTS_DIR / "its_pipeline.R"),
            "--input_dir",  input_dir,
            "--output_dir", output_dir,
            "--its_region", region,
            "--maxEE_f",    str(params.maxEE_F),
            "--maxEE_r",    str(params.maxEE_R),
            "--threads",    str(params.nThreads),
            "--job_name",   params.job_name or job_id,
            "--topN",       str(params.topN),
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        cmd += ["--error_rate", str(params.cutadaptErrorRate),
                "--min_overlap", str(params.cutadaptOverlap)]
        if params.discardUntrimmed:
            cmd += ["--discard_untrimmed", "TRUE"]
        if params.metadata:
            cmd += ["--metadata", metadata_path]
        if db_paths_json:
            cmd += ["--db_paths", db_paths_json]

    elif marker == "COX1":
        # ── COX1 metabarcoding pipeline ────────────────────────────────────
        cmd = [
            "Rscript", str(R_SCRIPTS_DIR / "cox1_pipeline.R"),
            "--input_dir",   input_dir,
            "--output_dir",  output_dir,
            "--truncLen_f",  str(params.truncLen_cox1_f),
            "--truncLen_r",  str(params.truncLen_cox1_r),
            "--maxEE_f",     str(params.maxEE_F),
            "--maxEE_r",     str(params.maxEE_R),
            "--codon_table", str(params.codon_table),
            "--min_length",  str(params.cox1_min_len),
            "--max_length",  str(params.cox1_max_len),
            "--threads",     str(params.nThreads),
            "--job_name",    params.job_name or job_id,
            "--topN",        str(params.topN),
            "--lulu",        str(params.run_lulu).upper(),
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        cmd += ["--error_rate", str(params.cutadaptErrorRate),
                "--min_overlap", str(params.cutadaptOverlap)]
        if params.discardUntrimmed:
            cmd += ["--discard_untrimmed", "TRUE"]
        if params.metadata:
            cmd += ["--metadata", metadata_path]
        if db_paths_json:
            cmd += ["--db_paths", db_paths_json]

    elif marker in ("ONT-16S", "ONT16S", "ONT"):
        # ── ONT 16S pipeline via Emu ────────────────────────────────────────
        emu_script = R_SCRIPTS_DIR.parent / "emu_pipeline.py"
        # Lookup emu db from db_paths.json if ont_db_path not specified
        # Try: ont_db_path → emu_db_mar2026 → emu_silva → any key starting with "emu_"
        _emu_db_auto = ""
        try:
            _dp = json.loads(Path(db_paths_json).read_text()) if Path(db_paths_json).exists() else {}
            for _key in ("emu_db_mar2026", "emu_silva"):
                _v = _dp.get(_key, "")
                if _v and Path(_v).exists():
                    _emu_db_auto = _v
                    break
            if not _emu_db_auto:
                for _key, _v in _dp.items():
                    if _key.startswith("emu_") and _v and Path(_v).exists():
                        _emu_db_auto = _v
                        break
        except Exception:
            pass
        ont_db = params.ont_db_path or _emu_db_auto or db_path
        # Validate: if ont_db points to a FILE (e.g. species_taxid.fasta was selected),
        # use the parent directory instead
        if ont_db and Path(ont_db).is_file():
            print(f"[emu] WARNING: db_path pointed to a file, using parent dir: {ont_db}")
            ont_db = str(Path(ont_db).parent)
        # Also validate taxonomy.tsv exists in the resolved db dir
        if ont_db and not (Path(ont_db) / "taxonomy.tsv").exists():
            print(f"[emu] WARNING: taxonomy.tsv not found in {ont_db}")
        cmd = [
            sys.executable, str(emu_script),
            "--input",         input_dir,
            "--output",        output_dir,
            "--db_path",       ont_db,
            "--threads",       str(params.nThreads),
            "--region",        params.ont_region,
            "--min_abundance", str(params.ont_min_abundance),
            "--topN",          str(params.topN),
            "--marker",        params.marker,
            "--job_name",      params.job_name or job_id,
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        if params.metadata:
            cmd += ["--metadata", metadata_path]

    elif marker == "PACBIO":
        # ── PacBio CCS long-read 16S pipeline ──────────────────────────────
        cmd = [
            "Rscript", str(R_SCRIPTS_DIR / "pacbio_pipeline.R"),
            "--input_dir",  input_dir,
            "--output_dir", output_dir,
            "--min_length", str(params.pb_min_len),
            "--max_length", str(params.pb_max_len),
            "--maxEE",      str(params.pb_maxEE),
            "--threads",    str(params.nThreads),
            "--region",     params.pb_region,
            "--job_name",   params.job_name or job_id,
            "--pool",       params.pool,
            "--dbPath",     db_path,
            "--minBoot",    str(params.minBoot),
            "--topN",       str(params.topN),
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        if params.metadata:
            cmd += ["--metadata", metadata_path]
        if db_paths_json:
            cmd += ["--db_paths_json", db_paths_json]

    elif params.sequencerType == "qiime2_vsearch":
        # ── 16S / 12S / 18S-nema — QIIME2/VSEARCH OTU-picking pipeline ─────
        # Alternative engine to DADA2: mirrors a manually-run QIIME2/VSEARCH
        # OTU-clustering workflow. Calls vsearch/cutadapt/chopper directly
        # (no qiime CLI / conda env required). Reuses the same Emu-format
        # reference database (sequences.fasta + seq2taxid.tsv + taxonomy.tsv)
        # as the ONT-16S (Emu) marker for taxonomy lookup.
        vsearch_script = R_SCRIPTS_DIR.parent / "qiime2_vsearch_pipeline.py"
        _vs_db_auto = ""
        try:
            _dp = json.loads(Path(db_paths_json).read_text()) if Path(db_paths_json).exists() else {}
            for _key in ("emu_db_mar2026", "emu_silva"):
                _v = _dp.get(_key, "")
                if _v and Path(_v).exists():
                    _vs_db_auto = _v
                    break
            if not _vs_db_auto:
                for _key, _v in _dp.items():
                    if _key.startswith("emu_") and _v and Path(_v).exists():
                        _vs_db_auto = _v
                        break
        except Exception:
            pass
        vsearch_db = params.ont_db_path or _vs_db_auto or db_path
        if vsearch_db and Path(vsearch_db).is_file():
            vsearch_db = str(Path(vsearch_db).parent)

        cmd = [
            sys.executable, str(vsearch_script),
            "--input",          input_dir,
            "--output",         output_dir,
            "--db_path",        vsearch_db,
            "--threads",        str(params.nThreads),
            "--otu_similarity", str(params.otuSimilarity),
            "--min_len",        str(params.ontMinLen),
            "--max_len",        str(params.ontMaxLen),
            "--topN",           str(params.topN),
            "--marker",         params.marker,
            "--job_name",       params.job_name or job_id,
            "--cutadapt_error_rate", str(params.cutadaptErrorRate),
            "--cutadapt_overlap",    str(params.cutadaptOverlap),
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        if params.metadata:
            cmd += ["--metadata", metadata_path]

    else:
        # ── 16S / 12S / 18S-nema — standard DADA2 pipeline ─────────────────
        cmd = [
            "Rscript", str(R_SCRIPTS_DIR / "dada2_pipeline.R"),
            "--input",         input_dir,      "--output",        output_dir,
            "--marker",        params.marker,  "--truncLen_F",    str(params.truncLen_F),
            "--truncLen_R",    str(params.truncLen_R),
            "--maxEE_F",       str(params.maxEE_F),
            "--maxEE_R",       str(params.maxEE_R),
            "--trimLeft_F",    str(params.trimLeft_F),
            "--trimLeft_R",    str(params.trimLeft_R),
            "--nbases",        str(int(params.nbases)),
            "--pool",          params.pool,
            "--chimeraMethod", params.chimeraMethod,
            "--taxDatabase",   params.taxDatabase,
            "--dbPath",        db_path,
            "--minBoot",       str(params.minBoot),
            "--topN",          str(params.topN),
            "--threads",       str(params.nThreads),
            "--metadata",      metadata_path,
            "--db_paths_json", db_paths_json,
        ]
        if params.primer_f:
            cmd += ["--primer_f", params.primer_f]
        if params.primer_r:
            cmd += ["--primer_r", params.primer_r]
        cmd += ["--error_rate", str(params.cutadaptErrorRate),
                "--min_overlap", str(params.cutadaptOverlap)]
        if params.discardUntrimmed:
            cmd += ["--discard_untrimmed", "TRUE"]
        if params.run_tax4fun:
            cmd += ["--tax4fun", "TRUE"]
            if db_paths_json:
                cmd += ["--tax4fun_ref", ""]  # will be auto-detected from db_paths.json
        if params.sequencerType == "ont":
            cmd += ["--single_end", "TRUE"]
            cmd += ["--ont_minLen", str(params.ontMinLen),
                    "--ont_maxLen", str(params.ontMaxLen)]

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        with jobs_lock:
            jobs[job_id]["pid"] = proc.pid
        log_lines: list[str] = []
        save_counter = 0

        # ── Step timing tracking ──────────────────────────────────────
        step_history: list[dict] = []
        cur_step_label: str | None = None
        cur_step_start: float | None = None

        with open(log_file, "w") as lf:
            for raw in proc.stdout:  # type: ignore
                line = raw.rstrip()
                if not line:
                    continue
                lf.write(line + "\n"); lf.flush()

                if line.startswith("CHECKPOINT:"):
                    # R script hit a checkpoint — pause and show warning to user
                    try:
                        body = line[11:]  # "low_merge|{json}"
                        chk_type, chk_json = body.split("|", 1) if "|" in body else (body, "{}")
                        chk_data = json.loads(chk_json)
                        with jobs_lock:
                            jobs[job_id]["status"]          = "waiting_checkpoint"
                            jobs[job_id]["step_label"]      = "⚠️  Low merge rate — waiting for user decision"
                            jobs[job_id]["checkpoint_type"] = chk_type.strip()
                            jobs[job_id]["checkpoint_data"] = chk_data
                            log_lines.append(f"⚠️ Checkpoint: merged={chk_data.get('merged_pct',0):.1f}%, nonchim={chk_data.get('nonchim_pct',0):.1f}%")
                            jobs[job_id]["log_lines"] = log_lines[-500:]
                        save_jobs()
                    except Exception as e:
                        print(f"[checkpoint] parse error: {e}")

                elif line.startswith("PROGRESS:"):
                    try:
                        body = line[9:]
                        pct_str, label = body.split("|", 1) if "|" in body else (body, "")
                        new_label = label.strip()

                        # Detect step transition — new label = new step
                        if new_label and new_label != cur_step_label:
                            if cur_step_label is not None:
                                step_history.append({
                                    "label":      cur_step_label,
                                    "started_at": cur_step_start,
                                    "ended_at":   time.time(),
                                    "status":     "completed",
                                })
                            cur_step_label = new_label
                            cur_step_start = time.time()

                        live_steps = step_history + (
                            [{"label": cur_step_label, "started_at": cur_step_start,
                              "ended_at": None, "status": "running"}]
                            if cur_step_label else []
                        )

                        with jobs_lock:
                            # Don't overwrite waiting_checkpoint status
                            if jobs[job_id].get("status") != "waiting_checkpoint":
                                jobs[job_id]["status"] = "running"
                            jobs[job_id]["progress"]     = int(pct_str)
                            jobs[job_id]["step_label"]   = new_label
                            jobs[job_id]["step_history"] = live_steps
                            if new_label:
                                log_lines.append(new_label)
                                jobs[job_id]["log_lines"] = log_lines[-500:]
                    except Exception:
                        pass
                else:
                    # Keep every non-empty output line (previously only lines matching
                    # a keyword shortlist were kept, which silently hid most of the
                    # pipeline's actual output from the UI). Capped at 500 most-recent
                    # lines per job to bound memory — full untruncated output still
                    # always goes to the on-disk log_file above regardless of this cap.
                    log_lines.append(line)
                    with jobs_lock:
                        jobs[job_id]["log_lines"] = log_lines[-500:]

                # Save to disk every ~20 lines
                save_counter += 1
                if save_counter % 20 == 0:
                    save_jobs()

        proc.wait(timeout=7200)
        now = time.time()

        # Close the last step
        if cur_step_label is not None:
            step_history.append({
                "label":      cur_step_label,
                "started_at": cur_step_start,
                "ended_at":   now,
                "status":     "completed" if proc.returncode == 0 else "error",
            })

        # ── PICRUSt2 functional prediction (post-pipeline) ────────────────
        if proc.returncode == 0 and params.run_picrust2:
            try:
                from picrust2_pipeline import run_picrust2
                with jobs_lock:
                    jobs[job_id]["step_label"] = "PICRUSt2 functional prediction..."
                    jobs[job_id]["progress"]   = 95

                asv_fasta = str(Path(output_dir) / "asvs.fasta")
                asv_table = str(Path(output_dir) / "asv_table.csv")
                pic_out   = str(Path(output_dir) / "PICRUSt2")

                if Path(asv_fasta).exists() and Path(asv_table).exists():
                    pic_result = run_picrust2(
                        job_id          = job_id,
                        asv_fasta       = asv_fasta,
                        asv_table_csv   = asv_table,
                        output_dir      = pic_out,
                        threads         = 4,
                    )
                    with jobs_lock:
                        jobs[job_id]["picrust2"] = pic_result
                else:
                    with jobs_lock:
                        jobs[job_id]["picrust2"] = {
                            "success": False,
                            "message": "ASV FASTA or table not found — PICRUSt2 skipped"
                        }
            except Exception as e:
                with jobs_lock:
                    jobs[job_id]["picrust2"] = {"success": False, "message": str(e)}

        # ── Auto-run viz_pipeline.R for DADA2 pipelines ──────────────────────
        # Generates r_tables/*.csv and r_plots/*.pdf used by Edit Charts.
        _viz_script = R_SCRIPTS_DIR / "viz_pipeline.R"
        if proc.returncode == 0 and _is_dada2_marker and _viz_script.exists():
            try:
                with jobs_lock:
                    jobs[job_id]["step_label"] = "Building interactive charts (r_tables)..."
                    jobs[job_id]["progress"]   = 97
                _viz_cmd = [
                    shutil.which("Rscript") or "Rscript",
                    str(_viz_script),
                    "--output_dir", output_dir,
                    "--marker",     params.marker,
                    "--threads",    str(params.nThreads),
                ]
                if metadata_path and Path(metadata_path).exists():
                    _viz_cmd += ["--metadata", metadata_path]
                _viz_proc = subprocess.run(
                    _viz_cmd, capture_output=True, text=True, timeout=600
                )
                with open(log_file, "a") as _lf:
                    _lf.write("\n--- viz_pipeline.R ---\n")
                    if _viz_proc.stdout:
                        _lf.write(_viz_proc.stdout)
                    if _viz_proc.returncode != 0 and _viz_proc.stderr:
                        _lf.write("\n[viz ERROR]\n" + _viz_proc.stderr[-1000:])
                if _viz_proc.returncode == 0:
                    print(f"[viz] r_tables/ populated for job {job_id}")
                else:
                    print(f"[viz] WARNING: viz_pipeline.R exited {_viz_proc.returncode} for job {job_id}")
            except subprocess.TimeoutExpired:
                print(f"[viz] WARNING: viz_pipeline.R timed out for job {job_id}")
            except Exception as _ve:
                print(f"[viz] WARNING: viz_pipeline.R failed: {_ve}")

        with jobs_lock:
            if proc.returncode == 0:
                jobs[job_id]["status"]   = "completed"
                jobs[job_id]["progress"] = 100
            else:
                jobs[job_id]["status"] = "error"
                with open(log_file) as lf:
                    jobs[job_id]["error"] = "".join(lf.readlines()[-30:])
            jobs[job_id]["finished_at"]  = now
            jobs[job_id]["step_history"] = step_history
            save_jobs()

    except subprocess.TimeoutExpired:
        with jobs_lock:
            jobs[job_id]["status"]      = "error"
            jobs[job_id]["error"]       = "Timeout — exceeded 2 hours"
            jobs[job_id]["finished_at"] = time.time()
            save_jobs()
    except Exception as e:
        with jobs_lock:
            jobs[job_id]["status"]      = "error"
            jobs[job_id]["error"]       = str(e)
            jobs[job_id]["finished_at"] = time.time()
            save_jobs()

# ── 3. Status ─────────────────────────────────────────────────────────────────
@app.get("/status/{job_id}")
def get_status(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    return jobs[job_id]

# ── 4. Progress + logs ────────────────────────────────────────────────────────
@app.get("/progress/{job_id}")
def get_progress(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    cpu_pct:     float | None = None
    ram_mb:      float | None = None
    ram_total_mb: float | None = None

    if jobs[job_id].get("status") == "running" and HAS_PSUTIL:
        try:
            pid = jobs[job_id].get("pid")

            # ── CPU: use system-wide measurement ───────────────────────────────
            # Per-process cpu_percent(interval=None) returns 0.0 on first call
            # for every newly-forked child (multithread R forks many short-lived
            # children). System-wide avoids this — and for single-job workloads
            # it accurately reflects what the pipeline is consuming.
            cpu_pct = round(psutil.cpu_percent(interval=None), 1)

            # ── RAM: sum parent + all children (recursive) ──────────────────────
            if pid:
                if pid not in _proc_cache:
                    _proc_cache[pid] = psutil.Process(pid)

                p = _proc_cache[pid]
                try:
                    children = p.children(recursive=True)
                except psutil.NoSuchProcess:
                    children = []

                # Cache new children so future cpu_percent calls have a baseline
                for child in children:
                    if child.pid not in _proc_cache:
                        try:
                            _proc_cache[child.pid] = child
                            child.cpu_percent(interval=None)   # establish baseline
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass

                all_procs = [p] + [_proc_cache.get(c.pid, c) for c in children]
                ram_bytes = 0
                for pp in all_procs:
                    try:
                        if pp.is_running():
                            ram_bytes += pp.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                ram_mb = round(ram_bytes / 1024 / 1024, 1)
            else:
                # No PID yet — use system RAM used as fallback
                vm = psutil.virtual_memory()
                ram_mb = round(vm.used / 1024 / 1024, 1)

            # Always report machine totals for the UI
            vm_total     = psutil.virtual_memory()
            ram_total_mb = round(vm_total.total / 1024 / 1024, 1)

        except Exception:
            # Process may have ended — clean up cache
            if pid and pid in _proc_cache:
                del _proc_cache[pid]

    # Collect static machine info (available even when not running)
    cpu_info: dict = {}
    if HAS_PSUTIL:
        try:
            cpu_info = {
                "logical":  psutil.cpu_count(logical=True)  or 1,
                "physical": psutil.cpu_count(logical=False) or 1,
                "ram_total_mb": round(
                    psutil.virtual_memory().total / 1024 / 1024, 1),
            }
        except Exception:
            pass

    return {
        "percent":      jobs[job_id].get("progress", 0),
        "logs":         jobs[job_id].get("log_lines", []),
        "step_label":   jobs[job_id].get("step_label", ""),
        "cpu_pct":      cpu_pct,
        "ram_mb":       ram_mb,
        "ram_total_mb": ram_total_mb,
        "started_at":   jobs[job_id].get("started_at"),
        **cpu_info,
    }

# ── 5. Job detail (step breakdown) ───────────────────────────────────────────
@app.get("/detail/{job_id}")
def get_detail(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Not found"})
    j = jobs[job_id]
    started  = j.get("started_at")
    finished = j.get("finished_at")
    total_secs = round(finished - started, 1) if started and finished else None
    return {
        "status":       j.get("status"),
        "started_at":   started,
        "finished_at":  finished,
        "total_secs":   total_secs,
        "step_history": j.get("step_history", []),
        "error":        j.get("error", ""),
        "params":       j.get("params", {}),
        "run_at":       j.get("run_at"),
    }

# ── 5b. All jobs dashboard ────────────────────────────────────────────────────
@app.get("/jobs")
def list_jobs():
    summary = []
    for jid, j in jobs.items():
        summary.append({
            "job_id":     jid,
            "job_name":   j.get("job_name", ""),
            "status":     j.get("status"),
            "progress":   j.get("progress", 0),
            "marker":     j.get("marker", "—"),
            "database":   j.get("database", "—"),
            "files":      j.get("files", []),
            "step_label": j.get("step_label", ""),
            "error":      j.get("error", ""),
            "started_at": j.get("started_at"),
        })
    return {"jobs": summary, "max_workers": MAX_WORKERS}

# ── 6. Cancel job ─────────────────────────────────────────────────────────────
@app.delete("/jobs/{job_id}")
def cancel_job(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    with jobs_lock:
        jobs[job_id]["status"]     = "cancelled"
        jobs[job_id]["step_label"] = "Cancelled by user"
        save_jobs()
    return {"job_id": job_id, "status": "cancelled"}

# ── 6b. Checkpoint continue/abort ─────────────────────────────────────────────
@app.post("/jobs/{job_id}/checkpoint/continue")
def checkpoint_continue(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    output_dir = jobs[job_id].get("output_dir", "")
    if not output_dir:
        return JSONResponse(status_code=400, content={"error": "No output dir"})
    # Write signal file — R script is polling for this
    signal_file = Path(output_dir) / "checkpoint_signal"
    signal_file.write_text("continue")
    with jobs_lock:
        jobs[job_id]["status"]     = "running"
        jobs[job_id]["step_label"] = "Continuing pipeline after checkpoint..."
        jobs[job_id].pop("checkpoint_data", None)
        jobs[job_id].pop("checkpoint_type", None)
        save_jobs()
    return {"ok": True, "message": "Pipeline continuing"}

@app.post("/jobs/{job_id}/checkpoint/abort")
def checkpoint_abort(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    output_dir = jobs[job_id].get("output_dir", "")
    if output_dir:
        signal_file = Path(output_dir) / "checkpoint_signal"
        signal_file.write_text("abort")
    with jobs_lock:
        jobs[job_id]["status"]     = "cancelled"
        jobs[job_id]["step_label"] = "Aborted at checkpoint — adjust settings and re-run"
        jobs[job_id].pop("checkpoint_data", None)
        jobs[job_id].pop("checkpoint_type", None)
        save_jobs()
    return {"ok": True, "message": "Pipeline aborted"}

# ── 6c. Reset job to "uploaded" so it can be re-run with new params -----------
@app.post("/jobs/{job_id}/reset")
def reset_job(job_id: str):
    """
    Reset a paused/aborted/error job back to 'uploaded' so the user can
    re-submit it with adjusted parameters without re-uploading files.
    Clears results folder so the pipeline starts fresh.
    """
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    j = jobs[job_id]
    allowed = ("waiting_checkpoint", "cancelled", "error", "completed")
    if j.get("status") not in allowed:
        return JSONResponse(status_code=400,
                            content={"error": f"Cannot reset a job with status '{j.get('status')}'"})

    # Write abort signal if pipeline is still waiting at checkpoint
    output_dir = j.get("output_dir", "")
    if output_dir:
        signal_file = Path(output_dir) / "checkpoint_signal"
        try:
            signal_file.write_text("abort")
        except Exception:
            pass
        # Remove results so R starts fresh
        try:
            shutil.rmtree(output_dir, ignore_errors=True)
        except Exception:
            pass

    with jobs_lock:
        jobs[job_id]["status"]          = "uploaded"
        jobs[job_id]["progress"]        = 0
        jobs[job_id]["step_label"]      = "Waiting to re-run..."
        jobs[job_id]["step_history"]    = []
        jobs[job_id]["log_lines"]       = []
        jobs[job_id].pop("checkpoint_data", None)
        jobs[job_id].pop("checkpoint_type", None)
        jobs[job_id].pop("output_dir", None)
        jobs[job_id].pop("error", None)
        jobs[job_id].pop("pid", None)
        save_jobs()

    return {"ok": True, "job_id": job_id, "status": "uploaded"}

# ── 6d. Get checkpoint data ────────────────────────────────────────────────────
@app.get("/jobs/{job_id}/checkpoint")
def get_checkpoint(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    j = jobs[job_id]
    if j.get("status") != "waiting_checkpoint":
        return JSONResponse(status_code=400, content={"error": "No active checkpoint"})
    return {
        "checkpoint_type": j.get("checkpoint_type"),
        "checkpoint_data": j.get("checkpoint_data", {}),
    }

# ── 7. Replot with custom colours ─────────────────────────────────────────────
class ReplotBody(BaseModel):
    colors: dict   # {"Phylum": {"Bacillota": "#ff0000", ...}, "Class": {...}, ...}

@app.post("/replot/{job_id}")
def replot_job(job_id: str, body: ReplotBody):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    j = jobs[job_id]
    if j.get("status") not in ("completed", "error"):
        return JSONResponse(status_code=400, content={"error": "Job not completed yet"})

    output_dir = j.get("output_dir", "")
    if not output_dir or not os.path.isdir(output_dir):
        return JSONResponse(status_code=400, content={"error": "Output directory not found"})

    # Save colours to a JSON file the R script will read
    colors_file = os.path.join(output_dir, "custom_colors.json")
    with open(colors_file, "w") as f:
        json.dump(body.colors, f)

    # Also save a copy of the metadata file path so replot.R can load it
    meta_path = j.get("metadata_path", "")

    # Build replot.R command
    replot_script = str(R_SCRIPTS_DIR / "replot.R")
    cmd = [
        "Rscript", replot_script,
        "--output",     output_dir,
        "--colorsFile", colors_file,
    ]
    if meta_path and os.path.isfile(meta_path):
        import shutil
        dest = os.path.join(output_dir, "metadata.csv")
        # Only copy if the source and destination are different files
        if os.path.abspath(meta_path) != os.path.abspath(dest):
            shutil.copy(meta_path, dest)

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            return JSONResponse(status_code=500, content={
                "error": "Replot failed",
                "detail": result.stderr[-2000:] if result.stderr else ""
            })
        return {"ok": True, "message": "Plots regenerated with custom colours"}
    except subprocess.TimeoutExpired:
        return JSONResponse(status_code=504, content={"error": "Replot timed out (>5 min)"})
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

# ── helpers: delete files on disk for a job ──────────────────────────────────
def _delete_job_files(job_id: str):
    """Delete upload folder and results folder for a job (best-effort)."""
    for folder in (UPLOAD_DIR / job_id, RESULTS_DIR / job_id):
        if folder.exists():
            try:
                shutil.rmtree(folder)
                print(f"[delete] Removed {folder}")
            except Exception as e:
                print(f"[delete] Could not remove {folder}: {e}")

# ── 8. Delete single job from history ─────────────────────────────────────────
@app.delete("/history/{job_id}")
def delete_history(job_id: str):
    """Remove a finished job from history and delete its files on disk."""
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    if jobs[job_id].get("status") == "running":
        return JSONResponse(status_code=400, content={"error": "Cannot delete a running job"})
    _delete_job_files(job_id)
    with jobs_lock:
        del jobs[job_id]
        save_jobs()
    return {"deleted": job_id}

# -- 9a. Clear all finished history ------------------------------------------
@app.delete("/history")
def clear_history():
    with jobs_lock:
        to_del = [jid for jid, j in jobs.items()
                  if j.get("status") in ("completed", "cancelled", "error")]
        for jid in to_del:
            _delete_job_files(jid)
            del jobs[jid]
        save_jobs()
    return {"cleared": len(to_del)}

# -- 9b. Cleanup orphaned result folders on disk -------------------------------
@app.delete("/results/cleanup")
def cleanup_orphan_results():
    """Remove result/upload folders on disk that have no matching job in history."""
    known_ids = set(jobs.keys())
    removed = []
    errors = []
    for folder_root in (RESULTS_DIR, UPLOAD_DIR):
        if not folder_root.exists():
            continue
        for subfolder in folder_root.iterdir():
            if subfolder.is_dir() and subfolder.name not in known_ids:
                try:
                    shutil.rmtree(subfolder)
                    removed.append(subfolder.name)
                    print(f"[cleanup] Removed orphan: {subfolder}")
                except Exception as e:
                    errors.append(str(e))
                    print(f"[cleanup] Could not remove {subfolder}: {e}")
    return {"removed": len(removed), "paths": removed, "errors": errors}

# -- 9c. Database paths -------------------------------------------------------
@app.get("/databases")
def get_databases():
    """Return available taxonomy databases and their paths."""
    db_dir = BASE_DIR / "databases"
    db_paths_file = db_dir / "db_paths.json"
    result = {"db_dir": str(db_dir), "databases": {}, "available": []}
    if db_paths_file.exists():
        try:
            import json as _json
            db_paths = _json.loads(db_paths_file.read_text())
            result["databases"] = db_paths
            result["available"] = [k for k, v in db_paths.items() if v and Path(v).exists()]
        except Exception as e:
            result["error"] = str(e)
    return result

# -- 9b2. Database file browser -----------------------------------------------
@app.get("/databases/browse")
def browse_databases():
    """List all database files inside the databases/ folder for the browser picker."""
    db_dir = BASE_DIR / "databases"
    db_dir.mkdir(parents=True, exist_ok=True)

    DB_EXTS = {".gz", ".fasta", ".fa", ".fastq"}
    files = []
    for f in sorted(db_dir.rglob("*")):
        if f.is_file() and f.suffix in DB_EXTS:
            rel = f.relative_to(db_dir)
            files.append({
                "name":     f.name,
                "path":     str(f),
                "rel_path": str(rel),
                "size_mb":  round(f.stat().st_size / 1_048_576, 1),
            })

    return {
        "db_dir": str(db_dir),
        "files":  files,
    }

# -- 9c. Worker config --------------------------------------------------------
@app.get("/config")
def get_config():
    return {"max_workers": MAX_WORKERS}

@app.post("/config")
def set_config(cfg: WorkerConfig):
    global MAX_WORKERS, executor
    MAX_WORKERS = max(1, min(cfg.max_workers, 16))
    executor.shutdown(wait=False)
    executor = ThreadPoolExecutor(max_workers=MAX_WORKERS)
    return {"max_workers": MAX_WORKERS}

# -- 10. Download ZIP ---------------------------------------------------------
@app.get("/download/{job_id}")
def download_zip(job_id: str):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    output_dir = jobs[job_id].get("output_dir", str(RESULTS_DIR / job_id))
    if not Path(output_dir).exists():
        return JSONResponse(status_code=404, content={"error": "Results not found"})
    import zipfile, io
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath in Path(output_dir).rglob("*"):
            if fpath.is_file():
                zf.write(fpath, fpath.relative_to(output_dir))
    buf.seek(0)
    from fastapi.responses import StreamingResponse
    marker = jobs[job_id].get("marker", "")
    safe_name = (jobs[job_id].get("job_name") or job_id).replace(" ", "_")
    filename = f"{safe_name}_{marker}_results.zip"
    return StreamingResponse(buf, media_type="application/zip",
                             headers={"Content-Disposition": f'attachment; filename="{filename}"'})

# -- 12. Taxonomy for colour picker -------------------------------------------
@app.get("/taxonomy/{job_id}")
def get_taxonomy(job_id: str, level: str = "Phylum"):
    if job_id not in jobs:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    output_dir = jobs[job_id].get("output_dir", "")
    tax_file = Path(output_dir) / "taxonomy_table.csv" if output_dir else None
    if not tax_file or not tax_file.exists():
        tax_file = RESULTS_DIR / job_id / "taxonomy_table.csv"
    if not tax_file or not tax_file.exists():
        return JSONResponse(status_code=404, content={"error": "taxonomy_table.csv not found"})
    import csv
    taxa = set()
    with open(tax_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = row.get(level, "").strip()
            if val and val not in ("NA", ""):
                # Strip SILVA/UNITE prefixes like p__, c__, etc.
                val = val.lstrip("dpcofgs").lstrip("_").strip()
                if val:
                    taxa.add(val)
    return {"level": level, "taxa": sorted(taxa)}

# -- 13. System stats ---------------------------------------------------------
@app.get("/system-stats")
def system_stats():
    stats: dict = {}
    if HAS_PSUTIL:
        vm = psutil.virtual_memory()
        stats["cpu_pct"]    = psutil.cpu_percent(interval=0.2)
        stats["ram_used_mb"]  = round(vm.used / 1024**2)
        stats["ram_total_mb"] = round(vm.total / 1024**2)
        stats["ram_pct"]      = vm.percent
    stats["active_jobs"] = sum(1 for j in jobs.values() if j.get("status") == "running")
    stats["queued_jobs"] = sum(1 for j in jobs.values() if j.get("status") == "queued")
    return stats

# -- 14. Update endpoints -----------------------------------------------------
try:
    from updater import check_update, apply_update, save_token, get_token

    @app.get("/update/check")
    def update_check():
        return check_update()

    @app.post("/update/apply")
    def update_apply():
        return apply_update()

    @app.post("/update/save-token")
    def update_save_token(body: dict):
        token = body.get("token", "")
        if not token:
            return JSONResponse(status_code=400, content={"error": "token required"})
        save_token(token)
        return {"ok": True}

    @app.get("/update/token-status")
    def update_token_status():
        tok = get_token()
        return {"has_token": bool(tok), "configured": bool(tok)}

except ImportError:
    pass

# -- 15a. QIIME2 pipeline router -----------------------------------------------
try:
    from qiime2_pipeline import router as _q2_router, set_shared_stores as _q2_set_stores
    _q2_set_stores(jobs, log_queues)
    app.include_router(_q2_router)
    print("[qiime2] QIIME2 pipeline router registered at /qiime2/*")
except ImportError as _q2_err:
    print(f"[qiime2] qiime2_pipeline.py not found — QIIME2 endpoints disabled: {_q2_err}")
except Exception as _q2_err:
    print(f"[qiime2] Failed to load qiime2_pipeline: {_q2_err}")

# -- 14b. Result Preview endpoints -------------------------------------------
import re as _re
from fastapi import HTTPException
from fastapi.responses import Response

@app.get("/results/{job_id}/tables")
def preview_tables(job_id: str):
    """List available CSVs and PDFs for a job.
    Priority: r_tables/ and r_plots/ subdirs first (viz_pipeline output).
    Fallback: scan root dir for DADA2-generated taxonomy_*.csv and related files.
    """
    job_dir      = RESULTS_DIR / job_id
    tables_dir   = job_dir / "r_tables"
    plots_dir    = job_dir / "r_plots"
    summary_file = job_dir / "summary.json"

    tables = sorted([f.name for f in tables_dir.glob("*.csv")]) if tables_dir.exists() else []
    plots  = sorted([f.name for f in plots_dir.glob("*.pdf")])  if plots_dir.exists()  else []

    # ── DADA2 fallback: root-level CSVs (taxonomy_*.csv, alpha_diversity.csv, etc.) ──
    # dada2_pipeline.R writes directly to job root; expose them as if they're in r_tables/
    if not tables and job_dir.exists():
        _dada2_csv_patterns = (
            "taxonomy_phylum.csv", "taxonomy_class.csv", "taxonomy_order.csv",
            "taxonomy_family.csv", "taxonomy_genus.csv", "taxonomy_species.csv",
            "alpha_diversity.csv", "bray_curtis_distance_matrix.csv",
            "read_tracking.csv",
            "rarefaction.csv", "shannon_rarefaction.csv", "rank_abundance.csv",
            "specaccum.csv", "faiths_pd.csv", "pca_scores.csv",
            "nmds_bray.csv", "nmds_jaccard.csv", "beta_heatmap.csv",
            "jaccard_heatmap.csv", "pca_scree.csv",
            # extra viz CSVs (generated by dada2_extra_viz.R)
            "asv_lengths.csv",
        )
        tables = [p for p in _dada2_csv_patterns if (job_dir / p).exists()]

    # ── DADA2 fallback: root-level PDFs for r_plots ──────────────────────────
    if not plots and job_dir.exists():
        _dada2_pdf_patterns = (
            # viz_pipeline.R naming (r_plots-style)
            "05_beta_UPGMA.pdf", "12_upgma_jaccard.pdf",
            # DADA2 root naming (dada2_pipeline.R direct output)
            "alpha_diversity.pdf", "observed_asvs.pdf",
            "asv_length_distribution.pdf", "rarefaction_curves.pdf",
            "prevalence_abundance.pdf", "qc_readcount_boxplot.pdf",
            "read_tracking_plot.pdf",
            "beta_pcoa.pdf", "beta_upgma.pdf",
            "beta_heatmap.pdf", "beta_heatmap_jaccard.pdf",
            "beta_nmds_jaccard.pdf",
            "taxonomy_heatmap_phylum.pdf", "taxonomy_heatmap_family.pdf",
            "taxonomy_heatmap_genus.pdf",
        )
        plots = [p for p in _dada2_pdf_patterns if (job_dir / p).exists()]

    return {
        "tables":  tables,
        "plots":   plots,
        "summary": json.loads(summary_file.read_text()) if summary_file.exists() else None,
    }

@app.get("/results/{job_id}/table/{filename}")
def preview_table(job_id: str, filename: str):
    """Serve a single CSV table as text/csv.
    Looks in r_tables/ first, then falls back to job root (DADA2 output).
    """
    if not _re.match(r'^[\w\-]+\.csv$', filename):
        raise HTTPException(status_code=400, detail="Invalid filename")
    job_dir = RESULTS_DIR / job_id
    # Try r_tables/ first, then root
    for search_dir in (job_dir / "r_tables", job_dir):
        path = search_dir / filename
        if path.exists():
            return Response(content=path.read_text(encoding="utf-8"), media_type="text/csv")
    raise HTTPException(status_code=404, detail="Table not found")

@app.get("/results/{job_id}/plot/{filename}")
def preview_plot_pdf(job_id: str, filename: str):
    """Serve a pre-generated PDF from r_plots/ or job root."""
    if not _re.match(r'^[\w\-\.]+\.pdf$', filename):
        raise HTTPException(status_code=400, detail="Invalid filename")
    job_dir = RESULTS_DIR / job_id
    # Try r_plots/ first, then root
    for search_dir in (job_dir / "r_plots", job_dir):
        path = search_dir / filename
        if path.exists():
            return FileResponse(str(path), media_type="application/pdf",
                                headers={"Content-Disposition": f'inline; filename="{filename}"'})
    raise HTTPException(status_code=404, detail="Plot not found")

@app.get("/results/{job_id}/settings")
def get_preview_settings(job_id: str):
    """Return saved preview customization settings — legacy endpoint (checks both locations)."""
    # Try new edit_charts/ location first, then legacy root
    for candidate in (
        RESULTS_DIR / job_id / "edit_charts" / "settings.json",
        RESULTS_DIR / job_id / "preview_settings.json",
    ):
        if candidate.exists():
            try:
                return json.loads(candidate.read_text(encoding="utf-8"))
            except Exception:
                pass
    return {}

@app.post("/results/{job_id}/settings")
async def save_preview_settings(job_id: str, request: Request):
    """Legacy save — redirects to edit_charts/settings.json."""
    body = await request.json()
    ec_dir = RESULTS_DIR / job_id / "edit_charts"
    ec_dir.mkdir(parents=True, exist_ok=True)
    (ec_dir / "settings.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"ok": True}

# ── Edit Charts folder endpoints ──────────────────────────────────────────────

@app.get("/results/{job_id}/edit_charts/settings")
def get_edit_charts_settings(job_id: str):
    """Return saved Edit Charts settings from edit_charts/settings.json."""
    path = RESULTS_DIR / job_id / "edit_charts" / "settings.json"
    if not path.exists():
        # Fall back to legacy location
        legacy = RESULTS_DIR / job_id / "preview_settings.json"
        if legacy.exists():
            try:
                return json.loads(legacy.read_text(encoding="utf-8"))
            except Exception:
                pass
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}

@app.post("/results/{job_id}/edit_charts/settings")
async def save_edit_charts_settings(job_id: str, request: Request):
    """Save Edit Charts customization to edit_charts/settings.json (included in ZIP export)."""
    body = await request.json()
    ec_dir = RESULTS_DIR / job_id / "edit_charts"
    ec_dir.mkdir(parents=True, exist_ok=True)
    (ec_dir / "settings.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"ok": True}

@app.post("/results/{job_id}/edit_charts/png")
async def save_edit_charts_png(job_id: str, request: Request):
    """Save a chart PNG to edit_charts/{tab_id}.png (included in ZIP export)."""
    import base64 as _b64
    body = await request.json()
    tab_id = _re.sub(r'[^\w\-]', '_', body.get("tab_id", "chart"))
    b64    = body.get("b64", "")
    if not b64:
        raise HTTPException(status_code=400, detail="Missing b64 image data")
    ec_dir = RESULTS_DIR / job_id / "edit_charts"
    ec_dir.mkdir(parents=True, exist_ok=True)
    png_bytes = _b64.b64decode(b64)
    (ec_dir / f"{tab_id}.png").write_bytes(png_bytes)
    return {"ok": True, "path": f"edit_charts/{tab_id}.png"}

# ── Metadata endpoints ────────────────────────────────────────────────────────

@app.get("/results/{job_id}/metadata")
def get_metadata(job_id: str):
    """Return metadata.csv for a job as JSON array of {SampleID, Group, ...}."""
    path = RESULTS_DIR / job_id / "metadata.csv"
    if not path.exists():
        return {"rows": [], "columns": []}
    try:
        import csv as _csv
        with open(path, newline="", encoding="utf-8") as f:
            reader = _csv.DictReader(f)
            rows   = list(reader)
            cols   = reader.fieldnames or []
        return {"rows": rows, "columns": list(cols)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/results/{job_id}/metadata")
async def save_metadata(job_id: str, request: Request):
    """Save metadata.csv (CSV text) to the job results directory."""
    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    body = await request.body()
    # Accept both JSON {"csv": "..."} and raw CSV text
    try:
        parsed = json.loads(body)
        csv_text = parsed.get("csv", "")
    except Exception:
        csv_text = body.decode("utf-8", errors="replace")
    if not csv_text.strip():
        raise HTTPException(status_code=400, detail="Empty metadata")
    (job_dir / "metadata.csv").write_text(csv_text, encoding="utf-8")
    return {"ok": True, "path": "metadata.csv"}

@app.post("/results/{job_id}/metadata/json")
async def save_metadata_json(job_id: str, request: Request):
    """Save metadata from JSON array [{SampleID, Group, ...}] → metadata.csv."""
    import csv as _csv, io as _io
    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    body = await request.json()
    rows = body.get("rows", [])
    if not rows:
        raise HTTPException(status_code=400, detail="No rows provided")
    cols = list(rows[0].keys())
    buf  = _io.StringIO()
    w    = _csv.DictWriter(buf, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
    (job_dir / "metadata.csv").write_text(buf.getvalue(), encoding="utf-8")
    return {"ok": True, "rows": len(rows)}

@app.post("/results/{job_id}/rerun-viz")
def rerun_visualization(job_id: str):
    """Re-run viz_pipeline.R on an existing job's data to generate r_tables/ and r_plots/."""
    import subprocess as _sp
    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")

    # Read job summary to get marker
    summary_file = job_dir / "summary.json"
    marker = "16S"
    if summary_file.exists():
        try:
            marker = json.loads(summary_file.read_text()).get("marker", "16S")
        except Exception:
            pass

    rscript = shutil.which("Rscript") or "Rscript"
    viz_script = _BACKEND_DIR / "r_scripts" / "viz_pipeline.R"
    cmd = [rscript, str(viz_script),
           "--output_dir", str(job_dir),
           "--marker", marker,
           "--threads", "4"]

    try:
        result = _sp.run(cmd, capture_output=True, text=True, timeout=300)
    except _sp.TimeoutExpired:
        raise HTTPException(status_code=504, detail="Visualization timed out (>5 min)")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    # ── Also run dada2_extra_viz.R directly for any job with asv_table.csv ──
    # This generates rarefaction.csv, pca_scores.csv, etc. for DADA2 16S jobs
    # even if viz_pipeline.R exited early (no QIIME2 artifacts present).
    extra_viz_script = _BACKEND_DIR / "r_scripts" / "run_extra_viz.R"
    asv_table        = job_dir / "asv_table.csv"
    if extra_viz_script.exists() and asv_table.exists():
        try:
            _sp.run(
                [rscript, str(extra_viz_script), str(job_dir)],
                capture_output=True, text=True, timeout=180
            )
        except Exception:
            pass  # non-fatal — extra viz is best-effort

    # ── Run advanced visualizations (clustering heatmap, Venn, phylo tree) ──
    adv_viz_script = _BACKEND_DIR / "r_scripts" / "run_advanced_viz.R"
    if adv_viz_script.exists():
        try:
            _sp.run(
                [rscript, str(adv_viz_script), str(job_dir)],
                capture_output=True, text=True, timeout=300
            )
        except Exception:
            pass  # non-fatal

    return {"ok": True, "stdout": result.stdout[-3000:], "stderr": result.stderr[-2000:]}

@app.post("/results/{job_id}/preview-charts")
async def save_preview_charts(job_id: str, request: Request):
    """Receive rendered PNG charts from frontend and save into job directory."""
    import base64 as _b64
    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Results not found")
    charts_dir = job_dir / "preview_charts"
    charts_dir.mkdir(exist_ok=True)
    body = await request.json()
    saved = 0
    for chart in body.get("charts", []):
        name = chart.get("name", "")
        b64  = chart.get("b64", "")
        if not name or not b64:
            continue
        # Sanitise filename — allow only safe characters
        safe_name = _re.sub(r'[^\w\-\.]', '_', name)
        (charts_dir / safe_name).write_bytes(_b64.b64decode(b64))
        saved += 1
    return {"saved": saved}


@app.get("/results/{job_id}/download")
def download_results(job_id: str):
    """Zip entire results directory and serve — all files included.
    Pre-compressed formats (gz, pdf, png …) use ZIP_STORED to avoid
    slow recompression; text files use ZIP_DEFLATED.
    """
    import zipfile, io as _io
    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Results not found")

    # Store (no re-compression) for already-compressed formats
    store_exts = {
        ".gz", ".bz2", ".zip", ".qza", ".qzv", ".biom",
        ".fastq", ".fq", ".pdf", ".png", ".jpg", ".jpeg",
        ".bam", ".cram",
    }
    # Skip only true temp/lock files
    skip_exts = {".tmp", ".lock", ".part"}

    buf = _io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for f in sorted(job_dir.rglob("*")):
            if not f.is_file():
                continue
            rel = f.relative_to(job_dir)
            ext = f.suffix.lower()
            if ext in skip_exts:
                continue
            compress = zipfile.ZIP_STORED if ext in store_exts else zipfile.ZIP_DEFLATED
            zf.write(f, rel, compress_type=compress)

    buf.seek(0)
    safe_id = _re.sub(r'[^\w\-]', '_', job_id)
    return Response(
        content=buf.read(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{safe_id}_results.zip"'},
    )

# -- 15. Serve frontend static files (SPA) ------------------------------------
# Looks for dist/ next to the backend/ directory (i.e. ~/r16s-app/frontend/dist)
from fastapi.staticfiles import StaticFiles

_FRONTEND_DIRS = [
    BASE_DIR.parent / "frontend" / "dist",   # Vite build output
    BASE_DIR.parent / "frontend" / "build",  # CRA build output (legacy)
]
_FRONTEND_DIR = next((d for d in _FRONTEND_DIRS if d.is_dir()), None)

if _FRONTEND_DIR:
    # Mount static assets (JS/CSS chunks) — must come before the catch-all
    app.mount("/assets", StaticFiles(directory=_FRONTEND_DIR / "assets"), name="assets")

    @app.get("/{full_path:path}", include_in_schema=False)
    def serve_spa(full_path: str):
        """Return index.html for any unknown path so React Router works."""
        index = _FRONTEND_DIR / "index.html"
        if index.exists():
            return FileResponse(str(index))
        return JSONResponse(status_code=404, content={"error": "Frontend not built"})

# -- 15. License endpoints ----------------------------------------------------
class ActivateLicenseRequest(BaseModel):
    license_key: str

try:
    from license import check_license, activate_license, deactivate_license

    @app.get("/license/status")
    def license_status():
        return check_license()

    @app.post("/license/activate")
    def license_activate(req: ActivateLicenseRequest):
        return activate_license(req.license_key)

    @app.delete("/license/deactivate")
    def license_deactivate():
        return deactivate_license()

    @app.post("/license/deactivate")
    def license_deactivate_post():
        return deactivate_license()

    from license import get_machine_id as _get_mid
    @app.get("/license/machine-id")
    def license_machine_id():
        return {"machine_id": _get_mid()}

except Exception as _lic_err:
    print(f"[license] Could not load license module: {_lic_err}")
    import traceback; traceback.print_exc()
