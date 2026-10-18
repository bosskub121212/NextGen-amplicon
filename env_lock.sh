#!/usr/bin/env bash
# =============================================================================
#  env_lock.sh — pin the exact software this machine analyses with
#
#    bash ~/r16s-app/env_lock.sh export     write environment/lock/ from what is installed now
#    bash ~/r16s-app/env_lock.sh check      compare what is installed now with the lock
#    bash ~/r16s-app/env_lock.sh restore    rebuild a MISSING conda env / venv from the lock
#
#  What is locked:
#    conda envs  emu, picrust2, ngamp-wgs, qiime2-amplicon*   conda list --explicit --md5
#                (exact package files + checksums: the env can be recreated bit for bit,
#                 offline from the conda cache, with: conda create -n X --file lock)
#    venvs       ~/r16s-app/venv, ~/cutadapt-env               pip freeze
#    R           R version + every installed package version   (dada2, phyloseq, maaslin3…)
#
#  check exits 1 when anything drifted, so it can gate a run or a deploy.
#  restore never touches an env that exists: a working env is not replaced
#  behind your back. To rebuild one, rename it first (conda rename) and re-run.
# =============================================================================
set -uo pipefail
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCK="$APP/environment/lock"
CMD="${1:-check}"

G='\033[0;32m'; Y='\033[1;33m'; R='\033[0;31m'; B='\033[1m'; N='\033[0m'
ok(){ echo -e "  ${G}✓${N} $*"; }; warn(){ echo -e "  ${Y}!${N} $*"; }; bad(){ echo -e "  ${R}✗${N} $*"; }

# conda / mamba / micromamba, whichever exists
CONDA=""
for c in conda mamba micromamba; do command -v "$c" >/dev/null 2>&1 && { CONDA="$c"; break; }; done
if [[ -z "$CONDA" ]]; then
  for d in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/anaconda3" "$HOME/mambaforge"; do
    [[ -x "$d/bin/conda" ]] && { CONDA="$d/bin/conda"; break; }
  done
fi

conda_envs(){   # names of the envs this app uses that exist here
  [[ -z "$CONDA" ]] && return
  "$CONDA" env list 2>/dev/null | awk '!/^#/ && NF {print $1}' |
    grep -E '^(emu|picrust2|ngamp-wgs|qiime2-amplicon.*)$' || true
}
VENVS=("$APP/venv" "$HOME/cutadapt-env")
venv_tag(){ basename "$1"; }

r_versions(){   # "pkg version" per line, plus R itself
  command -v Rscript >/dev/null 2>&1 || return 0
  Rscript --vanilla -e 'ip <- installed.packages()[, c("Package","Version")];
    ip <- ip[!duplicated(ip[, 1]), , drop = FALSE]; ip <- ip[order(ip[, 1]), , drop = FALSE];
    cat(sprintf("R %s\n", paste(R.version$major, R.version$minor, sep = ".")));
    cat(sprintf("%s %s\n", ip[, 1], ip[, 2]), sep = "")' 2>/dev/null
}

snapshot(){     # write a full snapshot into $1
  local out="$1"; mkdir -p "$out"
  local e
  for e in $(conda_envs); do
    "$CONDA" list -n "$e" --explicit --md5 > "$out/conda-$e.txt" 2>/dev/null \
      && ok "conda env $e ($(grep -c '^http' "$out/conda-$e.txt") packages)" \
      || warn "conda env $e could not be listed"
  done
  local v
  for v in "${VENVS[@]}"; do
    [[ -x "$v/bin/python" ]] || continue
    "$v/bin/python" -m pip freeze --all 2>/dev/null | sort -f > "$out/pip-$(venv_tag "$v").txt" \
      && ok "venv $(venv_tag "$v") ($(wc -l < "$out/pip-$(venv_tag "$v").txt") packages)"
  done
  r_versions > "$out/R-packages.txt"
  [[ -s "$out/R-packages.txt" ]] && ok "R $(head -1 "$out/R-packages.txt" | cut -d' ' -f2) ($(($(wc -l < "$out/R-packages.txt") - 1)) packages)" \
    || { rm -f "$out/R-packages.txt"; warn "Rscript not found — R not locked"; }
  python3 - "$out" "$APP" <<'PY'
import json, sys, socket, datetime, pathlib
out, app = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
try: ver = json.loads((app / "version.json").read_text()).get("version")
except Exception: ver = None
(out / "lock.json").write_text(json.dumps({
    "host": socket.gethostname(), "app_version": ver,
    "written": datetime.datetime.now().isoformat(timespec="seconds"),
    "files": sorted(p.name for p in out.iterdir() if p.name != "lock.json")}, indent=1))
PY
}

compare(){      # compare lock dir $1 with snapshot dir $2; exit status 1 on drift
  python3 - "$1" "$2" <<'PY'
import sys, re, pathlib
L, C = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
def parse(p):
    d = {}
    if not p.exists(): return None
    for ln in p.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith(("#", "@")): continue
        if ln.startswith("http"):                       # conda explicit: .../name-ver-build.conda#md5
            f = ln.split("/")[-1].split("#")[0]
            f = re.sub(r"\.(conda|tar\.bz2)$", "", f)
            parts = f.rsplit("-", 2)
            if len(parts) == 3: d[parts[0]] = f"{parts[1]}-{parts[2]}"
            else: d[f] = f
        elif "==" in ln:                                # pip freeze
            n, v = ln.split("==", 1); d[n.lower()] = v
        elif " @ " in ln:
            n, v = ln.split(" @ ", 1); d[n.lower()] = v
        elif p.name.startswith("R-"):                  # R: "pkg version"
            parts = ln.split()
            if len(parts) == 2: d[parts[0]] = parts[1]
    return d
drift = 0
names = sorted({p.name for p in L.glob("*.txt")} | {p.name for p in C.glob("*.txt")})
for n in names:
    a, b = parse(L / n), parse(C / n)
    label = n[:-4]
    if a is None:
        print(f"  ! {label}: installed here, not in the lock"); drift += 1; continue
    if b is None:
        print(f"  ✗ {label}: in the lock, MISSING on this machine"); drift += 1; continue
    chg = [(k, a[k], b[k]) for k in sorted(a) if k in b and a[k] != b[k]]
    gone = [k for k in sorted(a) if k not in b]
    new = [k for k in sorted(b) if k not in a]
    if not (chg or gone or new):
        print(f"  ✓ {label}: identical ({len(a)} packages)"); continue
    drift += 1
    print(f"  ✗ {label}: {len(chg)} changed, {len(gone)} missing, {len(new)} added")
    for k, x, y in chg[:25]: print(f"      {k}: {x} → {y}")
    if len(chg) > 25: print(f"      … {len(chg) - 25} more")
    for k in gone[:10]: print(f"      − {k} {a[k]}")
    for k in new[:10]:  print(f"      + {k} {b[k]}")
sys.exit(1 if drift else 0)
PY
}

case "$CMD" in
  export)
    echo -e "${B}Locking the environment of $(hostname) → environment/lock/${N}"
    if [[ -d "$LOCK" ]]; then
      bk="$APP/environment/lock.$(date +%Y%m%d-%H%M%S)"
      mv "$LOCK" "$bk" && warn "previous lock kept as $(basename "$bk")"
    fi
    snapshot "$LOCK"
    echo -e "\n  Done. Check it later with:  ${B}bash $APP/env_lock.sh check${N}"
    ;;
  check)
    [[ -d "$LOCK" ]] || { bad "no lock yet — run: bash $APP/env_lock.sh export"; exit 2; }
    echo -e "${B}Installed now vs. environment/lock/${N} ($(python3 -c "import json;d=json.load(open('$LOCK/lock.json'));print(d['host'],d['written'],'v'+str(d['app_version']))" 2>/dev/null))"
    tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
    snapshot "$tmp" >/dev/null
    if compare "$LOCK" "$tmp"; then
      echo -e "\n  ${G}${B}No drift.${N} Results from this machine are reproducible with the locked versions."
    else
      echo -e "\n  ${Y}${B}Drift found.${N} Re-run export if the change was intended (an update), or restore."
      exit 1
    fi
    ;;
  restore)
    [[ -d "$LOCK" ]] || { bad "no lock to restore from"; exit 2; }
    [[ -n "$CONDA" ]] || warn "conda not found — conda envs cannot be restored"
    existing="$(conda_envs)"
    for f in "$LOCK"/conda-*.txt; do
      [[ -e "$f" && -n "$CONDA" ]] || continue
      e="$(basename "$f" .txt)"; e="${e#conda-}"
      if grep -qx "$e" <<<"$existing"; then
        ok "conda env $e exists — left as it is (check shows any drift)"
      else
        echo "  … creating conda env $e from the lock"
        "$CONDA" create -n "$e" --file "$f" -y >/dev/null && ok "conda env $e restored" || bad "conda env $e failed"
      fi
    done
    for v in "${VENVS[@]}"; do
      f="$LOCK/pip-$(venv_tag "$v").txt"
      [[ -f "$f" ]] || continue
      if [[ -x "$v/bin/python" ]]; then
        ok "venv $(venv_tag "$v") exists — left as it is"
      else
        python3 -m venv "$v" && "$v/bin/python" -m pip install -q -r "$f" \
          && ok "venv $(venv_tag "$v") restored" || bad "venv $(venv_tag "$v") failed"
      fi
    done
    if [[ -f "$LOCK/R-packages.txt" ]]; then
      warn "R packages are not reinstalled automatically (CRAN/Bioconductor serve one version per R"
      warn "release). 'check' lists any that differ; pin one with:"
      echo "      Rscript -e 'remotes::install_version(\"PKG\", version = \"X.Y.Z\")'"
    fi
    ;;
  *)
    echo "usage: bash env_lock.sh export | check | restore"; exit 2 ;;
esac
