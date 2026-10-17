#!/usr/bin/env bash
# Mock-community tests on the analysis machine (needs R, DADA2 and the databases).
#   bash ~/r16s-app/tests/run_mock_tests.sh                 # every pipeline
#   bash ~/r16s-app/tests/run_mock_tests.sh dada2 cox1      # some
#   bash ~/r16s-app/tests/run_mock_tests.sh --keep ~/mock_test
# See tests/run_mock_tests.py for what is checked. 5-20 minutes.
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$APP/venv/bin/python"; [[ -x "$PY" ]] || PY="python3"
# vsearch / cutadapt / chopper may live in a conda env rather than on PATH
for b in miniconda3 anaconda3 miniforge3 mambaforge; do
  for e in emu ngamp-wgs; do
    d="$HOME/$b/envs/$e/bin"; [[ -d "$d" ]] && PATH="$PATH:$d"
  done
done
export PATH
exec "$PY" "$APP/tests/run_mock_tests.py" "$@"
