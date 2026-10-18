#!/usr/bin/env bash
# Shotgun metagenome mock test (beta): simulated Illumina + Nanopore reads of the
# 8 Zymo D6300 bacteria through shotgun_pipeline.py. Needs setup_shotgun.sh and
# setup_typestrains.sh. A few minutes.
#   bash ~/r16s-app/tests/run_shotgun_mock.sh            # both platforms
#   bash ~/r16s-app/tests/run_shotgun_mock.sh illumina --keep ~/sg_mock
APP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$APP/venv/bin/python"; [[ -x "$PY" ]] || PY="python3"
exec "$PY" "$APP/tests/run_shotgun_mock.py" "$@"
