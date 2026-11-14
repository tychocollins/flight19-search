#!/usr/bin/env bash
# Reproduce the final Flight 19 probability map from the repo root. Takes a few minutes.
set -euo pipefail
mkdir -p work
MAGVAR=1 TAGX=_m python3 code/model.py era5 1 --island --record --baxter --radio --fuelearly            | tee work/evidence_east.json
MAGVAR=1 TAGX=_m python3 code/model.py era5 1 --island --record --baxter --radio --fuelearly --zigzag   | tee work/evidence_zigzag.json
python3 code/combine.py
MODE=v3 python3 code/search1945.py mid
python3 code/export_pkg.py
python3 code/plot.py
python3 code/plot_simple.py
