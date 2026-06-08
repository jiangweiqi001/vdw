#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if [[ -z "${PYTHON_BIN:-}" && -x ".venv/bin/python" ]]; then
  PYTHON_BIN=".venv/bin/python"
else
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi
SCAN_OUTPUT="${SCAN_OUTPUT:-results/semicore_c6_q2_candidate_scan_wsl.csv}"
TARGET_OUTPUT="${TARGET_OUTPUT:-results/semicore_c6_workflow_targets_wsl.csv}"

echo "[semicore-sr-smoke] python: $("$PYTHON_BIN" -c 'import sys; print(sys.executable)')"
"$PYTHON_BIN" -c 'import pyscf; print("[semicore-sr-smoke] pyscf:", pyscf.__version__)'

mkdir -p results

echo "[semicore-sr-smoke] scanning Sr q2 large-core candidates"
"$PYTHON_BIN" probe_large_core_q2_candidates.py \
  --atom Sr \
  --pseudo-file external_data/cp2k/GTH_POTENTIALS \
  --pseudo-file external_data/cp2k/POTENTIAL_UZH_CASR_Q2 \
  --candidate-basis-csv external_data/cp2k/large_core_q2_basis_candidates.csv \
  --output "$SCAN_OUTPUT"

echo "[semicore-sr-smoke] summarizing Sr workflow target readiness"
"$PYTHON_BIN" run_semicore_c6_workflow.py \
  --candidate-scan "$SCAN_OUTPUT" \
  --atom Sr \
  --output "$TARGET_OUTPUT"

echo "[semicore-sr-smoke] wrote $SCAN_OUTPUT"
echo "[semicore-sr-smoke] wrote $TARGET_OUTPUT"
