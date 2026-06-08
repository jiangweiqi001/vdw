#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if [[ -z "${PYTHON_BIN:-}" && -x ".venv/bin/python" ]]; then
  PYTHON_BIN=".venv/bin/python"
else
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi

TARGETS="${TARGETS:-results/semicore_c6_workflow_targets_wsl.csv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-results/semicore_sr_validation}"

echo "[semicore-sr-validation] refreshing Sr smoke targets"
PYTHON_BIN="$PYTHON_BIN" TARGET_OUTPUT="$TARGETS" bash scripts/run_semicore_sr_smoke.sh

echo "[semicore-sr-validation] running Sr PSP + core Sternheimer validation"
"$PYTHON_BIN" run_semicore_sr_validation.py \
  --targets "$TARGETS" \
  --output-root "$OUTPUT_ROOT"

echo "[semicore-sr-validation] wrote $OUTPUT_ROOT"
