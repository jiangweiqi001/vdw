#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

BOOTSTRAP_PYTHON="${BOOTSTRAP_PYTHON:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"

if ! "$BOOTSTRAP_PYTHON" -m venv "$VENV_DIR"; then
  echo "[setup-wsl-pyscf] Could not create venv with $BOOTSTRAP_PYTHON -m venv." >&2
  echo "[setup-wsl-pyscf] On Ubuntu, install venv support with:" >&2
  echo "  sudo apt update && sudo apt install -y python3-venv python3-pip" >&2
  exit 1
fi

"$VENV_DIR/bin/python" -m pip install --upgrade pip
"$VENV_DIR/bin/python" -m pip install -r requirements.txt
"$VENV_DIR/bin/python" -c 'import pyscf; print("[setup-wsl-pyscf] pyscf:", pyscf.__version__)'

echo "[setup-wsl-pyscf] ready: $VENV_DIR/bin/python"
echo "[setup-wsl-pyscf] next: bash scripts/run_semicore_sr_smoke.sh"
