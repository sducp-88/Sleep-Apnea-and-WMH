#!/usr/bin/env bash
set -euo pipefail

# Run the exported WMH analysis steps in the same order as WMH_Analysis.ipynb.
# Usage examples:
#   bash scripts/run_notebook_steps.sh             # run every step
#   START_STEP=2 END_STEP=5 bash scripts/run_notebook_steps.sh
#   PYTHON=python3 bash scripts/run_notebook_steps.sh

PYTHON="${PYTHON:-python}"
START_STEP="${START_STEP:-1}"
END_STEP="${END_STEP:-99}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STEP_DIR="${SCRIPT_DIR}/notebook_steps"

# The notebook cells use repository-root-relative paths, so make that explicit.
cd "${REPO_ROOT}"

for step_file in "${STEP_DIR}"/[0-9][0-9]_*.py; do
  step_num="$(basename "${step_file}" | cut -d_ -f1)"
  # Force base-10 interpretation for zero-padded numbers.
  step_num=$((10#${step_num}))
  if (( step_num < START_STEP || step_num > END_STEP )); then
    continue
  fi
  echo "==> Running step ${step_num}: ${step_file}"
  "${PYTHON}" "${step_file}"
done
