#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$project_root"
python_runtime="${PYTHON:-$project_root/backend/.venv/bin/python}"
target="${1:-all}"
python3 scripts/check_private_data.py workspace
if [[ "$target" == backend || "$target" == all ]]; then
  "$python_runtime" -m pytest -q --import-mode=importlib -o addopts='' backend/tests Integration/tests tests
  "$python_runtime" -m ruff check backend/src backend/tests Integration/crispr_integration scripts tests
fi
if [[ "$target" == frontend || "$target" == all ]]; then
  (cd frontend && npm run check:source && npm test && VITE_API_BASE_URL="${VITE_API_BASE_URL:-https://analysis.example.org}" npm run build)
fi
if [[ "$target" != backend && "$target" != frontend && "$target" != all ]]; then
  echo 'Usage: scripts/check.sh [backend|frontend|all]' >&2
  exit 2
fi
