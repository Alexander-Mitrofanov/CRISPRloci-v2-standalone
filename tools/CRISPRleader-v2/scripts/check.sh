#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$project_root"
python_runtime="${PYTHON:-python}"
python_runtime="$($python_runtime -c 'import sys; print(sys.executable)')"
"$python_runtime" -m pytest -q
"$python_runtime" -m ruff check src tests scripts
"$python_runtime" -m mypy src
"$python_runtime" -m build --no-isolation
"$python_runtime" scripts/check_wheel.py
