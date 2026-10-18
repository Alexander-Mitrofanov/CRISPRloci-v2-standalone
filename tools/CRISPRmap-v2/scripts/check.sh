#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python_runtime="${PYTHON:-python3}"
"$python_runtime" scripts/check_private_data.py workspace
"$python_runtime" -m pytest -q
"$python_runtime" -m ruff check src tests scripts benchmarks
