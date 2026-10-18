#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$project_root"
python3 scripts/install_git_guards.py --allow-source-archive
python_runtime="${PYTHON:-python3.12}"
"$python_runtime" -c 'import sys; assert sys.version_info[:2] == (3, 12), "Use Python 3.12"'
if [[ ! -x backend/.venv/bin/python ]]; then
  "$python_runtime" -m venv backend/.venv
fi
backend/.venv/bin/python -m pip install -e './backend[test]' -e ./Integration build playwright
backend/.venv/bin/python -m pip check
(cd frontend && npm ci --ignore-scripts)
python3 scripts/sync_tools.py
printf '%s\n' 'Development environment ready. Continue with scripts/provision_tools.py and the release workflow in README.md.'
