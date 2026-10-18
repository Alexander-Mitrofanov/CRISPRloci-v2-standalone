#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if ! command -v python3 >/dev/null 2>&1; then
  echo 'Python 3.10 or newer is needed to open the installer. Install python3, then run bash install.sh again.' >&2
  exit 1
fi
exec python3 "$project_root/scripts/install.py" "$@"
