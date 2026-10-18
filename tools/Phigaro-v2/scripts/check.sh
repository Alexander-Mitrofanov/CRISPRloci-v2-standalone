#!/bin/sh
set -eu
python -m pytest -q tests
python -m pip check
python scripts/check_private_data.py workspace
