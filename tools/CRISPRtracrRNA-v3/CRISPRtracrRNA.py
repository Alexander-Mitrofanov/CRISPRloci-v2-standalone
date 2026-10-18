#!/usr/bin/env python3
"""Compatibility launcher for source checkouts; prefer the installed command."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from crisprtracrrna.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
