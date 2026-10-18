#!/usr/bin/env python3
"""Open the local graphical CRISPRloci installer; requires only Python's stdlib."""
import sys

if sys.version_info < (3, 10):
    raise SystemExit("The installer needs Python 3.10 or newer.")

from installer.server import main

if __name__ == "__main__":
    raise SystemExit(main())
