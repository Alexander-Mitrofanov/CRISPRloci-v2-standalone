"""Make the source checkout importable without requiring installation."""

from pathlib import Path
import sys


SOURCE = Path(__file__).resolve().parent / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))
