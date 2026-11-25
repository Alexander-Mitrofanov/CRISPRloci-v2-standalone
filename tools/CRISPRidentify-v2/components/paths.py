"""Checkout-relative assets shared by every package, independent of cwd."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MODEL_ROOT = PROJECT_ROOT / "trained_models"
TOOLS_ROOT = PROJECT_ROOT / "tools"
NATIVE_LIBRARY = PROJECT_ROOT / "components" / "native" / "libliteral_fuzzy.so"


def read_release_version():
    """Read the release label, including in incomplete source distributions."""
    try:
        return (PROJECT_ROOT / "VERSION").read_text(encoding="utf-8").strip()
    except OSError:
        return "unknown"
