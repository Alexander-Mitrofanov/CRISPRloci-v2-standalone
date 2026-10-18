"""RepeatTyper: preserved classifier with explicit uncertainty and provenance."""
__version__ = "2.0.1"

from .inference import classify
from .input import Record

__all__ = ["Record", "__version__", "classify"]
