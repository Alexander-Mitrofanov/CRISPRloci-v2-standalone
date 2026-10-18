"""Compatibility imports for mode-specific summary modules."""

from .adapter import _adapter_summary
from .detection import _detection_summary
from .pipeline import _pipeline_summary
from .result_tables import _orientation_summary, _reconstruction_summary

__all__ = [
    "_adapter_summary",
    "_detection_summary",
    "_orientation_summary",
    "_pipeline_summary",
    "_reconstruction_summary",
]
