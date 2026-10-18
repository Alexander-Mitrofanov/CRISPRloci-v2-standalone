"""Compatibility façade for bounded public-summary construction."""

from .summaries.builder import build_summary
from .summaries.common import _MAX_PUBLIC_SUMMARY_BYTES, sanitize_public_summary

__all__ = [
    "_MAX_PUBLIC_SUMMARY_BYTES",
    "build_summary",
    "sanitize_public_summary",
]
