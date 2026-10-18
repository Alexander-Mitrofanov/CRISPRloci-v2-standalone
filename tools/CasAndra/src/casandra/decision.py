"""Deployment-safe Cas detection decision rule."""

from __future__ import annotations


def is_cas_hit(scan: dict[str, object], threshold: float) -> bool:
    """Require real positive-profile evidence as well as the calibrated margin."""
    return scan.get("best_positive_profile") is not None and float(scan["margin"]) >= threshold
