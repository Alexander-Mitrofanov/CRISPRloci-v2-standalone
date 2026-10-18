"""Frozen experimental admission thresholds shared by filtering and reporting."""

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class ShortRescuePolicy:
    min_query_span: int = 12
    min_identity: float = 0.8
    min_anchored_paired_bases: int = 10
    min_anchored_paired_fraction: float = 0.5
    energy_strictly_below: int = 0
    cm_evalue_max: float = 0.01
    cm_duplex_overlap_strictly_below: int = 10
    terminator_evalue_max: float = 0.001

    def public(self) -> dict[str, int | float]:
        return {
            "min_query_span": self.min_query_span,
            "min_identity": self.min_identity,
            "min_anchored_paired_bases": self.min_anchored_paired_bases,
            "min_anchored_paired_fraction": self.min_anchored_paired_fraction,
            "energy_strictly_below": self.energy_strictly_below,
            "cm_evalue_max": self.cm_evalue_max,
            "cm_duplex_overlap_strictly_below": self.cm_duplex_overlap_strictly_below,
            "terminator_evalue_max": self.terminator_evalue_max,
        }


SHORT_RESCUE_POLICY: Final = ShortRescuePolicy()
