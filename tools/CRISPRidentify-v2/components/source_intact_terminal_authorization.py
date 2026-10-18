"""Explicit terminal-evidence authorization for source-intact recovery.

The ordinary structural filter can evaluate spacer-terminal dominance only
when at least five non-empty spacers are present.  A five-repeat CRISPR array
has exactly four spacers, so requiring the literal status ``passed`` makes a
declared five-repeat recovery minimum unreachable.

This module supplies a narrow source-only exception.  Four-spacer candidates
remain eligible only when the complete source-intact existence gate passes,
all four spacers are non-empty and distinct, neither spacer edge is unanimous,
and every other structural predicate passes.  It does not perform category
mutation, classification, boundary certification, or the mandatory coding-
tandem veto.
"""

from dataclasses import asdict, dataclass
import math
from numbers import Real
from typing import Optional

from components.components_detection_refinement import (
    AdvancedFuzzySearchFilterEvidence,
    SpacerTerminalDominanceEvidence,
)
from components.source_intact_existence_rescue import (
    SourceIntactExistenceDecision,
)


STATISTICALLY_PASSED = "statistically_passed"
GUARDED_FOUR_SPACER = "guarded_four_spacer_source_diversity"
REJECTED = "rejected"


@dataclass(frozen=True)
class SourceIntactTerminalAuthorizationRule:
    """Frozen policy for the only accepted insufficient-sample branch."""

    name: str = "source-intact-terminal-authorization"
    version: str = "1.0.0"
    small_sample_repeat_count: int = 5
    small_sample_spacer_count: int = 4
    structural_minimum_terminal_sample_size: int = 5
    maximum_small_sample_edge_dominance: float = 0.75

    def __post_init__(self):
        for field_name in ("name", "version"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value:
                raise ValueError(
                    "{} must be a non-empty string".format(field_name)
                )
        for field_name in (
            "small_sample_repeat_count",
            "small_sample_spacer_count",
            "structural_minimum_terminal_sample_size",
        ):
            value = getattr(self, field_name)
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 1
            ):
                raise ValueError(
                    "{} must be a positive integer".format(field_name)
                )
        if self.small_sample_repeat_count != self.small_sample_spacer_count + 1:
            raise ValueError(
                "small-sample repeat/spacer cardinality is inconsistent"
            )
        if (
            self.structural_minimum_terminal_sample_size
            <= self.small_sample_spacer_count
        ):
            raise ValueError(
                "small-sample branch must remain below the statistical "
                "minimum"
            )
        dominance = self.maximum_small_sample_edge_dominance
        if (
            not isinstance(dominance, Real)
            or isinstance(dominance, bool)
            or not math.isfinite(float(dominance))
            or not 0.0 <= float(dominance) < 1.0
        ):
            raise ValueError(
                "maximum_small_sample_edge_dominance must be finite in [0, 1)"
            )

    def as_dict(self):
        return asdict(self)


SOURCE_INTACT_TERMINAL_AUTHORIZATION_RULE = (
    SourceIntactTerminalAuthorizationRule()
)


@dataclass(frozen=True)
class SourceIntactTerminalAuthorization:
    """One immutable authorization decision with explicit observed evidence."""

    authorized: bool
    mode: str
    reason_codes: tuple
    observed_terminal_status: Optional[str]
    nonempty_spacer_count: Optional[int]
    begin_dominance: Optional[float]
    end_dominance: Optional[float]
    source_intact_existence_eligible: bool
    structural_filter_passed: bool
    rule: SourceIntactTerminalAuthorizationRule

    def as_dict(self):
        return asdict(self)


def _finite_fraction(value):
    return (
        isinstance(value, Real)
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and 0.0 <= float(value) <= 1.0
    )


def authorize_source_intact_terminal_evidence(
    existence,
    structural,
    rule=SOURCE_INTACT_TERMINAL_AUTHORIZATION_RULE,
):
    """Authorize ordinary or narrowly guarded small-sample terminal evidence.

    The returned decision is not sufficient to admit a candidate.  Callers
    must still require complete family/context evidence, fresh source-CDS
    coding-tandem clearance, source reconstruction, and transactional conflict
    handling.
    """

    if not isinstance(rule, SourceIntactTerminalAuthorizationRule):
        raise TypeError(
            "rule must be a SourceIntactTerminalAuthorizationRule"
        )
    if not isinstance(existence, SourceIntactExistenceDecision):
        raise TypeError("existence must be a SourceIntactExistenceDecision")
    if not isinstance(structural, AdvancedFuzzySearchFilterEvidence):
        raise TypeError(
            "structural must be an AdvancedFuzzySearchFilterEvidence"
        )

    reasons = []
    if not existence.eligible:
        reasons.append("source_intact_existence_not_eligible")
    if structural.decision != "passed":
        reasons.append("structural_filter_not_passed")
    terminal = structural.terminal_dominance
    if not isinstance(terminal, SpacerTerminalDominanceEvidence):
        reasons.append("terminal_evidence_missing_or_invalid")
        return SourceIntactTerminalAuthorization(
            authorized=False,
            mode=REJECTED,
            reason_codes=tuple(sorted(set(reasons))),
            observed_terminal_status=None,
            nonempty_spacer_count=None,
            begin_dominance=None,
            end_dominance=None,
            source_intact_existence_eligible=existence.eligible,
            structural_filter_passed=structural.decision == "passed",
            rule=rule,
        )

    if terminal.status == "passed":
        authorized = not reasons
        return SourceIntactTerminalAuthorization(
            authorized=authorized,
            mode=STATISTICALLY_PASSED if authorized else REJECTED,
            reason_codes=tuple(sorted(set(reasons))),
            observed_terminal_status=terminal.status,
            nonempty_spacer_count=terminal.nonempty_spacer_count,
            begin_dominance=terminal.begin_dominance,
            end_dominance=terminal.end_dominance,
            source_intact_existence_eligible=existence.eligible,
            structural_filter_passed=structural.decision == "passed",
            rule=rule,
        )

    metrics = existence.metrics
    if terminal.status != "insufficient_sample":
        reasons.append("terminal_status_not_authorized")
    if terminal.nonempty_spacer_count != rule.small_sample_spacer_count:
        reasons.append("small_sample_terminal_spacer_count_differs")
    if (
        terminal.minimum_sample_size
        != rule.structural_minimum_terminal_sample_size
    ):
        reasons.append("terminal_statistical_minimum_differs")
    if metrics.repeat_count != rule.small_sample_repeat_count:
        reasons.append("small_sample_repeat_count_differs")
    if metrics.supplied_spacer_count != rule.small_sample_spacer_count:
        reasons.append("small_sample_supplied_spacer_count_differs")
    if metrics.nonempty_supplied_spacer_count != rule.small_sample_spacer_count:
        reasons.append("small_sample_nonempty_spacer_count_differs")
    if metrics.distinct_spacer_count != rule.small_sample_spacer_count:
        reasons.append("small_sample_spacers_not_all_distinct")
    if metrics.exact_duplicate_spacer_pair_count != 0:
        reasons.append("small_sample_duplicate_spacer_pair_present")
    for side, value in (
        ("begin", terminal.begin_dominance),
        ("end", terminal.end_dominance),
    ):
        if not _finite_fraction(value):
            reasons.append("small_sample_{}_dominance_invalid".format(side))
        elif value > rule.maximum_small_sample_edge_dominance:
            reasons.append(
                "small_sample_{}_dominance_above_maximum".format(side)
            )

    authorized = not reasons
    return SourceIntactTerminalAuthorization(
        authorized=authorized,
        mode=GUARDED_FOUR_SPACER if authorized else REJECTED,
        reason_codes=tuple(sorted(set(reasons))),
        observed_terminal_status=terminal.status,
        nonempty_spacer_count=terminal.nonempty_spacer_count,
        begin_dominance=terminal.begin_dominance,
        end_dominance=terminal.end_dominance,
        source_intact_existence_eligible=existence.eligible,
        structural_filter_passed=structural.decision == "passed",
        rule=rule,
    )


__all__ = (
    "GUARDED_FOUR_SPACER",
    "REJECTED",
    "SOURCE_INTACT_TERMINAL_AUTHORIZATION_RULE",
    "STATISTICALLY_PASSED",
    "SourceIntactTerminalAuthorization",
    "SourceIntactTerminalAuthorizationRule",
    "authorize_source_intact_terminal_evidence",
)
