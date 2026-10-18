"""Conservative low-spacer-diversity evidence for tandem-like arrays.

The rule in this module deliberately addresses only the small-array regime.
It does not treat one duplicated spacer inside an otherwise diverse array as
tandem evidence.  Every spacer must be non-empty canonical DNA and every pair
must meet the same full-length global edit-identity threshold.
"""

from dataclasses import asdict, dataclass
from fractions import Fraction
import math
from numbers import Real


CATEGORY_NAME = "likely_tandem"
RULE_NAME = "complete-link-low-spacer-diversity"
RULE_VERSION = "1.0.0"
COORDINATE_SYSTEM = "source_forward_strand"


@dataclass(frozen=True)
class LikelyTandemPolicy:
    """Fixed, bounded policy for the requested small-array screen."""

    minimum_spacer_count: int = 2
    maximum_spacer_count: int = 8
    minimum_pairwise_edit_identity: float = 0.95
    canonical_alphabet: str = "ACGT"

    def __post_init__(self):
        if (
                isinstance(self.minimum_spacer_count, bool)
                or not isinstance(self.minimum_spacer_count, int)
                or self.minimum_spacer_count < 2):
            raise ValueError("minimum_spacer_count must be an integer >= 2")
        if (
                isinstance(self.maximum_spacer_count, bool)
                or not isinstance(self.maximum_spacer_count, int)
                or self.maximum_spacer_count < self.minimum_spacer_count):
            raise ValueError(
                "maximum_spacer_count must be >= minimum_spacer_count"
            )
        threshold = self.minimum_pairwise_edit_identity
        if (
                isinstance(threshold, bool)
                or not isinstance(threshold, Real)
                or not math.isfinite(float(threshold))
                or not 0.0 <= float(threshold) <= 1.0):
            raise ValueError(
                "minimum_pairwise_edit_identity must be between 0 and 1"
            )
        if (
                not isinstance(self.canonical_alphabet, str)
                or not self.canonical_alphabet
                or len(set(self.canonical_alphabet))
                != len(self.canonical_alphabet)
                or self.canonical_alphabet != self.canonical_alphabet.upper()):
            raise ValueError(
                "canonical_alphabet must contain unique uppercase symbols"
            )

    def as_dict(self):
        return asdict(self)


LIKELY_TANDEM_POLICY = LikelyTandemPolicy()


@dataclass(frozen=True)
class LikelyTandemEvidence:
    """Complete audit evidence for one spacer collection."""

    decision: str
    category: object
    reason_code: str
    rule_name: str
    rule_version: str
    coordinate_system: str
    policy: LikelyTandemPolicy
    spacer_count: int
    nonempty_spacer_count: int
    canonical_spacer_count: int
    distinct_spacer_count: int
    exact_duplicate_excess_count: int
    shortest_spacer_length: object
    longest_spacer_length: object
    pairwise_comparison_count: int
    minimum_pairwise_edit_identity: object
    maximum_pairwise_edit_identity: object
    maximum_pairwise_edit_distance: object
    minimum_identity_pair_1based: object
    minimum_identity_pair_edit_distance: object
    minimum_identity_pair_maximum_length: object
    complete_link_threshold_satisfied: object

    @property
    def likely_tandem(self):
        return self.decision == CATEGORY_NAME

    def as_dict(self):
        value = asdict(self)
        pair = value["minimum_identity_pair_1based"]
        if pair is not None:
            value["minimum_identity_pair_1based"] = list(pair)
        return value


def global_edit_distance(first, second):
    """Return Levenshtein distance with full-sequence end-gap penalties."""

    if not isinstance(first, str) or not isinstance(second, str):
        raise TypeError("edit-distance inputs must be strings")
    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for first_index, first_base in enumerate(first, 1):
        current = [first_index]
        for second_index, second_base in enumerate(second, 1):
            current.append(min(
                current[-1] + 1,
                previous[second_index] + 1,
                previous[second_index - 1]
                + (first_base != second_base),
            ))
        previous = current
    return previous[-1]


def global_edit_identity(first, second):
    """Return ``1 - Levenshtein/max(lengths)`` for two full sequences."""

    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    return 1.0 - global_edit_distance(first, second) / float(denominator)


def _unassessed_evidence(
        policy,
        spacers,
        normalized,
        reason_code,
        canonical_count):
    lengths = tuple(len(spacer) for spacer in normalized)
    distinct_count = len(set(normalized))
    return LikelyTandemEvidence(
        decision="retained",
        category=None,
        reason_code=reason_code,
        rule_name=RULE_NAME,
        rule_version=RULE_VERSION,
        coordinate_system=COORDINATE_SYSTEM,
        policy=policy,
        spacer_count=len(spacers),
        nonempty_spacer_count=sum(bool(spacer) for spacer in normalized),
        canonical_spacer_count=canonical_count,
        distinct_spacer_count=distinct_count,
        exact_duplicate_excess_count=len(normalized) - distinct_count,
        shortest_spacer_length=min(lengths) if lengths else None,
        longest_spacer_length=max(lengths) if lengths else None,
        pairwise_comparison_count=0,
        minimum_pairwise_edit_identity=None,
        maximum_pairwise_edit_identity=None,
        maximum_pairwise_edit_distance=None,
        minimum_identity_pair_1based=None,
        minimum_identity_pair_edit_distance=None,
        minimum_identity_pair_maximum_length=None,
        complete_link_threshold_satisfied=None,
    )


def _identity_meets_threshold(distance, maximum_length, threshold):
    """Compare edit identity with a decimal policy threshold exactly."""

    threshold_fraction = Fraction(str(float(threshold)))
    identity_numerator = maximum_length - distance
    return (
        identity_numerator * threshold_fraction.denominator
        >= maximum_length * threshold_fraction.numerator
    )


def evaluate_likely_tandem(spacers, policy=LIKELY_TANDEM_POLICY):
    """Classify a spacer sequence as a conservative tandem-like quarantine.

    Similarity is complete-link: every pair, not merely a connected chain of
    near neighbours, must meet the threshold.  Spacer reverse complements are
    intentionally not canonicalized; all comparisons use the array's emitted
    source orientation.
    """

    if not isinstance(policy, LikelyTandemPolicy):
        raise TypeError("policy must be a LikelyTandemPolicy")
    if isinstance(spacers, (str, bytes)) or spacers is None:
        raise TypeError("spacers must be a sequence of strings")
    try:
        spacers = tuple(spacers)
    except TypeError as error:
        raise TypeError("spacers must be a sequence of strings") from error
    if any(not isinstance(spacer, str) for spacer in spacers):
        raise TypeError("every spacer must be a string")

    normalized = tuple(spacer.upper() for spacer in spacers)
    alphabet = set(policy.canonical_alphabet)
    canonical_count = sum(
        bool(spacer) and set(spacer).issubset(alphabet)
        for spacer in normalized
    )
    if len(spacers) < policy.minimum_spacer_count:
        return _unassessed_evidence(
            policy, spacers, normalized,
            "spacer_count_below_policy_minimum", canonical_count,
        )
    if len(spacers) > policy.maximum_spacer_count:
        return _unassessed_evidence(
            policy, spacers, normalized,
            "spacer_count_above_policy_maximum", canonical_count,
        )
    if any(not spacer for spacer in normalized):
        return _unassessed_evidence(
            policy, spacers, normalized,
            "empty_spacer_requires_separate_review", canonical_count,
        )
    if canonical_count != len(normalized):
        return _unassessed_evidence(
            policy, spacers, normalized,
            "noncanonical_spacer_requires_separate_review", canonical_count,
        )

    identities = []
    distances = []
    minimum_pair = None
    for left_index, left in enumerate(normalized[:-1]):
        for right_index in range(left_index + 1, len(normalized)):
            right = normalized[right_index]
            distance = global_edit_distance(left, right)
            maximum_length = max(len(left), len(right))
            identity_numerator = maximum_length - distance
            identity = identity_numerator / float(maximum_length)
            identities.append(identity)
            distances.append(distance)
            if (
                    minimum_pair is None
                    or (
                        identity_numerator * minimum_pair[1]
                        < minimum_pair[0] * maximum_length
                    )
                    or (
                        identity_numerator * minimum_pair[1]
                        == minimum_pair[0] * maximum_length
                        and (left_index, right_index)
                        < minimum_pair[2]
                    )):
                minimum_pair = (
                    identity_numerator,
                    maximum_length,
                    (left_index, right_index),
                    distance,
                )

    minimum_identity = minimum_pair[0] / float(minimum_pair[1])
    threshold_satisfied = _identity_meets_threshold(
        minimum_pair[3],
        minimum_pair[1],
        policy.minimum_pairwise_edit_identity,
    )
    distinct_count = len(set(normalized))
    lengths = tuple(len(spacer) for spacer in normalized)
    return LikelyTandemEvidence(
        decision=CATEGORY_NAME if threshold_satisfied else "retained",
        category=CATEGORY_NAME if threshold_satisfied else None,
        reason_code=(
            "all_spacer_pairs_at_or_above_identity_threshold"
            if threshold_satisfied
            else "at_least_one_spacer_pair_below_identity_threshold"
        ),
        rule_name=RULE_NAME,
        rule_version=RULE_VERSION,
        coordinate_system=COORDINATE_SYSTEM,
        policy=policy,
        spacer_count=len(normalized),
        nonempty_spacer_count=len(normalized),
        canonical_spacer_count=len(normalized),
        distinct_spacer_count=distinct_count,
        exact_duplicate_excess_count=len(normalized) - distinct_count,
        shortest_spacer_length=min(lengths),
        longest_spacer_length=max(lengths),
        pairwise_comparison_count=len(identities),
        minimum_pairwise_edit_identity=minimum_identity,
        maximum_pairwise_edit_identity=max(identities),
        maximum_pairwise_edit_distance=max(distances),
        minimum_identity_pair_1based=(
            minimum_pair[2][0] + 1,
            minimum_pair[2][1] + 1,
        ),
        minimum_identity_pair_edit_distance=minimum_pair[3],
        minimum_identity_pair_maximum_length=minimum_pair[1],
        complete_link_threshold_satisfied=threshold_satisfied,
    )


__all__ = [
    "CATEGORY_NAME",
    "COORDINATE_SYSTEM",
    "LIKELY_TANDEM_POLICY",
    "LikelyTandemEvidence",
    "LikelyTandemPolicy",
    "RULE_NAME",
    "RULE_VERSION",
    "evaluate_likely_tandem",
    "global_edit_distance",
    "global_edit_identity",
]
