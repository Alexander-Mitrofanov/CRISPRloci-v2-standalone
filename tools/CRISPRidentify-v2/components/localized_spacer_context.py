"""Source-bound pending observations for a localized duplicate-spacer block.

This module only recognizes and captures an existing complete geometry. It
does not change structural decisions, generate cropped arrays, or classify.
"""
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass, fields, is_dataclass
import hashlib
import json

from components.candidate_origins import candidate_origin_evidence
from components.reconstruction_validator import validate_candidate_reconstruction


RULE = "localized-identical-spacer-context-v1"
MINIMUM_OWN_SCORE = 0.50
SCORE_TIE_TOLERANCE = 1e-12


def plain(value):
    if is_dataclass(value):
        return {field.name: plain(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, tuple) and hasattr(value, "_fields"):
        return {name: plain(getattr(value, name)) for name in value._fields}
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    return value


def encoded(value):
    return json.dumps(plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def localized_block(candidate, structural_filter):
    """Return one fully localized, internal maximal run, or None.

    Flanks include their adjacent block boundary repeat, but no duplicated
    block spacer. They are diagnostic checks, not independent votes.
    """
    spacers = candidate.list_spacers
    consecutive = structural_filter.max_allowed_consecutive_spacers
    total = structural_filter.max_allowed_same_spacers
    minimum = max(3, structural_filter.min_number_repeats)
    if (type(consecutive) is not int or type(total) is not int
            or min(consecutive, total) < 1 or not spacers
            or any(not isinstance(s, str) or not s for s in spacers)):
        return None
    runs, start = [], 0
    while start < len(spacers):
        end = start + 1
        while end < len(spacers) and spacers[end] == spacers[start]:
            end += 1
        runs.append((start, end, spacers[start]))
        start = end
    violating = [(a, b, s) for a, b, s in runs if b-a > consecutive]
    if len(violating) != 1:
        return None
    a, b, sequence = violating[0]
    counts = Counter(spacers)
    if counts[sequence] != b-a or any(count > total for s, count in counts.items() if s != sequence):
        return None
    count = len(candidate.list_repeats)
    if a == 0 or b == len(spacers) or min(a+1, count-b) < minimum:
        return None
    return {
        "rule": RULE, "spacer_run": [a, b], "block_repeat_indexes": [a, b+1],
        "spacer_sequence": sequence, "multiplicity": b-a,
        "left_repeat_indexes": [0, a+1], "right_repeat_indexes": [b, count],
        "minimum_flank_repeats": minimum,
        "maximum_consecutive_identical_spacers": consecutive,
        "maximum_identical_spacers": total,
        "flank_alignment_recipe": "inherited_columns_remove_all_gap_columns_then_own_consensus",
    }


def candidate_capture(candidate, source_hash):
    """Value record of the current unscored observation, never birth history."""
    return plain({
        "repeat_intervals": [[s, s+len(r)] for s, r in
                             zip(candidate.list_repeat_starts, candidate.list_repeats)],
        "repeats": candidate.list_repeats, "spacers": candidate.list_spacers,
        "gapped_repeats": candidate.list_repeats_gaped, "gapped_repeats_alias": candidate.list_gaped_repeats,
        "consensus": candidate.consensus,
        "eligible": getattr(candidate, "source_path_classification_eligible", True),
        "eligibility_reason": getattr(candidate, "source_path_classification_reason", None),
        "origins": candidate_origin_evidence(candidate, source_hash),
        "source_path_provenance": getattr(candidate, "source_path_provenance", ()),
        "restriction_reasons": getattr(candidate, "source_path_restriction_reasons", ()),
        "boundary_repair_records": getattr(candidate, "boundary_repair_records", ()),
        "score_status": "unscored_at_refinement_capture",
    })


@dataclass(frozen=True)
class PendingSpacerContext:
    locus_key: tuple
    candidate: object
    source_sha256: str
    observation_json: str
    block_json: str

    @property
    def observation(self):
        return json.loads(self.observation_json)

    @property
    def block(self):
        return json.loads(self.block_json)

    def validate_current(self, source_hash):
        if self.source_sha256 != source_hash:
            raise ValueError("pending spacer context source differs")
        if encoded(candidate_capture(self.candidate, source_hash)) != self.observation_json:
            raise ValueError("pending spacer context representation or provenance changed")


def capture_pending_context(source, source_hash, locus_key, candidate, structural_filter):
    if getattr(candidate, "source_path_classification_eligible", True) is not True:
        return None
    block = localized_block(candidate, structural_filter)
    if block is None:
        return None
    if not validate_candidate_reconstruction(source, candidate).valid:
        return None
    # Other failing predicates are intentionally retained as possible unsafe
    # rank peers. Capture is not a structural exception or admission.
    owned = deepcopy(candidate)
    return PendingSpacerContext(tuple(locus_key), owned, source_hash,
                                encoded(candidate_capture(owned, source_hash)), encoded(block))
