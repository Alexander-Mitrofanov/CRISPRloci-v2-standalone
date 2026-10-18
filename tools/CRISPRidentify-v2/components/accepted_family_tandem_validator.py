"""Source-only coding-tandem validation for one accepted array family.

This module is deliberately independent of classifier features and pipeline
category mutation.  Every supplied candidate is reconstructed from the
original source before its coordinates are considered.  A tandem decision is
made only in the strand and reading frame of one supplied source CDS; there is
no six-frame maximization.
"""

from dataclasses import asdict, dataclass
import hashlib
from itertools import combinations
from typing import Optional, Tuple

from components.architecture_rescue_coding_context import SourceCDS
from components.reconstruction_validator import (
    ReconstructionValidationError,
    validate_candidate_reconstruction,
)
from components.sequence_utils import reverse_complement


RULE_NAME = "accepted-family-source-cds-tandem"
RULE_VERSION = "1.0.0"
COORDINATE_SYSTEM = "zero_based_half_open"

MINIMUM_FAMILY_ENVELOPE_COVERAGE = 0.80
MINIMUM_CDS_ENVELOPE_COVERAGE = 0.80
MINIMUM_CDS_CONFIDENCE = 80.0
MINIMUM_ADJACENT_PERIOD_COUNT = 5
MINIMUM_PERIOD = 30
MAXIMUM_PERIOD = 240
MINIMUM_ADJACENT_PEPTIDE_IDENTITY = 0.60
MINIMUM_MEAN_ADJACENT_PEPTIDE_IDENTITY = 0.70
MINIMUM_MEAN_ALL_PAIR_PEPTIDE_IDENTITY = 0.65
MINIMUM_MEAN_ADJACENT_DNA_IDENTITY = 0.75


_STANDARD_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class AcceptedFamilyTandemDecision:
    """Frozen, runtime-object-free evidence for one direct tandem decision."""

    supported: bool
    decision: str
    reason: str
    rule_name: str
    rule_version: str
    coordinate_system: str
    accepted_interval: Tuple[int, int]
    family_evidence_interval: Optional[Tuple[int, int]]
    family_envelope_coverage: float
    selected_cds_interval: Optional[Tuple[int, int]]
    selected_cds_strand: Optional[str]
    selected_cds_confidence: Optional[float]
    selected_cds_envelope_coverage: float
    period: Optional[int]
    period_unit_intervals: Tuple[Tuple[int, int], ...]
    translation_frame_offset: Optional[int]
    peptide_lengths: Tuple[int, ...]
    adjacent_peptide_identities: Tuple[float, ...]
    minimum_adjacent_peptide_identity: Optional[float]
    mean_adjacent_peptide_identity: Optional[float]
    mean_all_pair_peptide_identity: Optional[float]
    adjacent_dna_identities: Tuple[float, ...]
    mean_adjacent_dna_identity: Optional[float]
    family_unique_path_count: int
    source_sha256: str
    accepted_source_slice_sha256: str
    family_source_slice_sha256: Optional[str]
    selected_cds_source_slice_sha256: Optional[str]
    period_unit_source_slice_sha256s: Tuple[str, ...]

    def as_dict(self):
        """Return a JSON-serializable report record with the frozen rule."""

        record = asdict(self)
        record["rule"] = {
            "name": RULE_NAME,
            "version": RULE_VERSION,
            "minimum_family_envelope_coverage": (
                MINIMUM_FAMILY_ENVELOPE_COVERAGE
            ),
            "minimum_cds_envelope_coverage": (
                MINIMUM_CDS_ENVELOPE_COVERAGE
            ),
            "minimum_cds_confidence": MINIMUM_CDS_CONFIDENCE,
            "minimum_adjacent_period_count": (
                MINIMUM_ADJACENT_PERIOD_COUNT
            ),
            "minimum_period": MINIMUM_PERIOD,
            "maximum_period": MAXIMUM_PERIOD,
            "period_must_be_divisible_by_three": True,
            "minimum_adjacent_peptide_identity": (
                MINIMUM_ADJACENT_PEPTIDE_IDENTITY
            ),
            "minimum_mean_adjacent_peptide_identity": (
                MINIMUM_MEAN_ADJACENT_PEPTIDE_IDENTITY
            ),
            "minimum_mean_all_pair_peptide_identity": (
                MINIMUM_MEAN_ALL_PAIR_PEPTIDE_IDENTITY
            ),
            "minimum_mean_adjacent_dna_identity": (
                MINIMUM_MEAN_ADJACENT_DNA_IDENTITY
            ),
            "translation_scope": "selected_source_cds_strand_and_frame_only",
            "translation_must_be_stop_and_x_free": True,
            "all_full_period_units_must_be_inside_selected_cds": True,
            "classifier_features_used": False,
            "six_frame_maximization_used": False,
        }
        return record


@dataclass(frozen=True)
class _CandidatePath:
    starts: Tuple[int, ...]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    envelope: Tuple[int, int]


@dataclass(frozen=True)
class _Witness:
    path: _CandidatePath
    family_coverage: float
    cds: SourceCDS
    cds_coverage: float
    period: int
    units: Tuple[Tuple[int, int], ...]
    frame_offset: int
    peptides: Tuple[str, ...]
    adjacent_peptide_identities: Tuple[float, ...]
    minimum_adjacent_peptide_identity: float
    mean_adjacent_peptide_identity: float
    mean_all_pair_peptide_identity: float
    adjacent_dna_identities: Tuple[float, ...]
    mean_adjacent_dna_identity: float
    translation_valid: bool
    metric_reason: Optional[str]

    @property
    def supported(self):
        return self.translation_valid and self.metric_reason is None


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _overlap_fraction(interval, covering_interval):
    overlap = max(
        0,
        min(interval[1], covering_interval[1])
        - max(interval[0], covering_interval[0]),
    )
    return overlap / float(interval[1] - interval[0])


def _candidate_path(source, candidate, label):
    try:
        repeats = tuple(candidate.list_repeats)
        starts = tuple(candidate.list_repeat_starts)
        tuple(candidate.list_spacers)
    except AttributeError as error:
        raise TypeError(
            "{} does not expose standard CrisprCandidate fields".format(label)
        ) from error
    reconstruction = validate_candidate_reconstruction(source, candidate)
    if not reconstruction.valid:
        raise ReconstructionValidationError(reconstruction)
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    return _CandidatePath(
        starts=starts,
        repeat_intervals=repeat_intervals,
        envelope=(repeat_intervals[0][0], repeat_intervals[-1][1]),
    )


def _path_key(path):
    return path.starts, path.repeat_intervals


def _validate_source(source):
    if not isinstance(source, str):
        raise TypeError("source must be a string")
    if not source:
        raise ValueError("source must be non-empty")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source must contain ASCII DNA symbols") from error


def _validate_source_cds(source, source_cds):
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be a tuple of SourceCDS objects")
    validated = []
    for index, cds in enumerate(source_cds):
        if not isinstance(cds, SourceCDS):
            raise TypeError(
                "source_cds[{}] is not a SourceCDS".format(index)
            )
        if cds.end > len(source):
            raise ValueError(
                "source_cds[{}] lies outside the source".format(index)
            )
        validated.append(cds)
    return tuple(sorted(
        set(validated),
        key=lambda cds: (
            cds.start,
            cds.end,
            cds.is_complement,
            cds.confidence,
        ),
    ))


def _equal_period_runs(starts):
    deltas = tuple(
        right - left for left, right in zip(starts, starts[1:])
    )
    runs = []
    index = 0
    while index < len(deltas):
        period = deltas[index]
        end = index + 1
        while end < len(deltas) and deltas[end] == period:
            end += 1
        if (
                end - index >= MINIMUM_ADJACENT_PERIOD_COUNT
                and MINIMUM_PERIOD <= period <= MAXIMUM_PERIOD
                and period % 3 == 0):
            runs.append((
                period,
                tuple(
                    (starts[unit_index], starts[unit_index + 1])
                    for unit_index in range(index, end)
                ),
            ))
        index = end
    return tuple(runs)


def _contained_unit_runs(units, cds):
    contained = tuple(
        cds.start <= start and end <= cds.end for start, end in units
    )
    runs = []
    index = 0
    while index < len(units):
        if not contained[index]:
            index += 1
            continue
        end = index + 1
        while end < len(units) and contained[end]:
            end += 1
        if end - index >= MINIMUM_ADJACENT_PERIOD_COUNT:
            runs.append(units[index:end])
        index = end
    return tuple(runs)


def _translate_unit_in_cds_frame(source, unit, cds):
    start, end = unit
    if cds.is_complement:
        oriented = reverse_complement(source[start:end]).upper()
        frame_offset = (end - cds.end) % 3
    else:
        oriented = source[start:end].upper()
        frame_offset = (cds.start - start) % 3
    peptide = "".join(
        _STANDARD_GENETIC_CODE.get(oriented[offset:offset + 3], "X")
        for offset in range(frame_offset, len(oriented) - 2, 3)
    )
    return peptide, frame_offset


def _sequence_identity(first, second):
    if len(first) != len(second):
        raise ValueError("identity requires equal-length source units")
    if not first:
        return 0.0
    return sum(
        left == right for left, right in zip(first, second)
    ) / float(len(first))


def _mean(values):
    return sum(values) / float(len(values))


def _metric_gate_reason(
        minimum_adjacent_peptide_identity,
        mean_adjacent_peptide_identity,
        mean_all_pair_peptide_identity,
        mean_adjacent_dna_identity):
    """Return the first failed calibrated similarity gate, inclusively."""

    if (
            minimum_adjacent_peptide_identity
            < MINIMUM_ADJACENT_PEPTIDE_IDENTITY):
        return "adjacent_peptide_minimum_below_threshold"
    if (
            mean_adjacent_peptide_identity
            < MINIMUM_MEAN_ADJACENT_PEPTIDE_IDENTITY):
        return "adjacent_peptide_mean_below_threshold"
    if (
            mean_all_pair_peptide_identity
            < MINIMUM_MEAN_ALL_PAIR_PEPTIDE_IDENTITY):
        return "all_pair_peptide_mean_below_threshold"
    if mean_adjacent_dna_identity < MINIMUM_MEAN_ADJACENT_DNA_IDENTITY:
        return "adjacent_dna_mean_below_threshold"
    return None


def _build_witness(source, path, family_coverage, cds, cds_coverage,
                   period, units):
    translated = tuple(
        _translate_unit_in_cds_frame(source, unit, cds) for unit in units
    )
    peptides = tuple(item[0] for item in translated)
    frame_offsets = {item[1] for item in translated}
    if len(frame_offsets) != 1:
        raise RuntimeError("equal in-frame period units changed CDS phase")
    frame_offset, = frame_offsets
    adjacent_peptide = tuple(
        _sequence_identity(first, second)
        for first, second in zip(peptides, peptides[1:])
    )
    all_pair_peptide = tuple(
        _sequence_identity(first, second)
        for first, second in combinations(peptides, 2)
    )
    dna_units = tuple(source[start:end].upper() for start, end in units)
    adjacent_dna = tuple(
        _sequence_identity(first, second)
        for first, second in zip(dna_units, dna_units[1:])
    )
    minimum_adjacent_peptide = min(adjacent_peptide)
    mean_adjacent_peptide = _mean(adjacent_peptide)
    mean_all_pair_peptide = _mean(all_pair_peptide)
    mean_adjacent_dna = _mean(adjacent_dna)
    translation_valid = bool(peptides) and all(
        peptide and "*" not in peptide and "X" not in peptide
        for peptide in peptides
    )
    return _Witness(
        path=path,
        family_coverage=family_coverage,
        cds=cds,
        cds_coverage=cds_coverage,
        period=period,
        units=units,
        frame_offset=frame_offset,
        peptides=peptides,
        adjacent_peptide_identities=adjacent_peptide,
        minimum_adjacent_peptide_identity=minimum_adjacent_peptide,
        mean_adjacent_peptide_identity=mean_adjacent_peptide,
        mean_all_pair_peptide_identity=mean_all_pair_peptide,
        adjacent_dna_identities=adjacent_dna,
        mean_adjacent_dna_identity=mean_adjacent_dna,
        translation_valid=translation_valid,
        metric_reason=_metric_gate_reason(
            minimum_adjacent_peptide,
            mean_adjacent_peptide,
            mean_all_pair_peptide,
            mean_adjacent_dna,
        ),
    )


def _witness_sort_key(witness):
    return (
        not witness.supported,
        not witness.translation_valid,
        witness.metric_reason is not None,
        -len(witness.units),
        -witness.minimum_adjacent_peptide_identity,
        -witness.mean_adjacent_peptide_identity,
        -witness.mean_all_pair_peptide_identity,
        -witness.mean_adjacent_dna_identity,
        witness.path.envelope,
        witness.cds.start,
        witness.cds.end,
        witness.cds.is_complement,
        -witness.cds.confidence,
        witness.period,
        witness.units,
    )


def _decision(source, accepted, unique_path_count, reason, witness=None):
    accepted_slice = source[accepted.envelope[0]:accepted.envelope[1]]
    if witness is None:
        return AcceptedFamilyTandemDecision(
            supported=False,
            decision="not_supported",
            reason=reason,
            rule_name=RULE_NAME,
            rule_version=RULE_VERSION,
            coordinate_system=COORDINATE_SYSTEM,
            accepted_interval=accepted.envelope,
            family_evidence_interval=None,
            family_envelope_coverage=0.0,
            selected_cds_interval=None,
            selected_cds_strand=None,
            selected_cds_confidence=None,
            selected_cds_envelope_coverage=0.0,
            period=None,
            period_unit_intervals=(),
            translation_frame_offset=None,
            peptide_lengths=(),
            adjacent_peptide_identities=(),
            minimum_adjacent_peptide_identity=None,
            mean_adjacent_peptide_identity=None,
            mean_all_pair_peptide_identity=None,
            adjacent_dna_identities=(),
            mean_adjacent_dna_identity=None,
            family_unique_path_count=unique_path_count,
            source_sha256=_sha256(source),
            accepted_source_slice_sha256=_sha256(accepted_slice),
            family_source_slice_sha256=None,
            selected_cds_source_slice_sha256=None,
            period_unit_source_slice_sha256s=(),
        )

    family_slice = source[
        witness.path.envelope[0]:witness.path.envelope[1]
    ]
    cds_slice = source[witness.cds.start:witness.cds.end]
    supported = witness.supported
    return AcceptedFamilyTandemDecision(
        supported=supported,
        decision=(
            "tandem_coding_repeat" if supported else "not_supported"
        ),
        reason=(
            "tandem_coding_repeat_supported"
            if supported
            else reason
        ),
        rule_name=RULE_NAME,
        rule_version=RULE_VERSION,
        coordinate_system=COORDINATE_SYSTEM,
        accepted_interval=accepted.envelope,
        family_evidence_interval=witness.path.envelope,
        family_envelope_coverage=witness.family_coverage,
        selected_cds_interval=witness.cds.interval,
        selected_cds_strand=witness.cds.strand,
        selected_cds_confidence=witness.cds.confidence,
        selected_cds_envelope_coverage=witness.cds_coverage,
        period=witness.period,
        period_unit_intervals=witness.units,
        translation_frame_offset=witness.frame_offset,
        peptide_lengths=tuple(len(peptide) for peptide in witness.peptides),
        adjacent_peptide_identities=(
            witness.adjacent_peptide_identities
        ),
        minimum_adjacent_peptide_identity=(
            witness.minimum_adjacent_peptide_identity
        ),
        mean_adjacent_peptide_identity=(
            witness.mean_adjacent_peptide_identity
        ),
        mean_all_pair_peptide_identity=(
            witness.mean_all_pair_peptide_identity
        ),
        adjacent_dna_identities=witness.adjacent_dna_identities,
        mean_adjacent_dna_identity=witness.mean_adjacent_dna_identity,
        family_unique_path_count=unique_path_count,
        source_sha256=_sha256(source),
        accepted_source_slice_sha256=_sha256(accepted_slice),
        family_source_slice_sha256=_sha256(family_slice),
        selected_cds_source_slice_sha256=_sha256(cds_slice),
        period_unit_source_slice_sha256s=tuple(
            _sha256(source[start:end]) for start, end in witness.units
        ),
    )


def validate_accepted_family_tandem(
        source, accepted_candidate, family_candidates, source_cds):
    """Validate direct coding-tandem evidence for one accepted array family.

    ``family_candidates`` may contain the accepted object, omit it, repeat
    objects, or arrive in any order.  All supplied objects are independently
    reconstructed before geometrically duplicate paths are collapsed.
    ``source_cds`` must be the tuple of source-coordinate :class:`SourceCDS`
    records already produced for this source.
    """

    _validate_source(source)
    accepted = _candidate_path(source, accepted_candidate, "accepted_candidate")
    if isinstance(family_candidates, (str, bytes)) or family_candidates is None:
        raise TypeError("family_candidates must be an iterable of candidates")
    try:
        family_candidates = tuple(family_candidates)
    except TypeError as error:
        raise TypeError(
            "family_candidates must be an iterable of candidates"
        ) from error
    reconstructed = []
    for index, candidate in enumerate(family_candidates):
        reconstructed.append(_candidate_path(
            source,
            candidate,
            "family_candidates[{}]".format(index),
        ))
    cds_records = _validate_source_cds(source, source_cds)

    unique_by_key = {_path_key(accepted): accepted}
    for path in reconstructed:
        unique_by_key.setdefault(_path_key(path), path)
    unique_paths = tuple(
        unique_by_key[key] for key in sorted(unique_by_key)
    )

    covered_paths = tuple(
        (path, _overlap_fraction(accepted.envelope, path.envelope))
        for path in unique_paths
        if _overlap_fraction(accepted.envelope, path.envelope)
        >= MINIMUM_FAMILY_ENVELOPE_COVERAGE
    )
    if not covered_paths:
        return _decision(
            source,
            accepted,
            len(unique_paths),
            "family_envelope_coverage_below_minimum",
        )

    eligible_cds = tuple(
        (cds, _overlap_fraction(accepted.envelope, cds.interval))
        for cds in cds_records
        if (
            cds.confidence >= MINIMUM_CDS_CONFIDENCE
            and _overlap_fraction(accepted.envelope, cds.interval)
            >= MINIMUM_CDS_ENVELOPE_COVERAGE
        )
    )
    if not eligible_cds:
        return _decision(
            source,
            accepted,
            len(unique_paths),
            "eligible_source_cds_not_found",
        )

    periodic_paths = tuple(
        (path, coverage, _equal_period_runs(path.starts))
        for path, coverage in covered_paths
        if _equal_period_runs(path.starts)
    )
    if not periodic_paths:
        return _decision(
            source,
            accepted,
            len(unique_paths),
            "equal_in_frame_period_run_not_found",
        )

    witnesses = []
    for path, family_coverage, period_runs in periodic_paths:
        for cds, cds_coverage in eligible_cds:
            for period, units in period_runs:
                for contained_units in _contained_unit_runs(units, cds):
                    witnesses.append(_build_witness(
                        source,
                        path,
                        family_coverage,
                        cds,
                        cds_coverage,
                        period,
                        contained_units,
                    ))
    if not witnesses:
        return _decision(
            source,
            accepted,
            len(unique_paths),
            "full_period_units_not_contained_in_one_cds",
        )

    selected = min(witnesses, key=_witness_sort_key)
    if not selected.translation_valid:
        reason = "selected_cds_translation_contains_stop_or_unknown"
    elif selected.metric_reason is not None:
        reason = selected.metric_reason
    else:
        reason = "tandem_coding_repeat_supported"
    return _decision(
        source,
        accepted,
        len(unique_paths),
        reason,
        selected,
    )


__all__ = (
    "AcceptedFamilyTandemDecision",
    "validate_accepted_family_tandem",
)
