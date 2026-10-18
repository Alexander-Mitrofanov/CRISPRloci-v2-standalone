"""Source-wide Prodigal evidence for the architecture-rescue coding veto.

The architecture policy requires an overlap fraction and CDS score from the
same source-context CDS.  This module computes that pair without changing any
candidate, score, category, or pipeline decision.  Prodigal is run once on the
complete source sequence, regardless of the number of candidates.

All public intervals are zero-based and end-exclusive.  Prodigal's GenBank
locations are one-based and inclusive and are converted exactly once, at the
parser boundary.
"""

from dataclasses import dataclass
import math
import os
import tempfile
from typing import Optional, Sequence, Tuple

from components.architecture_rescue_policy import (
    AUDITED_ARCHITECTURE_RESCUE_POLICY,
    ArchitectureRescuePolicy,
)
from components.external_tools import (
    ExternalToolError,
    require_nonempty_file,
    run_external_tool,
    parse_prodigal_scores,
)
from components.reconstruction_validator import (
    IUPAC_DNA_ALPHABET,
    ReconstructionValidationError,
    validate_candidate_reconstruction,
)


PRODIGAL_MODE = "meta"
_SOURCE_RECORD_IDENTIFIER = "candidate_0"


@dataclass(frozen=True)
class SourceCDS:
    """One validated source-coordinate CDS predicted by Prodigal."""

    start: int
    end: int
    is_complement: bool
    confidence: float

    def __post_init__(self):
        if (
            not isinstance(self.start, int)
            or isinstance(self.start, bool)
            or not isinstance(self.end, int)
            or isinstance(self.end, bool)
            or self.start < 0
            or self.end <= self.start
        ):
            raise ValueError("CDS interval must be zero-based and non-empty")
        if not isinstance(self.is_complement, bool):
            raise TypeError("is_complement must be a boolean")
        if (
            isinstance(self.confidence, bool)
            or not isinstance(self.confidence, (int, float))
            or not math.isfinite(float(self.confidence))
            or not 0.0 <= float(self.confidence) <= 100.0
        ):
            raise ValueError("CDS confidence must be finite and in [0, 100]")
        object.__setattr__(self, "confidence", float(self.confidence))

    @property
    def interval(self):
        return self.start, self.end

    @property
    def strand(self):
        return "-" if self.is_complement else "+"


@dataclass(frozen=True)
class CandidateCodingContextEvidence:
    """The single-CDS evidence pair selected for one source-valid candidate."""

    candidate_index: int
    candidate_start: int
    candidate_end: int
    coding_overlap_fraction: float
    context_cds_score: float
    selected_cds: Optional[SourceCDS]

    @property
    def candidate_interval(self):
        return self.candidate_start, self.candidate_end

    @property
    def selected_cds_interval(self):
        if self.selected_cds is None:
            return None
        return self.selected_cds.interval

    def as_policy_kwargs(self):
        """Return the paired keyword arguments consumed by the rescue policy."""

        return {
            "coding_overlap_fraction": self.coding_overlap_fraction,
            "context_cds_score": self.context_cds_score,
        }


@dataclass(frozen=True)
class BatchSourceCodingContextEvidence:
    """One source-wide Prodigal result and all candidate-specific evidence."""

    source_length: int
    cds: Tuple[SourceCDS, ...]
    candidates: Tuple[CandidateCodingContextEvidence, ...]


class CandidateSourceValidationError(ReconstructionValidationError):
    """A batch member failed independent reconstruction from source DNA."""

    def __init__(self, candidate_index, result):
        self.candidate_index = candidate_index
        super().__init__(result)
        detail = self.args[0]
        self.args = ("candidate {}: {}".format(candidate_index, detail),)


def _normalize_source_sequence(source_sequence):
    if not isinstance(source_sequence, str):
        raise TypeError("source_sequence must be a string")
    if not source_sequence:
        raise ValueError("source_sequence must be non-empty")
    invalid = sorted(set(source_sequence).difference(IUPAC_DNA_ALPHABET))
    if invalid:
        raise ValueError(
            "source_sequence contains non-IUPAC DNA symbols: {}".format(
                ", ".join(repr(symbol) for symbol in invalid)
            )
        )
    return source_sequence.upper()


def _write_source_fasta(source_sequence, output_path):
    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write(">{}\n".format(_SOURCE_RECORD_IDENTIFIER))
        for offset in range(0, len(source_sequence), 80):
            output_file.write(source_sequence[offset:offset + 80] + "\n")


def _predict_normalized_source_cds(source_sequence, runner):
    with tempfile.TemporaryDirectory(
        prefix="crispridentify-source-prodigal-"
    ) as temporary_directory:
        input_path = os.path.join(temporary_directory, "source.fna")
        output_path = os.path.join(temporary_directory, "source.gbk")
        _write_source_fasta(source_sequence, input_path)
        runner(
            "prodigal",
            [
                "-i",
                input_path,
                "-o",
                output_path,
                "-p",
                PRODIGAL_MODE,
            ],
        )
        require_nonempty_file(output_path, "source Prodigal output")
        parsed = parse_prodigal_scores(output_path, [len(source_sequence)])

    cds_records = []
    for (start_one_based, end_inclusive, is_complement), confidence in (
        parsed[0].items()
    ):
        try:
            cds = SourceCDS(
                start=start_one_based - 1,
                end=end_inclusive,
                is_complement=is_complement,
                confidence=confidence,
            )
        except (TypeError, ValueError) as error:
            raise ExternalToolError(
                "Prodigal returned an invalid source CDS: {}".format(error)
            ) from error
        if cds.end > len(source_sequence):
            raise ExternalToolError(
                "source Prodigal CDS lies outside the source sequence: "
                "{}..{} of {}".format(
                    cds.start,
                    cds.end,
                    len(source_sequence),
                )
            )
        cds_records.append(cds)

    return tuple(sorted(
        cds_records,
        key=lambda cds: (
            cds.start,
            cds.end,
            cds.is_complement,
            cds.confidence,
        ),
    ))


def predict_source_cds(source_sequence, runner=None):
    """Run Prodigal once and return deterministic source-coordinate CDS calls."""

    normalized_source = _normalize_source_sequence(source_sequence)
    if runner is None:
        runner = run_external_tool
    if not callable(runner):
        raise TypeError("runner must be callable")
    return _predict_normalized_source_cds(normalized_source, runner)


def _candidate_intervals(source_sequence, candidates):
    if isinstance(candidates, (str, bytes)) or candidates is None:
        raise TypeError("candidates must be an iterable of candidate objects")
    try:
        candidate_records = tuple(candidates)
    except TypeError as error:
        raise TypeError(
            "candidates must be an iterable of candidate objects"
        ) from error

    intervals = []
    for candidate_index, candidate in enumerate(candidate_records):
        try:
            reconstruction = validate_candidate_reconstruction(
                source_sequence,
                candidate,
            )
        except AttributeError as error:
            raise TypeError(
                "candidate {} does not expose standard CRISPRidentify "
                "fields".format(candidate_index)
            ) from error
        if not reconstruction.valid:
            raise CandidateSourceValidationError(
                candidate_index,
                reconstruction,
            )
        if (
            reconstruction.array_start is None
            or reconstruction.array_end is None
        ):
            raise RuntimeError(
                "valid candidate {} has no source interval".format(
                    candidate_index
                )
            )
        intervals.append((
            reconstruction.array_start,
            reconstruction.array_end,
        ))
    return tuple(intervals)


def _overlap_fraction(candidate_interval, cds):
    overlap = max(
        0,
        min(candidate_interval[1], cds.end)
        - max(candidate_interval[0], cds.start),
    )
    return overlap / float(candidate_interval[1] - candidate_interval[0])


def _selected_cds_evidence(candidate_interval, source_cds, policy):
    overlapping = []
    for cds in source_cds:
        overlap_fraction = _overlap_fraction(candidate_interval, cds)
        if overlap_fraction > 0.0:
            triggers_veto = (
                overlap_fraction >= policy.coding_overlap_veto_fraction
                and cds.confidence
                >= policy.context_cds_score_veto_minimum
            )
            overlapping.append((triggers_veto, overlap_fraction, cds))

    if not overlapping:
        return 0.0, 0.0, None

    _, overlap_fraction, selected_cds = max(
        overlapping,
        key=lambda evidence: (
            evidence[0],
            evidence[1],
            evidence[2].confidence,
            -evidence[2].start,
            -evidence[2].end,
            not evidence[2].is_complement,
        ),
    )
    return overlap_fraction, selected_cds.confidence, selected_cds


def batch_source_coding_context_evidence(
    source_sequence,
    candidates: Sequence[object],
    runner=None,
    policy: ArchitectureRescuePolicy = AUDITED_ARCHITECTURE_RESCUE_POLICY,
):
    """Compute paired single-CDS coding evidence for a candidate batch.

    Every candidate is source-validated before Prodigal is invoked.  Selection
    first prefers a CDS whose *own* overlap and confidence meet both veto
    thresholds.  Within the triggering and non-triggering groups, larger
    overlap, then confidence, then stable source coordinates break ties.  This
    prevents independent maxima from different CDS calls from either creating
    or hiding a veto.
    """

    normalized_source = _normalize_source_sequence(source_sequence)
    intervals = _candidate_intervals(normalized_source, candidates)
    if not isinstance(policy, ArchitectureRescuePolicy):
        raise TypeError("policy must be an ArchitectureRescuePolicy")
    if runner is None:
        runner = run_external_tool
    if not callable(runner):
        raise TypeError("runner must be callable")

    source_cds = _predict_normalized_source_cds(normalized_source, runner)
    evidence_records = []
    for candidate_index, candidate_interval in enumerate(intervals):
        overlap_fraction, context_cds_score, selected_cds = (
            _selected_cds_evidence(candidate_interval, source_cds, policy)
        )
        evidence_records.append(CandidateCodingContextEvidence(
            candidate_index=candidate_index,
            candidate_start=candidate_interval[0],
            candidate_end=candidate_interval[1],
            coding_overlap_fraction=overlap_fraction,
            context_cds_score=context_cds_score,
            selected_cds=selected_cds,
        ))

    return BatchSourceCodingContextEvidence(
        source_length=len(normalized_source),
        cds=source_cds,
        candidates=tuple(evidence_records),
    )
