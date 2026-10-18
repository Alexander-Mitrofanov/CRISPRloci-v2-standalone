"""Generic bounded PWM occurrence scanning, independent of leader decisions.

Matrices hold positive base probabilities. Log2 likelihood ratios are rounded
per position to the nearest integer ``score_scale`` unit (Python ties-to-even).
Dynamic programming enumerates the *complete* null distribution of that integer
score under the stated independent-base background. Tail probabilities are
exact for this discrete statistic, subject only to floating-point summation;
they are not probabilities of a biological function or a leader annotation.

BH q-values cover every eligible position, strand, model, and candidate passed
to one call, including ambiguous windows conservatively assigned p=1. Reporting
thresholds apply only after adjustment and never shrink the testing universe.
"""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math
from typing import Any, Literal

from crisprleader.errors import InputError, ReferenceError
from crisprleader.models import FlankCandidate, MotifHit
from crisprleader.references import BASES, DNA, validate_motif

MAX_NULL_STATES = 100000
DEFAULT_MAX_TESTS = 1000000
COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASE_INDEX = {base: index for index, base in enumerate(BASES)}
COMP_INDEX = (3, 2, 1, 0)


@dataclass(frozen=True)
class _ScoreTable:
    integer_rows: tuple[tuple[int, ...], ...]
    real_rows: tuple[tuple[float, ...], ...]
    scores: tuple[int, ...]
    tails: tuple[float, ...]

    def pvalue(self, score: int) -> float:
        index = bisect_left(self.scores, score)
        return self.tails[index] if index < len(self.scores) else 0.0


@lru_cache(maxsize=32)
def _score_table(
    probability_rows: tuple[tuple[float, ...], ...],
    background: tuple[float, ...],
    scale: int,
) -> _ScoreTable:
    if len(probability_rows) * math.log10(min(background)) < -300:
        raise ReferenceError(
            "motif/background combination exceeds reliable floating tail precision; "
            "the scanner will not report underflowed zero p-values"
        )
    real_rows = tuple(
        tuple(math.log2(probability) - math.log2(background[index])
              for index, probability in enumerate(row))
        for row in probability_rows
    )
    integer_rows = tuple(tuple(round(value * scale) for value in row) for row in real_rows)
    distribution = {0: 1.0}
    for row in integer_rows:
        next_distribution: dict[int, list[float]] = defaultdict(list)
        for score, mass in distribution.items():
            for index, increment in enumerate(row):
                next_distribution[score + increment].append(mass * background[index])
                if len(next_distribution) > MAX_NULL_STATES:
                    raise ReferenceError(
                        "motif null-score distribution exceeds the state bound; "
                        "use a lower declared score_scale or shorter matrix"
                    )
        distribution = {score: math.fsum(masses) for score, masses in next_distribution.items()}
    ordered = sorted(distribution)
    tails_reversed: list[float] = []
    running = 0.0
    for score in reversed(ordered):
        running = math.fsum((running, distribution[score]))
        tails_reversed.append(min(1.0, max(0.0, running)))
    # Sum-of-products floating error is not evidence for a tail below one at
    # the minimum score; its null event includes the entire sample space.
    tails_reversed[-1] = 1.0
    return _ScoreTable(integer_rows, real_rows, tuple(ordered), tuple(reversed(tails_reversed)))


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    """BH adjusted values including ties; no filtering of the testing family."""
    if any(type(value) not in (int, float) or not math.isfinite(value)
           or value < 0 or value > 1 for value in pvalues):
        raise ValueError("p-values must be finite numbers in [0,1]")
    order = sorted(range(len(pvalues)), key=lambda index: pvalues[index])
    result = [1.0] * len(order)
    previous = 1.0
    for rank_minus_one in range(len(order) - 1, -1, -1):
        index = order[rank_minus_one]
        adjusted = pvalues[index] * len(order) / (rank_minus_one + 1)
        previous = min(previous, adjusted, 1.0)
        result[index] = previous
    return result


def _positions(candidate: FlankCandidate, model: dict[str, Any]) -> range:
    width = len(model["pwm"])
    boundary = len(candidate.sequence)
    total = boundary + len(candidate.observed_repeat)
    domain = model["domain"]
    if domain == "leader":
        first, last = 0, boundary - width
    elif domain == "repeat":
        first, last = boundary, total - width
    else:
        # Junction windows must actually contain bases on both sides.
        first, last = max(0, boundary - width + 1), min(boundary - 1, total - width)
    if "junction_end_range" in model:
        minimum, maximum = model["junction_end_range"]
        first = max(first, boundary + minimum - width)
        last = min(last, boundary + maximum - width)
    return range(first, max(first, last + 1))


def _matrix(model: dict[str, Any], strand: str) -> tuple[tuple[float, ...], ...]:
    matrix = tuple(tuple(float(row[base]) for base in BASES) for row in model["pwm"])
    if strand == "-":
        return tuple(tuple(row[index] for index in COMP_INDEX) for row in reversed(matrix))
    return matrix


def _validate_candidates(candidates: Sequence[FlankCandidate]) -> None:
    seen = set()
    for candidate in candidates:
        if candidate.id in seen:
            raise InputError(f"motif scope contains duplicate candidate id: {candidate.id}")
        seen.add(candidate.id)
        if not candidate.observed_repeat:
            raise InputError("motif scanning requires the actual observed terminal-repeat sequence")
        if not set(candidate.junction_sequence) <= DNA:
            raise InputError("motif candidates must contain uppercase IUPAC DNA")


def scan_motifs(
    candidates: Sequence[FlankCandidate],
    motifs: Sequence[dict[str, Any]],
    *,
    max_tests: int = DEFAULT_MAX_TESTS,
) -> list[MotifHit]:
    """Report PWM occurrences with a complete, explicit per-call test scope.

    Candidate sequences are already oriented into each proposed array end.
    Local coordinates index ``candidate.sequence + candidate.observed_repeat``.
    The hit strand is relative to that coordinate frame. Reverse matrices are
    evaluated against the same forward background, so asymmetric backgrounds
    receive the correct strand-specific null distribution.

    Ambiguous IUPAC windows are counted in the testing family with p=1 and are
    never emitted as occurrences. No testing or reporting cap silently drops a
    window: exceeding ``max_tests`` raises an actionable error before scoring.
    A model's optional ``junction_end_range`` is an inclusive, predeclared range
    for the interval end offset from the junction, on the candidate frame.
    """
    if type(max_tests) is not int or max_tests < 1:
        raise ValueError("max_tests must be a positive integer")
    _validate_candidates(candidates)
    models = [validate_motif(motif) for motif in motifs]
    if len({model["id"] for model in models}) != len(models):
        raise ReferenceError("motif testing scope contains duplicate model ids")
    total_tests = sum(
        len(_positions(candidate, model)) * len(model["strands"])
        for candidate in candidates for model in models
    )
    if total_tests > max_tests:
        raise InputError(
            f"motif testing scope has {total_tests} tests, exceeds bound {max_tests}; "
            "reduce the explicit analysis scope or increase the operator-selected bound"
        )
    scope: dict[str, Any] = {
        "adjustment": "Benjamini-Hochberg",
        "family": "all eligible model/candidate/position/strand tests in this scan call",
        "number_tests": total_tests,
        "background_assumption": "fixed independent identically distributed bases per model",
        "ambiguous_window_policy": "count with p=1; do not report an occurrence",
        "candidates": [
            {"id": candidate.id, "sequence_sha256": hashlib.sha256(
                candidate.junction_sequence.encode("ascii")).hexdigest(),
             "flank_length": len(candidate.sequence),
             "observed_repeat_length": len(candidate.observed_repeat)}
            for candidate in candidates
        ],
        "models": [{"id": model["id"], "domain": model["domain"],
                    "strands": model["strands"],
                    "model_sha256": hashlib.sha256(json.dumps(
                        model, sort_keys=True, separators=(",", ":"), allow_nan=False
                    ).encode()).hexdigest()} for model in models],
    }
    # Tuples are candidate, model, local start, strand, integer score, raw
    # floating log-odds, p-value; ambiguous positions remain None with p=1.
    tests: list[tuple[FlankCandidate, dict[str, Any], int, str, int, float, float] | None] = []
    pvalues: list[float] = []
    # Prepare matrix inputs once; null distributions retain their bounded cache.
    matrices = {(model["id"], strand): _matrix(model, strand)
                for model in models for strand in model["strands"]} if candidates else {}
    backgrounds = {model["id"]: tuple(float(model["background"][base]) for base in BASES)
                   for model in models} if candidates else {}
    ambiguous = 0
    for candidate in candidates:
        indices = tuple(BASE_INDEX.get(base, -1) for base in candidate.junction_sequence)
        # A prefix count recognizes ambiguous windows without rescanning bases.
        ambiguous_prefix = [0]
        for index in indices:
            ambiguous_prefix.append(ambiguous_prefix[-1] + (index < 0))
        for model in models:
            width = len(model["pwm"])
            positions = _positions(candidate, model)
            for strand in model["strands"]:
                table = _score_table(
                    matrices[model["id"], strand], backgrounds[model["id"]],
                    model.get("score_scale", 100),
                )
                for start in positions:
                    if ambiguous_prefix[start + width] != ambiguous_prefix[start]:
                        tests.append(None)
                        pvalues.append(1.0)
                        ambiguous += 1
                        continue
                    concrete = indices[start:start + width]
                    integer_score = sum(
                        table.integer_rows[position][index]
                        for position, index in enumerate(concrete)
                    )
                    raw_score = math.fsum(
                        table.real_rows[position][index]
                        for position, index in enumerate(concrete)
                    )
                    pvalue = table.pvalue(integer_score)
                    pvalues.append(pvalue)
                    tests.append((candidate, model, start, strand, integer_score, raw_score, pvalue))
    adjusted = benjamini_hochberg(pvalues)
    scope["number_ambiguous_tests"] = ambiguous
    scope_digest = hashlib.sha256(json.dumps(
        scope, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()).hexdigest()
    hits = []
    for test, qvalue in zip(tests, adjusted, strict=True):
        if test is None:
            continue
        candidate, model, start, strand, integer_score, raw_score, pvalue = test
        if pvalue > model.get("report_pvalue", 0.001):
            continue
        end = start + len(model["pwm"])
        boundary = len(candidate.sequence)
        # Slice the two frame components without rebuilding the complete context per hit.
        matched = (candidate.sequence[start:min(end, boundary)]
                   + candidate.observed_repeat[max(0, start - boundary):max(0, end - boundary)])
        model_strand: Literal["+", "-"] = "+" if strand == "+" else "-"
        hits.append(MotifHit(
            candidate_id=candidate.id, model_id=model["id"], domain=model["domain"],
            start=start, end=end, strand=model_strand,
            score=integer_score / model.get("score_scale", 100), pvalue=pvalue, qvalue=qvalue,
            details={
                "score_units": "quantized_log2_likelihood_ratio_bits",
                "unquantized_log2_likelihood_ratio_bits": raw_score,
                "score_scale": model.get("score_scale", 100),
                "integer_score": integer_score,
                "pvalue_statistic": "exact tail of rounded integer score under declared IID background",
                "rounding": "nearest integer per position, ties to even",
                "matched_sequence": matched,
                "motif_oriented_sequence": matched if strand == "+" else matched.translate(
                    COMPLEMENT)[::-1],
                "junction_end_offset": end - len(candidate.sequence),
                "background": dict(model["background"]),
                "testing_scope_id": scope_digest,
                "testing_scope": scope,
                "citations": list(model["citations"]),
                "interpretation": "sequence occurrence; not a probability of biological function",
            },
        ))
    return hits


__all__ = ["benjamini_hochberg", "scan_motifs"]
