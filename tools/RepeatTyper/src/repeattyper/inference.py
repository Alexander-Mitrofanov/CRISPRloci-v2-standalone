"""Validated inference with deduplication and occurrence-preserving projection."""
import hashlib
import math

import numpy as np
import xgboost as xgb

from . import __version__
from .errors import InputError, ModelError
from .features import FEATURE_NAMES, matrix
from .input import DNA, Record, validate_records
from .model import load_model


def classify(records: list[Record], *, input_sha256: str, threshold: float = 0.75) -> dict:
    """Return versioned JSON-compatible results; scores are closed-set model scores."""
    validate_records(records)
    if (not isinstance(threshold, (int, float)) or isinstance(threshold, bool)
            or not 0 <= threshold <= 1 or not math.isfinite(threshold)):
        raise InputError("Threshold must be a finite number from zero to one")
    if (not isinstance(input_sha256, str) or len(input_sha256) != 64
            or any(char not in "0123456789abcdef" for char in input_sha256)):
        raise InputError("Input SHA-256 is required")
    booster, provenance = load_model()
    unique_sequences = dict.fromkeys(record.sequence for record in records)
    sequences = [sequence for sequence in unique_sequences if DNA.issuperset(sequence)]
    sequence_digests = {sequence: hashlib.sha256(sequence.encode("ascii")).hexdigest()
                        for sequence in unique_sequences}
    predictions: dict[str, tuple[str, float]] = {}
    if sequences:
        try:
            scores = booster.predict(xgb.DMatrix(matrix(sequences), feature_names=FEATURE_NAMES),
                                     iteration_range=(0, 41))
        except xgb.core.XGBoostError as exc:
            raise ModelError("Classifier inference failed") from exc
        if (not isinstance(scores, np.ndarray) or scores.dtype.kind != "f"
                or scores.shape != (len(sequences) + 1, 37) or not np.isfinite(scores).all()
                or (scores < 0).any() or (scores > 1).any()
                or not np.allclose(scores.sum(axis=1), 1, rtol=0, atol=1e-5)):
            raise ModelError("Classifier emitted invalid score matrix")
        scores = scores[1:]
        for sequence, row in zip(sequences, scores, strict=True):
            winner_index = int(row.argmax())
            predictions[sequence] = provenance["labels"][winner_index], float(row[winner_index])
    results = []
    counts = dict.fromkeys(("predicted", "uncertain", "unsupported"), 0)
    for record in records:
        winner, score = predictions.get(record.sequence, (None, None))
        reason: str | None
        if score is None:
            status, reason = "unsupported", "ambiguous_iupac"
        elif score < threshold:
            status, reason = "uncertain", "below_threshold"
        else:
            status, reason = "predicted", None
        counts[status] += 1
        results.append({"id": record.id, "sequence": record.sequence,
                        "sequence_sha256": sequence_digests[record.sequence],
                        "status": status, "winner": winner, "native_score": score,
                        "prediction": winner if status == "predicted" else None,
                        "reason": reason})
    return {"schema": {"name": "repeattyper.classification", "version": "1.0.0"},
            "tool": {"name": "RepeatTyper", "version": __version__}, "execution_completed": True,
            "input": {"sha256": input_sha256, "record_count": len(records)},
            "model": provenance, "policy": {"threshold": threshold, "kmer": 4, "iteration_range": [0, 41],
            "score_semantics": "closed_set_model_score_not_calibrated_confidence", "orientation": "canonical_kmers_no_strand_inference"},
            "counts": {"records": len(records), "distinct_classified_sequences": len(sequences),
                       **counts},
            "results": results}
