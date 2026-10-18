"""Single-model scoring for CRISPRidentify v2."""

import numpy as np

from components.compact_feature_extraction import extract_compact_feature_rows
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.compact_model import predict_compact_scores


def evaluate_compact_model(evaluation):
    """Score every candidate using only the compact ten-feature model."""

    rows = extract_compact_feature_rows(
        evaluation.dict_crispr_array_candidates
    )
    for key in evaluation.dict_crispr_array_candidates:
        evaluation.dict_scored_result[key] = []
        evaluation.dict_scored_result_with_all_vectors[key] = []
    if not rows:
        return

    matrix = np.vstack([row for _key, _candidate, row in rows])
    if matrix.shape[1] != len(COMPACT_FEATURE_SCHEMA):
        raise RuntimeError("compact feature matrix has the wrong width")
    scores = predict_compact_scores(matrix)
    for score, (key, candidate, feature_row) in zip(scores, rows):
        feature_vector = feature_row.reshape(1, -1)
        scored_row = [float(score), candidate, [feature_vector]]
        evaluation.dict_scored_result[key].append(scored_row)
        evaluation.dict_scored_result_with_all_vectors[key].append(
            (float(score), candidate, [feature_vector])
        )
