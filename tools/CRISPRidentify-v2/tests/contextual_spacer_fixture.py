"""Actual NC013407 sequence with stored vectors and synthetic local origins.

Features here are frozen feasibility outputs. The scorer replays the portable
model on those vectors; it does not independently extract features. Translating
coordinates does not confer the original detector's historical lineage.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from components.compact_model import predict_compact_scores
from components.components_detection_refinement import CrisprCandidate
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.localized_spacer_context import capture_pending_context
from components.module_evaluation import build_advanced_fuzzy_search_filter
from tests.test_candidate_evidence import PARAMETERS


def load_fixture():
    data = json.loads((Path(__file__).parent / "fixtures/nc013407_spacer_context.json").read_text())
    source = data["source_window"]["sequence"]
    if hashlib.sha256(source.encode()).hexdigest() != data["source_window"]["sha256"]:
        raise ValueError("context test source hash differs")
    return data


def candidate_key(candidate):
    return (tuple(candidate.list_repeat_starts), tuple(candidate.list_repeats),
            tuple(candidate.list_repeats_gaped), tuple(candidate.list_spacers))


def local_candidate(assessment, offset):
    raw = assessment["candidate"]
    candidate = CrisprCandidate(
        list(raw["repeats"]), list(raw["gapped_repeats"]), list(raw["spacers"]),
        [start - offset for start, _end in raw["source_intervals"]],
    )
    if candidate.consensus != raw["consensus"]:
        raise ValueError("fixture consensus differs after local construction")
    candidate.source_path_classification_eligible = assessment[
        "classification_eligible_at_observed_stage"]
    if candidate.source_path_classification_eligible is not True:
        candidate.source_path_classification_reason = "recorded_ineligible_fixture_path"
    return candidate


def fixture(*, boundary_tie=False, ordinals=None):
    data = load_fixture()
    source = data["source_window"]["sequence"]
    source_hash = hashlib.sha256(source.encode()).hexdigest()
    offset = data["source_window"]["absolute_start"]
    ordinals = list(ordinals) if ordinals is not None else ([258, 260] if boundary_tie else [260])
    structural = build_advanced_fuzzy_search_filter(PARAMETERS)
    vectors, pending = {}, []
    for ordinal in ordinals:
        result = data["cases"][str(ordinal)]["observed_feasibility_result"]
        for role, assessment in result["assessments"].items():
            candidate = local_candidate(assessment, offset)
            vector = np.array([[assessment["feature_values"][name]
                                for name in COMPACT_FEATURE_SCHEMA]])
            score = float(predict_compact_scores(vector)[0])
            if score != assessment["own_score"]:
                raise ValueError("fixture vector-to-model score differs")
            key = candidate_key(candidate)
            if key in vectors and not np.array_equal(vectors[key], vector):
                raise ValueError("fixture has conflicting vectors for one representation")
            vectors[key] = vector
            if role == "full":
                capture = capture_pending_context(
                    source, source_hash, (0, len(source)), candidate, structural)
                if capture is not None:
                    pending.append(capture)

    def scorer(requests):
        return {key: [[float(predict_compact_scores(vectors[candidate_key(candidate)])[0]),
                       candidate, [vectors[candidate_key(candidate)].copy()]]
                      for candidate in candidates] for key, candidates in requests.items()}

    return source, [{}, {}, {}, {}, {}], pending, scorer
