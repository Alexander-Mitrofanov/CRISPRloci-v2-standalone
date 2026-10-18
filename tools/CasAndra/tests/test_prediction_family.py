import json
from pathlib import Path

import pytest

from casandra.api import DEFAULT_MODEL_BUNDLE
from casandra.model import ProteinModel
from casandra.prediction import predict_query_records, profile_cas_family


@pytest.mark.parametrize(
    ("profile_id", "expected"),
    [
        ("C25_Cas2_1", "Cas2"),
        ("C25_Cas4_1", "Cas4"),
        ("C25_Csx32_1", "Csx32"),
    ],
)
def test_primary_family_is_available_without_system_classification(profile_id, expected):
    metadata = json.loads(
        (DEFAULT_MODEL_BUNDLE / "protein" / "profile_metadata.json").read_text(
            encoding="utf-8"
        )
    )

    assert metadata[profile_id]["labels"] == {"class": [], "subtype": [], "type": []}
    assert profile_cas_family(metadata, profile_id) == expected


def test_profile_family_requires_curated_primary_metadata():
    assert profile_cas_family({}, "missing") is None
    assert profile_cas_family({"profile": {"families": []}}, "profile") is None


def test_positive_prediction_emits_family_when_system_type_is_unavailable(monkeypatch):
    model = ProteinModel(
        directory=Path("unused"),
        manifest={"model_id": "test-model"},
        card={"report_evalue": 1e-3, "detection_threshold": 0.0},
        metadata={
            "C25_Cas2_1": {
                "families": ["cas2", "cas"],
                "intended_use": "primary_search",
                "labels": {"class": [], "type": [], "subtype": []},
            }
        },
        hmms=(),
    )
    monkeypatch.setattr(
        "casandra.prediction.scan_profiles_fast",
        lambda *_args, **_kwargs: [
            {
                "positive_score": 50.0,
                "negative_score": 0.0,
                "margin": 50.0,
                "best_positive_profile": "C25_Cas2_1",
                "classification": {"class": None, "type": None, "subtype": None},
                "hit_count": 1,
            }
        ],
    )

    row = predict_query_records(
        model, [{"sequence_id": "protein", "sequence": "MKTW"}], threads=1
    )[0]

    assert row["is_cas"] is True
    assert row["cas_family"] == "Cas2"
    assert row["result"] == "Cas2"
    assert row["best_positive_profile"] == "C25_Cas2_1"
    assert row["classification"] == {"class": None, "type": None, "subtype": None}
