"""Scientific record serialization preserves coordinates, uncertainty and snapshots."""
from crisprleader.models import ArrayAnnotation, RepeatAnchor, Segment, json_record


def test_record_snapshot_preserves_nested_values_without_sharing_containers() -> None:
    first = RepeatAnchor((Segment(4, 6),), "AC")
    last = RepeatAnchor((Segment(10, 12),), "GT")
    record = ArrayAnnotation(
        "array", "source", (Segment(4, 12),), (first, last), "AC",
        metadata={"score": None, "weights": (0.0, 1e-20), "audit": [{"certainty": "unknown"}]},
    )
    encoded = json_record(record)
    assert encoded == {
        "id": "array", "source_id": "source", "segments": [{"start": 4, "end": 12}],
        "terminal_repeats": [{"segments": [{"start": 4, "end": 6}], "sequence": "AC"},
                             {"segments": [{"start": 10, "end": 12}], "sequence": "GT"}],
        "consensus_repeat": "AC", "strand": "?", "strand_constraint": False, "category": None,
        "metadata": {"score": None, "weights": [0.0, 1e-20], "audit": [{"certainty": "unknown"}]},
    }
    encoded["metadata"]["audit"][0]["certainty"] = "edited export"
    encoded["terminal_repeats"][0]["segments"][0]["start"] = 0
    assert record.metadata["audit"] == [{"certainty": "unknown"}]
    assert first.segments == (Segment(4, 6),)
