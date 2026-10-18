import hashlib
import json

import pytest


def record(sequence, name="record", assembly=None, topology="linear"):
    return {"id": name, "sequence": sequence, "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "length": len(sequence), "assembly_id": assembly, "topology": topology}


def report(source, repeat_intervals=((3, 4), (9, 10)), category="Bona-fide", strand=None):
    repeats, spacers, deletions = [], [], []
    for ordinal, (start, end) in enumerate(repeat_intervals, 1):
        repeats.append({"index": ordinal, "sequence": source["sequence"][start-1:end],
                        "source_interval": {"start": start, "end": end, "length": end-start+1}})
        if ordinal < len(repeat_intervals):
            right = repeat_intervals[ordinal][0]-1
            if right == end:
                deletions.append({"index": ordinal, "sequence": "", "length": 0,
                                  "source_boundary": {"position": end, "coordinate_system": "0-based interbase", "reference": "input_forward_strand"}})
            else:
                spacers.append({"index": ordinal, "sequence": source["sequence"][end:right],
                                "source_interval": {"start": end+1, "end": right, "length": right-end}})
    start, end = repeat_intervals[0][0], repeat_intervals[-1][1]
    native = {"id": "arrayA", "source_interval": {"start": start, "end": end, "length": end-start+1},
              "source_sequence": source["sequence"][start-1:end], "category": category, "orientation": {"strand": strand},
              "repeat_count": len(repeats), "spacer_count": len(spacers)+len(deletions),
              "non_deleted_spacer_count": len(spacers), "spacer_deletion_count": len(deletions),
              "repeats": repeats, "spacers": spacers, "spacer_deletions": deletions}
    return {"schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
            "coordinate_system": {"indexing": 1, "end_inclusive": True, "orientation_is_separate": True, "reference": "input_forward_strand"},
            "source": {"id": source["id"], "length": source["length"], "sha256": source["sha256"], "sha256_scope": "uppercase_sequence"},
            "arrays": [native]}


@pytest.fixture
def native_fixture(tmp_path):
    source = record("TTACGATTACGGGGGATTTT")
    data = report(source)
    path = tmp_path / "report.json"
    path.write_text(json.dumps(data))
    fasta = tmp_path / "source.fasta"
    fasta.write_text(">record\n"+source["sequence"]+"\n")
    return source, data, path, fasta
