import json
import tempfile
from dataclasses import asdict
from itertools import product
from pathlib import Path

from crisprtracrrna.adapters import load_json, nearest_cas
from crisprtracrrna.config import Settings
from crisprtracrrna.evidence import short_rescue_admission
from crisprtracrrna.models import AntiRepeat, Array, Candidate, Hit, Interval, Sequence
from crisprtracrrna.search import parse_alignment_table

result = {"alignment": [], "fraction": [], "cas_ties": [], "finite_json": []}
source = Sequence("genome", "genome", "ACGT" * 100)
array = Array("array", Interval(250, 300), "T" * 40, "+", "Possible", 0.8)
settings = Settings()
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "table.tsv"
    for inverted, qreverse, sreverse, span, identity, short in product(
        (False, True),
        (False, True),
        (False, True),
        (11, 12, 13, 23, 24, 25, 40),
        (69.99999999999999, 70.0, 70.00000000000001, 79.99999999999999, 80.0, 80.00000000000001),
        (False, True),
    ):
        qstart, qend = (span, 1) if qreverse else (1, span)
        sstart, send = (100 + span - 1, 100) if sreverse else (100, 100 + span - 1)
        if inverted:
            ids = "genome repeat"
            coords = f"{sstart} {send} {qstart} {qend}"
        else:
            ids = "repeat genome"
            coords = f"{qstart} {qend} {sstart} {send}"
        path.write_text(f"{ids} {identity} {span} 0 0 {coords} 1 5\n")
        hits = parse_alignment_table(
            path,
            {"repeat": (array, array.repeat, "+")},
            source,
            [array],
            settings,
            genome_is_query=inverted,
            short=short,
        )
        valid = (
            (span / 40 < 0.6 and span >= 12 and identity / 100 >= 0.8)
            if short
            else (span / 40 >= 0.6 and identity / 100 >= 0.7)
        )
        assert len(hits) == (int(valid) if short else 2 * int(valid))
        if hits:
            assert hits[0].strand == ("-" if qreverse == sreverse else "+")
        result["alignment"].append([asdict(hit) for hit in hits])

    for anchor_length, paired, strand in product((19, 20, 21), (9, 10, 11), ("+", "-")):
        anti = AntiRepeat(
            array,
            array.repeat,
            "+",
            Interval(100, 100 + anchor_length),
            strand,
            0.9,
            anchor_length / 40,
            True,
        )
        candidate = Candidate(
            "candidate", anti, Interval(95, 180), "A" * 85, Interval(5, 5 + anchor_length)
        )
        interaction = {
            "target_interval": Interval(5, 5 + anchor_length).public(),
            "energy_kcal_mol": -0.0 if paired == 9 else -1e-300,
            "target_paired_positions": list(range(6, 6 + paired)),
            "repeat_paired_positions": list(range(1, 1 + paired)),
        }
        model = Hit("candidate", Interval(5 + anchor_length, 60), "+", 5, 0.01)
        admission = short_rescue_admission(candidate, [model], [], interaction, settings)
        assert admission["accepted"] is (paired >= 10 and paired / anchor_length >= 0.5)
        if paired == 9:
            assert admission["rejection_reasons"][0] == "nonnegative_energy"
        result["fraction"].append(admission)

    for rows in (["z", "a"], ["a", "z"]):
        cas = {
            "cassettes": [
                {
                    "cassette_id": identifier,
                    "start_1based": 1,
                    "end_1based_inclusive": 10,
                    "final_classification": {"type": "II"},
                }
                for identifier in rows
            ]
        }
        nearest = nearest_cas(Interval(20, 30), cas, "II")
        assert nearest["id"] == "a"
        result["cas_ties"].append(nearest)

    for literal in (
        "1.7976931348623157e308",
        "-1.7976931348623157e308",
        "5e-324",
        "-0.0",
        "1e-999",
    ):
        path.write_text('{"nested": [{"number": ' + literal + "}]}")
        result["finite_json"].append(load_json(path))

print(json.dumps(result, separators=(",", ":")))
