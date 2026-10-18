"""Compare valid scientific decisions using isolated pre/post source imports."""

import copy
import json
import tempfile
from dataclasses import asdict, replace
from itertools import product
from pathlib import Path

from crisprtracrrna.adapters import nearest_cas, parse_casandra, parse_identify, sha256_file
from crisprtracrrna.config import Settings
from crisprtracrrna.evidence import (
    candidate_evidence,
    model_rescue_predictions,
    short_rescue_admission,
)
from crisprtracrrna.models import AntiRepeat, Array, Candidate, Hit, Interval, Sequence
from crisprtracrrna.search import parse_alignment_table
from crisprtracrrna.sequence import write_fasta

results = {
    "identify": [],
    "casandra": [],
    "alignment": [],
    "standard": [],
    "experimental": [],
    "model_rescue": [],
}
with tempfile.TemporaryDirectory() as directory_name:
    directory = Path(directory_name)
    source = Sequence("seq", "seq", "ACGT" * 100)
    fasta = directory / "input.fa"
    write_fasta(fasta, [source])
    for strand in ("+", "-", None):
        report = {
            "schema": {"name": "CRISPRidentify-report", "version": "1.9.0"},
            "coordinate_system": {
                "indexing": 1,
                "end_inclusive": True,
                "reference": "input_forward_strand",
            },
            "source": {"id": source.id, "length": len(source.bases), "sha256": source.sha256},
            "validation": {"status": "valid"},
            "optional": {"future": [1e100, "résultat"]},
            "arrays": [
                {
                    "id": "array1",
                    "category": "Bona-fide",
                    "certainty_score": 0.9,
                    "source_interval": {"start": 1, "end": 4},
                    "source_sequence": "ACGT",
                    "consensus_repeat": "acgu",
                    "orientation": {"strand": strand},
                    "validation": {"status": "valid"},
                },
                {"category": "Alternative"},
            ],
        }
        path = directory / "identify.json"
        path.write_text(json.dumps(report), encoding="utf-8")
        arrays, preserved = parse_identify(path, source)
        results["identify"].append({"arrays": [asdict(a) for a in arrays], "report": preserved})
    genes = [
        {
            "contig_id": "seq",
            "start_1based": 1,
            "end_1based_inclusive": 9,
            "prediction": {"is_cas": True},
            "optional": [1e100, "第二"],
        },
        {
            "contig_id": "seq",
            "start_1based": 20,
            "end_1based_inclusive": 29,
            "prediction": {"is_cas": False},
        },
    ]
    cassettes = [
        {
            "cassette_id": f"c{i}",
            "contig_id": "seq",
            "start_1based": 1 + i * 30,
            "end_1based_inclusive": 12 + i * 30,
            "final_classification": {"type": kind},
            "optional": 1e-100,
        }
        for i, kind in enumerate(("II", "II", "V"))
    ]
    for filename, rows in (("proteins.jsonl", genes), ("cassettes.jsonl", cassettes)):
        (directory / filename).write_text(
            "\n\n".join(json.dumps(r, ensure_ascii=False) for r in rows), encoding="utf-8"
        )
    run = {
        "schema_version": 5,
        "program": "CasAndra",
        "genes": 2,
        "cas_cassettes": 3,
        "inputs": [
            {"kind": "protein_fasta", "sha256": "x"},
            {"kind": "genome_fasta", "sha256": sha256_file(fasta)},
        ],
    }
    (directory / "run.json").write_text(json.dumps(run))
    manifest = {
        "schema_version": 5,
        "files": {
            name: {"sha256": sha256_file(directory / name)}
            for name in ("proteins.jsonl", "cassettes.jsonl")
        },
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))
    cas = parse_casandra(directory, fasta, source)
    results["casandra"].append(cas)
    for kind in ("II", "V", "other"):
        results["casandra"].append(nearest_cas(Interval(15, 25), cas, kind))

    array = Array("array", Interval(250, 300), "T" * 40, "+", "Possible", 0.8)
    lines = []
    for length in (11, 12, 13, 23, 24, 25, 40):
        for identity in (69, 70, 79, 80, 90, 100):
            for reverse in (False, True):
                start, end = (30 + length - 1, 30) if reverse else (30, 30 + length - 1)
                lines.append(f"q seq {identity} {length} 0 0 1 {length} {start} {end} 1 5\n")
    lines.append("q seq 100 40 0 0 1 40 250 289 1 5\n")
    align = directory / "alignment.tsv"
    align.write_text("".join(lines))
    for settings in (Settings(), Settings(coverage=0.7, similarity=0.9)):
        for short in (False, True):
            hits = parse_alignment_table(
                align, {"q": (array, array.repeat, "+")}, source, [array], settings, short=short
            )
            results["alignment"].append([asdict(h) for h in hits])

    for strand in ("+", "-"):
        anti = AntiRepeat(array, array.repeat, "+", Interval(100, 120), strand, 0.9, 0.5, True)
        candidate = Candidate("c", anti, Interval(95, 180), "A" * 85, Interval(5, 25))
        model = Hit("c", Interval(25, 60), "+", 5, 0.01)
        term = Hit("c", Interval(25, 50), "+", 5, 0.001)
        duplex = {
            "target_interval": Interval(5, 25).public(),
            "energy_kcal_mol": -4,
            "target_paired_positions": list(range(6, 26)),
            "repeat_paired_positions": list(range(1, 21)),
        }
        for models, terms in (
            ([], []),
            ([model], []),
            ([], [term]),
            ([model, replace(model, strand="-")], [term, replace(term, strand="-")]),
        ):
            for interaction in (None, duplex):
                results["standard"].append(
                    candidate_evidence(
                        candidate, source, models, terms, interaction, cas, Settings()
                    )
                )
        for settings in (Settings(), Settings(short_anti_repeat_rescue=True)):
            scenarios = product(
                (9, 10, 11, 20),
                (-0.0001, 0, 0.0001),
                (0.01, 0.01001),
                (15, 16, 25),
                (0.001, 0.001001),
                (24, 25),
            )
            for paired, energy, model_evalue, model_start, term_evalue, term_start in scenarios:
                interaction = copy.deepcopy(duplex)
                interaction["energy_kcal_mol"] = energy
                interaction["target_paired_positions"] = list(range(6, 6 + paired))
                interaction["repeat_paired_positions"] = list(range(1, paired + 1))
                results["experimental"].append(
                    short_rescue_admission(
                        candidate,
                        [replace(model, evalue=model_evalue, interval=Interval(model_start, 60))],
                        [replace(term, evalue=term_evalue, interval=Interval(term_start, 50))],
                        interaction,
                        settings,
                    )
                )
            results["experimental"].append(
                short_rescue_admission(candidate, [model], [term], None, settings)
            )
        for hits in ([], [Hit("seq", Interval(100, 130), strand, 20, 1e-8, "custom")]):
            results["model_rescue"].append(
                model_rescue_predictions(source, hits, [array], [], cas, Settings())
            )
print(json.dumps(results, ensure_ascii=False, separators=(",", ":")))
