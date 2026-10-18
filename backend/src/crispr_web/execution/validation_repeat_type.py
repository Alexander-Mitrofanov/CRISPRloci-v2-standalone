"""Independent backend boundary for RepeatTyper; no Integration/tool imports."""
import hashlib
import json
import math

from .validation_annotations import _payload, _strict_json, _fasta

MODEL = {'id': 'crisprcastyper-1.9.0-repeat-model',
 'upstream_repository': 'https://github.com/Russel88/CRISPRCasTyper',
 'upstream_version': '1.9.0',
 'upstream_commit': 'a2d51831d4802533bb113bac23cf9df2853063b0',
 'upstream_wheel_sha256': '9d1113ef7f6e375253986bc60740942a404b0439e1333da6ed2520ac0224cca7',
 'files': {'xgb_repeats.json': 'f9170b07805d1bb89fa1c02cfb861836e42498738e458f5e9e6db181dfcbe57a',
           'type_dict.tab': 'b01468938d868e038a14c7a42f3187bd234e2d63d49dfc4c04f92094cb4112c8'},
 'labels': ['I-A',
            'I-B',
            'I-C',
            'I-D',
            'I-E',
            'I-F',
            'I-F_T',
            'I-G',
            'II-A',
            'II-B',
            'II-C',
            'III-A',
            'III-B',
            'III-C',
            'III-D',
            'III-E',
            'III-F',
            'IV-A1',
            'IV-A2',
            'IV-A3',
            'IV-D',
            'IV-E',
            'V-A',
            'V-B1',
            'V-E',
            'V-F1',
            'V-F2',
            'V-F3',
            'V-G',
            'V-I',
            'V-J',
            'V-K',
            'VI-A',
            'VI-B1',
            'VI-B2',
            'VI-C',
            'VI-D'],
 'license': 'MIT',
 'feature_contract': 'canonical-4mer-counts-GC-fraction-length',
 'iteration_range': [0, 41]}


def sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def same(a, b):
    return json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False)


def validate_repeat_type_pipeline(claimed, root, expected, integration, input_root, checkpoint=None):
    try:
        return _validate(claimed, root, expected, integration, input_root, checkpoint)
    except (OSError, ValueError, TypeError, KeyError, AttributeError, IndexError, RecursionError, UnicodeError):
        return False


def _validate(claimed, root, expected, integration, input_root, checkpoint):
    if input_root is None or claimed.molecule != "DNA" or integration.get("workflow") != "repeat_type":
        return False
    source = {}
    for path in sorted(input_root.iterdir()):
        key, sequence = _fasta(_payload(input_root, path, checkpoint))
        if key in source or expected.get(key) != (len(sequence), sha(sequence)) or not 1 <= len(sequence) <= 200:
            return False
        source[key] = sequence
    if set(source) != set(expected) or not 1 <= len(source) <= 1000:
        return False
    doc = _strict_json(_payload(root, root / "repeat-type/report.json", checkpoint))
    threshold = float(claimed.tool_options.get("repeattyper", {}).get("--threshold", 0.75))
    if (doc["schema"] != {"name": "repeattyper.classification", "version": "1.0.0"}
            or doc["tool"] not in ({"name": "RepeatTyper", "version": "2.0.0"},
                                   {"name": "RepeatTyper", "version": "2.0.1"})
            or doc["execution_completed"] is not True or not same(doc["model"], MODEL)
            or not same(doc["input"], {"sha256": sha("".join(f">{k}\n{v}\n" for k, v in source.items())), "record_count": len(source)})
            or not same(doc["policy"], {"threshold": threshold, "kmer": 4, "iteration_range": [0, 41],
                "score_semantics": "closed_set_model_score_not_calibrated_confidence", "orientation": "canonical_kmers_no_strand_inference"})
            or len(doc["results"]) != len(source)):
        return False
    seen, canonical = {}, {}
    counts = {"records": len(source), "distinct_classified_sequences": len({s for s in source.values() if set(s) <= set("ACGT")}),
              "predicted": 0, "uncertain": 0, "unsupported": 0}
    for row, (key, sequence) in zip(doc["results"], source.items(), strict=True):
        if checkpoint:
            checkpoint()
        if row["id"] != key or row["sequence"] != sequence or row["sequence_sha256"] != sha(sequence):
            return False
        score, winner = row["native_score"], row["winner"]
        if set(sequence) - set("ACGT"):
            if score is not None or winner is not None:
                return False
            status, reason = "unsupported", "ambiguous_iupac"
        else:
            if type(score) not in (float, int) or not math.isfinite(score) or not 0 <= score <= 1 or winner not in MODEL["labels"]:
                return False
            status = "predicted" if score >= threshold else "uncertain"
            reason = "below_threshold" if status == "uncertain" else None
        if row["status"] != status or row["reason"] != reason or row["prediction"] != (winner if status == "predicted" else None):
            return False
        if sequence in seen and seen[sequence] != (winner, score):
            return False
        seen[sequence] = winner, score
        if set(sequence) <= set("ACGT"):
            identity = min(sequence, sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])
            if identity in canonical and canonical[identity] != (winner, score):
                return False
            canonical[identity] = winner, score
        counts[status] += 1
    stage = integration["stages"]["repeattyper"]
    return (same(doc["counts"], counts) and same(stage["counts"], counts) and stage["status"] == "completed"
            and stage["model_id"] == MODEL["id"] and stage["model_sha256"] == MODEL["files"]["xgb_repeats.json"])
