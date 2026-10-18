"""Independent source/TSV publication gate for the ISEScan worker stage."""

import json
import hashlib
from .validation_annotations import _fasta, _payload, _strict_json
from .validation_insertion_evidence import document

MODELS = {
    "clusters.faa.hmm": "f5de1bf92059e39f66490340f1d4e585df3c4b63c79e92c74cf387d23b0007db",
    "clusters.single.faa": "da7b2b67f57394cca9c8952c21cfbf5077aaf3aabf38cefd285f0d702e234295",
}


def validate_insertion_context_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    try:
        if (
            input_root is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "insertion_context"
        ):
            return False
        records = {}
        for path in sorted(input_root.iterdir()):
            ident, seq = _fasta(_payload(input_root, path, checkpoint))
            if ident in records:
                return False
            records[ident] = seq
        if {
            key: (len(seq), hashlib.sha256(seq.encode()).hexdigest())
            for key, seq in records.items()
        } != expected:
            return False
        payload = "".join(f">{key}\n{seq}\n" for key, seq in records.items()).encode(
            "ascii"
        )
        if _payload(root, root / "insertion-input/source.fasta", checkpoint) != payload:
            return False
        report = _strict_json(_payload(root, root / "isescan/isescan.json", checkpoint))
        options = report["options"]
        if (
            set(options) != {"threads", "include_partial"}
            or type(options["threads"]) is not int
            or not 1 <= options["threads"] <= 8
            or options["include_partial"] is not True
        ):
            return False
        computed = document(
            payload,
            _payload(root, root / "isescan/native.tsv", checkpoint),
            options,
            MODELS,
        )
        if json.dumps(computed, sort_keys=True, allow_nan=False) != json.dumps(
            report, sort_keys=True, allow_nan=False
        ):
            return False
        from .validation_context_annotations import validate_cas_context
        validate_cas_context(root, records, integration["stages"]["casandra"], checkpoint)
        elements = report["elements"]
        counts = {
            "sources": len(records),
            "elements": len(elements),
            "complete": sum(e["native_type"] == "c" for e in elements),
            "partial": sum(e["native_type"] == "p" for e in elements),
            "with_terminal_inverted_repeats": sum(
                e["terminal_inverted_repeats"] is not None for e in elements
            ),
        }
        stage = integration["stages"]["isescan"]
        return stage["status"] == "completed" and stage["counts"] == counts
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        AttributeError,
        IndexError,
        RecursionError,
        UnicodeError,
    ):
        return False
