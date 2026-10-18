"""Independent CCTK publication gate: recompute identities, set metrics and components."""

import hashlib
import json
from ..array_exchange import native_input
from .validation_annotations import _fasta, _payload, _strict_json

INTERPRETATION = "Observed sharing only; source-forward order is not inferred acquisition order or a phylogeny."


def sha(value):
    return hashlib.sha256(
        value.encode() if isinstance(value, str) else value
    ).hexdigest()


def recompute(native, options, digest):
    arrays = []
    memberships = {}
    original = native["arrays"]
    if options["order"] == "id":
        original = sorted(original, key=lambda a: a["id"])
    for array in original:
        slots = []
        n = 0
        members = set()
        for occurrence in array["spacers"]:
            sequence = occurrence["sequence"]
            symbol = None
            if sequence is not None:
                if options["identity"] == "reverse-complement":
                    reverse = "".join(
                        {"A": "T", "C": "G", "G": "C", "T": "A"}[base]
                        for base in reversed(sequence)
                    )
                    sequence = min(sequence, reverse)
                symbol = "s_" + sha(sequence)
                members.add(symbol)
                n += 1
            slots.append({**occurrence, "symbol": symbol})
        arrays.append(
            {
                **array,
                "spacers": slots,
                "observed_spacer_count": n,
                "evidence_status": "observed" if n else "no_observed_spacers",
            }
        )
        memberships[array["id"]] = members
    ids = list(memberships)
    edges = []
    partition = {key: {key} for key in ids}
    for index, a in enumerate(ids):
        for b in ids[index + 1 :]:
            common = memberships[a].intersection(memberships[b])
            if len(common) < options["min_shared"]:
                continue
            edges.append(
                {
                    "a": a,
                    "b": b,
                    "nshared": len(common),
                    "jaccard": len(common) / len(memberships[a].union(memberships[b])),
                    "shared_symbols": sorted(common),
                }
            )
            merged = partition[a] | partition[b]
            for node in merged:
                partition[node] = merged
    components = []
    done = set()
    for key in ids:
        if key not in done:
            group = [node for node in ids if node in partition[key]]
            components.append(group)
            done.update(group)
    return {
        "schema": {"name": "cctk.comparison", "version": "1.0.0"},
        "tool": {"name": "CCTK", "version": "2.0.1"},
        "options": options,
        "arrays": arrays,
        "edges": edges,
        "components": components,
        "isolates": [g[0] for g in components if len(g) == 1],
        "interpretation": INTERPRETATION,
        "input_sha256": digest,
    }


def validate_array_compare_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    try:
        if (
            input_root is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "array_compare"
        ):
            return False
        job_root = input_root.parent.parent
        raw = _payload(job_root, job_root / "array-map.json", checkpoint)
        if sha(raw) != claimed.array_map_sha256:
            return False
        mapping = _strict_json(raw)
        if set(mapping) != {"schema", "submitted_sha256", "document"} or mapping[
            "schema"
        ] != {"name": "crisprloci.array-map", "version": "1.0.0"}:
            return False
        native, records = native_input(mapping["document"])
        source = {}
        for path in sorted(input_root.iterdir()):
            key, sequence = _fasta(_payload(input_root, path, checkpoint))
            source[key] = sequence
        if (
            records != source
            or {k: (len(v), sha(v)) for k, v in source.items()} != expected
        ):
            return False
        # Reject rewritten native input and annotations even if a report is forged.
        if (
            _strict_json(
                _payload(root, root / "array-compare-input/map.json", checkpoint)
            )
            != mapping
        ):
            return False
        serialized = (
            json.dumps(native, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        if (
            _payload(root, root / "array-compare-input/arrays.json", checkpoint)
            != serialized.encode()
        ):
            return False
        opts = claimed.tool_options.get("cctk", {})
        options = {
            "identity": opts.get("--identity", "exact"),
            "min_shared": opts.get("--min-shared", 1),
            "order": opts.get("--order", "supplied"),
        }
        report = _strict_json(
            _payload(root, root / "array-compare/comparison.json", checkpoint)
        )
        if report.get("tool") not in ({"name": "CCTK", "version": "2.0.1"},
                                      {"name": "CCTK", "version": "2.0.2"}):
            return False
        computed = recompute(native, options, sha(serialized))
        # Both supported patch releases use the same observed-comparison schema
        # and metric definitions; independently reconstruct every data field.
        computed["tool"] = report["tool"]
        if json.dumps(report, sort_keys=True, allow_nan=False) != json.dumps(
            computed, sort_keys=True, allow_nan=False
        ):
            return False
        arrays = computed["arrays"]
        counts = {
            "arrays": len(arrays),
            "spacer_occurrences": sum(a["observed_spacer_count"] for a in arrays),
            "deletions": sum(
                len(a["spacers"]) - a["observed_spacer_count"] for a in arrays
            ),
            "edges": len(computed["edges"]),
            "components": len(computed["components"]),
            "isolates": len(computed["isolates"]),
        }
        stage = integration["stages"]["cctk"]
        return (
            stage["status"] == "completed"
            and stage["array_map_sha256"] == claimed.array_map_sha256
            and json.dumps(stage["counts"], sort_keys=True)
            == json.dumps(counts, sort_keys=True)
        )
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
