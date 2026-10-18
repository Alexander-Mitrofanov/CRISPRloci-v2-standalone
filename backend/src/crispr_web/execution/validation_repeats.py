"""Independent repeat-output validation against immutable job inputs.

No imports from Integration or CRISPRrepeat: source slicing, orientation,
pair/motif arithmetic and effective settings are checked independently here.
"""

from __future__ import annotations

import hashlib
import math
import re

from ..service import CATEGORY_POLICY_CATEGORIES
from .validation_annotations import _COMPLEMENT, _json, _payload
from .validation_identify import _validate_identify_reports


def _hash(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def _number(value):
    return type(value) in {int, float} and math.isfinite(value)


def _close(actual, expected):
    return _number(actual) and abs(actual - expected) <= 1e-8


def _pairs(rows, sequence, checkpoint, *, global_model=False, span=None):
    length = len(sequence)
    if not isinstance(rows, list) or len(rows) > 250000:
        raise ValueError("Invalid pair count")
    pairs, sums = {}, [0.0] * length
    for number, row in enumerate(rows):
        if checkpoint and number % 1024 == 0:
            checkpoint()
        i, j, p = row
        if (
            type(i) is not int
            or type(j) is not int
            or not 0 <= i < j < length
            or not _number(p)
            or not 0 < p <= 1
            or j - i < 4
            or sequence[i] + sequence[j] not in {"AU", "UA", "GC", "CG", "GU", "UG"}
            or (i, j) in pairs
            or (span is not None and j - i + 1 > span)
        ):
            raise ValueError("Invalid pair value or coordinate")
        pairs[i, j] = p
        sums[i] += p
        sums[j] += p
    if global_model and any(x > 1 + 1e-8 for x in sums):
        raise ValueError("Global pair marginals exceed one")
    return pairs, sums


def _unpaired(values, n, sums=None):
    if not isinstance(values, list) or len(values) != n:
        raise ValueError("Incomplete accessibility")
    for i, value in enumerate(values):
        if not _number(value) or not 0 <= value <= 1:
            raise ValueError("Invalid accessibility")
        if sums is not None and not _close(value, max(0, 1 - sums[i])):
            raise ValueError("Inconsistent global unpaired probability")


def _isolated(result, sequence, checkpoint):
    if set(sequence) - set("ACGU"):
        if result != {"status": "unsupported", "reason": "ambiguous_bases"}:
            raise ValueError("Ambiguous folding must be unavailable")
        return None, None
    if (
        result.get("status") != "completed"
        or result.get("method") != "global_partition_function"
    ):
        raise ValueError("Missing isolated fold")
    pairs, sums = _pairs(result["pairs"], sequence, checkpoint, global_model=True)
    _unpaired(result["unpaired"], len(sequence), sums)
    if not all(_number(result.get(k)) for k in ("mfe_kcal_mol", "ensemble_kcal_mol")):
        raise ValueError("Invalid energy")
    structure = result.get("structure")
    if not isinstance(structure, str) or len(structure) != len(sequence):
        raise ValueError("Invalid structure length")
    stack, motif = [], []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")" and stack:
            i = stack.pop()
            if index - i < 4 or sequence[i] + sequence[index] not in {
                "AU",
                "UA",
                "GC",
                "CG",
                "GU",
                "UG",
            }:
                raise ValueError("MFE pair violates model")
            motif.append([i, index])
        elif char != ".":
            raise ValueError("Unbalanced structure")
    if stack or sorted(motif) != result["motif_pairs"]:
        raise ValueError("MFE pair set mismatch")
    return pairs, sorted(motif)


def validate_repeat_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    if input_root is None:
        return False
    try:
        molecule = claimed.molecule
        options = claimed.tool_options.get("crisprrepeat", {})
        requested_window = options.get("--window", 150)
        requested_span = options.get("--span", 80)
        array_mode = claimed.mode == "repeat_context"
        if molecule not in {"DNA", "RNA"} or (array_mode and molecule != "DNA"):
            return False
        if (
            integration.get("workflow") != claimed.mode
            or integration.get("status") != "completed"
        ):
            return False
        sources = {}
        for path in sorted(input_root.iterdir()):
            lines = _payload(input_root, path, checkpoint).decode("ascii").splitlines()
            if (
                not lines
                or not lines[0].startswith(">")
                or any(x.startswith(">") for x in lines[1:])
            ):
                return False
            identity, seq = lines[0][1:], "".join(lines[1:])
            alphabet = set(
                "ACGURYSWKMBDHVN" if molecule == "RNA" else "ACGTRYSWKMBDHVN"
            )
            if (
                identity in sources
                or not seq
                or set(seq) - alphabet
                or expected.get(identity) != (len(seq), _hash(seq))
            ):
                return False
            sources[identity] = seq
        if set(sources) != set(expected):
            return False
        doc = _json(root, root / "repeat/repeat-analysis.json", checkpoint)
        if (
            doc.get("schema") != {"name": "CRISPRrepeat-report", "version": "1.0.0"}
            or doc.get("tool") not in ({"name": "CRISPRrepeat", "version": "1.0.0"},
                                       {"name": "CRISPRrepeat", "version": "1.0.1"})
            or doc.get("status") != "completed"
            or doc.get("molecule") != molecule
            or doc.get("input_kind") != ("arrays" if array_mode else "repeats")
            or doc.get("coordinate_system") != "zero-based half-open, source forward"
            or doc.get("orientation_policy")
            != ("forward" if molecule == "RNA" else "both")
        ):
            return False
        model = doc["model"]
        if any(
            model.get(k) != v
            for k, v in {
                "engine": "ViennaRNA",
                "version": "2.7.0",
                "parameters": "RNA_Turner2004",
                "dangles": 2,
                "salt_molar": 1.021,
                "temperature": options.get("--temperature", 37.0),
                "window": requested_window,
                "span": requested_span,
                "no_lonely_pairs": options.get("--no-lonely-pairs", False),
                "numeric_cutoff": 0.0,
                "minimum_loop_size": 3,
                "allow_GU": True,
                "pseudoknots": False,
            }.items()
        ) or not re.fullmatch(r"[0-9a-f]{64}", model.get("parameter_sha256", "")):
            return False
        if doc["sources"] != [
            {"id": k, "length": len(v), "sha256": _hash(v)} for k, v in sources.items()
        ]:
            return False
        categories = list(CATEGORY_POLICY_CATEGORIES[claimed.category_policy])
        strands = ["forward"] if molecule == "RNA" else ["forward", "reverse"]
        arrays, desired = {}, []
        if not array_mode and (
            doc["arrays"] or doc["contexts"] or doc["consensus_folds"]
        ):
            return False
        if doc["categories"] != (categories if array_mode else []):
            return False
        if array_mode:
            if (
                integration.get("stages", {}).get("crispridentify", {}).get("status")
                != "completed"
            ):
                return False
            if not _validate_identify_reports(root, expected, checkpoint):
                return False
            raw_arrays = []
            for path in sorted((root / "crispridentify").rglob("report.json")):
                report = _json(root, path, checkpoint)
                raw_arrays.extend(
                    (report["source"]["id"], row)
                    for row in report["arrays"]
                    if row["category"] in categories
                )
            if len(raw_arrays) != len(doc["arrays"]):
                return False
            for number, (array, (identity, raw)) in enumerate(
                zip(doc["arrays"], raw_arrays, strict=True), 1
            ):
                seq = sources[identity]
                start, end = (
                    raw["source_interval"]["start"] - 1,
                    raw["source_interval"]["end"],
                )
                if (
                    type(raw["source_interval"]["start"]) is not int
                    or type(end) is not int
                    or not 0 <= start < end <= len(seq)
                    or end - start > 5000
                    or raw["source_sequence"] != seq[start:end]
                ):
                    return False
                repeats, previous = [], start
                consensus = raw.get("consensus_repeat")
                if (
                    not isinstance(consensus, str)
                    or not 1 <= len(consensus) <= 200
                    or set(consensus) - set("ACGTRYSWKMBDHVN")
                    or not isinstance(raw.get("repeats"), list)
                    or type(raw.get("repeat_count")) is not int
                ):
                    return False
                for index, repeat in enumerate(raw["repeats"], 1):
                    rs, re_ = (
                        repeat["source_interval"]["start"] - 1,
                        repeat["source_interval"]["end"],
                    )
                    if (
                        type(repeat.get("index")) is not int
                        or repeat.get("index") != index
                        or type(repeat["source_interval"]["start"]) is not int
                        or type(re_) is not int
                        or not previous <= rs < re_ <= end
                        or re_ - rs > 200
                        or repeat["sequence"] != seq[rs:re_]
                    ):
                        return False
                    repeats.append(
                        {
                            "index": index,
                            "start": rs,
                            "end": re_,
                            "sequence": seq[rs:re_],
                        }
                    )
                    previous = re_
                orientation = (raw.get("orientation") or {}).get("strand")
                bound = {
                    "id": f"array_{number:06d}",
                    "native_id": raw["id"],
                    "source_id": identity,
                    "start": start,
                    "end": end,
                    "category": raw["category"],
                    "source_sequence": seq[start:end],
                    "consensus_repeat": raw["consensus_repeat"],
                    "repeats": repeats,
                    "detector_strand": orientation
                    if orientation in {"+", "-"}
                    else None,
                }
                if array != bound or not repeats or raw["repeat_count"] != len(repeats):
                    return False
                arrays[array["id"]] = array
                for strand in strands:
                    desired.extend(
                        (identity, array, repeat, strand) for repeat in repeats
                    )
        else:
            desired = [
                (identity, None, None, strand)
                for identity in sources
                for strand in strands
            ]
        if sum(len(array["source_sequence"]) for array in arrays.values()) > 20000:
            return False
        if (
            len(desired) > 2000
            or len(doc["instances"]) != len(desired)
            or len(doc["contexts"]) != len(arrays) * len(strands)
        ):
            return False
        contexts, context_keys = {}, set()
        for number, context in enumerate(doc["contexts"], 1):
            key = context["array_id"], context["orientation"]
            if (
                key in context_keys
                or key[1] not in strands
                or context["id"] != f"context_{number:06d}"
            ):
                return False
            context_keys.add(key)
            seq = arrays[key[0]]["source_sequence"]
            if key[1] == "reverse":
                seq = seq.translate(_COMPLEMENT)[::-1]
            seq = seq.replace("T", "U")
            if context["rna_sequence"] != seq or context["rna_sha256"] != _hash(seq):
                return False
            window, span = min(len(seq), requested_window), min(len(seq), requested_span)
            if (
                context.get("window"),
                context.get("span"),
                context.get("requested_window"),
                context.get("requested_span"),
            ) != (window, span, requested_window, requested_span):
                return False
            if context.get("method") != "local_window_pair_averages":
                return False
            if set(seq) - set("ACGU"):
                if (
                    context.get("status") != "unsupported"
                    or context.get("reason") != "ambiguous_bases"
                ):
                    return False
                q = {}
            elif len(seq) < 5:
                if (
                    context.get("reason") != "context_too_short"
                    or context.get("status") != "unsupported"
                ):
                    return False
                q = {}
            else:
                if context.get("status") != "completed":
                    return False
                q, _ = _pairs(context["pairs"], context["rna_sequence"], checkpoint, span=span)
                _unpaired(context["unpaired"], len(seq))  # No local row-sum constraint.
            if context["status"] != "completed" and (
                "pairs" in context or "unpaired" in context
            ):
                return False
            contexts[context["id"]] = context, q
        if len(doc["consensus_folds"]) != len(context_keys):
            return False
        consensus_seen = set()
        for row in doc["consensus_folds"]:
            key = row["array_id"], row["orientation"]
            if key not in context_keys or key in consensus_seen:
                return False
            consensus_seen.add(key)
            seq = arrays[key[0]]["consensus_repeat"]
            seq = (
                seq.translate(_COMPLEMENT)[::-1] if key[1] == "reverse" else seq
            ).replace("T", "U")
            if row["rna_sequence"] != seq:
                return False
            _isolated(row["isolated"], seq, checkpoint)
        for number, (row, (identity, array, repeat, strand)) in enumerate(
            zip(doc["instances"], desired, strict=True), 1
        ):
            if type(row.get("terminal")) is not bool or any(
                row.get(key) is not None and type(row[key]) is not int
                for key in ("repeat_index", "start", "end", "context_start")
            ):
                return False
            sequence = repeat["sequence"] if repeat else sources[identity]
            rna = (
                sequence.translate(_COMPLEMENT)[::-1]
                if strand == "reverse"
                else sequence
            ).replace("T", "U")
            offset = (
                (
                    repeat["start"] - array["start"]
                    if strand == "forward"
                    else array["end"] - repeat["end"]
                )
                if repeat
                else None
            )
            required = {
                "id": f"repeat_{number:06d}",
                "source_id": identity,
                "array_id": array["id"] if array else None,
                "repeat_index": repeat["index"] if repeat else None,
                "orientation": strand,
                "start": repeat["start"] if repeat else None,
                "end": repeat["end"] if repeat else None,
                "source_sequence": sequence,
                "rna_sequence": rna,
                "rna_sha256": _hash(rna),
                "context_start": offset,
                "terminal": bool(
                    repeat and repeat["index"] in {1, len(array["repeats"])}
                ),
            }
            if any(row.get(k) != v for k, v in required.items()) or len(rna) > 200:
                return False
            context, q = contexts[row["context_id"]] if array else (None, {})
            if array and (
                context["array_id"] != array["id"] or context["orientation"] != strand
            ):
                return False
            if not array and row["context_id"] is not None:
                return False
            p, motif = _isolated(row["isolated"], rna, checkpoint)
            comp = row["comparison"]
            reason, iso, ctx = None, None, None
            if p is None:
                reason = "ambiguous_bases"
            elif not motif:
                reason = "no_reference_pairs"
            else:
                iso = sum(p.get(tuple(pair), 0) for pair in motif) / len(motif)
                if context is None:
                    reason = "context_not_supplied"
                elif context["status"] != "completed":
                    reason = context["reason"]
                elif any(j - i + 1 > context["span"] for i, j in motif):
                    reason = "motif_exceeds_context_span"
                else:
                    ctx = sum(
                        q.get((i + offset, j + offset), 0) for i, j in motif
                    ) / len(motif)
            if (
                comp.get("motif_origin") != "isolated_mfe_hypothesis"
                or comp.get("motif_pairs") != (motif or [])
                or comp.get("reason") != reason
                or comp.get("status") != ("unsupported" if reason else "completed")
            ):
                return False
            for key, value in {
                "isolated_pair_support": iso,
                "context_pair_support": ctx,
                "support_change": ctx - iso if ctx is not None else None,
            }.items():
                if (value is None and comp.get(key) is not None) or (
                    value is not None and not _close(comp.get(key), value)
                ):
                    return False
            # The initial public workflow has no configured reference/classifier.
            if row["references"] != {
                "status": "unsupported",
                "reason": "reference_not_configured",
                "hits": [],
            }:
                return False
        if doc["classification"] != {
            key: {"status": "unsupported", "reason": f"no_{label}_{key}_model"}
            for key, label in [
                ("strand", "calibrated"),
                ("subtype", "validated"),
                ("family", "validated"),
            ]
        }:
            return False
        counts = {
            "sources": len(sources),
            "arrays": len(arrays),
            "instances": len(desired),
            "contexts": len(contexts),
        }
        stage = integration["stages"]["crisprrepeat"]
        return (
            all(type(value) is int for value in doc["counts"].values())
            and doc["counts"] == counts
            and stage.get("status") == "completed"
            and all(
                type(stage.get(k)) is int and stage[k] == v for k, v in counts.items()
            )
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        IndexError,
        AttributeError,
        RecursionError,
        ZeroDivisionError,
    ):
        return False
