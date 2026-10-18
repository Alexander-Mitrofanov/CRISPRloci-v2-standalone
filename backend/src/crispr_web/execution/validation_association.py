"""Independent SpacePHARER publication gate; no Integration or tool imports."""

import hashlib
import json
import math
import re
import sqlite3

from .validation_association_reference import reference_data

from .validation_annotations import _fasta, _payload, _strict_json

_json = _strict_json


def _records(path):
    files = sorted(path.iterdir()) if path.is_dir() else [path]
    result = {}
    for file in files:
        key, sequence = _fasta(_payload(file.parent, file))
        if key in result:
            raise ValueError("Duplicate sequence identity")
        result[key] = sequence
    return result


SUPPORTED_VERSIONS = {"2.1.0", "2.1.1"}
SHA = re.compile(r"^[0-9a-f]{64}$")
ID = re.compile(r"^[A-Za-z0-9_.:-]{1,160}$")


def same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(
        right, sort_keys=True, allow_nan=False
    )


def sha(data):
    return hashlib.sha256(data.encode() if isinstance(data, str) else data).hexdigest()


def validated_groups(document, records):
    if (
        not isinstance(document, dict)
        or set(document) != {"schema_version", "groups"}
        or document["schema_version"] != "spacepharer.groups/1.0.0"
    ):
        raise ValueError("Invalid group map")
    groups = document["groups"]
    if not isinstance(groups, list) or not 1 <= len(groups) <= 100:
        raise ValueError("Association requires 1–100 explicit groups")
    members, ids = [], set()
    for group in groups:
        if (
            set(group) != {"id", "members"}
            or not ID.fullmatch(group["id"])
            or group["id"] in ids
            or not isinstance(group["members"], list)
            or not group["members"]
        ):
            raise ValueError("Invalid association group")
        ids.add(group["id"])
        members.extend(group["members"])
    if len(members) != len(set(members)) or set(members) != set(records):
        raise ValueError("Group map must partition all occurrences")
    return groups


def number(value, minimum=0, maximum=1):
    if (
        type(value) not in (int, float)
        or not math.isfinite(value)
        or not minimum <= value <= maximum
    ):
        raise ValueError("Invalid scientific value")
    return value


def interval(sequence, start, end):
    if (
        type(start) is not int
        or type(end) is not int
        or not 1 <= min(start, end) <= max(start, end) <= len(sequence)
    ):
        raise ValueError("Invalid association coordinate")
    result = sequence[min(start, end) - 1 : max(start, end)]
    return (
        result
        if start <= end
        else result.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    )


def validate_report(
    doc, records, group_document, reference, targets, input_sha, groups_sha, fdr
):
    groups = validated_groups(group_document, records)
    if (
        set(doc)
        != {
            "schema_version",
            "tool_version",
            "input_sha256",
            "groups_sha256",
            "requested_fdr",
            "calibration_status",
            "calibration",
            "warnings",
            "reference",
            "records",
            "groups",
            "associations",
            "interpretation",
        }
        or doc["schema_version"] != "spacepharer.association/2.0.0"
        or doc["tool_version"] not in SUPPORTED_VERSIONS
        or doc["input_sha256"] != input_sha
        or doc["groups_sha256"] != groups_sha
        or number(doc["requested_fdr"], 0.001, 0.5) != fdr
        or not same(doc["reference"], reference)
        or not same(
            doc["records"],
            [{"id": key, "sequence": seq} for key, seq in records.items()],
        )
    ):
        raise ValueError("Association input or reference identity mismatch")
    expected_groups = []
    membership = {}
    for group in groups:
        unsupported = [
            key
            for key in group["members"]
            if len(records[key]) < 27 or set(records[key]) - set("ACGT")
        ]
        supported = set(group["members"]) - set(unsupported)
        expected_groups.append(
            {
                **group,
                "supported_occurrence_count": len(supported),
                "unique_sequence_count": len({records[key] for key in supported}),
                "unsupported_occurrence_ids": unsupported,
            }
        )
        for key in group["members"]:
            membership[key] = group["id"]
    if not same(doc["groups"], expected_groups):
        raise ValueError("Association group/unsupported count mismatch")
    supported_count = sum(row["supported_occurrence_count"] for row in expected_groups)
    status = doc["calibration_status"]
    calibration = doc["calibration"]
    if not supported_count:
        if status != "not_run" or calibration is not None or doc["associations"]:
            raise ValueError("Unsupported inputs cannot become negative evidence")
    else:
        if (
            status not in {"estimated", "unavailable"}
            or not isinstance(calibration, dict)
            or set(calibration)
            != {"schema_version", "status", "target_score_count", "control_score_count"}
            or calibration["schema_version"] != "spacepharer.calibration/1.0.0"
            or calibration["status"] != status
        ):
            raise ValueError("Missing native calibration evidence")
        if any(
            type(calibration[k]) is not int or not 0 <= calibration[k] <= 10000000
            for k in ("target_score_count", "control_score_count")
        ) or (
            status == "estimated"
            and min(
                calibration["target_score_count"], calibration["control_score_count"]
            )
            < 2
        ):
            raise ValueError("Invalid native calibration support")
    rows = doc["associations"]
    if not isinstance(rows, list) or len(rows) > 5000:
        raise ValueError("Association count exceeds budget")
    seen, total_hits = set(), 0
    group_lookup = {row["id"]: row["members"] for row in groups}
    for row in rows:
        if set(row) != {
            "group_id",
            "target_id",
            "combined_score",
            "native_fdr",
            "status",
            "native_hit_count",
            "hits",
        }:
            raise ValueError("Unexpected association evidence fields")
        pair = row["group_id"], row["target_id"]
        if pair in seen or pair[0] not in group_lookup or pair[1] not in targets:
            raise ValueError("Association group or target mismatch")
        seen.add(pair)
        number(row["combined_score"], 0, 1e300)
        if status != "estimated":
            if (
                row["native_fdr"] is not None
                or row["status"] != "uncalibrated_candidate"
            ):
                raise ValueError("Uncalibrated results cannot be accepted associations")
        else:
            native = number(row["native_fdr"])
            if row["status"] != (
                "association" if native <= fdr else "above_cutoff_candidate"
            ):
                raise ValueError("Association cutoff semantics mismatch")
        hits = row["hits"]
        if (
            not isinstance(hits, list)
            or not hits
            or type(row["native_hit_count"]) is not int
            or row["native_hit_count"] != len(hits)
        ):
            raise ValueError("Association support count mismatch")
        hit_seen = set()
        for hit in hits:
            total_hits += 1
            if total_hits > 10000 or set(hit) != {
                "occurrence_ids",
                "query_sequence",
                "p_best_hit",
                "query_start",
                "query_end",
                "target_start",
                "target_end",
                "query_alignment",
                "target_alignment",
            }:
                raise ValueError("Invalid association support shape or limit")
            sequence = hit["query_sequence"]
            ids = [key for key in group_lookup[pair[0]] if records[key] == sequence]
            if (
                not ids
                or hit["occurrence_ids"] != ids
                or len(sequence) < 27
                or set(sequence) - set("ACGT")
            ):
                raise ValueError("Support occurrence binding mismatch")
            number(hit["p_best_hit"])
            qa, ta = hit["query_alignment"], hit["target_alignment"]
            if (
                not isinstance(qa, str)
                or not isinstance(ta, str)
                or not 1 <= len(qa) <= 400
                or len(qa) != len(ta)
                or set(qa + ta) - set("ACGTN-")
                or any(left == right == "-" for left, right in zip(qa, ta, strict=True))
            ):
                raise ValueError("Invalid alignment")
            if qa.replace("-", "") != interval(
                sequence, hit["query_start"], hit["query_end"]
            ) or ta.replace("-", "") != interval(
                targets[pair[1]], hit["target_start"], hit["target_end"]
            ):
                raise ValueError("Alignment does not reconstruct source sequences")
            signature = (
                sequence,
                hit["query_start"],
                hit["query_end"],
                hit["target_start"],
                hit["target_end"],
            )
            if signature in hit_seen:
                raise ValueError("Duplicated support alignment")
            hit_seen.add(signature)
    if (
        doc["interpretation"]
        != "Spacer similarity supports a group–target association hypothesis, not active infection, cleavage or a PAM prediction."
    ):
        raise ValueError("Invalid association interpretation")
    if (
        not isinstance(doc["warnings"], list)
        or len(doc["warnings"]) > 10
        or any(
            not isinstance(x, str) or len(x) > 500 or any(ord(c) < 32 for c in x)
            for x in doc["warnings"]
        )
        or not isinstance(doc["interpretation"], str)
        or len(doc["interpretation"]) > 500
    ):
        raise ValueError("Invalid association narrative")
    return {
        "records": len(records),
        "groups": len(groups),
        "associations": sum(row["status"] == "association" for row in rows),
        "candidates": sum(row["status"] != "association" for row in rows),
        "support_hits": total_hits,
        "unsupported": len(records) - supported_count,
    }


def validate_association_pipeline(
    claimed, root, expected, integration, input_root, reference_root, checkpoint=None
):
    targets = None
    try:
        if (
            input_root is None
            or reference_root is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "spacer_association"
        ):
            return False
        records = _records(input_root)
        if not 1 <= len(records) <= 1000 or any(
            not 1 <= len(seq) <= 200 for seq in records.values()
        ):
            return False
        if {key: (len(seq), sha(seq)) for key, seq in records.items()} != expected:
            return False
        job_root = input_root.parent.parent
        group_raw = _payload(job_root, job_root / "association-groups.json", checkpoint)
        if sha(group_raw) != claimed.association_groups_sha256:
            return False
        groups = _json(group_raw)
        reference, targets = reference_data(
            reference_root, claimed.association_reference_sha256, checkpoint
        )
        if checkpoint:
            checkpoint()
        doc = _json(_payload(root, root / "association/report.json", checkpoint))
        fdr = float(claimed.tool_options.get("spacepharer", {}).get("--fdr", 0.05))
        counts = validate_report(
            doc,
            records,
            groups,
            reference,
            targets,
            sha("".join(f">{key}\n{seq}\n" for key, seq in records.items())),
            sha(group_raw),
            fdr,
        )
        stage = integration["stages"]["spacepharer"]
        return (
            stage["status"] == "completed"
            and same(stage["counts"], counts)
            and stage["reference_sha256"] == claimed.association_reference_sha256
            and stage["calibration_status"] == doc["calibration_status"]
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
        sqlite3.Error,
    ):
        return False
    finally:
        if targets is not None:
            targets.close()
