"""Narrow public views of the validated CasAndra and tracrRNA native reports."""

from __future__ import annotations
from collections import Counter
import json

from .common import _json, _safe_file
from .detection import _array_view


def _rows(root, path, checkpoint):
    if _safe_file(root, path) is None:
        return []
    if path.stat().st_size > 20_000_000:
        raise ValueError("Annotation table exceeds summary read limit")
    rows = []
    with path.open() as stream:
        for line in stream:
            if checkpoint:
                checkpoint()
            if line.strip():
                rows.append(json.loads(line))
    return rows


def annotation_summary(
    root, *, max_rows, accepted_categories, warnings, checkpoint=None, max_arrays=None
):
    array_limit = max_rows if max_arrays is None else max_arrays
    cas = {"status": "not_requested", "genes": [], "cassettes": []}
    tracr = {"status": "not_requested", "candidates": [], "sources": []}
    detection = None
    try:
        directories = []
        source_mapping = {}
        report_path = root / "tracrrna/report.json"
        if _safe_file(root, report_path) is not None:
            report = _json(root, report_path, checkpoint=checkpoint)
            if report.get("status") != "complete":
                raise ValueError("tracrRNA workflow is incomplete")
            records = report["records"]
            mode = report["configuration"]["mode"]
            tracr.update(
                status="completed",
                mode=mode,
                model_type=report["configuration"]["model_type"],
                prediction_count=report["summary"]["prediction_count"],
                sequence_count=report["summary"]["sequence_count"],
                score_is_probability=False,
            )
            arrays, counts, selected_counts = [], Counter(), Counter()
            schemas = Counter()
            array_count = selected_count = 0
            for record in records:
                if checkpoint:
                    checkpoint()
                source = record["source"]
                source_mapping[record["key"]] = source["id"]
                tracr["sources"].append(
                    {
                        "id": source["id"],
                        "length": source["length"],
                        "sha256": source["sha256"],
                    }
                )
                for item in record["predictions"]:
                    if len(tracr["candidates"]) >= max_rows:
                        break
                    interval = item["interval"]
                    model = item.get("model_hit") or {}
                    anti = item.get("anti_repeat") or {}
                    interaction = item.get("interaction") or {}
                    array = item.get("array") or {}
                    tracr["candidates"].append(
                        {
                            "id": item["export_id"],
                            "source_id": source["id"],
                            "start": interval["start"],
                            "end": interval["end"],
                            "strand": item["strand"],
                            "sequence": item["sequence"],
                            "origin": item.get("candidate_origin", "covariance_model" if mode == "model" else None),
                            "evidence_class": item.get("evidence_class", "cm_only" if mode == "model" else None),
                            "boundary_kind": item.get("boundary_kind", "merged_model_region" if mode == "model" else None),
                            "score": (item.get("ranking") or {}).get(
                                "score", item.get("best_score")
                            ),
                            "score_kind": "ranking"
                            if item.get("ranking")
                            else "covariance_model_bits",
                            "evalue": model.get("evalue", item.get("best_evalue")),
                            "array_id": array.get("id"),
                            "anti_repeat_identity": anti.get("identity"),
                            "anti_repeat_coverage": anti.get("coverage"),
                            "interaction_energy": interaction.get("energy_kcal_mol"),
                        }
                    )
                if mode == "complete":
                    native = root / "tracrrna/work" / record["key"]
                    directories.append(native / "casandra")
                    for path in (native / "identify").rglob("report.json"):
                        if _safe_file(root, path) is None:
                            continue
                        identify = _json(root, path, checkpoint=checkpoint)
                        schemas[str(identify["schema"]["version"])] += 1
                        for item in identify["arrays"]:
                            view = _array_view(item, source["id"])
                            counts[view["category"]] += 1
                            array_count += 1
                            if view["category"] in accepted_categories:
                                selected_count += 1
                                selected_counts[view["category"]] += 1
                                if len(arrays) < array_limit:
                                    arrays.append(
                                        {**view, "selected_for_analysis": True}
                                    )
            tracr["truncated"] = len(tracr["candidates"]) < tracr["prediction_count"]
            if mode == "complete":
                detection = {
                    "source_count": len(records),
                    "report_count": len(records),
                    "array_count": array_count,
                    "selected_array_count": selected_count,
                    "unselected_array_count": array_count - selected_count,
                    "category_counts": dict(counts),
                    "selected_category_counts": dict(selected_counts),
                    "accepted_categories": sorted(accepted_categories),
                    "arrays": arrays,
                    "arrays_truncated": len(arrays) < selected_count,
                    "model_score_is_probability": False,
                    "report_schema_versions": dict(schemas),
                }
        elif (root / "casandra").is_dir():
            directories = sorted(
                path.parent for path in (root / "casandra").glob("*/run.json")
            )
        if directories:
            cas.update(
                status="completed", gene_count=0, cas_gene_count=0, cassette_count=0
            )
            for directory in directories:
                run = _json(root, directory / "run.json", checkpoint=checkpoint)
                cas["gene_count"] += run["genes"]
                cas["cas_gene_count"] += run["cas_proteins"]
                cas["cassette_count"] += run["cas_cassettes"]
                cas["model"] = run.get("bundle_id")
                for row in _rows(root, directory / "proteins.jsonl", checkpoint):
                    prediction = row["prediction"]
                    if not prediction["is_cas"] or len(cas["genes"]) >= max_rows:
                        continue
                    category = prediction["classification"]
                    cas["genes"].append(
                        {
                            "id": row["protein_id"],
                            "source_id": source_mapping.get(
                                row["contig_id"], row["contig_id"]
                            ),
                            "start": row["start_1based"],
                            "end": row["end_1based_inclusive"],
                            "strand": row["strand"],
                            "type": category.get("type"),
                            "subtype": category.get("subtype"),
                            "profile": prediction.get("best_positive_profile"),
                            "score_margin": prediction.get("score_margin"),
                        }
                    )
                for row in _rows(root, directory / "cassettes.jsonl", checkpoint):
                    if len(cas["cassettes"]) >= max_rows:
                        break
                    category = row["final_classification"]
                    cas["cassettes"].append(
                        {
                            "id": row["cassette_id"],
                            "source_id": source_mapping.get(
                                row["contig_id"], row["contig_id"]
                            ),
                            "start": row["start_1based"],
                            "end": row["end_1based_inclusive"],
                            "cas_gene_count": row["cas_gene_count"],
                            **{
                                key: category.get(key)
                                for key in [
                                    "class",
                                    "type",
                                    "subtype",
                                    "method",
                                    "confidence",
                                ]
                            },
                        }
                    )
            cas["truncated"] = (
                len(cas["genes"]) < cas["cas_gene_count"]
                or len(cas["cassettes"]) < cas["cassette_count"]
            )
    except (OSError, ValueError, KeyError, TypeError) as error:
        warnings.append(f"Annotation summaries unavailable: {type(error).__name__}")
    return cas, tracr, detection
