"""Cas, tracrRNA and flank evidence without unrequested module placeholders."""

from .common import array_data, fields, register, sources


def cas_data(read):
    doc = read.json("annotations.json")["casandra"]
    if doc["status"] != "completed" or doc["truncated"]:
        raise ValueError("Cas output is incomplete")
    return {
        "coordinate_system": "1-based-inclusive-source-forward",
        "score_semantics": "uncalibrated-model-evidence",
        "evaluated_gene_count": doc["gene_count"],
        "genes": [
            fields(
                row, "id source_id start end strand type subtype profile score_margin"
            )
            for row in doc["genes"]
        ],
        "cassettes": [
            fields(
                row,
                "id source_id start end cas_gene_count class type subtype method confidence",
            )
            for row in doc["cassettes"]
        ],
    }


def tracr_data(read):
    doc = read.json("annotations.json")["tracrrna"]
    if doc["status"] != "completed" or doc["truncated"]:
        raise ValueError("tracrRNA output is incomplete")
    return {
        "coordinate_system": "1-based-inclusive-source-forward",
        "score_semantics": "uncalibrated-model-evidence",
        "method": doc["mode"],
        "model_type": doc["model_type"],
        "sources": sources(doc["sources"]),
        "candidates": [
            fields(
                row,
                "id source_id start end strand sequence origin evidence_class boundary_kind score score_kind evalue array_id anti_repeat_identity anti_repeat_coverage interaction_energy",
            )
            for row in doc["candidates"]
        ],
    }


def leader_data(read):
    doc = read.json("leader-contexts.json")
    if not doc["complete"] or doc["truncated"]:
        raise ValueError("flank output is incomplete")
    arrays, contexts = {}, []
    for row in doc["contexts"]:
        read.step()
        array = {
            "id": row["array_id"],
            "source_id": row["source_id"],
            "category": row["category"],
            "segments": row["array_segments"],
        }
        register(arrays, (array["source_id"], array["id"]), array)
        contexts.append(
            fields(
                row,
                "id source_id array_id side strand segments sequence truncated_reason observed_repeat repeat_segments",
            )
        )
    return {
        "coordinate_system": "zero_based_half_open",
        "prediction_performed": False,
        "context_semantics": "analysis-window-not-biological-boundary",
        "flank_length": doc["flank_length"],
        "sources": sources(doc["sources"]),
        "arrays": list(arrays.values()),
        "contexts": contexts,
    }


def loci_data(read):
    components = {
        "detection": array_data(read),
        "cas": cas_data(read),
        "tracrrna": tracr_data(read),
        "flanks": leader_data(read),
    }
    combined = {}
    for component in components.values():
        for source in component.pop("sources", []):
            previous = combined.get(source["id"], {})
            if any(
                previous[key] != value
                for key, value in source.items()
                if key in previous
            ):
                raise ValueError("inconsistent locus source")
            combined[source["id"]] = {**previous, **source}
    return {"sources": list(combined.values()), **components}
