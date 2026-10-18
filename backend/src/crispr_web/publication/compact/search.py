"""Compact reference searches with shared query/reference identities."""

from .common import fields, register, sources


def protospacer_data(read):
    doc = read.json("protospacer/analysis.json")
    spacers, matches = {}, []
    for row in doc["matches"]:
        read.step()
        register(
            spacers,
            row["spacer_id"],
            {
                "id": row["spacer_id"],
                "sequence": row["spacer_sequence"],
                "sha256": row["spacer_sha256"],
            },
        )
        matches.append(
            fields(row, "spacer_id target_id target_start target_end relative_strand")
        )
    return {
        "coordinate_system": "1-based-inclusive-target-forward",
        "outcome": doc["outcome"],
        "search_complete": doc["search_complete"],
        "reference": fields(doc["reference"], "name source_release database_sha256"),
        "policy": fields(
            doc["policy"],
            "method substitutions gaps strands topology array_overlap pam",
        ),
        "targets": sources(doc["targets"]),
        "spacers": list(spacers.values()),
        "matches": matches,
    }


def viral_data(read):
    doc = read.json("viral/result.json")
    subjects, query_sequences, matches = {}, {}, []
    for row in doc["matches"]:
        read.step()
        register(query_sequences, row["query_occurrence_id"], row["query_sequence"])
        register(
            subjects,
            row["subject_accession"],
            {
                "accession": row["subject_accession"],
                "title": row["subject_title"],
                "length": row["subject_length"],
                "sha256": row["subject_sha256"],
                "taxid": row["taxid"],
            },
        )
        matches.append(
            {
                **fields(
                    row,
                    "id subject_accession subject_start subject_end relative_strand substitutions evalue bitscore subject_sequence_forward",
                ),
                "query_id": row["query_occurrence_id"],
            }
        )
    queries = []
    for row in doc["queries"]:
        read.step()
        queries.append(
            {
                "id": row["occurrence_id"],
                "source_id": row["original_id"],
                **fields(
                    row, "length sha256 outcome skip_reasons subject_limit_saturated"
                ),
                **(
                    {"sequence": query_sequences[row["occurrence_id"]]}
                    if row["occurrence_id"] in query_sequences
                    else {}
                ),
            }
        )
    return {
        "coordinate_system": "1-based-inclusive-subject-forward",
        "retrieval_completeness": doc["retrieval_completeness"],
        "outcome": doc["outcome"],
        "reference": fields(
            doc["reference"], "name source_release content_sha256 reference_sha256"
        ),
        "policy": fields(
            doc["policy"],
            "profile max_substitutions strands topology array_overlap host_inference pam",
        ),
        "queries": queries,
        "subjects": list(subjects.values()),
        "matches": matches,
    }


def repeat_map_data(read):
    doc = read.json("repeat-mapping.json")
    references, results, tie_policies = {}, [], set()
    for row in doc["results"]:
        read.step()
        hits = []
        for hit in row["hits"]:
            register(
                references,
                hit["reference_id"],
                fields(hit, "reference_id family motif accession"),
            )
            tie_policies.add(hit["alignment_tie_policy"])
            hits.append(
                fields(
                    hit,
                    "reference_id orientation aligned_orientation distance similarity query_aligned reference_aligned",
                )
            )
        results.append(
            {
                **fields(
                    row,
                    "id sequence status reason best_distance best_hit_count hits_truncated",
                ),
                "annotations": {
                    key: fields(value, "status values missing_reference_count")
                    for key, value in row["annotations"].items()
                    if key in {"family", "motif"}
                },
                "hits": hits,
            }
        )
    return {
        "reference": fields(doc["reference"], "namespace content_sha256 file_sha256"),
        "policy": {
            **fields(doc["policy"], "max_distance top_k orientation_scope"),
            "alignment_tie_policies": sorted(tie_policies),
        },
        "references": list(references.values()),
        "results": results,
    }
