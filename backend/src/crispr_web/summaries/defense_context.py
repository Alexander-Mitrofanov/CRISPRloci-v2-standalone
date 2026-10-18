"""Bounded defensive-system preview with explicit model and sequence scope."""
from .common import _json, _safe_file

GENE_FIELDS = ("id", "source_id", "replicon", "position", "start", "end", "strand", "translation_table",
               "partial_begin", "partial_end", "initiator_methionine", "cds_sha256", "protein", "protein_sha256")
SYSTEM_FIELDS = ("id", "source_id", "type", "subtype", "model", "gene_ids", "wholeness", "score", "start", "end")
COMPONENT_FIELDS = ("system_id", "gene_id", "profile", "status", "score", "i_evalue", "profile_coverage")
HIT_FIELDS = ("gene_id", "profile", "i_evalue", "score", "profile_coverage", "sequence_coverage", "protein_start", "protein_end")


def public_report(doc, root=None, checkpoint=None):
    from .context_annotations import cas_neighbors
    return {"cas_context": cas_neighbors(root, doc["systems"], checkpoint) if root is not None else [], "schema": {"name": "crisprloci.defense-context", "version": "1.0.0"}, "execution_completed": True,
            "model": {k: doc["model"][k] for k in ("id", "version", "commit", "license", "scope", "manifest_sha256")},
            "policy": {k: doc["policy"][k] for k in ("coordinates", "topology", "gene_caller", "coverage_profile", "scope", "functional_activity")},
            "counts": {k: doc["counts"][k] for k in ("sources", "genes", "systems", "components", "hits")},
            "sources": [{k: row[k] for k in ("id", "length", "sha256")} for row in doc["sources"]],
            **{key: [{k: row[k] for k in fields} for row in doc[key]] for key, fields in (
                ("genes", GENE_FIELDS), ("systems", SYSTEM_FIELDS), ("components", COMPONENT_FIELDS), ("hits", HIT_FIELDS))}}


def defense_context_summary(root, *, max_rows, checkpoint=None):
    path = root / "defense-context/report.json"
    if _safe_file(root, path) is None:
        return {"available": False, "systems": [], "components": [], "hits": [], "genes": []}
    doc = public_report(_json(root, path, max_bytes=20_000_000, checkpoint=checkpoint), root, checkpoint)
    limit = max(0, min(max_rows, 500))
    # Keep native components and hit identity separate; previews never imply global absence.
    for key in ("systems", "components", "hits"):
        doc[f"{key}_truncated"] = len(doc[key]) > limit
        doc[key] = doc[key][:limit]
    wanted = {gene for system in doc["systems"] for gene in system["gene_ids"]}
    wanted.update(row["gene_id"] for key in ("components", "hits") for row in doc[key])
    doc["genes"] = [{k: v for k, v in gene.items() if k != "protein"} for gene in doc["genes"] if gene["id"] in wanted]
    doc["cas_context"] = doc["cas_context"][:limit]
    doc["available"] = True
    return doc
