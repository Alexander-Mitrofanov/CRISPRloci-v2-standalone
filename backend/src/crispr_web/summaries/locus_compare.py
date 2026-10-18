"""Bounded gene-neighborhood preview with separate native Cas evidence."""

from .common import _json, _safe_file


def public_report(root, checkpoint=None):
    report = _json(
        root,
        root / "locus-compare/comparison.json",
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    from .context_annotations import cas_neighbors
    return {
        "cas_context": cas_neighbors(root, report["loci"], checkpoint),
        "schema": {"name": "crisprloci.locus-comparison", "version": "1.0.0"},
        "execution_completed": True,
        **{
            key: report[key]
            for key in (
                "tool",
                "coordinate_system",
                "input_sha256",
                "features_sha256",
                "options",
                "limits",
                "sources",
                "loci",
                "links",
                "evaluated_protein_pairs",
                "alignment_cells",
                "interpretation",
            )
        },
        "counts": {
            "loci": len(report["loci"]),
            "genes": sum(len(locus["genes"]) for locus in report["loci"]),
            "links": len(report["links"]),
            "evaluated_protein_pairs": report["evaluated_protein_pairs"],
        },
    }


def locus_compare_summary(root, *, max_rows, checkpoint=None):
    if _safe_file(root, root / "locus-compare/comparison.json") is None:
        return {"available": False, "loci": [], "links": []}
    report = public_report(root, checkpoint)
    limit = min(max_rows, 100)
    report.update(available=True, links_truncated=len(report["links"]) > limit)
    report["links"] = [
        {
            key: value
            for key, value in link.items()
            if key not in {"aligned_query", "aligned_target"}
        }
        for link in report["links"][:limit]
    ]
    for locus in report["loci"]:
        for gene in locus["genes"]:
            gene["protein_length"] = len(gene.pop("translation"))
    return report
