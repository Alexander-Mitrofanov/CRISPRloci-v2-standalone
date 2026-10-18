"""Independent source, native annotation, interval and spacer-search gate."""

from .validation_annotations import _fasta, _payload, _strict_json
from .validation_prophage_evidence import validate_report, sha
from .validation_source_spacer_evidence import context_document
from .validation_identify import _validate_identify_reports


def validate_prophage_context_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    try:
        if (
            input_root is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "prophage_context"
        ):
            return False
        records = {}
        for path in sorted(input_root.iterdir()):
            key, seq = _fasta(_payload(input_root, path, checkpoint))
            if key in records:
                return False
            records[key] = seq
        if {
            k: (len(s), sha(s)) for k, s in records.items()
        } != expected or not _validate_identify_reports(root, expected, checkpoint):
            return False
        payload = "".join(f">{k}\n{s}\n" for k, s in records.items()).encode("ascii")
        if _payload(root, root / "phigaro-input/source.fasta", checkpoint) != payload:
            return False
        report = _strict_json(_payload(root, root / "phigaro/report.json", checkpoint))
        proteins = (
            _payload(root, root / "phigaro/work/proteins.faa", checkpoint)
            if any(len(s) >= 20000 for s in records.values())
            else b""
        )
        hits = (
            _payload(root, root / "phigaro/work/hits.tsv", checkpoint)
            if report["genes"]
            else b""
        )
        counts = validate_report(report, records, payload, proteins, hits)
        reports = [
            _strict_json(_payload(root, path, checkpoint))
            for path in sorted((root / "crispridentify").rglob("report.json"))
        ]
        native = _strict_json(
            _payload(root, root / "source-spacer/analysis.json", checkpoint)
        )
        context = context_document(native, reports, records, report["regions"])
        if (
            _strict_json(_payload(root, root / "phigaro/context.json", checkpoint))
            != context
        ):
            return False
        stages = integration["stages"]
        return (
            stages["phigaro"]["status"] == "completed"
            and stages["phigaro"]["counts"] == counts
            and stages["source_spacer"]["status"] == "completed"
            and stages["source_spacer"]["counts"] == context["search_counts"]
            and stages["crispridentify"]["status"] == "completed"
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
