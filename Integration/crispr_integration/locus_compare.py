"""Isolated CasAndra all-CDS extraction followed by bounded clinker comparison."""

from pathlib import Path
from .annotations import _payload, _strict_json, validate_casandra, json_lines
from .protospacer import _records
from .locus_compare_validation import expected_loci, validate_report, sha


def feature_document(records, rows_by_source):
    loci = []
    for index, (source, sequence) in enumerate(records.items(), 1):
        genes = []
        for row in rows_by_source[source]:
            if row["contig_id"] != source:
                raise ValueError("CDS belongs to another source")
            genes.append(
                {
                    "id": row["protein_id"],
                    "start": row["start_1based"],
                    "end": row["end_1based_inclusive"],
                    "strand": row["strand"],
                    "phase": row["phase"],
                    "translation_table": row["translation_table"],
                    "partial_5prime": row["partial_5prime"],
                    "partial_3prime": row["partial_3prime"],
                    "translation": row["protein_sequence"],
                }
            )
        genes.sort(key=lambda g: (g["start"], g["end"], g["id"]))
        loci.append({"id": f"locus{index:03d}", "source_id": source, "start": 1, "end": len(sequence), "genes": genes})
    doc = {
        "schema": {"name": "clinker.features", "version": "1.0.0"},
        "coordinate_system": "1-based-inclusive-source-forward",
        "loci": loci,
    }
    expected_loci(doc, records)
    return doc


def run_locus_compare_pipeline(input_path, output_path, *, runner_config, cpu=1, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    records = _records(Path(input_path).absolute())
    if not 2 <= len(records) <= 8 or any(len(s) > 50000 for s in records.values()):
        raise ValueError("Supply 2–8 genomic-neighborhood records, each at most 50,000 nt")
    job = _empty_job_root(output_path)
    logs = job / "logs"
    logs.mkdir()
    inputs = job / "locus-compare-input"
    inputs.mkdir()
    (inputs / "sources").mkdir()
    (job / "casandra").mkdir()
    combined = inputs / "input.fasta"
    combined.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()), encoding="ascii")
    result = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "workflow": "locus_compare",
        "status": "running",
        "tool_options": runner_config.tool_options,
        "stages": {},
    }
    marker = job / "integration_result.json"
    _atomic_json(marker, result)
    try:
        rows = {}
        counts = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
        for index, (source, sequence) in enumerate(records.items(), 1):
            fasta = inputs / "sources" / f"record-{index:06d}.fasta"
            fasta.write_text(f">{source}\n{sequence}\n", encoding="ascii")
            output = job / "casandra" / f"record-{index:06d}"
            _run_process(
                f"casandra-{index:06d}",
                [
                    *runner_config.casandra,
                    "predict-genome",
                    "--genome",
                    str(fasta),
                    "--output",
                    str(output),
                    "--gene-mode",
                    "auto",
                    "--threads",
                    "1",
                ],
                logs,
                timeout_seconds=timeout_seconds,
            )
            for key, value in validate_casandra(job, output, fasta).items():
                counts[key] += value
            rows[source] = json_lines(job, output / "proteins.jsonl")
        result["stages"]["casandra"] = {"status": "completed", **counts}
        features = feature_document(records, rows)
        _atomic_json(inputs / "features.json", features)
        identity = runner_config.tool_options.get("clinker", {}).get("--identity", 0.3)
        stage = _run_process(
            "clinker",
            [
                *runner_config.clinker,
                "compare",
                "--input",
                str(combined),
                "--features",
                str(inputs / "features.json"),
                "--outdir",
                str(job / "locus-compare"),
                "--identity",
                str(identity),
            ],
            logs,
            timeout_seconds=timeout_seconds,
        )
        report = _strict_json(_payload(job, job / "locus-compare/comparison.json"))
        stage["counts"] = validate_report(
            report,
            features,
            records,
            sha(_payload(job, combined)),
            sha(_payload(job, inputs / "features.json")),
            identity,
        )
        result["stages"]["clinker"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(marker, result)
        raise StageError("clinker", str(exc)) from exc
    _atomic_json(marker, result)
    return result
