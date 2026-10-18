"""Source-bound Cas context shared by independent mobile/defense stages."""


def run_cas_context(job, records, config, timeout_seconds):
    from .annotations import validate_casandra
    from .runner import _run_process

    inputs = job / "context-input"
    inputs.mkdir()
    (job / "casandra").mkdir()
    counts = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
    for index, (source, sequence) in enumerate(records.items(), 1):
        path = inputs / f"record-{index:06d}.fasta"
        path.write_text(f">{source}\n{sequence}\n", encoding="ascii")
        output = job / "casandra" / f"record-{index:06d}"
        _run_process(
            f"casandra-{index:06d}",
            [*config.casandra, "predict-genome", "--genome", str(path),
             "--output", str(output), "--gene-mode", "auto", "--threads", "1"],
            job / "logs", timeout_seconds=timeout_seconds,
        )
        for key, value in validate_casandra(job, output, path).items():
            counts[key] += value
    return {"status": "completed", **counts}
