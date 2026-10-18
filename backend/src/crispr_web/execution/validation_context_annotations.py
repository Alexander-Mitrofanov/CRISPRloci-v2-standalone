"""Independent check of the Cas context produced alongside a new stage."""

from .validation_annotations import _cas, _payload


def validate_cas_context(root, records, stage, checkpoint=None):
    expected_directories = {f"record-{index:06d}" for index in range(1, len(records) + 1)}
    if ({path.name for path in (root / "casandra").iterdir()} != expected_directories
            or {path.name for path in (root / "context-input").iterdir()}
            != {name + ".fasta" for name in expected_directories}):
        raise ValueError("Unexpected Cas context source inventory")
    counts = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
    for index, (source, sequence) in enumerate(records.items(), 1):
        if checkpoint:
            checkpoint()
        payload = f">{source}\n{sequence}\n".encode("ascii")
        if _payload(root, root / "context-input" / f"record-{index:06d}.fasta", checkpoint) != payload:
            raise ValueError("Cas context input differs from the submitted sequence")
        for key, value in _cas(root, root / "casandra" / f"record-{index:06d}", payload, source, sequence, checkpoint).items():
            counts[key] += value
    if (stage != {"status": "completed", **counts}
            or any(type(stage[key]) is not int for key in counts)):
        raise ValueError("Cas context counts disagree")
