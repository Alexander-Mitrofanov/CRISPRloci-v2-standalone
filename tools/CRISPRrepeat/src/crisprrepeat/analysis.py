"""Repeat evidence composition with explicit orientation alternatives."""

from __future__ import annotations

from . import __version__
from .folding import Folder, Settings, comparison
from .inputs import load_arrays, read_fasta, reverse_complement, sha
from .references import References

LIMITATIONS = [
    "Pair support is a descriptive thermodynamic model comparison, not processing or immune activity.",
    "RNAplfold pair values use different windows; their sums are not partner probabilities.",
    "Array boundaries are analysis boundaries; transcription ends and protein binding are not modeled.",
    "Isolated MFE motifs are hypotheses; non-hairpin and protein-stabilized repeats remain valid.",
    "Different instance motifs are not automatically homologous; do not pool them as one shared motif.",
]


def analyze(
    fasta,
    *,
    reports=None,
    molecule="DNA",
    orientation="both",
    settings=Settings(),
    categories=None,
    reference=None,
):
    kind = "arrays" if reports is not None else "repeats"
    if kind == "arrays" and molecule != "DNA":
        raise ValueError("Genomic array context requires DNA")
    if orientation not in {"both", "forward", "reverse"}:
        raise ValueError("Invalid orientation")
    if molecule == "RNA" and orientation != "forward":
        raise ValueError("RNA is supplied in transcribed 5-prime to 3-prime orientation")
    categories = categories or ["Bona-fide", "Possible"]
    sources = read_fasta(fasta, molecule, repeats=kind == "repeats")
    folder = Folder(settings)
    references = References(reference) if reference else None
    arrays = load_arrays(reports, sources, categories) if reports else []
    orientations = ["forward", "reverse"] if orientation == "both" else [orientation]
    result = {
        "schema": {"name": "CRISPRrepeat-report", "version": "1.0.0"},
        "tool": {"name": "CRISPRrepeat", "version": __version__},
        "status": "completed",
        "input_kind": kind,
        "molecule": molecule,
        "orientation_policy": orientation,
        "coordinate_system": "zero-based half-open, source forward",
        "normalization": "ASCII uppercase; DNA T to RNA U for folding; original molecular identity retained",
        "model": folder.provenance,
        "sources": [{"id": k, "length": len(v), "sha256": sha(v)} for k, v in sources.items()],
        "categories": categories if kind == "arrays" else [],
        "arrays": arrays,
        "contexts": [],
        "consensus_folds": [],
        "instances": [],
        "limitations": LIMITATIONS,
        "classification": {
            "strand": {"status": "unsupported", "reason": "no_calibrated_strand_model"},
            "subtype": {"status": "unsupported", "reason": "no_validated_subtype_model"},
            "family": {"status": "unsupported", "reason": "no_validated_family_model"},
        },
    }

    def add_instance(source_id, sequence, strand, *, array=None, repeat=None, context=None, lookup=None):
        oriented = reverse_complement(sequence) if strand == "reverse" else sequence
        rna = oriented.replace("T", "U")
        isolated = folder.isolated(rna)
        offset = None
        if array:
            offset = repeat["start"] - array["start"] if strand == "forward" else array["end"] - repeat["end"]
        result["instances"].append(
            {
                "id": f"repeat_{len(result['instances']) + 1:06d}",
                "source_id": source_id,
                "array_id": array["id"] if array else None,
                "repeat_index": repeat["index"] if repeat else None,
                "start": repeat["start"] if repeat else None,
                "end": repeat["end"] if repeat else None,
                "terminal": bool(array and repeat["index"] in {1, len(array["repeats"])}),
                "orientation": strand,
                "source_sequence": sequence,
                "rna_sequence": rna,
                "rna_sha256": sha(rna),
                "context_id": context["id"] if context else None,
                "context_start": offset,
                "isolated": isolated,
                "comparison": comparison(isolated, context, offset or 0, lookup),
                "references": references.search(rna)
                if references
                else {"status": "unsupported", "reason": "reference_not_configured", "hits": []},
            }
        )

    if kind == "repeats":
        for source_id, sequence in sources.items():
            for strand in orientations:
                add_instance(source_id, sequence, strand)
    for array in arrays:
        for strand in orientations:
            sequence = array["source_sequence"]
            if strand == "reverse":
                sequence = reverse_complement(sequence)
            rna = sequence.replace("T", "U")
            context = {
                "id": f"context_{len(result['contexts']) + 1:06d}",
                "array_id": array["id"],
                "orientation": strand,
                "rna_sequence": rna,
                "rna_sha256": sha(rna),
                **folder.context(rna),
            }
            result["contexts"].append(context)
            lookup = {(i, j): p for i, j, p in context.get("pairs", [])}
            consensus = array["consensus_repeat"]
            if strand == "reverse":
                consensus = reverse_complement(consensus)
            consensus = consensus.replace("T", "U")
            result["consensus_folds"].append(
                {
                    "array_id": array["id"],
                    "orientation": strand,
                    "rna_sequence": consensus,
                    "isolated": folder.isolated(consensus),
                }
            )
            for repeat in array["repeats"]:
                add_instance(
                    array["source_id"],
                    repeat["sequence"],
                    strand,
                    array=array,
                    repeat=repeat,
                    context=context,
                    lookup=lookup,
                )
    result["counts"] = {
        "sources": len(sources),
        "arrays": len(arrays),
        "instances": len(result["instances"]),
        "contexts": len(result["contexts"]),
    }
    return result
