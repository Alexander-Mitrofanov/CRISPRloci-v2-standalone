"""Convert CRISPRidentify reports into the shared CCDB directory contract."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Iterable

from .parser import ArrayRecord, load_identified_arrays, normalize_category


_COMPLEMENT = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
})


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical_dna(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def _safe_name(value: str, *, fallback: str, max_length: int = 80) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return (name or fallback)[:max_length]


def _portable_report_path(report_path: Path, input_root: Path) -> str:
    """Return a stable report identifier without exposing a host filesystem path."""

    try:
        return report_path.resolve().relative_to(input_root).as_posix()
    except ValueError:
        # This is a defensive fallback: discovered reports normally live under
        # ``input_root``. Keep provenance useful without publishing a path
        # outside the submitted result tree if that invariant ever changes.
        basename = _safe_name(report_path.name, fallback="report.json", max_length=120)
        return f"external/{basename}"


def _array_token(array: ArrayRecord, *, input_root: Path | None = None) -> str:
    report_reference = (
        _portable_report_path(array.report_path, input_root)
        if input_root is not None
        else _safe_name(array.report_path.name, fallback="report.json", max_length=120)
    )
    identity = "\0".join((
        report_reference, array.source_id, array.array_id,
        array.consensus_repeat, *array.spacers,
    )).encode("utf-8")
    digest = hashlib.sha256(identity).hexdigest()[:10]
    source = _safe_name(array.source_id, fallback="source", max_length=36)
    array_id = _safe_name(array.array_id, fallback="array", max_length=36)
    return f"{source}__{array_id}__{digest}"


class _DisjointSet:
    def __init__(self, size: int):
        self.parents = list(range(size))

    def find(self, item: int) -> int:
        while self.parents[item] != item:
            self.parents[item] = self.parents[self.parents[item]]
            item = self.parents[item]
        return item

    def union(self, first: int, second: int) -> None:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root != second_root:
            self.parents[second_root] = first_root


def _overlap_components(
    arrays: list[ArrayRecord], *, input_root: Path | None = None
) -> list[list[ArrayRecord]]:
    """Find transitive array components sharing orientation-neutral spacers."""

    disjoint = _DisjointSet(len(arrays))
    first_by_spacer: dict[str, int] = {}
    for array_index, array in enumerate(arrays):
        for spacer in {canonical_dna(sequence) for sequence in array.spacers}:
            previous = first_by_spacer.setdefault(spacer, array_index)
            disjoint.union(previous, array_index)
    grouped: dict[int, list[ArrayRecord]] = defaultdict(list)
    for index, array in enumerate(arrays):
        grouped[disjoint.find(index)].append(array)
    return sorted(
        (
            sorted(component, key=lambda item: _array_token(item, input_root=input_root))
            for component in grouped.values()
        ),
        key=lambda component: _array_token(component[0], input_root=input_root),
    )


def group_arrays(
    arrays: Iterable[ArrayRecord], *, input_root: Path | None = None
) -> list[tuple[str, list[ArrayRecord]]]:
    """Group arrays by orientation-neutral repeat and transitive spacer overlap."""

    by_repeat: dict[str, list[ArrayRecord]] = defaultdict(list)
    for array in arrays:
        by_repeat[canonical_dna(array.consensus_repeat)].append(array)

    groups: list[tuple[str, list[ArrayRecord]]] = []
    for repeat in sorted(by_repeat):
        digest = hashlib.sha256(repeat.encode("ascii")).hexdigest()[:12]
        components = _overlap_components(by_repeat[repeat], input_root=input_root)
        for component_index, component in enumerate(components, 1):
            groups.append((f"repeat_{digest}_group_{component_index:03d}", component))
    return groups


def _atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    os.replace(temporary, path)


def _ensure_empty_output(path: Path) -> None:
    if path.exists() and not path.is_dir():
        raise FileExistsError(f"adapter output exists and is not a directory: {path}")
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"adapter output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _manifest_array(
    array: ArrayRecord, *, input_root: Path, reason: str | None = None
) -> dict[str, object]:
    record: dict[str, object] = {
        "array_id": array.array_id,
        "category": array.category,
        "consensus_repeat": array.consensus_repeat,
        "input_sequence_orientation": array.sequence_orientation,
        "omitted_empty_spacers": array.omitted_empty_spacers,
        "report_path": _portable_report_path(array.report_path, input_root),
        "source_id": array.source_id,
        "spacer_count": len(array.spacers),
        "strand": array.strand,
    }
    if reason is not None:
        record["reason"] = reason
    return record


def adapt_reports(
    report_path: str | Path,
    output_path: str | Path,
    *,
    categories: Iterable[str] = ("Bona-fide", "Possible"),
    minimum_arrays: int = 2,
    expected_fasta_path: str | Path | None = None,
) -> dict[str, object]:
    """Write a deterministic CCDB tree and return its machine manifest."""

    if minimum_arrays < 2:
        raise ValueError("minimum_arrays must be at least 2 for evolutionary reconstruction")
    accepted_categories = {normalize_category(value) for value in categories}
    if not accepted_categories:
        raise ValueError("at least one CRISPRidentify category must be accepted")

    source = Path(report_path).expanduser().resolve()
    input_is_directory = source.is_dir()
    input_root = source if input_is_directory else source.parent
    input_reference = "." if input_is_directory else source.name

    destination = Path(output_path).expanduser().resolve()
    _ensure_empty_output(destination)
    ccdb_path = destination / "ccdb"
    ccdb_path.mkdir()

    all_arrays = load_identified_arrays(
        report_path,
        expected_fasta_path=expected_fasta_path,
    )
    selected: list[ArrayRecord] = []
    skipped: list[dict[str, object]] = []
    for array in all_arrays:
        if array.category not in accepted_categories:
            skipped.append(_manifest_array(
                array, input_root=input_root, reason="category_not_selected"
            ))
        elif not array.spacers:
            skipped.append(_manifest_array(
                array, input_root=input_root, reason="no_non_empty_spacers"
            ))
        elif array.sequence_orientation == "source" and array.strand is None:
            # CRISPR-evOr compares one orientation for the complete group with
            # its global reversal. It cannot repair one independently
            # misoriented member, so assuming forward here can corrupt both the
            # phylogeny and the likelihood comparison. Explicitly biologically
            # oriented sequences remain safe without a strand label because
            # their sequence frame has already been resolved by the producer.
            skipped.append(_manifest_array(
                array, input_root=input_root, reason="orientation_not_determined"
            ))
        else:
            selected.append(array)

    emitted_groups: list[dict[str, object]] = []
    emitted_array_keys: set[tuple[Path, str, str]] = set()
    for group_name, arrays in group_arrays(selected, input_root=input_root):
        if len(arrays) < minimum_arrays:
            skipped.extend(
                _manifest_array(
                    array, input_root=input_root, reason="undersized_repeat_overlap_group"
                )
                for array in arrays
            )
            continue

        group_path = ccdb_path / group_name
        group_manifest: dict[str, object] = {
            "name": group_name,
            "array_count": len(arrays),
            "repeat_key": canonical_dna(arrays[0].consensus_repeat),
            "arrays": [],
        }
        for array in arrays:
            if array.sequence_orientation == "biological":
                output_strand = "+"
                output_spacers = array.spacers
            else:
                output_strand = "-" if array.strand == "-" else "+"
                output_spacers = array.spacers
            strand_folder = group_path / ("neg_strand" if output_strand == "-" else "pos_strand")
            strand_folder.mkdir(parents=True, exist_ok=True)
            token = _array_token(array, input_root=input_root)
            array_path = strand_folder / f"{token}.fasta"
            with array_path.open("x", encoding="utf-8") as stream:
                for spacer_index, spacer in enumerate(output_spacers, 1):
                    stream.write(f">{token}_spacer_{spacer_index}\n{spacer}\n")
            array_entry = _manifest_array(array, input_root=input_root)
            array_entry.update({
                "ccdb_file": str(array_path.relative_to(destination)),
                "ccdb_strand": output_strand,
            })
            group_manifest["arrays"].append(array_entry)
            emitted_array_keys.add((array.report_path, array.source_id, array.array_id))
        emitted_groups.append(group_manifest)

    manifest: dict[str, object] = {
        "schema": {"name": "crispr-tool-integration-adapter", "version": "1.1.0"},
        "input": input_reference,
        "input_kind": "directory" if input_is_directory else "file",
        "path_contract": {
            "base": "adapter_input_root",
            "input_root": ".",
            "report_paths": (
                "POSIX paths relative to the adapter input root; "
                "external/<safe-basename> is used defensively for out-of-root reports"
            ),
        },
        "output_format": "ccdb",
        "source_binding": {
            "status": "verified" if expected_fasta_path is not None else "not_requested",
            "identity": ["source.id", "source.length", "source.sha256"],
            "sha256_scope": "uppercase_sequence",
        },
        "orientation_contract": {
            "source_sequences": "stored unchanged; reverse arrays are placed in neg_strand",
            "biological_sequences": "stored in pos_strand to prevent double reversal",
            "unknown_source_orientation": "excluded; a group-wide reversal cannot repair one misoriented member",
        },
        "grouping": {
            "repeat": "exact IUPAC sequence after reverse-complement canonicalization",
            "component": "transitive overlap of reverse-complement-canonicalized spacer sequences",
            "minimum_arrays": minimum_arrays,
        },
        "accepted_categories": sorted(accepted_categories),
        "input_array_count": len(all_arrays),
        "selected_array_count": len(selected),
        "emitted_array_count": len(emitted_array_keys),
        "emitted_group_count": len(emitted_groups),
        "groups": emitted_groups,
        "skipped_arrays": sorted(
            skipped,
            key=lambda item: (
                str(item["report_path"]), str(item["source_id"]), str(item["array_id"]), str(item["reason"])
            ),
        ),
    }
    _atomic_json(destination / "manifest.json", manifest)
    return manifest
