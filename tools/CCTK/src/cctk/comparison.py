"""Observed spacer sharing: occurrence order and set membership remain distinct."""

from __future__ import annotations
from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any

from . import __version__
from .errors import CCTKError
from .io import (
    MAX_ARRAYS, MAX_OCCURRENCES, MAX_SPACERS_PER_ARRAY, write_json, safe_open,
    publish_directory,
)
from .sequence_operations import build_network, rev_comp

DNA = frozenset("ACGT")


@dataclass(frozen=True)
class ObservedArray:
    id: str
    source_id: str
    array_id: str
    strand: str | None
    occurrences: tuple[dict[str, Any], ...]
    spacers: tuple[str, ...]


def identifier(value: Any) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 256
        or any(ord(c) < 32 or ord(c) == 127 or 0xD800 <= ord(c) <= 0xDFFF for c in value)
    ):
        raise CCTKError("invalid identifier")
    return value


def parse_document(document: Any, identity: str = "exact") -> list[ObservedArray]:
    if not isinstance(identity, str) or identity not in {"exact", "reverse-complement"}:
        raise CCTKError("unknown sequence identity policy")
    if (
        not isinstance(document, dict)
        or set(document) != {"schema", "arrays"}
        or document["schema"] != {"name": "cctk.arrays", "version": "1.0.0"}
    ):
        raise CCTKError("expected cctk.arrays 1.0.0")
    raw_arrays = document["arrays"]
    if not isinstance(raw_arrays, list) or not 1 <= len(raw_arrays) <= MAX_ARRAYS:
        raise CCTKError("unsupported array count")
    ids: set[str] = set()
    scoped: set[tuple[str, str]] = set()
    output = []
    total = 0
    symbols_by_sequence: dict[str, str] = {}
    for raw in raw_arrays:
        if not isinstance(raw, dict) or set(raw) != {
            "id",
            "source_id",
            "array_id",
            "strand",
            "spacers",
        }:
            raise CCTKError("invalid array fields")
        aid, source, original = (identifier(raw[key]) for key in ("id", "source_id", "array_id"))
        if aid in ids or (source, original) in scoped or raw["strand"] not in (None, "+", "-"):
            raise CCTKError("duplicate array identity or unsupported strand")
        ids.add(aid)
        scoped.add((source, original))
        spacers = raw["spacers"]
        if not isinstance(spacers, list) or not 1 <= len(spacers) <= MAX_SPACERS_PER_ARRAY:
            raise CCTKError("unsupported spacer count")
        total += len(spacers)
        if total > MAX_OCCURRENCES:
            raise CCTKError("too many spacer occurrences")
        occurrences, symbols = [], []
        seen: set[str] = set()
        last = 0
        for spacer in spacers:
            if not isinstance(spacer, dict) or set(spacer) != {
                "id",
                "ordinal",
                "sequence",
                "status",
            }:
                raise CCTKError("invalid spacer fields")
            ident = identifier(spacer["id"])
            ordinal = spacer["ordinal"]
            if ident in seen or type(ordinal) is not int or ordinal != last + 1:
                raise CCTKError("duplicate occurrence or noncontiguous ordinal")
            seen.add(ident)
            last = ordinal
            sequence = spacer["sequence"]
            symbol = None
            if spacer["status"] == "deletion" and sequence is None:
                pass
            elif (
                spacer["status"] == "observed"
                and isinstance(sequence, str)
                and 1 <= len(sequence) <= 1000
                and set(sequence) <= DNA
            ):
                if sequence not in symbols_by_sequence:
                    canonical = (
                        min(sequence, rev_comp(sequence))
                        if identity == "reverse-complement" else sequence
                    )
                    symbols_by_sequence[sequence] = (
                        "s_" + hashlib.sha256(canonical.encode()).hexdigest()
                    )
                symbol = symbols_by_sequence[sequence]
                symbols.append(symbol)
            else:
                raise CCTKError(
                    "spacers require observed A/C/G/T DNA or an explicit deletion; ambiguity is unsupported"
                )
            occurrences.append({**spacer, "symbol": symbol})
        output.append(
            ObservedArray(aid, source, original, raw["strand"], tuple(occurrences), tuple(symbols))
        )
    return output


def components(ids: list[str], edges: list[dict[str, Any]]) -> list[list[str]]:
    neighbours: dict[str, set[str]] = {ident: set() for ident in ids}
    for edge in edges:
        neighbours[edge["a"]].add(edge["b"])
        neighbours[edge["b"]].add(edge["a"])
    visited: set[str] = set()
    result = []
    for ident in ids:
        if ident in visited:
            continue
        stack = [ident]
        group: set[str] = set()
        while stack:
            node = stack.pop()
            if node in group:
                continue
            group.add(node)
            visited.add(node)
            stack.extend(neighbours[node] - group)
        result.append([node for node in ids if node in group])
    return result


def compare(
    document: Any, *, identity: str = "exact", min_shared: int = 1, order: str = "supplied"
) -> dict[str, Any]:
    if (
        type(min_shared) is not int
        or not 1 <= min_shared <= MAX_SPACERS_PER_ARRAY
        or not isinstance(order, str)
        or order not in {"supplied", "id"}
    ):
        raise CCTKError("invalid comparison options")
    arrays = parse_document(document, identity)
    if order == "id":
        arrays.sort(key=lambda array: array.id)
    by_id = {array.id: array for array in arrays}
    edges = []
    for edge in build_network(arrays):
        if edge.nshared >= min_shared:
            edges.append(
                {
                    "a": edge.a,
                    "b": edge.b,
                    "nshared": edge.nshared,
                    "jaccard": edge.jaccard,
                    "shared_symbols": sorted(edge.shared_spacers),
                }
            )
    groups = components(list(by_id), edges)
    return {
        "schema": {"name": "cctk.comparison", "version": "1.0.0"},
        "tool": {"name": "CCTK", "version": __version__},
        "options": {"identity": identity, "min_shared": min_shared, "order": order},
        "arrays": [
            {
                "id": a.id,
                "source_id": a.source_id,
                "array_id": a.array_id,
                "strand": a.strand,
                "spacers": list(a.occurrences),
                "observed_spacer_count": len(a.spacers),
                "evidence_status": "observed" if a.spacers else "no_observed_spacers",
            }
            for a in arrays
        ],
        "edges": edges,
        "components": groups,
        "isolates": [group[0] for group in groups if len(group) == 1],
        "interpretation": "Observed sharing only; source-forward order is not inferred acquisition order or a phylogeny.",
    }


def run_compare(
    input_path: Path,
    outdir: Path,
    *,
    identity: str,
    min_shared: int,
    order: str,
    plot: bool = False,
    seed: int = 2,
) -> dict[str, Any]:
    import json
    import shutil
    import tempfile
    from .io import reject_duplicate_keys

    with safe_open(input_path, "rb") as handle:
        snapshot = handle.read()
    digest = hashlib.sha256(snapshot).hexdigest()
    try:
        document = json.loads(
            snapshot,
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=lambda _: (_ for _ in ()).throw(CCTKError("nonfinite JSON")),
        )
    except (json.JSONDecodeError, UnicodeError, RecursionError) as error:
        raise CCTKError("invalid JSON input") from error
    result = compare(document, identity=identity, min_shared=min_shared, order=order)
    result["input_sha256"] = digest
    if outdir.exists() or outdir.is_symlink() or any(p.is_symlink() for p in outdir.parents):
        raise CCTKError("output directory must be new and nonsymlinked")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".cctk-", dir=outdir.parent))
    try:
        write_json(staging / "comparison.json", result)
        with safe_open(staging / "network.tsv", "w") as handle:
            handle.write("array_a\tarray_b\tshared_spacers\tjaccard\n")
            for edge in result["edges"]:
                handle.write(
                    f"{table_cell(edge['a'])}\t{table_cell(edge['b'])}\t{edge['nshared']}\t{edge['jaccard']:.17g}\n"
                )
        if plot:
            plot_comparison(result, staging / "comparison.png", seed)
        publish_directory(staging, outdir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return result


def plot_comparison(result: dict[str, Any], path: Path, seed: int = 2) -> None:
    """Use the native CRISPRdiff geometry; deletion slots are separate gap symbols."""
    from . import plotting, colour_schemes

    arrays = {
        a["id"]: [
            s["symbol"] if s["symbol"] is not None else f"__gap_{index}_{s['ordinal']}"
            for s in a["spacers"]
        ]
        for index, a in enumerate(result["arrays"])
    }
    gap_symbols = {
        symbol for symbols in arrays.values() for symbol in symbols if symbol.startswith("__gap_")
    }
    all_symbols: dict[str, int] = {}
    for symbols in arrays.values():
        for symbol in dict.fromkeys(symbols):
            all_symbols[symbol] = all_symbols.get(symbol, 0) + 1
    shared = [symbol for symbol, count in all_symbols.items() if count > 1]
    colours = dict(zip(shared, colour_schemes.choose_col_scheme(len(shared), seed)))
    plotting.plot_diffplot(
        arrays,
        list(arrays)[::-1],
        shared,
        colours,
        plot_width=8,
        plot_height=min(16, max(3, len(arrays) * 0.3)),
        dpi=120,
        outfile=str(path),
        gap_symbols=gap_symbols,
        observed_view=True,
    )


def table_cell(value: str) -> str:
    return ("'" + value) if value.lstrip().startswith(("=", "+", "-", "@")) else value
