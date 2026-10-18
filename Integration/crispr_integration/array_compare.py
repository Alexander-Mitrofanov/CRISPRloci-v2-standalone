"""Isolated CCTK observed comparisons; no detector/evolution eligibility gate."""

from pathlib import Path
import hashlib
import json
from .protospacer import _records, _payload, _json
from .scientific_options import tool_arguments


def native_document(mapping, records):
    if mapping.get("schema") != {"name": "crisprloci.array-map", "version": "1.0.0"}:
        raise ValueError("Invalid array map")
    doc = mapping["document"]
    if (
        doc.get("schema") != {"name": "crisprloci.arrays", "version": "1.0.0"}
        or doc.get("coordinate_system") != "1-based-inclusive-source-forward"
    ):
        raise ValueError("Invalid observed arrays")
    arrays = doc["arrays"]
    output = []
    expected = {}
    scopes = set()
    if not isinstance(arrays, list) or not 1 <= len(arrays) <= 200 or sum(len(a["spacers"]) for a in arrays) > 10000:
        raise ValueError("Array count limit")
    for i, array in enumerate(arrays):
        scope = (array["source_id"], array["id"])
        if scope in scopes or array["strand"] not in (None, "+", "-"):
            raise ValueError("Invalid scoped identity/strand")
        scopes.add(scope)
        slots = []
        ids = set()
        if not 1 <= len(array["spacers"]) <= 500:
            raise ValueError("Spacer slot limit")
        for n, unit in enumerate(array["spacers"], 1):
            if unit["ordinal"] != n or type(unit["ordinal"]) is not int:
                raise ValueError("Invalid ordinal")
            if unit.get("status") == "deletion" and unit.get("sequence") is None:
                ident = f"gap{n:06d}"
                while ident in {u.get("id") for u in array["spacers"]}:
                    ident += "_"
                row = {"id": unit.get("id") or ident, "ordinal": n, "sequence": None, "status": "deletion"}
            elif (
                unit.get("status") is None
                and isinstance(unit.get("sequence"), str)
                and 1 <= len(unit["sequence"]) <= 1000
                and set(unit["sequence"]) <= set("ACGT")
            ):
                expected[f"occ{i:06d}_{n:06d}"] = unit["sequence"]
                row = {"id": unit["id"], "ordinal": n, "sequence": unit["sequence"], "status": "observed"}
            else:
                raise ValueError("Incomplete or ambiguous spacer evidence")
            if row["id"] in ids:
                raise ValueError("Duplicate occurrence")
            ids.add(row["id"])
            slots.append(row)
        output.append(
            {
                "id": f"array{i:06d}",
                "source_id": array["source_id"],
                "array_id": array["id"],
                "strand": array["strand"],
                "spacers": slots,
            }
        )
    if expected != records or not expected:
        raise ValueError("Array map differs from source DNA")
    return {"schema": {"name": "cctk.arrays", "version": "1.0.0"}, "arrays": output}


def validate_report(report, native, digest, options):
    if (
        set(report)
        != {"schema", "tool", "options", "arrays", "edges", "components", "isolates", "interpretation", "input_sha256"}
        or report.get("interpretation")
        != "Observed sharing only; source-forward order is not inferred acquisition order or a phylogeny."
    ):
        raise ValueError("CCTK schema or interpretation differs")
    if (
        report["schema"] != {"name": "cctk.comparison", "version": "1.0.0"}
        or report["tool"] not in ({"name": "CCTK", "version": "2.0.1"},
                                 {"name": "CCTK", "version": "2.0.2"})
        or report["input_sha256"] != digest
        or report["options"] != options
    ):
        raise ValueError("CCTK identity or options differ")
    arrays = native["arrays"] if options["order"] == "supplied" else sorted(native["arrays"], key=lambda a: a["id"])
    symbols = {}
    expected_arrays = []
    for a in arrays:
        slots = []
        observed = 0
        for s in a["spacers"]:
            sequence = s["sequence"]
            symbol = None
            if sequence is not None:
                if options["identity"] == "reverse-complement":
                    sequence = min(sequence, sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])
                symbol = "s_" + hashlib.sha256(sequence.encode()).hexdigest()
                observed += 1
            slots.append({**s, "symbol": symbol})
        symbols[a["id"]] = {s["symbol"] for s in slots if s["symbol"] is not None}
        expected_arrays.append(
            {
                **a,
                "spacers": slots,
                "observed_spacer_count": observed,
                "evidence_status": "observed" if observed else "no_observed_spacers",
            }
        )
    edges = []
    ids = list(symbols)
    neighbors = {i: set() for i in ids}
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            shared = symbols[a] & symbols[b]
            if len(shared) >= options["min_shared"]:
                edges.append(
                    {
                        "a": a,
                        "b": b,
                        "nshared": len(shared),
                        "jaccard": len(shared) / len(symbols[a] | symbols[b]),
                        "shared_symbols": sorted(shared),
                    }
                )
                neighbors[a].add(b)
                neighbors[b].add(a)
    groups = []
    seen = set()
    for a in ids:
        if a in seen:
            continue
        pending = [a]
        group = set()
        while pending:
            node = pending.pop()
            if node in group:
                continue
            group.add(node)
            pending.extend(neighbors[node] - group)
        seen.update(group)
        groups.append([x for x in ids if x in group])
    expected = {
        **report,
        "arrays": expected_arrays,
        "edges": edges,
        "components": groups,
        "isolates": [g[0] for g in groups if len(g) == 1],
    }
    if json.dumps(expected, sort_keys=True, allow_nan=False) != json.dumps(report, sort_keys=True, allow_nan=False):
        raise ValueError("CCTK observed counts, occurrences or components differ")
    return {
        "arrays": len(arrays),
        "spacer_occurrences": sum(a["observed_spacer_count"] for a in expected_arrays),
        "deletions": sum(len(a["spacers"]) - a["observed_spacer_count"] for a in expected_arrays),
        "edges": len(edges),
        "components": len(groups),
        "isolates": sum(len(g) == 1 for g in groups),
    }


def run_array_compare_pipeline(input_path, output_path, *, runner_config, map_path, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    records = _records(Path(input_path).absolute())
    raw = _payload(Path(map_path).parent, Path(map_path))
    mapping = _json(raw)
    native = native_document(mapping, records)
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "array-compare-input").mkdir()
    _atomic_json(job / "array-compare-input/arrays.json", native)
    _atomic_json(job / "array-compare-input/map.json", mapping)
    native_hash = hashlib.sha256((job / "array-compare-input/arrays.json").read_bytes()).hexdigest()
    supplied = runner_config.tool_options.get("cctk", {})
    options = {
        "identity": supplied.get("--identity", "exact"),
        "min_shared": supplied.get("--min-shared", 1),
        "order": supplied.get("--order", "supplied"),
    }
    result = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "workflow": "array_compare",
        "status": "running",
        "stages": {},
        "tool_options": runner_config.tool_options,
    }
    _atomic_json(job / "integration_result.json", result)
    try:
        stage = _run_process(
            "cctk",
            [
                *runner_config.cctk,
                "compare",
                "--input",
                str(job / "array-compare-input/arrays.json"),
                "--outdir",
                str(job / "array-compare"),
                "--plot",
                *tool_arguments(runner_config, "cctk"),
            ],
            job / "logs",
            timeout_seconds=timeout_seconds,
        )
        report = _json(_payload(job, job / "array-compare/comparison.json"))
        stage.update(
            counts=validate_report(report, native, native_hash, options),
            array_map_sha256=hashlib.sha256(raw).hexdigest(),
        )
        result["stages"]["cctk"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("cctk", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
