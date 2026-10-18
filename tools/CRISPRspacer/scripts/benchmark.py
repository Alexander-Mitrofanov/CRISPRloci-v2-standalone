"""Seeded synthetic benchmark; no biological performance claims."""

from pathlib import Path
import json
import platform
import random
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crisprspacer import __version__  # noqa: E402
from crisprspacer.inputs import from_spacers  # noqa: E402
from crisprspacer.io import digest  # noqa: E402
from crisprspacer.search import Policy, analyze  # noqa: E402


def record(seq, name, topology="linear"):
    return {"id": name, "sequence": seq, "length": len(seq), "sha256": digest(seq),
            "topology": topology, "assembly_id": None}


def main():
    rng = random.Random(20260928)
    text = "".join(rng.choice("ACGT") for _ in range(100_000))
    queries = [record(text[i:i+32], f"spacer-{i}") for i in range(0, 100_000, 1000)]
    inventory = from_spacers(queries, digest(b"seed=20260928;synthetic-benchmark-v1"))
    rows = []
    for topology in ("linear", "circular"):
        for substitutions in (0, 1, 2):
            start = time.perf_counter()
            report = analyze(inventory, [record(text, "synthetic", topology)], Policy(substitutions=substitutions))
            rows.append({"topology": topology, "substitutions": substitutions,
                         "seconds": round(time.perf_counter()-start, 4), "counts": report["counts"]})
    print(json.dumps({"benchmark": "synthetic-v1", "tool_version": __version__,
                      "python": platform.python_version(), "target_bases": len(text),
                      "query_occurrences": len(queries), "query_length": 32, "seed": 20260928, "runs": rows}, indent=2))


if __name__ == "__main__":
    main()
