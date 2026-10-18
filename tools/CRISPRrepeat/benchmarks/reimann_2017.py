"""Run two attributed experimental substrates without deriving cleavage predictions."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time

from crisprrepeat import __version__
from crisprrepeat.folding import Folder, Settings, comparison


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = ROOT / "tests" / "fixtures" / "reimann_2017.json"
LIMITATIONS = [
    "Two related substrates from one study are a reproduction exercise, not held-out biological validation.",
    "Current pair sets come from isolated MFE hypotheses, not a recovered curated motif from the paper.",
    "ViennaRNA 2.7.0 and these scores differ from the published model; exact published scores are not reproduced.",
    "Experimental cleavage labels and published model scores are retained separately from current predictions.",
    "No cleavage threshold, classification or universal accuracy estimate is produced.",
    "Changing the window can materially change support even for the same sequence and maximum pair span.",
    "Elapsed time describes this small run and is not a scalability benchmark.",
]


def digest(value):
    return hashlib.sha256(value).hexdigest()


def read_fixture(path):
    raw = Path(path).read_bytes()
    fixture = json.loads(raw)
    for case in fixture["cases"]:
        sequence = case["rna_sequence"]
        if (len(sequence) != case["length"] or set(sequence) - set("ACGU")
                or digest(sequence.encode("ascii")) != case["rna_sha256"]
                or not sequence.startswith(case["synthetic_5prime_prefix"])):
            raise ValueError("Published transcript sequence/checksum mismatch")
        for repeat in case["repeats"]:
            a, b = repeat["start"], repeat["end"]
            if not 0 <= a < b <= len(sequence) or sequence[a:b] != fixture["repeat_sequence"]:
                raise ValueError("Published repeat interval does not reconstruct")
    return fixture, digest(raw)


def run_benchmark(fixture_path=DEFAULT_FIXTURE):
    started = time.perf_counter()
    fixture, fixture_sha = read_fixture(fixture_path)
    cases = []
    for case in fixture["cases"]:
        # The third preset changes only W relative to the second; keep that
        # contrast separate from the default's different lonely-pair policy.
        presets = [
            ("default_local", Settings()),
            ("local_no_lonely_pairs", Settings(no_lonely_pairs=True)),
            ("whole_fragment_no_lonely_pairs", Settings(window=case["length"], no_lonely_pairs=True)),
        ]
        runs = []
        for name, settings in presets:
            folder = Folder(settings)
            isolated = folder.isolated(fixture["repeat_sequence"])
            context = folder.context(case["rna_sequence"])
            repeats = []
            for repeat in case["repeats"]:
                repeats.append({
                    "id": repeat["id"],
                    "start": repeat["start"],
                    "end": repeat["end"],
                    "experimental_cleavage": repeat["experimental_cleavage"],
                    "published_model_score": repeat["published_model_score"],
                    "computed_comparison": comparison(isolated, context, repeat["start"]),
                })
            runs.append({
                "preset": name,
                "model": folder.provenance,
                "effective_context": {key: context[key] for key in (
                    "status", "method", "window", "span", "requested_window", "requested_span")},
                "isolated_mfe_structure": isolated["structure"],
                "isolated_mfe_kcal_mol": isolated["mfe_kcal_mol"],
                "repeats": repeats,
            })
        cases.append({
            "id": case["id"], "length": case["length"], "rna_sequence": case["rna_sequence"],
            "rna_sha256": case["rna_sha256"], "purpose": case["purpose"], "runs": runs,
        })
    return {
        "schema": {"name": "CRISPRrepeat-published-benchmark", "version": "1.0.0"},
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool": {"name": "CRISPRrepeat", "version": __version__},
        "runtime": {"python": platform.python_version(), "platform": platform.platform()},
        "fixture_sha256": fixture_sha,
        "script_sha256": digest(Path(__file__).read_bytes()),
        "source": fixture["source"],
        "molecule": fixture["molecule"],
        "orientation": fixture["orientation"],
        "coordinate_system": fixture["coordinate_system"],
        "repeat_sequence": fixture["repeat_sequence"],
        "cases": cases,
        "limitations": LIMITATIONS,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--output", type=Path, help="Write JSON here; omit to print JSON to stdout")
    args = parser.parse_args()
    output = json.dumps(run_benchmark(args.fixture), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
