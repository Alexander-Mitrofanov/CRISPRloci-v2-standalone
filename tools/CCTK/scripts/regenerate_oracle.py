"""Regenerate scientific goldens from a hash-verified upstream checkout; never download."""

from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(__doc__)
parser.add_argument("upstream", type=Path, help="Extracted upstream repository root")
parser.add_argument("output", type=Path, help="New directory for regenerated fixtures")
args = parser.parse_args()
manifest = json.loads((root / "UPSTREAM.json").read_text())
for name, digest in manifest["files"].items():
    if name.startswith("CRISPR_comparison_toolkit/"):
        if hashlib.sha256((args.upstream / name).read_bytes()).hexdigest() != digest:
            raise SystemExit(f"Upstream source differs: {name}")
args.output.mkdir()
sys.path.insert(0, str(args.upstream / "CRISPR_comparison_toolkit"))
from cctkpkg import sequence_operations as seq, CRISPRdiff as diff, array_parsimony as pars  # noqa: E402
from cctkpkg import colour_schemes as colors, tree_operations as trees  # noqa: E402

rng = random.Random(193)
cases = []
for _ in range(50):
    a = [str(rng.randrange(1, 10)) for _ in range(rng.randrange(1, 10))]
    b = [str(rng.randrange(1, 10)) for _ in range(rng.randrange(1, 10))]
    first, second = pars.find_modules(pars.Array("a", a), pars.Array("b", b))

    def modules(arr):
        return {
            "aligned": arr.aligned,
            "modules": [
                {"type": m.type, "spacers": m.spacers, "indices": m.indices} for m in arr.modules
            ],
        }

    cases.append(
        {
            "a": a,
            "b": b,
            "alignment": seq.needle(a, b),
            "modules": [modules(first), modules(second)],
        }
    )
arrays = {
    "A": ["1", "2", "2", "3", "4"],
    "B": ["5", "2", "3", "4"],
    "C": ["6", "5"],
    "isolated": ["9", "10"],
}
out = {
    "upstream_commit": manifest["commit"],
    "pairs": cases,
    "order": diff.decide_array_order_global_best(arrays),
    "score": diff.get_list_score(arrays, list(arrays)),
    "colors": colors.choose_col_scheme(12, 2),
    "nodes": trees.create_internal_node_ids(40),
}
(args.output / "oracle.json").write_text(json.dumps(out, indent=2) + "\n")
for name in ["arrays.txt", "mapping.txt", "constrain.nwk"]:
    shutil.copyfile(root / "tests/fixtures" / name, args.output / name)
cli = args.upstream / "CRISPR_comparison_toolkit/cctk"
env = {**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(args.output / "mpl-cache")}


def run(*argv):
    result = subprocess.run(
        [sys.executable, str(cli), *map(str, argv)],
        capture_output=True,
        text=True,
        check=True,
        env=env,
        timeout=120,
    )
    (args.output / (str(argv[0]) + ".stdout")).write_text(result.stdout)
    (args.output / (str(argv[0]) + ".stderr")).write_text(result.stderr)


a = args.output / "arrays.txt"
run("network", "-i", a, "-o", args.output)
run(
    "crisprdiff",
    "-a",
    a,
    "-o",
    args.output / "diff.png",
    "--preordered",
    "--seed",
    "2",
    "--colour-scheme-outfile",
    args.output / "colours.json",
)
run(
    "crisprtree",
    "-a",
    a,
    "--output-arrays",
    args.output / "tree_arrays.txt",
    "--seed",
    "3",
    "-r",
    "2",
    "-x",
    "A",
    "B",
    "C",
)
run("evolve", "-n", "12", "-s", "3", "-o", args.output)
shutil.copyfile(next(args.output.glob("evolved_arrays*")), args.output / "evolution_arrays.txt")
shutil.copyfile(next(args.output.glob("evolved_tree_file*")), args.output / "evolution_tree.nwk")
run(
    "constrain",
    "-a",
    a,
    "-t",
    args.output / "constrain.nwk",
    "-g",
    args.output / "mapping.txt",
    "-o",
    args.output / "constrain.png",
)
receipt = {
    "upstream_commit": manifest["commit"],
    "python": platform.python_version(),
    "packages": {
        name: importlib.metadata.version(name)
        for name in ["numpy", "matplotlib", "DendroPy", "pillow"]
    },
    "fixture_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(args.output.iterdir())
        if p.is_file()
    },
}
(args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
