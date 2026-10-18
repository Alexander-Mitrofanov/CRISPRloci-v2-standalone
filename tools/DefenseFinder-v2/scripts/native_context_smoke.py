"""Run the real pinned HMM workflow and verify a known upstream fragment."""
import argparse
from pathlib import Path
import tempfile

from defense_context.pipeline import run

parser = argparse.ArgumentParser()
parser.add_argument("--models-dir", required=True, type=Path)
parser.add_argument("--hmmer", required=True, type=Path)
args = parser.parse_args()
with tempfile.TemporaryDirectory() as temporary:
    report = run(Path(__file__).resolve().parents[1] / "tests/data/context_fragment.fna",
                 Path(temporary) / "result", args.models_dir.absolute(), args.hmmer.absolute(), 1)
    assert report["counts"] == {"sources": 1, "genes": 42, "systems": 1, "components": 2, "hits": 5}
    assert report["systems"][0]["subtype"] == "RM_Type_II"
    print("Native DefenseFinder fixture passed: 42 proteins, 1 system, 2 components, 5 profile hits")
