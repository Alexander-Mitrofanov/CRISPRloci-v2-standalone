"""Frozen pre-refinement native bytes for indexed, interleaved alignment rows."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

NATIVE = os.environ.get("SPACEPHARER_NATIVE")
FROZEN = json.loads((Path(__file__).parent / "fixtures/indexed-records.json").read_text())
pytestmark = pytest.mark.skipif(not NATIVE, reason="candidate native executable required")


def database(path: Path, text: str, kind: int) -> str:
    data = text.encode() + b"\0"
    path.write_bytes(data)
    path.with_suffix(".index").write_text(f"0\t0\t{len(data)}\n")
    path.with_suffix(".dbtype").write_bytes(kind.to_bytes(4, "little"))
    return str(path)


@pytest.mark.parametrize("threads", [1, 2])
def test_native_interleaved_combination_bytes(tmp_path, threads):
    """Repeated target keys retain the final nucleotide row and protein order."""
    output = tmp_path / "combined"
    subprocess.run(
        [
            NATIVE,
            "combineprotnuclaln",
            database(tmp_path / "protein", FROZEN["protein"], 5),
            database(tmp_path / "nucleotide", FROZEN["nucleotide"], 5),
            str(output),
            "--threads",
            str(threads),
            "-v",
            "0",
        ],
        check=True,
        timeout=30,
    )
    expected = FROZEN["combined"].encode()
    # MMseqs preserves multiple unmerged thread shards for result databases.
    data = (
        output.read_bytes()
        if output.is_file()
        else b"".join(part.read_bytes() for part in sorted(output.parent.glob(output.name + ".[0-9]")))
    )
    assert data == expected
    assert output.with_suffix(".index").read_bytes() == f"0\t0\t{len(expected)}\n".encode()
    assert output.with_suffix(".dbtype").read_bytes() == (5).to_bytes(4, "little")


@pytest.mark.parametrize("fmt", [0, 1, 2])
@pytest.mark.parametrize("report_fdr", [0, 1])
@pytest.mark.parametrize("threads", [1, 2])
@pytest.mark.parametrize("missing_position", [0, 1, 2])
def test_native_interleaved_summary_bytes(tmp_path, fmt, report_fdr, threads, missing_position):
    """Keep target/support order, gaps, descending coordinates and PAM columns."""
    output = tmp_path / "summary.tsv"
    # The absent HSP key must advance at the beginning, middle and end of a match list.
    matches = FROZEN["matches"].splitlines(keepends=True)
    missing = matches.pop()
    assert missing.split("\t", 1)[0] not in {
        line.split("\t", 1)[0] for line in FROZEN["alignments"].splitlines()
    }
    matches.insert(missing_position, missing)
    subprocess.run(
        [
            NATIVE,
            "summarizeresults",
            database(tmp_path / "target", "ACGT\n", 1),
            database(tmp_path / "matches", "".join(matches), 12),
            database(tmp_path / "alignments", FROZEN["alignments"], 12),
            str(output),
            "--fmt",
            str(fmt),
            "--report-fdr",
            str(report_fdr),
            "--threads",
            str(threads),
            "-v",
            "0",
        ],
        check=True,
        timeout=30,
    )
    assert output.read_bytes() == FROZEN["summaries"][f"{fmt}-{report_fdr}"].encode()
    assert not output.with_suffix(".tsv.index").exists()
