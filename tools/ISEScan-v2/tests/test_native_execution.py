from pathlib import Path
import sys
from unittest.mock import patch
import pytest

NATIVE = Path(__file__).resolve().parents[1] / "isescan_v2/native"
sys.path.insert(0, str(NATIVE))
import is_analysis  # noqa: E402
import tools as native_tools  # noqa: E402


def test_fraggenescan_uses_direct_argv_and_preserves_shell_characters():
    with (
        patch.object(is_analysis.shutil, "which", return_value="/safe/FragGeneScan"),
        patch.object(is_analysis.subprocess, "call", return_value=0) as call,
    ):
        assert (
            is_analysis.translate_genome_dna_v3(
                ("/a space/x;echo.fna", "/out space", "0", "illumina_5", 1)
            )
            == 0
        )
        assert call.call_args.args[0] == [
            "/safe/FragGeneScan",
            "-s",
            "/a space/x;echo.fna",
            "-o",
            "/out space",
            "-w",
            "0",
            "-t",
            "illumina_5",
            "-p",
            "1",
        ]
        assert call.call_args.kwargs["shell"] is False


def test_blast_failure_cannot_be_interpreted_as_no_hit():
    with patch.object(native_tools.subprocess, "Popen") as popen:
        popen.return_value.communicate.return_value = ("", "failed")
        popen.return_value.returncode = 2
        with pytest.raises(RuntimeError, match="BLAST execution failed"):
            native_tools.doBlastnOnStream(">a\nAAAA\n", "/db space")
        assert "/db space" in popen.call_args.args[0]


def test_hmmer_paths_preserved():
    with patch.object(is_analysis.subprocess, "call", return_value=0) as call:
        is_analysis.is_hmmsearch_v2(("/m space", "/d space", "/o space", 1))
        assert call.call_args.args[0][-2:] == ["/m space", "/d space"]
