"""Kernel publication must preserve concurrent outputs and their identity."""
import pytest

from casandra.utils import _rename_new_directory


@pytest.mark.parametrize("kind", ["directory", "file", "symlink"])
def test_publication_never_replaces_a_concurrent_output(tmp_path, kind):
    staging, target = tmp_path / "staging", tmp_path / "target"
    staging.mkdir()
    (staging / "result.json").write_text("complete")
    if kind == "directory":
        target.mkdir()
    elif kind == "file":
        target.write_text("preserve")
    else:
        target.symlink_to(tmp_path / "unrelated")
    identity = target.lstat().st_ino
    with pytest.raises((FileExistsError, OSError)):
        _rename_new_directory(staging, target)
    assert target.lstat().st_ino == identity
    assert (staging / "result.json").read_text() == "complete"
    if kind == "directory":
        assert list(target.iterdir()) == []
    elif kind == "file":
        assert target.read_text() == "preserve"
    else:
        assert target.is_symlink()


def test_publication_moves_a_complete_new_directory(tmp_path):
    staging, target = tmp_path / "staging", tmp_path / "target"
    staging.mkdir()
    (staging / "result.json").write_text("complete")
    _rename_new_directory(staging, target)
    assert not staging.exists()
    assert (target / "result.json").read_text() == "complete"
