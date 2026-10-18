"""Annotation outputs fail closed before any scientific interpretation."""

import hashlib

import pytest

from crispr_integration import annotations


@pytest.mark.parametrize("payload", [
    '{"status":"failed","status":"complete"}',
    '{"nested":{"score":0,"score":1}}',
    '{"score":NaN}', '{"score":Infinity}', '{"score":1e999}',
])
@pytest.mark.parametrize("reader", [annotations.read_object, annotations.json_lines])
def test_annotation_json_rejects_ambiguous_and_nonfinite_evidence(tmp_path, payload, reader):
    path = tmp_path / "report.json"
    path.write_text(payload)
    with pytest.raises(ValueError):
        reader(tmp_path, path)


@pytest.mark.parametrize("reader", [annotations.read_object, annotations.json_lines])
def test_annotation_reader_rejects_parent_symlink(tmp_path, reader):
    source = tmp_path / "source"
    source.mkdir()
    (source / "report.json").write_text('{"score":null}')
    (tmp_path / "link").symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic links"):
        reader(tmp_path, tmp_path / "link/report.json")


def test_annotation_checksum_is_bounded(tmp_path, monkeypatch):
    path = tmp_path / "report.json"
    path.write_bytes(b"12345")
    assert annotations.checksum(path) == hashlib.sha256(b"12345").hexdigest()
    monkeypatch.setattr(annotations, "MAX_REPORT_BYTES", 4)
    with pytest.raises(ValueError, match="oversized"):
        annotations.checksum(path)


def test_annotation_readers_preserve_missing_evidence(tmp_path):
    path = tmp_path / "report.json"
    path.write_text('{"score":null,"other":1.5}')
    assert annotations.read_object(tmp_path, path) == {"score": None, "other": 1.5}
    assert annotations.json_lines(tmp_path, path) == [{"score": None, "other": 1.5}]
