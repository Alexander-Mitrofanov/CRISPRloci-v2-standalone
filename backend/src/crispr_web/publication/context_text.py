"""Atomic, bounded text exports without following tool-created destinations."""

from contextlib import contextmanager
import os
from pathlib import Path


class _BoundedText:
    def __init__(self, stream, encoding, limit, checkpoint):
        self.stream = stream
        self.encoding = encoding
        self.limit = limit
        self.checkpoint = checkpoint
        self.size = 0

    def write(self, value):
        if self.checkpoint:
            self.checkpoint()
        self.size += len(value.encode(self.encoding))
        if self.size > self.limit:
            raise ValueError("Context text export exceeds byte limit")
        return self.stream.write(value)


@contextmanager
def atomic_text(root, name, *, encoding="utf-8", max_bytes=20_000_000, checkpoint=None):
    """Replace a fixed export name atomically; never open its existing inode."""
    if Path(name).name != name or name in {"", ".", ".."}:
        raise ValueError("Context export must have a simple filename")
    temporary = root / f".{name}.tmp"
    created = False
    try:
        with temporary.open("x", encoding=encoding, newline="") as stream:
            created = True
            yield _BoundedText(stream, encoding, max_bytes, checkpoint)
            stream.flush()
            os.fsync(stream.fileno())
            os.fchmod(stream.fileno(), 0o600)
        os.replace(temporary, root / name)
    finally:
        if created:
            temporary.unlink(missing_ok=True)
