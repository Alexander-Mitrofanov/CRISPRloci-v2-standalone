"""Publish a complete bounded result directory via a same-filesystem rename."""
import csv
import ctypes
import errno
import io
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from .errors import OutputError


def _rename_new_directory(source: Path, destination: Path) -> None:
    """Atomically publish without replacing a concurrently created directory."""
    if os.name == "nt":
        # Windows rename already refuses to replace an existing destination.
        source.rename(destination)
        return
    library = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith("linux") and hasattr(library, "renameat2"):
        rename = library.renameat2
        rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        # AT_FDCWD, RENAME_NOREPLACE: the existence check occurs inside rename.
        result = rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    elif sys.platform == "darwin" and hasattr(library, "renamex_np"):
        rename = library.renamex_np
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        result = rename(os.fsencode(source), os.fsencode(destination), 4)  # RENAME_EXCL
    else:
        raise OutputError("This platform does not support atomic directory publication without replacement")
    if result:
        error = ctypes.get_errno()
        if error in (errno.EEXIST, errno.ENOTEMPTY):
            raise OutputError("Output must be a new non-symlink directory")
        raise OSError(error, os.strerror(error), str(destination))


def publish(report: dict, destination: Path) -> None:
    if destination.exists() or any(p.is_symlink() for p in (destination, *destination.parents)):
        raise OutputError("Output must be a new non-symlink directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(report, indent=2, allow_nan=False).encode() + b"\n"
    if len(encoded) > 32 * 1024 * 1024:
        raise OutputError("Result exceeds 32 MiB")
    temporary = Path(tempfile.mkdtemp(prefix=".repeattyper-", dir=destination.parent))
    try:
        (temporary / "report.json").write_bytes(encoded)
        table = io.StringIO(newline="")
        fields = ["id", "sequence", "status", "prediction", "winner", "native_score", "reason"]
        writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["results"])
        (temporary / "predictions.tsv").write_text(table.getvalue(), encoding="ascii")
        _rename_new_directory(temporary, destination)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
