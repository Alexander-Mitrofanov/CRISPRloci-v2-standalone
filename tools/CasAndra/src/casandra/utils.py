"""Small dependency-free utilities shared by the CasAndra runtime."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import os
import sys
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of *path* without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


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
        raise FileExistsError("This platform does not support atomic directory publication without replacement")
    if result:
        error = ctypes.get_errno()
        if error in (errno.EEXIST, errno.ENOTEMPTY):
            raise FileExistsError("Output must be a new non-symlink directory")
        raise OSError(error, os.strerror(error), str(destination))
