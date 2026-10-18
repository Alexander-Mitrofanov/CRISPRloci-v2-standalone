"""Bounded, checksum-verified operator downloads with real byte progress."""
from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
from urllib.parse import urlparse
from urllib.request import Request, urlopen


def file_digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1024 * 1024):
            value.update(block)
    return value.hexdigest()


def download(url, destination, expected, report, cancelled=lambda: False, *, limit=256 * 1024**2):
    destination = Path(destination)
    if urlparse(url).scheme != "https":
        raise ValueError("Installer downloads require HTTPS")
    if destination.is_symlink():
        raise ValueError("Download target must not be a symlink")
    if destination.exists():
        if file_digest(destination) != expected:
            raise ValueError("Cached download checksum differs; remove that cache file before retrying")
        size = destination.stat().st_size
        report(destination.name, size, size, "bytes")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".partial")
    if temporary.is_symlink():
        raise ValueError("Download staging file must not be a symlink")
    temporary.unlink(missing_ok=True)
    digest = hashlib.sha256()
    try:
        request = Request(url, headers={"User-Agent": "CRISPRloci-standalone-installer", "Accept-Encoding": "identity"})
        with urlopen(request, timeout=60) as response, temporary.open("xb") as output:
            if urlparse(response.geturl()).scheme != "https":
                raise ValueError("Download redirected away from HTTPS")
            raw_length = response.headers.get("Content-Length", "")
            length = int(raw_length) if raw_length.isdecimal() else None
            if length is not None and length > limit:
                raise ValueError("Download is larger than the supported limit")
            report(destination.name, 0, length, "bytes")
            size = 0
            while block := response.read(1024 * 1024):
                if cancelled():
                    raise InterruptedError("Installation cancelled")
                size += len(block)
                if size > limit or shutil.disk_usage(destination.parent).free < 1024**3:
                    raise ValueError("Download exceeds its byte or free-disk budget")
                digest.update(block)
                output.write(block)
                report(destination.name, size, length, "bytes")
        if not size or (length is not None and size != length) or digest.hexdigest() != expected:
            raise ValueError("Download is incomplete or differs from its pinned checksum")
        temporary.replace(destination)
        report(destination.name, size, size, "bytes")
        return destination
    finally:
        temporary.unlink(missing_ok=True)
