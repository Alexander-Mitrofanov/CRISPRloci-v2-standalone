"""Explicit operator-only acquisition of pinned NCBI viral RefSeq source files.

No job or search calls this module. Plans permit only the official release URLs;
downloads verify the pinned official catalog and every compressed source file.
"""
from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
from html.parser import HTMLParser
from http.client import HTTPException
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener
import zlib

from . import __version__
from .io import InputError, LimitError, atomic_directory, canonical, decode_json, digest, integer, read_bytes, write_json

BASE = "https://ftp.ncbi.nlm.nih.gov/refseq/release/"
VIRAL_URL = BASE + "viral/"
RELEASE_URL = BASE + "RELEASE_NUMBER"
PLAN_SCHEMA = {"name": "crisprspacer.ncbi-viral-plan", "version": "1.0.0"}
DOWNLOAD_SCHEMA = {"name": "crisprspacer.ncbi-viral-download", "version": "1.0.0"}
MAX_METADATA = 8 * 1024**2
MAX_SOURCE_BYTES = 8 * 1024**3
MAX_TOTAL_BYTES = 16 * 1024**3
MAX_UNCOMPRESSED_BYTES = 20 * 1024**3
MIN_FREE_BYTES = 20 * 1024**3
FILE_NAME = re.compile(r"viral\.[1-9][0-9]*\.[1-9][0-9]*\.genomic\.fna\.gz\Z")
RELEASE = re.compile(r"[1-9][0-9]{0,5}\Z")
MD5 = re.compile(r"[0-9a-f]{32}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise InputError("NCBI download redirects are not permitted")


def _allowed_url(url):
    if url in {VIRAL_URL, RELEASE_URL}:
        return True
    if url.startswith(VIRAL_URL) and FILE_NAME.fullmatch(url[len(VIRAL_URL):]):
        return True
    return bool(re.fullmatch(re.escape(BASE) + r"release-catalog/release[1-9][0-9]{0,5}\.files\.installed", url))


def _open(url):
    if not isinstance(url, str) or not _allowed_url(url):
        raise InputError("Only fixed official NCBI viral RefSeq URLs are permitted")
    request = Request(url, headers={"User-Agent": f"CRISPRspacer/{__version__} operator-reference-provisioning",
                                    "Accept-Encoding": "identity"})
    try:
        response = build_opener(_NoRedirect()).open(request, timeout=30)
    except (HTTPError, URLError, HTTPException, TimeoutError, OSError) as exc:
        raise InputError("Cannot retrieve the pinned NCBI resource") from exc
    if response.geturl() != url or response.status != 200:
        response.close()
        raise InputError("Unexpected NCBI response URL/status")
    return response


def _fetch(url, limit=MAX_METADATA):
    deadline = time.monotonic() + 300
    result = bytearray()
    try:
        with _open(url) as response:
            while chunk := response.read(min(1024 * 1024, limit + 1 - len(result))):
                result.extend(chunk)
                if len(result) > limit:
                    raise LimitError("NCBI metadata exceeds its byte budget")
                if time.monotonic() > deadline:
                    raise LimitError("NCBI metadata download deadline exceeded")
    except (OSError, HTTPException, TimeoutError) as exc:
        raise InputError("Interrupted NCBI metadata download") from exc
    return bytes(result)


def _release(data):
    try:
        value = data.decode("ascii").strip()
    except UnicodeError as exc:
        raise InputError("Malformed NCBI release number") from exc
    if not RELEASE.fullmatch(value):
        raise InputError("Malformed NCBI release number")
    return value


class _Listing(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.names = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for key, value in attrs:
                if key == "href" and value is not None and FILE_NAME.fullmatch(value):
                    self.names.append(value)


def _listing(data):
    parser = _Listing()
    try:
        parser.feed(data.decode("utf-8"))
        parser.close()
    except (UnicodeError, ValueError) as exc:
        raise InputError("Malformed NCBI viral listing") from exc
    if not parser.names or len(parser.names) != len(set(parser.names)) or len(parser.names) > 1000:
        raise InputError("NCBI viral listing has an empty, duplicate, or excessive source inventory")
    return sorted(parser.names)


def _catalog(data):
    """Parse official MD5 rows; catalog format is deliberately fail-closed."""
    try:
        lines = data.decode("ascii").splitlines()
    except UnicodeError as exc:
        raise InputError("Malformed NCBI release catalog") from exc
    result = {}
    for line in lines:
        tokens = line.split()
        names = [token.removeprefix("*").removeprefix("viral/") for token in tokens]
        names = [name for name in names if FILE_NAME.fullmatch(name)]
        if not names:
            continue
        if (len(tokens) != 2 or len(names) != 1 or not MD5.fullmatch(tokens[0].lower())
                or tokens[1] != names[0] or names[0] in result):
            raise InputError("Malformed or duplicate viral source checksum in NCBI catalog")
        result[names[0]] = {"name": names[0], "url": VIRAL_URL + names[0], "md5": tokens[0].lower(), "bytes": None}
    if not result or len(result) > 1000:
        raise InputError("NCBI catalog has no supported viral genomic FASTA sources")
    return [result[name] for name in sorted(result)]


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _atomic_plan(destination, value):
    destination = Path(destination).absolute()
    if any(p.is_symlink() for p in (destination, *destination.parents)):
        raise InputError("Output paths must not traverse symlinks")
    if destination.exists():
        raise InputError("Plan output already exists; choose a new destination")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".crisprspacer-plan-", dir=destination.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(canonical(value))
            handle.flush()
            os.fsync(handle.fileno())
        # Hard link publishes complete bytes without replacing a concurrent writer.
        os.link(temporary, destination)
        directory_fd = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except FileExistsError as exc:
        raise InputError("Plan output appeared during execution; refusing replacement") from exc
    finally:
        temporary.unlink()


def plan_ncbi_viral(destination: Path):
    """Read bounded metadata only and atomically save the complete release plan."""
    release_bytes = _fetch(RELEASE_URL, 128)
    release = _release(release_bytes)
    listing_bytes = _fetch(VIRAL_URL)
    names = _listing(listing_bytes)
    catalog_url = BASE + f"release-catalog/release{release}.files.installed"
    catalog_bytes = _fetch(catalog_url)
    files = _catalog(catalog_bytes)
    if names != [row["name"] for row in files]:
        raise InputError("NCBI viral listing and release catalog disagree; retry after release publication")
    if (_fetch(catalog_url) != catalog_bytes or _fetch(VIRAL_URL) != listing_bytes
            or _fetch(RELEASE_URL, 128) != release_bytes):
        raise InputError("NCBI release metadata changed during planning; no plan published")
    plan = {"schema": PLAN_SCHEMA, "created_at_utc": _now(), "source_release": release,
            "source_url": VIRAL_URL, "files": files,
            "catalog": {"url": catalog_url, "sha256": digest(catalog_bytes), "bytes": len(catalog_bytes)},
            "listing": {"url": VIRAL_URL, "sha256": digest(listing_bytes), "bytes": len(listing_bytes)},
            "release_number": {"url": RELEASE_URL, "sha256": digest(release_bytes), "bytes": len(release_bytes)}}
    _atomic_plan(destination, plan)
    return plan


def _validate_plan(plan):
    expected = {"schema", "created_at_utc", "source_release", "source_url", "files", "catalog", "listing", "release_number"}
    if not isinstance(plan, dict) or set(plan) != expected or plan["schema"] != PLAN_SCHEMA:
        raise InputError("Unsupported NCBI viral plan schema")
    release = plan["source_release"]
    if not isinstance(release, str) or not RELEASE.fullmatch(release) or plan["source_url"] != VIRAL_URL:
        raise InputError("Invalid NCBI release identity")
    if not isinstance(plan["created_at_utc"], str) or len(plan["created_at_utc"]) > 64:
        raise InputError("Invalid NCBI plan timestamp")
    urls = {"catalog": BASE + f"release-catalog/release{release}.files.installed",
            "listing": VIRAL_URL, "release_number": RELEASE_URL}
    for key, url in urls.items():
        row = plan[key]
        if not isinstance(row, dict) or set(row) != {"url", "sha256", "bytes"} or row["url"] != url:
            raise InputError("Invalid NCBI plan metadata URL/fields")
        if not isinstance(row["sha256"], str) or not SHA256.fullmatch(row["sha256"]):
            raise InputError("Invalid NCBI metadata SHA256")
        integer(row["bytes"], "metadata bytes", 1, 128 if key == "release_number" else MAX_METADATA)
    files = plan["files"]
    if not isinstance(files, list) or not 1 <= len(files) <= 1000:
        raise InputError("Invalid NCBI source file inventory")
    names = []
    for row in files:
        if not isinstance(row, dict) or set(row) != {"name", "url", "md5", "bytes"}:
            raise InputError("Invalid NCBI source file fields")
        name = row["name"]
        if not isinstance(name, str) or not FILE_NAME.fullmatch(name) or row["url"] != VIRAL_URL + name:
            raise InputError("Only named official NCBI viral genomic FASTA files are permitted")
        if not isinstance(row["md5"], str) or not MD5.fullmatch(row["md5"]):
            raise InputError("Invalid NCBI source MD5")
        if row["bytes"] is not None:
            integer(row["bytes"], "source bytes", 1, MAX_SOURCE_BYTES)
        names.append(name)
    if names != sorted(set(names)):
        raise InputError("NCBI source file inventory must be sorted and unique")
    return plan


def _verify_gzip(path, deadline, remaining):
    total = 0
    try:
        with gzip.open(path, "rb") as handle:
            while chunk := handle.read(1024 * 1024):
                total += len(chunk)
                if total > remaining:
                    raise LimitError("NCBI source exceeds decompressed byte budget")
                if time.monotonic() > deadline:
                    raise LimitError("NCBI source validation deadline exceeded")
    except (OSError, EOFError, zlib.error) as exc:
        raise InputError("NCBI source gzip integrity check failed") from exc
    if not total:
        raise InputError("NCBI source gzip is empty")
    return total


def download_ncbi_viral(plan: Path, destination: Path, min_free_bytes=MIN_FREE_BYTES):
    """Download only the pinned plan into a new directory, with no automatic update."""
    integer(min_free_bytes, "minimum free bytes", 0, 1024**5)
    plan_bytes = read_bytes(plan, MAX_METADATA)
    pinned = _validate_plan(decode_json(plan_bytes))
    catalog_bytes = _fetch(pinned["catalog"]["url"])
    if (digest(catalog_bytes) != pinned["catalog"]["sha256"]
            or len(catalog_bytes) != pinned["catalog"]["bytes"]
            or _catalog(catalog_bytes) != pinned["files"]):
        raise InputError("Pinned NCBI catalog or source inventory changed; refusing download")
    deadline = time.monotonic() + 24 * 60 * 60
    with atomic_directory(destination) as work:
        if shutil.disk_usage(work).free < min_free_bytes:
            raise LimitError("Insufficient free disk space for viral reference acquisition")
        total = 0
        decompressed = 0
        files = []
        for row in pinned["files"]:
            sha256 = hashlib.sha256()
            md5 = hashlib.md5(usedforsecurity=False)
            size = 0
            path = work / row["name"]
            try:
                with _open(row["url"]) as response, path.open("xb") as output:
                    while chunk := response.read(1024 * 1024):
                        size += len(chunk)
                        total += len(chunk)
                        if size > MAX_SOURCE_BYTES or total > MAX_TOTAL_BYTES:
                            raise LimitError("NCBI source download exceeds its byte budget")
                        if time.monotonic() > deadline:
                            raise LimitError("NCBI source download deadline exceeded")
                        if shutil.disk_usage(work).free < min(256 * 1024**2, min_free_bytes):
                            raise LimitError("Insufficient free disk space during NCBI download")
                        md5.update(chunk)
                        sha256.update(chunk)
                        output.write(chunk)
            except (OSError, HTTPException, TimeoutError) as exc:
                raise InputError("Interrupted NCBI source download") from exc
            if not size or md5.hexdigest() != row["md5"] or row["bytes"] not in {None, size}:
                raise InputError("NCBI source size/official MD5 integrity check failed")
            decompressed += _verify_gzip(path, deadline, MAX_UNCOMPRESSED_BYTES - decompressed)
            files.append({**row, "bytes": size, "sha256": sha256.hexdigest()})
        # A changed historical catalog must not silently bless a mixed snapshot.
        if _fetch(pinned["catalog"]["url"]) != catalog_bytes:
            raise InputError("NCBI catalog changed during download; no source directory published")
        (work / "plan.json").write_bytes(plan_bytes)
        (work / "release.files.installed").write_bytes(catalog_bytes)
        receipt = {"schema": DOWNLOAD_SCHEMA, "created_at_utc": _now(),
                   "source_release": pinned["source_release"], "source_url": VIRAL_URL,
                   "plan_file": "plan.json", "plan_sha256": digest(plan_bytes),
                   "catalog": pinned["catalog"], "files": files}
        write_json(work / "download.json", receipt)
    return receipt
