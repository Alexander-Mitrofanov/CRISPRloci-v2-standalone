"""Offline, bounded validation of immutable, non-executable reference data.

The checksums detect corruption and identify content; they do not authenticate
the publisher. Supply bundles only from a provenance source you intend to use.
"""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
from importlib.resources import files as package_files
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any

from crisprleader.errors import ReferenceError

SCHEMA_VERSION = "1.0.0"
MAX_MANIFEST_BYTES = 8 * 1024 * 1024
MAX_FILE_BYTES = 64 * 1024 * 1024
MAX_BUNDLE_BYTES = 256 * 1024 * 1024
MAX_FILES = 4096
MAX_FAMILIES = 4096
MAX_MOTIFS = 256
MAX_MOTIF_LENGTH = 128
DNA = frozenset("ACGTRYSWKMBDHVN")
BASES = "ACGT"
IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
MANIFEST_KEYS = frozenset({
    "schema_version", "bundle_id", "version", "families", "motifs", "provenance",
    "files", "policy",
})


def _fail(message: str) -> None:
    raise ReferenceError(message)


def _object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ReferenceError(f"{label} must be an object")
    return value


def _text(value: Any, label: str, *, maximum: int = 4096) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum or "\x00" in value:
        raise ReferenceError(f"{label} must be a nonempty string of at most {maximum} characters")
    return value


def _identifier(value: Any, label: str) -> str:
    text = _text(value, label, maximum=200)
    if not IDENTIFIER.fullmatch(text):
        raise ReferenceError(f"{label} must be a portable identifier")
    return text


def _number(value: Any, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ReferenceError(f"{label} must be a finite number")
    return float(value)


def _strings(value: Any, label: str, *, nonempty: bool = True) -> list[str]:
    if not isinstance(value, list) or (nonempty and not value):
        raise ReferenceError(f"{label} must be {'a nonempty' if nonempty else 'a'} list")
    result = [_text(item, label) for item in value]
    if len(result) != len(set(result)):
        raise ReferenceError(f"{label} must not contain duplicates")
    return result


def _dna_strings(value: Any, label: str) -> None:
    sequences = _strings(value, label)
    if len(sequences) > 10000:
        _fail(f"{label} exceeds the 10000-sequence limit")
    for sequence in sequences:
        if len(sequence) > 10000 or not set(sequence) <= DNA:
            _fail(f"{label} requires uppercase IUPAC DNA, at most 10000 bases per sequence")


def _keys(record: dict[str, Any], required: set[str], optional: set[str], label: str) -> None:
    missing = required - record.keys()
    unknown = record.keys() - required - optional
    if missing or unknown:
        _fail(f"{label} keys invalid: missing={sorted(missing)}, unknown={sorted(unknown)}")


def _json_tree(value: Any, depth: int = 0) -> None:
    if depth > 64:
        _fail("reference JSON nesting exceeds 64 levels")
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                _fail("reference JSON object keys must be strings")
            _json_tree(item, depth + 1)
    elif isinstance(value, list):
        for item in value:
            _json_tree(item, depth + 1)
    elif type(value) is float and not math.isfinite(value):
        _fail("nonfinite reference JSON numbers are forbidden")
    elif value is not None and type(value) not in (str, int, float, bool):
        _fail("reference content must use JSON-compatible values")


def validate_motif(value: Any) -> dict[str, Any]:
    """Check a generic position-probability matrix and its explicit scan scope."""
    _json_tree(value)
    motif = _object(value, "motif")
    _keys(
        motif,
        {"id", "domain", "pwm", "background", "strands", "citations", "metadata"},
        {"score_scale", "report_pvalue", "junction_end_range"},
        "motif",
    )
    _identifier(motif["id"], "motif.id")
    if motif["domain"] not in ("leader", "junction", "repeat"):
        _fail("motif.domain must be leader, junction or repeat")
    rows = motif["pwm"]
    if not isinstance(rows, list) or not 1 <= len(rows) <= MAX_MOTIF_LENGTH:
        _fail(f"motif.pwm must have 1..{MAX_MOTIF_LENGTH} probability rows")
    for label, row in [("background", motif["background"])] + [
        (f"pwm[{index}]", row) for index, row in enumerate(rows)
    ]:
        row = _object(row, f"motif.{label}")
        if set(row) != set(BASES):
            _fail(f"motif.{label} requires exactly A,C,G,T")
        probabilities = [_number(row[base], f"motif.{label}.{base}") for base in BASES]
        if any(number <= 0 or number > 1 for number in probabilities):
            _fail(f"motif.{label} probabilities must be in (0,1]; zero requires explicit smoothing")
        if not math.isclose(math.fsum(probabilities), 1.0, rel_tol=0, abs_tol=1e-10):
            _fail(f"motif.{label} probabilities must sum to 1")
    strands = motif["strands"]
    if (
        not isinstance(strands, list) or not strands
        or any(not isinstance(strand, str) for strand in strands)
        or len(strands) != len(set(strands))
    ):
        _fail("motif.strands must be a nonempty unique list")
    if any(strand not in ("+", "-") for strand in strands):
        _fail("motif.strands accepts only + and -")
    scale = motif.get("score_scale", 100)
    if type(scale) is not int or not 1 <= scale <= 1000:
        _fail("motif.score_scale must be an integer in 1..1000")
    threshold = _number(motif.get("report_pvalue", 0.001), "motif.report_pvalue")
    if not 0 <= threshold <= 1:
        _fail("motif.report_pvalue must be in [0,1]")
    if "junction_end_range" in motif:
        bounds = motif["junction_end_range"]
        if (
            not isinstance(bounds, list) or len(bounds) != 2
            or any(type(bound) is not int for bound in bounds)
            or bounds[0] > bounds[1] or any(abs(bound) > 100000 for bound in bounds)
        ):
            _fail("motif.junction_end_range requires ordered inclusive integer offsets")
    _strings(motif["citations"], "motif.citations")
    _object(motif["metadata"], "motif.metadata")
    return motif


def validate_manifest(value: Any) -> dict[str, Any]:
    """Validate manifest content without touching paths or executing model files."""
    _json_tree(value)
    manifest = _object(value, "reference manifest")
    if set(manifest) != MANIFEST_KEYS:
        _fail(
            "reference manifest keys invalid: "
            f"missing={sorted(MANIFEST_KEYS - manifest.keys())}, "
            f"unknown={sorted(manifest.keys() - MANIFEST_KEYS)}"
        )
    if manifest["schema_version"] != SCHEMA_VERSION:
        _fail(f"unsupported reference schema; expected {SCHEMA_VERSION}")
    _identifier(manifest["bundle_id"], "bundle_id")
    _identifier(manifest["version"], "version")
    families = manifest["families"]
    if not isinstance(families, list) or len(families) > MAX_FAMILIES:
        _fail(f"families must be a list with at most {MAX_FAMILIES} entries")
    seen: set[str] = set()
    for family in families:
        family = _object(family, "family")
        _keys(
            family,
            {"id", "repeat_sequences", "core_sequences", "boundary_kind",
             "training_locus_ids", "citations", "metadata"},
            {"description", "model_file", "model_format", "calibration"},
            "family",
        )
        identifier = _identifier(family["id"], "family.id")
        if identifier in seen:
            _fail(f"duplicate family id {identifier}")
        seen.add(identifier)
        _dna_strings(family["repeat_sequences"], "family.repeat_sequences")
        _dna_strings(family["core_sequences"], "family.core_sequences")
        if family["boundary_kind"] not in (
            "comparative_conserved_segment", "experimentally_interrogated_segment", "unresolved"
        ):
            _fail("family.boundary_kind is unsupported")
        _strings(family["training_locus_ids"], "family.training_locus_ids")
        _strings(family["citations"], "family.citations")
        _object(family["metadata"], "family.metadata")
        if "description" in family:
            _text(family["description"], "family.description")
        if ("model_file" in family) != ("model_format" in family):
            _fail("family.model_file and model_format must be supplied together")
        if "model_file" in family:
            member = _safe_member(family["model_file"])
            if member not in _object(manifest["files"], "files"):
                _fail(f"family model_file is not checksum-listed: {member}")
            if family["model_format"] not in ("hmmer3-dna", "json"):
                _fail("family.model_format must be hmmer3-dna or json; executable models forbidden")
        if "calibration" in family:
            _object(family["calibration"], "family.calibration")
    motifs = manifest["motifs"]
    if not isinstance(motifs, list) or len(motifs) > MAX_MOTIFS:
        _fail(f"motifs must be a list with at most {MAX_MOTIFS} entries")
    seen.clear()
    for motif in motifs:
        motif = validate_motif(motif)
        if motif["id"] in seen:
            _fail(f"duplicate motif id {motif['id']}")
        seen.add(motif["id"])
    provenance = _object(manifest["provenance"], "provenance")
    for key in ("build", "sources", "rights", "limitations"):
        if key not in provenance:
            _fail(f"provenance.{key} is required")
    _object(provenance["build"], "provenance.build")
    if not isinstance(provenance["sources"], list) or not isinstance(provenance["rights"], list):
        _fail("provenance.sources and provenance.rights must be lists")
    _strings(provenance["limitations"], "provenance.limitations", nonempty=False)
    _object(manifest["policy"], "policy")
    checksums = _object(manifest["files"], "files")
    if len(checksums) > MAX_FILES:
        _fail(f"files exceeds {MAX_FILES} entries")
    for name, digest in checksums.items():
        _safe_member(name)
        if not isinstance(digest, str) or not SHA256.fullmatch(digest):
            _fail(f"files[{name}] requires a lowercase SHA256 digest")
    return manifest


def _safe_member(value: Any) -> str:
    name = _text(value, "reference member path", maximum=1024)
    path = PurePosixPath(name)
    if (
        path.is_absolute() or "\\" in name or ":" in name or "//" in name
        or any(part in (".", "..", "") for part in name.split("/"))
        or name == "manifest.json" or any(ord(character) < 32 for character in name)
    ):
        _fail(f"unsafe reference member path: {name!r}")
    if path.suffix.lower() in {".pkl", ".pickle", ".py", ".pyc", ".so", ".dll", ".sh"}:
        _fail(f"executable reference member format is forbidden: {name}")
    return name


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    _fail(f"nonfinite JSON constant is forbidden: {value}")


def _read_regular(path: Path, limit: int, *, root_fd: int | None = None) -> bytes:
    directory_fds: list[int] = []
    try:
        flags = os.O_RDONLY | os.O_NONBLOCK | getattr(os, "O_NOFOLLOW", 0)
        if root_fd is None:
            descriptor = os.open(path, flags)
        else:
            parent_fd = root_fd
            for part in path.parts[:-1]:
                parent_fd = os.open(
                    part, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=parent_fd,
                )
                directory_fds.append(parent_fd)
            descriptor = os.open(path.name, flags, dir_fd=parent_fd)
        with os.fdopen(descriptor, "rb") as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode):
                _fail(f"reference member is not a regular file: {path.name}")
            if before.st_size > limit:
                _fail(f"reference member exceeds byte limit {limit}: {path.name}")
            content = stream.read(limit + 1)
            after = os.fstat(stream.fileno())
            if len(content) > limit:
                _fail(f"reference member exceeds byte limit {limit}: {path.name}")
            if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                after.st_size, after.st_mtime_ns, after.st_ctime_ns
            ):
                _fail(f"reference member changed during verification: {path.name}")
            return content
    except OSError as error:
        raise ReferenceError(f"cannot read reference member {path.name}: {error.strerror}") from error
    finally:
        for directory_fd in reversed(directory_fds):
            os.close(directory_fd)


def _inventory(root: Path) -> set[str]:
    result: set[str] = set()
    directories = 0
    for current, directory_names, filenames in os.walk(root, followlinks=False):
        directories += 1
        if directories > MAX_FILES:
            _fail("reference bundle directory limit exceeded")
        for name in directory_names + filenames:
            member = Path(current) / name
            if member.is_symlink():
                _fail(f"symlinks are forbidden in reference bundles: {member.relative_to(root)}")
        for name in filenames:
            member = Path(current) / name
            if not member.is_file():
                _fail(f"nonregular reference member: {member.relative_to(root)}")
            relative = member.relative_to(root).as_posix()
            if relative != "manifest.json":
                _safe_member(relative)
                result.add(relative)
            if len(result) > MAX_FILES:
                _fail("reference bundle file limit exceeded")
    return result


def load_bundle(path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """Verify every member and return a plain manifest mapping plus its receipt.

    ``None`` addresses the installed ``data/default`` directory. A missing
    default is a setup error, never an unsupported biological-family result.
    No network, deserialization of executable objects, or reference updates
    occur. The verified member bytes are retained, avoiding later path races.
    """
    if path is None:
        root = Path(str(package_files("crisprleader").joinpath("data", "default")))
    else:
        root = Path(path)
    if root.name == "manifest.json":
        root = root.parent
    if root.is_symlink() or not root.is_dir():
        raise ReferenceError(
            f"reference bundle directory unavailable: {root}; "
            "install a verified biological reference bundle or provide --references"
        )
    root = root.absolute()
    manifest_bytes = _read_regular(root / "manifest.json", MAX_MANIFEST_BYTES)
    try:
        value = json.loads(
            manifest_bytes.decode("utf-8"), object_pairs_hook=_no_duplicates,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, ValueError, RecursionError) as error:
        raise ReferenceError(f"invalid reference manifest JSON: {error}") from error
    manifest = validate_manifest(value)
    inventory = _inventory(root)
    expected = set(manifest["files"])
    if inventory != expected:
        _fail(
            "reference file inventory mismatch: "
            f"missing={sorted(expected - inventory)[:10]}, "
            f"unlisted={sorted(inventory - expected)[:10]}"
        )
    members: dict[str, bytes] = {}
    total = len(manifest_bytes)
    try:
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ReferenceError(f"cannot open reference root: {error.strerror}") from error
    try:
        for name, expected_digest in sorted(manifest["files"].items()):
            # Open every path component relative to one root descriptor. A
            # concurrent replacement with a symlink cannot redirect reads.
            content = _read_regular(Path(name), MAX_FILE_BYTES, root_fd=root_fd)
            total += len(content)
            if total > MAX_BUNDLE_BYTES:
                _fail(f"reference bundle exceeds byte limit {MAX_BUNDLE_BYTES}")
            digest = hashlib.sha256(content).hexdigest()
            if digest != expected_digest:
                _fail(f"reference checksum mismatch: {name}")
            members[name] = content
    finally:
        os.close(root_fd)
    manifest["_receipt"] = {
        "schema_version": SCHEMA_VERSION,
        "bundle_id": manifest["bundle_id"],
        "version": manifest["version"],
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "files": dict(manifest["files"]),
        "verified_bytes": total,
    }
    manifest["_members"] = members
    return manifest


def reference_receipt(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Return a JSON-safe copy of the content-verification receipt."""
    receipt = bundle.get("_receipt")
    if not isinstance(receipt, dict):
        raise ReferenceError("bundle has no verification receipt; load it with load_bundle")
    return json.loads(json.dumps(receipt))


__all__ = ["load_bundle", "reference_receipt", "validate_manifest", "validate_motif"]
