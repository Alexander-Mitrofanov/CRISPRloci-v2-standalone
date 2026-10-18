"""Strict operator-owned configuration.

No setting in this module is populated from an anonymous job request.  In
particular, executable paths, resource budgets, categories, and timeouts are
fixed by the service operator.
"""

from __future__ import annotations

from . import __version__

from dataclasses import dataclass, field
import hashlib
import hmac
import ipaddress
import json
import stat
import os
import re
from pathlib import Path
from typing import Iterable
from urllib.parse import urlsplit


def _integer(name: str, default: int, minimum: int, maximum: int) -> int:
    raw = os.getenv(name)
    try:
        value = default if raw is None else int(raw)
    except ValueError as error:
        raise ValueError(f"{name} must be an integer") from error
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return value


_SHA256_HEX = re.compile(r"[0-9a-f]{64}")


def _sha256_identifier(name: str, value: object) -> str:
    if value is None:
        raise ValueError(f"{name} is required")
    if not isinstance(value, str) or _SHA256_HEX.fullmatch(value) is None:
        raise ValueError(f"{name} must be exactly 64 lowercase hexadecimal characters")
    return value


def _boolean(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    normalized = raw.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{name} must be true or false")


_DNS_LABEL = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?$")


def _canonical_origin(value: str) -> str:
    if not value or any(
        ord(character) <= 32 or ord(character) == 127 for character in value
    ):
        raise ValueError("CORS origins cannot contain whitespace or control characters")
    if "*" in value:
        raise ValueError("CORS origins cannot contain wildcards")
    parsed = urlsplit(value)
    try:
        hostname = parsed.hostname
        port = parsed.port
    except ValueError as error:
        raise ValueError(f"invalid CORS origin: {value!r}") from error
    if (
        hostname is None
        or not parsed.netloc
        or parsed.path
        or parsed.query
        or parsed.fragment
        or parsed.username
        or parsed.password
    ):
        raise ValueError(f"invalid exact CORS origin: {value!r}")

    address: ipaddress.IPv4Address | ipaddress.IPv6Address | None = None
    try:
        address = ipaddress.ip_address(hostname)
    except ValueError:
        if len(hostname) > 253 or any(
            not _DNS_LABEL.fullmatch(label) for label in hostname.split(".")
        ):
            raise ValueError(f"invalid CORS origin hostname: {hostname!r}") from None

    loopback = hostname == "localhost" or bool(address and address.is_loopback)
    if parsed.scheme != "https" and not (parsed.scheme == "http" and loopback):
        raise ValueError(
            "CORS origins must use HTTPS; HTTP is allowed only for exact "
            "loopback development origins"
        )

    if address is not None:
        rendered_host = (
            f"[{address.compressed}]"
            if isinstance(address, ipaddress.IPv6Address)
            else address.compressed
        )
    else:
        rendered_host = hostname.lower()
    default_port = 443 if parsed.scheme == "https" else 80
    rendered_port = "" if port is None else f":{port}"
    if port == default_port:
        raise ValueError("CORS origins must omit their scheme's default port")
    canonical = f"{parsed.scheme}://{rendered_host}{rendered_port}"
    if value != canonical:
        raise ValueError(f"CORS origin must be canonical: {canonical!r}")
    return canonical


def _origins(raw: str) -> tuple[str, ...]:
    origins: list[str] = []
    for item in raw.split(","):
        if not item.strip():
            continue
        if item != item.strip():
            raise ValueError("CRISPR_WEB_CORS_ORIGINS entries cannot have padding")
        origins.append(_canonical_origin(item))
    if len(origins) != len(set(origins)):
        raise ValueError("CRISPR_WEB_CORS_ORIGINS contains duplicates")
    return tuple(origins)


def _networks(raw: str) -> tuple[ipaddress.IPv4Network | ipaddress.IPv6Network, ...]:
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(ipaddress.ip_network(item, strict=False))
    return tuple(values)


def _command(raw: str) -> tuple[str, ...]:
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError(
            "CRISPR_WEB_INTEGRATION_COMMAND must be a JSON array"
        ) from error
    if (
        not isinstance(value, list)
        or not value
        or not all(
            isinstance(part, str) and part and "\0" not in part for part in value
        )
    ):
        raise ValueError(
            "CRISPR_WEB_INTEGRATION_COMMAND must be a non-empty JSON string array"
        )
    return tuple(value)


def _safe_root(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved == Path(resolved.anchor):
        raise ValueError("CRISPR_WEB_DATA_ROOT cannot be a filesystem root")
    return resolved


_SCIENTIFIC_ARTIFACTS = {
    "crispridentify": ("source_tree", "source_tree_v1"),
    "spacerplacer": ("wheel", "file_bytes"),
    "crispr_evor": ("wheel", "file_bytes"),
    "integration": ("wheel", "file_bytes"),
}
_SCIENTIFIC_ARTIFACTS_V2 = {
    **{name: ("git_source", "source_inventory_v1") for name in
       ("crispridentify", "spacerplacer", "crispr_evor", "casandra", "crispr_tracrrna")},
    "integration": ("wheel", "file_bytes"),
}
_MANIFEST_IDENTITY_FIELDS = frozenset(
    {"name", "display_name", "version", "kind", "path", "sha256", "sha256_scope"}
)
_MANIFEST_MAX_BYTES = 1_000_000


def _strict_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"scientific release manifest repeats key {key!r}")
        result[key] = value
    return result


def _manifest_string(
    mapping: dict[str, object], key: str, *, maximum: int = 500
) -> str:
    value = mapping.get(key)
    if (
        not isinstance(value, str)
        or not value
        or len(value) > maximum
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError(f"scientific release manifest field {key!r} is invalid")
    return value


def _manifest_path(mapping: dict[str, object], key: str = "path") -> str:
    value = _manifest_string(mapping, key, maximum=4_096)
    if not Path(value).is_absolute():
        raise ValueError(f"scientific release manifest field {key!r} must be absolute")
    return value


def _validated_scientific_manifest(
    document: object,
) -> dict[str, dict[str, str]]:
    required_fields = {
        "schema",
        "schema_version",
        "release_id",
        "scientific_artifacts",
    }
    if not isinstance(document, dict) or not required_fields.issubset(document):
        raise ValueError("scientific release manifest is missing required fields")
    if document.get("schema") != "crispr-web-scientific-release":
        raise ValueError("scientific release manifest schema is unsupported")
    if document.get("schema_version") not in {"1.0.0", "2.0.0", "3.0.0", "4.0.0", "5.0.0", "6.0.0", "7.0.0", "8.0.0", "9.0.0", "10.0.0"}:
        raise ValueError("scientific release manifest schema version is unsupported")
    _manifest_string(document, "release_id")

    expected_artifacts = _SCIENTIFIC_ARTIFACTS_V2 if document["schema_version"] == "2.0.0" else _SCIENTIFIC_ARTIFACTS
    if document["schema_version"] == "3.0.0":
        expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "crisprleader": ("git_source", "source_inventory_v1")}
    artifacts = document.get("scientific_artifacts")
    if document["schema_version"] == "4.0.0":
        expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "crisprrepeat": ("git_source", "source_inventory_v1")}
        if isinstance(artifacts, dict) and "crisprleader" in artifacts:
            expected_artifacts["crisprleader"] = ("git_source", "source_inventory_v1")
    if document["schema_version"] in {"5.0.0", "6.0.0", "7.0.0", "8.0.0"}:
        expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "crisprspacer": ("git_source", "source_inventory_v1")}
        if document["schema_version"] in {"6.0.0", "7.0.0", "8.0.0"}:
            expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "crisprmap": ("git_source", "source_inventory_v1")}
            if isinstance(artifacts, dict) and "crisprspacer" in artifacts:
                expected_artifacts["crisprspacer"] = ("git_source", "source_inventory_v1")
        if document["schema_version"] == "7.0.0":
            expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "repeattyper": ("git_source", "source_inventory_v1")}
            for optional in ("crisprspacer", "crisprmap"):
                if isinstance(artifacts, dict) and optional in artifacts:
                    expected_artifacts[optional] = ("git_source", "source_inventory_v1")
        for optional in ("crisprleader", "crisprrepeat"):
            if isinstance(artifacts, dict) and optional in artifacts:
                expected_artifacts[optional] = ("git_source", "source_inventory_v1")
    if document["schema_version"] == "8.0.0":
        expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "spacepharer": ("git_source", "source_inventory_v1")}
        for optional in ("crisprleader", "crisprrepeat", "crisprspacer", "crisprmap", "repeattyper"):
            if isinstance(artifacts,dict) and optional in artifacts:
                expected_artifacts[optional] = ("git_source", "source_inventory_v1")
    if document["schema_version"] == "9.0.0":
        expected_artifacts = {**_SCIENTIFIC_ARTIFACTS_V2, "cctk": ("git_source", "source_inventory_v1")}
        for optional in ("crisprleader", "crisprrepeat", "crisprspacer", "crisprmap", "repeattyper", "spacepharer"):
            if isinstance(artifacts, dict) and optional in artifacts:
                expected_artifacts[optional] = ("git_source", "source_inventory_v1")
    if document["schema_version"] == "10.0.0":
        additions = ("clinker", "phigaro", "isescan", "defensefinder", "rnaalifold")
        if not isinstance(artifacts, dict) or not any(name in artifacts for name in additions):
            raise ValueError("Context release requires at least one context tool")
        expected_artifacts = dict(_SCIENTIFIC_ARTIFACTS_V2)
        for optional in ("crisprleader", "crisprrepeat", "crisprspacer", "crisprmap", "repeattyper", "spacepharer", "cctk", *additions):
            if optional in artifacts:
                expected_artifacts[optional] = ("git_source", "source_inventory_v1")
    if not isinstance(artifacts, dict) or set(artifacts) != set(expected_artifacts):
        raise ValueError("scientific release manifest artifact set is invalid")
    public_identities: dict[str, dict[str, str]] = {}
    for key, (expected_kind, expected_scope) in expected_artifacts.items():
        identity = artifacts.get(key)
        expected_fields = _MANIFEST_IDENTITY_FIELDS | ({"repository", "commit"} if expected_kind == "git_source" else set())
        if not isinstance(identity, dict) or set(identity) != expected_fields:
            raise ValueError(f"scientific release manifest identity {key!r} is invalid")
        if identity.get("name") != key:
            raise ValueError(
                f"scientific release manifest identity {key!r} is misnamed"
            )
        if identity.get("kind") != expected_kind:
            raise ValueError(
                f"scientific release manifest identity {key!r} has wrong kind"
            )
        if identity.get("sha256_scope") != expected_scope:
            raise ValueError(
                f"scientific release manifest identity {key!r} has wrong hash scope"
            )
        _manifest_path(identity)
        digest = _sha256_identifier(
            f"scientific_artifacts.{key}.sha256", identity.get("sha256")
        )
        public_identities[key] = {
            "name": key,
            "display_name": _manifest_string(identity, "display_name"),
            "version": _manifest_string(identity, "version"),
            "kind": expected_kind,
            "sha256": digest,
            "sha256_scope": expected_scope,
        }

        if expected_kind == "git_source":
            repository = _manifest_string(identity, "repository")
            commit = _manifest_string(identity, "commit")
            if not re.fullmatch(r"https://github\.com/[A-Za-z0-9_-]+/[A-Za-z0-9_.-]+\.git", repository) or not re.fullmatch(r"[0-9a-f]{40}", commit):
                raise ValueError("Invalid scientific repository identity")
            public_identities[key].update(repository=repository, commit=commit)

    return public_identities


def _load_scientific_release_manifest(
    path: Path, expected_sha256: str
) -> dict[str, dict[str, str]]:
    if not path.is_absolute():
        raise ValueError("CRISPR_RELEASE_MANIFEST must be an absolute path")
    try:
        resolved_path = path.resolve(strict=True)
    except OSError as error:
        raise ValueError("scientific release manifest is unavailable") from error
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(resolved_path, flags)
    except OSError as error:
        raise ValueError("scientific release manifest is unavailable") from error
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_size <= 0
            or before.st_size > _MANIFEST_MAX_BYTES
        ):
            raise ValueError(
                "scientific release manifest must be a bounded regular file"
            )
        payload = bytearray()
        while block := os.read(descriptor, min(65_536, _MANIFEST_MAX_BYTES + 1)):
            payload.extend(block)
            if len(payload) > _MANIFEST_MAX_BYTES:
                raise ValueError("scientific release manifest exceeds its size limit")
        after = os.fstat(descriptor)
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise ValueError("scientific release manifest changed while being read")
    finally:
        os.close(descriptor)
    actual_sha256 = hashlib.sha256(payload).hexdigest()
    if not hmac.compare_digest(actual_sha256, expected_sha256):
        raise ValueError(
            "scientific release manifest SHA-256 does not match deployment"
        )
    try:
        document = json.loads(
            payload.decode("utf-8"), object_pairs_hook=_strict_json_object
        )
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(
            "scientific release manifest is not strict UTF-8 JSON"
        ) from error
    return _validated_scientific_manifest(document)


@dataclass(frozen=True, slots=True)
class Settings:
    """All public limits and operator-controlled worker settings."""

    data_root: Path
    database_path: Path
    runner_config: Path
    integration_command: tuple[str, ...]
    cors_origins: tuple[str, ...]
    release_id: str
    backend_wheel_sha256: str
    scientific_release_manifest: Path
    scientific_release_manifest_sha256: str
    _scientific_release_identities_json: str = field(init=False, repr=False)
    trusted_proxy_networks: tuple[
        ipaddress.IPv4Network | ipaddress.IPv6Network, ...
    ] = ()
    service_version: str = __version__
    max_request_bytes: int = 4_500_000
    max_total_bases: int = 2_000_000
    max_record_bases: int = 2_000_000
    max_records: int = 20
    max_header_characters: int = 200
    max_queued_jobs: int = 8
    max_active_jobs: int = 8
    max_active_jobs_per_client: int = 2
    readiness_cache_seconds: int = 5
    max_job_attempts: int = 2
    rate_limit_submissions: int = 4
    rate_limit_window_seconds: int = 3_600
    retention_seconds: int = 86_400
    worker_cpu: int = 8
    identify_folder_jobs: int = 1
    stage_timeout_seconds: int = 7_200
    job_timeout_seconds: int = 21_900
    worker_poll_seconds: int = 2
    worker_lease_seconds: int = 30
    worker_stale_seconds: int = 90
    disable_bias_corrections: bool = False
    max_summary_arrays: int = 1_000
    max_summary_rows: int = 1_000
    max_bundle_files: int = 5_000
    max_artifact_file_bytes: int = 256_000_000
    max_archive_bytes: int = 512_000_000
    max_job_disk_bytes: int = 2_000_000_000
    max_job_entries: int = 50_000
    max_launcher_log_bytes: int = 20_000_000
    min_free_bytes: int = 2_000_000_000
    min_free_inodes: int = 100_000
    job_disk_check_seconds: int = 1
    global_rate_limit_submissions: int = 30
    global_rate_limit_window_seconds: int = 60
    hsts: bool = True
    association_reference_path: Path | None = None
    association_reference_sha256: str | None = None
    repeat_reference_path: Path | None = None
    repeat_reference_sha256: str | None = None
    viral_reference_path: Path | None = None
    viral_reference_sha256: str | None = None
    spacer_database_path: Path | None = None
    spacer_database_sha256: str | None = None
    catalog_database_path: Path | None = None
    catalog_search_index_path: Path | None = None

    def __post_init__(self) -> None:
        _sha256_identifier("release_id", self.release_id)
        _sha256_identifier("backend_wheel_sha256", self.backend_wheel_sha256)
        scientific_digest = _sha256_identifier(
            "scientific_release_manifest_sha256",
            self.scientific_release_manifest_sha256,
        )
        manifest_path = self.scientific_release_manifest.expanduser()
        object.__setattr__(self, "scientific_release_manifest", manifest_path)
        identities = _load_scientific_release_manifest(
            manifest_path,
            scientific_digest,
        )
        object.__setattr__(
            self,
            "_scientific_release_identities_json",
            json.dumps(identities, sort_keys=True, separators=(",", ":")),
        )
        object.__setattr__(self, "data_root", _safe_root(self.data_root))
        object.__setattr__(
            self, "database_path", self.database_path.expanduser().resolve()
        )
        if (self.spacer_database_path is None) != (self.spacer_database_sha256 is None):
            raise ValueError("Spacer database path and SHA-256 must be configured together")
        if self.spacer_database_path is not None:
            _sha256_identifier("spacer_database_sha256", self.spacer_database_sha256)
            path = self.spacer_database_path.expanduser()
            if not path.is_absolute() or not path.is_dir() or any(p.is_symlink() for p in (path, *path.parents)):
                raise ValueError("Spacer database must be an existing absolute directory without symlinks")
            if "crisprspacer" not in identities:
                raise ValueError("Configured spacer database requires a release-pinned CRISPRspacer tool")
            object.__setattr__(self, "spacer_database_path", path)
        if (self.viral_reference_path is None) != (self.viral_reference_sha256 is None):
            raise ValueError("Viral reference path and SHA-256 must be configured together")
        if self.viral_reference_path is not None:
            _sha256_identifier("viral_reference_sha256", self.viral_reference_sha256)
            path = self.viral_reference_path.expanduser()
            if not path.is_absolute() or not path.is_dir() or any(p.is_symlink() for p in (path, *path.parents)):
                raise ValueError("Viral reference must be an existing absolute directory without symlinks")
            if identities.get("crisprspacer", {}).get("version") not in {"0.3.0", "0.3.1"}:
                raise ValueError("Viral reference requires release-pinned CRISPRspacer 0.3.0 or 0.3.1")
            object.__setattr__(self, "viral_reference_path", path)
        if (self.association_reference_path is None) != (self.association_reference_sha256 is None):
            raise ValueError("Association panel path and digest must be configured together")
        if self.association_reference_path is not None:
            _sha256_identifier("association_reference_sha256", self.association_reference_sha256)
            path = self.association_reference_path.expanduser()
            if not path.is_absolute() or not path.is_dir() or any(p.is_symlink() for p in (path,*path.parents)):
                raise ValueError("Association panel must be an absolute directory without symlinks")
            if identities.get("spacepharer", {}).get("version") not in {"2.1.0", "2.1.1"}:
                raise ValueError("Association panel requires release-pinned SpacePHARER 2.1.0 or 2.1.1")
            object.__setattr__(self,"association_reference_path",path)
        if (self.repeat_reference_path is None) != (self.repeat_reference_sha256 is None):
            raise ValueError("Repeat reference path and SHA-256 must be configured together")
        if self.repeat_reference_path is not None:
            _sha256_identifier("repeat_reference_sha256", self.repeat_reference_sha256)
            path = self.repeat_reference_path.expanduser()
            if not path.is_absolute() or not path.is_file() or any(p.is_symlink() for p in (path, *path.parents)):
                raise ValueError("Repeat reference must be an absolute regular file without symlinks")
            if "crisprmap" not in identities:
                raise ValueError("Repeat reference requires release-pinned CRISPRmap")
            object.__setattr__(self, "repeat_reference_path", path)
        if self.catalog_database_path is not None:
            catalog_path = self.catalog_database_path.expanduser().absolute()
            if catalog_path.resolve() == self.database_path:
                raise ValueError("catalogue database must be separate from the job database")
            object.__setattr__(self, "catalog_database_path", catalog_path)
        if self.catalog_search_index_path is not None:
            object.__setattr__(self, "catalog_search_index_path", self.catalog_search_index_path.expanduser().absolute())
        object.__setattr__(
            self, "runner_config", self.runner_config.expanduser().resolve()
        )
        if self.worker_stale_seconds <= self.worker_lease_seconds:
            raise ValueError("worker_stale_seconds must exceed worker_lease_seconds")
        if self.job_timeout_seconds <= self.stage_timeout_seconds:
            raise ValueError("job_timeout_seconds must exceed stage_timeout_seconds")
        if self.identify_folder_jobs > self.worker_cpu:
            raise ValueError("identify_folder_jobs cannot exceed worker_cpu")
        if self.max_record_bases > self.max_total_bases:
            raise ValueError("max_record_bases cannot exceed max_total_bases")
        if self.max_archive_bytes > self.max_job_disk_bytes:
            raise ValueError("max_archive_bytes cannot exceed max_job_disk_bytes")
        if self.max_launcher_log_bytes > self.max_job_disk_bytes:
            raise ValueError("max_launcher_log_bytes cannot exceed max_job_disk_bytes")
        if not self.integration_command:
            raise ValueError("integration_command cannot be empty")

    def verify_scientific_release(self) -> dict[str, dict[str, str]]:
        """Re-read and verify the exact deployed scientific manifest bytes."""

        if self.association_reference_path is not None:
            path = self.association_reference_path / "reference.json"
            if any(p.is_symlink() for p in (path, *path.parents)) or not path.is_file():
                raise ValueError("Association reference manifest is unavailable")
            with path.open("rb") as stream:
                payload = stream.read(2 * 1024 * 1024 + 1)
            if len(payload) > 2 * 1024 * 1024 or hashlib.sha256(payload).hexdigest() != self.association_reference_sha256:
                raise ValueError("Pinned association reference manifest changed")
        return _load_scientific_release_manifest(
            self.scientific_release_manifest,
            self.scientific_release_manifest_sha256,
        )

    @property
    def scientific_release_identities(self) -> dict[str, dict[str, str]]:
        identities = json.loads(self._scientific_release_identities_json)
        if not isinstance(identities, dict):
            raise RuntimeError("cached scientific release identities are invalid")
        return identities

    def scientific_release_ok(self) -> bool:
        try:
            self.verify_scientific_release()
        except (OSError, ValueError):
            return False
        return True

    @property
    def jobs_root(self) -> Path:
        return self.data_root / "jobs"

    @property
    def submission_admission_free_bytes(self) -> int:
        """Keep the emergency reserve plus every active input materialization free."""

        return self.min_free_bytes + (2 * self.max_request_bytes * self.max_active_jobs)

    @property
    def submission_admission_free_inodes(self) -> int:
        """Reserve the inode envelope of all concurrently materialized inputs."""

        # Job/input/records directories, normalized FASTA and submission lock,
        # plus at most two mode-specific maps (handoff + association groups).
        # Observed comparison needs one array map within this same envelope.
        input_entries = self.max_records + 7
        return self.min_free_inodes + input_entries * self.max_active_jobs

    @property
    def worker_admission_free_bytes(self) -> int:
        """Start work only when a complete job budget fits above the reserve."""

        return self.min_free_bytes + self.max_job_disk_bytes

    @property
    def worker_admission_free_inodes(self) -> int:
        """Start work only when a complete entry budget fits above the reserve."""

        return self.min_free_inodes + self.max_job_entries

    @classmethod
    def from_env(cls) -> "Settings":
        data_root = _safe_root(
            Path(os.getenv("CRISPR_WEB_DATA_ROOT", "/srv/crispr/jobs/web"))
        )
        return cls(
            data_root=data_root,
            database_path=Path(
                os.getenv("CRISPR_WEB_DATABASE", str(data_root / "queue.sqlite3"))
            ),
            catalog_database_path=(
                Path(os.environ["CRISPR_WEB_CATALOG_DATABASE"])
                if os.getenv("CRISPR_WEB_CATALOG_DATABASE") else None
            ),
            catalog_search_index_path=(
                Path(os.environ["CRISPR_WEB_CATALOG_SEARCH_INDEX"])
                if os.getenv("CRISPR_WEB_CATALOG_SEARCH_INDEX") else None
            ),
            runner_config=Path(
                os.getenv("CRISPR_WEB_RUNNER_CONFIG", "/etc/crispr-tools/runners.json")
            ),
            integration_command=_command(
                os.getenv(
                    "CRISPR_WEB_INTEGRATION_COMMAND",
                    '["/opt/crispr-workers/integration/bin/crispr-tools"]',
                )
            ),
            cors_origins=_origins(os.getenv("CRISPR_WEB_CORS_ORIGINS", "")),
            release_id=_sha256_identifier(
                "CRISPR_WEB_RELEASE_ID", os.getenv("CRISPR_WEB_RELEASE_ID")
            ),
            backend_wheel_sha256=_sha256_identifier(
                "CRISPR_WEB_BACKEND_WHEEL_SHA256",
                os.getenv("CRISPR_WEB_BACKEND_WHEEL_SHA256"),
            ),
            scientific_release_manifest=Path(
                os.getenv(
                    "CRISPR_RELEASE_MANIFEST",
                    "/etc/crispr-tools/release-manifest.json",
                )
            ),
            scientific_release_manifest_sha256=_sha256_identifier(
                "CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256",
                os.getenv("CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256"),
            ),
            trusted_proxy_networks=_networks(
                os.getenv("CRISPR_WEB_TRUSTED_PROXY_CIDRS", "")
            ),
            service_version=os.getenv("CRISPR_WEB_RELEASE", __version__),
            max_request_bytes=_integer(
                "CRISPR_WEB_MAX_REQUEST_BYTES", 4_500_000, 1_024, 50_000_000
            ),
            max_total_bases=_integer(
                "CRISPR_WEB_MAX_TOTAL_BASES", 2_000_000, 100, 20_000_000
            ),
            max_record_bases=_integer(
                "CRISPR_WEB_MAX_RECORD_BASES", 2_000_000, 100, 20_000_000
            ),
            max_records=_integer("CRISPR_WEB_MAX_RECORDS", 20, 1, 100),
            max_header_characters=_integer(
                "CRISPR_WEB_MAX_HEADER_CHARACTERS", 200, 20, 1_000
            ),
            max_queued_jobs=_integer("CRISPR_WEB_MAX_QUEUED_JOBS", 8, 1, 100),
            max_active_jobs=_integer("CRISPR_WEB_MAX_ACTIVE_JOBS", 8, 1, 100),
            max_active_jobs_per_client=_integer(
                "CRISPR_WEB_MAX_ACTIVE_JOBS_PER_CLIENT", 2, 1, 100
            ),
            readiness_cache_seconds=_integer(
                "CRISPR_WEB_READINESS_CACHE_SECONDS", 5, 1, 60
            ),
            max_job_attempts=_integer("CRISPR_WEB_MAX_JOB_ATTEMPTS", 2, 1, 10),
            rate_limit_submissions=_integer(
                "CRISPR_WEB_RATE_LIMIT_SUBMISSIONS", 4, 1, 1_000
            ),
            rate_limit_window_seconds=_integer(
                "CRISPR_WEB_RATE_LIMIT_WINDOW_SECONDS", 3_600, 1, 86_400
            ),
            retention_seconds=_integer(
                "CRISPR_WEB_RETENTION_SECONDS", 86_400, 300, 2_592_000
            ),
            worker_cpu=_integer("CRISPR_WEB_WORKER_CPU", 8, 1, 256),
            identify_folder_jobs=_integer("CRISPR_WEB_IDENTIFY_FOLDER_JOBS", 1, 1, 256),
            stage_timeout_seconds=_integer(
                "CRISPR_WEB_STAGE_TIMEOUT_SECONDS", 7_200, 30, 172_800
            ),
            job_timeout_seconds=_integer(
                "CRISPR_WEB_JOB_TIMEOUT_SECONDS", 21_900, 60, 604_800
            ),
            worker_poll_seconds=_integer("CRISPR_WEB_WORKER_POLL_SECONDS", 2, 1, 60),
            worker_lease_seconds=_integer(
                "CRISPR_WEB_WORKER_LEASE_SECONDS", 30, 10, 600
            ),
            worker_stale_seconds=_integer(
                "CRISPR_WEB_WORKER_STALE_SECONDS", 90, 20, 3_600
            ),
            disable_bias_corrections=_boolean(
                "CRISPR_WEB_DISABLE_BIAS_CORRECTIONS", False
            ),
            max_summary_arrays=_integer(
                "CRISPR_WEB_MAX_SUMMARY_ARRAYS", 1_000, 1, 10_000
            ),
            max_summary_rows=_integer("CRISPR_WEB_MAX_SUMMARY_ROWS", 1_000, 1, 10_000),
            max_bundle_files=_integer("CRISPR_WEB_MAX_BUNDLE_FILES", 5_000, 10, 10_000),
            max_artifact_file_bytes=_integer(
                "CRISPR_WEB_MAX_ARTIFACT_FILE_BYTES",
                256_000_000,
                1_024,
                2_000_000_000,
            ),
            max_archive_bytes=_integer(
                "CRISPR_WEB_MAX_ARCHIVE_BYTES",
                512_000_000,
                1_024,
                4_000_000_000,
            ),
            max_job_disk_bytes=_integer(
                "CRISPR_WEB_MAX_JOB_DISK_BYTES",
                2_000_000_000,
                10_000_000,
                100_000_000_000,
            ),
            max_job_entries=_integer(
                "CRISPR_WEB_MAX_JOB_ENTRIES", 50_000, 100, 1_000_000
            ),
            max_launcher_log_bytes=_integer(
                "CRISPR_WEB_MAX_LAUNCHER_LOG_BYTES",
                20_000_000,
                1_024,
                1_000_000_000,
            ),
            min_free_bytes=_integer(
                "CRISPR_WEB_MIN_FREE_BYTES", 2_000_000_000, 0, 100_000_000_000
            ),
            min_free_inodes=_integer(
                "CRISPR_WEB_MIN_FREE_INODES", 100_000, 0, 1_000_000_000
            ),
            job_disk_check_seconds=_integer(
                "CRISPR_WEB_JOB_DISK_CHECK_SECONDS", 1, 1, 60
            ),
            global_rate_limit_submissions=_integer(
                "CRISPR_WEB_GLOBAL_RATE_LIMIT_SUBMISSIONS", 30, 1, 100_000
            ),
            global_rate_limit_window_seconds=_integer(
                "CRISPR_WEB_GLOBAL_RATE_LIMIT_WINDOW_SECONDS", 60, 1, 86_400
            ),
            spacer_database_path=Path(os.environ["CRISPR_WEB_SPACER_DATABASE"]) if os.getenv("CRISPR_WEB_SPACER_DATABASE") else None,
            association_reference_path=Path(os.environ["CRISPR_WEB_ASSOCIATION_REFERENCE"]) if os.getenv("CRISPR_WEB_ASSOCIATION_REFERENCE") else None,
            association_reference_sha256=os.getenv("CRISPR_WEB_ASSOCIATION_REFERENCE_SHA256") or None,
            repeat_reference_path=Path(os.environ["CRISPR_WEB_REPEAT_REFERENCE"]) if os.getenv("CRISPR_WEB_REPEAT_REFERENCE") else None,
            repeat_reference_sha256=os.getenv("CRISPR_WEB_REPEAT_REFERENCE_SHA256") or None,
            viral_reference_path=Path(os.environ["CRISPR_WEB_VIRAL_REFERENCE"]) if os.getenv("CRISPR_WEB_VIRAL_REFERENCE") else None,
            viral_reference_sha256=os.getenv("CRISPR_WEB_VIRAL_REFERENCE_SHA256") or None,
            spacer_database_sha256=os.getenv("CRISPR_WEB_SPACER_DATABASE_SHA256") or None,
            hsts=_boolean("CRISPR_WEB_HSTS", True),
        )


def origin_strings(origins: Iterable[str]) -> list[str]:
    """Return an ordinary list for FastAPI/Starlette configuration."""

    return list(origins)
