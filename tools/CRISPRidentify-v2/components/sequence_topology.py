"""Strict, source-bound topology evidence for sequence records.

Topology is never inferred from a FASTA header, assembly status, component
type, or sequence content.  A source is circular or linear only after a pinned
local INSDC flatfile has been replayed against that exact source.  Every other
source is explicitly ``unknown`` and therefore cannot authorize a circular
join.

Manifest coordinates and source hashes refer to the complete uppercase source
sequence.  Public provenance retains a canonical manifest-relative artifact
locator so reports are portable; a separate private absolute path is used only
for I/O.  Frozen flatfiles are deliberately opened only when evidence for the
corresponding source is requested.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import date
import hashlib
import json
from pathlib import Path
import re
from typing import Optional, Tuple

from components.reconstruction_validator import IUPAC_DNA_ALPHABET


MANIFEST_SCHEMA_NAME = "CRISPRidentify-sequence-topology-manifest"
MANIFEST_SCHEMA_VERSION = "1.0.0"
EVIDENCE_SCHEMA_NAME = "CRISPRidentify-sequence-topology-evidence"
EVIDENCE_SCHEMA_VERSION = "1.0.0"
EXTENSION_KEY = "sequence_topology"

TOPOLOGY_UNKNOWN = "unknown"
TOPOLOGY_LINEAR = "linear"
TOPOLOGY_CIRCULAR = "circular"
TOPOLOGIES = frozenset((
    TOPOLOGY_UNKNOWN,
    TOPOLOGY_LINEAR,
    TOPOLOGY_CIRCULAR,
))
KNOWN_TOPOLOGIES = frozenset((TOPOLOGY_LINEAR, TOPOLOGY_CIRCULAR))

METHOD_FROZEN_INSDC_FLATFILE = "frozen_insdc_flatfile"
METHOD_UNKNOWN_DEFAULT = "unknown_default"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_VERSION_PATTERN = re.compile(r"^VERSION\s+(\S+)(?:\s+.*)?$")
_LOCUS_PATTERN = re.compile(r"^LOCUS\s+")
_LOCATOR_SEGMENT_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class SequenceTopologyError(ValueError):
    """Base class for topology evidence that cannot be trusted."""


class TopologyManifestError(SequenceTopologyError):
    """A topology manifest is malformed, duplicated, conflicting, or stale."""


class TopologyEvidenceError(SequenceTopologyError):
    """A frozen evidence artifact cannot support its topology assertion."""


class SourceTopologyBindingError(SequenceTopologyError):
    """Evidence is being used with a different source record."""


def _sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def _canonical_sha256(value):
    encoded = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return _sha256_bytes(encoded)


def _require_token(value, label, error_type=SequenceTopologyError):
    if (
            not isinstance(value, str)
            or not value
            or value.strip() != value
            or any(character.isspace() for character in value)):
        raise error_type("{} must be one non-empty whitespace-free string".format(
            label
        ))
    return value


def _require_sha256(value, label, error_type=SequenceTopologyError):
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise error_type(
            "{} must be a lowercase hexadecimal SHA-256".format(label)
        )
    return value


def _require_date(value, label, error_type=SequenceTopologyError):
    if not isinstance(value, str):
        raise error_type("{} must be an ISO calendar date".format(label))
    try:
        date.fromisoformat(value)
    except ValueError as error:
        raise error_type("{} must be an ISO calendar date".format(label)) from error
    return value


def _require_manifest_relative_locator(
        value, label, error_type=SequenceTopologyError):
    if not isinstance(value, str) or not value:
        raise error_type(
            "{} must be a canonical manifest-relative path".format(label)
        )
    segments = value.split("/")
    if (
            "\\" in value
            or any(segment in ("", ".", "..") for segment in segments)
            or any(
                _LOCATOR_SEGMENT_PATTERN.fullmatch(segment) is None
                for segment in segments
            )):
        raise error_type(
            "{} must be a canonical manifest-relative path".format(label)
        )
    return value


def _is_plain_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _normalize_source_sequence(sequence):
    if not isinstance(sequence, str) or not sequence:
        raise SourceTopologyBindingError(
            "source sequence must be one non-empty DNA string"
        )
    try:
        sequence.encode("ascii")
    except UnicodeEncodeError as error:
        raise SourceTopologyBindingError(
            "source sequence must contain ASCII DNA symbols"
        ) from error
    normalized = sequence.upper()
    if set(normalized).difference(IUPAC_DNA_ALPHABET):
        raise SourceTopologyBindingError(
            "source sequence contains non-IUPAC DNA symbols"
        )
    return normalized


@dataclass(frozen=True)
class SourceTopologyBinding:
    """Identity of one complete source after uppercase normalization."""

    source_id: str
    source_length: int
    source_sha256: str

    def __post_init__(self):
        _require_token(self.source_id, "source_id", SourceTopologyBindingError)
        if not _is_plain_int(self.source_length) or self.source_length <= 0:
            raise SourceTopologyBindingError(
                "source_length must be a positive integer"
            )
        _require_sha256(
            self.source_sha256,
            "source_sha256",
            SourceTopologyBindingError,
        )

    @classmethod
    def from_source(cls, source_id, source_sequence):
        source_id = _require_token(
            source_id,
            "source_id",
            SourceTopologyBindingError,
        )
        normalized = _normalize_source_sequence(source_sequence)
        return cls(
            source_id=source_id,
            source_length=len(normalized),
            source_sha256=_sha256_bytes(normalized.encode("ascii")),
        )

    def as_dict(self):
        return {
            "source_id": self.source_id,
            "source_length": self.source_length,
            "source_sha256": self.source_sha256,
            "hash_normalization": "uppercase_ascii_iupac_without_reformatting",
        }


@dataclass(frozen=True)
class TopologyProvenance:
    """Explicit provenance fields retained in every report extension."""

    method: str
    accession: str
    locator: Optional[str]
    date: Optional[str]
    evidence_sha256: Optional[str]

    def __post_init__(self):
        if self.method not in (
                METHOD_FROZEN_INSDC_FLATFILE,
                METHOD_UNKNOWN_DEFAULT):
            raise SequenceTopologyError("unsupported topology evidence method")
        _require_token(self.accession, "provenance accession")
        if self.method == METHOD_FROZEN_INSDC_FLATFILE:
            _require_manifest_relative_locator(
                self.locator,
                "frozen evidence locator",
            )
            _require_date(self.date, "provenance date")
            _require_sha256(self.evidence_sha256, "evidence_sha256")
        elif any(value is not None for value in (
                self.locator,
                self.date,
                self.evidence_sha256)):
            raise SequenceTopologyError(
                "unknown-default provenance cannot claim external evidence"
            )

    def as_dict(self):
        return {
            "method": self.method,
            "accession": self.accession,
            "locator": self.locator,
            "date": self.date,
            "evidence_sha256": self.evidence_sha256,
        }


@dataclass(frozen=True)
class SequenceTopologyEvidence:
    """Immutable topology decision bound to exactly one source record."""

    binding: SourceTopologyBinding
    topology: str
    provenance: TopologyProvenance
    evidence_complete: bool
    evidence_ambiguous: bool

    def __post_init__(self):
        if not isinstance(self.binding, SourceTopologyBinding):
            raise SequenceTopologyError("binding has the wrong type")
        if self.topology not in TOPOLOGIES:
            raise SequenceTopologyError(
                "topology must be circular, linear, or unknown"
            )
        if not isinstance(self.provenance, TopologyProvenance):
            raise SequenceTopologyError("provenance has the wrong type")
        if not isinstance(self.evidence_complete, bool):
            raise SequenceTopologyError("evidence_complete must be a boolean")
        if not isinstance(self.evidence_ambiguous, bool):
            raise SequenceTopologyError("evidence_ambiguous must be a boolean")
        if self.provenance.accession != self.binding.source_id:
            raise SequenceTopologyError(
                "provenance accession must equal the exact source_id"
            )

        if self.topology == TOPOLOGY_UNKNOWN:
            if (
                    self.provenance.method != METHOD_UNKNOWN_DEFAULT
                    or self.evidence_complete
                    or not self.evidence_ambiguous):
                raise SequenceTopologyError(
                    "unknown topology must be incomplete, ambiguous, and defaulted"
                )
        elif (
                self.provenance.method != METHOD_FROZEN_INSDC_FLATFILE
                or not self.evidence_complete
                or self.evidence_ambiguous):
            raise SequenceTopologyError(
                "known topology requires complete unambiguous frozen evidence"
            )

    @property
    def source_id(self):
        return self.binding.source_id

    @property
    def source_length(self):
        return self.binding.source_length

    @property
    def source_sha256(self):
        return self.binding.source_sha256

    def assert_matches(self, source_id, source_sequence):
        supplied = SourceTopologyBinding.from_source(source_id, source_sequence)
        if supplied != self.binding:
            raise SourceTopologyBindingError(
                "topology evidence does not match the supplied source record"
            )

    def as_dict(self):
        return {
            "schema_name": EVIDENCE_SCHEMA_NAME,
            "schema_version": EVIDENCE_SCHEMA_VERSION,
            "source": self.binding.as_dict(),
            "topology": self.topology,
            "provenance": self.provenance.as_dict(),
            "evidence_complete": self.evidence_complete,
            "evidence_ambiguous": self.evidence_ambiguous,
            "inference_used": False,
        }

    def as_extension(self):
        """Return a fresh JSON-safe report extension."""

        return {EXTENSION_KEY: self.as_dict()}


@dataclass(frozen=True)
class _ManifestEntry:
    binding: SourceTopologyBinding
    topology: str
    provenance: TopologyProvenance
    artifact_path: str

    def __post_init__(self):
        if (
                not isinstance(self.artifact_path, str)
                or not Path(self.artifact_path).is_absolute()):
            raise TopologyManifestError(
                "internal artifact_path must be an absolute path"
            )


@dataclass(frozen=True)
class SequenceTopologyManifest:
    """Validated manifest metadata with lazy flatfile verification."""

    manifest_path: str
    entries: Tuple[_ManifestEntry, ...]

    def __post_init__(self):
        if not isinstance(self.manifest_path, str) or not self.manifest_path:
            raise TopologyManifestError("manifest_path must be a non-empty string")
        if not isinstance(self.entries, tuple):
            raise TopologyManifestError("manifest entries must be immutable")
        if any(not isinstance(entry, _ManifestEntry) for entry in self.entries):
            raise TopologyManifestError("manifest contains an invalid entry")

    @property
    def source_ids(self):
        return tuple(entry.binding.source_id for entry in self.entries)

    def _entry_for(self, source_id):
        for entry in self.entries:
            if entry.binding.source_id == source_id:
                return entry
        return None

    def validate_source_binding(self, source_id, source_sequence):
        binding = SourceTopologyBinding.from_source(source_id, source_sequence)
        entry = self._entry_for(binding.source_id)
        if entry is not None and entry.binding != binding:
            raise TopologyManifestError(
                "manifest source binding is stale for {!r}".format(source_id)
            )
        return binding

    def evidence_for(self, source_id, source_sequence):
        """Resolve evidence for one source, reading its artifact only now."""

        binding = self.validate_source_binding(source_id, source_sequence)
        entry = self._entry_for(binding.source_id)
        if entry is None:
            return _unknown_evidence(binding)
        return _verify_frozen_insdc_entry(entry, source_sequence)


def _unknown_evidence(binding):
    return SequenceTopologyEvidence(
        binding=binding,
        topology=TOPOLOGY_UNKNOWN,
        provenance=TopologyProvenance(
            method=METHOD_UNKNOWN_DEFAULT,
            accession=binding.source_id,
            locator=None,
            date=None,
            evidence_sha256=None,
        ),
        evidence_complete=False,
        evidence_ambiguous=True,
    )


def unknown_sequence_topology(source_id, source_sequence):
    """Return the mandatory fail-closed default for an unproven source."""

    return _unknown_evidence(
        SourceTopologyBinding.from_source(source_id, source_sequence)
    )


def _exact_keys(value, expected, label):
    if not isinstance(value, dict):
        raise TopologyManifestError("{} must be a JSON object".format(label))
    actual = frozenset(value)
    expected = frozenset(expected)
    if actual != expected:
        missing = sorted(expected.difference(actual))
        extra = sorted(actual.difference(expected))
        raise TopologyManifestError(
            "{} has incorrect keys (missing={}, extra={})".format(
                label,
                missing,
                extra,
            )
        )


def _unique_json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise TopologyManifestError(
                "manifest contains duplicate JSON key {!r}".format(key)
            )
        result[key] = value
    return result


def _reject_json_constant(value):
    raise TopologyManifestError(
        "manifest contains non-finite JSON value {!r}".format(value)
    )


def _read_manifest_json(path):
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise TopologyManifestError(
            "cannot read topology manifest {!s}".format(path)
        ) from error
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise TopologyManifestError("topology manifest is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_json_constant,
        )
    except TopologyManifestError:
        raise
    except (json.JSONDecodeError, RecursionError) as error:
        raise TopologyManifestError("topology manifest is not valid JSON") from error


def _parse_manifest_entry(raw, index, manifest_directory):
    label = "sources[{}]".format(index)
    _exact_keys(
        raw,
        (
            "source_id",
            "source_length",
            "source_sha256",
            "topology",
            "provenance",
        ),
        label,
    )
    source_id = _require_token(
        raw["source_id"],
        label + ".source_id",
        TopologyManifestError,
    )
    source_length = raw["source_length"]
    if not _is_plain_int(source_length) or source_length <= 0:
        raise TopologyManifestError(
            label + ".source_length must be a positive integer"
        )
    source_hash = _require_sha256(
        raw["source_sha256"],
        label + ".source_sha256",
        TopologyManifestError,
    )
    topology = raw["topology"]
    if topology not in KNOWN_TOPOLOGIES:
        raise TopologyManifestError(
            label + ".topology must be explicitly linear or circular"
        )

    raw_provenance = raw["provenance"]
    _exact_keys(
        raw_provenance,
        ("method", "accession", "locator", "date", "evidence_sha256"),
        label + ".provenance",
    )
    if raw_provenance["method"] != METHOD_FROZEN_INSDC_FLATFILE:
        raise TopologyManifestError(
            label + ".provenance.method is unsupported"
        )
    accession = _require_token(
        raw_provenance["accession"],
        label + ".provenance.accession",
        TopologyManifestError,
    )
    if accession != source_id:
        raise TopologyManifestError(
            label + " accession conflicts with its exact source_id"
        )
    locator = _require_manifest_relative_locator(
        raw_provenance["locator"],
        label + ".provenance.locator",
        TopologyManifestError,
    )
    manifest_directory = manifest_directory.resolve(strict=False)
    unresolved_artifact_path = manifest_directory / locator
    artifact_path = unresolved_artifact_path.resolve(strict=False)
    try:
        artifact_path.relative_to(manifest_directory)
    except ValueError as error:
        raise TopologyManifestError(
            label + ".provenance.locator resolves outside the manifest directory"
        ) from error
    if artifact_path != unresolved_artifact_path:
        raise TopologyManifestError(
            label + ".provenance.locator must not resolve through a symlink"
        )

    evidence_date = _require_date(
        raw_provenance["date"],
        label + ".provenance.date",
        TopologyManifestError,
    )
    evidence_hash = _require_sha256(
        raw_provenance["evidence_sha256"],
        label + ".provenance.evidence_sha256",
        TopologyManifestError,
    )
    try:
        provenance = TopologyProvenance(
            method=METHOD_FROZEN_INSDC_FLATFILE,
            accession=accession,
            locator=locator,
            date=evidence_date,
            evidence_sha256=evidence_hash,
        )
        binding = SourceTopologyBinding(
            source_id=source_id,
            source_length=source_length,
            source_sha256=source_hash,
        )
    except SequenceTopologyError as error:
        raise TopologyManifestError(str(error)) from error
    return _ManifestEntry(
        binding=binding,
        topology=topology,
        provenance=provenance,
        artifact_path=str(artifact_path),
    )


def load_sequence_topology_manifest(path, source_records=None):
    """Load strict metadata without opening any frozen evidence artifact.

    ``source_records`` may be a mapping from exact source IDs to unformatted
    sequence strings.  Matching manifest entries are checked for staleness at
    load time; entries not supplied in that mapping remain available for later
    per-source resolution.
    """

    manifest_path = Path(path).expanduser().resolve(strict=False)
    document = _read_manifest_json(manifest_path)
    _exact_keys(
        document,
        ("schema_name", "schema_version", "sources"),
        "manifest",
    )
    if document["schema_name"] != MANIFEST_SCHEMA_NAME:
        raise TopologyManifestError("unsupported topology manifest schema_name")
    if document["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise TopologyManifestError("unsupported topology manifest schema_version")
    raw_sources = document["sources"]
    if not isinstance(raw_sources, list):
        raise TopologyManifestError("manifest.sources must be a JSON array")

    entries = []
    source_ids = set()
    accessions = set()
    locators = set()
    for index, raw_entry in enumerate(raw_sources):
        entry = _parse_manifest_entry(raw_entry, index, manifest_path.parent)
        source_id = entry.binding.source_id
        if source_id in source_ids:
            raise TopologyManifestError(
                "manifest contains duplicate source_id {!r}".format(source_id)
            )
        accession = entry.provenance.accession
        if accession in accessions:
            raise TopologyManifestError(
                "manifest contains conflicting accession {!r}".format(accession)
            )
        locator = entry.artifact_path
        if locator in locators:
            raise TopologyManifestError(
                "one-record evidence locator is assigned to multiple sources"
            )
        source_ids.add(source_id)
        accessions.add(accession)
        locators.add(locator)
        entries.append(entry)

    manifest = SequenceTopologyManifest(
        manifest_path=str(manifest_path),
        entries=tuple(entries),
    )
    if source_records is not None:
        if not isinstance(source_records, Mapping):
            raise TopologyManifestError(
                "source_records must map exact source IDs to sequences"
            )
        for source_id, source_sequence in source_records.items():
            try:
                manifest.validate_source_binding(source_id, source_sequence)
            except SourceTopologyBindingError as error:
                raise TopologyManifestError(str(error)) from error
    return manifest


def _parse_frozen_insdc_flatfile(raw, locator):
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as error:
        raise TopologyEvidenceError(
            "frozen INSDC flatfile is not ASCII: {}".format(locator)
        ) from error
    lines = text.splitlines()
    nonempty = [index for index, line in enumerate(lines) if line.strip()]
    if not nonempty or _LOCUS_PATTERN.match(lines[nonempty[0]]) is None:
        raise TopologyEvidenceError(
            "frozen INSDC flatfile must start with one LOCUS record"
        )
    locus_indexes = [
        index
        for index, line in enumerate(lines)
        if _LOCUS_PATTERN.match(line) is not None
    ]
    origin_indexes = [
        index for index, line in enumerate(lines) if line.strip() == "ORIGIN"
    ]
    terminator_indexes = [
        index for index, line in enumerate(lines) if line.strip() == "//"
    ]
    if len(locus_indexes) != 1 or len(origin_indexes) != 1:
        raise TopologyEvidenceError(
            "frozen INSDC flatfile must contain exactly one LOCUS and ORIGIN"
        )
    if len(terminator_indexes) != 1:
        raise TopologyEvidenceError(
            "frozen INSDC flatfile must contain exactly one record terminator"
        )
    locus_index = locus_indexes[0]
    origin_index = origin_indexes[0]
    terminator_index = terminator_indexes[0]
    if not locus_index < origin_index < terminator_index:
        raise TopologyEvidenceError("frozen INSDC record sections are out of order")
    if any(line.strip() for line in lines[terminator_index + 1:]):
        raise TopologyEvidenceError(
            "frozen INSDC flatfile contains data after its one record"
        )

    locus_tokens = lines[locus_index].split()
    topology_tokens = [
        token.lower()
        for token in locus_tokens
        if token.lower() in KNOWN_TOPOLOGIES
    ]
    if len(topology_tokens) != 1:
        raise TopologyEvidenceError(
            "LOCUS must contain exactly one explicit linear/circular token"
        )
    bp_indexes = [
        index for index, token in enumerate(locus_tokens)
        if token.lower() == "bp"
    ]
    if (
            len(bp_indexes) != 1
            or bp_indexes[0] == 0
            or not locus_tokens[bp_indexes[0] - 1].isdigit()):
        raise TopologyEvidenceError("LOCUS must contain one explicit base count")
    locus_length = int(locus_tokens[bp_indexes[0] - 1])

    version_values = []
    for line in lines[locus_index + 1:origin_index]:
        if not line.startswith("VERSION"):
            continue
        match = _VERSION_PATTERN.fullmatch(line)
        if match is None:
            raise TopologyEvidenceError("VERSION line is malformed")
        version_values.append(match.group(1))
    if len(version_values) != 1:
        raise TopologyEvidenceError(
            "frozen INSDC record must contain exactly one VERSION accession"
        )

    origin_characters = []
    for line in lines[origin_index + 1:terminator_index]:
        for character in line:
            if character.isspace() or character.isdigit():
                continue
            if not character.isalpha():
                raise TopologyEvidenceError(
                    "ORIGIN contains a non-sequence character"
                )
            origin_characters.append(character.upper())
    origin = "".join(origin_characters)
    if not origin:
        raise TopologyEvidenceError("ORIGIN sequence is empty")
    if set(origin).difference(IUPAC_DNA_ALPHABET):
        raise TopologyEvidenceError("ORIGIN contains non-IUPAC DNA symbols")
    if locus_length != len(origin):
        raise TopologyEvidenceError("LOCUS base count disagrees with ORIGIN")
    return topology_tokens[0], version_values[0], origin


def _verify_frozen_insdc_entry(entry, source_sequence):
    locator = Path(entry.artifact_path)
    if locator.resolve(strict=False) != locator:
        raise TopologyEvidenceError(
            "frozen topology evidence path changed or resolves through a symlink"
        )
    try:
        raw = locator.read_bytes()
    except OSError as error:
        raise TopologyEvidenceError(
            "cannot read frozen topology evidence {!s}".format(locator)
        ) from error
    artifact_hash = _sha256_bytes(raw)
    if artifact_hash != entry.provenance.evidence_sha256:
        raise TopologyEvidenceError(
            "frozen topology evidence SHA-256 does not match the manifest"
        )
    topology, version_accession, origin = _parse_frozen_insdc_flatfile(
        raw,
        locator,
    )
    if version_accession != entry.provenance.accession:
        raise TopologyEvidenceError(
            "VERSION accession does not match the manifest/source ID"
        )
    if topology != entry.topology:
        raise TopologyEvidenceError(
            "LOCUS topology conflicts with the manifest assertion"
        )
    source = _normalize_source_sequence(source_sequence)
    if len(origin) != entry.binding.source_length or origin != source:
        raise TopologyEvidenceError(
            "normalized ORIGIN does not equal the complete source record"
        )
    if _sha256_bytes(origin.encode("ascii")) != entry.binding.source_sha256:
        raise TopologyEvidenceError(
            "normalized ORIGIN hash does not match the source binding"
        )
    return SequenceTopologyEvidence(
        binding=entry.binding,
        topology=entry.topology,
        provenance=entry.provenance,
        evidence_complete=True,
        evidence_ambiguous=False,
    )


def to_circular_topology_evidence(
        evidence,
        source_id,
        source_sequence):
    """Adapt validated evidence for :mod:`circular_array_join`.

    The source binding is checked again at the trust boundary.  The unknown
    adapter uses a deterministic internal decision digest solely to satisfy the
    target data type; it remains incomplete and ambiguous and therefore cannot
    authorize a join.
    """

    if not isinstance(evidence, SequenceTopologyEvidence):
        raise TypeError("evidence has the wrong type")
    evidence.assert_matches(source_id, source_sequence)
    from components.circular_array_join import CircularTopologyEvidence

    if evidence.topology == TOPOLOGY_UNKNOWN:
        provenance = {
            "method": METHOD_UNKNOWN_DEFAULT,
            "source": evidence.binding.as_dict(),
            "reason": "no_verified_external_topology_evidence",
        }
        evidence_method = METHOD_UNKNOWN_DEFAULT
        evidence_date = "0001-01-01"
        evidence_sha256 = _canonical_sha256(provenance)
        accession = evidence.source_id
    else:
        evidence_method = evidence.provenance.method
        evidence_date = evidence.provenance.date
        evidence_sha256 = evidence.provenance.evidence_sha256
        accession = evidence.provenance.accession
    return CircularTopologyEvidence(
        topology=evidence.topology,
        accession=accession,
        source_sha256=evidence.source_sha256,
        evidence_method=evidence_method,
        evidence_date=evidence_date,
        evidence_sha256=evidence_sha256,
        evidence_complete=evidence.evidence_complete,
        evidence_ambiguous=evidence.evidence_ambiguous,
    )


# Short alias for callers that already namespace the operation by module.
load_topology_manifest = load_sequence_topology_manifest
