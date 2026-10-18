"""Strict subprocess adapter for the bundled CasAndra Cas annotator.

CasAndra intentionally runs behind a process boundary.  CRISPRidentify only
accepts its frozen schema-v4 output after validating the installed model,
output manifest, input identity, source coordinates, and cross-file counts.
No legacy Cas backend is selected when this boundary fails.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from itertools import product
from pathlib import Path, PurePosixPath
import sys
import tempfile

from components.bounded_subprocess import (
    BoundedSubprocessTimeoutError,
    BoundedSubprocessTruncationError,
    run_bounded,
)
from components.external_tools import ExternalToolError


CASANDRA_VERSION = "0.1.0.dev0+crispridentify.3"
CASANDRA_SCHEMA_VERSION = 4
CASANDRA_BUNDLE_ID = "casandra-cas-only-cpu-bundle-v4-fast-final"
CASANDRA_BUNDLE_ROLE = "deployment_refit"
CASANDRA_MODEL_ID = "casandra-production-profile-hmm-type-v-v4-fast"
CASANDRA_MODEL_MANIFEST_SHA256 = (
    "8d722e68a720d95b3267ec9c83ff3ec8be8bd61e1e2adee9ce7a670e35a70432"
)
CASANDRA_PROFILE_METADATA_SHA256 = (
    "72e8b0d62d2c369f6a61db62d12cd554940fa61bb9961755efa546d136b211d1"
)
CASANDRA_SOURCE_TREE_SHA256 = (
    "97b895b349cddf9cf73d0d31e1b7d1669c013304ed9bacd77c4c2dc202727796"
)
TYPE_VII_EXPERIMENTAL_REASON = (
    "CasAndra v4 Type VII calls lack independent whole-genome cassette "
    "validation"
)
_CAS_CLASSES = frozenset(("1", "2"))
_CAS_TYPES = frozenset(("I", "II", "III", "IV", "V", "VI", "VII"))
_CAS_TYPE_TO_CLASS = {
    "I": "1", "III": "1", "IV": "1", "VII": "1",
    "II": "2", "V": "2", "VI": "2",
}
_CAS_SUBTYPES = frozenset((
    "I-A", "I-B", "I-C", "I-D", "I-E", "I-F1", "I-F2", "I-G",
    "II-A", "II-B", "II-C", "II-D",
    "III-A", "III-B", "III-C", "III-D", "III-E", "III-F", "III-G",
    "III-H", "III-I",
    "IV-A", "IV-B", "IV-C",
    "V-A", "V-B", "V-C", "V-D", "V-E", "V-F", "V-G", "V-H",
    "V-I", "V-J", "V-K", "V-L", "V-M", "V-N", "V-O", "V-P", "V-Q",
    "VI-A", "VI-B", "VI-C", "VI-D", "VI-E", "VI-F",
))
_CASSETTE_MAX_NON_CAS_GENES = 5
_CASSETTE_MAX_GAP_NT = 5000
_CASSETTE_MIN_SUBTYPE_VOTE = 0.5
_CASANDRA_STDOUT_LIMIT = 1024 * 1024
_CASANDRA_STDERR_LIMIT = 64 * 1024
_BASES = "TCAG"
_AMINO_ACIDS_TABLE_11 = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
_CODON_TABLE_11 = {
    first + second + third: amino_acid
    for amino_acid, (first, second, third) in zip(
        _AMINO_ACIDS_TABLE_11,
        (
            (first, second, third)
            for first in _BASES
            for second in _BASES
            for third in _BASES
        ),
    )
}
_IUPAC_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


_REQUIRED_OUTPUTS = frozenset((
    "run.json",
    "proteins.jsonl",
    "protein_predictions.jsonl",
    "cassettes.jsonl",
    "rejected_cassette_candidates.jsonl",
    "extraction/audit.json",
    "extraction/manifest.json",
))
_IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)
_AMBIGUOUS_TO_N = str.maketrans({base: "N" for base in "RYSWKMBDHV"})


@dataclass(frozen=True)
class CasAndraInstallation:
    """Validated metadata for the CasAndra installation in this interpreter."""

    package_root: Path
    version: str
    bundle_id: str
    bundle_role: str
    model_id: str
    model_manifest_sha256: str
    profile_metadata: dict


class _DuplicateJsonKey(ValueError):
    pass


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise _DuplicateJsonKey("duplicate JSON key {!r}".format(key))
        value[key] = item
    return value


def _reject_json_constant(value):
    raise ValueError("non-finite JSON number {}".format(value))


def _loads_json(value, description):
    try:
        return json.loads(
            value,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_json_constant,
        )
    except (json.JSONDecodeError, _DuplicateJsonKey, ValueError) as error:
        raise ExternalToolError(
            "CasAndra returned invalid {}: {}".format(description, error)
        ) from error


def _read_json(path, description):
    try:
        value = Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise ExternalToolError(
            "could not read CasAndra {}: {}".format(description, error)
        ) from error
    result = _loads_json(value, description)
    if not isinstance(result, dict):
        raise ExternalToolError(
            "CasAndra {} must contain a JSON object".format(description)
        )
    return result


def _read_jsonl(path, description):
    rows = []
    try:
        with Path(path).open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                if not raw_line.strip():
                    raise ExternalToolError(
                        "CasAndra {} contains a blank line at {}".format(
                            description, line_number
                        )
                    )
                row = _loads_json(
                    raw_line,
                    "{} line {}".format(description, line_number),
                )
                if not isinstance(row, dict):
                    raise ExternalToolError(
                        "CasAndra {} line {} must be an object".format(
                            description, line_number
                        )
                    )
                rows.append(row)
    except ExternalToolError:
        raise
    except (OSError, UnicodeError) as error:
        raise ExternalToolError(
            "could not read CasAndra {}: {}".format(description, error)
        ) from error
    return rows


def _sha256(path):
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
    except OSError as error:
        raise ExternalToolError(
            "could not hash required CasAndra file {}: {}".format(path, error)
        ) from error
    return digest.hexdigest()


def _required_string(mapping, key, description):
    value = mapping.get(key) if isinstance(mapping, dict) else None
    if not isinstance(value, str) or not value:
        raise ExternalToolError(
            "CasAndra {} has invalid {}".format(description, key)
        )
    return value


def _required_integer(mapping, key, description, minimum=0):
    value = mapping.get(key) if isinstance(mapping, dict) else None
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ExternalToolError(
            "CasAndra {} has invalid {}".format(description, key)
        )
    return value


def _required_number(mapping, key, description):
    value = mapping.get(key) if isinstance(mapping, dict) else None
    if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)):
        raise ExternalToolError(
            "CasAndra {} has invalid {}".format(description, key)
        )
    return float(value)


def _require_value(mapping, key, expected, description):
    observed = mapping.get(key) if isinstance(mapping, dict) else None
    if observed != expected:
        raise ExternalToolError(
            "CasAndra {} {} is {!r}, expected {!r}".format(
                description, key, observed, expected
            )
        )


def _validate_artifact_manifest(root, manifest):
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise ExternalToolError("CasAndra model manifest has no artifacts")
    for relative_name, record in sorted(artifacts.items()):
        if not isinstance(relative_name, str):
            raise ExternalToolError("CasAndra model manifest path is invalid")
        relative = PurePosixPath(relative_name)
        if (
                relative.is_absolute()
                or ".." in relative.parts
                or "\\" in relative_name):
            raise ExternalToolError(
                "CasAndra model manifest path is unsafe: {!r}".format(
                    relative_name
                )
            )
        path = root.joinpath(*relative.parts)
        if path.is_symlink() or not path.is_file():
            raise ExternalToolError(
                "CasAndra model artifact is missing: {}".format(path)
            )
        if not isinstance(record, dict):
            raise ExternalToolError(
                "CasAndra model artifact record is malformed: {}".format(
                    relative_name
                )
            )
        size = _required_integer(record, "size", "model artifact", minimum=1)
        digest = _required_string(record, "sha256", "model artifact")
        if path.stat().st_size != size or _sha256(path) != digest:
            raise ExternalToolError(
                "CasAndra model artifact failed integrity validation: {}"
                .format(path)
            )


def _casandra_source_aggregate(package_root):
    lines = []
    for path in Path(package_root).glob("*.py"):
        if path.is_symlink() or not path.is_file():
            raise ExternalToolError(
                "CasAndra source module is not a regular file: {}".format(path)
            )
        lines.append("{}  {}\n".format(_sha256(path), path.name))
    if not lines:
        raise ExternalToolError("the installed CasAndra package has no modules")
    return hashlib.sha256(
        "".join(sorted(lines)).encode("utf-8")
    ).hexdigest()


def _validate_casandra_source(package_root):
    observed = _casandra_source_aggregate(package_root)
    if observed != CASANDRA_SOURCE_TREE_SHA256:
        raise ExternalToolError(
            "the installed CasAndra source does not match the pinned "
            "CRISPRidentify snapshot"
        )


def inspect_casandra_installation():
    """Validate the exact CasAndra package and deployment model in use."""

    try:
        distribution = importlib.metadata.distribution("CasAndra")
    except importlib.metadata.PackageNotFoundError as error:
        raise ExternalToolError(
            "CasAndra is not installed in {}".format(sys.executable)
        ) from error
    version = distribution.version
    if version != CASANDRA_VERSION:
        raise ExternalToolError(
            "CasAndra version {} is installed; this release requires {}"
            .format(version, CASANDRA_VERSION)
        )

    try:
        spec = importlib.util.find_spec("casandra")
    except (ImportError, AttributeError, ValueError) as error:
        raise ExternalToolError(
            "the installed CasAndra package cannot be located"
        ) from error
    if spec is None or spec.origin is None:
        raise ExternalToolError("the installed CasAndra package cannot be located")
    package_root = Path(spec.origin).resolve().parent
    metadata_package_root = Path(
        distribution.locate_file("casandra")
    ).resolve()
    if package_root != metadata_package_root:
        raise ExternalToolError(
            "CasAndra import path disagrees with installed distribution metadata"
        )
    if not (package_root / "__main__.py").is_file():
        raise ExternalToolError(
            "the installed CasAndra package has no python -m entry point"
        )

    _validate_casandra_source(package_root)
    model_root = package_root / "models"
    model_manifest_path = model_root / "manifest.json"
    if (
            not model_manifest_path.is_file()
            or _sha256(model_manifest_path) != CASANDRA_MODEL_MANIFEST_SHA256):
        raise ExternalToolError(
            "the installed CasAndra model manifest does not match the pinned "
            "CRISPRidentify deployment bundle"
        )
    model_manifest = _read_json(model_manifest_path, "model manifest")
    _require_value(
        model_manifest, "schema_version", CASANDRA_SCHEMA_VERSION,
        "model manifest",
    )
    _require_value(
        model_manifest, "bundle_id", CASANDRA_BUNDLE_ID, "model manifest"
    )
    _validate_artifact_manifest(model_root, model_manifest)

    config = _read_json(model_root / "config.json", "model config")
    _require_value(config, "schema_version", CASANDRA_SCHEMA_VERSION, "model config")
    _require_value(config, "bundle_id", CASANDRA_BUNDLE_ID, "model config")
    _require_value(config, "bundle_role", CASANDRA_BUNDLE_ROLE, "model config")
    _require_value(config, "cpu_only", True, "model config")
    _require_value(config, "offline_inference", True, "model config")
    _require_value(config, "crispr_array_prediction", False, "model config")

    protein_manifest = _read_json(
        model_root / "protein" / "manifest.json", "protein model manifest"
    )
    _require_value(
        protein_manifest, "schema_version", CASANDRA_SCHEMA_VERSION,
        "protein model manifest",
    )
    _require_value(
        protein_manifest, "model_id", CASANDRA_MODEL_ID,
        "protein model manifest",
    )

    profile_path = model_root / "protein" / "profile_metadata.json"
    if _sha256(profile_path) != CASANDRA_PROFILE_METADATA_SHA256:
        raise ExternalToolError(
            "the installed CasAndra profile metadata failed integrity validation"
        )
    profiles = _read_json(profile_path, "profile metadata")
    if not profiles:
        raise ExternalToolError("CasAndra profile metadata is empty")
    return CasAndraInstallation(
        package_root=package_root,
        version=version,
        bundle_id=CASANDRA_BUNDLE_ID,
        bundle_role=CASANDRA_BUNDLE_ROLE,
        model_id=CASANDRA_MODEL_ID,
        model_manifest_sha256=CASANDRA_MODEL_MANIFEST_SHA256,
        profile_metadata=profiles,
    )


def casandra_installation_error():
    """Return an actionable preflight error, or ``None`` when ready."""

    try:
        inspect_casandra_installation()
    except (ExternalToolError, OSError) as error:
        return (
            "--cas True requires the pinned CasAndra {} backend: {}; run "
            "scripts/setup_local.sh and launch with "
            "scripts/run_crispridentify.sh".format(CASANDRA_VERSION, error)
        )
    return None


def _validate_input_fasta(path, expected_source_id, expected_dna):
    header = None
    sequence = []
    record_count = 0
    try:
        with Path(path).open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    record_count += 1
                    if record_count > 1:
                        raise ExternalToolError(
                            "CasAndra input must contain exactly one FASTA record"
                        )
                    header = line[1:].strip()
                    continue
                if header is None:
                    raise ExternalToolError(
                        "CasAndra input has sequence before its FASTA header "
                        "at line {}".format(line_number)
                    )
                if any(character.isspace() for character in line):
                    raise ExternalToolError(
                        "CasAndra input contains embedded sequence whitespace"
                    )
                sequence.append(line.upper())
    except ExternalToolError:
        raise
    except (OSError, UnicodeError) as error:
        raise ExternalToolError(
            "could not read CasAndra input FASTA: {}".format(error)
        ) from error
    if record_count != 1 or not header or not sequence:
        raise ExternalToolError(
            "CasAndra input must contain exactly one non-empty FASTA record"
        )
    source_id = header.split()[0]
    if source_id != expected_source_id:
        raise ExternalToolError(
            "CasAndra input source ID {!r} disagrees with {!r}".format(
                source_id, expected_source_id
            )
        )
    observed_dna = "".join(sequence)
    if observed_dna != expected_dna:
        raise ExternalToolError(
            "CasAndra input sequence disagrees with the validated pipeline source"
        )


def _safe_output_path(output_root, relative_name):
    if not isinstance(relative_name, str):
        raise ExternalToolError("CasAndra output manifest path is invalid")
    relative = PurePosixPath(relative_name)
    if (
            relative.is_absolute()
            or ".." in relative.parts
            or "\\" in relative_name):
        raise ExternalToolError(
            "CasAndra output manifest path is unsafe: {!r}".format(relative_name)
        )
    return output_root.joinpath(*relative.parts)


def _validate_output_manifest(output_root, installation):
    manifest_path = output_root / "manifest.json"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise ExternalToolError("CasAndra did not produce manifest.json")
    manifest = _read_json(manifest_path, "output manifest")
    _require_value(
        manifest, "schema_version", CASANDRA_SCHEMA_VERSION, "output manifest"
    )
    _require_value(
        manifest, "bundle_manifest_sha256",
        installation.model_manifest_sha256, "output manifest",
    )
    files = manifest.get("files")
    if not isinstance(files, dict):
        raise ExternalToolError("CasAndra output manifest files are malformed")
    missing = sorted(_REQUIRED_OUTPUTS.difference(files))
    if missing:
        raise ExternalToolError(
            "CasAndra output manifest is missing: {}".format(", ".join(missing))
        )

    manifest_names = set()
    for relative_name, record in sorted(files.items()):
        path = _safe_output_path(output_root, relative_name)
        manifest_names.add(relative_name)
        if path.is_symlink() or not path.is_file():
            raise ExternalToolError(
                "CasAndra manifested output is missing: {}".format(relative_name)
            )
        if not isinstance(record, dict):
            raise ExternalToolError(
                "CasAndra output manifest record is malformed: {}".format(
                    relative_name
                )
            )
        size = _required_integer(record, "size", "output manifest", minimum=0)
        digest = _required_string(record, "sha256", "output manifest")
        if path.stat().st_size != size or _sha256(path) != digest:
            raise ExternalToolError(
                "CasAndra output failed manifest validation: {}".format(
                    relative_name
                )
            )

    actual_names = set()
    for path in output_root.rglob("*"):
        if path.is_symlink():
            raise ExternalToolError(
                "CasAndra output contains a symbolic link: {}".format(path)
            )
        if path.is_file():
            actual_names.add(path.relative_to(output_root).as_posix())
    expected_names = manifest_names.union(("manifest.json",))
    if actual_names != expected_names:
        raise ExternalToolError(
            "CasAndra output files disagree with its manifest"
        )
    return manifest


def _reverse_complement(sequence):
    return sequence.translate(_IUPAC_COMPLEMENT)[::-1]


def _coordinate_nucleotide(source_dna, start, end, strand):
    sequence = source_dna[start - 1:end]
    return _reverse_complement(sequence) if strand == "-" else sequence


def _translate_table_11(nucleotide):
    if len(nucleotide) % 3:
        raise ExternalToolError(
            "CasAndra coding interval length is not divisible by three"
        )
    amino_acids = []
    for offset in range(0, len(nucleotide), 3):
        codon = nucleotide[offset:offset + 3].upper()
        try:
            expansions = product(*(_IUPAC_BASES[base] for base in codon))
        except KeyError as error:
            raise ExternalToolError(
                "CasAndra coding interval contains invalid base {!r}".format(
                    error.args[0]
                )
            ) from error
        translations = {
            _CODON_TABLE_11["".join(expansion)]
            for expansion in expansions
        }
        amino_acids.append(
            next(iter(translations)) if len(translations) == 1 else "X"
        )
    protein = "".join(amino_acids)
    return protein[:-1] if protein.endswith("*") else protein


def _profile_annotation(profile, reported_families, installation):
    metadata = installation.profile_metadata.get(profile)
    if not isinstance(metadata, dict):
        raise ExternalToolError(
            "CasAndra prediction references unknown profile {!r}".format(profile)
        )
    if metadata.get("intended_use") != "primary_search":
        raise ExternalToolError(
            "CasAndra positive prediction references non-primary profile {!r}"
            .format(profile)
        )
    metadata_families = metadata.get("families")
    if not isinstance(metadata_families, list):
        raise ExternalToolError(
            "CasAndra profile {!r} families are malformed".format(profile)
        )
    expected_families = []
    for value in metadata_families:
        if not isinstance(value, str):
            raise ExternalToolError(
                "CasAndra profile {!r} family is malformed".format(profile)
            )
        family = value.strip()
        if family and family not in expected_families:
            expected_families.append(family)
    if reported_families != expected_families:
        raise ExternalToolError(
            "CasAndra prediction families disagree with profile {!r}"
            .format(profile)
        )
    return (
        reported_families[0] if reported_families else profile,
        list(reported_families),
    )


def _classification_values(
        classification, description, *, enforce_hierarchy=True):
    if not isinstance(classification, dict):
        raise ExternalToolError(
            "CasAndra {} classification is malformed".format(description)
        )
    values = []
    for key in ("class", "type", "subtype"):
        value = classification.get(key)
        if value is not None and (not isinstance(value, str) or not value):
            raise ExternalToolError(
                "CasAndra {} classification {} is invalid".format(
                    description, key
                )
            )
        if key == "class" and value is not None and value not in _CAS_CLASSES:
            raise ExternalToolError(
                "CasAndra {} class is outside the frozen enum".format(
                    description
                )
            )
        if key == "type" and value is not None and value not in _CAS_TYPES:
            raise ExternalToolError(
                "CasAndra {} type is outside the frozen enum".format(
                    description
                )
            )
        if key == "subtype" and value is not None:
            if value not in _CAS_SUBTYPES:
                raise ExternalToolError(
                    "CasAndra {} subtype is outside the frozen enum".format(
                        description
                    )
                )
        values.append(value)
    class_name, cas_type, subtype = values
    if enforce_hierarchy and (
            class_name is not None
            and cas_type is not None
            and class_name != _CAS_TYPE_TO_CLASS[cas_type]):
        raise ExternalToolError(
            "CasAndra {} class and type are inconsistent".format(
                description
            )
        )
    if enforce_hierarchy and subtype is not None and (
            cas_type is None or not subtype.startswith(cas_type + "-")):
        raise ExternalToolError(
            "CasAndra {} type and subtype are inconsistent".format(
                description
            )
        )
    return tuple(values)


def _validate_prediction(prediction, protein_id, protein_length, installation):
    description = "prediction for {}".format(protein_id)
    _require_value(prediction, "sequence_id", protein_id, description)
    _require_value(prediction, "sequence_length", protein_length, description)
    is_cas = prediction.get("is_cas")
    if not isinstance(is_cas, bool):
        raise ExternalToolError("CasAndra {} has invalid is_cas".format(description))
    positive = _required_number(prediction, "positive_profile_score", description)
    negative = _required_number(
        prediction, "hard_negative_profile_score", description
    )
    margin = _required_number(prediction, "score_margin", description)
    if not math.isclose(margin, positive - negative, rel_tol=1e-10, abs_tol=1e-8):
        raise ExternalToolError(
            "CasAndra {} score margin is inconsistent".format(description)
        )
    classification = prediction.get("classification")
    class_name, cas_type, subtype = _classification_values(
        classification, description
    )
    profile = prediction.get("best_positive_profile")
    if profile is not None and (not isinstance(profile, str) or not profile):
        raise ExternalToolError(
            "CasAndra {} profile is invalid".format(description)
        )
    if is_cas and profile is None:
        raise ExternalToolError(
            "CasAndra positive {} has no positive profile".format(description)
        )
    profile_families = prediction.get("profile_families", [])
    if (
            not isinstance(profile_families, list)
            or len(profile_families) != len(set(profile_families))
            or any(
                not isinstance(family, str) or not family
                for family in profile_families
            )):
        raise ExternalToolError(
            "CasAndra {} profile families are malformed".format(description)
        )
    if not is_cas and profile_families:
        raise ExternalToolError(
            "CasAndra negative {} has profile families".format(description)
        )
    if not is_cas and any(
            value is not None for value in (class_name, cas_type, subtype)):
        raise ExternalToolError(
            "CasAndra negative {} has a Cas classification".format(description)
        )

    evidence = prediction.get("evidence")
    if not isinstance(evidence, dict):
        raise ExternalToolError(
            "CasAndra {} evidence is malformed".format(description)
        )
    _require_value(evidence, "model_id", installation.model_id, description)
    _require_value(evidence, "positive_profile_required", True, description)
    _required_integer(evidence, "profile_hits", description, minimum=0)
    threshold = _required_number(
        evidence, "decision_threshold", description
    )
    report_evalue = _required_number(evidence, "report_evalue", description)
    _require_semantic_equal(threshold, 0.0, "protein decision threshold")
    _require_semantic_equal(
        report_evalue, 1000.0, "protein profile reporting E-value"
    )
    _require_value(
        evidence, "type_v_subtype_masked", False, description
    )
    raw_classification = evidence.get("raw_classification")
    if (
            not isinstance(raw_classification, dict)
            or set(raw_classification) != {"class", "type", "subtype"}):
        raise ExternalToolError(
            "CasAndra raw {} classification is malformed".format(
                description
            )
        )
    raw_class_name, raw_type, raw_subtype = _classification_values(
        raw_classification,
        "raw {}".format(description),
        enforce_hierarchy=False,
    )
    expected_classification = (
        _CAS_TYPE_TO_CLASS[raw_type]
        if raw_type is not None else raw_class_name,
        raw_type,
        (
            raw_subtype
            if (
                raw_type is not None
                and raw_subtype is not None
                and raw_subtype.startswith(raw_type + "-")
            )
            else None
        ),
    )
    observed_classification = (class_name, cas_type, subtype)
    if observed_classification != expected_classification:
        raise ExternalToolError(
            "CasAndra {} classification disagrees with its hierarchy-"
            "normalized raw classification".format(description)
        )
    expected_conflict_masked = (
        (raw_class_name, raw_type, raw_subtype) != expected_classification
    )
    conflict_masked = evidence.get("hierarchy_conflict_masked")
    if not isinstance(conflict_masked, bool):
        raise ExternalToolError(
            "CasAndra {} has invalid hierarchy_conflict_masked".format(
                description
            )
        )
    if conflict_masked != expected_conflict_masked:
        raise ExternalToolError(
            "CasAndra {} hierarchy_conflict_masked disagrees with its raw "
            "classification".format(description)
        )
    if not is_cas and any(
            value is not None
            for value in (raw_class_name, raw_type, raw_subtype)):
        raise ExternalToolError(
            "CasAndra negative {} has a raw Cas classification".format(
                description
            )
        )
    expected_is_cas = profile is not None and margin >= threshold
    if is_cas != expected_is_cas:
        raise ExternalToolError(
            "CasAndra {} is_cas disagrees with its frozen protein decision"
            .format(description)
        )
    return {
        "classification_evidence": {
            "raw_classification": dict(raw_classification),
            "hierarchy_conflict_masked": conflict_masked,
        },
        "is_cas": is_cas,
        "profile": profile,
        "profile_families": profile_families,
        "score": positive,
        "score_kind": "CasAndra positive-profile score",
        "classification": {
            "class": class_name,
            "type": cas_type,
            "subtype": subtype,
        },
    }


def _validate_proteins(
        proteins, predictions, source_id, source_dna, source_sha256,
        gene_mode, translation_table, installation):
    predictions_by_id = {}
    for prediction in predictions:
        protein_id = _required_string(prediction, "sequence_id", "prediction")
        if protein_id in predictions_by_id:
            raise ExternalToolError(
                "CasAndra returned duplicate prediction {}".format(protein_id)
            )
        predictions_by_id[protein_id] = prediction

    proteins_by_id = {}
    normalized = {}
    previous = None
    caller_versions = set()
    selected_modes = set()
    for protein in proteins:
        protein_id = _required_string(protein, "protein_id", "protein")
        if protein_id in proteins_by_id:
            raise ExternalToolError(
                "CasAndra returned duplicate protein {}".format(protein_id)
            )
        proteins_by_id[protein_id] = protein
        _require_value(protein, "contig_id", source_id, "protein")
        _require_value(
            protein, "genome_source_sha256", source_sha256, "protein"
        )
        start = _required_integer(protein, "start_1based", "protein", minimum=1)
        end = _required_integer(
            protein, "end_1based_inclusive", "protein", minimum=start
        )
        if end > len(source_dna):
            raise ExternalToolError(
                "CasAndra protein {} lies outside source {}".format(
                    protein_id, source_id
                )
            )
        strand = protein.get("strand")
        if strand not in ("+", "-"):
            raise ExternalToolError(
                "CasAndra protein {} has invalid strand".format(protein_id)
            )
        ordering = (start, end, protein_id)
        if previous is not None and ordering < previous:
            raise ExternalToolError("CasAndra proteins are not source ordered")
        previous = ordering

        nucleotide = _required_string(
            protein, "nucleotide_sequence", "protein {}".format(protein_id)
        ).upper()
        expected_nucleotide = _coordinate_nucleotide(
            source_dna, start, end, strand
        )
        if (
                nucleotide != expected_nucleotide
                and nucleotide != expected_nucleotide.translate(_AMBIGUOUS_TO_N)):
            raise ExternalToolError(
                "CasAndra coordinates do not reconstruct protein {} from source {}"
                .format(protein_id, source_id)
            )
        _require_value(
            protein, "nucleotide_sha256",
            hashlib.sha256(nucleotide.encode("ascii")).hexdigest(),
            "protein {}".format(protein_id),
        )
        protein_sequence = _required_string(
            protein, "protein_sequence", "protein {}".format(protein_id)
        )
        protein_length = _required_integer(
            protein, "protein_length", "protein {}".format(protein_id),
            minimum=1,
        )
        if len(protein_sequence) != protein_length:
            raise ExternalToolError(
                "CasAndra protein {} length is inconsistent".format(protein_id)
            )
        try:
            protein_bytes = protein_sequence.encode("ascii")
        except UnicodeEncodeError as error:
            raise ExternalToolError(
                "CasAndra protein {} sequence is not ASCII".format(protein_id)
            ) from error
        _require_value(
            protein, "protein_sha256",
            hashlib.sha256(protein_bytes).hexdigest(),
            "protein {}".format(protein_id),
        )
        translated = _translate_table_11(nucleotide)
        if translated != protein_sequence:
            raise ExternalToolError(
                "CasAndra reported protein {} does not match the table-11 "
                "translation of its source interval".format(protein_id)
            )
        _require_value(
            protein, "translation_table", translation_table,
            "protein {}".format(protein_id),
        )
        caller = protein.get("caller")
        if not isinstance(caller, dict):
            raise ExternalToolError(
                "CasAndra protein {} caller is malformed".format(protein_id)
            )
        _require_value(caller, "program", "Pyrodigal", "gene caller")
        _require_value(caller, "requested_mode", gene_mode, "gene caller")
        _require_value(
            caller, "translation_table", translation_table, "gene caller"
        )
        caller_versions.add(_required_string(caller, "version", "gene caller"))
        selected_mode = _required_string(caller, "mode", "gene caller")
        if selected_mode not in ("single", "meta"):
            raise ExternalToolError("CasAndra selected an invalid gene mode")
        selected_modes.add(selected_mode)

        prediction = predictions_by_id.get(protein_id)
        if prediction is None or protein.get("prediction") != prediction:
            raise ExternalToolError(
                "CasAndra protein and prediction records disagree for {}"
                .format(protein_id)
            )
        prediction_data = _validate_prediction(
            prediction, protein_id, protein_length, installation
        )
        if prediction_data["is_cas"]:
            annotation, profile_families = _profile_annotation(
                prediction_data["profile"],
                prediction_data["profile_families"],
                installation,
            )
        else:
            annotation, profile_families = "unknown", []
        normalized[protein_id] = {
            "annotation": annotation,
            "classification": prediction_data["classification"],
            "classification_evidence": prediction_data[
                "classification_evidence"
            ],
            "end": end,
            "is_cas": prediction_data["is_cas"],
            "profile": prediction_data["profile"],
            "profile_families": profile_families,
            "protein_id": protein_id,
            "protein_length": protein_length,
            "protein_sequence": protein_sequence,
            "score": prediction_data["score"],
            "score_kind": prediction_data["score_kind"],
            "source_id": source_id,
            "start": start,
            "strand": strand,
            "translation_table": translation_table,
        }

    if set(predictions_by_id) != set(proteins_by_id):
        raise ExternalToolError(
            "CasAndra protein and prediction identifiers disagree"
        )
    if len(caller_versions) > 1 or len(selected_modes) > 1:
        raise ExternalToolError(
            "CasAndra returned inconsistent gene-caller provenance"
        )
    return normalized, caller_versions, selected_modes


def _weighted_label(genes, head):
    votes = {}
    supporters = {}
    for gene in genes:
        label = gene["classification"].get(head)
        if label is None:
            continue
        weight = max(1.0, float(gene["score"]))
        votes[label] = votes.get(label, 0.0) + weight
        supporters.setdefault(label, []).append(gene["protein_id"])
    if not votes:
        return None, 0.0, {}, {}
    ordered = sorted(votes, key=lambda label: (-votes[label], label))
    winner = ordered[0]
    total = sum(votes.values())
    return (
        winner,
        votes[winner] / total,
        dict(sorted(votes.items())),
        supporters,
    )


def _expected_direct_classification(cas_genes):
    class_name, class_fraction, class_votes, _ = _weighted_label(
        cas_genes, "class"
    )
    cas_type, type_fraction, type_votes, _ = _weighted_label(
        cas_genes, "type"
    )
    raw_subtype, subtype_fraction, subtype_votes, supporters = (
        _weighted_label(cas_genes, "subtype")
    )
    subtype = (
        raw_subtype
        if subtype_fraction >= _CASSETTE_MIN_SUBTYPE_VOTE else None
    )
    if (
            subtype is not None
            and cas_type is not None
            and not subtype.startswith(cas_type + "-")):
        subtype = None
    return {
        "abstained_subtype": (
            raw_subtype is not None and subtype is None
        ),
        "class": class_name,
        "confidence": {
            "class_vote_fraction": class_fraction,
            "subtype_vote_fraction": subtype_fraction,
            "type_vote_fraction": type_fraction,
        },
        "raw_subtype_winner": raw_subtype,
        "subtype": subtype,
        "subtype_supporting_proteins": (
            supporters.get(raw_subtype, []) if raw_subtype else []
        ),
        "type": cas_type,
        "votes": {
            "class": class_votes,
            "subtype": subtype_votes,
            "type": type_votes,
        },
    }


def _semantic_equal(observed, expected):
    if (
            isinstance(observed, (int, float))
            and not isinstance(observed, bool)
            and isinstance(expected, (int, float))
            and not isinstance(expected, bool)):
        return (
            math.isfinite(observed)
            and math.isfinite(expected)
            and math.isclose(
                float(observed), float(expected),
                rel_tol=1e-12, abs_tol=1e-12,
            )
        )
    if isinstance(expected, dict):
        return (
            isinstance(observed, dict)
            and set(observed) == set(expected)
            and all(
                _semantic_equal(observed[key], value)
                for key, value in expected.items()
            )
        )
    if isinstance(expected, list):
        return (
            isinstance(observed, list)
            and len(observed) == len(expected)
            and all(
                _semantic_equal(item, expected_item)
                for item, expected_item in zip(observed, expected)
            )
        )
    return type(observed) is type(expected) and observed == expected


def _require_semantic_equal(observed, expected, description):
    if not _semantic_equal(observed, expected):
        raise ExternalToolError(
            "CasAndra {} disagrees with recomputed frozen-v4 semantics".format(
                description
            )
        )


def _expected_candidate_claims(genes_by_id):
    ordered = sorted(
        genes_by_id.values(),
        key=lambda gene: (
            gene["start"], gene["end"], gene["protein_id"]
        ),
    )
    index_by_id = {
        gene["protein_id"]: index for index, gene in enumerate(ordered)
    }
    positives = [gene for gene in ordered if gene["is_cas"]]
    groups = []
    for gene in positives:
        if not groups:
            groups.append([gene])
            continue
        previous = groups[-1][-1]
        intervening = (
            index_by_id[gene["protein_id"]]
            - index_by_id[previous["protein_id"]]
            - 1
        )
        gap_nt = gene["start"] - previous["end"] - 1
        if (
                intervening <= _CASSETTE_MAX_NON_CAS_GENES
                and gap_nt <= _CASSETTE_MAX_GAP_NT):
            groups[-1].append(gene)
        else:
            groups.append([gene])

    claims = []
    for ordinal, group in enumerate(groups, 1):
        first_index = index_by_id[group[0]["protein_id"]]
        last_index = index_by_id[group[-1]["protein_id"]]
        context = ordered[first_index:last_index + 1]
        source_id = group[0]["source_id"]
        start = min(gene["start"] for gene in group)
        end = max(gene["end"] for gene in group)
        claims.append({
            "cas_gene_count": len(group),
            "cas_protein_ids": [gene["protein_id"] for gene in group],
            "cassette_id": (
                "casandra|contig={}|cassette={:04d}|loc={}-{}"
                .format(source_id, ordinal, start, end)
            ),
            "classification": _expected_direct_classification(group),
            "construction": {
                "crispr_array_evidence_used": False,
                "max_cas_gap_nt": _CASSETTE_MAX_GAP_NT,
                "max_non_cas_genes": _CASSETTE_MAX_NON_CAS_GENES,
            },
            "context_gene_count": len(context),
            "context_protein_ids": [
                gene["protein_id"] for gene in context
            ],
            "contig_id": source_id,
            "end_1based_inclusive": end,
            "start_1based": start,
        })
    return claims


def _validate_final_route(row, expected_claim):
    direct = expected_claim["classification"]
    final = row.get("final_classification")
    if not isinstance(final, dict):
        raise ExternalToolError(
            "CasAndra cassette final classification is malformed"
        )
    required_fields = {
        "class", "confidence", "direct_profile_result", "evidence",
        "method", "subtype", "type",
    }
    if set(final) != required_fields:
        raise ExternalToolError(
            "CasAndra cassette final classification fields are malformed"
        )
    _require_semantic_equal(
        final["direct_profile_result"],
        direct,
        "final direct-profile result",
    )

    evidence = final.get("evidence")
    if (
            not isinstance(evidence, dict)
            or set(evidence) != {
                "raw_classification", "hierarchy_conflict_masked"}):
        raise ExternalToolError(
            "CasAndra final cassette classification evidence is malformed"
        )
    raw_classification = evidence.get("raw_classification")
    if (
            not isinstance(raw_classification, dict)
            or set(raw_classification) != {"class", "type", "subtype"}):
        raise ExternalToolError(
            "CasAndra raw final cassette classification is malformed"
        )
    raw_class_name, raw_type, raw_subtype = _classification_values(
        raw_classification,
        "raw final cassette",
        enforce_hierarchy=False,
    )
    expected_classification = (
        _CAS_TYPE_TO_CLASS[raw_type]
        if raw_type is not None else raw_class_name,
        raw_type,
        (
            raw_subtype
            if (
                raw_type is not None
                and raw_subtype is not None
                and raw_subtype.startswith(raw_type + "-")
            )
            else None
        ),
    )
    class_name, cas_type, subtype = _classification_values(
        final, "final cassette"
    )
    if (class_name, cas_type, subtype) != expected_classification:
        raise ExternalToolError(
            "CasAndra final cassette classification disagrees with its "
            "hierarchy-normalized raw route"
        )
    expected_conflict_masked = (
        (raw_class_name, raw_type, raw_subtype) != expected_classification
    )
    conflict_masked = evidence.get("hierarchy_conflict_masked")
    if not isinstance(conflict_masked, bool):
        raise ExternalToolError(
            "CasAndra final cassette hierarchy_conflict_masked is invalid"
        )
    if conflict_masked != expected_conflict_masked:
        raise ExternalToolError(
            "CasAndra final cassette hierarchy_conflict_masked disagrees "
            "with its raw classification"
        )
    if raw_class_name != direct["class"] or raw_type != direct["type"]:
        raise ExternalToolError(
            "CasAndra raw final class/type disagrees with direct profile "
            "votes"
        )
    method = _required_string(final, "method", "final cassette")
    confidence = _required_number(final, "confidence", "final cassette")
    if confidence < 0.0 or confidence > 1.0:
        raise ExternalToolError(
            "CasAndra cassette confidence is outside [0, 1]"
        )

    if direct["type"] == "III":
        if method != "type_iii_ordered_profile_architecture_extratrees":
            raise ExternalToolError(
                "CasAndra Type III cassette did not use architecture routing"
            )
        if raw_subtype is None or not raw_subtype.startswith("III-"):
            raise ExternalToolError(
                "CasAndra Type III architecture result has invalid subtype"
            )
    else:
        has_subtype = direct["subtype"] is not None
        expected_method = (
            "direct_profile_aggregation"
            if has_subtype else "direct_profile_type_aggregation"
        )
        expected_confidence_key = (
            "subtype_vote_fraction"
            if has_subtype else "type_vote_fraction"
        )
        if method != expected_method or raw_subtype != direct["subtype"]:
            raise ExternalToolError(
                "CasAndra direct cassette route disagrees with profile votes"
            )
        _require_semantic_equal(
            confidence,
            direct["confidence"][expected_confidence_key],
            "direct-route confidence",
        )
    return final


def _expected_genome_gate(expected_claim, final, genes_by_id):
    if final["subtype"] is None and final["type"] != "VII":
        return {"accepted": False, "rule": "no_supported_subtype"}
    count = expected_claim["cas_gene_count"]
    if final["class"] == "1":
        return {
            "accepted": count >= 2,
            "rule": "class1_requires_two_cas_genes",
        }
    if final["class"] == "2":
        if count >= 2:
            return {"accepted": True, "rule": "class2_multi_gene"}
        protein_id = expected_claim["cas_protein_ids"][0]
        gene = genes_by_id[protein_id]
        accepted = gene["protein_length"] >= 400 and gene["score"] >= 50.0
        return {
            "accepted": accepted,
            "rule": "class2_single_effector_requires_length400_score50",
        }
    return {
        "accepted": count >= 2,
        "rule": "unresolved_class_requires_two_cas_genes",
    }


def _validate_candidate_partition(accepted, rejected, genes_by_id):
    expected_claims = _expected_candidate_claims(genes_by_id)
    expected_by_id = {
        claim["cassette_id"]: claim for claim in expected_claims
    }
    actual_rows = list(accepted) + list(rejected)
    actual_by_id = {}
    for row in actual_rows:
        cassette_id = _required_string(row, "cassette_id", "cassette")
        if cassette_id in actual_by_id:
            raise ExternalToolError(
                "CasAndra returned duplicate cassette {}".format(cassette_id)
            )
        actual_by_id[cassette_id] = row
    if set(actual_by_id) != set(expected_by_id):
        raise ExternalToolError(
            "CasAndra accepted/rejected cassettes do not partition all "
            "recomputed Cas-positive gene groups"
        )

    accepted_ids = {
        _required_string(row, "cassette_id", "accepted cassette")
        for row in accepted
    }
    expected_claim_fields = set(next(iter(expected_claims), {}))
    for cassette_id, expected_claim in expected_by_id.items():
        row = actual_by_id[cassette_id]
        expected_fields = expected_claim_fields.union((
            "final_classification", "genome_evidence_gate",
        ))
        if set(row) != expected_fields:
            raise ExternalToolError(
                "CasAndra cassette fields disagree with frozen schema v4"
            )
        for key, expected in expected_claim.items():
            _require_semantic_equal(
                row.get(key),
                expected,
                "cassette {}".format(key),
            )
        final = _validate_final_route(row, expected_claim)
        expected_gate = _expected_genome_gate(
            expected_claim, final, genes_by_id
        )
        _require_semantic_equal(
            row.get("genome_evidence_gate"),
            expected_gate,
            "genome evidence gate",
        )
        if (cassette_id in accepted_ids) != expected_gate["accepted"]:
            raise ExternalToolError(
                "CasAndra accepted/rejected partition disagrees with its "
                "frozen genome evidence gate"
            )


def _validate_rejected_cassettes(
        rows, source_id, source_length, proteins_by_id, seen_ids):
    for row in rows:
        cassette_id = _required_string(row, "cassette_id", "rejected cassette")
        if cassette_id in seen_ids:
            raise ExternalToolError(
                "CasAndra returned duplicate cassette {}".format(cassette_id)
            )
        seen_ids.add(cassette_id)
        _require_value(row, "contig_id", source_id, "rejected cassette")
        start = _required_integer(
            row, "start_1based", "rejected cassette", minimum=1
        )
        end = _required_integer(
            row, "end_1based_inclusive", "rejected cassette", minimum=start
        )
        if end > source_length:
            raise ExternalToolError(
                "CasAndra rejected cassette lies outside its source"
            )
        context = row.get("context_protein_ids")
        if (
                not isinstance(context, list)
                or len(context) != len(set(context))
                or any(item not in proteins_by_id for item in context)):
            raise ExternalToolError(
                "CasAndra rejected cassette context is malformed"
            )
        gate = row.get("genome_evidence_gate")
        if not isinstance(gate, dict) or gate.get("accepted") is not False:
            raise ExternalToolError(
                "CasAndra rejected cassette lacks rejection evidence"
            )


def _validate_cassettes(rows, source_id, source_length, genes_by_id):
    systems = []
    seen_ids = set()
    previous = None
    for cassette_number, row in enumerate(rows, 1):
        provider_id = _required_string(row, "cassette_id", "cassette")
        if provider_id in seen_ids:
            raise ExternalToolError(
                "CasAndra returned duplicate cassette {}".format(provider_id)
            )
        seen_ids.add(provider_id)
        _require_value(row, "contig_id", source_id, "cassette")
        start = _required_integer(row, "start_1based", "cassette", minimum=1)
        end = _required_integer(
            row, "end_1based_inclusive", "cassette", minimum=start
        )
        if end > source_length:
            raise ExternalToolError("CasAndra cassette lies outside its source")
        ordering = (start, end, provider_id)
        if previous is not None and ordering < previous:
            raise ExternalToolError("CasAndra cassettes are not source ordered")
        previous = ordering

        context_ids = row.get("context_protein_ids")
        cas_ids = row.get("cas_protein_ids")
        if (
                not isinstance(context_ids, list)
                or not context_ids
                or len(context_ids) != len(set(context_ids))
                or any(item not in genes_by_id for item in context_ids)):
            raise ExternalToolError("CasAndra cassette context is malformed")
        if (
                not isinstance(cas_ids, list)
                or not cas_ids
                or len(cas_ids) != len(set(cas_ids))
                or any(item not in context_ids for item in cas_ids)):
            raise ExternalToolError("CasAndra cassette Cas genes are malformed")
        _require_value(row, "context_gene_count", len(context_ids), "cassette")
        _require_value(row, "cas_gene_count", len(cas_ids), "cassette")
        observed_cas = {
            protein_id for protein_id in context_ids
            if genes_by_id[protein_id]["is_cas"]
        }
        if observed_cas != set(cas_ids):
            raise ExternalToolError(
                "CasAndra cassette Cas gene membership is inconsistent"
            )
        context_genes = [genes_by_id[protein_id] for protein_id in context_ids]
        if context_genes != sorted(
                context_genes,
                key=lambda gene: (gene["start"], gene["end"], gene["protein_id"])):
            raise ExternalToolError(
                "CasAndra cassette context genes are not source ordered"
            )
        if (
                start != min(gene["start"] for gene in context_genes)
                or end != max(gene["end"] for gene in context_genes)):
            raise ExternalToolError(
                "CasAndra cassette boundaries disagree with its context genes"
            )

        final = row.get("final_classification")
        class_name, cas_type, subtype = _classification_values(
            final, "final cassette"
        )
        if class_name is None or cas_type is None:
            raise ExternalToolError(
                "CasAndra accepted cassette lacks a supported type"
            )
        if subtype is None and cas_type != "VII":
            raise ExternalToolError(
                "CasAndra accepted non-Type-VII cassette lacks a subtype"
            )
        if subtype is not None and not subtype.startswith(cas_type + "-"):
            raise ExternalToolError(
                "CasAndra cassette type and subtype are inconsistent"
            )
        method = _required_string(final, "method", "final cassette")
        confidence_kinds = {
            "direct_profile_aggregation": (
                "positive_bitscore_weighted_subtype_vote_fraction"
            ),
            "direct_profile_type_aggregation": (
                "positive_bitscore_weighted_type_vote_fraction"
            ),
            "type_iii_ordered_profile_architecture_extratrees": (
                "architecture_extratrees_winner_probability"
            ),
        }
        if method not in confidence_kinds:
            raise ExternalToolError(
                "CasAndra cassette has unsupported classification method {!r}"
                .format(method)
            )
        confidence = _required_number(final, "confidence", "final cassette")
        if confidence < 0.0 or confidence > 1.0:
            raise ExternalToolError(
                "CasAndra cassette confidence is outside [0, 1]"
            )
        gate = row.get("genome_evidence_gate")
        if not isinstance(gate, dict) or gate.get("accepted") is not True:
            raise ExternalToolError(
                "CasAndra accepted cassette lacks acceptance evidence"
            )
        rule = _required_string(gate, "rule", "genome evidence gate")

        experimental = subtype is None and cas_type == "VII"
        experimental_reason = (
            TYPE_VII_EXPERIMENTAL_REASON if experimental else None
        )
        system_genes = []
        for gene in context_genes:
            system_gene = dict(gene)
            system_gene.pop("protein_sequence", None)
            system_genes.append(system_gene)
        systems.append({
            "cassette_id": cassette_number,
            "classification": {
                "class": class_name,
                "confidence": confidence,
                "confidence_kind": confidence_kinds[method],
                "method": method,
                "subtype": subtype,
                "type": cas_type,
            },
            "end": end,
            "evidence": {
                "genome_evidence_gate": rule,
                "experimental_reason": experimental_reason,
                "provider_direct_profile_result": final.get(
                    "direct_profile_result"
                ),
                "provider_final_classification_evidence": {
                    "raw_classification": dict(
                        final["evidence"]["raw_classification"]
                    ),
                    "hierarchy_conflict_masked": final["evidence"][
                        "hierarchy_conflict_masked"],
                },
                "type_vii_experimental": experimental,
            },
            "experimental": experimental,
            "experimental_reason": experimental_reason,
            "genes": system_genes,
            "provider_cassette_id": provider_id,
            "source_id": source_id,
            "start": start,
            "type": "CAS-{}".format(subtype or cas_type),
        })
    return systems, seen_ids


class CasAndraPredictor:
    """Run CasAndra for one validated FASTA source and normalize its output."""

    def __init__(
            self, input_fasta, source_id, source_dna, cpu_count=1,
            gene_mode="auto", translation_table=11, installation=None,
            command_runner=None, timeout_seconds=3600):
        if gene_mode not in ("auto", "single", "meta"):
            raise ValueError("CasAndra gene mode must be auto, single, or meta")
        if translation_table != 11:
            raise ValueError("CasAndra integration currently requires table 11")
        if (
                isinstance(timeout_seconds, bool)
                or not isinstance(timeout_seconds, (int, float))
                or not math.isfinite(timeout_seconds)
                or timeout_seconds <= 0):
            raise ValueError(
                "CasAndra timeout_seconds must be a positive finite number"
            )
        self.input_fasta = Path(input_fasta).resolve()
        self.source_id = source_id
        self.source_dna = source_dna.upper()
        self.cpu_count = max(1, int(cpu_count))
        self.gene_mode = gene_mode
        self.translation_table = translation_table
        self.timeout_seconds = float(timeout_seconds)
        self.installation = installation
        self.command_runner = command_runner or run_bounded

    def run(self):
        _validate_input_fasta(
            self.input_fasta, self.source_id, self.source_dna
        )
        installation = self.installation or inspect_casandra_installation()
        source_sha256 = _sha256(self.input_fasta)

        with tempfile.TemporaryDirectory(
                prefix="crispridentify_casandra_") as temporary_directory:
            output_root = Path(temporary_directory) / "output"
            command = [
                sys.executable,
                "-P",
                "-m",
                "casandra",
                "predict-genome",
                "--genome",
                str(self.input_fasta),
                "--output",
                str(output_root),
                "--gene-mode",
                self.gene_mode,
                "--translation-table",
                str(self.translation_table),
                "--threads",
                str(self.cpu_count),
            ]
            child_environment = os.environ.copy()
            child_environment["PYTHONPYCACHEPREFIX"] = str(
                Path(temporary_directory) / "pycache"
            )
            try:
                completed = self.command_runner(
                    command,
                    env=child_environment,
                    errors="strict",
                    stderr_limit=_CASANDRA_STDERR_LIMIT,
                    stdout_limit=_CASANDRA_STDOUT_LIMIT,
                    timeout=self.timeout_seconds,
                )
            except BoundedSubprocessTimeoutError as error:
                raise ExternalToolError(
                    "CasAndra exceeded its {:.3f} second timeout; no legacy "
                    "Cas backend was used".format(self.timeout_seconds)
                ) from error
            except BoundedSubprocessTruncationError as error:
                raise ExternalToolError(
                    "CasAndra stdout exceeded its {} byte safety limit; no "
                    "legacy Cas backend was used".format(
                        _CASANDRA_STDOUT_LIMIT
                    )
                ) from error
            except UnicodeError as error:
                raise ExternalToolError(
                    "CasAndra output is not strict UTF-8; no legacy Cas "
                    "backend was used"
                ) from error
            except OSError as error:
                raise ExternalToolError(
                    "could not execute CasAndra with {}: {}; no legacy Cas "
                    "backend was used".format(
                        sys.executable, error,
                    )
                ) from error
            if completed.returncode != 0:
                detail = (completed.stderr or "").strip()
                if len(detail) > 2000:
                    detail = detail[-2000:]
                raise ExternalToolError(
                    "CasAndra failed with exit status {}{}; no legacy Cas "
                    "backend was used".format(
                        completed.returncode,
                        ": {}".format(detail) if detail else "",
                    )
                )

            _validate_output_manifest(output_root, installation)
            run = _read_json(output_root / "run.json", "run metadata")
            stdout_run = _loads_json(completed.stdout or "", "CLI summary")
            if stdout_run != run:
                raise ExternalToolError(
                    "CasAndra CLI summary disagrees with run.json"
                )
            self._validate_run(run, source_sha256, installation)

            proteins = _read_jsonl(
                output_root / "proteins.jsonl", "proteins.jsonl"
            )
            predictions = _read_jsonl(
                output_root / "protein_predictions.jsonl",
                "protein_predictions.jsonl",
            )
            cassettes = _read_jsonl(
                output_root / "cassettes.jsonl", "cassettes.jsonl"
            )
            rejected = _read_jsonl(
                output_root / "rejected_cassette_candidates.jsonl",
                "rejected_cassette_candidates.jsonl",
            )
            genes_by_id, caller_versions, selected_modes = _validate_proteins(
                proteins,
                predictions,
                self.source_id,
                self.source_dna,
                source_sha256,
                self.gene_mode,
                self.translation_table,
                installation,
            )
            _validate_candidate_partition(
                cassettes, rejected, genes_by_id
            )
            systems, seen_cassette_ids = _validate_cassettes(
                cassettes,
                self.source_id,
                len(self.source_dna),
                genes_by_id,
            )
            _validate_rejected_cassettes(
                rejected,
                self.source_id,
                len(self.source_dna),
                genes_by_id,
                seen_cassette_ids,
            )
            self._validate_counts(
                run, proteins, genes_by_id, systems, cassettes, rejected
            )

            if caller_versions:
                caller_version = next(iter(caller_versions))
                selected_mode = next(iter(selected_modes))
            else:
                extraction_manifest = _read_json(
                    output_root / "extraction" / "manifest.json",
                    "extraction manifest",
                )
                caller = extraction_manifest.get("caller")
                if not isinstance(caller, dict):
                    raise ExternalToolError(
                        "CasAndra extraction manifest has no gene caller"
                    )
                _require_value(caller, "program", "Pyrodigal", "gene caller")
                _require_value(
                    caller, "requested_mode", self.gene_mode, "gene caller"
                )
                caller_version = _required_string(
                    caller, "version", "gene caller"
                )
                selected_mode = _required_string(caller, "mode", "gene caller")

            positive_genes = []
            for gene in genes_by_id.values():
                if gene["is_cas"]:
                    positive_genes.append(dict(gene))
            hierarchy_conflict_masked_proteins = sum(
                gene["classification_evidence"]["hierarchy_conflict_masked"]
                for gene in positive_genes
            )
            return {
                "coordinate_system": {
                    "end_inclusive": True,
                    "indexing": 1,
                    "reference": "input_forward_strand",
                },
                "gene_count": len(positive_genes),
                "genes": positive_genes,
                "method": {
                    "bundle_id": installation.bundle_id,
                    "bundle_role": installation.bundle_role,
                    "gene_caller": "Pyrodigal {}".format(caller_version),
                    "gene_mode": self.gene_mode,
                    "model_manifest_sha256": (
                        installation.model_manifest_sha256
                    ),
                    "pipeline_schema_version": CASANDRA_SCHEMA_VERSION,
                    "provider": "CasAndra",
                    "provider_version": installation.version,
                    "selected_gene_mode": selected_mode,
                    "timeout_seconds": self.timeout_seconds,
                    "translation_table": self.translation_table,
                },
                "protein_count": len(proteins),
                "run_statistics": {
                    "hierarchy_conflict_masked_cassettes": run[
                        "hierarchy_conflict_masked_cassettes"
                    ],
                    "hierarchy_conflict_masked_proteins": (
                        hierarchy_conflict_masked_proteins
                    ),
                    "raw_cassette_candidates": run[
                        "raw_cassette_candidates"
                    ],
                    "rejected_cassette_candidates": run[
                        "rejected_cassette_candidates"
                    ],
                },
                "status": "computed",
                "system_count": len(systems),
                "systems": systems,
            }

    def _validate_run(self, run, source_sha256, installation):
        _require_value(run, "schema_version", CASANDRA_SCHEMA_VERSION, "run")
        _require_value(run, "program", "CasAndra", "run")
        _require_value(run, "bundle_id", installation.bundle_id, "run")
        _require_value(run, "bundle_role", installation.bundle_role, "run")
        _require_value(run, "cpu_threads", self.cpu_count, "run")
        _require_value(run, "offline_inference", True, "run")
        _require_value(run, "crispr_array_prediction", False, "run")
        input_data = run.get("input")
        if not isinstance(input_data, dict):
            raise ExternalToolError("CasAndra run input provenance is malformed")
        _require_value(input_data, "sha256", source_sha256, "run input")
        _require_value(
            input_data, "mode", "pyrodigal_de_novo", "run input"
        )
        for key in (
                "genes", "cas_proteins", "raw_cassette_candidates",
                "rejected_cassette_candidates",
                "hierarchy_conflict_masked_cassettes", "cas_cassettes",
                "classified_cassettes"):
            _required_integer(run, key, "run", minimum=0)

    @staticmethod
    def _validate_counts(
            run, proteins, genes_by_id, systems, accepted, rejected):
        positive_count = sum(gene["is_cas"] for gene in genes_by_id.values())
        expected = {
            "genes": len(proteins),
            "cas_proteins": positive_count,
            "cas_cassettes": len(systems),
            "classified_cassettes": len(systems),
            "rejected_cassette_candidates": len(rejected),
            "raw_cassette_candidates": len(systems) + len(rejected),
            "hierarchy_conflict_masked_cassettes": sum(
                cassette["final_classification"]["evidence"][
                    "hierarchy_conflict_masked"
                ]
                for cassette in accepted
            ),
        }
        for key, value in expected.items():
            _require_value(run, key, value, "run")
