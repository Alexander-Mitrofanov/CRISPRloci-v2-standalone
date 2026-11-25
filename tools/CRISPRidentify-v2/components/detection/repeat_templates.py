"""Strict, immutable supplied repeat panels; deliberately standard-library-only."""

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path

POLICY = "supplied-literal-discovery-v1"
NORMALIZATION = "acgt-min-reverse-complement-v1"
PANEL_OPTIONS = ("--repeat-template-fasta", "--repeat_template_fasta")


class TemplateDiscoveryError(ValueError):
    """A declared template inventory could not be completely and safely evaluated."""


@dataclass(frozen=True)
class TemplatePolicy:
    raw_bytes: int = 65536
    records: int = 256
    canonical_motifs: int = 32
    oriented_patterns: int = 64
    source_length: int = 50_000_000
    base_pattern_positions: int = 1_000_000_000
    hits_per_pattern: int = 8192
    total_hits: int = 65536
    tasks: int = 4096
    hits_per_task: int = 256
    window_span: int = 50000

    def __post_init__(self):
        if any(type(v) is not int or v < 1 for v in asdict(self).values()):
            raise TemplateDiscoveryError("template limits must be positive integers")


DEFAULT_POLICY = TemplatePolicy()


def policy_for_engine(engine):
    # Pinned regex is much more expensive per source position than the native
    # mask scan. These are operational bounds, never acceptance thresholds.
    if engine == "full":
        return replace(DEFAULT_POLICY, source_length=1_000_000,
                       base_pattern_positions=2_000_000)
    if engine != "native":
        raise TemplateDiscoveryError("unsupported template engine: " + str(engine))
    return DEFAULT_POLICY


def validate_template_source(source, panel, policy):
    check_limit("source_length", len(source), policy)
    check_limit("oriented_patterns", len(panel.patterns), policy)
    check_limit("base_pattern_positions", len(source) * len(panel.patterns), policy)


def check_limit(name, count, policy=DEFAULT_POLICY):
    if count > getattr(policy, name):
        raise TemplateDiscoveryError(
            "{}: {} limit exceeded ({} > {})".format(POLICY, name, count, getattr(policy, name)))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode("ascii")).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@dataclass(frozen=True)
class TemplateAlias:
    id: str
    sequence: str
    canonical: str
    orientation: str


@dataclass(frozen=True)
class RepeatTemplatePanel:
    raw_bytes: bytes
    raw_sha256: str
    normalized_sha256: str
    canonical_sequences: tuple
    patterns: tuple
    aliases: tuple


def parse_repeat_template_panel(raw, policy=DEFAULT_POLICY):
    if type(raw) is not bytes:
        raise TemplateDiscoveryError("panel must be a captured byte snapshot")
    check_limit("raw_bytes", len(raw), policy)
    if any(c < 32 and c not in (9, 10, 13) or c == 127 for c in raw):
        raise TemplateDiscoveryError("repeat panel contains invalid control bytes")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise TemplateDiscoveryError("repeat panel must be ASCII") from error
    records, names, current, chunks = [], set(), None, []

    def finish():
        if current is None:
            return
        sequence = "".join(chunks).upper()
        if not 15 <= len(sequence) <= 63 or set(sequence) - set("ACGT"):
            raise TemplateDiscoveryError("template {} requires 15..63 literal A/C/G/T bases".format(current))
        canonical = min(sequence, reverse_complement(sequence))
        records.append(TemplateAlias(current, sequence, canonical,
                                     "+" if sequence == canonical else "-"))
        check_limit("records", len(records), policy)

    for line in lines:
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header = line[1:].split()
            if not header or any(ord(c) < 33 or ord(c) > 126 for c in header[0]):
                raise TemplateDiscoveryError("template header needs a nonempty ASCII identifier")
            current = header[0]
            if current in names:
                raise TemplateDiscoveryError("duplicate template ID: " + current)
            names.add(current)
            chunks = []
        else:
            if current is None:
                raise TemplateDiscoveryError("template sequence precedes its header")
            if any(c not in "ACGTacgt" for c in line):
                raise TemplateDiscoveryError("template sequence contains whitespace or non-ACGT bases")
            chunks.append(line)
    finish()
    if not records:
        raise TemplateDiscoveryError("repeat panel is empty")
    canonical = tuple(sorted({r.canonical for r in records}))
    patterns = tuple(sorted({s for c in canonical for s in (c, reverse_complement(c))}))
    check_limit("canonical_motifs", len(canonical), policy)
    check_limit("oriented_patterns", len(patterns), policy)
    return RepeatTemplatePanel(raw, hashlib.sha256(raw).hexdigest(),
                               digest([NORMALIZATION, canonical]), canonical, patterns,
                               tuple(sorted(records, key=lambda r: r.id)))


def load_repeat_template_panel(path, policy=DEFAULT_POLICY):
    try:
        with open(path, "rb") as handle:
            raw = handle.read(policy.raw_bytes + 1)
    except OSError as error:
        raise TemplateDiscoveryError("cannot read repeat panel: {}".format(error)) from error
    return parse_repeat_template_panel(raw, policy)


def panel_manifest(panel):
    return {"normalization": NORMALIZATION, "raw_sha256": panel.raw_sha256,
            "normalized_sha256": panel.normalized_sha256,
            "raw_bytes": len(panel.raw_bytes),
            "canonical_sequences": list(panel.canonical_sequences),
            "patterns": list(panel.patterns), "aliases": [asdict(a) for a in panel.aliases]}


def panel_configuration(panel, engine):
    return {"policy": POLICY, "engine": engine, "limits": asdict(policy_for_engine(engine)),
            "panel": panel_manifest(panel)}


def validate_template_workflow(preset, accepted_only, engine):
    if preset not in ("balanced", "full") or accepted_only or engine not in ("native", "full"):
        raise TemplateDiscoveryError(
            "repeat templates require balanced/full, accepted_only=False and fuzzy_search_mode=native/full")


def snapshot_repeat_template_panel(panel, destination):
    destination = Path(destination).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_bytes() != panel.raw_bytes:
            raise TemplateDiscoveryError("repeat panel snapshot content differs")
        return str(destination)
    fd, temporary = tempfile.mkstemp(prefix=".repeat-panel-", dir=str(destination.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(panel.raw_bytes)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return str(destination)
