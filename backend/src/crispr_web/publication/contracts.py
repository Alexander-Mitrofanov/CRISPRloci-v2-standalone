"""Artifact allowlists and publication limits."""

from __future__ import annotations

import re


_INDIVIDUAL_EXACT = frozenset(
    {
        "locus-comparison.json", "locus-alignments.tsv", "locus-proteins.faa",
        "prophage-context.json", "prophage-regions.tsv", "prophage-spacer-matches.tsv",
        "insertion-context.json", "insertion-elements.tsv",
        "defense-context.json", "defense-systems.tsv",
        "rna-family.json", "rna-family-alignment.fasta", "rna-family-consensus.dbn", "rna-family-pairs.tsv",
        "result.json",
        "input-sequences.json",
        "arrays.json",
        "repeats.json",
        "spacers.json",
        "spacer-associations.json",
        "array-comparison.json",
        "array-sharing.tsv",
        "spacer-associations.tsv",
        "repeat-typing.json",
        "repeat-typing.tsv",
        "repeat-mapping.json",
        "repeat-mapping.tsv",
        "viral-matches.json",
        "viral-matches.tsv",
        "protospacer-analysis.json",
        "protospacer-matches.tsv",
        "repeat-analysis.json",
        "repeat-instances.tsv",
        "repeat-analysis.fasta",
        "repeat-structures.txt",
        "leader-contexts.json",
        "leader-contexts.tsv",
        "leader-contexts.fasta",
        "leader-contexts.gff3",
        "leader-contexts.bed",
        "annotations.json",
        "array-evidence.json",
        "cas-genes.tsv",
        "cas-cassettes.tsv",
        "tracr-candidates.tsv",
        "tracr-predictions.fasta",
        "spacerplacer/0_results.csv",
        "spacerplacer/0_protocol_skipped.csv",
        "spacerplacer/summary.txt",
        "crispr_evor/0_results.csv",
        "crispr_evor/0_protocol_skipped.csv",
        "crispr_evor/summary.txt",
        "crispr_evor/crispr_evor_results.json",
    }
)
_IDENTIFY_FILENAMES = frozenset(
    {"report.json", "Summary.csv", "Complete_summary.csv", "Complete_json_summary.json"}
)
_BUNDLE_SUFFIXES = frozenset(
    {".dbn", ".faa", ".bed", ".csv", ".tsv", ".gff3", ".json", ".txt", ".log", ".nwk", ".tree", ".fa", ".fasta"}
)
_DENIED_PARTS = frozenset({"locus-compare", "locus-compare-input", "phigaro", "phigaro-input", "source-spacer", "insertion-input", "isescan", "defense-context", "defense-input", "rna-family", "rna-family-input", "context-input", "work_folder", "__pycache__", "work", "leader", "repeat", "protospacer", "protospacer-input", "viral", "viral-input", "repeat-map", "repeat-map-input", "repeat-type", "repeat-type-input", "association", "association-input", "array-compare", "array-compare-input"})
_MEDIA_TYPES = {
    ".dbn": "text/plain",
    ".faa": "text/plain",
    ".bed": "text/plain",
    ".json": "application/json",
    ".csv": "text/csv",
    ".tsv": "text/tab-separated-values",
    ".gff3": "text/plain",
    ".txt": "text/plain",
    ".log": "text/plain",
    ".nwk": "text/plain",
    ".tree": "text/plain",
    ".fa": "text/plain",
    ".fasta": "text/plain",
    ".zip": "application/zip",
}
_INTERNAL_JSON = frozenset({"integration_result.json", "adapter/manifest.json"})
_STATIC_PRIVATE_ROOTS = (
    "/srv",
    "/opt",
    "/etc",
    "/home",
    "/root",
    "/tmp",
    "/var/tmp",
    "/var/lib",
)
_PRIVATE_MANIFEST_KEYS = frozenset(
    {
        "argv",
        "ccdb_file",
        "command",
        "commands",
        "executable",
        "input",
        "report_path",
        "runner_config",
        "stderr",
        "stdout",
    }
)
_OMIT = object()
_MAX_PUBLIC_JSON_BYTES = 20_000_000
_EMBEDDED_ABSOLUTE_PATH = re.compile(
    r"(?<![A-Za-z0-9:/])/(?:[A-Za-z0-9._-]+/)+[A-Za-z0-9._-]+"
)
