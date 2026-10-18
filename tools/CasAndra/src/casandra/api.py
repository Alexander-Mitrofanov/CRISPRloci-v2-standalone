"""Stable Python facade for CLI, workers, and future web services."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from casandra.architecture_runtime import (
    ArchitectureModel,
    load_architecture_model,
)
from casandra.cassette_input import classify_protein_cassette
from casandra.model import ProteinModel, normalize_bundle_config, validate_threads, verify_bundle
from casandra.pipeline import predict_genome
from casandra.prediction import predict_proteins as predict_protein_fasta
from casandra.protein_annotation import annotate_protein_fasta

DEFAULT_MODEL_BUNDLE = Path(__file__).resolve().parent / "models"


@dataclass(frozen=True)
class GenomeInput:
    genome_fasta: Path | None = None
    genbank: Path | None = None
    gff3: Path | None = None


@dataclass(frozen=True)
class AnalysisOptions:
    gene_mode: str = "auto"
    translation_table: int = 11
    threads: int = 6


@dataclass(frozen=True)
class RunResult:
    output_dir: Path
    summary: dict[str, object]


class Analyzer:
    """A verified, preloaded CasAndra model suitable for a long-lived worker."""

    def __init__(
        self,
        bundle_dir: Path,
        manifest: dict[str, Any],
        config: dict[str, Any],
        protein_model: ProteinModel,
        architecture_model: ArchitectureModel,
    ) -> None:
        self.bundle_dir = Path(bundle_dir)
        self.manifest = manifest
        self.config = config
        self.protein_model = protein_model
        self.architecture_model = architecture_model

    @classmethod
    def from_bundle(
        cls, bundle_dir: Path = DEFAULT_MODEL_BUNDLE, *, verify: bool = True
    ) -> Analyzer:
        bundle_dir = Path(bundle_dir)
        if verify:
            manifest, config = verify_bundle(bundle_dir)
        else:
            import json

            manifest = json.loads((bundle_dir / "manifest.json").read_text(encoding="utf-8"))
            config = normalize_bundle_config(json.loads((bundle_dir / "config.json").read_text(encoding="utf-8")))
        protein_model = ProteinModel.from_directory(
            bundle_dir / config["protein_model_directory"], verify=False
        )
        if "specificity_runtime_contract" in config and protein_model.specificity is None:
            raise ValueError("This bundle requires a loaded specificity overlay")
        architecture_model = load_architecture_model(
            bundle_dir / config["cassette_architecture_directory"]
        )
        return cls(bundle_dir, manifest, config, protein_model, architecture_model)

    def inspect(self) -> dict[str, object]:
        result = {
            "schema_version": self.manifest.get("schema_version"),
            "bundle_id": self.config["bundle_id"],
            "bundle_role": self.config.get("bundle_role", self.manifest.get("role")),
            "protein_model_id": self.protein_model.model_id,
            "profile_hmms": len(self.protein_model.hmms),
            "cpu_only": bool(self.config.get("cpu_only")),
            "offline_inference": bool(self.config.get("offline_inference")),
            "integrity": "verified",
        }
        if self.protein_model.specificity is not None:
            overlay = self.protein_model.specificity
            result['specificity_overlay'] = {
                'id': overlay.manifest['id'], 'mode': overlay.manifest['mode'],
                'search_spaces': overlay.manifest['search_spaces'],
                'additional_banks': [{'policy_id':b.card['id'],'profile_hmms':len(b.hmms),
                    'Z':b.card['auxiliary_Z'],'manifest_sha256':b.manifest_sha256} for b in overlay.extensions],
                'deployment_status': 'local_candidate_not_live'}
        return result

    def predict_proteins(
        self, input_fasta: Path, output_jsonl: Path, *, threads: int = 6
    ) -> dict[str, object]:
        return predict_protein_fasta(
            self.protein_model,
            Path(input_fasta),
            Path(output_jsonl),
            threads=validate_threads(threads),
        )

    def annotate_proteins(
        self, input_fasta: Path, output_dir: Path, *, threads: int = 6
    ) -> RunResult:
        """Annotate supplied proteins into an atomic, provenance-bound result directory."""

        summary = annotate_protein_fasta(
            self.bundle_dir,
            Path(input_fasta),
            Path(output_dir),
            protein_model=self.protein_model,
            bundle_manifest=self.manifest,
            bundle_config=self.config,
            threads=validate_threads(threads),
        )
        return RunResult(Path(output_dir), summary)

    def classify_cassette(
        self, input_fasta: Path, output_dir: Path, *, threads: int = 6
    ) -> RunResult:
        """Classify one ordered set of supplied proteins as a Cas cassette."""

        summary = classify_protein_cassette(
            self.bundle_dir,
            Path(input_fasta),
            Path(output_dir),
            protein_model=self.protein_model,
            architecture_model=self.architecture_model,
            bundle_manifest=self.manifest,
            bundle_config=self.config,
            threads=validate_threads(threads),
        )
        return RunResult(Path(output_dir), summary)

    def analyze_genome(
        self,
        genome_input: GenomeInput,
        output_dir: Path,
        options: AnalysisOptions | None = None,
    ) -> RunResult:
        options = options or AnalysisOptions()
        summary = predict_genome(
            self.bundle_dir,
            Path(output_dir),
            genome_fasta=genome_input.genome_fasta,
            genbank=genome_input.genbank,
            gff3=genome_input.gff3,
            gene_mode=options.gene_mode,
            translation_table=options.translation_table,
            threads=validate_threads(options.threads),
            protein_model=self.protein_model,
            architecture_model=self.architecture_model,
            verified_bundle=(self.manifest, self.config),
        )
        return RunResult(Path(output_dir), summary)
