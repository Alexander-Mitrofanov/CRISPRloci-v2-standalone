"""Explicit orchestration with isolated per-record working directories."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import time
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from . import __version__
from .adapters import load_json, run_casandra, run_identify, sha256_file
from .config import DATA, Settings
from .errors import TracrError
from .evidence import (
    candidate_evidence,
    evidence_class,
    merge_model_hits,
    model_rescue_predictions,
    short_rescue_admission,
)
from .models import AntiRepeat, Hit, Sequence
from .process import Runner
from .reporting import _assign_export_ids
from .reporting import atomic_json as atomic_json
from .reporting import write_exports as write_exports
from .search import interactions, make_candidates, scan_models, search_anti_repeats, terminators
from .sequence import input_files, read_fasta, write_fasta

LOG = logging.getLogger(__name__)


def doctor(settings: Settings) -> dict[str, Any]:
    required = {"cmscan": settings.cmscan}
    if settings.mode == "complete":
        required.update(
            {
                "identify": settings.identify_command[0],
                "casandra": settings.casandra_command[0],
                "IntaRNA": settings.intarna,
                "erpin": settings.erpin,
                settings.search_backend: getattr(settings, settings.search_backend),
            }
        )
    search_path = settings.environment.get("PATH", os.environ.get("PATH"))
    executables = {name: shutil.which(value, path=search_path) for name, value in required.items()}
    missing = [
        f"executable: {name} ({required[name]})" for name, path in executables.items() if not path
    ]
    resources = list(settings.model_paths())
    for model in settings.model_paths():
        resources.extend(Path(str(model) + suffix) for suffix in (".i1f", ".i1i", ".i1m", ".i1p"))
    if settings.mode == "complete":
        resources.append(Path(settings.erpin_model))
    missing.extend(f"resource: {path}" for path in resources if not path.is_file())
    manifest = load_json(DATA / "manifest.json")
    for relative, expected in manifest.items():
        resource = DATA / relative
        if not resource.is_file() or sha256_file(resource) != expected:
            missing.append(f"Bundled resource checksum mismatch: {relative}")
    return {
        "ready": not missing,
        "executables": executables,
        "errors": missing,
        "note": (
            "Checks paths and bundled resources; upstream runtimes are verified during execution."
        ),
    }


def command_provenance(
    command: tuple[str, ...], environment: dict[str, str] | None = None
) -> dict[str, Any]:
    """Record checkout revision when a command directly references a repository."""
    result: dict[str, Any] = {"command": command}
    for part in command:
        path = Path(part)
        if path.is_file():
            result.setdefault("files", {})[str(path.resolve())] = sha256_file(path)
            try:
                revision = subprocess.run(
                    ["git", "-C", str(path.resolve().parent), "rev-parse", "HEAD"],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False,
                )
                if revision.returncode == 0:
                    result["git_commit"] = revision.stdout.strip()
            except (OSError, subprocess.TimeoutExpired):
                pass
    # A Python interpreter may live in an unrelated checkout. Its revision is
    # not the revision of a package loaded with `python -m package`.
    if "-m" in command:
        index = command.index("-m")
        if index + 1 < len(command):
            module = command[index + 1]
            result["module"] = module
            if "git_commit" in result:
                result["interpreter_checkout_git_commit"] = result.pop("git_commit")
            search = (environment or {}).get("PYTHONPATH", os.environ.get("PYTHONPATH", ""))
            result["module_resolution"] = "unresolved; no import executed"
            for root in search.split(os.pathsep):
                if not root:
                    continue
                module_path = Path(root).joinpath(*module.split("."))
                choices = [module_path / "__init__.py", module_path.with_suffix(".py")]
                found = next((path for path in choices if path.is_file()), None)
                if found is None:
                    continue
                result["module_source"] = str(found.resolve())
                result["module_source_sha256"] = sha256_file(found)
                result["module_resolution"] = (
                    "static explicit effective PYTHONPATH; no import executed"
                )
                try:
                    revision = subprocess.run(
                        ["git", "-C", str(found.resolve().parent), "rev-parse", "HEAD"],
                        capture_output=True,
                        text=True,
                        timeout=10,
                        check=False,
                    )
                    if revision.returncode == 0:
                        result["git_commit"] = revision.stdout.strip()
                except (OSError, subprocess.TimeoutExpired):
                    pass
                break
    return result


def analyze_record(source: Sequence, work: Path, settings: Settings) -> dict[str, Any]:
    work.mkdir()
    fasta = work / "input.fasta"
    write_fasta(fasta, [source])
    runner = Runner(work / "logs", settings.timeout, settings.environment)
    if settings.mode == "model":
        LOG.info("Searching covariance models in %s (%s nt)", source.id, len(source.bases))
        hits = scan_models([source], work, settings, runner)
        return {
            "predictions": merge_model_hits(
                source, hits, settings.merge_gap, settings.type_v_motifs
            ),
            "raw_model_hits": [h.public() for h in hits],
            "arrays": [],
            "cas": {"status": "not_requested"},
        }
    LOG.info("Running CRISPRidentify for %s", source.id)
    arrays, identify_report = run_identify(fasta, source, work, settings, runner)
    LOG.info("Running standalone CasAndra for %s", source.id)
    cas = run_casandra(fasta, source, work, settings, runner)
    LOG.info("Searching anti-repeats for %s (%s accepted arrays)", source.id, len(arrays))
    short_proposals: list[AntiRepeat] | None = [] if settings.short_anti_repeat_rescue else None
    anti = search_anti_repeats(source, fasta, arrays, work, settings, runner, short_proposals)
    candidates = make_candidates(source, anti, settings)
    LOG.info("Scoring %s transcript hypotheses for %s", len(candidates), source.id)
    cm_hits = scan_models(
        [Sequence(c.id, c.id, c.sequence) for c in candidates], work, settings, runner
    )
    by_query: dict[str, list[Hit]] = defaultdict(list)
    for hit in cm_hits:
        by_query[hit.query].append(hit)
    terms = terminators(candidates, work, settings, runner)
    predictions = []
    interaction_map = interactions(candidates, work, settings, runner)
    for candidate in candidates:
        rna_interaction = interaction_map[candidate.id]
        predictions.append(
            candidate_evidence(
                candidate,
                source,
                by_query[candidate.id],
                terms.get(candidate.id, []),
                rna_interaction,
                cas,
                settings,
            )
        )
    genome_hits = []
    if settings.model_rescue:
        # Scan the genome independently: scanning only anti-repeat windows makes
        # complete mode blind to strong CM loci missed by the sequence branch.
        model_work = work / "genomic-models"
        model_work.mkdir()
        genome_hits = scan_models([source], model_work, settings, runner)
        predictions.extend(
            model_rescue_predictions(source, genome_hits, arrays, predictions, cas, settings)
        )
    # Standard CM-only rescue is evaluated before experimental insertions, so
    # enabling the branch never removes or renumbers an existing candidate.
    proposal_records = []
    if short_proposals:
        short_work = work / "short-rescue"
        short_work.mkdir()
        short_candidates = make_candidates(source, short_proposals, settings, "short_candidate")
        short_models: dict[str, list[Hit]] = defaultdict(list)
        for hit in scan_models(
            [Sequence(c.id, c.id, c.sequence) for c in short_candidates],
            short_work,
            settings,
            runner,
        ):
            short_models[hit.query].append(hit)
        short_terms = terminators(short_candidates, short_work, settings, runner)
        short_interactions = interactions(
            short_candidates, short_work, settings, runner, anchored=True
        )
        for candidate in short_candidates:
            models = short_models[candidate.id]
            terms_for_candidate = short_terms.get(candidate.id, [])
            duplex = short_interactions[candidate.id]
            admission = short_rescue_admission(
                candidate, models, terms_for_candidate, duplex, settings
            )
            row = candidate_evidence(
                candidate, source, models, terms_for_candidate, duplex, cas, settings
            )
            row["admission_branch"] = "short_rescue"
            row["short_rescue_admission"] = admission
            row["evidence_class"] = evidence_class(
                duplex is not None,
                bool(admission["supporting_model_hits"]),
                bool(admission["supporting_terminator_hits"]),
            )
            proposal_records.append(row)
            if admission["accepted"]:
                predictions.append(row)
    predictions.sort(key=lambda p: (-p["ranking"]["score"], p["id"]))
    return {
        "predictions": predictions,
        "raw_genomic_model_hits": [h.public() for h in genome_hits],
        "short_rescue_proposals": proposal_records,
        "arrays": [
            a["id"] for a in identify_report["arrays"] if a["category"] in {"Bona-fide", "Possible"}
        ],
        "identify": {"schema": identify_report["schema"], "run": identify_report["run"]},
        "cas": {
            "status": cas["status"],
            "run": cas["run"],
            "cas_gene_count": sum(g["prediction"]["is_cas"] for g in cas["genes"]),
            "cassettes": cas["cassettes"],
        },
    }


def run(input_path: Path, output: Path, settings: Settings) -> dict[str, Any]:
    files = input_files(input_path)
    readiness = doctor(settings)
    if not readiness["ready"]:
        raise TracrError("Preflight failed:\n" + "\n".join(readiness["errors"]))
    output = output.resolve()
    if output.exists():
        raise TracrError(f"Refusing to overwrite existing output: {output}")
    output.mkdir(parents=True)
    started = time.monotonic()
    records: list[dict[str, Any]] = []
    try:
        work = output / "work"
        work.mkdir()
        config = settings.public()
        atomic_json(output / "configuration.json", config)
        for file_index, path in enumerate(files, 1):
            input_digest = sha256_file(path)
            for sequence_index, original in enumerate(read_fasta(path), 1):
                key = f"f{file_index:06d}.s{sequence_index:06d}"
                # Tools receive safe opaque IDs; the original ID and header are retained.
                source = Sequence(key, key, original.bases)
                result = analyze_record(source, work / key, settings)
                record = {
                    "key": key,
                    "source": {
                        "file": str(path),
                        "file_sha256": input_digest,
                        "id": original.id,
                        "header": original.description,
                        "length": len(source.bases),
                        "sha256": source.sha256,
                        "normalization": "uppercase; U converted to T",
                    },
                    **result,
                }
                _assign_export_ids(record)
                atomic_json(work / key / "record.json", record)
                records.append(record)
        write_exports(output, records)
        report = {
            "schema": {"name": "CRISPRtracrRNA-report", "version": "3.0.0"},
            "program_version": __version__,
            "status": "complete",
            "created_at": datetime.now(UTC).isoformat(),
            "wall_seconds": round(time.monotonic() - started, 6),
            "coordinate_system": {
                "indexing": 1,
                "end_inclusive": True,
                "reference": "input_forward_strand",
                "strand_is_separate": True,
            },
            "local_coordinate_reference": "oriented_prediction_sequence",
            "configuration": config,
            "tool_provenance": {
                "identify": command_provenance(settings.identify_command, settings.environment),
                "casandra": command_provenance(settings.casandra_command, settings.environment),
            },
            "model_sha256": {str(p): sha256_file(p) for p in settings.model_paths()},
            "summary": {
                "sequence_count": len(records),
                "prediction_count": sum(len(r["predictions"]) for r in records),
            },
            "records": records,
        }
        # Hash all retained evidence before publishing the completion marker last.
        manifest = {
            str(p.relative_to(output)): sha256_file(p)
            for p in sorted(output.rglob("*"))
            if p.is_file()
        }
        atomic_json(output / "manifest.json", {"algorithm": "sha256", "files": manifest})
        atomic_json(output / "report.json", report)
        return report
    except BaseException as exc:
        try:
            atomic_json(
                output / "failure.json",
                {"status": "failed", "error": str(exc), "completed_sequences": len(records)},
            )
        except BaseException:
            LOG.exception("Could not publish failure marker at %s", output / "failure.json")
        raise
