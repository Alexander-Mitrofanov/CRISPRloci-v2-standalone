"""Validated, provider-neutral tabular output for Cas predictions."""

import copy
import csv
import importlib.util
import math
from pathlib import Path

try:
    from components import output_transaction as _output_publication
except ImportError as error:
    # Some distribution tests load this module directly by file path.
    if "output_transaction" not in str(error):
        raise
    _SPEC = importlib.util.spec_from_file_location(
        "v2_output_transaction_for_cas",
        Path(__file__).with_name("output_transaction.py"),
    )
    _output_publication = importlib.util.module_from_spec(_SPEC)
    _SPEC.loader.exec_module(_output_publication)


OutputTransaction = _output_publication.OutputTransaction


SYSTEM_FIELDS = (
    "cassette_id",
    "provider_cassette_id",
    "source_id",
    "start",
    "end",
    "type",
    "classification_class",
    "classification_type",
    "classification_subtype",
    "classification_method",
    "confidence",
    "confidence_kind",
    "experimental",
    "experimental_reason",
    "cas_gene_count",
    "context_gene_count",
    "non_cas_context_gene_count",
)
GENE_FIELDS = (
    "cassette_id",
    "provider_cassette_id",
    "source_id",
    "protein_id",
    "start",
    "end",
    "strand",
    "is_cas",
    "annotation",
    "profile",
    "score",
    "score_kind",
)

_CAS_TYPE_TO_CLASS = {
    "I": "1", "III": "1", "IV": "1", "VII": "1",
    "II": "2", "V": "2", "VI": "2",
}


class CasReportError(ValueError):
    """A Cas result cannot be represented without changing its claims."""


def _non_empty_text(value, label):
    if not isinstance(value, str) or not value:
        raise CasReportError("{} must be a non-empty string".format(label))
    return value


def _optional_text(value, label):
    if value is None:
        return None
    return _non_empty_text(value, label)


def _integer(value, label, minimum=None):
    if isinstance(value, bool) or not isinstance(value, int):
        raise CasReportError("{} must be an integer".format(label))
    if minimum is not None and value < minimum:
        raise CasReportError(
            "{} must be at least {}".format(label, minimum)
        )
    return value


def _finite_number(value, label, minimum=None, maximum=None):
    if value is None:
        return None
    if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)):
        raise CasReportError("{} must be a finite number or null".format(label))
    value = float(value)
    if minimum is not None and value < minimum:
        raise CasReportError("{} is below {}".format(label, minimum))
    if maximum is not None and value > maximum:
        raise CasReportError("{} is above {}".format(label, maximum))
    return value


def _normalize_classification_evidence(raw, label):
    required_fields = {
        "raw_classification",
        "hierarchy_conflict_masked",
    }
    if not isinstance(raw, dict) or set(raw) != required_fields:
        raise CasReportError("{} must contain the frozen evidence fields".format(
            label
        ))
    raw_classification = raw["raw_classification"]
    classification_fields = {"class", "type", "subtype"}
    if (
            not isinstance(raw_classification, dict)
            or set(raw_classification) != classification_fields):
        raise CasReportError(
            "{}.raw_classification must contain class, type, and subtype"
            .format(label)
        )
    normalized_raw = {}
    for field in ("class", "type", "subtype"):
        normalized_raw[field] = _optional_text(
            raw_classification[field],
            "{}.raw_classification.{}".format(label, field),
        )
    conflict_masked = raw["hierarchy_conflict_masked"]
    if not isinstance(conflict_masked, bool):
        raise CasReportError(
            "{}.hierarchy_conflict_masked must be a boolean".format(label)
        )
    return {
        "raw_classification": normalized_raw,
        "hierarchy_conflict_masked": conflict_masked,
    }


def _normalize_final_classification_evidence(raw, classification, label):
    evidence = _normalize_classification_evidence(raw, label)
    raw_classification = evidence["raw_classification"]
    raw_class = raw_classification["class"]
    raw_type = raw_classification["type"]
    raw_subtype = raw_classification["subtype"]
    if raw_class not in ("1", "2"):
        raise CasReportError(
            "{}.raw_classification.class is outside the frozen enum"
            .format(label)
        )
    if raw_type not in _CAS_TYPE_TO_CLASS:
        raise CasReportError(
            "{}.raw_classification.type is outside the frozen enum"
            .format(label)
        )
    expected = {
        "class": _CAS_TYPE_TO_CLASS[raw_type],
        "type": raw_type,
        "subtype": (
            raw_subtype
            if (
                raw_subtype is not None
                and raw_subtype.startswith(raw_type + "-")
            )
            else None
        ),
    }
    observed = {
        key: classification[key] for key in ("class", "type", "subtype")
    }
    if observed != expected:
        raise CasReportError(
            "{} disagrees with the hierarchy-normalized raw provider claim"
            .format(label)
        )
    expected_masked = raw_classification != expected
    if evidence["hierarchy_conflict_masked"] != expected_masked:
        raise CasReportError(
            "{}.hierarchy_conflict_masked disagrees with raw classification"
            .format(label)
        )
    return evidence


def _normalize_gene(
        raw_gene, *, label, source_id=None, source_length=None,
        container_start=None, container_end=None, require_evidence=False):
    if not isinstance(raw_gene, dict):
        raise CasReportError("{} must be an object".format(label))
    gene = copy.deepcopy(raw_gene)
    protein_id = _non_empty_text(gene.get("protein_id"), label + ".protein_id")
    gene_source = _non_empty_text(gene.get("source_id"), label + ".source_id")
    if source_id is not None and gene_source != source_id:
        raise CasReportError("{} source ID disagrees with input".format(label))

    start = _integer(gene.get("start"), label + ".start", 1)
    end = _integer(gene.get("end"), label + ".end", start)
    if source_length is not None and end > source_length:
        raise CasReportError("{} lies outside the input sequence".format(label))
    if (
            container_start is not None
            and (start < container_start or end > container_end)):
        raise CasReportError("{} lies outside its interval".format(label))
    strand = gene.get("strand")
    if strand not in ("+", "-"):
        raise CasReportError("{}.strand must be + or -".format(label))

    annotation = _non_empty_text(
        gene.get("annotation", "unknown"), label + ".annotation"
    )
    is_cas = gene.get("is_cas")
    if is_cas is None:
        # Compatibility with the former legacy projection. New providers
        # must make this biological claim explicit.
        if require_evidence:
            raise CasReportError("{}.is_cas must be a boolean".format(label))
        is_cas = annotation != "unknown"
    if not isinstance(is_cas, bool):
        raise CasReportError("{}.is_cas must be a boolean".format(label))

    profile = gene.get("profile")
    score = gene.get("score")
    score_kind = gene.get("score_kind")
    if "profile" not in gene and not require_evidence:
        profile = annotation if is_cas and annotation != "unknown" else None
    if "score" not in gene and "bitscore" in gene:
        score = gene.get("bitscore")
    if score_kind is None and score is not None and not require_evidence:
        score_kind = "HMM bit score" if "bitscore" in gene else "provider score"
    profile = _optional_text(profile, label + ".profile")
    score = _finite_number(score, label + ".score")
    score_kind = _optional_text(score_kind, label + ".score_kind")
    if require_evidence and score is not None and score_kind is None:
        raise CasReportError(
            "{}.score_kind is required when score is reported".format(label)
        )
    if require_evidence:
        gene["classification_evidence"] = _normalize_classification_evidence(
            gene.get("classification_evidence"), label + ".classification_evidence"
        )

    gene.update({
        "annotation": annotation,
        "end": end,
        "is_cas": is_cas,
        "profile": profile,
        "protein_id": protein_id,
        "score": score,
        "score_kind": score_kind,
        "source_id": gene_source,
        "start": start,
        "strand": strand,
    })
    return gene


def _legacy_probabilities(system, label):
    probabilities = system.get("probabilities", [])
    if not isinstance(probabilities, list):
        raise CasReportError("{}.probabilities must be a list".format(label))
    normalized = []
    seen = set()
    for index, raw_score in enumerate(probabilities, 1):
        score_label = "{}.probabilities[{}]".format(label, index)
        if not isinstance(raw_score, dict):
            raise CasReportError("{} must be an object".format(score_label))
        cas_type = _non_empty_text(
            raw_score.get("cas_type"), score_label + ".cas_type"
        )
        if cas_type in seen:
            raise CasReportError("{} repeats a Cas type".format(score_label))
        seen.add(cas_type)
        value = _finite_number(
            raw_score.get("probability"),
            score_label + ".probability",
            minimum=0.0,
            maximum=1.0,
        )
        normalized.append({"cas_type": cas_type, "probability": value})
    return normalized


def _legacy_classification(system, probabilities, method):
    display_type = system["type"]
    provider_label = (
        display_type[4:] if display_type.startswith("CAS-") else display_type
    )
    parts = provider_label.split("-", 1)
    confidence = probabilities[0]["probability"] if probabilities else None
    classifier = method.get("classifier", "legacy classifier")
    return {
        "class": None,
        "type": parts[0] if parts else None,
        "subtype": provider_label,
        "method": str(classifier),
        "confidence": confidence,
        "confidence_kind": (
            "legacy_classifier_score" if confidence is not None else None
        ),
    }


def _normalize_classification(raw, label):
    if not isinstance(raw, dict):
        raise CasReportError("{} must be an object".format(label))
    classification = copy.deepcopy(raw)
    for field in ("class", "type", "subtype"):
        classification[field] = _optional_text(
            classification.get(field), "{}.{}".format(label, field)
        )
    classification["method"] = _non_empty_text(
        classification.get("method"), label + ".method"
    )
    classification["confidence"] = _finite_number(
        classification.get("confidence"),
        label + ".confidence",
        minimum=0.0,
        maximum=1.0,
    )
    classification["confidence_kind"] = _optional_text(
        classification.get("confidence_kind"), label + ".confidence_kind"
    )
    if (
            classification["confidence"] is not None
            and classification["confidence_kind"] is None):
        raise CasReportError(
            "{}.confidence_kind is required when confidence is reported"
            .format(label)
        )
    return classification


def _normalize_experimental_claim(
        system, classification, label, require_casandra):
    if not require_casandra:
        return (
            system.get("experimental", False),
            system.get("experimental_reason"),
            system.get("evidence"),
        )

    expected = (
        classification["type"] == "VII"
        and classification["subtype"] is None
    )
    experimental = system.get("experimental")
    if not isinstance(experimental, bool):
        raise CasReportError(
            "{}.experimental must be a boolean".format(label)
        )
    if experimental != expected:
        raise CasReportError(
            "{} experimental maturity disagrees with classification"
            .format(label)
        )
    reason = _optional_text(
        system.get("experimental_reason"),
        label + ".experimental_reason",
    )
    if expected and reason is None:
        raise CasReportError(
            "{} experimental Type VII claim requires a reason".format(label)
        )
    if not expected and reason is not None:
        raise CasReportError(
            "{} non-experimental system cannot carry an experimental reason"
            .format(label)
        )

    evidence = system.get("evidence")
    if not isinstance(evidence, dict):
        raise CasReportError("{}.evidence must be an object".format(label))
    required_evidence_fields = {
        "experimental_reason",
        "genome_evidence_gate",
        "provider_direct_profile_result",
        "provider_final_classification_evidence",
        "type_vii_experimental",
    }
    if set(evidence) != required_evidence_fields:
        raise CasReportError(
            "{}.evidence must contain the frozen CasAndra fields"
            .format(label)
        )
    evidence = copy.deepcopy(evidence)
    evidence["genome_evidence_gate"] = _non_empty_text(
        evidence["genome_evidence_gate"],
        label + ".evidence.genome_evidence_gate",
    )
    if not isinstance(evidence["provider_direct_profile_result"], dict):
        raise CasReportError(
            "{}.evidence.provider_direct_profile_result must be an object"
            .format(label)
        )
    evidence["provider_final_classification_evidence"] = (
        _normalize_final_classification_evidence(
            evidence["provider_final_classification_evidence"],
            classification,
            label + ".evidence.provider_final_classification_evidence",
        )
    )
    direct = evidence["provider_direct_profile_result"]
    direct_claim = {
        field: _optional_text(
            direct.get(field),
            "{}.evidence.provider_direct_profile_result.{}".format(
                label, field
            ),
        )
        for field in ("class", "type", "subtype")
    }
    raw_final = evidence[
        "provider_final_classification_evidence"
    ]["raw_classification"]
    if (
            raw_final["class"] != direct_claim["class"]
            or raw_final["type"] != direct_claim["type"]):
        raise CasReportError(
            "{} raw final class/type disagrees with direct provider evidence"
            .format(label)
        )
    route = classification["method"]
    if route in (
            "direct_profile_aggregation",
            "direct_profile_type_aggregation"):
        expected_route = (
            "direct_profile_aggregation"
            if direct_claim["subtype"] is not None
            else "direct_profile_type_aggregation"
        )
        if (
                route != expected_route
                or raw_final["subtype"] != direct_claim["subtype"]):
            raise CasReportError(
                "{} direct route disagrees with provider evidence"
                .format(label)
            )
    elif route == "type_iii_ordered_profile_architecture_extratrees":
        if (
                direct_claim["type"] != "III"
                or raw_final["subtype"] is None
                or not raw_final["subtype"].startswith("III-")):
            raise CasReportError(
                "{} Type III architecture route is inconsistent"
                .format(label)
            )
    else:
        raise CasReportError(
            "{} uses an unsupported CasAndra classification route"
            .format(label)
        )
    evidence_experimental = evidence.get("type_vii_experimental")
    if not isinstance(evidence_experimental, bool):
        raise CasReportError(
            "{}.evidence.type_vii_experimental must be a boolean"
            .format(label)
        )
    if evidence_experimental != experimental:
        raise CasReportError(
            "{} evidence maturity disagrees with system maturity".format(label)
        )
    evidence_reason = _optional_text(
        evidence.get("experimental_reason"),
        label + ".evidence.experimental_reason",
    )
    if evidence_reason != reason:
        raise CasReportError(
            "{} evidence experimental reason disagrees with system"
            .format(label)
        )
    return experimental, reason, evidence


def _normalize_method(raw_method, *, require_casandra):
    if not isinstance(raw_method, dict):
        raise CasReportError("computed Cas method must be an object")
    method = copy.deepcopy(raw_method)
    if not require_casandra:
        provider = method.get("provider")
        if provider is None:
            method["provider"] = "CRISPRcasIdentifier (legacy)"
        else:
            _non_empty_text(provider, "Cas method.provider")
        return method

    text_fields = (
        "provider_version",
        "bundle_id",
        "bundle_role",
        "gene_caller",
        "gene_mode",
        "model_manifest_sha256",
    )
    for field in text_fields:
        method[field] = _non_empty_text(
            method.get(field), "Cas method.{}".format(field)
        )
    if method["provider"] != "CasAndra":
        raise CasReportError("CasAndra provider name must be 'CasAndra'")
    if method.get("pipeline_schema_version") != 4:
        raise CasReportError("CasAndra pipeline schema version must be 4")
    _integer(method.get("translation_table"), "Cas method.translation_table", 1)
    digest = method["model_manifest_sha256"]
    if (
            len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)):
        raise CasReportError(
            "Cas method.model_manifest_sha256 must be lowercase SHA-256"
        )
    return method


def _normalize_run_statistics(raw, top_genes, system_count):
    required_fields = {
        "raw_cassette_candidates",
        "rejected_cassette_candidates",
        "hierarchy_conflict_masked_cassettes",
        "hierarchy_conflict_masked_proteins",
    }
    if not isinstance(raw, dict) or set(raw) != required_fields:
        raise CasReportError(
            "CasAndra run_statistics must contain the frozen fields"
        )
    statistics = copy.deepcopy(raw)
    for field in required_fields:
        statistics[field] = _integer(
            statistics[field],
            "Cas run_statistics.{}".format(field),
            0,
        )
    if (
            statistics["raw_cassette_candidates"]
            != system_count + statistics["rejected_cassette_candidates"]):
        raise CasReportError(
            "CasAndra cassette candidate statistics are inconsistent"
        )
    if statistics["hierarchy_conflict_masked_cassettes"] > system_count:
        raise CasReportError(
            "CasAndra masked cassette count exceeds reported systems"
        )
    observed_conflicts = sum(
        gene["classification_evidence"]["hierarchy_conflict_masked"]
        for gene in top_genes
    )
    if (
            statistics["hierarchy_conflict_masked_proteins"]
            != observed_conflicts):
        raise CasReportError(
            "CasAndra hierarchy conflict count disagrees with gene evidence"
        )
    return statistics


def normalize_cas_data(cas_data, source_id=None, source_length=None):
    """Return one validated canonical Cas report object.

    CasAndra projections are checked strictly. The previous in-process
    CRISPRcasIdentifier projection remains readable and is upgraded with
    provider-neutral evidence fields without relabelling its scores as
    calibrated probabilities.
    """

    if not isinstance(cas_data, dict):
        raise CasReportError("Cas report data must be an object")
    status = cas_data.get("status")
    if status == "not_requested":
        systems = cas_data.get("systems", [])
        if systems != []:
            raise CasReportError("not-requested Cas data cannot contain systems")
        result = copy.deepcopy(cas_data)
        result["systems"] = []
        return result
    if status != "computed":
        raise CasReportError("Cas report status must be computed or not_requested")

    provider = (
        cas_data.get("method", {}).get("provider")
        if isinstance(cas_data.get("method"), dict)
        else None
    )
    require_casandra = provider == "CasAndra"
    method = _normalize_method(
        cas_data.get("method", {}), require_casandra=require_casandra
    )
    coordinate_system = cas_data.get("coordinate_system")
    if require_casandra:
        if coordinate_system != {
                "end_inclusive": True,
                "indexing": 1,
                "reference": "input_forward_strand"}:
            raise CasReportError("CasAndra coordinate-system contract is invalid")

    systems = cas_data.get("systems")
    if not isinstance(systems, list):
        raise CasReportError("computed Cas systems must be a list")
    if cas_data.get("system_count") != len(systems):
        raise CasReportError("Cas system count disagrees with systems")

    top_genes_raw = cas_data.get("genes", [])
    if not isinstance(top_genes_raw, list):
        raise CasReportError("computed Cas genes must be a list")
    top_genes = []
    top_by_id = {}
    for index, raw_gene in enumerate(top_genes_raw, 1):
        gene = _normalize_gene(
            raw_gene,
            label="Cas genes[{}]".format(index),
            source_id=source_id,
            source_length=source_length,
            require_evidence=require_casandra,
        )
        if gene["protein_id"] in top_by_id:
            raise CasReportError("Cas top-level protein IDs must be unique")
        if require_casandra and not gene["is_cas"]:
            raise CasReportError("CasAndra top-level genes must be Cas-positive")
        top_by_id[gene["protein_id"]] = gene
        top_genes.append(gene)
    if "gene_count" in cas_data and cas_data["gene_count"] != len(top_genes):
        raise CasReportError("Cas gene count disagrees with genes")
    if require_casandra and "gene_count" not in cas_data:
        raise CasReportError("CasAndra result is missing gene_count")
    if "protein_count" in cas_data:
        protein_count = _integer(
            cas_data["protein_count"], "Cas protein_count", 0
        )
        if protein_count < len(top_genes):
            raise CasReportError("Cas protein_count is smaller than gene_count")
    run_statistics = None
    if require_casandra:
        run_statistics = _normalize_run_statistics(
            cas_data.get("run_statistics"),
            top_genes,
            len(systems),
        )

    normalized_systems = []
    provider_ids = set()
    previous_interval = None
    for position, raw_system in enumerate(systems, 1):
        label = "Cas systems[{}]".format(position)
        if not isinstance(raw_system, dict):
            raise CasReportError("{} must be an object".format(label))
        system = copy.deepcopy(raw_system)
        cassette_id = _integer(
            system.get("cassette_id"), label + ".cassette_id", 1
        )
        if cassette_id != position:
            raise CasReportError("Cas cassette display IDs must be sequential")
        provider_cassette_id = system.get("provider_cassette_id")
        if provider_cassette_id is None and not require_casandra:
            provider_cassette_id = str(cassette_id)
        provider_cassette_id = _non_empty_text(
            provider_cassette_id, label + ".provider_cassette_id"
        )
        if provider_cassette_id in provider_ids:
            raise CasReportError("Cas provider cassette IDs must be unique")
        provider_ids.add(provider_cassette_id)

        system_source = _non_empty_text(
            system.get("source_id"), label + ".source_id"
        )
        if source_id is not None and system_source != source_id:
            raise CasReportError("{} source ID disagrees with input".format(label))
        start = _integer(system.get("start"), label + ".start", 1)
        end = _integer(system.get("end"), label + ".end", start)
        if source_length is not None and end > source_length:
            raise CasReportError("{} lies outside the input sequence".format(label))
        interval_key = (start, end, provider_cassette_id)
        if previous_interval is not None and interval_key < previous_interval:
            raise CasReportError("Cas systems must be in source-coordinate order")
        previous_interval = interval_key
        display_type = _non_empty_text(system.get("type"), label + ".type")

        probabilities = _legacy_probabilities(system, label)
        raw_classification = system.get("classification")
        if raw_classification is None:
            if require_casandra:
                raise CasReportError("{} is missing classification".format(label))
            classification = _legacy_classification(system, probabilities, method)
        else:
            classification = _normalize_classification(
                raw_classification, label + ".classification"
            )
        if require_casandra:
            expected_type = classification["subtype"] or classification["type"]
            if expected_type is None or display_type != "CAS-" + expected_type:
                raise CasReportError(
                    "{} display type disagrees with CasAndra classification"
                    .format(label)
                )

        experimental, experimental_reason, evidence = (
            _normalize_experimental_claim(
                system, classification, label, require_casandra
            )
        )

        for duplicate_field in ("confidence", "confidence_kind"):
            if (
                    duplicate_field in system
                    and system[duplicate_field] != classification[duplicate_field]):
                raise CasReportError(
                    "{}.{} disagrees with classification".format(
                        label, duplicate_field
                    )
                )

        raw_genes = system.get("genes")
        if not isinstance(raw_genes, list) or not raw_genes:
            raise CasReportError("{} must contain context genes".format(label))
        genes = []
        context_ids = set()
        previous_gene_key = None
        for gene_position, raw_gene in enumerate(raw_genes, 1):
            gene = _normalize_gene(
                raw_gene,
                label="{}.genes[{}]".format(label, gene_position),
                source_id=system_source,
                source_length=source_length,
                container_start=start,
                container_end=end,
                require_evidence=require_casandra,
            )
            gene_key = (gene["start"], gene["end"], gene["protein_id"])
            if previous_gene_key is not None and gene_key < previous_gene_key:
                raise CasReportError("{} genes must be source ordered".format(label))
            previous_gene_key = gene_key
            if gene["protein_id"] in context_ids:
                raise CasReportError("{} repeats a protein ID".format(label))
            context_ids.add(gene["protein_id"])
            if require_casandra and gene["is_cas"]:
                top_gene = top_by_id.get(gene["protein_id"])
                if top_gene is None:
                    raise CasReportError(
                        "{} Cas-positive context gene is absent from top-level genes"
                        .format(label)
                    )
                for field in (
                        "source_id", "start", "end", "strand", "is_cas",
                        "annotation", "classification_evidence", "profile",
                        "score", "score_kind"):
                    if gene[field] != top_gene[field]:
                        raise CasReportError(
                            "{} disagrees with top-level gene {}".format(
                                label, gene["protein_id"]
                            )
                        )
            genes.append(gene)

        non_cas_count = sum(not gene["is_cas"] for gene in genes)
        if (
                "missing_gene_count" in system
                and system["missing_gene_count"] != non_cas_count):
            raise CasReportError(
                "{}.missing_gene_count must count non-Cas context genes"
                .format(label)
            )
        system.update({
            "cassette_id": cassette_id,
            "classification": classification,
            "confidence": classification["confidence"],
            "confidence_kind": classification["confidence_kind"],
            "end": end,
            "experimental": experimental,
            "experimental_reason": experimental_reason,
            "evidence": evidence,
            "genes": genes,
            "missing_gene_count": non_cas_count,
            "provider_cassette_id": provider_cassette_id,
            "source_id": system_source,
            "start": start,
            "type": display_type,
        })
        if "probabilities" in raw_system or not require_casandra:
            # Preserve the legacy field, but never manufacture it for a new
            # provider projection.
            system["probabilities"] = probabilities
        normalized_systems.append(system)

    if require_casandra:
        observed_masked_cassettes = sum(
            system["evidence"][
                "provider_final_classification_evidence"
            ]["hierarchy_conflict_masked"]
            for system in normalized_systems
        )
        if run_statistics["hierarchy_conflict_masked_cassettes"] != (
                observed_masked_cassettes):
            raise CasReportError(
                "CasAndra masked cassette count disagrees with system evidence"
            )

    result = copy.deepcopy(cas_data)
    result.update({
        "genes": top_genes,
        "gene_count": len(top_genes),
        "method": method,
        "system_count": len(normalized_systems),
        "systems": normalized_systems,
    })
    if require_casandra:
        result["run_statistics"] = run_statistics
    return result


def _stage_csv(publication, path, fields, rows):
    def write(stream):
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    publication.stage_writer(path, write, newline="")


def stage_cas_outputs(
        publication, result_path, cas_data, *,
        source_id=None, source_length=None):
    result_path = Path(result_path)
    systems_path = result_path / "Cas_systems.csv"
    genes_path = result_path / "Cas_genes.csv"
    normalized = normalize_cas_data(
        cas_data, source_id=source_id, source_length=source_length
    )
    if normalized["status"] == "not_requested":
        for path in (systems_path, genes_path):
            publication.delete(path)
        return normalized

    system_rows = []
    gene_rows = []
    seen_context_proteins = set()
    for system in normalized["systems"]:
        classification = system["classification"]
        cas_gene_count = sum(gene["is_cas"] for gene in system["genes"])
        system_rows.append({
            "cassette_id": system["cassette_id"],
            "provider_cassette_id": system["provider_cassette_id"],
            "source_id": system["source_id"],
            "start": system["start"],
            "end": system["end"],
            "type": system["type"],
            "classification_class": classification["class"],
            "classification_type": classification["type"],
            "classification_subtype": classification["subtype"],
            "classification_method": classification["method"],
            "confidence": classification["confidence"],
            "confidence_kind": classification["confidence_kind"],
            "experimental": system["experimental"],
            "experimental_reason": system["experimental_reason"],
            "cas_gene_count": cas_gene_count,
            "context_gene_count": len(system["genes"]),
            "non_cas_context_gene_count": len(system["genes"]) - cas_gene_count,
        })
        for gene in system["genes"]:
            seen_context_proteins.add(gene["protein_id"])
            gene_rows.append({
                "cassette_id": system["cassette_id"],
                "provider_cassette_id": system["provider_cassette_id"],
                "source_id": gene["source_id"],
                "protein_id": gene["protein_id"],
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "is_cas": gene["is_cas"],
                "annotation": gene["annotation"],
                "profile": gene["profile"],
                "score": gene["score"],
                "score_kind": gene["score_kind"],
            })

    # CasAndra reports Cas-positive proteins that did not enter an accepted
    # cassette. Preserve them instead of silently losing that evidence.
    for gene in normalized["genes"]:
        if gene["protein_id"] in seen_context_proteins:
            continue
        gene_rows.append({
            "cassette_id": None,
            "provider_cassette_id": None,
            "source_id": gene["source_id"],
            "protein_id": gene["protein_id"],
            "start": gene["start"],
            "end": gene["end"],
            "strand": gene["strand"],
            "is_cas": gene["is_cas"],
            "annotation": gene["annotation"],
            "profile": gene["profile"],
            "score": gene["score"],
            "score_kind": gene["score_kind"],
        })

    _stage_csv(publication, systems_path, SYSTEM_FIELDS, system_rows)
    _stage_csv(publication, genes_path, GENE_FIELDS, gene_rows)
    return normalized


def write_cas_outputs(
        result_path, cas_data, *, source_id=None, source_length=None):
    """Validate and transactionally publish the standalone Cas CSV pair."""

    with OutputTransaction() as publication:
        normalized = stage_cas_outputs(
            publication,
            result_path,
            cas_data,
            source_id=source_id,
            source_length=source_length,
        )
        publication.publish()
    return normalized
