from __future__ import annotations

import json
from pathlib import Path

from crispr_web.api import create_app
from crispr_web.models import CategoryPolicy, JobStatus, ScientificMode


CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "frontend"
    / "src"
    / "contracts"
    / "public-api-v1.json"
)
HTTP_METHODS = frozenset({"get", "post", "put", "patch", "delete", "options", "head"})


def _schema_name(schema: dict[str, object]) -> str | None:
    reference = schema.get("$ref")
    return str(reference).rsplit("/", maxsplit=1)[-1] if reference else None


def _schema_signature(schema: dict[str, object]) -> str:
    reference = _schema_name(schema)
    if reference:
        return f"ref:{reference}"
    if "anyOf" in schema:
        variants = (_schema_signature(item) for item in schema["anyOf"])
        return "|".join(sorted(variants))
    schema_type = str(schema.get("type"))
    if schema_type == "array":
        return f"array:{_schema_signature(schema['items'])}"
    if schema_type == "object" and isinstance(schema.get("additionalProperties"), dict):
        return f"map:{_schema_signature(schema['additionalProperties'])}"
    return schema_type


def _schema_references(value: object) -> set[str]:
    if isinstance(value, dict):
        references = {
            name
            for candidate in [value]
            if (name := _schema_name(candidate)) is not None
        }
        return references | {
            reference
            for child in value.values()
            for reference in _schema_references(child)
        }
    if isinstance(value, list):
        return {reference for child in value for reference in _schema_references(child)}
    return set()


def test_frontend_contract_matches_backend_openapi_and_enums(settings_factory) -> None:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    schema = create_app(settings_factory()).openapi()

    expected_routes = {
        (route["path"], route["method"].lower())
        for route in contract["routes"].values()
    }
    actual_routes = {
        (path, method)
        for path, operations in schema["paths"].items()
        for method in operations
        if method in HTTP_METHODS
    }

    assert contract["api_version"] == "v1"
    parameters = schema["paths"]["/api/v1/catalog/repeats"]["get"]["parameters"]
    match_parameter = next(item for item in parameters if item["name"] == "match")
    assert match_parameter["schema"]["enum"] == contract["enums"]["catalog_match_modes"]
    assert actual_routes == expected_routes
    assert set(contract["enums"]["analysis_modes"]) == {
        value.value for value in ScientificMode
    }
    assert set(contract["enums"]["category_policies"]) == {
        value.value for value in CategoryPolicy
    }
    assert set(contract["enums"]["job_statuses"]) == {
        value.value for value in JobStatus
    }
    assert set(contract["enums"]["terminal_statuses"]) == {
        JobStatus.completed.value,
        JobStatus.completed_no_eligible_groups.value,
        JobStatus.failed.value,
        JobStatus.cancelled.value,
        "expired",
    }

    openapi_contract = contract["openapi"]
    scheme_name = openapi_contract["security_scheme"]
    assert schema["components"]["securitySchemes"][scheme_name] == {
        "type": "http",
        "scheme": "bearer",
    }

    for route_name, expected in openapi_contract["operations"].items():
        route = contract["routes"][route_name]
        operation = schema["paths"][route["path"]][route["method"].lower()]
        response = operation["responses"][str(expected["success_status"])]
        assert set(response.get("content", {})) == set(expected["response_media_types"])
        response_schema = expected["response_schema"]
        if response_schema:
            assert (
                _schema_name(response["content"]["application/json"]["schema"])
                == response_schema
            )
        if request_schema := expected.get("request_schema"):
            request_body = operation["requestBody"]
            assert request_body["required"] is True
            request_media_type = expected["request_media_type"]
            assert set(request_body["content"]) == {request_media_type}
            assert (
                _schema_name(request_body["content"][request_media_type]["schema"])
                == request_schema
            )
        expected_security = [{scheme_name: []}] if expected["protected"] else []
        assert operation.get("security", []) == expected_security

    component_schemas = schema["components"]["schemas"]
    roots = {
        model_name
        for operation in openapi_contract["operations"].values()
        for model_name in (
            operation.get("request_schema"),
            operation.get("response_schema"),
        )
        if model_name
    }
    reachable = set(roots)
    pending = list(roots)
    while pending:
        model_name = pending.pop()
        for reference in _schema_references(component_schemas[model_name]):
            if reference not in reachable:
                reachable.add(reference)
                pending.append(reference)
    reachable_objects = {
        model_name
        for model_name in reachable
        if component_schemas[model_name].get("type") == "object"
    }
    assert set(openapi_contract["models"]) == reachable_objects

    for model_name, expected in openapi_contract["models"].items():
        model = component_schemas[model_name]
        assert set(model["properties"]) == set(expected["properties"])
        assert set(model.get("required", [])) == set(expected["required"])
        assert model.get("additionalProperties") == expected["additional_properties"]
        signatures = {
            name: _schema_signature(property_schema)
            for name, property_schema in model["properties"].items()
        }
        assert signatures == expected["types"]
        nullable = {
            name
            for name, signature in signatures.items()
            if "null" in signature.split("|")
        }
        assert nullable == set(expected["nullable"])
