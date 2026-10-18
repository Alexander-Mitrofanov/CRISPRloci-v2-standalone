"""Bias-correction policy and warning normalization."""

from __future__ import annotations

from typing import Any, Mapping

from .contracts import (
    _BIAS_FALLBACK_SCOPES,
    _BIAS_FALLBACK_WARNING_CODE,
    _BIAS_RHO_FAILURE_POLICIES,
    _BIAS_STAGE_LABELS,
    _BIAS_WARNING_ORIENTATIONS,
    _CREDENTIAL_BEARING_TEXT,
    _MAX_BIAS_GROUP_CHARACTERS,
    _MAX_BIAS_MESSAGE_CHARACTERS,
    _MAX_BIAS_POLICY_CHARACTERS,
    _OPAQUE_TOKEN_TEXT,
)
from .sanitize import _path_bearing_string, _summary_path_markers


def _bounded_bias_text(value: object, max_characters: int) -> str | None:
    if not isinstance(value, str):
        return None
    rendered = value.strip()
    if (
        not rendered
        or len(rendered) > max_characters
        or any(ord(character) < 32 for character in rendered)
        or _path_bearing_string(rendered, _summary_path_markers(()))
        or "../" in rendered.replace("\\", "/")
        or _CREDENTIAL_BEARING_TEXT.search(rendered)
        or _OPAQUE_TOKEN_TEXT.search(rendered)
    ):
        return None
    return rendered


def _bias_policy(stage_name: str, stage: Mapping[str, Any]) -> dict[str, object] | None:
    raw_policy = stage.get("bias_correction_policy")
    if not isinstance(raw_policy, Mapping):
        return None
    alpha_requested = raw_policy.get("alpha_requested")
    rho_requested = raw_policy.get("rho_requested")
    rho_failure = raw_policy.get("rho_failure")
    fallback_scope = raw_policy.get("fallback_scope")
    orientation_effect = _bounded_bias_text(
        raw_policy.get("orientation_likelihood_effect"),
        _MAX_BIAS_POLICY_CHARACTERS,
    )
    if (
        not isinstance(alpha_requested, bool)
        or not isinstance(rho_requested, bool)
        or not isinstance(rho_failure, str)
        or rho_failure not in _BIAS_RHO_FAILURE_POLICIES
        or fallback_scope != _BIAS_FALLBACK_SCOPES[stage_name]
        or orientation_effect is None
    ):
        return None
    return {
        "alpha_requested": alpha_requested,
        "rho_requested": rho_requested,
        "rho_failure": str(rho_failure),
        "fallback_scope": str(fallback_scope),
        "orientation_likelihood_effect": orientation_effect,
    }


def _bias_warning(stage_name: str, value: object) -> dict[str, str] | None:
    if not isinstance(value, Mapping):
        return None
    if value.get("code") != _BIAS_FALLBACK_WARNING_CODE:
        return None
    group = _bounded_bias_text(value.get("group"), _MAX_BIAS_GROUP_CHARACTERS)
    message = _bounded_bias_text(value.get("message"), _MAX_BIAS_MESSAGE_CHARACTERS)
    if group is None or message is None:
        return None
    warning = {
        "group": group,
        "code": _BIAS_FALLBACK_WARNING_CODE,
        "message": message,
    }
    orientation = value.get("orientation")
    if (
        not isinstance(orientation, str)
        or orientation not in _BIAS_WARNING_ORIENTATIONS[stage_name]
    ):
        return None
    warning["orientation"] = orientation
    return warning


def _bias_fallback_metadata(
    stage_name: str,
    stage: Mapping[str, Any],
    policy: Mapping[str, object],
    *,
    max_warnings: int,
) -> tuple[int, list[dict[str, str]]] | None:
    reported_count = stage.get("bias_fallback_count")
    raw_warnings = stage.get("bias_warnings")
    if (
        isinstance(reported_count, bool)
        or not isinstance(reported_count, int)
        or reported_count < 0
        or reported_count > max_warnings
        or not isinstance(raw_warnings, list)
        or len(raw_warnings) != reported_count
    ):
        return None
    if reported_count > 0 and (
        policy.get("rho_requested") is not True
        or policy.get("rho_failure") != "fallback_uncorrected"
    ):
        return None
    normalized: list[dict[str, str]] = []
    for raw_warning in raw_warnings:
        warning = _bias_warning(stage_name, raw_warning)
        if warning is None:
            return None
        normalized.append(warning)
    return reported_count, normalized


def _public_bias_warning(stage_name: str, warning: Mapping[str, str]) -> dict[str, str]:
    orientation = warning.get("orientation")
    direction = f" ({orientation.replace('_', ' ')})" if orientation is not None else ""
    return {
        "stage": stage_name,
        **warning,
        "title": (
            f"{_BIAS_STAGE_LABELS[stage_name]} bias fallback for "
            f"{warning['group']}{direction}"
        ),
    }


