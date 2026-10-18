"""Independent published-estimate and effective-correction checks."""
import math
from collections.abc import Mapping

_FINITE = {"ln_lh_bdm", "ln_lh_idm", "reversed_ln_lh_bdm", "reversed_ln_lh_idm",
           "ln_lh_bdm - reversed_ln_lh_bdm", "test_statistic (-2*ln_lh_ratio)"}
_POSITIVE = {"deletion_rate_bdm", "insertion_rate_bdm", "per_spacer_deletion_rate_bdm", "alpha_bdm",
             "deletion_rate_idm", "insertion_rate_idm", "avg array length", "tree length", "tree height",
             "run_time", "chi2_quantile", "significance level"}
_WARNING = "rho_bias_fit_unavailable_used_uncorrected"


COUNTS = {'nb of reconstructed deletions', 'nb of leafs (after combining non-uniques)', 'nb of reconstructed independent gains', 'nb of unique spacers', 'nb of reconstructed rearrangements', 'nb of unique spacer arrays', 'nb of reconstructed other dup. events', 'nb of leafs (before combining non-uniques)', 'nb of reconstructed insertions', 'nb of reconstructed ind. acquisitions, not duplicates', 'nb of spacers in alignment', 'nb of reconstructed dup/rearr candidates', 'nb of reconstructed reacquisitions', 'nb of reconstructed duplications'}

def valid_metrics(row):
    for key, raw in row.items():
        is_count = key in COUNTS
        if key not in _FINITE and key not in _POSITIVE and not is_count:
            continue
        if raw is None or not raw.strip():
            continue
        try:
            number = float(raw)
        except (ValueError, TypeError):
            return False
        if not math.isfinite(number) or ((is_count or key in _POSITIVE) and number < 0):
            return False
        if is_count and not number.is_integer() or key == "significance level" and number > 1:
            return False
    return row.get("test result", "").strip().lower() in {"", "true", "false"}


def valid_policy(policy, requested, failure, scope):
    return (isinstance(policy, Mapping) and policy.get("alpha_requested") is requested
            and policy.get("rho_requested") is requested and policy.get("rho_failure") == failure
            and policy.get("fallback_scope") == scope
            and isinstance(policy.get("orientation_likelihood_effect"), str)
            and bool(policy["orientation_likelihood_effect"].strip()))


def correction(effective, requested, failure, *, allow_not_applicable=False):
    if not isinstance(effective, Mapping):
        return False
    alpha = {"applied"} if requested else {"disabled"}
    rho = {"applied", "fallback_uncorrected"} if requested else {"disabled"}
    if allow_not_applicable:
        alpha.add("not_applicable")
        rho.add("not_applicable")
    if effective.get("alpha") not in alpha or effective.get("rho") not in rho:
        return False
    code, warning = effective.get("warning_code"), effective.get("warning")
    if effective["rho"] == "fallback_uncorrected":
        return failure == "fallback_uncorrected" and code == _WARNING and isinstance(warning, str) and bool(warning.strip())
    return code is None and warning is None


def csv_correction(row):
    fields = ("alpha_bias_correction_status", "rho_bias_correction_status",
              "rho_bias_correction_warning_code", "rho_bias_correction_warning", "rho_bias_failure_policy")
    if any(field not in row for field in fields):
        return None
    return {"alpha": row[fields[0]].strip(), "rho": row[fields[1]].strip(),
            "warning_code": row[fields[2]].strip() or None, "warning": row[fields[3]].strip() or None}


def stage_corrections(stage, policy, warnings):
    return (stage.get("bias_correction_policy") == policy and stage.get("bias_fallback_count") == len(warnings)
            and type(stage.get("bias_fallback_count")) is int and stage.get("bias_warnings") == warnings)
