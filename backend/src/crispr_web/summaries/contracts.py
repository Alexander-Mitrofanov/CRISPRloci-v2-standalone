"""Limits, allowlists, and public-summary contract constants."""

from __future__ import annotations

import re


_MAX_SUMMARY_INPUT_BYTES = 20_000_000
_MAX_PUBLIC_SUMMARY_BYTES = 5_000_000
_MAX_RECONSTRUCTION_INPUT_BYTES = 5_000_000
_MAX_RECONSTRUCTION_TOTAL_INPUT_BYTES = 40_000_000
_MAX_RECONSTRUCTION_SPACERS = 10_000
_MAX_RECONSTRUCTION_NODES = 2_000
_MAX_RECONSTRUCTION_STATE_CELLS = 1_000_000
_MAX_RECONSTRUCTION_EVENT_REFERENCES = 100_000
_MAX_RECONSTRUCTION_IDENTITY_CHARACTERS = 2_048
_MAX_RECONSTRUCTION_NAME_CHARACTERS = 240
_MAX_RECONSTRUCTION_WARNINGS = 100
_RECONSTRUCTION_EVENT_FIELDS = {
    "rec_contra_dict": "contradictions",
    "rec_duplications_dict": "duplications",
    "rec_rearrangements_dict": "rearrangements",
    "rec_reacquisition_dict": "reacquisitions",
    "rec_indep_gain_dict": "independent_gains",
    "rec_other_dup_events_dict": "other_duplication_events",
}
_BIAS_FALLBACK_WARNING_CODE = "rho_bias_fit_unavailable_used_uncorrected"
_BIAS_STAGE_NAMES = frozenset({"crispr_evor", "spacerplacer"})
_BIAS_STAGE_LABELS = {
    "crispr_evor": "CRISPR-evOr",
    "spacerplacer": "SpacerPlacer",
}
_MAX_BIAS_WARNINGS = 2_000
_MAX_BIAS_GROUP_CHARACTERS = 240
_MAX_BIAS_MESSAGE_CHARACTERS = 2_000
_MAX_BIAS_POLICY_CHARACTERS = 1_000
_BIAS_RHO_FAILURE_POLICIES = frozenset({"fail", "fallback_uncorrected"})
_BIAS_FALLBACK_SCOPES = {
    "crispr_evor": "per_group_per_orientation",
    "spacerplacer": "per_group",
}
_BIAS_WARNING_ORIENTATIONS = {
    "crispr_evor": frozenset({"forward", "reverse"}),
    "spacerplacer": frozenset({"input_order"}),
}
_CREDENTIAL_BEARING_TEXT = re.compile(
    r"(?i)(?:\b(?:access[ _-]?token|api[ _-]?key|authorization|bearer|secret)\b"
    r"|\bgh[pousr]_[A-Za-z0-9]{20,}\b|\bsk-[A-Za-z0-9_-]{20,}\b)"
)
_OPAQUE_TOKEN_TEXT = re.compile(
    r"(?<![A-Za-z0-9_-])[A-Za-z0-9_-]{40,}(?![A-Za-z0-9_-])"
)
_PUBLIC_RESULT_FIELDS = frozenset(
    {
        "name",
        "deletion_rate_bdm",
        "alpha_bdm",
        "per_spacer_deletion_rate_bdm",
        "insertion_rate_bdm",
        "deletion_rate_idm",
        "insertion_rate_idm",
        "ln_lh_idm",
        "ln_lh_bdm",
        "test_statistic (-2*ln_lh_ratio)",
        "test result",
        "Deletion model preferred by LRT",
        "chi2_quantile",
        "significance level",
        "avg array length",
        "nb of spacers in alignment",
        "nb of unique spacers",
        "nb of reconstructed insertions",
        "nb of reconstructed deletions",
        "nb of reconstructed dup/rearr candidates",
        "nb of reconstructed rearrangements",
        "nb of reconstructed duplications",
        "nb of reconstructed reacquisitions",
        "nb of reconstructed independent gains",
        "nb of reconstructed other dup. events",
        "nb of reconstructed ind. acquisitions, not duplicates",
        "tree height",
        "tree length",
        "min/max tree branch lengths",
        "nb of leafs (before combining non-uniques)",
        "nb of leafs (after combining non-uniques)",
        "nb of unique spacer arrays",
        "run_time",
        "reversed_ln_lh_idm",
        "reversed_ln_lh_bdm",
        "ln_lh_bdm - reversed_ln_lh_bdm",
        "predicted orientation",
        "recommend reversing array",
        "confidence threshold",
    }
)


