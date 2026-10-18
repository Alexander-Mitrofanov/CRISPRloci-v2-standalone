"""Typed public evolutionary estimates; no scientific computation or fitting."""
import math

SIGNED = {"ln_lh_bdm", "ln_lh_idm", "reversed_ln_lh_bdm", "reversed_ln_lh_idm",
          "ln_lh_bdm - reversed_ln_lh_bdm", "test_statistic (-2*ln_lh_ratio)"}
NONNEGATIVE = {"deletion_rate_bdm", "insertion_rate_bdm", "per_spacer_deletion_rate_bdm",
               "alpha_bdm", "deletion_rate_idm", "insertion_rate_idm", "avg array length",
               "tree length", "tree height", "run_time", "chi2_quantile", "significance level"}


COUNTS = {'nb of reconstructed deletions', 'nb of leafs (after combining non-uniques)', 'nb of reconstructed independent gains', 'nb of unique spacers', 'nb of reconstructed rearrangements', 'nb of unique spacer arrays', 'nb of reconstructed other dup. events', 'nb of leafs (before combining non-uniques)', 'nb of reconstructed insertions', 'nb of reconstructed ind. acquisitions, not duplicates', 'nb of spacers in alignment', 'nb of reconstructed dup/rearr candidates', 'nb of reconstructed reacquisitions', 'nb of reconstructed duplications'}

def validate_metrics(row):
    for field, raw in row.items():
        count = field in COUNTS
        if field not in SIGNED | NONNEGATIVE and not count:
            continue
        if raw is None or not raw.strip():
            continue  # Optional absence remains missing; required fields are checked by the caller.
        try:
            value = float(raw)
        except (TypeError, ValueError) as error:
            raise ValueError(f"result CSV field {field!r} is invalid") from error
        if not math.isfinite(value):
            raise ValueError(f"result CSV field {field!r} is non-finite")
        if count and (value < 0 or not value.is_integer()):
            raise ValueError(f"result CSV field {field!r} must be a non-negative integer")
        if field in NONNEGATIVE and value < 0:
            raise ValueError(f"result CSV field {field!r} must be non-negative")
        if field == "significance level" and value > 1:
            raise ValueError("Invalid significance level")
    if row.get("test result", "").strip().lower() not in {"", "true", "false"}:
        raise ValueError("Invalid evolutionary test result")
