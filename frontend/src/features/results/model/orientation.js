import { normalizeReconstructionRow, normalizeTree } from "./reconstruction.js";
import { asList, firstDefined, firstRecord } from "./values.js";

function normalizeComparison(value, index) {
  const comparison = firstRecord(value) || {};
  return {
    ...comparison,
    confidence_threshold: firstDefined(comparison, "confidence_threshold", "threshold"),
    forward_ln_likelihood_bdm: firstDefined(
      comparison,
      "forward_ln_likelihood_bdm",
      "forward_ln_likelihood",
    ),
    forward_minus_reverse_ln_likelihood_bdm: firstDefined(
      comparison,
      "forward_minus_reverse_ln_likelihood_bdm",
      "delta_ln_likelihood",
      "delta_lnL",
      "delta_log_likelihood",
    ),
    group: String(
      firstDefined(comparison, "group", "name", "group_id", "id") || `group_${index + 1}`,
    ),
    reverse_ln_likelihood_bdm: firstDefined(
      comparison,
      "reverse_ln_likelihood_bdm",
      "reverse_ln_likelihood",
    ),
  };
}

function normalizeHistory(value, index) {
  const history = firstRecord(value) || {};
  return {
    ...history,
    group: String(firstDefined(history, "group", "name", "group_id", "id") || `group_${index + 1}`),
  };
}

export function normalizeOrientation(summary) {
  const source = firstRecord(summary?.orientation, summary?.orientation_evidence);
  if (!source) return null;
  const comparisons = asList(source.comparisons).length
    ? asList(source.comparisons)
    : asList(firstDefined(source, "groups") ?? summary?.orientation_groups);
  return {
    ...source,
    comparisons: comparisons.map(normalizeComparison),
    confidence_threshold: firstDefined(source, "confidence_threshold", "threshold"),
    forward_ln_likelihood_bdm: firstDefined(
      source,
      "forward_ln_likelihood_bdm",
      "forward_ln_likelihood",
    ),
    forward_minus_reverse_ln_likelihood_bdm: firstDefined(
      source,
      "forward_minus_reverse_ln_likelihood_bdm",
      "delta_ln_likelihood",
      "delta_lnL",
      "delta_log_likelihood",
    ),
    reconstructions: asList(source.reconstructions).map(normalizeHistory),
    reverse_ln_likelihood_bdm: firstDefined(
      source,
      "reverse_ln_likelihood_bdm",
      "reverse_ln_likelihood",
    ),
    selected_reconstructions: asList(source.selected_reconstructions).map(
      normalizeReconstructionRow,
    ),
    trees: asList(source.trees).map(normalizeTree),
  };
}
