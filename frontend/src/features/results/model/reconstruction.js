import { asList, firstDefined, firstRecord } from "./values.js";

const METRIC_ALIASES = {
  acquisitions: [
    "acquisitions",
    "nb of reconstructed insertions",
    "gains",
    "insertions",
    "gain_events",
  ],
  aligned_spacers: [
    "nb of spacers in alignment",
    "nb of spacers in model matrix",
    "alignment_spacers",
    "aligned_spacers",
  ],
  bdm_log_likelihood: ["bdm_log_likelihood", "ln_lh_bdm", "log_likelihood", "ln_likelihood", "lnL"],
  deletion_rate_bdm: ["deletion_rate_bdm", "deletion_rate", "loss_rate"],
  deletions: ["nb of reconstructed deletions", "deletions", "losses", "deletion_events"],
  duplications: ["nb of reconstructed duplications", "duplications"],
  idm_log_likelihood: ["ln_lh_idm", "idm_log_likelihood"],
  independent_gains: ["nb of reconstructed independent gains", "independent_gains"],
  leaf_count: ["nb of leafs (after combining non-uniques)", "leaf_count"],
  likelihood_ratio_statistic: ["test_statistic (-2*ln_lh_ratio)", "likelihood_ratio_statistic"],
  model_selection_cutoff: ["chi2_quantile", "model_selection_cutoff"],
  preferred_model: ["Deletion model preferred by LRT", "preferred_model", "model_name", "model"],
  reacquisitions: ["nb of reconstructed reacquisitions", "reacquisitions"],
  rearrangements: ["nb of reconstructed rearrangements", "rearrangements"],
  runtime_seconds: ["run_time", "runtime_seconds", "duration_seconds"],
  unique_arrays: ["nb of unique spacer arrays", "unique_arrays"],
  unique_spacers: ["nb of unique spacers", "unique_spacers"],
};

export function normalizeReconstructionRow(value, index = 0) {
  const row = firstRecord(value) || {};
  const normalized = {
    ...row,
    name: String(firstDefined(row, "name", "group", "group_id", "id") || `group_${index + 1}`),
  };
  for (const [field, aliases] of Object.entries(METRIC_ALIASES)) {
    normalized[field] = firstDefined(row, ...aliases);
  }
  return normalized;
}

export function normalizeTree(value, index = 0) {
  const tree = firstRecord(value) || {};
  return {
    ...tree,
    group: String(firstDefined(tree, "group", "name", "group_id", "id") || `group_${index + 1}`),
  };
}

export function normalizeReconstruction(summary) {
  const source = firstRecord(summary?.reconstruction, summary?.spacerplacer);
  if (!source) return null;
  return {
    ...source,
    results: asList(source.results).map(normalizeReconstructionRow),
    selected_model: source.selected_model
      ? normalizeReconstructionRow(source.selected_model)
      : null,
    trees: asList(source.trees).map(normalizeTree),
  };
}
