import { finiteMetric } from "./formatting.js";

export const reconstructionMetric = (row, name) => {
  const field = {
    acquisitions: "acquisitions",
    aligned: "aligned_spacers",
    bdmLikelihood: "bdm_log_likelihood",
    cutoff: "model_selection_cutoff",
    deletions: "deletions",
    duplications: "duplications",
    idmLikelihood: "idm_log_likelihood",
    independentGains: "independent_gains",
    leaves: "leaf_count",
    lrt: "likelihood_ratio_statistic",
    patterns: "unique_arrays",
    reacquisitions: "reacquisitions",
    rearrangements: "rearrangements",
    unique: "unique_spacers",
  }[name];
  return finiteMetric(field ? row?.[field] : null);
};

export const preferredDeletionModel = (row) => String(row?.preferred_model || "Not reported");
