import { selectedToolOptions } from "./toolOptions.js";

export function buildSubmission({ sequence, filename, mode, options }) {
  const toolOptions = selectedToolOptions(mode, options.toolOptions);
  return {
    sequence,
    filename,
    mode,
    ...(Object.keys(toolOptions).length ? { tool_options: toolOptions } : {}),
    ...(mode === "viral_search" ? { viral_max_mismatches: options.viralMaxMismatches ?? 2 } : {}),
    ...(mode === "spacer_association"
      ? { association_grouping: options.associationGrouping || null }
      : {}),
    ...(["repeats", "rna_family"].includes(mode) ? { molecule: options.molecule || "DNA" } : {}),
    ...(["loci", "tracrrna"].includes(mode)
      ? { tracr_model_type: options.tracrModelType || "II" }
      : {}),
    ...(["loci", "leader"].includes(mode)
      ? { leader_flank_length: options.leaderFlankLength ?? 500 }
      : {}),
    category_policy: options.categoryPolicy,
    spacer_distance: options.spacerEditDistance,
    bias_corrections: options.biasCorrection,
  };
}
