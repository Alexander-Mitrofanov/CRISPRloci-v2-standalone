import { asList, asRecord } from "./values.js";

const text = (value) => (typeof value === "string" ? value : null);
const number = (value) => (typeof value === "number" && Number.isFinite(value) ? value : null);
const count = (value) => (Number.isInteger(value) && value >= 0 ? value : null);
const probability = (value) => (number(value) != null && value >= 0 && value <= 1 ? value : null);
const stage = (value) => {
  const record = asRecord(value) || {};
  return { status: text(record.status) || "not_reported", reason: text(record.reason) };
};
const pairs = (value, length) =>
  asList(value).filter(
    (row) =>
      Array.isArray(row) &&
      row.length === 3 &&
      count(row[0]) != null &&
      count(row[1]) != null &&
      row[0] < row[1] &&
      row[1] < length &&
      probability(row[2]) != null,
  );
const motif = (value, length) =>
  asList(value).filter(
    (row) =>
      Array.isArray(row) &&
      row.length === 2 &&
      count(row[0]) != null &&
      count(row[1]) != null &&
      row[0] < row[1] &&
      row[1] < length,
  );

function normalizeInstance(value, index) {
  const row = asRecord(value) || {};
  const rna = text(row.rna_sequence);
  const length = rna?.length || 0;
  const isolated = asRecord(row.isolated) || {};
  const context = asRecord(row.context) || {};
  const comparison = asRecord(row.comparison) || {};
  const references = asRecord(row.references) || {};
  return {
    id: text(row.id) || `preview_${index + 1}`,
    source_id: text(row.source_id),
    array_id: text(row.array_id),
    repeat_index: count(row.repeat_index),
    start: count(row.start),
    end: count(row.end),
    terminal: row.terminal === true,
    orientation: ["forward", "reverse"].includes(row.orientation) ? row.orientation : null,
    source_sequence: text(row.source_sequence),
    rna_sequence: rna,
    context_start: count(row.context_start),
    isolated: {
      ...stage(isolated),
      structure: text(isolated.structure),
      mfe_kcal_mol: number(isolated.mfe_kcal_mol),
      ensemble_kcal_mol: number(isolated.ensemble_kcal_mol),
      pairs: pairs(isolated.pairs, length),
      pairs_available:
        Array.isArray(isolated.pairs) &&
        pairs(isolated.pairs, length).length === isolated.pairs.length,
      unpaired: asList(isolated.unpaired).map(probability),
      motif_pairs: motif(isolated.motif_pairs, length),
    },
    context: {
      ...stage(context),
      intrarepeat_pairs: pairs(context.intrarepeat_pairs, length),
      pairs_available:
        Array.isArray(context.intrarepeat_pairs) &&
        pairs(context.intrarepeat_pairs, length).length === context.intrarepeat_pairs.length,
      unpaired_probabilities: asList(context.unpaired_probabilities).map(probability),
      effective_window: count(context.effective_window),
      effective_span: count(context.effective_span),
      competitor_count: count(context.competitor_count),
      competitors_truncated: context.competitors_truncated === true,
      competitors: asList(context.competitors).filter(
        (pair) =>
          count(pair?.i) != null &&
          count(pair?.j) != null &&
          pair.i < pair.j &&
          probability(pair.probability) != null,
      ),
    },
    comparison: {
      ...stage(comparison),
      motif_origin: text(comparison.motif_origin),
      motif_pairs: motif(comparison.motif_pairs, length),
      motif_pair_count: Array.isArray(comparison.motif_pairs)
        ? motif(comparison.motif_pairs, length).length
        : null,
      isolated_pair_support: probability(comparison.isolated_pair_support),
      context_pair_support: probability(comparison.context_pair_support),
      support_change: number(comparison.support_change),
    },
    references: {
      ...stage(references),
      release: text(references.release),
      reported_hit_count: count(references.reported_hit_count),
      total_hit_count: count(references.total_hit_count),
    },
  };
}

export function normalizeRepeats(value) {
  const report = asRecord(value) || {};
  const model = asRecord(report.model) || {};
  const counts = asRecord(report.counts) || {};
  const classification = asRecord(report.classification) || {};
  return {
    ...stage(report),
    molecule: text(report.molecule),
    model: {
      engine: text(model.engine),
      version: text(model.version),
      parameters: text(model.parameters),
      temperature: number(model.temperature),
      window: count(model.window),
      span: count(model.span),
    },
    counts: {
      sources: count(counts.sources),
      arrays: count(counts.arrays),
      instances: count(counts.instances),
      contexts: count(counts.contexts),
    },
    classification: {
      strand: stage(classification.strand),
      subtype: stage(classification.subtype),
      family: stage(classification.family),
    },
    interpretation: text(report.interpretation),
    instances: asList(report.instances).filter(asRecord).slice(0, 20).map(normalizeInstance),
    truncated:
      report.truncated === true ||
      asList(report.instances).length > 20 ||
      (count(counts.instances) != null && counts.instances > asList(report.instances).length),
  };
}

export function repeatStageText(value) {
  const reasons = {
    no_reference_pairs:
      "No reference pairs: the isolated MFE is unpaired. Folding still completed.",
    context_not_supplied: "Array context was not supplied.",
    ambiguous_bases: "Ambiguous bases prevent this folding analysis.",
    context_too_short: "The array context is too short for local folding.",
    motif_exceeds_context_span:
      "The motif exceeds the local pairing span; the comparison is unavailable.",
    reference_not_configured: "No repeat reference release was configured.",
    no_calibrated_strand_model: "No calibrated strand model is available.",
    no_validated_subtype_model: "No validated repeat subtype model is available.",
    no_validated_family_model: "No validated repeat family model is available.",
  };
  return (
    reasons[value?.reason] ||
    value?.reason?.replaceAll("_", " ") ||
    {
      completed: "Completed",
      no_hit: "No supported match in this reference release.",
      unsupported: "Unsupported",
      failed: "Analysis failed",
      truncated: "Limited results",
      not_requested: "Not requested",
      not_reported: "Not reported",
    }[value?.status] ||
    "Not reported"
  );
}
