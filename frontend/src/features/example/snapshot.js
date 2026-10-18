import {
  DNA_ONLY,
  EXAMPLE_SCHEMA_VERSION,
  FORBIDDEN_IDENTITY,
  FORBIDDEN_KEYS,
  MASKED_RECORD_ID,
  SHA256_HEX,
  fail,
  finite,
  object,
  positiveInteger,
} from "./contract.js";
import { validateEvolutionaryReconstructions } from "./reconstruction.js";
import { validateModuleSnapshot } from "./moduleSnapshot.js";

function scanPublicSnapshot(value) {
  if (typeof value === "string") {
    if ((value.length >= 40 && DNA_ONLY.test(value)) || FORBIDDEN_IDENTITY.test(value)) fail();
    return;
  }
  if (Array.isArray(value)) {
    value.forEach(scanPublicSnapshot);
    return;
  }
  if (!object(value)) return;
  for (const [key, item] of Object.entries(value)) {
    if (FORBIDDEN_KEYS.has(key.toLowerCase())) fail();
    scanPublicSnapshot(item);
  }
}

function validateRecord(record) {
  const allowedKeys = new Set([
    "record_id",
    "sequence_length",
    "source_array_orientation",
    "expected_spacer_count",
  ]);
  if (
    !object(record) ||
    Object.keys(record).some((key) => !allowedKeys.has(key)) ||
    !MASKED_RECORD_ID.test(record.record_id) ||
    !positiveInteger(record.sequence_length) ||
    !["pos", "neg"].includes(record.source_array_orientation) ||
    !positiveInteger(record.expected_spacer_count)
  )
    fail();
}

function validateComparison(comparison) {
  if (!object(comparison)) fail();
  for (const key of [
    "confidence_threshold",
    "forward_ln_likelihood_bdm",
    "reverse_ln_likelihood_bdm",
    "forward_minus_reverse_ln_likelihood_bdm",
  ]) {
    if (!finite(comparison[key])) fail();
  }
  const delta =
    Number(comparison.forward_ln_likelihood_bdm) - Number(comparison.reverse_ln_likelihood_bdm);
  if (Math.abs(delta - Number(comparison.forward_minus_reverse_ln_likelihood_bdm)) > 1e-9) fail();
  const expected =
    delta > Number(comparison.confidence_threshold)
      ? "Forward"
      : delta < -Number(comparison.confidence_threshold)
        ? "Reverse"
        : "ND";
  if (comparison.prediction !== expected) fail();
}

function validateTeachingClaims(example, job) {
  const findings = example.findings;
  const summary = job.summary;
  const detection = summary.detection;
  const adapter = summary.adapter;
  const comparisons = summary.orientation.comparisons;
  const reconstructions = summary.orientation.selected_reconstructions;
  const detectedSpacerCounts = detection.arrays.map((item) => Number(item.spacer_count));
  const detectedSources = new Set(detection.arrays.map((item) => String(item.source_id)));
  const observedSpacerRange = [
    Math.min(...detectedSpacerCounts),
    Math.max(...detectedSpacerCounts),
  ];
  if (
    findings.detection.arrays !== detection.array_count ||
    findings.detection.bona_fide !== (detection.category_counts?.["Bona-fide"] ?? 0) ||
    findings.detection.possible !== (detection.category_counts?.Possible ?? 0) ||
    findings.detection.spacer_count_range?.[0] !== observedSpacerRange[0] ||
    findings.detection.spacer_count_range?.[1] !== observedSpacerRange[1] ||
    findings.preflight.modeled_arrays !== adapter.emitted_array_count ||
    findings.preflight.eligible_groups !== adapter.emitted_group_count ||
    findings.preflight.excluded_arrays !== adapter.skipped_array_count ||
    !Array.isArray(comparisons) ||
    comparisons.length < 1 ||
    !Array.isArray(reconstructions) ||
    reconstructions.length < 1 ||
    !Array.isArray(adapter.groups)
  )
    fail();
  comparisons.forEach(validateComparison);
  const mainComparison = comparisons.find((item) => item.group === findings.orientation.group);
  const mainReconstruction = reconstructions.find(
    (item) => item.name === findings.reconstruction.group,
  );
  const mainGroup = adapter.groups.find((item) => item.name === findings.orientation.group);
  const delta = Number(mainComparison?.forward_minus_reverse_ln_likelihood_bdm);
  const threshold = Number(mainComparison?.confidence_threshold);
  const expectedDecision =
    delta > threshold
      ? "Input order supported"
      : delta < -threshold
        ? "Reverse input order supported"
        : "Unresolved";
  if (
    !mainComparison ||
    !mainReconstruction ||
    !mainGroup ||
    findings.reconstruction.group !== findings.orientation.group ||
    findings.orientation.delta_ln_likelihood !==
      mainComparison.forward_minus_reverse_ln_likelihood_bdm ||
    findings.orientation.confidence_threshold !== mainComparison.confidence_threshold ||
    findings.orientation.decision !== expectedDecision ||
    expectedDecision === "Unresolved" ||
    mainComparison.decisive !== true ||
    Math.abs(delta) <= threshold ||
    !Array.isArray(mainGroup.arrays) ||
    mainGroup.arrays.length !== mainGroup.array_count ||
    mainGroup.arrays.length !== findings.reconstruction.array_count ||
    mainGroup.arrays.some((item) => !detectedSources.has(String(item.source_id))) ||
    mainGroup.arrays_truncated !== false ||
    typeof mainGroup.repeat_key !== "string" ||
    !mainGroup.repeat_key ||
    findings.reconstruction.unique_spacers !== mainReconstruction["nb of unique spacers"] ||
    findings.reconstruction.insertions !== mainReconstruction["nb of reconstructed insertions"] ||
    findings.reconstruction.deletions !== mainReconstruction["nb of reconstructed deletions"] ||
    findings.reconstruction.duplications !==
      mainReconstruction["nb of reconstructed duplications"] ||
    findings.reconstruction.preferred_deletion_model !==
      mainReconstruction["Deletion model preferred by LRT"]
  )
    fail();
  validateEvolutionaryReconstructions(summary, mainGroup, mainComparison, mainReconstruction);
}

export function validateExampleSnapshot(value) {
  if (value?.schema?.name === "crisprloci.web-example") return validateModuleSnapshot(value);
  scanPublicSnapshot(value);
  const example = value?.example;
  const source = example?.source;
  const input = example?.input;
  const records = example?.records;
  const job = value?.job;
  const summary = job?.summary;
  if (
    value?.schema?.name !== "crispr-evolution-web-example" ||
    value?.schema?.version !== EXAMPLE_SCHEMA_VERSION ||
    !object(example) ||
    typeof example.analysis_question !== "string" ||
    typeof example.analysis_takeaway !== "string" ||
    !object(source) ||
    typeof source.masking_policy !== "string" ||
    typeof source.provenance_note !== "string" ||
    !object(input) ||
    !Array.isArray(records) ||
    records.length < 3 ||
    records.length !== input.record_count ||
    !positiveInteger(input.base_count) ||
    typeof input.filename !== "string" ||
    !/\.fasta$/i.test(input.filename) ||
    !SHA256_HEX.test(input.file_sha256) ||
    !SHA256_HEX.test(input.normalized_sha256) ||
    source.displayed_locus_bases !== input.base_count ||
    job?.status !== "completed" ||
    job?.mode !== "orientation" ||
    !object(summary?.detection) ||
    !object(summary?.adapter) ||
    !object(summary?.orientation) ||
    summary.pipeline_status !== "completed" ||
    !Array.isArray(job.artifacts) ||
    job.artifacts.length !== 0 ||
    job.options?.category_policy !== "bona_fide_possible" ||
    job.options?.spacer_distance !== 1 ||
    job.options?.bias_corrections_requested !== true ||
    job.options?.bias_corrections_effective !== true
  )
    fail();
  records.forEach(validateRecord);
  if (
    new Set(records.map((record) => record.record_id)).size !== records.length ||
    records.reduce((total, record) => total + record.sequence_length, 0) !== input.base_count
  )
    fail();
  validateTeachingClaims(example, job);
  return value;
}
