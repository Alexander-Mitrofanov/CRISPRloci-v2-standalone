import contract from "../../contracts/public-api-v1.json" with { type: "json" };
import { fail, object, positiveInteger, SHA256_HEX, FORBIDDEN_KEYS } from "./contract.js";

// Scientific result rows legitimately include sequences and null evidence.
// Reject capabilities and local paths, while retaining exact scientific output.
function scan(value) {
  if (typeof value === "string") {
    if (
      /(?:^|[^A-Za-z0-9])(?:CP|FR|LN|LR|AP)\d{6}/i.test(value) ||
      /(?:\/home\/|\/tmp\/|\.ops-private\/|\.runtime\/)/.test(value)
    )
      fail();
  } else if (Array.isArray(value)) value.forEach(scan);
  else if (object(value))
    for (const [key, row] of Object.entries(value)) {
      if (FORBIDDEN_KEYS.has(key.toLowerCase()) && key !== "accession") fail();
      if (key === "accession" && row !== null && !/^SYNTHETIC_\d+\.\d+$/.test(row)) fail();
      scan(row);
    }
}

export function validateModuleSnapshot(value) {
  scan(value);
  const { example, job } = value;
  const input = example?.input;
  const verification = example?.verification;
  const mode = job?.mode;
  const section =
    {
      cas: "casandra",
      leader: "crisprleader",
      loci: "casandra",
      detection: "detection",
      repeat_context: "repeats",
    }[mode] || mode;
  if (
    value.schema?.version !== "1.0.0" ||
    !contract.enums.analysis_modes.includes(mode) ||
    job.status !== "completed" ||
    job.summary?.pipeline_status !== "completed" ||
    !object(job.summary?.[section]) ||
    !object(job.options) ||
    !object(example?.submission_options) ||
    typeof example?.source_note !== "string" ||
    !object(input) ||
    !positiveInteger(input.record_count) ||
    !positiveInteger(input.base_count) ||
    !["input.fasta", "input.json"].includes(input.filename) ||
    !SHA256_HEX.test(input.file_sha256) ||
    !SHA256_HEX.test(input.normalized_sha256) ||
    input.normalized_sha256 !== job.input?.sha256 ||
    input.record_count !== job.input?.record_count ||
    input.base_count !== job.input?.base_count ||
    verification?.independent_runs !== 2 ||
    verification.scientific_results_identical !== true ||
    !SHA256_HEX.test(verification.release_id) ||
    !Array.isArray(job.artifacts) ||
    job.artifacts.length
  )
    fail();
  return value;
}
