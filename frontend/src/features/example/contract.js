export const EXAMPLE_FASTA_PATH = "example-input.fasta";
export const EXAMPLE_RESULT_PATH = "example-result.json";
export const EXAMPLE_SCHEMA_VERSION = "1.3.0";

export const SHA256_HEX = /^[0-9a-f]{64}$/;
export const MASKED_RECORD_ID = /^example_record_\d{2}$/;
export const FORBIDDEN_KEYS = new Set([
  "access_token",
  "artifact_url",
  "download_url",
  "job_id",
  "organism",
  "strain",
  "accession",
  "ncbi_url",
  "region_start_1based",
  "region_end_1based",
  "token",
  "token_digest",
]);
export const DNA_ONLY = /^[ACGTRYSWKMBDHVN]+$/i;
export const FORBIDDEN_IDENTITY = /(?:CP|FR|LN|LR|AP)\d{6}/i;

export function fail(
  message = "The example result is incomplete or incompatible with this interface.",
) {
  throw new Error(message);
}

export function object(value) {
  return value != null && typeof value === "object" && !Array.isArray(value);
}

export function finite(value) {
  return Number.isFinite(Number(value));
}

export function positiveInteger(value) {
  return Number.isInteger(Number(value)) && Number(value) > 0;
}

export function canonicalInteger(value) {
  return typeof value === "number" && Number.isInteger(value) && value > 0;
}

export function nonNegativeInteger(value) {
  return typeof value === "number" && Number.isInteger(value) && value >= 0;
}
