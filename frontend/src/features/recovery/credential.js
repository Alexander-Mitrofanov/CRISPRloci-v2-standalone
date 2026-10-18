const JOB_ID_PATTERN = /^[0-9a-f]{32}$/;
const TOKEN_PATTERN = /^[A-Za-z0-9_-]{43}$/;

export function normalizeJobCredential(value, now = Date.now()) {
  if (!value || typeof value !== "object" || Array.isArray(value)) {
    throw new Error("Recovery credential must be an object.");
  }
  const expectedKeys = ["accessToken", "expiresAt", "jobId"];
  const actualKeys = Object.keys(value).sort();
  if (
    actualKeys.length !== expectedKeys.length ||
    actualKeys.some((key, index) => key !== expectedKeys[index])
  ) {
    throw new Error("Recovery credential must contain exactly the supported fields.");
  }
  const { jobId, accessToken, expiresAt } = value;
  if (
    typeof jobId !== "string" ||
    typeof accessToken !== "string" ||
    !JOB_ID_PATTERN.test(jobId) ||
    !TOKEN_PATTERN.test(accessToken)
  ) {
    throw new Error("Recovery credential has an invalid job ID or access token.");
  }
  if (expiresAt != null) {
    const expiry = Date.parse(String(expiresAt));
    if (!Number.isFinite(expiry)) throw new Error("Recovery credential has an invalid expiry.");
    if (expiry <= now) throw new Error("Recovery credential has expired.");
  }
  return { jobId, accessToken, expiresAt: expiresAt == null ? null : String(expiresAt) };
}
