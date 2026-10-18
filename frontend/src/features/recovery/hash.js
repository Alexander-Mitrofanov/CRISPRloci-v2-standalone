import { normalizeJobCredential } from "./credential.js";

const RECOVERY_PREFIX = "#job=";
const RECOVERY_PAYLOAD_PATTERN = /^([0-9a-f]{32})\.([A-Za-z0-9_-]{43})$/;

export function isRecoveryHash(value) {
  return typeof value === "string" && value.startsWith(RECOVERY_PREFIX);
}

export function serializeRecoveryHash(value) {
  const credential = normalizeJobCredential(value, 0);
  return `${RECOVERY_PREFIX}${credential.jobId}.${credential.accessToken}`;
}

export function parseRecoveryHash(value) {
  if (!isRecoveryHash(value))
    throw new Error("Recovery link does not contain a supported job hash.");
  const match = RECOVERY_PAYLOAD_PATTERN.exec(value.slice(RECOVERY_PREFIX.length));
  if (!match) throw new Error("Recovery link contains an invalid job hash.");
  return normalizeJobCredential({ jobId: match[1], accessToken: match[2], expiresAt: null });
}

export function buildRecoveryUrl(value, currentUrl = window.location.href) {
  const url = new URL(currentUrl);
  url.hash = serializeRecoveryHash(value).slice(1);
  return url.href;
}
