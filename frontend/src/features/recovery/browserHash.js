import { isRecoveryHash, parseRecoveryHash, serializeRecoveryHash } from "./hash.js";

export function readBrowserRecovery(browser = window) {
  if (!isRecoveryHash(browser.location.hash)) return { credential: null, error: "" };
  try {
    return { credential: parseRecoveryHash(browser.location.hash), error: "" };
  } catch {
    return { credential: null, error: "This recovery link is invalid or no longer supported." };
  }
}

export function replaceBrowserRecovery(credential, browser = window) {
  const url = new URL(browser.location.href);
  if (credential) url.hash = serializeRecoveryHash(credential).slice(1);
  else if (isRecoveryHash(url.hash)) url.hash = "";
  browser.history.replaceState(browser.history.state, "", url);
}
