import { onBeforeUnmount, ref, watch } from "vue";

import { api } from "../api.js";
import {
  isRecoveryHash,
  installRecoverySectionNavigation,
  normalizeJobCredential,
  readBrowserRecovery,
  replaceBrowserRecovery,
} from "../features/recovery/index.js";
import { TERMINAL_STATUSES } from "../science.js";
import { revealSection } from "../utils/dom.js";

export function useJobSession(client = api, browser = window) {
  const recovered = readBrowserRecovery(browser);
  const credential = ref(recovered.credential);
  const job = ref(null);
  const exampleSnapshot = ref(null);
  const pollError = ref(recovered.error);
  const cancelling = ref(false);
  let cancellingCredential;
  let pollTimer;
  let pollController;
  let revealTimer;

  function stopPolling() {
    window.clearTimeout(pollTimer);
    pollController?.abort();
    pollController = undefined;
  }

  async function poll() {
    const current = credential.value;
    if (!current) return;
    pollController = new AbortController();
    try {
      const latest = await client.getJob(current.jobId, current.accessToken, {
        signal: pollController.signal,
      });
      if (credential.value !== current) return;
      job.value = latest;
      pollError.value = "";
      if (latest.expires_at && latest.expires_at !== current.expiresAt) {
        credential.value = normalizeJobCredential({ ...current, expiresAt: latest.expires_at });
        return;
      }
      if (!TERMINAL_STATUSES.has(latest.status)) pollTimer = window.setTimeout(poll, 2_500);
    } catch (error) {
      if (error.name === "AbortError" || credential.value !== current) return;
      if ([401, 403, 404, 410].includes(error.status)) {
        credential.value = null;
        job.value = null;
      }
      pollError.value = error.message || "Job status could not be refreshed.";
      if (![401, 403, 404, 410].includes(error.status)) pollTimer = window.setTimeout(poll, 5_000);
    }
  }

  watch(
    credential,
    (next) => {
      stopPolling();
      if (cancellingCredential && !sameJobCapability(cancellingCredential, next)) {
        cancellingCredential = undefined;
        cancelling.value = false;
      }
      replaceBrowserRecovery(next, browser);
      if (next) void poll();
    },
    { flush: "sync", immediate: true },
  );

  function onHashChange() {
    if (isRecoveryHash(browser.location.hash)) {
      const next = readBrowserRecovery(browser);
      if (next.credential) {
        credential.value = next.credential;
        job.value = null;
        pollError.value = "";
      } else {
        credential.value = null;
        job.value = null;
        pollError.value = next.error;
        replaceBrowserRecovery(null, browser);
      }
      return;
    }
    if (credential.value) replaceBrowserRecovery(credential.value, browser);
  }

  browser.addEventListener("hashchange", onHashChange);
  const removeRecoveryNavigation = installRecoverySectionNavigation(
    () => Boolean(credential.value),
    browser,
  );

  function onSubmitted(nextCredential, initialJob) {
    exampleSnapshot.value = null;
    credential.value = nextCredential;
    job.value = initialJob;
    pollError.value = "";
    window.clearTimeout(revealTimer);
    revealTimer = window.setTimeout(() => revealSection("job-status", "#job-heading"), 50);
  }

  function onExampleLoaded(snapshot) {
    exampleSnapshot.value = snapshot;
    window.clearTimeout(revealTimer);
    if (snapshot)
      revealTimer = window.setTimeout(
        () => revealSection("example-result", "#results-heading"),
        50,
      );
  }

  async function cancel() {
    const current = credential.value;
    if (!current || sameJobCapability(cancellingCredential, current)) return;
    cancellingCredential = current;
    cancelling.value = true;
    pollError.value = "";
    try {
      const response = await client.cancelJob(current.jobId, current.accessToken);
      if (!sameJobCapability(credential.value, current)) return;
      job.value = response?.job || response;
    } catch (error) {
      if (!sameJobCapability(credential.value, current)) return;
      pollError.value = error.message || "The cancellation request failed.";
    } finally {
      if (sameJobCapability(cancellingCredential, current)) {
        cancellingCredential = undefined;
        cancelling.value = false;
      }
    }
  }

  function forget() {
    window.clearTimeout(revealTimer);
    credential.value = null;
    job.value = null;
    pollError.value = "";
  }

  onBeforeUnmount(() => {
    window.clearTimeout(revealTimer);
    stopPolling();
    browser.removeEventListener("hashchange", onHashChange);
    removeRecoveryNavigation();
  });
  return {
    credential,
    job,
    exampleSnapshot,
    pollError,
    cancelling,
    onSubmitted,
    onExampleLoaded,
    cancel,
    forget,
  };
}

function sameJobCapability(left, right) {
  return Boolean(
    left && right && left.jobId === right.jobId && left.accessToken === right.accessToken,
  );
}
