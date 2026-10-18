import { onBeforeUnmount, onMounted, ref } from "vue";

import { api } from "../api.js";

const DEFAULT_LIMITS = Object.freeze({
  maxBases: 0,
  maxRecordBases: 0,
  maxRecords: 0,
  maxRequestBytes: 0,
  maxArchiveBytes: 0,
  maxHeaderCharacters: 200,
});
const RETRY_INTERVAL_MS = 10_000;

export function useServiceConfig(client = api) {
  const service = ref({ state: "checking", message: "Checking analysis service" });
  const limits = ref({ ...DEFAULT_LIMITS });
  let controller;
  let retryTimer;

  async function refresh() {
    window.clearTimeout(retryTimer);
    controller?.abort();
    if (!client.configured) {
      service.value = { state: "offline", message: "Analysis endpoint not configured" };
      return;
    }
    const request = new AbortController();
    controller = request;
    try {
      const [health, config] = await Promise.all([
        client.health({ signal: request.signal }),
        client.config({ signal: request.signal }),
      ]);
      if (request.signal.aborted) return;
      service.value = {
        state: "online",
        message: "Analysis service ready",
        version: health?.version || config?.api_version,
        modes: Array.isArray(config?.modes)
          ? config.modes.filter((mode) => typeof mode === "string")
          : [],
        expiresHours:
          config?.retention_hours ||
          health?.retention_hours ||
          (config?.retention_seconds ? Math.round(config.retention_seconds / 3600) : null),
      };
      limits.value = {
        maxBases:
          config?.max_total_bases || config?.max_sequence_bases || health?.max_sequence_bases || 0,
        maxRecordBases: config?.max_record_bases || 0,
        maxRecords: config?.max_records || health?.max_records || 0,
        maxRequestBytes: config?.max_request_bytes || 0,
        maxArchiveBytes: config?.max_archive_bytes || 0,
        maxHeaderCharacters: config?.max_header_characters || 200,
      };
    } catch (error) {
      if (!request.signal.aborted && error.name !== "AbortError") {
        // Cancel the other read before scheduling another pair of requests.
        request.abort();
        service.value = {
          state: "offline",
          message: error.message || "The analysis API could not be reached.",
        };
        retryTimer = window.setTimeout(refresh, RETRY_INTERVAL_MS);
      }
    } finally {
      if (controller === request) controller = undefined;
    }
  }

  function retryWhenAvailable() {
    if (service.value.state === "offline" && !controller) void refresh();
  }

  window.addEventListener("focus", retryWhenAvailable);
  window.addEventListener("online", retryWhenAvailable);
  onMounted(refresh);
  onBeforeUnmount(() => {
    window.clearTimeout(retryTimer);
    controller?.abort();
    window.removeEventListener("focus", retryWhenAvailable);
    window.removeEventListener("online", retryWhenAvailable);
  });
  return { service, limits, refresh };
}
