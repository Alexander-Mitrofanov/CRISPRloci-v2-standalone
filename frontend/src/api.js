import publicApiContract from "./contracts/public-api-v1.json";

const DEFAULT_API_BASE = String(import.meta.env.VITE_API_BASE_URL || "").trim();
const PRIVATE_REQUEST_POLICY = Object.freeze({
  cache: "no-store",
  credentials: "omit",
  referrerPolicy: "no-referrer",
});

export class ApiError extends Error {
  constructor(message, status = 0, code = "request_failed") {
    super(message);
    this.name = "ApiError";
    this.status = status;
    this.code = code;
  }
}

export function normalizeApiBase(value) {
  const base = String(value || "")
    .trim()
    .replace(/\/+$/, "");
  if (!base) return "";
  let parsed;
  try {
    parsed = new URL(base);
  } catch {
    throw new ApiError("The analysis API address is not a valid URL.", 0, "invalid_api_base");
  }
  // URL.hostname retains brackets around IPv6 literals in browsers and Node.
  const local = ["localhost", "127.0.0.1", "[::1]"].includes(parsed.hostname);
  if (
    (parsed.protocol !== "https:" && !(local && parsed.protocol === "http:")) ||
    parsed.username ||
    parsed.password ||
    parsed.pathname !== "/" ||
    parsed.search ||
    parsed.hash
  ) {
    throw new ApiError("The analysis API must use HTTPS.", 0, "insecure_api_base");
  }
  return parsed.origin;
}

function createUrl(base, path) {
  const normalized = normalizeApiBase(base);
  if (!normalized) {
    throw new ApiError(
      "The site operator has not configured the analysis API.",
      0,
      "api_not_configured",
    );
  }
  return `${normalized}${path.startsWith("/") ? path : `/${path}`}`;
}

function contractPath(name, values = {}) {
  let path = publicApiContract.routes[name].path;
  for (const [key, value] of Object.entries(values)) {
    path = path.replace(`{${key}}`, encodeURIComponent(value));
  }
  if (path.includes("{")) {
    throw new ApiError("The analysis API route is incomplete.", 0, "invalid_api_route");
  }
  return path;
}

async function parseJsonResponse(response) {
  let payload = null;
  try {
    payload = await response.json();
  } catch {
    // A controlled status message is safer and clearer than a JSON parser detail.
  }
  if (!response.ok) {
    const detail = payload?.detail;
    const message = typeof detail === "string" ? detail : detail?.message;
    throw new ApiError(
      message || `The analysis service returned status ${response.status}.`,
      response.status,
      detail?.code,
    );
  }
  return payload;
}

async function boundedText(response, maxBytes) {
  const tooLarge = () => new ApiError("This sequence file exceeds the input byte limit.");
  if (Number(response.headers.get("content-length")) > maxBytes) {
    await response.body?.cancel();
    throw tooLarge();
  }
  if (!response.body) throw new ApiError("The sequence file is empty.");
  const reader = response.body.getReader();
  const decoder = new TextDecoder("utf-8", { fatal: true });
  const parts = [];
  let bytes = 0;
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      bytes += value.byteLength;
      if (bytes > maxBytes) throw tooLarge();
      parts.push(decoder.decode(value, { stream: true }));
    }
    parts.push(decoder.decode());
    return parts.join("");
  } finally {
    await reader.cancel();
    reader.releaseLock();
  }
}

export function createApiClient(baseUrl = DEFAULT_API_BASE, fetchImpl = globalThis.fetch) {
  const endpoint = (path) => createUrl(baseUrl, path);
  const jobHeaders = (accessToken, extra = {}) => ({
    Accept: "application/json",
    Authorization: `Bearer ${accessToken}`,
    ...extra,
  });

  return {
    configured: Boolean(String(baseUrl || "").trim()),
    displayBase: String(baseUrl || "")
      .trim()
      .replace(/\/+$/, ""),

    async catalogSummary({ signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("catalog_summary")), {
          ...PRIVATE_REQUEST_POLICY,
          headers: { Accept: "application/json" },
          signal,
        }),
      );
    },

    async catalogPage(entity, filters = {}, { signal } = {}) {
      const routes = {
        assemblies: "catalog_assemblies",
        arrays: "catalog_arrays",
        repeats: "catalog_repeats",
        spacers: "catalog_spacers",
        cas_genes: "catalog_cas_genes",
        cas_systems: "catalog_cas_systems",
      };
      if (!routes[entity])
        throw new ApiError("Unknown database record type.", 0, "invalid_catalog_entity");
      const query = new URLSearchParams();
      for (const key of ["q", "assembly_accession", "array_id", "cursor", "limit", "match"]) {
        if (filters[key] !== undefined && filters[key] !== null && filters[key] !== "") {
          query.set(key, String(filters[key]));
        }
      }
      return parseJsonResponse(
        await fetchImpl(endpoint(`${contractPath(routes[entity])}?${query}`), {
          ...PRIVATE_REQUEST_POLICY,
          headers: { Accept: "application/json" },
          signal,
        }),
      );
    },

    async catalogArray(arrayId, { signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("catalog_array", { array_id: arrayId })), {
          ...PRIVATE_REQUEST_POLICY,
          headers: { Accept: "application/json" },
          signal,
        }),
      );
    },

    async health({ signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("health")), {
          ...PRIVATE_REQUEST_POLICY,
          headers: { Accept: "application/json" },
          signal,
        }),
      );
    },

    async config({ signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("config")), {
          ...PRIVATE_REQUEST_POLICY,
          headers: { Accept: "application/json" },
          signal,
        }),
      );
    },

    async submit(payload, { signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("submit_job")), {
          ...PRIVATE_REQUEST_POLICY,
          method: publicApiContract.routes.submit_job.method,
          headers: { Accept: "application/json", "Content-Type": "application/json" },
          body: JSON.stringify(payload),
          signal,
        }),
      );
    },

    async getJob(jobId, accessToken, { signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("job_status", { job_id: jobId })), {
          ...PRIVATE_REQUEST_POLICY,
          headers: jobHeaders(accessToken),
          signal,
        }),
      );
    },

    async cancelJob(jobId, accessToken, { signal } = {}) {
      return parseJsonResponse(
        await fetchImpl(endpoint(contractPath("cancel_job", { job_id: jobId })), {
          ...PRIVATE_REQUEST_POLICY,
          method: publicApiContract.routes.cancel_job.method,
          headers: jobHeaders(accessToken),
          signal,
        }),
      );
    },

    async downloadBundle(jobId, accessToken, { signal } = {}) {
      const response = await fetchImpl(endpoint(contractPath("result_bundle", { job_id: jobId })), {
        ...PRIVATE_REQUEST_POLICY,
        headers: jobHeaders(accessToken, { Accept: "application/zip, application/octet-stream" }),
        signal,
      });
      if (!response.ok) return parseJsonResponse(response);
      return response.blob();
    },

    async downloadArtifact(jobId, artifactId, accessToken, { signal } = {}) {
      const response = await fetchImpl(
        endpoint(contractPath("artifact", { job_id: jobId, artifact_id: artifactId })),
        {
          ...PRIVATE_REQUEST_POLICY,
          headers: jobHeaders(accessToken, { Accept: "application/octet-stream" }),
          signal,
        },
      );
      if (!response.ok) return parseJsonResponse(response);
      return response.blob();
    },

    async readSequenceArtifact(
      jobId,
      artifactId,
      accessToken,
      { signal, maxBytes = 20_000_000 } = {},
    ) {
      const response = await fetchImpl(
        endpoint(contractPath("artifact", { job_id: jobId, artifact_id: artifactId })),
        {
          ...PRIVATE_REQUEST_POLICY,
          headers: jobHeaders(accessToken),
          signal,
        },
      );
      if (!response.ok) return parseJsonResponse(response);
      return boundedText(response, Math.min(maxBytes, 20_000_000));
    },
  };
}

export const api = createApiClient();
