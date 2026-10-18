import { defineConfig, loadEnv } from "vite";
import vue from "@vitejs/plugin-vue";
import { igvStyleHashes } from "./scripts/igv-csp.mjs";
import { auditIgvVendor } from "./scripts/igv-vendor.mjs";

function normalizeBase(value) {
  const trimmed = String(value || "/").trim();
  if (!trimmed || trimmed === "/") return "/";
  return `/${trimmed.replace(/^\/+|\/+$/g, "")}/`;
}

function productionCsp(command, value) {
  if (command !== "build") return null;
  const raw = String(value || "").trim();
  if (!raw) throw new Error("VITE_API_BASE_URL is required for a production build.");

  let parsed;
  try {
    parsed = new URL(raw);
  } catch {
    throw new Error("VITE_API_BASE_URL must be a valid HTTPS origin.");
  }
  if (
    parsed.protocol !== "https:" ||
    parsed.username ||
    parsed.password ||
    parsed.pathname !== "/" ||
    parsed.search ||
    parsed.hash
  ) {
    throw new Error(
      "VITE_API_BASE_URL must be an exact HTTPS origin without credentials, path, query, or fragment.",
    );
  }

  const policy = [
    "default-src 'none'",
    "base-uri 'none'",
    "object-src 'none'",
    "script-src 'self'",
    "script-src-attr 'none'",
    `style-src 'self' ${igvStyleHashes.join(" ")}`,
    "style-src-attr 'none'",
    "img-src 'self' data:",
    "font-src 'self'",
    `connect-src 'self' ${parsed.origin}`,
    "media-src 'none'",
    "frame-src 'none'",
    "child-src 'none'",
    "worker-src 'none'",
    "manifest-src 'self'",
    "form-action 'none'",
    "upgrade-insecure-requests",
  ].join("; ");

  return {
    name: "production-content-security-policy",
    enforce: "post",
    transformIndexHtml() {
      return [
        {
          tag: "meta",
          attrs: { "http-equiv": "Content-Security-Policy", content: policy },
          injectTo: "head-prepend",
        },
      ];
    },
  };
}

export default defineConfig(({ command, mode }) => {
  const env = loadEnv(mode, process.cwd(), "VITE_");
  const csp = productionCsp(command, process.env.VITE_API_BASE_URL || env.VITE_API_BASE_URL);
  return {
    base: normalizeBase(process.env.VITE_BASE_PATH || env.VITE_BASE_PATH),
    define:
      command === "serve" &&
      mode === "development" &&
      !env.VITE_API_BASE_URL &&
      !process.env.VITE_API_BASE_URL
        ? { "import.meta.env.VITE_API_BASE_URL": JSON.stringify("http://127.0.0.1:8000") }
        : {},
    plugins: [vue(), auditIgvVendor(), ...(csp ? [csp] : [])],
    test: {
      environment: "jsdom",
      globals: true,
      setupFiles: "./tests/setup.js",
      css: true,
    },
  };
});
