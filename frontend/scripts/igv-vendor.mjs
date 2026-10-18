import { createHash } from "node:crypto";
import { readFileSync, existsSync } from "node:fs";
import { resolve } from "node:path";

const sourcePath = new URL("../node_modules/igv/dist/igv.esm.js", import.meta.url);
const sourceHash = "baa5bf1bd9fb4c414b506e5a2fe383cfa5b8ad19e55e44b7b01bec5dec9122fc";
const manifestPath = "igv-vendor.json";
const hash = (value) => createHash("sha256").update(value).digest("hex");

// IGV contains public reference-assembly identifiers. Exempt only its isolated,
// checksum-pinned vendor chunk from the private-example identity scanner, never
// application code or data. A dependency upgrade requires explicit re-audit.
export function auditIgvVendor() {
  return {
    name: "audit-igv-vendor",
    apply: "build",
    generateBundle(_options, bundle) {
      if (hash(readFileSync(sourcePath)) !== sourceHash)
        throw new Error("Unreviewed IGV distribution");
      const chunks = Object.values(bundle).filter(
        (chunk) =>
          chunk.type === "chunk" &&
          Object.keys(chunk.modules).some((id) => id.endsWith("/node_modules/igv/dist/igv.esm.js")),
      );
      if (chunks.length !== 1 || Object.keys(chunks[0].modules).length !== 1) {
        throw new Error("IGV must remain in an isolated vendor chunk");
      }
      const chunk = chunks[0];
      this.emitFile({
        type: "asset",
        fileName: manifestPath,
        source: JSON.stringify({
          package: "igv",
          version: "3.8.10",
          source_sha256: sourceHash,
          file: chunk.fileName,
          sha256: hash(chunk.code),
        }),
      });
    },
  };
}

export function verifiedIgvBundle(root) {
  const manifest = resolve(root, manifestPath);
  if (!existsSync(manifest)) return null;
  const info = JSON.parse(readFileSync(manifest, "utf8"));
  if (
    info.package !== "igv" ||
    info.version !== "3.8.10" ||
    info.source_sha256 !== sourceHash ||
    !/^assets\/igv(?:\.esm)?-[A-Za-z0-9_-]+\.js$/.test(info.file) ||
    hash(readFileSync(sourcePath)) !== sourceHash
  ) {
    throw new Error("Invalid IGV vendor audit");
  }
  const file = resolve(root, info.file);
  if (hash(readFileSync(file)) !== info.sha256)
    throw new Error("IGV vendor bundle changed after audit");
  return file;
}
