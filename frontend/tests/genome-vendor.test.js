// @vitest-environment node
import { mkdtempSync, mkdirSync, writeFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { auditIgvVendor, verifiedIgvBundle } from "../scripts/igv-vendor.mjs";
import { hasForbiddenPublicIdentity } from "../scripts/public-identity-policy.mjs";

describe("genome vendor public-export boundary", () => {
  it("rejects a vendor chunk containing application modules", () => {
    expect(() =>
      auditIgvVendor().generateBundle(
        {},
        {
          vendor: {
            type: "chunk",
            modules: { "/node_modules/igv/dist/igv.esm.js": {}, "/src/app.js": {} },
          },
        },
      ),
    ).toThrow(/isolated vendor/);
  });

  it("accepts only the exact audited output and rejects post-build changes", () => {
    const root = mkdtempSync(join(tmpdir(), "genome-vendor-test-"));
    try {
      const fileName = "assets/igv.esm-example.js";
      const code = 'const builtInReference="GCF_000001405.40";';
      let manifest;
      auditIgvVendor().generateBundle.call(
        {
          emitFile: (value) => {
            manifest = value;
          },
        },
        {},
        {
          vendor: {
            type: "chunk",
            modules: { "/node_modules/igv/dist/igv.esm.js": {} },
            fileName,
            code,
          },
        },
      );
      mkdirSync(join(root, "assets"));
      writeFileSync(join(root, fileName), code);
      writeFileSync(join(root, manifest.fileName), manifest.source);
      expect(verifiedIgvBundle(root)).toBe(join(root, fileName));
      expect(hasForbiddenPublicIdentity(code, { compiledBundle: true })).toBe(true);
      writeFileSync(join(root, fileName), code + 'const source="CP000001";');
      expect(() => verifiedIgvBundle(root)).toThrow(/changed after audit/);
    } finally {
      rmSync(root, { recursive: true, force: true });
    }
  });
});
