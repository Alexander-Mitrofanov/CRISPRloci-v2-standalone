import { describe, expect, it } from "vitest";
import { hasForbiddenPublicIdentity } from "../scripts/public-identity-policy.mjs";

describe("public export identity boundary", () => {
  it("preserves CCTK deletion-slot identifiers without confusing them with AP accessions", () => {
    expect(hasForbiddenPublicIdentity('{"id":"gap000003"}')).toBe(false);
    expect(hasForbiddenPublicIdentity('{"id":"source_AP009389.1"}')).toBe(true);
  });
  it("allows catalog schema vocabulary only in compiled transport code", () => {
    const schema = '{"properties":["id","accession"],"types":{"accession":"string"}}';
    expect(hasForbiddenPublicIdentity(schema)).toBe(true);
    expect(hasForbiddenPublicIdentity(schema, { compiledBundle: true })).toBe(false);
    expect(hasForbiddenPublicIdentity('organism="private source"', { compiledBundle: true })).toBe(
      true,
    );
  });

  it.each(["CP000001", "AP009389.1", "GCA_000000001.1", "GCF_000000001.2", "NC_000001.1"])(
    "rejects embedded genomic identity %s even in compiled JavaScript",
    (accession) => {
      expect(
        hasForbiddenPublicIdentity(`const record={accession:"${accession}"}`, {
          compiledBundle: true,
        }),
      ).toBe(true);
    },
  );
});
