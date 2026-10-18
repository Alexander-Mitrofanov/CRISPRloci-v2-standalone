import { describe, expect, it } from "vitest";

import {
  buildRecoveryUrl,
  isRecoveryHash,
  normalizeJobCredential,
  parseRecoveryHash,
  serializeRecoveryHash,
} from "../src/jobStore.js";

const credential = {
  jobId: "0123456789abcdef0123456789abcdef",
  accessToken: "a".repeat(43),
  expiresAt: "2099-01-01T00:00:00Z",
};

describe("job recovery links", () => {
  it("round-trips a strict fragment capability without browser storage", () => {
    const hash = serializeRecoveryHash(credential);
    expect(hash).toBe(`#job=${credential.jobId}.${credential.accessToken}`);
    expect(parseRecoveryHash(hash)).toEqual({ ...credential, expiresAt: null });
    expect(isRecoveryHash(hash)).toBe(true);
    expect(window.localStorage.length).toBe(0);
    expect(window.sessionStorage.length).toBe(0);
  });

  it("builds a same-page recovery URL without placing the capability in its query", () => {
    const value = buildRecoveryUrl(
      credential,
      "https://example.test/crispr-evolution-web/?language=en#workflow",
    );
    const url = new URL(value);
    expect(url.pathname).toBe("/crispr-evolution-web/");
    expect(url.search).toBe("?language=en");
    expect(url.search).not.toContain(credential.accessToken);
    expect(url.hash).toBe(serializeRecoveryHash(credential));
  });

  it("keeps queued credentials valid before terminal retention is known", () => {
    expect(normalizeJobCredential({ ...credential, expiresAt: null })).toEqual({
      ...credential,
      expiresAt: null,
    });
  });

  it("rejects expired credentials and malformed or extended recovery hashes", () => {
    expect(() =>
      normalizeJobCredential(
        { ...credential, expiresAt: "2020-01-01T00:00:00Z" },
        Date.parse("2021-01-01"),
      ),
    ).toThrow(/expired/i);
    expect(() => parseRecoveryHash("#workflow")).toThrow(/supported job hash/i);
    expect(() => parseRecoveryHash("#job=bad")).toThrow(/invalid job hash/i);
    expect(() => parseRecoveryHash(`${serializeRecoveryHash(credential)}&extra=1`)).toThrow(
      /invalid job hash/i,
    );
    expect(() => normalizeJobCredential({ ...credential, accessToken: "short" })).toThrow(
      /invalid job ID or access token/i,
    );
    expect(() => normalizeJobCredential({ ...credential, unexpected: true })).toThrow(/exactly/i);
  });
});
