import { describe, expect, it, vi } from "vitest";

import { ApiError, createApiClient, normalizeApiBase } from "../src/api.js";

function jsonResponse(payload, status = 200) {
  return new Response(JSON.stringify(payload), {
    status,
    headers: { "Content-Type": "application/json" },
  });
}

describe("API client", () => {
  it("reads sequence artifacts with bounded bytes and header-only authentication", async () => {
    const fetchMock = vi.fn(async () => new Response('{"records":[]}'));
    const client = createApiClient("https://analysis.example.org", fetchMock);
    expect(
      await client.readSequenceArtifact("job", "port", "private-token", { maxBytes: 100 }),
    ).toBe('{"records":[]}');
    const [url, options] = fetchMock.mock.calls[0];
    expect(url).not.toContain("private-token");
    expect(options).toMatchObject({
      cache: "no-store",
      credentials: "omit",
      referrerPolicy: "no-referrer",
      headers: { Authorization: "Bearer private-token" },
    });
    await expect(
      client.readSequenceArtifact("job", "port", "private-token", { maxBytes: 3 }),
    ).rejects.toThrow(/byte limit/);
  });

  it("rejects a declared oversized sequence file before reading it", async () => {
    const cancel = vi.fn();
    const fetchMock = vi.fn(async () => ({
      ok: true,
      headers: new Headers({ "content-length": "20000001" }),
      body: { cancel },
    }));
    const client = createApiClient("https://analysis.example.org", fetchMock);
    await expect(client.readSequenceArtifact("job", "port", "token")).rejects.toThrow(/byte limit/);
    expect(cancel).toHaveBeenCalledOnce();
  });
  it("requires HTTPS except on a local development host", () => {
    expect(normalizeApiBase("https://analysis.example.org/")).toBe("https://analysis.example.org");
    expect(normalizeApiBase("http://localhost:8080")).toBe("http://localhost:8080");
    expect(normalizeApiBase("http://[::1]:8080")).toBe("http://[::1]:8080");
    expect(() => normalizeApiBase("http://analysis.example.org")).toThrow(ApiError);
    expect(() => normalizeApiBase("https://user@example.org")).toThrow(ApiError);
    expect(() => normalizeApiBase("https://analysis.example.org/api")).toThrow(ApiError);
    expect(() => normalizeApiBase("https://analysis.example.org?x=1")).toThrow(ApiError);
    expect(() => normalizeApiBase("https://analysis.example.org:bad")).toThrow(ApiError);
  });

  it("submits JSON without placing a job token in the request", async () => {
    const fetchMock = vi
      .fn()
      .mockResolvedValue(jsonResponse({ job_id: "abc", access_token: "secret" }, 202));
    const client = createApiClient("https://analysis.example.org", fetchMock);
    const payload = { sequence: ">a\nACGT\n", mode: "detection" };

    await client.submit(payload);

    const [url, request] = fetchMock.mock.calls[0];
    expect(url).toBe("https://analysis.example.org/api/v1/jobs");
    expect(request.method).toBe("POST");
    expect(request.headers.Authorization).toBeUndefined();
    expect(JSON.parse(request.body)).toEqual(payload);
  });

  it("authenticates job polling with a Bearer header and never a URL token", async () => {
    const fetchMock = vi.fn().mockResolvedValue(jsonResponse({ id: "job-1", status: "queued" }));
    const client = createApiClient("https://analysis.example.org", fetchMock);

    await client.getJob("job/one", "private-token");

    const [url, request] = fetchMock.mock.calls[0];
    expect(url).toBe("https://analysis.example.org/api/v1/jobs/job%2Fone");
    expect(url).not.toContain("private-token");
    expect(request.headers.Authorization).toBe("Bearer private-token");
  });

  it("authenticates cancellation with the same credential", async () => {
    const fetchMock = vi.fn().mockResolvedValue(jsonResponse({ id: "job-1", status: "cancelled" }));
    const client = createApiClient("https://analysis.example.org", fetchMock);

    await client.cancelJob("job-1", "private-token");

    const [, request] = fetchMock.mock.calls[0];
    expect(request.method).toBe("DELETE");
    expect(request.headers.Authorization).toBe("Bearer private-token");
  });

  it("downloads artifacts as authenticated blobs without token-bearing links", async () => {
    const blob = new Blob(["tree"], { type: "text/plain" });
    const fetchMock = vi.fn().mockResolvedValue(new Response(blob));
    const client = createApiClient("https://analysis.example.org", fetchMock);

    const result = await client.downloadArtifact("job-1", "tree.nwk", "private-token");

    const [url, request] = fetchMock.mock.calls[0];
    expect(url).toBe("https://analysis.example.org/api/v1/jobs/job-1/artifacts/tree.nwk");
    expect(url).not.toContain("private-token");
    expect(request.headers.Authorization).toBe("Bearer private-token");
    expect(result).toBeInstanceOf(Blob);
  });

  it("omits ambient cookies and cache state from every service request", async () => {
    const fetchMock = vi.fn(async () => jsonResponse({ ok: true }));
    const client = createApiClient("https://analysis.example.org", fetchMock);

    await client.health();
    await client.config();
    await client.submit({ sequence: ">a\nACGT\n", mode: "detection" });
    await client.getJob("a".repeat(32), "private-token");
    await client.cancelJob("a".repeat(32), "private-token");
    await client.downloadBundle("a".repeat(32), "private-token");
    await client.downloadArtifact("a".repeat(32), "b".repeat(24), "private-token");

    expect(fetchMock).toHaveBeenCalledTimes(7);
    for (const [, request] of fetchMock.mock.calls) {
      expect(request).toMatchObject({
        cache: "no-store",
        credentials: "omit",
        referrerPolicy: "no-referrer",
      });
    }
  });
});
