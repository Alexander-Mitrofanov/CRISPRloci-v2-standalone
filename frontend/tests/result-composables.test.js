import { effectScope, ref } from "vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import { useAdapterMembership, useArtifactDownloads } from "../src/features/results/index.js";

const credential = {
  jobId: "0123456789abcdef0123456789abcdef",
  accessToken: "a".repeat(43),
};

const scopes = [];

function inScope(factory) {
  const scope = effectScope();
  scopes.push(scope);
  return scope.run(factory);
}

afterEach(() => {
  scopes.splice(0).forEach((scope) => scope.stop());
});

describe("result request composables", () => {
  it("puts the compact result before sequence ports and downloads its exact contents", async () => {
    const blob = new Blob([JSON.stringify({ complete: true, data: { matches: [] } })]);
    const client = { downloadArtifact: vi.fn().mockResolvedValue(blob) };
    const save = vi.fn();
    const artifacts = [
      { artifact_id: "native", name: "annotations.json" },
      { artifact_id: "spacers", name: "spacers.json" },
      { artifact_id: "compact", filename: "result.json" },
      { artifact_id: "repeats", name: "repeats.json" },
      { artifact_id: "arrays", name: "arrays.json" },
    ];
    const result = inScope(() =>
      useArtifactDownloads(ref({ artifacts }), ref(credential), client, save),
    );
    expect(result.individual.value.map((artifact) => artifact.artifact_id)).toEqual([
      "compact",
      "arrays",
      "repeats",
      "spacers",
      "native",
    ]);
    await result.download(result.individual.value[0]);
    expect(client.downloadArtifact).toHaveBeenCalledWith(
      credential.jobId,
      "compact",
      credential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
    expect(save).toHaveBeenCalledWith(blob, "result.json");
    expect(artifacts[0].artifact_id).toBe("native");
  });

  it("selects a bundle, filters empty artifacts, and saves authenticated downloads", async () => {
    const blob = new Blob(["zip"]);
    const client = { downloadArtifact: vi.fn().mockResolvedValue(blob) };
    const save = vi.fn();
    const job = ref({
      artifacts: [
        { artifact_id: "empty", filename: "empty.txt", size_bytes: 0 },
        { artifact_id: "bundle-1", filename: "results.zip", size_bytes: 3 },
        { artifact_id: "tree-1", filename: "tree / unsafe.nwk", size_bytes: 4 },
      ],
    });
    const result = inScope(() => useArtifactDownloads(job, ref(credential), client, save));

    expect(result.bundle.value.artifact_id).toBe("bundle-1");
    expect(result.individual.value.map((artifact) => artifact.artifact_id)).toEqual(["tree-1"]);
    await result.download(result.individual.value[0]);

    expect(client.downloadArtifact).toHaveBeenCalledWith(
      credential.jobId,
      "tree-1",
      credential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
    expect(save).toHaveBeenCalledWith(blob, "_unsafe.nwk");
    expect(result.downloading.value).toBe("");
    expect(result.error.value).toBe("");
  });

  it("uses the synthesized bundle endpoint when no archive artifact exists", async () => {
    const blob = new Blob(["zip"]);
    const client = { downloadBundle: vi.fn().mockResolvedValue(blob) };
    const save = vi.fn();
    const result = inScope(() =>
      useArtifactDownloads(ref({ artifacts: [] }), ref(credential), client, save),
    );

    await result.download();

    expect(client.downloadBundle).toHaveBeenCalledWith(
      credential.jobId,
      credential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
    expect(save).toHaveBeenCalledWith(blob, `crispr-analysis-${credential.jobId}.zip`);
  });

  it("loads and sanitizes legacy adapter membership behind one request boundary", async () => {
    const group = {
      name: "repeat_group_1",
      arrays: [{ source_id: "record", array_id: "array-1", category: "Bona-fide" }],
    };
    const client = {
      downloadArtifact: vi
        .fn()
        .mockResolvedValue(
          new Blob([JSON.stringify({ groups: [group] })], { type: "application/json" }),
        ),
    };
    const result = inScope(() =>
      useAdapterMembership(
        ref({ adapter: { groups: [] } }),
        ref({ artifacts: [{ artifact_id: "manifest-1", name: "adapter/manifest.json" }] }),
        ref(credential),
        client,
      ),
    );

    await vi.waitFor(() => expect(result.membershipStatus.value).toBe("loaded"));

    expect(result.artifactGroups.value[0]).toMatchObject({ name: "repeat_group_1" });
    expect(client.downloadArtifact).toHaveBeenCalledWith(
      credential.jobId,
      "manifest-1",
      credential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
  });

  it("does not commit a manifest from an obsolete watcher generation", async () => {
    let resolveOld;
    let resolveNew;
    const requests = {
      old: new Promise((resolve) => {
        resolveOld = resolve;
      }),
      new: new Promise((resolve) => {
        resolveNew = resolve;
      }),
    };
    const client = {
      downloadArtifact: vi.fn((_jobId, artifact) => requests[artifact]),
    };
    const summary = ref({ adapter: { groups: [] } });
    const job = ref({ artifacts: [{ artifact_id: "old", name: "adapter/manifest.json" }] });
    const result = inScope(() => useAdapterMembership(summary, job, ref(credential), client));
    await vi.waitFor(() => expect(client.downloadArtifact).toHaveBeenCalledTimes(1));

    job.value = { artifacts: [{ artifact_id: "new", name: "adapter/manifest.json" }] };
    await vi.waitFor(() => expect(client.downloadArtifact).toHaveBeenCalledTimes(2));
    resolveNew(
      new Blob([
        JSON.stringify({
          groups: [{ name: "new-group", arrays: [{ source_id: "new", array_id: "new-1" }] }],
        }),
      ]),
    );
    await vi.waitFor(() => expect(result.artifactGroups.value?.[0]?.name).toBe("new-group"));

    resolveOld(
      new Blob([
        JSON.stringify({
          groups: [{ name: "old-group", arrays: [{ source_id: "old", array_id: "old-1" }] }],
        }),
      ]),
    );
    await new Promise((resolve) => window.setTimeout(resolve, 0));

    expect(result.artifactGroups.value[0].name).toBe("new-group");
  });

  it("finishes a download while the same job and capability refresh their metadata", async () => {
    let resolveDownload;
    const client = {
      downloadArtifact: vi.fn(
        () =>
          new Promise((resolve) => {
            resolveDownload = resolve;
          }),
      ),
    };
    const save = vi.fn();
    const job = ref({
      job_id: credential.jobId,
      artifacts: [{ artifact_id: "compact", filename: "result.json" }],
    });
    const currentCredential = ref(credential);
    const result = inScope(() => useArtifactDownloads(job, currentCredential, client, save));
    const pending = result.download(result.individual.value[0]);
    const requestOptions = client.downloadArtifact.mock.calls[0][3];

    currentCredential.value = { ...credential, expiresAt: "2099-01-01T00:00:00Z" };
    job.value = { ...job.value, status: "completed", artifacts: [...job.value.artifacts] };
    expect(requestOptions.signal.aborted).toBe(false);
    expect(result.downloading.value).toBe("compact");
    const blob = new Blob(["{}"]);
    resolveDownload(blob);
    await pending;
    expect(save).toHaveBeenCalledWith(blob, "result.json");
  });

  it.each(["jobId", "accessToken", "job", "clear"])(
    "does not save a download after %s changes",
    async (change) => {
      let resolveDownload;
      const client = {
        downloadArtifact: vi.fn(
          () =>
            new Promise((resolve) => {
              resolveDownload = resolve;
            }),
        ),
      };
      const save = vi.fn();
      const job = ref({
        job_id: credential.jobId,
        artifacts: [{ artifact_id: "tree-1", filename: "tree.nwk" }],
      });
      const currentCredential = ref(credential);
      const result = inScope(() => useArtifactDownloads(job, currentCredential, client, save));

      const pending = result.download(result.individual.value[0]);
      await vi.waitFor(() => expect(client.downloadArtifact).toHaveBeenCalledTimes(1));
      const requestOptions = client.downloadArtifact.mock.calls[0][3];
      if (change === "jobId") {
        currentCredential.value = { ...credential, jobId: "fedcba9876543210fedcba9876543210" };
      } else if (change === "accessToken") {
        currentCredential.value = { ...credential, accessToken: "b".repeat(43) };
      } else if (change === "job") {
        job.value = { ...job.value, job_id: "fedcba9876543210fedcba9876543210" };
      } else {
        job.value = null;
      }
      expect(requestOptions.signal.aborted).toBe(true);
      resolveDownload(new Blob(["tree"]));
      await pending;

      expect(save).not.toHaveBeenCalled();
      expect(result.error.value).toBe("");
      expect(result.downloading.value).toBe("");
    },
  );
});
