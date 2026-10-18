import { computed, onScopeDispose, ref, watch } from "vue";

import { api } from "../../api.js";
import { saveBlob } from "../../utils/download.js";
import { asArray, downloadName } from "../../utils/formatting.js";

function artifactId(artifact) {
  return String(artifact?.artifact_id || artifact?.id || "");
}

export function useArtifactDownloads(job, credential, client = api, save = saveBlob) {
  const downloading = ref("");
  const error = ref("");
  let downloadLatch = false;
  let downloadController;

  const artifacts = computed(() =>
    asArray(job.value?.artifacts || job.value?.summary?.artifacts).filter(
      (artifact) =>
        !Number.isFinite(Number(artifact?.size_bytes)) || Number(artifact.size_bytes) > 0,
    ),
  );
  const bundle = computed(() =>
    artifacts.value.find((artifact) =>
      `${artifact?.kind || ""} ${artifact?.filename || artifact?.name || ""}`
        .toLowerCase()
        .match(/bundle|archive|results\.zip/),
    ),
  );
  const individual = computed(() =>
    artifacts.value
      .filter((artifact) => artifact !== bundle.value)
      .sort((a, b) => {
        const priority = [
          "result.json",
          "arrays.json",
          "repeats.json",
          "spacers.json",
          "input-sequences.json",
        ];
        const rank = (artifact) =>
          priority.includes(artifact.filename || artifact.name)
            ? priority.indexOf(artifact.filename || artifact.name)
            : priority.length;
        return rank(a) - rank(b);
      }),
  );

  function abortDownload() {
    downloadController?.abort();
    downloadController = undefined;
    downloadLatch = false;
    downloading.value = "";
  }

  // Completion and expiry refreshes replace these objects without changing the
  // job capability. Only leaving that job should interrupt its download.
  watch(
    [
      () => job.value?.job_id,
      () => Boolean(job.value),
      () => credential.value?.jobId,
      () => credential.value?.accessToken,
    ],
    abortDownload,
    { flush: "sync" },
  );
  onScopeDispose(abortDownload);

  async function download(artifact = null) {
    if (downloadLatch) return;
    downloadLatch = true;
    const id = artifact ? artifactId(artifact) : "bundle";
    downloading.value = id;
    error.value = "";
    const controller = new AbortController();
    downloadController = controller;
    const currentCredential = credential.value;
    const currentBundle = bundle.value;
    try {
      const blob = artifact
        ? await client.downloadArtifact(
            currentCredential.jobId,
            id,
            currentCredential.accessToken,
            { signal: controller.signal },
          )
        : currentBundle
          ? await client.downloadArtifact(
              currentCredential.jobId,
              artifactId(currentBundle),
              currentCredential.accessToken,
              { signal: controller.signal },
            )
          : await client.downloadBundle(currentCredential.jobId, currentCredential.accessToken, {
              signal: controller.signal,
            });
      if (controller.signal.aborted) return;
      save(
        blob,
        artifact
          ? downloadName(artifact.filename || artifact.name, `${id}.dat`)
          : `crispr-analysis-${currentCredential.jobId}.zip`,
      );
    } catch (downloadError) {
      if (!controller.signal.aborted && downloadError.name !== "AbortError") {
        error.value = downloadError.message || "Download failed.";
      }
    } finally {
      if (downloadController === controller) {
        downloadController = undefined;
        downloadLatch = false;
        downloading.value = "";
      }
    }
  }

  return { artifacts, bundle, individual, downloading, error, download };
}
