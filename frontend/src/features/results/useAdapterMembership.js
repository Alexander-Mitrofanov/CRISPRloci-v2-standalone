import { ref, watch } from "vue";

import { api } from "../../api.js";
import { asArray } from "../../utils/formatting.js";
import {
  MAX_ADAPTER_MANIFEST_BYTES,
  adapterHasMembership,
  sanitizeAdapterMembership,
} from "../../utils/results.js";

export function useAdapterMembership(summary, job, credential, client = api) {
  const artifactGroups = ref(null);
  const membershipStatus = ref("idle");

  watch(
    [summary, job, credential],
    async ([currentSummary, currentJob, currentCredential], _old, onCleanup) => {
      if (adapterHasMembership(currentSummary)) {
        artifactGroups.value = null;
        membershipStatus.value = "inline";
        return;
      }
      const artifact = asArray(currentJob?.artifacts || currentSummary?.artifacts).find(
        (item) => String(item?.name || item?.filename || "") === "adapter/manifest.json",
      );
      const artifactId = artifact ? String(artifact.artifact_id || artifact.id || "") : "";
      if (!artifactId || !currentCredential?.jobId || !currentCredential?.accessToken) {
        artifactGroups.value = null;
        membershipStatus.value = "unavailable";
        return;
      }
      const controller = new AbortController();
      onCleanup(() => controller.abort());
      artifactGroups.value = null;
      membershipStatus.value = "loading";
      try {
        const blob = await client.downloadArtifact(
          currentCredential.jobId,
          artifactId,
          currentCredential.accessToken,
          { signal: controller.signal },
        );
        if (controller.signal.aborted) return;
        if (Number(blob?.size) > MAX_ADAPTER_MANIFEST_BYTES) {
          throw new Error("Adapter manifest is too large.");
        }
        const text = await blob.text();
        if (controller.signal.aborted) return;
        if (new TextEncoder().encode(text).byteLength > MAX_ADAPTER_MANIFEST_BYTES) {
          throw new Error("Adapter manifest is too large.");
        }
        const groups = sanitizeAdapterMembership(JSON.parse(text));
        if (controller.signal.aborted) return;
        if (!groups.some((group) => group.arrays.length)) {
          throw new Error("Adapter manifest has no valid membership.");
        }
        artifactGroups.value = groups;
        membershipStatus.value = "loaded";
      } catch (error) {
        if (!controller.signal.aborted && error.name !== "AbortError") {
          membershipStatus.value = "unavailable";
        }
      }
    },
    { immediate: true },
  );

  return { artifactGroups, membershipStatus };
}
