import { onScopeDispose, ref, shallowRef, watch } from "vue";
import { api } from "../../api.js";
import { genomeArtifacts, requireValue } from "./contract.js";
import { buildGenomeModel } from "./model.js";

export function useGenomeData(props, client = api) {
  const model = shallowRef(null);
  const loading = ref(false);
  const error = ref("");
  let controller;

  function reset() {
    controller?.abort();
    controller = undefined;
    model.value = null;
    error.value = "";
    loading.value = false;
  }

  async function load() {
    if (loading.value || model.value || !props.active) return;
    const request = new AbortController();
    controller = request;
    loading.value = true;
    error.value = "";
    try {
      requireValue(props.jobId && props.accessToken, "This job's recovery capability is required.");
      const artifacts = genomeArtifacts(props.artifacts);
      const documents = [];
      for (const artifact of [artifacts.result, artifacts.input]) {
        requireValue(
          artifact?.artifact_id && artifact.size_bytes <= 20_000_000,
          "The complete genome artifacts are unavailable within the viewing limits.",
        );
        const text = await client.readSequenceArtifact(
          props.jobId,
          artifact.artifact_id,
          props.accessToken,
          { signal: request.signal, maxBytes: 20_000_000 },
        );
        if (request.signal.aborted) return;
        documents.push(JSON.parse(text));
      }
      const value = await buildGenomeModel(documents[0], documents[1], props.mode);
      if (!request.signal.aborted) model.value = value;
    } catch (failure) {
      if (!request.signal.aborted)
        error.value = failure.message || "The genome evidence could not be loaded.";
    } finally {
      if (controller === request) {
        loading.value = false;
        controller = undefined;
      }
    }
  }

  watch(
    () => [props.jobId, props.accessToken, props.mode],
    () => {
      reset();
      void load();
    },
    { immediate: true, flush: "sync" },
  );
  watch(
    () => props.active,
    (active) => {
      if (active && !error.value) void load();
    },
  );
  onScopeDispose(reset);
  return { model, loading, error, load };
}
