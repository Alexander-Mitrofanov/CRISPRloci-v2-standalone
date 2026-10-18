import { computed, onScopeDispose, ref, watch } from "vue";
import { api } from "../../api.js";
import ports from "../../contracts/module-ports-v1.json";
import { ANALYSIS_MODES, analysisModeAvailable } from "../../science.js";
import { inspectSequenceInput } from "../submission/index.js";

const SUCCESS = new Set(["completed", "completed_no_eligible_groups"]);

export function useSequenceHandoff(props, emit, client = api) {
  const selected = ref("");
  const loading = ref(false);
  const error = ref("");
  let controller;
  const choices = computed(() => {
    if (!props.credential || !SUCCESS.has(props.job?.status)) return [];
    return (props.job.artifacts || []).flatMap((artifact) => {
      const kind = {
        "arrays.json": "array",
        "repeats.json": "repeat",
        "spacers.json": "spacer",
        "input-sequences.json":
          props.job.mode === "array_compare" ? "spacer" : ports[props.job.mode],
      }[artifact.name];
      if (!kind || artifact.size_bytes > 20_000_000) return [];
      return ANALYSIS_MODES.filter(
        (mode) => ports[mode.id] === kind && analysisModeAvailable(mode.id, props.service),
      ).map((mode) => ({
        id: `${artifact.artifact_id}:${mode.id}`,
        artifact,
        mode: mode.id,
        label: `${artifact.name} → ${mode.title}`,
      }));
    });
  });

  function abort() {
    controller?.abort();
    controller = undefined;
    loading.value = false;
    error.value = "";
    selected.value = "";
  }
  watch(
    [() => props.job?.job_id, () => props.credential?.jobId, () => props.credential?.accessToken],
    abort,
    { flush: "sync" },
  );
  onScopeDispose(abort);

  async function prepare() {
    const choice = choices.value.find((item) => item.id === selected.value);
    if (!choice || loading.value) return;
    const request = new AbortController();
    controller = request;
    loading.value = true;
    error.value = "";
    try {
      const text = await client.readSequenceArtifact(
        props.credential.jobId,
        choice.artifact.artifact_id,
        props.credential.accessToken,
        { signal: request.signal, maxBytes: props.limits.maxRequestBytes || 20_000_000 },
      );
      if (request.signal.aborted) return;
      const packet = JSON.parse(text.replace(/^\uFEFF/, ""));
      const molecule = packet.molecule || "DNA";
      if (molecule === "RNA" && !["repeats", "rna_family"].includes(choice.mode))
        throw new Error("This analysis does not accept RNA input.");
      const inspection = inspectSequenceInput(text, { mode: choice.mode, molecule });
      if (!inspection.valid) throw new Error(inspection.errors.join(" "));
      emit("continue-analysis", {
        mode: choice.mode,
        sequence: text,
        filename: choice.artifact.name,
        molecule,
        sourceJobId: props.job.job_id,
      });
    } catch (failure) {
      if (!request.signal.aborted) error.value = failure.message || "Could not load this input.";
    } finally {
      if (controller === request) {
        controller = undefined;
        loading.value = false;
      }
    }
  }
  return { choices, selected, loading, error, prepare };
}
