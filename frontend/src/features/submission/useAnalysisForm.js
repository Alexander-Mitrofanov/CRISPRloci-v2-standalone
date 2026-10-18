import { validAssociationGrouping } from "./association.js";
import { computed, ref } from "vue";

import { ApiError, api } from "../../api.js";
import {
  exampleAssetsFor,
  exampleSubmissionMatches,
  readExampleResponse,
  validateExampleInput,
} from "../../example.js";
import { inspectSequenceInput } from "./sequenceInput.js";
import { normalizeJobCredential } from "../../jobStore.js";
import { ANALYSIS_MODES, DEFAULT_CATEGORY_POLICY, analysisModeAvailable } from "../../science.js";
import { buildSubmission } from "../../submission.js";
import { toolOptionCount, validToolOptions } from "../../toolOptions.js";

const INITIAL_OPTIONS = Object.freeze({
  categoryPolicy: DEFAULT_CATEGORY_POLICY,
  spacerEditDistance: 1,
  biasCorrection: true,
  tracrModelType: "II",
  leaderFlankLength: 500,
  molecule: "DNA",
  associationGrouping: "",
  viralMaxMismatches: 2,
  toolOptions: {},
});

export function useAnalysisForm(props, emit, client = api, fetcher = globalThis.fetch) {
  const mode = ref(props.initialMode || "loci");
  const sequence = ref(props.initialInput?.sequence || "");
  const filename = ref(props.initialInput?.filename || "input.fasta");
  const options = ref({ ...INITIAL_OPTIONS, molecule: props.initialInput?.molecule || "DNA" });
  const submitting = ref(false);
  const loadingExample = ref(false);
  const error = ref("");
  const preparedExample = ref(null);
  const preparedExampleSequence = ref("");
  let submittingLatch = false;
  let exampleLatch = false;

  const molecule = computed(() =>
    ["repeats", "rna_family"].includes(mode.value) ? options.value.molecule : "DNA",
  );
  const hasExample = computed(() => Boolean(exampleAssetsFor(mode.value, options.value)));
  const inputLimits = computed(() =>
    ["locus_compare", "rna_family"].includes(mode.value)
      ? {
          ...props.limits,
          maxRecords: Math.min(
            props.limits.maxRecords || 50,
            mode.value === "locus_compare" ? 8 : 50,
          ),
          maxRecordBases: Math.min(
            props.limits.maxRecordBases || 50000,
            mode.value === "locus_compare" ? 50000 : 500,
          ),
        }
      : ["prophage_context", "insertion_context", "defense_context"].includes(mode.value)
        ? {
            ...props.limits,
            maxRecords: Math.min(props.limits.maxRecords || 20, 20),
            maxBases: Math.min(props.limits.maxBases || 2_000_000, 2_000_000),
            maxRecordBases: Math.min(props.limits.maxRecordBases || 2_000_000, 2_000_000),
          }
        : mode.value === "viral_search"
          ? { ...props.limits, maxRecordBases: Math.min(props.limits.maxRecordBases || 80, 80) }
          : ["repeats", "repeat_map", "repeat_type", "spacer_association"].includes(mode.value)
            ? {
                ...props.limits,
                maxRecords: Math.min(
                  props.limits.maxRecords || 1000,
                  mode.value === "repeat_map" ? 100 : 1000,
                ),
                maxRecordBases: Math.min(props.limits.maxRecordBases || 200, 200),
              }
            : mode.value === "protospacer"
              ? {
                  ...props.limits,
                  maxBases: Math.min(props.limits.maxBases || 5_000_000, 5_000_000),
                  maxRecordBases: Math.min(props.limits.maxRecordBases || 5_000_000, 5_000_000),
                }
              : mode.value === "repeat_context"
                ? { ...props.limits, maxRecords: Math.min(props.limits.maxRecords || 1000, 1000) }
                : props.limits,
  );

  const inspection = computed(() =>
    inspectSequenceInput(sequence.value, {
      mode: mode.value,
      maxHeaderCharacters: props.limits.maxHeaderCharacters || 200,
      molecule: molecule.value,
    }),
  );
  const selectedMode = computed(() => ANALYSIS_MODES.find((item) => item.id === mode.value));
  const submission = computed(() =>
    buildSubmission({
      sequence: sequence.value,
      filename: filename.value,
      mode: mode.value,
      options: options.value,
    }),
  );
  const requestBytes = computed(
    () => new TextEncoder().encode(JSON.stringify(submission.value)).byteLength,
  );
  const precomputedPolicyMatches = computed(() => {
    const recorded = preparedExample.value?.job?.options;
    return Boolean(
      preparedExample.value &&
      sequence.value === preparedExampleSequence.value &&
      mode.value === preparedExample.value.job?.mode &&
      (preparedExample.value.schema.name === "crisprloci.web-example"
        ? exampleSubmissionMatches(preparedExample.value, submission.value)
        : toolOptionCount(mode.value, options.value.toolOptions) === 0 &&
          submission.value.category_policy === recorded?.category_policy &&
          submission.value.spacer_distance === recorded?.spacer_distance &&
          submission.value.bias_corrections === recorded?.bias_corrections_requested),
    );
  });
  const withinLimits = computed(
    () =>
      (!inputLimits.value.maxRecords ||
        inspection.value.recordCount <= inputLimits.value.maxRecords) &&
      (!inputLimits.value.maxBases || inspection.value.baseCount <= inputLimits.value.maxBases) &&
      (!inputLimits.value.maxRecordBases ||
        inspection.value.records.every(
          (record) => record.sequence.length <= inputLimits.value.maxRecordBases,
        )),
  );
  const withinRequest = computed(
    () => !props.limits.maxRequestBytes || requestBytes.value <= props.limits.maxRequestBytes,
  );
  const ready = computed(
    () =>
      inspection.value.valid &&
      inspection.value.recordCount >= selectedMode.value.minimumRecords &&
      (mode.value !== "rna_family" ||
        (inspection.value.records.every(
          (record) => record.sequence.length >= 30 && /^[ACGTU]+$/.test(record.sequence),
        ) &&
          new Set(inspection.value.records.map((record) => record.sequence.replaceAll("T", "U")))
            .size >= 2)) &&
      (precomputedPolicyMatches.value || withinLimits.value) &&
      (precomputedPolicyMatches.value || withinRequest.value) &&
      validToolOptions(mode.value, options.value.toolOptions) &&
      (mode.value !== "spacer_association" ||
        (validAssociationGrouping(sequence.value, options.value.associationGrouping) &&
          (options.value.associationGrouping !== "per_record" ||
            inspection.value.recordCount <= 100))) &&
      (!["loci", "leader"].includes(mode.value) ||
        (Number.isInteger(options.value.leaderFlankLength) &&
          options.value.leaderFlankLength >= 1 &&
          options.value.leaderFlankLength <= 5000)) &&
      (props.service.state === "online" || precomputedPolicyMatches.value) &&
      (precomputedPolicyMatches.value || analysisModeAvailable(mode.value, props.service)) &&
      !props.hasActiveJob,
  );

  async function loadExample() {
    const assets = exampleAssetsFor(mode.value, options.value);
    if (exampleLatch || !assets || props.hasActiveJob) return;
    exampleLatch = true;
    loadingExample.value = true;
    error.value = "";
    try {
      const [inputResponse, resultResponse] = await Promise.all([
        fetcher(`${import.meta.env.BASE_URL}${assets.input}`, {
          cache: "no-store",
          credentials: "omit",
          referrerPolicy: "no-referrer",
        }),
        fetcher(`${import.meta.env.BASE_URL}${assets.result}`, {
          cache: "no-store",
          credentials: "omit",
          referrerPolicy: "no-referrer",
        }),
      ]);
      if (!inputResponse.ok || !resultResponse.ok) {
        throw new Error("The example demonstration could not be loaded.");
      }
      const [exampleSequence, snapshotText] = await Promise.all([
        readExampleResponse(inputResponse),
        readExampleResponse(resultResponse),
      ]);
      const rawSnapshot = JSON.parse(snapshotText);
      const { snapshot } = await validateExampleInput(rawSnapshot, exampleSequence, {
        maxHeaderCharacters: props.limits.maxHeaderCharacters || 200,
      });
      if (snapshot.job.mode !== mode.value)
        throw new Error("The example result belongs to a different method.");
      const recorded = snapshot.example.submission_options || {};
      options.value = {
        ...INITIAL_OPTIONS,
        categoryPolicy: recorded.category_policy || INITIAL_OPTIONS.categoryPolicy,
        spacerEditDistance: recorded.spacer_distance ?? INITIAL_OPTIONS.spacerEditDistance,
        biasCorrection: recorded.bias_corrections ?? INITIAL_OPTIONS.biasCorrection,
        molecule: recorded.molecule || "DNA",
        tracrModelType: recorded.tracr_model_type || "II",
        leaderFlankLength: recorded.leader_flank_length ?? 500,
        viralMaxMismatches: recorded.viral_max_mismatches ?? 2,
        associationGrouping: recorded.association_grouping || "",
        toolOptions: recorded.tool_options || {},
      };
      sequence.value = exampleSequence;
      filename.value = snapshot.example.input.filename;
      preparedExample.value = snapshot;
      preparedExampleSequence.value = exampleSequence;
      emit("example-loaded", null);
    } catch (loadError) {
      error.value = loadError.message || "The example demonstration could not be loaded.";
    } finally {
      exampleLatch = false;
      loadingExample.value = false;
    }
  }

  async function submit() {
    if (!ready.value || submittingLatch) return;
    submittingLatch = true;
    submitting.value = true;
    error.value = "";
    try {
      if (precomputedPolicyMatches.value) {
        const snapshot = preparedExample.value;
        const exampleSequence = sequence.value;
        await validateExampleInput(snapshot, exampleSequence, {
          maxHeaderCharacters: props.limits.maxHeaderCharacters || 200,
        });
        if (
          !precomputedPolicyMatches.value ||
          preparedExample.value !== snapshot ||
          sequence.value !== exampleSequence ||
          props.hasActiveJob
        )
          throw new Error("The example input or options changed. Review them before continuing.");
        emit("example-loaded", snapshot);
        return;
      }
      const response = await client.submit(submission.value);
      const submittedJob = response?.job || response;
      const jobId = response?.job_id || submittedJob?.job_id || submittedJob?.id;
      if (!jobId || !response?.access_token) {
        throw new ApiError("The service returned an incomplete job credential.");
      }
      const credential = normalizeJobCredential({
        jobId,
        accessToken: response.access_token,
        expiresAt: response?.expires_at || submittedJob?.expires_at,
      });
      emit("submitted", credential, { ...submittedJob, mode: mode.value });
    } catch (submitError) {
      error.value = submitError.message || "The analysis could not be submitted.";
    } finally {
      submittingLatch = false;
      submitting.value = false;
    }
  }

  return {
    mode,
    molecule,
    inputLimits,
    sequence,
    filename,
    options,
    submitting,
    loadingExample,
    hasExample,
    error,
    inspection,
    selectedMode,
    requestBytes,
    precomputedPolicyMatches,
    ready,
    loadExample,
    submit,
  };
}
