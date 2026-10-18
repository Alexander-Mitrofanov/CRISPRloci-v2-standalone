<script setup>
import { useAnalysisForm } from "../../features/submission/index.js";
import AppIcon from "../common/AppIcon.vue";
import AdvancedOptions from "./AdvancedOptions.vue";
import FastaInput from "./FastaInput.vue";
import { ANALYSIS_MODES } from "../../science.js";
import AssociationGrouping from "./AssociationGrouping.vue";
import SpacerSearchControls from "./SpacerSearchControls.vue";
import ReadinessPanel from "./ReadinessPanel.vue";

const props = defineProps({
  service: { type: Object, required: true },
  limits: { type: Object, required: true },
  hasActiveJob: Boolean,
  initialInput: { type: Object, default: null },
  initialMode: {
    type: String,
    default: "loci",
    validator: (value) => ANALYSIS_MODES.some((item) => item.id === value),
  },
});
const emit = defineEmits(["submitted", "example-loaded", "back", "change-mode"]);
const {
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
} = useAnalysisForm(props, emit);
</script>

<template>
  <section
    id="workflow"
    class="workflow"
    aria-labelledby="workflow-title"
  >
    <button
      class="back-to-methods"
      type="button"
      @click="$emit('back')"
    >
      All methods
    </button>
    <div class="section-intro">
      <h1
        id="workflow-title"
        tabindex="-1"
      >
        {{
          ["protospacer", "viral_search"].includes(mode) ? "Spacer searches" : selectedMode.title
        }}
      </h1>
      <p
        v-if="!['protospacer', 'viral_search'].includes(mode)"
        class="method-description"
      >
        {{ selectedMode.description }}
      </p>
    </div>
    <p
      v-if="initialInput"
      role="note"
    >
      Loaded {{ initialInput.filename }} from a completed analysis. Review the input and options.
      Save the previous job's recovery link before submitting; the new job has its own recovery
      link.
    </p>
    <SpacerSearchControls
      v-if="['protospacer', 'viral_search'].includes(mode)"
      :mode="mode"
      v-model:mismatches="options.viralMaxMismatches"
      :service="service"
      @direction="emit('change-mode', $event)"
    />
    <form
      id="analysis-form"
      novalidate
      @submit.prevent="submit"
    >
      <section
        class="input-section"
        aria-labelledby="input-step-heading"
      >
        <h3
          id="input-step-heading"
          class="section-title"
        >
          Sequences
        </h3>
        <fieldset
          v-if="['repeats', 'rna_family'].includes(mode)"
          class="repeat-molecule"
        >
          <legend>Input molecule</legend>
          <label
            ><input
              v-model="options.molecule"
              type="radio"
              value="DNA"
              name="repeat-molecule"
            />
            {{
              mode === "rna_family"
                ? "DNA — family members in matching orientation"
                : "DNA repeats — both orientation hypotheses"
            }}</label
          >
          <label
            ><input
              v-model="options.molecule"
              type="radio"
              value="RNA"
              name="repeat-molecule"
            />
            {{
              mode === "rna_family"
                ? "RNA — family members in transcribed 5′→3′ order"
                : "RNA repeats — transcribed 5′→3′ sequence"
            }}</label
          >
          <p v-if="mode === 'repeats'">
            Up to 1,000 repeats, 200 nt per record, within the service limits. Use T for DNA and U
            for RNA; mixed T/U input is rejected.
          </p>
          <p v-else>
            2–50 unambiguous sequences, 30–500 nt each, within service limits. At least two distinct
            sequences are required. Exact duplicates receive one folding weight. Supply comparable
            boundaries; homology and orientation are not inferred.
          </p>
        </fieldset>
        <AssociationGrouping
          v-if="mode === 'spacer_association'"
          v-model="options.associationGrouping"
        />
        <div class="input-layout">
          <FastaInput
            v-model:sequence="sequence"
            v-model:filename="filename"
            :inspection="inspection"
            :array-input="mode === 'array_compare'"
            :repeat-input="['repeats', 'repeat_map', 'repeat_type'].includes(mode)"
            :family-input="mode === 'rna_family'"
            :locus-input="mode === 'locus_compare'"
            :spacer-input="['viral_search', 'spacer_association'].includes(mode)"
            :molecule="molecule"
            :loading-example="loadingExample"
            :show-example="hasExample"
            :example-disabled="hasActiveJob"
            :max-request-bytes="limits.maxRequestBytes"
            @load-example="loadExample"
          /><ReadinessPanel
            :inspection="inspection"
            :selected-mode="selectedMode"
            :limits="inputLimits"
            :molecule="molecule"
            :service="service"
            :request-bytes="requestBytes"
            :precomputed="precomputedPolicyMatches"
          />
        </div>
      </section>
      <section
        v-if="
          ![
            'locus_compare',
            'prophage_context',
            'insertion_context',
            'defense_context',
            'rna_family',
          ].includes(mode)
        "
        class="policy-section"
        aria-labelledby="policy-step-heading"
      >
        <h3
          id="policy-step-heading"
          class="sr-only"
        >
          Options
        </h3>
        <AdvancedOptions
          v-model="options"
          :mode="mode"
        />
      </section>
      <div
        v-if="error"
        class="submit-error"
        role="alert"
      >
        <AppIcon name="warning" /><span>{{ error }}</span>
      </div>
      <div
        v-if="hasActiveJob"
        class="active-job-lock"
        role="status"
      >
        <AppIcon name="info" /><span
          ><strong>Another job is open.</strong> Copy its recovery link, then cancel it or leave it
          after completion.</span
        >
      </div>
      <div
        class="privacy-notice"
        role="note"
        aria-label="Sequence privacy and retention notice"
      >
        <AppIcon name="shield" />
        <p>
          <strong>Non-sensitive research data only.</strong> Sequences are sent to the service
          operator. Results are available as soon as the analysis completes. Download a copy to
          keep.
        </p>
        <details class="privacy-details">
          <summary>Data privacy</summary>
          <p>
            The job token protects retrieval; it is not end-to-end encryption from the operator. Do
            not submit personal, clinical, controlled, or unpublished sensitive sequences; use an
            institutionally approved private route. Results are stored temporarily; each job shows
            its expiry time. The bundled example is matched locally and is never submitted while its
            input and options match the stored result.
          </p>
        </details>
      </div>
      <div class="submit-bar">
        <div>
          <span>{{ selectedMode.tools.join(" → ") }}</span>
        </div>
        <button
          class="primary-button"
          type="submit"
          :disabled="!ready || submitting"
        >
          {{
            submitting
              ? precomputedPolicyMatches
                ? "Loading result…"
                : "Submitting…"
              : hasActiveJob
                ? "Current job still open"
                : precomputedPolicyMatches
                  ? "View precomputed result"
                  : ["protospacer", "viral_search"].includes(mode)
                    ? "Search references"
                    : "Compute"
          }}<AppIcon name="arrow" />
        </button>
      </div>
    </form>
  </section>
</template>
