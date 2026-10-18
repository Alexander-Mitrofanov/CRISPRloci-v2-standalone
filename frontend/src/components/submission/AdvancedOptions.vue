<script setup>
import { CATEGORY_POLICIES } from "../../science.js";
import { computed } from "vue";
import { toolOptionCount } from "../../toolOptions.js";
import AppIcon from "../common/AppIcon.vue";
import ToolOptions from "./ToolOptions.vue";

const props = defineProps({
  mode: { type: String, required: true },
  modelValue: { type: Object, required: true },
});
const emit = defineEmits(["update:modelValue"]);
const update = (values) => emit("update:modelValue", { ...props.modelValue, ...values });
const inclusivePolicy = CATEGORY_POLICIES[0];
const strictPolicy = CATEGORY_POLICIES[1];
const changedCount = computed(() => toolOptionCount(props.mode, props.modelValue.toolOptions));
</script>

<template>
  <details class="advanced">
    <summary>
      <span class="advanced-summary-copy"
        ><strong>Analysis options</strong
        ><small
          v-if="changedCount"
          class="advanced-changed-count"
          >{{ changedCount }} tool {{ changedCount === 1 ? "flag" : "flags" }} changed</small
        ></span
      ><span
        class="advanced-summary-action"
        aria-hidden="true"
        ><span class="advanced-action-closed">Show options</span
        ><span class="advanced-action-open">Hide options</span></span
      >
    </summary>
    <ToolOptions
      :mode="mode"
      :model-value="modelValue.toolOptions || {}"
      @update:model-value="update({ toolOptions: $event })"
    />
    <div
      v-if="
        ![
          'spacer_association',
          'array_compare',
          'repeat_type',
          'repeat_map',
          'repeats',
          'protospacer',
          'viral_search',
          'cas',
        ].includes(mode)
      "
      class="advanced-body"
    >
      <h4 class="workflow-options-title">Workflow settings</h4>
      <div
        v-if="['leader', 'loci'].includes(mode)"
        class="option-column"
      >
        <label for="leader-flank-length"
          ><strong>Leader context window (nt per side)</strong><code>--flank-length</code></label
        >
        <input
          id="leader-flank-length"
          type="number"
          min="1"
          max="5000"
          step="1"
          required
          :value="modelValue.leaderFlankLength"
          @input="update({ leaderFlankLength: Number($event.target.value) })"
        />
        <small
          >1–5,000 nt; default 500. Both sides of Bona-fide and Possible arrays are retained. Leader
          prediction is unavailable; this window is not an inferred leader boundary.</small
        >
      </div>
      <fieldset v-if="['loci', 'tracrrna'].includes(mode)">
        <legend>tracrRNA covariance models <code>--model-type</code></legend>
        <label class="radio-line"
          ><input
            type="radio"
            name="tracr-model"
            value="II"
            :checked="modelValue.tracrModelType === 'II'"
            @change="update({ tracrModelType: 'II' })"
          /><span
            ><strong>Type II</strong><small>Type II-associated tracrRNA models</small></span
          ></label
        >
        <label class="radio-line"
          ><input
            type="radio"
            name="tracr-model"
            value="V"
            :checked="modelValue.tracrModelType === 'V'"
            @change="update({ tracrModelType: 'V' })"
          /><span
            ><strong>Type V-K</strong
            ><small>Four retained Type V-K model families; not all Cas12 subtypes</small></span
          ></label
        >
      </fieldset>
      <fieldset
        v-if="['detection', 'reconstruction', 'orientation', 'repeat_context'].includes(mode)"
      >
        <legend>Eligible CRISPRidentify categories</legend>
        <label class="radio-line"
          ><input
            type="radio"
            name="category-policy"
            :value="inclusivePolicy"
            :checked="modelValue.categoryPolicy === inclusivePolicy"
            @change="update({ categoryPolicy: inclusivePolicy })"
          /><span
            ><strong>Bona-fide + Possible</strong
            ><small>Default · preserves plausible arrays for group preflight</small></span
          ></label
        ><label class="radio-line"
          ><input
            type="radio"
            name="category-policy"
            :value="strictPolicy"
            :checked="modelValue.categoryPolicy === strictPolicy"
            @change="update({ categoryPolicy: strictPolicy })"
          /><span
            ><strong>Strict Bona-fide only</strong
            ><small>Higher specificity, but may leave fewer comparable arrays</small></span
          ></label
        >
      </fieldset>
      <div
        v-if="['reconstruction', 'orientation'].includes(mode)"
        class="option-column"
      >
        <label for="edit-distance"
          ><strong>Spacer edit distance</strong><code>--cluster_spacers_max_distance</code
          ><small>Maximum sequence edits when matching homologous spacers</small></label
        >
        <div class="number-control">
          <button
            type="button"
            aria-label="Decrease spacer edit distance"
            @click="update({ spacerEditDistance: Math.max(0, modelValue.spacerEditDistance - 1) })"
          >
            <AppIcon
              name="minus"
              :size="16"
            /></button
          ><output
            id="edit-distance"
            aria-live="polite"
            >{{ modelValue.spacerEditDistance }}</output
          ><button
            type="button"
            aria-label="Increase spacer edit distance"
            @click="update({ spacerEditDistance: Math.min(2, modelValue.spacerEditDistance + 1) })"
          >
            <AppIcon
              name="plus"
              :size="16"
            />
          </button>
        </div>
      </div>
      <div
        v-if="['reconstruction', 'orientation'].includes(mode)"
        class="option-column"
      >
        <span
          ><strong>Deletion-parameter bias corrections</strong
          ><small
            >Refine α/ρ deletion-parameter estimates; orientation ΔlnL is unchanged. On this public
            service, if a group/direction-specific ρ correction has no finite fit, that fit alone
            uses uncorrected ρ estimates and records a warning instead of failing the job. The
            affected deletion-rate estimates remain uncorrected. Optimizer likelihoods and the
            IDM/BDM LRT are unchanged by ρ correction; strict all-corrections-or-error is a separate
            CLI policy.</small
          ></span
        >
        <div
          class="segmented"
          role="group"
          aria-label="Deletion-parameter bias corrections"
        >
          <button
            type="button"
            :aria-pressed="modelValue.biasCorrection"
            @click="update({ biasCorrection: true })"
          >
            Enabled</button
          ><button
            type="button"
            :aria-pressed="!modelValue.biasCorrection"
            @click="update({ biasCorrection: false })"
          >
            Disabled
          </button>
        </div>
      </div>
      <div
        v-if="['detection', 'reconstruction', 'orientation'].includes(mode)"
        class="policy-note"
      >
        <strong>{{ mode === "detection" ? "Detection only" : "Tree policy" }}</strong
        ><span>{{
          mode === "detection"
            ? "Spacer matching, tree estimation, and evolutionary modeling are not run."
            : "Trees are estimated separately for input and reversed order. The result reports the exact policy used."
        }}</span>
      </div>
    </div>
  </details>
</template>
