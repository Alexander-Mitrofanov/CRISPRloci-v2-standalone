<script setup>
import { useSequenceHandoff } from "../../../features/results/index.js";

const props = defineProps({
  job: { type: Object, required: true },
  credential: { type: Object, required: true },
  service: { type: Object, required: true },
  limits: { type: Object, required: true },
});
const emit = defineEmits(["continue-analysis"]);
const { choices, selected, loading, error, prepare } = useSequenceHandoff(props, emit);
</script>

<template>
  <section
    v-if="choices.length"
    class="result-section sequence-handoff"
    aria-labelledby="continue-heading"
  >
    <h3 id="continue-heading">Continue with these sequences</h3>
    <p>
      Load a complete output into another analysis. Review its input limits and options before
      submitting.
    </p>
    <label for="next-analysis">Sequence file and next analysis</label>
    <select
      id="next-analysis"
      v-model="selected"
      :disabled="loading"
    >
      <option value="">Choose a connection</option>
      <option
        v-for="choice in choices"
        :key="choice.id"
        :value="choice.id"
      >
        {{ choice.label }}
      </option>
    </select>
    <button
      type="button"
      class="secondary-button"
      :disabled="!selected || loading"
      @click="prepare"
    >
      {{ loading ? "Loading input…" : "Use as analysis input" }}
    </button>
    <p
      v-if="error"
      role="alert"
    >
      {{ error }}
    </p>
  </section>
</template>

<style scoped>
.sequence-handoff label {
  display: block;
  margin-block: 16px 8px;
  font-weight: 600;
  font-size: 0.875rem;
}
.sequence-handoff select {
  display: block;
  width: 100%;
  min-width: 0;
  min-height: 44px;
  padding: 10px;
  border: 1px solid var(--line-dark);
  border-radius: var(--radius-sm);
  background: var(--paper);
  color: var(--ink);
  font: inherit;
  font-size: 0.875rem;
}
.sequence-handoff button {
  margin-top: 12px;
  padding: 10px 16px;
  border: 1px solid var(--action);
  border-radius: var(--radius-sm);
  font-size: 0.875rem;
  font-weight: 600;
}
.sequence-handoff button:enabled {
  background: var(--action-soft);
  color: var(--action-text);
}
.sequence-handoff button:enabled:hover {
  background: var(--action-soft-hover);
}
.sequence-handoff [role="alert"] {
  color: var(--danger);
}
</style>
