<script setup>
import { computed, nextTick, ref, useId } from "vue";
import { catalogInterval, sequenceSearchIssue } from "../../../features/catalog/index.js";
import SequenceSearchPanel from "../../catalog/SequenceSearchPanel.vue";

const props = defineProps({
  repeats: { type: Array, required: true },
  spacers: { type: Array, required: true },
  label: { type: String, required: true },
  truncated: { type: Boolean, default: false },
});
const selected = ref(null);
const panelId = useId();
let trigger;
const sections = computed(() => [
  { kind: "repeats", label: "Repeat", rows: props.repeats },
  { kind: "spacers", label: "Spacer", rows: props.spacers },
]);
function search(section, row, event) {
  trigger = event.currentTarget;
  selected.value = {
    kind: section.kind,
    sequence: row.sequence,
    label: `${section.label} ${row.ordinal} in ${props.label}`,
  };
}
async function close() {
  selected.value = null;
  await nextTick();
  trigger?.focus();
}
</script>

<template>
  <details class="result-details array-sequences">
    <summary>Repeats &amp; spacers — {{ label }}</summary>
    <p
      v-if="truncated"
      role="status"
    >
      Some sequences were omitted from this summary. The complete results are in Files &amp;
      methods.
    </p>
    <section
      v-for="section in sections"
      :key="section.kind"
      :aria-label="`${section.label} sequences in ${label}`"
    >
      <h4>{{ section.kind === "repeats" ? "Observed repeats" : "Observed spacers" }}</h4>
      <div
        v-if="section.rows.length"
        class="table-wrap"
        tabindex="0"
        role="region"
        :aria-label="`${section.label} sequences`"
      >
        <table class="sequence-responsive-table">
          <thead>
            <tr>
              <th>Occurrence</th>
              <th>Coordinates</th>
              <th>Sequence</th>
              <th>Database</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="(row, index) in section.rows"
              :key="index"
            >
              <td data-label="Occurrence">{{ section.label }} {{ row.ordinal }}</td>
              <td data-label="Coordinates">{{ catalogInterval(row) }}</td>
              <td data-label="Sequence">
                <code class="sequence-search-candidate">{{ row.sequence || "Not reported" }}</code>
              </td>
              <td data-label="Database">
                <button
                  type="button"
                  class="catalog-button"
                  :disabled="Boolean(sequenceSearchIssue(row.sequence, row.sequence_status))"
                  :aria-label="`Search in DB: ${section.label.toLowerCase()} ${row.ordinal} in ${label}`"
                  :aria-controls="selected ? panelId : undefined"
                  @click="search(section, row, $event)"
                >
                  Search in DB
                </button>
                <small v-if="sequenceSearchIssue(row.sequence, row.sequence_status)">{{
                  sequenceSearchIssue(row.sequence, row.sequence_status)
                }}</small>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
      <p v-else>Individual {{ section.kind }} were not reported in this summary.</p>
    </section>
    <SequenceSearchPanel
      v-if="selected"
      :id="panelId"
      :key="selected.label"
      v-bind="selected"
      @close="close"
    />
  </details>
</template>
