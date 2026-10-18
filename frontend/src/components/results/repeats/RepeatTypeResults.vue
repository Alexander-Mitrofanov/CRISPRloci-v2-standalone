<script setup>
import { computed, ref } from "vue";
const props = defineProps({
  available: Boolean,
  modelId: { type: String, default: null },
  modelSha256: { type: String, default: null },
  labels: { type: Array, required: true },
  threshold: { type: Number, default: null },
  counts: { type: Object, required: true },
  results: { type: Array, required: true },
  truncated: Boolean,
  compact: Boolean,
});
const filter = ref("");
const rows = computed(() =>
  props.results.filter((row) =>
    `${row.id} ${row.sequence} ${row.status} ${row.winner ?? ""}`
      .toLowerCase()
      .includes(filter.value.toLowerCase()),
  ),
);
const display = (value) => value ?? "Not reported";
</script>
<template>
  <section
    class="result-section repeat-type"
    aria-label="Repeat subtype evidence"
  >
    <h3>Repeat subtype evidence</h3>
    <p v-if="!available">A completed repeat classification was not reported.</p>
    <template v-else>
      <p>
        {{ display(counts.predicted) }} predicted · {{ display(counts.uncertain) }} uncertain ·
        {{ display(counts.unsupported) }} unsupported
      </p>
      <p>
        Native scores describe a closed set of {{ labels.length }} labels; they are not calibrated
        confidence. Unfamiliar families can receive high scores. Repeat evidence does not establish
        a Cas system or transcription strand.
      </p>
      <p>
        Minimum model score: {{ display(threshold) }}. Ambiguous IUPAC bases remain unsupported.
      </p>
      <p v-if="compact">Open Repeat subtypes for individual occurrences and model provenance.</p>
      <template v-else>
        <label
          >Filter displayed repeats
          <input
            v-model="filter"
            type="search"
        /></label>
        <p v-if="truncated">
          This preview is limited. Download complete JSON or TSV from Files &amp; methods.
        </p>
        <p v-if="!rows.length">No displayed repeats contain that text.</p>
        <article
          v-for="row in rows"
          :key="row.id"
          class="repeat-instance"
        >
          <h4>{{ row.id }} · {{ row.status ?? "Not reported" }}</h4>
          <p>
            <code class="sequence-text">{{ row.sequence }}</code>
          </p>
          <p>
            Prediction: {{ row.prediction ?? "Unknown" }} · Model winner:
            {{ row.winner ?? "Not assessed" }} · Native score:
            {{ row.nativeScore == null ? "Not reported" : row.nativeScore.toFixed(4) }}
          </p>
          <p v-if="row.reason === 'below_threshold'">
            The winning label falls below the selected score threshold.
          </p>
          <p v-if="row.reason === 'ambiguous_iupac'">
            This sequence contains ambiguous bases; no prediction was made.
          </p>
        </article>
        <details>
          <summary>Model provenance and supported labels</summary>
          <p>{{ display(modelId) }}</p>
          <p class="sequence-text">SHA-256: {{ display(modelSha256) }}</p>
          <p>{{ labels.join(", ") }}</p>
        </details>
      </template>
    </template>
  </section>
</template>
<style scoped>
.sequence-text {
  overflow-wrap: anywhere;
}
</style>
