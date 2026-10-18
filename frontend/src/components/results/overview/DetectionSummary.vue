<script setup>
import { computed, ref } from "vue";

import { categoryClass } from "../../../science.js";
import { formatNumber } from "../../../utils/formatting.js";
import AppIcon from "../../common/AppIcon.vue";
import ArraySequences from "./ArraySequences.vue";

const props = defineProps({
  summary: { type: Object, required: true },
  arrays: { type: Array, default: () => [] },
});
const showScores = ref(false);
const categories = computed(() => {
  const reported = props.summary.category_counts;
  if (Object.keys(reported).length) return reported;
  return props.arrays.reduce((counts, row) => {
    const label = row.category || "Unclassified";
    counts[label] = (counts[label] || 0) + 1;
    return counts;
  }, {});
});
const ordered = computed(() =>
  [
    ...new Set([
      "Bona-fide",
      "Possible",
      "Possible discarded",
      "Low score",
      ...Object.keys(categories.value),
    ]),
  ].filter((key) => categories.value[key] != null),
);
</script>

<template>
  <section
    class="result-section category-section"
    aria-labelledby="category-heading"
  >
    <div class="result-heading">
      <div>
        <h3 id="category-heading">CRISPR arrays</h3>
      </div>
      <p>CRISPRidentify categories</p>
    </div>
    <div
      v-if="ordered.length"
      class="category-grid"
    >
      <div
        v-for="label in ordered"
        :key="label"
        :class="['category-card', `category-${categoryClass(label)}`]"
      >
        <span>{{ label }}</span
        ><strong>{{ formatNumber(categories[label]) }}</strong
        ><small>arrays</small>
      </div>
    </div>
    <div
      v-else
      class="empty-result"
    >
      No CRISPR array calls were reported.
    </div>
    <div
      v-if="arrays.length"
      class="table-wrap"
      role="region"
      tabindex="0"
      aria-label="Scrollable CRISPRidentify array results"
    >
      <table>
        <thead>
          <tr>
            <th>Record</th>
            <th>Coordinates</th>
            <th>Category</th>
            <th>Strand</th>
            <th>Spacers</th>
            <th v-if="showScores">Raw CRISPRidentify Model score</th>
          </tr>
        </thead>
        <tbody>
          <tr
            v-for="(row, index) in arrays"
            :key="`${row.source_id}:${row.array_id}:${index}`"
          >
            <td>
              <strong>{{ row.source_id || "—" }}</strong
              ><small v-if="showScores">{{ row.array_id || "" }}</small>
            </td>
            <td>
              {{
                row.start != null && row.end != null
                  ? `${formatNumber(row.start)}–${formatNumber(row.end)}`
                  : "—"
              }}
            </td>
            <td>
              <span :class="['category-pill', `category-${categoryClass(row.category)}`]">{{
                row.category || "Unclassified"
              }}</span>
            </td>
            <td>{{ row.strand || "Unknown" }}</td>
            <td>{{ formatNumber(row.spacer_count) }}</td>
            <td v-if="showScores">
              <span class="raw-score">{{ formatNumber(row.model_score, 4) }}</span>
            </td>
          </tr>
        </tbody>
      </table>
    </div>
    <label
      v-if="arrays.length"
      class="result-score-toggle"
      ><input
        v-model="showScores"
        type="checkbox"
      />
      Show model scores</label
    >
    <ArraySequences
      v-for="(row, index) in arrays"
      :key="`${row.source_id}:${row.array_id}:${index}`"
      :label="`${row.source_id || 'Record'} / ${row.array_id || `array ${index + 1}`}`"
      :repeats="row.repeats || []"
      :spacers="row.spacers || []"
      :truncated="row.sequences_truncated"
    />
    <div
      v-if="showScores"
      class="interpretation-note"
    >
      <AppIcon
        name="info"
        :size="18"
      />
      <p>
        The raw Model score is <strong>not a calibrated probability</strong>. Interpret it with the
        category and array context.
      </p>
    </div>
  </section>
</template>
