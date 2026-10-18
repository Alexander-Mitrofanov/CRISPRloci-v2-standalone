<script setup>
import { computed } from "vue";

import { asArray, finiteMetric, formatNumber } from "../../../utils/formatting.js";
import { comparisonDecisionFor } from "../../../utils/results.js";

const props = defineProps({
  summary: { type: Object, required: true },
  exampleSnapshot: { type: Object, default: null },
  noEligible: Boolean,
  hasEvolution: { type: Boolean, default: true },
  hasOrientation: { type: Boolean, default: true },
});
const cards = computed(() => {
  const detection = props.summary.detection;
  const orientation = props.summary.orientation || {};
  const comparison = asArray(orientation.comparisons)[0];
  const reconstruction = (
    asArray(orientation.selected_reconstructions).length
      ? asArray(orientation.selected_reconstructions)
      : asArray(props.summary?.reconstruction?.results)
  )[0];
  const categories = detection.category_counts;
  const arrayCount =
    finiteMetric(detection.array_count) ??
    (Object.keys(categories).length
      ? Object.values(categories).reduce((total, value) => total + (finiteMetric(value) ?? 0), 0)
      : asArray(detection.arrays).length || null);
  const decision = comparison ? comparisonDecisionFor(comparison, orientation).label : null;
  const acquisitions = finiteMetric(reconstruction?.acquisitions);
  const deletions = finiteMetric(reconstruction?.deletions);
  return [
    {
      label: "Detected arrays",
      value: formatNumber(arrayCount),
      detail: `${formatNumber(categories["Bona-fide"])} Bona-fide`,
    },
    {
      label: "Orientation",
      visible: props.hasOrientation && !props.noEligible,
      value: decision || "Not evaluated",
      detail: "First reported group",
    },
    {
      label: "Reported history",
      visible: props.hasEvolution && !props.noEligible,
      value:
        acquisitions == null || deletions == null
          ? "Not reconstructed"
          : `${formatNumber(acquisitions)} / ${formatNumber(deletions)}`,
      detail: "acquisitions / deletions",
    },
  ].filter((card) => card.visible !== false);
});
</script>

<template>
  <section
    :class="['result-synopsis', { 'precomputed-synopsis': exampleSnapshot?.example }]"
    aria-labelledby="synopsis-heading"
  >
    <div class="synopsis-copy">
      <h3 id="synopsis-heading">
        {{ exampleSnapshot ? "Five related CRISPR arrays" : "Run summary" }}
      </h3>
      <p v-if="exampleSnapshot">Precomputed example. No sequences were submitted.</p>
      <p v-else-if="hasEvolution && !noEligible">
        Evolutionary estimates summarize the first reported group.
      </p>
    </div>
    <div class="synopsis-cards">
      <div
        v-for="card in cards"
        :key="card.label"
      >
        <span>{{ card.label }}</span
        ><strong>{{ card.value }}</strong
        ><small>{{ card.detail }}</small>
      </div>
    </div>
    <details
      v-if="exampleSnapshot?.example"
      class="example-context"
    >
      <summary>About this example</summary>
      <p>{{ exampleSnapshot.example.analysis_question }}</p>
      <p>{{ exampleSnapshot.example.analysis_takeaway }}</p>
    </details>
  </section>
</template>
