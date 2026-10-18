<script setup>
import { computed } from "vue";

import { asArray, formatNumber } from "../../../utils/formatting.js";
import { reconstructionMetric } from "../../../utils/reconstruction.js";
import { comparisonDecisionFor, groupIdentity } from "../../../utils/results.js";
import AncestralHistoryExplorer from "../history/AncestralHistoryExplorer.vue";
import DeletionModelEvidence from "./DeletionModelEvidence.vue";
import ReconstructionEvents from "./ReconstructionEvents.vue";
import ReconstructionVerdict from "./ReconstructionVerdict.vue";
import SpacerInventory from "./SpacerInventory.vue";

const props = defineProps({
  row: { type: Object, required: true },
  index: { type: Number, required: true },
  summary: { type: Object, required: true },
});
const group = computed(() => String(props.row.name || `Group ${props.index + 1}`));
const orientation = computed(() => props.summary.orientation || {});
const comparison = computed(() =>
  asArray(orientation.value.comparisons).find(
    (item, index) => groupIdentity(item, index) === group.value,
  ),
);
const decision = computed(() =>
  comparison.value
    ? comparisonDecisionFor(comparison.value, orientation.value).label
    : "Input orientation; not compared",
);
const reportedByDefault = computed(() => decision.value === "Unresolved");
const tree = computed(() => {
  const entry = [
    ...asArray(orientation.value.trees),
    ...asArray(props.summary.reconstruction?.trees),
  ].find((item) => String(item?.group) === group.value);
  return entry?.selected_newick || entry?.newick || "";
});
</script>

<template>
  <article :class="['reconstruction-story', `group-tone-${index % 4}`]">
    <div class="reconstruction-story-heading">
      <div>
        <small>Reconstructed group {{ index + 1 }}</small>
        <h4>{{ group }}</h4>
        <p v-if="!comparison">
          Input orientation; not compared. Tree inferred from arrays with an uncalibrated scale.
        </p>
      </div>
      <span>{{ formatNumber(reconstructionMetric(row, "leaves")) }} leaves</span>
    </div>
    <AncestralHistoryExplorer
      :summary="summary"
      :group="group"
      :fallback-tree="tree"
      :reported-by-default="reportedByDefault"
    />
    <details class="model-details">
      <summary>Events &amp; model estimates</summary>
      <ReconstructionVerdict :row="row" />
      <div class="reconstruction-visual-grid">
        <ReconstructionEvents :row="row" /><DeletionModelEvidence :row="row" />
      </div>
      <SpacerInventory :row="row" />
    </details>
  </article>
</template>
