<script setup>
import { computed } from "vue";

import { formatNumber } from "../../../utils/formatting.js";
import { reconstructionMetric } from "../../../utils/reconstruction.js";

const props = defineProps({ row: { type: Object, required: true } });
const acquisitions = computed(() => reconstructionMetric(props.row, "acquisitions"));
const deletions = computed(() => reconstructionMetric(props.row, "deletions"));
const hasTally = computed(() => acquisitions.value != null && deletions.value != null);
const total = computed(() =>
  hasTally.value ? Math.max(0, acquisitions.value) + Math.max(0, deletions.value) : null,
);
const acquisitionShare = computed(() =>
  total.value > 0 ? (Math.max(0, acquisitions.value) / total.value) * 100 : null,
);
const deletionShare = computed(() =>
  acquisitionShare.value == null ? null : 100 - acquisitionShare.value,
);
const specialDefinitions = [
  ["duplication", "Duplication candidates", "duplications"],
  ["rearrangement", "Rearrangement candidates", "rearrangements"],
  ["reacquisition", "Reacquisition candidates", "reacquisitions"],
  ["independent", "Independent-gain candidates", "independentGains"],
];
const specialMetrics = computed(() =>
  specialDefinitions.map(([type, label, key]) => ({
    type,
    label,
    value: reconstructionMetric(props.row, key),
  })),
);
const specials = computed(() =>
  specialMetrics.value.filter((item) => item.value != null && item.value > 0),
);
const allSpecialsReported = computed(() =>
  specialMetrics.value.every((item) => item.value != null),
);
</script>

<template>
  <div class="event-graphic">
    <div class="graphic-label">
      <span>Reconstructed event totals</span><small>Reported ancestral history</small>
    </div>
    <div
      v-if="total > 0"
      class="event-ribbon"
      role="img"
      :aria-label="`${formatNumber(acquisitions)} acquisitions and ${formatNumber(deletions)} deletions`"
    >
      <svg
        class="event-ribbon-chart"
        viewBox="0 0 100 100"
        preserveAspectRatio="none"
        aria-hidden="true"
      >
        <rect
          class="event-ribbon-gains"
          x="0"
          y="0"
          :width="acquisitionShare"
          height="100"
        />
        <rect
          class="event-ribbon-losses"
          :x="acquisitionShare"
          y="0"
          :width="deletionShare"
          height="100"
        />
      </svg>
      <div
        class="event-ribbon-labels"
        aria-hidden="true"
      >
        <span class="event-ribbon-label-gains"
          ><b>{{ formatNumber(acquisitions) }}</b
          ><small>acquisitions</small></span
        ><span class="event-ribbon-label-losses"
          ><b>{{ formatNumber(deletions) }}</b
          ><small>deletions</small></span
        >
      </div>
    </div>
    <p
      v-else-if="hasTally"
      class="model-unavailable"
    >
      No acquisition or deletion events were reconstructed.
    </p>
    <p
      v-else
      class="model-unavailable"
    >
      Acquisition and deletion totals were not reported.
    </p>
    <div
      v-if="specials.length"
      class="special-event-grid"
    >
      <span
        v-for="item in specials"
        :key="item.type"
        ><i
          :class="['event-glyph', `event-glyph-${item.type}`]"
          aria-hidden="true"
        /><b>{{ formatNumber(item.value) }}</b
        ><small>{{ item.label }}</small></span
      >
    </div>
    <p
      v-else
      class="no-special-events"
    >
      {{
        allSpecialsReported
          ? "No duplication, rearrangement, reacquisition, or independent-gain candidates were reported."
          : "Special-event candidate totals were not reported."
      }}
    </p>
    <p>
      The visual key is adapted from SpacerPlacer: green denotes acquisitions, red denotes
      deletions, and distinct shapes flag special-event candidates when present.
    </p>
  </div>
</template>
