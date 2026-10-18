<script setup>
import { computed } from "vue";

import { formatNumber } from "../../../utils/formatting.js";
import { preferredDeletionModel, reconstructionMetric } from "../../../utils/reconstruction.js";

const props = defineProps({ row: { type: Object, required: true } });
const reportedStatistic = computed(() => reconstructionMetric(props.row, "lrt"));
const statistic = computed(() => Math.max(0, reportedStatistic.value ?? 0));
const cutoff = computed(() => reconstructionMetric(props.row, "cutoff"));
const preferred = computed(() => preferredDeletionModel(props.row));
const normalizedPreferred = computed(() => preferred.value.toUpperCase());
const extent = computed(() =>
  Math.max(1, cutoff.value == null ? 0 : cutoff.value * 1.45, statistic.value * 1.18),
);
const valuePosition = computed(() => Math.min(100, (statistic.value / extent.value) * 100));
const cutoffPosition = computed(() =>
  cutoff.value == null ? null : Math.min(100, (cutoff.value / extent.value) * 100),
);
const aria = computed(() =>
  [
    `Deletion model likelihood-ratio statistic ${formatNumber(statistic.value, 3)}`,
    cutoff.value == null ? "" : `, cutoff ${formatNumber(cutoff.value, 3)}`,
    `. Preferred model ${preferred.value}.`,
  ].join(""),
);
</script>

<template>
  <div class="model-selection">
    <div class="graphic-label">
      <span>Deletion-pattern model</span><small>IDM versus BDM</small>
    </div>
    <div class="model-call">
      <strong>{{ preferred }}</strong
      ><span>{{
        normalizedPreferred === "BDM"
          ? "Block deletion model favored by this LRT"
          : normalizedPreferred === "IDM"
            ? "BDM not favored by this LRT"
            : "Reported model"
      }}</span>
    </div>
    <template v-if="reportedStatistic != null">
      <div
        class="model-gauge"
        role="img"
        :aria-label="aria"
      >
        <svg
          viewBox="0 0 100 24"
          preserveAspectRatio="none"
          aria-hidden="true"
        >
          <rect
            class="model-gauge-track"
            x="0"
            y="0"
            width="100"
            height="24"
          />
          <rect
            class="model-gauge-fill"
            x="0"
            y="0"
            :width="valuePosition"
            height="24"
          />
          <line
            v-if="cutoffPosition != null"
            class="model-cutoff"
            :x1="cutoffPosition"
            :x2="cutoffPosition"
            y1="0"
            y2="24"
          />
          <line
            class="model-value"
            :x1="valuePosition"
            :x2="valuePosition"
            y1="0"
            y2="24"
          />
        </svg>
      </div>
      <div class="model-axis"><span>IDM retained</span><span>Evidence for BDM</span></div>
    </template>
    <p
      v-else
      class="model-unavailable"
    >
      Likelihood-ratio statistic not reported.
    </p>
    <div class="model-likelihoods">
      <span
        >LRT statistic
        <b>{{ formatNumber(reportedStatistic == null ? null : statistic, 3) }}</b></span
      ><span
        >Cutoff <b>{{ formatNumber(cutoff, 3) }}</b></span
      ><span
        >IDM lnL <b>{{ formatNumber(reconstructionMetric(row, "idmLikelihood"), 3) }}</b></span
      ><span
        >BDM lnL <b>{{ formatNumber(reconstructionMetric(row, "bdmLikelihood"), 3) }}</b></span
      >
    </div>
    <p class="model-interpretation">
      Retaining IDM is not proof that deletions are biologically independent. Rates are conditional
      on this model tree and its branch scale, not per-generation measurements.
    </p>
  </div>
</template>
