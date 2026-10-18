<script setup>
import { computed } from "vue";

import { formatNumber } from "../../../utils/formatting.js";

const props = defineProps({ summary: { type: Object, required: true } });
const adapter = computed(() => props.summary.adapter);
const reasons = computed(() => {
  const value = adapter.value.skipped_by_reason;
  return Array.isArray(value)
    ? value.map((item) => [item.reason || item.code || "Other", item.count ?? 1])
    : Object.entries(value || {});
});
const metrics = computed(() => [
  ["Retained arrays", adapter.value.emitted_array_count],
  ["Excluded arrays", adapter.value.skipped_array_count],
  ["Eligible groups", adapter.value.emitted_group_count],
  ["Unknown strand excluded", adapter.value.unknown_strand_excluded_count],
]);
</script>

<template>
  <section
    v-if="Object.keys(adapter).length"
    class="result-section"
    aria-labelledby="preflight-heading"
  >
    <div class="result-heading">
      <div>
        <h3 id="preflight-heading">What reached the model</h3>
      </div>
    </div>
    <div class="metric-grid">
      <div
        v-for="([label, value], index) in metrics"
        :key="label"
        :class="{ attention: index === 3 && Number(value) > 0 }"
      >
        <span>{{ label }}</span
        ><strong>{{ formatNumber(value) }}</strong>
      </div>
    </div>
    <div class="reason-block">
      <h4>Exclusion and skip reasons</h4>
      <p
        v-if="!reasons.length"
        class="muted"
      >
        No exclusions were reported.
      </p>
      <ul
        v-else
        class="reason-list"
      >
        <li
          v-for="[reason, count] in reasons"
          :key="reason"
        >
          <span>{{ String(reason).replaceAll("_", " ") }}</span
          ><b>{{ formatNumber(count) }}</b>
        </li>
      </ul>
    </div>
  </section>
</template>
