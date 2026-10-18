<script setup>
import { computed } from "vue";

import { formatNumber } from "../../../utils/formatting.js";
import { reconstructionMetric } from "../../../utils/reconstruction.js";

const props = defineProps({ row: { type: Object, required: true } });
const unique = computed(() => reconstructionMetric(props.row, "unique"));
const aligned = computed(() => reconstructionMetric(props.row, "aligned"));
const leaves = computed(() => reconstructionMetric(props.row, "leaves"));
const patterns = computed(() => reconstructionMetric(props.row, "patterns"));
const visible = computed(() => Math.min(18, Math.max(0, Math.round(unique.value || 0))));
const hasInventory = computed(() => unique.value != null || aligned.value != null);
</script>

<template>
  <div
    v-if="hasInventory"
    class="spacer-inventory"
  >
    <div class="graphic-label">
      <span>Spacer inventory</span><small>Count view, not branch placement</small>
    </div>
    <div
      class="spacer-blocks"
      role="img"
      :aria-label="`${unique == null ? 'Unknown number of' : formatNumber(unique)} unique spacers in ${formatNumber(aligned)} aligned spacer positions; count view only`"
    >
      <i
        v-for="index in visible"
        :key="index"
        class="spacer-block spacer-block-neutral"
        ><span>{{ index }}</span></i
      ><b v-if="unique > visible">+{{ formatNumber(unique - visible) }}</b>
    </div>
    <div class="inventory-counts">
      <span
        ><strong>{{ formatNumber(unique) }}</strong> unique spacers</span
      ><span
        ><strong>{{ formatNumber(aligned) }}</strong> aligned positions</span
      ><span
        ><strong>{{ formatNumber(leaves) }}</strong> tree leaves</span
      ><span
        ><strong>{{ formatNumber(patterns) }}</strong> unique array patterns</span
      >
    </div>
  </div>
</template>
