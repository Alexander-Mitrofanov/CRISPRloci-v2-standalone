<script setup>
import { computed } from "vue";

import { finiteMetric, formatNumber, signedNumber } from "../../../utils/formatting.js";
import { groupIdentity } from "../../../utils/results.js";

const props = defineProps({
  group: { type: Object, required: true },
  index: { type: Number, required: true },
  decision: { type: Object, required: true },
});
const forward = computed(() => finiteMetric(props.group.forward_ln_likelihood_bdm));
const reverse = computed(() => finiteMetric(props.group.reverse_ln_likelihood_bdm));
const gap = computed(() =>
  props.decision.delta == null
    ? null
    : Math.max(0, props.decision.threshold - Math.abs(props.decision.delta)),
);
const surplus = computed(() =>
  props.decision.delta == null
    ? null
    : Math.max(0, Math.abs(props.decision.delta) - props.decision.threshold),
);
</script>

<template>
  <article :class="['hypothesis-comparison', `group-tone-${index % 4}`]">
    <div class="hypothesis-heading">
      <div>
        <small>Array group {{ index + 1 }}</small
        ><strong>{{ groupIdentity(group, index) }}</strong>
      </div>
      <span class="orientation-chip">{{ decision.label }}</span>
    </div>
    <div
      class="hypothesis-pair"
      :aria-label="`Likelihood comparison for ${groupIdentity(group, index)}`"
    >
      <div
        :class="[
          'hypothesis-card',
          'hypothesis-input',
          { 'is-leading': forward != null && reverse != null && forward > reverse },
        ]"
      >
        <span>Input spacer order</span><strong>{{ formatNumber(forward, 3) }}</strong
        ><small>BDM log likelihood</small>
      </div>
      <div class="hypothesis-versus">
        <span>vs</span
        ><b>{{
          decision.delta == null ? "No delta" : `Delta ${signedNumber(decision.delta, 2)}`
        }}</b>
      </div>
      <div
        :class="[
          'hypothesis-card',
          'hypothesis-reverse',
          { 'is-leading': forward != null && reverse != null && reverse > forward },
        ]"
      >
        <span>Reversed spacer order</span><strong>{{ formatNumber(reverse, 3) }}</strong
        ><small>BDM log likelihood</small>
      </div>
    </div>
    <div
      :class="[
        'decision-distance',
        decision.label === 'Unresolved' ? 'is-unresolved' : 'is-decisive',
      ]"
    >
      <span>{{
        decision.label === "Unresolved" ? "Distance still needed" : "Boundary crossed by"
      }}</span
      ><strong>{{ formatNumber(decision.label === "Unresolved" ? gap : surplus, 2) }} ΔlnL</strong
      ><small>Decision boundary: ±{{ formatNumber(decision.threshold, 2) }}</small>
    </div>
  </article>
</template>
