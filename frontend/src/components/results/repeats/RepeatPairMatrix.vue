<script setup>
import { computed, ref, useId, watch } from "vue";

const props = defineProps({
  sequence: { type: String, required: true },
  isolatedPairs: { type: Array, required: true },
  contextPairs: { type: Array, required: true },
  motifPairs: { type: Array, required: true },
  isolatedAvailable: Boolean,
  contextAvailable: Boolean,
});
const heading = useId();
const first = ref(0);
const second = ref(1);
const n = computed(() => props.sequence.length);
const cell = computed(() => 560 / Math.max(n.value, 1));
const ticks = computed(() => {
  const stride = Math.max(1, Math.ceil(n.value / 8));
  return Array.from({ length: n.value }, (_, index) => index).filter(
    (index) => index % stride === 0 || index === n.value - 1,
  );
});
const key = (i, j) => `${Math.min(i, j)}:${Math.max(i, j)}`;
const isolated = computed(
  () => new Map(props.isolatedPairs.map(([i, j, value]) => [key(i, j), value])),
);
const context = computed(
  () => new Map(props.contextPairs.map(([i, j, value]) => [key(i, j), value])),
);
const validSelection = computed(
  () =>
    Number.isInteger(first.value) &&
    Number.isInteger(second.value) &&
    first.value >= 0 &&
    second.value >= 0 &&
    first.value < n.value &&
    second.value < n.value &&
    first.value !== second.value,
);
const selected = computed(() => key(first.value, second.value));
const display = (value) => (value == null ? "Not available" : String(value));
const isolatedValue = computed(() =>
  props.isolatedAvailable && validSelection.value
    ? (isolated.value.get(selected.value) ?? 0)
    : null,
);
const contextValue = computed(() =>
  props.contextAvailable && validSelection.value ? (context.value.get(selected.value) ?? 0) : null,
);
watch(
  () => props.sequence,
  () => {
    first.value = props.motifPairs[0]?.[0] ?? 0;
    second.value = props.motifPairs[0]?.[1] ?? Math.min(1, n.value - 1);
  },
  { immediate: true },
);
</script>

<template>
  <section
    class="repeat-matrix"
    :aria-labelledby="heading"
  >
    <h4 :id="heading">Pair support within this repeat</h4>
    <p class="repeat-matrix-key">
      <span class="repeat-isolated-key">Upper triangle: isolated RNA</span>
      <span class="repeat-context-key">Lower triangle: local array context</span>
      <span>Outlined cells: isolated MFE motif</span>
    </p>
    <svg
      viewBox="0 0 626 626"
      role="img"
      :aria-labelledby="`${heading}-title ${heading}-description`"
    >
      <title :id="`${heading}-title`">Isolated and array-context pair-support matrix</title>
      <desc :id="`${heading}-description`">
        Both axes are zero-based repeat positions. Color intensity is pair support from zero to one.
        Upper triangle shows isolated probabilities; lower triangle shows local context pair
        averages. Outlines mark the proposed MFE motif. Inspect exact values with the coordinate
        controls below.
      </desc>
      <rect
        x="42"
        y="22"
        width="560"
        height="560"
        fill="#f0f4f5"
      />
      <path
        v-if="!contextAvailable"
        d="M42 22 L42 582 L602 582 Z"
        fill="#e4e9ec"
      />
      <rect
        v-for="[i, j, p] in isolatedPairs"
        :key="`p:${i}:${j}`"
        :x="42 + j * cell"
        :y="22 + i * cell"
        :width="cell"
        :height="cell"
        fill="#087b78"
        :fill-opacity="p"
      >
        <title>Isolated ({{ i }}, {{ j }}): {{ p }}</title>
      </rect>
      <rect
        v-for="[i, j, p] in contextPairs"
        :key="`q:${i}:${j}`"
        :x="42 + i * cell"
        :y="22 + j * cell"
        :width="cell"
        :height="cell"
        fill="#3158a5"
        :fill-opacity="p"
      >
        <title>Context ({{ i }}, {{ j }}): {{ p }}</title>
      </rect>
      <template
        v-for="[i, j] in motifPairs"
        :key="`motif:${i}:${j}`"
      >
        <rect
          :x="42 + j * cell"
          :y="22 + i * cell"
          :width="cell"
          :height="cell"
          fill="none"
          stroke="#20323d"
          stroke-width="0.8"
        />
        <rect
          v-if="contextAvailable"
          :x="42 + i * cell"
          :y="22 + j * cell"
          :width="cell"
          :height="cell"
          fill="none"
          stroke="#20323d"
          stroke-width="0.8"
        />
      </template>
      <line
        x1="42"
        y1="22"
        x2="602"
        y2="582"
        stroke="#98a9b1"
        stroke-width="1"
      />
      <g
        v-for="index in ticks"
        :key="index"
        fill="#465a65"
        font-size="13"
      >
        <text
          :x="42 + (index + 0.5) * cell"
          y="602"
          text-anchor="middle"
        >
          {{ index }}
        </text>
        <text
          x="34"
          :y="22 + (index + 0.5) * cell + 4"
          text-anchor="end"
        >
          {{ index }}
        </text>
      </g>
      <text
        x="322"
        y="622"
        text-anchor="middle"
        fill="#465a65"
        font-size="13"
      >
        Repeat position (zero-based)
      </text>
    </svg>
    <p
      v-if="!contextAvailable"
      class="repeat-note"
    >
      The lower triangle is unavailable because no completed context result was reported.
    </p>
    <div class="repeat-pair-inspector">
      <label :for="`${heading}-first`"
        >First base (0–{{ n - 1 }})
        <input
          :id="`${heading}-first`"
          v-model.number="first"
          type="number"
          min="0"
          :max="n - 1"
          step="1"
        />
      </label>
      <label :for="`${heading}-second`"
        >Second base (0–{{ n - 1 }})
        <input
          :id="`${heading}-second`"
          v-model.number="second"
          type="number"
          min="0"
          :max="n - 1"
          step="1"
        />
      </label>
      <dl aria-live="polite">
        <div>
          <dt>Isolated pair</dt>
          <dd>{{ display(isolatedValue) }}</dd>
        </div>
        <div>
          <dt>Context pair</dt>
          <dd>{{ display(contextValue) }}</dd>
        </div>
      </dl>
    </div>
    <p class="repeat-note">
      Local context values average different windows. Do not sum them into a partner probability or
      compare them as experimental accuracy.
    </p>
  </section>
</template>
