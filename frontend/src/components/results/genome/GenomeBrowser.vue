<script setup>
import { computed, ref, useId, watch } from "vue";
import { useGenomeBrowser } from "../../../features/genome/index.js";

const props = defineProps({
  sources: { type: Array, required: true },
  tracks: { type: Array, required: true },
  features: { type: Array, required: true },
  notes: { type: Array, required: true },
  active: { type: Boolean, required: true },
});
const emit = defineEmits(["open-evidence", "ready"]);
const id = useId();
const host = ref(null);
const sourceId = ref(props.sources[0].id);
const source = computed(() => props.sources.find((row) => row.id === sourceId.value));
const from = ref(1);
const to = ref(source.value.length);
const query = ref("");
const selectedKey = ref("");
const overlapping = ref([]);
const byKey = new Map(props.features.map((feature) => [feature.key, feature]));
const selected = computed(() => byKey.get(selectedKey.value));
const filtered = computed(() =>
  props.features.filter(
    (feature) =>
      feature.sourceId === sourceId.value &&
      `${feature.label} ${feature.kind}`.toLowerCase().includes(query.value.toLowerCase()),
  ),
);
const choices = computed(() => filtered.value.slice(0, 100));
const dna = computed(() =>
  selected.value && !selected.value.point
    ? source.value.sequence.slice(
        selected.value.start,
        Math.min(selected.value.end, selected.value.start + 2000),
      )
    : "",
);

function picked(keys) {
  overlapping.value = keys.map((key) => byKey.get(key)).filter(Boolean);
  selectedKey.value = overlapping.value[0]?.key || "";
  if (selected.value) sourceId.value = selected.value.sourceId;
}
const { ready, error, navigate, zoom, exportFigure } = useGenomeBrowser(props, host, picked);
watch(ready, (value) => emit("ready", value));

async function wholeRecord() {
  from.value = 1;
  to.value = source.value.length;
  selectedKey.value = "";
  overlapping.value = [];
  await navigate(source.value);
}

async function focusFeature(key) {
  selectedKey.value = key;
  if (!selected.value) return;
  sourceId.value = selected.value.sourceId;
  const padding = Math.max(40, Math.ceil((selected.value.end - selected.value.start) / 2));
  from.value = Math.max(1, selected.value.start + 1 - padding);
  to.value = Math.min(source.value.length, selected.value.end + padding);
  await navigate(source.value, from.value, to.value);
}
</script>

<template>
  <div class="genome-explorer">
    <div class="genome-navigation">
      <label :for="`${id}-source`"
        >Sequence record
        <select
          :id="`${id}-source`"
          v-model="sourceId"
          aria-label="Sequence record"
          :disabled="!ready"
          @change="wholeRecord"
        >
          <option
            v-for="row in sources"
            :key="row.id"
            :value="row.id"
          >
            {{ row.inputId }} · {{ row.length.toLocaleString() }} nt
          </option>
        </select>
      </label>
      <form
        class="genome-interval"
        @submit.prevent="navigate(source, Number(from), Number(to))"
      >
        <label :for="`${id}-from`"
          >Start
          <input
            :id="`${id}-from`"
            v-model="from"
            type="number"
            min="1"
            :max="source.length"
            required
          />
        </label>
        <label :for="`${id}-to`"
          >End
          <input
            :id="`${id}-to`"
            v-model="to"
            type="number"
            min="1"
            :max="source.length"
            required
          />
        </label>
        <button
          type="submit"
          :disabled="!ready"
        >
          Go to interval
        </button>
      </form>
    </div>
    <div class="genome-actions">
      <button
        type="button"
        :disabled="!ready"
        @click="zoom('in')"
      >
        Zoom in
      </button>
      <button
        type="button"
        :disabled="!ready"
        @click="zoom('out')"
      >
        Zoom out
      </button>
      <button
        type="button"
        :disabled="!ready"
        @click="wholeRecord"
      >
        Whole record
      </button>
      <button
        type="button"
        :disabled="!ready"
        @click="exportFigure(source)"
      >
        Export SVG figure
      </button>
      <button
        type="button"
        @click="$emit('open-evidence')"
      >
        Open evidence tables
      </button>
    </div>
    <p
      v-if="!ready && !error"
      role="status"
    >
      Preparing genome tracks…
    </p>
    <p
      v-if="error"
      role="alert"
    >
      {{ error }}
    </p>
    <p class="genome-caption">
      {{ source.chr }}: {{ source.inputId }} · {{ source.length.toLocaleString() }} nt. Topology:
      {{ source.topology }}.
      <span v-if="source.assessment"
        >Assessment: {{ source.assessment.replaceAll("_", " ") }}.</span
      >
      Drag the tracks to pan; select a feature to inspect its evidence.
    </p>
    <div
      class="genome-canvas-scroll"
      tabindex="0"
      role="region"
      aria-label="Interactive genome tracks"
    >
      <div
        ref="host"
        class="genome-canvas"
        @contextmenu.capture.prevent.stop
      />
    </div>
    <ul
      class="genome-legend"
      aria-label="Genome track legend"
    >
      <li
        v-for="track in tracks"
        :key="track.id"
      >
        <svg
          width="14"
          height="14"
          aria-hidden="true"
        >
          <rect
            width="14"
            height="14"
            :fill="track.color"
          />
        </svg>
        {{ track.label }} ({{
          track.features.filter((feature) => feature.sourceId === source.id).length
        }})
      </li>
    </ul>
    <div class="genome-selection">
      <div>
        <label :for="`${id}-filter`"
          >Find a feature
          <input
            :id="`${id}-filter`"
            v-model="query"
            type="search"
            placeholder="Array, gene, spacer or feature name"
          />
        </label>
        <label :for="`${id}-feature`"
          >Feature
          <select
            :id="`${id}-feature`"
            :value="selectedKey"
            aria-label="Feature"
            :disabled="!ready"
            @change="
              overlapping = [];
              focusFeature($event.target.value);
            "
          >
            <option value="">Choose a reported feature</option>
            <option
              v-if="selected && !choices.some((feature) => feature.key === selected.key)"
              :value="selected.key"
            >
              {{ selected.label }}
            </option>
            <option
              v-for="feature in choices"
              :key="feature.key"
              :value="feature.key"
            >
              {{ feature.label }}
            </option>
          </select>
        </label>
        <p v-if="filtered.length > choices.length">
          {{ choices.length }} of {{ filtered.length }} matching features listed. Refine the search;
          all located features remain in the tracks.
        </p>
        <p v-if="!filtered.length">No reported features match this selection.</p>
      </div>
      <section
        v-if="selected"
        class="genome-feature"
        aria-label="Selected genomic feature"
      >
        <h4>{{ selected.label }}</h4>
        <p>
          {{ selected.kind }} ·
          {{
            selected.point
              ? `Interbase boundary ${selected.start}`
              : `${(selected.start + 1).toLocaleString()}–${selected.end.toLocaleString()} (${selected.end - selected.start} nt)`
          }}
          · Strand: {{ selected.strand ?? "Not reported" }}
        </p>
        <div v-if="overlapping.length > 1">
          <p>{{ overlapping.length }} features overlap this position:</p>
          <button
            v-for="feature in overlapping"
            :key="feature.key"
            type="button"
            @click="focusFeature(feature.key)"
          >
            {{ feature.label }}
          </button>
        </div>
        <dl>
          <template
            v-for="([label, value], index) in selected.details"
            :key="index"
          >
            <dt>{{ label }}</dt>
            <dd>{{ value }}</dd>
          </template>
        </dl>
        <details v-if="dna">
          <summary>
            Reference DNA at this interval (source-forward){{
              selected.end - selected.start > 2000 ? " (first 2,000 bases)" : ""
            }}
          </summary>
          <p class="genome-dna">{{ dna }}</p>
        </details>
      </section>
      <p v-else>
        Select a track feature or use the feature list to inspect coordinates and evidence.
      </p>
    </div>
    <details class="genome-notes">
      <summary>Scope and interpretation</summary>
      <p
        v-for="note in notes"
        :key="note"
      >
        {{ note }}
      </p>
      <p>Only this job's reported tracks are shown. Unknown topology never connects contig ends.</p>
      <p>
        Genome navigation uses
        <a
          href="https://igv.org/doc/igvjs/"
          target="_blank"
          rel="noopener noreferrer"
          >igv.js</a
        >. Sequences are loaded from this job and rendered in your browser.
      </p>
    </details>
  </div>
</template>

<style scoped>
.genome-explorer {
  margin-top: 1.2rem;
}
.genome-navigation,
.genome-interval,
.genome-actions,
.genome-legend {
  display: flex;
  flex-wrap: wrap;
  gap: 0.7rem;
  align-items: end;
}
.genome-navigation > label {
  flex: 1;
  min-width: 15rem;
}
label {
  display: grid;
  gap: 0.3rem;
  font-weight: 500;
}
input,
select {
  min-height: 2.6rem;
  padding: 0.45rem;
  border: 1px solid var(--line-dark);
  background: var(--paper);
  color: var(--ink);
}
select {
  max-width: 100%;
  width: 100%;
}
.genome-interval input {
  width: 8rem;
}
.genome-actions {
  margin: 0.8rem 0;
}
button {
  padding: 0.5rem 0.75rem;
  border: 1px solid var(--line-dark);
  border-radius: 3px;
  background: var(--paper);
  color: var(--ink);
}
button:hover:not(:disabled) {
  background: var(--surface-subtle);
}
button:disabled {
  opacity: 0.5;
  cursor: wait;
}
.genome-caption {
  color: var(--ink-soft);
  font-size: 0.9rem;
}
.genome-canvas-scroll {
  overflow-x: auto;
  border-block: 1px solid var(--line);
  padding-block: 0.5rem;
}
.genome-canvas {
  min-width: 540px;
  min-height: 160px;
}
.genome-legend {
  list-style: none;
  padding: 0;
  align-items: center;
  font-size: 0.85rem;
}
.genome-legend li {
  display: flex;
  align-items: center;
  gap: 0.35rem;
}
.genome-selection {
  display: grid;
  grid-template-columns: minmax(200px, 1fr) minmax(0, 2fr);
  gap: 1.4rem;
  margin: 1.5rem 0;
}
.genome-selection label + label {
  margin-top: 0.6rem;
}
.genome-feature {
  border-left: 3px solid var(--green);
  padding-left: 1rem;
  overflow-wrap: anywhere;
}
.genome-feature h4 {
  margin-top: 0;
}
dl {
  display: grid;
  grid-template-columns: minmax(120px, 1fr) 2fr;
  gap: 0.3rem 0.7rem;
}
dt {
  color: var(--ink-soft);
}
dd {
  margin: 0;
}
.genome-dna {
  font-family: monospace;
  overflow-wrap: anywhere;
}
.genome-notes {
  border-top: 1px solid var(--line);
  padding-top: 0.8rem;
}
.genome-notes p {
  max-width: 80ch;
}
@media (max-width: 680px) {
  .genome-selection {
    grid-template-columns: 1fr;
  }
  .genome-interval input {
    width: 6.5rem;
  }
}
</style>
