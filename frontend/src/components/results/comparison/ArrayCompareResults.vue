<script setup>
import { computed, ref } from "vue";
const props = defineProps({
  available: Boolean,
  arrays: { type: Array, required: true },
  edges: { type: Array, required: true },
  isolates: { type: Array, required: true },
  counts: { type: Object, required: true },
  identity: { type: String, default: null },
  minShared: { type: Number, default: null },
  truncated: Boolean,
  compact: Boolean,
});
const selected = ref(null);
const membership = computed(() => {
  const counts = new Map();
  for (const array of props.arrays)
    for (const symbol of new Set(array.spacers.map((s) => s.symbol).filter(Boolean)))
      counts.set(symbol, (counts.get(symbol) || 0) + 1);
  return counts;
});
const colour = (s) =>
  s.symbol && membership.value.get(s.symbol) > 1 ? `#${s.symbol.slice(2, 8)}` : "#202c31";
const nodes = computed(() =>
  props.arrays.map((a, i) => ({
    ...a,
    n: i + 1,
    x: 250 + 170 * Math.cos((2 * Math.PI * i) / props.arrays.length),
    y: 190 + 145 * Math.sin((2 * Math.PI * i) / props.arrays.length),
  })),
);
const nodeMap = computed(() => new Map(nodes.value.map((n) => [n.id, n])));
const visibleArrays = computed(() => (props.compact ? props.arrays.slice(0, 3) : props.arrays));
const name = (id) => {
  const a = nodeMap.value.get(id);
  return a ? `${a.sourceId} · ${a.arrayId}` : id;
};
const show = (v) => v ?? "Not reported";
</script>
<template>
  <section
    class="result-section array-comparison"
    aria-label="Observed array comparison"
  >
    <h3>Observed array comparison</h3>
    <p v-if="!available">A completed observed-array comparison was not reported.</p>
    <template v-else>
      <p>
        <strong>{{ show(counts.arrays) }} arrays</strong> · {{ show(counts.edges) }} sharing links ·
        {{ show(counts.isolates) }} isolated arrays · {{ show(counts.deletions) }} deletion slots
      </p>
      <p>
        Equality:
        {{
          identity === "reverse-complement"
            ? "explicit reverse-complement equivalence"
            : identity === "exact"
              ? "exact DNA sequence"
              : "not reported"
        }}. Links require at least {{ show(minShared) }} distinct shared spacers. Jaccard uses sets;
        repeated occurrences keep their positions.
      </p>
      <p>
        Order follows the supplied source-forward coordinates. Unknown orientation stays unknown.
        Sharing does not establish ancestry, horizontal transfer or an evolutionary tree.
      </p>
      <p>
        Uploaded source descriptors and coordinates remain annotations; this comparison does not
        independently verify the original genomes.
      </p>
      <p
        v-if="truncated"
        role="note"
      >
        <strong>Preview truncated.</strong> Up to 50 arrays, 50 slots per array and 100 links are
        shown. Download the complete JSON and sharing TSV from Files &amp; methods.
      </p>
      <p v-if="compact">
        Open Array comparison for the sharing graph, complete preview and pairwise metrics.
      </p>
      <div
        v-if="!compact && nodes.length"
        class="sharing-graph"
      >
        <svg
          viewBox="0 0 500 380"
          role="group"
          aria-label="Observed spacer sharing graph; numbered nodes correspond to arrays below"
        >
          <line
            v-for="edge in edges"
            :key="`${edge.a}/${edge.b}`"
            :x1="nodeMap.get(edge.a).x"
            :y1="nodeMap.get(edge.a).y"
            :x2="nodeMap.get(edge.b).x"
            :y2="nodeMap.get(edge.b).y"
            stroke="#9eafb6"
            stroke-width="1.5"
          />
          <g
            v-for="node in nodes"
            :key="node.id"
            role="button"
            tabindex="0"
            :aria-label="`Highlight ${name(node.id)}`"
            @click="selected = selected === node.id ? null : node.id"
            @keydown.enter="selected = selected === node.id ? null : node.id"
            @keydown.space.prevent="selected = selected === node.id ? null : node.id"
          >
            <title>{{ name(node.id) }}{{ isolates.includes(node.id) ? " · isolated" : "" }}</title>
            <circle
              :cx="node.x"
              :cy="node.y"
              r="12"
              :fill="
                selected === node.id
                  ? '#123e47'
                  : isolates.includes(node.id)
                    ? '#e5e9eb'
                    : '#bdd8d1'
              "
              stroke="#45666d"
            />
            <text
              :x="node.x"
              :y="node.y + 4"
              text-anchor="middle"
              font-size="10"
              :fill="selected === node.id ? '#ffffff' : '#143d43'"
            >
              {{ node.n }}
            </text>
          </g>
        </svg>
      </div>
      <p>
        Colored blocks share a DNA identity across shown slots in different arrays; dark blocks have
        no shared identity in this preview. × denotes a reported deletion, never a DNA sequence.
      </p>
      <article
        v-for="array in visibleArrays"
        :key="array.id"
        class="repeat-instance observed-array"
        :class="{ selected: selected === array.id }"
      >
        <h4>{{ nodeMap.get(array.id)?.n }}. {{ array.sourceId }} · {{ array.arrayId }}</h4>
        <p>
          Strand: <strong>{{ array.strand ?? "Unknown" }}</strong> ·
          {{ show(array.observedCount) }} observed occurrences{{
            isolates.includes(array.id) ? " · isolated" : ""
          }}{{
            array.evidenceStatus === "no_observed_spacers" ? " · no observed DNA evidence" : ""
          }}
        </p>
        <div
          class="spacer-strip"
          role="img"
          :aria-label="`Ordered spacer slots for ${array.arrayId}`"
        >
          <svg
            v-for="spacer in array.spacers"
            :key="spacer.ordinal"
            viewBox="0 0 24 32"
            width="24"
            height="32"
          >
            <title>
              {{ spacer.deletion ? "Deletion slot" : spacer.id }} · position {{ spacer.ordinal
              }}{{ spacer.sequence ? ` · ${spacer.sequence}` : "" }}
            </title>
            <text
              v-if="spacer.deletion"
              x="12"
              y="17"
              text-anchor="middle"
              font-size="20"
              fill="#666666"
            >
              ×
            </text>
            <rect
              v-else
              x="2"
              y="3"
              width="20"
              height="15"
              rx="2"
              :fill="colour(spacer)"
            />
            <text
              x="12"
              y="30"
              text-anchor="middle"
              font-size="9"
              fill="#36454d"
            >
              {{ spacer.ordinal }}
            </text>
          </svg>
        </div>
        <p v-if="array.truncated">
          First 50 slots shown; remaining occurrences are in the complete JSON.
        </p>
        <details v-if="!compact">
          <summary>Occurrence identities, coordinates and DNA</summary>
          <p
            v-for="spacer in array.spacers"
            :key="spacer.ordinal"
            class="occurrence-detail"
          >
            <strong>{{ spacer.ordinal }} · {{ spacer.id }}</strong> ·
            {{ spacer.deletion ? "Reported deletion" : `${show(spacer.start)}–${show(spacer.end)}`
            }}<br />{{ spacer.sequence ?? "No DNA sequence" }}
          </p>
        </details>
      </article>
      <div
        v-if="!compact"
        class="comparison-table"
      >
        <h4>Pairwise sharing</h4>
        <p v-if="!edges.length && counts.edges === 0">
          No links meet the selected threshold. Isolated arrays remain visible.
        </p>
        <p v-else-if="!edges.length">
          No sharing links in this preview; download the full comparison.
        </p>
        <table v-else>
          <thead>
            <tr>
              <th>First array</th>
              <th>Second array</th>
              <th>Shared</th>
              <th>Jaccard</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="edge in edges"
              :key="`${edge.a}/${edge.b}`"
            >
              <td>{{ name(edge.a) }}</td>
              <td>{{ name(edge.b) }}</td>
              <td>{{ edge.shared }}</td>
              <td>{{ edge.jaccard.toFixed(3) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
    </template>
  </section>
</template>
<style scoped>
.array-comparison,
.observed-array {
  min-width: 0;
  overflow-wrap: anywhere;
}
.sharing-graph {
  max-width: 34rem;
  margin-inline: auto;
}
.sharing-graph svg {
  display: block;
  width: 100%;
  height: auto;
}
.sharing-graph g[role="button"] {
  cursor: pointer;
}
.sharing-graph g:focus circle {
  stroke-width: 3;
}
.observed-array.selected {
  outline: 2px solid #45666d;
  outline-offset: 3px;
}
.spacer-strip {
  display: flex;
  flex-wrap: wrap;
  gap: 2px;
  padding-block: 0.5rem;
}
.spacer-strip svg {
  flex: 0 0 24px;
}
.occurrence-detail {
  font-family: var(--font-mono, monospace);
  font-size: 0.8rem;
}
.comparison-table {
  overflow-x: auto;
}
.comparison-table table {
  width: 100%;
}
.comparison-table th,
.comparison-table td {
  padding: 0.6rem;
  text-align: left;
}
@media (max-width: 480px) {
  .comparison-table table {
    table-layout: fixed;
    font-size: 0.75rem;
  }
  .comparison-table th,
  .comparison-table td {
    padding: 0.35rem 0.2rem;
  }
  .comparison-table th:nth-child(-n + 2) {
    width: 32%;
  }
}
</style>
