<script setup>
import { computed } from "vue";

const props = defineProps({
  available: Boolean,
  loci: { type: Array, default: () => [] },
  links: { type: Array, default: () => [] },
  counts: { type: Object, default: () => ({}) },
  identityCutoff: { type: Number, default: null },
  truncated: Boolean,
  compact: Boolean,
});
const shownLinks = computed(() => props.links.slice(0, props.compact ? 5 : 100));
const number = (value) => (value == null ? "Not reported" : value.toLocaleString());
const percent = (value) => (value == null ? "Not reported" : `${(100 * value).toFixed(1)}%`);
function arrow(gene, locus) {
  const span = locus.end - locus.start + 1;
  const left = 10 + (780 * (gene.start - locus.start)) / span;
  const right = Math.max(left + 2, 10 + (780 * (gene.end - locus.start + 1)) / span);
  const point = Math.min(5, (right - left) / 2);
  if (!gene.strand) return `${left},12 ${right},12 ${right},32 ${left},32`;
  return gene.strand === "-"
    ? `${left},22 ${left + point},12 ${right},12 ${right},32 ${left + point},32`
    : `${left},12 ${right - point},12 ${right},22 ${right - point},32 ${left},32`;
}
</script>

<template>
  <section
    class="locus-comparison"
    aria-label="Gene-neighborhood comparison"
  >
    <h3>Gene-neighborhood comparison</h3>
    <p>
      Protein similarity and observed gene order do not establish orthology, function or an organism
      phylogeny.
    </p>
    <p v-if="!available">Validated locus-comparison evidence is unavailable.</p>
    <template v-else>
      <p>
        {{ number(counts.loci) }} loci · {{ number(counts.genes) }} source-verified CDSs ·
        {{ number(counts.links) }} similarity links
      </p>
      <p>
        Global protein identity cutoff: {{ percent(identityCutoff) }}. Each neighborhood below has
        its own coordinate scale; arrows show source-forward positions and gene strand.
      </p>
      <article
        v-for="locus in loci"
        :key="locus.id"
        class="locus-row"
      >
        <h4>{{ locus.sourceId }} · {{ number(locus.start) }}–{{ number(locus.end) }}</h4>
        <svg
          v-if="locus.start !== null && locus.end >= locus.start"
          viewBox="0 0 800 44"
          role="img"
          :aria-label="`Gene map for ${locus.sourceId}`"
        >
          <line
            x1="10"
            y1="22"
            x2="790"
            y2="22"
            stroke="currentColor"
            stroke-opacity="0.3"
          />
          <template
            v-for="gene in locus.genes"
            :key="gene.id"
          >
            <polygon
              v-if="gene.start >= locus.start && gene.end >= gene.start && gene.end <= locus.end"
              :points="arrow(gene, locus)"
              :class="{ partial: gene.partial5 || gene.partial3 }"
            >
              <title>
                {{ gene.id }} · {{ gene.start }}–{{ gene.end }} · strand
                {{ gene.strand || "unknown" }}
              </title>
            </polygon>
          </template>
        </svg>
        <p v-if="!locus.genes.length">No CDSs are shown for this neighborhood.</p>
        <details v-if="!compact && locus.genes.length">
          <summary>Inspect {{ locus.genes.length }} CDSs and translations</summary>
          <div class="table-scroll">
            <table>
              <thead>
                <tr>
                  <th>Gene</th>
                  <th>Interval</th>
                  <th>Strand</th>
                  <th>Protein length</th>
                  <th>Genetic code</th>
                  <th>Partial call</th>
                </tr>
              </thead>
              <tbody>
                <tr
                  v-for="gene in locus.genes"
                  :key="gene.id"
                >
                  <td>{{ gene.id }}</td>
                  <td>{{ gene.start }}–{{ gene.end }}</td>
                  <td>{{ gene.strand || "Unknown" }}</td>
                  <td>{{ number(gene.proteinLength) }} aa</td>
                  <td>{{ number(gene.translationTable) }}</td>
                  <td>
                    {{
                      gene.partial5 === null || gene.partial3 === null
                        ? "Not reported"
                        : gene.partial5 || gene.partial3
                          ? "Partial"
                          : "No partial flag"
                    }}
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
        </details>
      </article>
      <p v-if="!links.length">
        {{
          truncated || counts.links !== 0
            ? "No similarity links in this preview."
            : "No protein alignments met the selected identity cutoff."
        }}
      </p>
      <div
        v-else
        class="table-scroll"
      >
        <table>
          <caption>
            Protein-similarity evidence
          </caption>
          <thead>
            <tr>
              <th>Query gene</th>
              <th>Target gene</th>
              <th>Identity</th>
              <th>Similarity</th>
              <th>Alignment score</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="link in shownLinks"
              :key="`${link.queryLocus}:${link.queryGene}:${link.targetLocus}:${link.targetGene}`"
            >
              <td>{{ link.queryLocus }} · {{ link.queryGene }}</td>
              <td>{{ link.targetLocus }} · {{ link.targetGene }}</td>
              <td>{{ percent(link.identity) }}</td>
              <td>{{ percent(link.similarity) }}</td>
              <td>{{ number(link.score) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p v-if="truncated || shownLinks.length < links.length">
        This preview is limited. Download the complete comparison and protein alignments for all
        evidence.
      </p>
      <p>
        Coordinates are 1-based and inclusive. Protein translations were checked against source DNA.
        Partial calls and source boundaries remain evidence limits.
      </p>
    </template>
  </section>
</template>

<style scoped>
.locus-comparison {
  display: grid;
  gap: 0.8rem;
}
.locus-comparison p {
  max-width: 90ch;
  margin: 0;
}
.locus-row {
  border-top: 1px solid var(--border-color, #ddd);
  padding-top: 0.75rem;
}
.locus-row h4 {
  overflow-wrap: anywhere;
}
svg {
  width: 100%;
  min-height: 44px;
}
polygon {
  fill: #245a75;
  stroke: #173c51;
}
polygon.partial {
  fill: #b98533;
}
.table-scroll {
  overflow-x: auto;
}
table {
  width: 100%;
  border-collapse: collapse;
}
th,
td {
  padding: 0.65rem;
  text-align: left;
  border-bottom: 1px solid var(--border-color, #ddd);
}
td {
  overflow-wrap: anywhere;
  max-width: 32ch;
}
caption {
  text-align: left;
  font-weight: 600;
}
</style>
