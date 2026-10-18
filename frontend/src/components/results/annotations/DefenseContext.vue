<script setup>
import { computed } from "vue";

const props = defineProps({
  available: { type: Boolean, default: false },
  systems: { type: Array, default: () => [] },
  components: { type: Array, default: () => [] },
  hits: { type: Array, default: () => [] },
  genes: { type: Array, default: () => [] },
  counts: { type: Object, default: () => ({}) },
  systemsTruncated: { type: Boolean, default: false },
  componentsTruncated: { type: Boolean, default: false },
  hitsTruncated: { type: Boolean, default: false },
});
const geneMap = computed(() => new Map(props.genes.map((gene) => [gene.id, gene])));
const bound = (value) => (Number.isInteger(value) && value >= 0 ? value : "Unknown");
const interval = (gene) =>
  gene && Number.isInteger(gene.start) && Number.isInteger(gene.end)
    ? `${gene.start.toLocaleString()}–${gene.end.toLocaleString()}`
    : "Unknown";
const score = (value) =>
  typeof value === "number" && Number.isFinite(value) ? value.toPrecision(3) : "Unknown";
</script>

<template>
  <section
    class="defense-context"
    aria-label="Defense context"
  >
    <h3>Predicted defense repertoire</h3>
    <p>
      DefenseFinder model release 1.2.3 assesses the submitted sequences. Newer defense models and
      CasFinder are not assessed. CRISPR classifications remain separate.
    </p>
    <p>
      System calls satisfy model rules; they do not establish activity or cooperative immunity. A
      fragment cannot establish genome-wide absence.
    </p>
    <p v-if="!available">Defense evidence is unavailable.</p>
    <template v-else>
      <p>
        {{ bound(counts.systems) }} system calls · {{ bound(counts.components) }} component hits ·
        {{ bound(counts.hits) }} profile hits · {{ bound(counts.genes) }} assessed proteins
      </p>
      <p v-if="systemsTruncated">
        System preview is limited. Download the complete evidence for all calls.
      </p>
      <p v-if="!systems.length">
        {{
          systemsTruncated || counts.systems !== 0
            ? "No system calls in this preview."
            : "No modeled systems detected in the assessed sequence."
        }}
      </p>
      <div
        v-else
        class="defense-table-wrap"
      >
        <table>
          <caption>
            System-rule predictions
          </caption>
          <thead>
            <tr>
              <th>Source</th>
              <th>System</th>
              <th>Interval (1-based)</th>
              <th>Components</th>
              <th>Model completeness</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="system in systems"
              :key="system.id"
            >
              <td>{{ system.source_id }}</td>
              <td>{{ system.subtype }}</td>
              <td>{{ interval(system) }}</td>
              <td>{{ Array.isArray(system.gene_ids) ? system.gene_ids.length : "Unknown" }}</td>
              <td>{{ score(system.wholeness) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <details v-if="components.length || componentsTruncated">
        <summary>Inspect component assignments</summary>
        <p v-if="componentsTruncated">
          Component preview is limited; the download contains all component assignments.
        </p>
        <ul>
          <li
            v-for="(item, index) in components"
            :key="`${item.system_id}:${item.gene_id}:${index}`"
          >
            {{ item.gene_id }} · {{ item.profile }} · {{ item.status }} ·
            {{ interval(geneMap.get(item.gene_id)) }}
          </li>
        </ul>
      </details>
      <details v-if="hits.length || hitsTruncated">
        <summary>Inspect profile hits, including hits outside system calls</summary>
        <p>Profile matches are separate evidence and are not automatically complete systems.</p>
        <p v-if="hitsTruncated">
          Profile-hit preview is limited; the download contains the complete table.
        </p>
        <ul>
          <li
            v-for="(hit, index) in hits"
            :key="`${hit.gene_id}:${hit.profile}:${index}`"
          >
            {{ hit.gene_id }} · {{ hit.profile }} · score {{ score(hit.score) }} ·
            {{ interval(geneMap.get(hit.gene_id)) }}
          </li>
        </ul>
      </details>
    </template>
  </section>
</template>

<style scoped>
.defense-context {
  display: grid;
  gap: 0.8rem;
}
.defense-context p {
  max-width: 85ch;
  margin: 0;
}
.defense-table-wrap {
  overflow-x: auto;
}
table {
  border-collapse: collapse;
  width: 100%;
}
th,
td {
  padding: 0.7rem;
  text-align: left;
  border-bottom: 1px solid var(--border-color, #ddd);
}
caption {
  text-align: left;
  font-weight: 600;
  margin-bottom: 0.5rem;
}
li {
  overflow-wrap: anywhere;
  margin: 0.4rem 0;
}
</style>
