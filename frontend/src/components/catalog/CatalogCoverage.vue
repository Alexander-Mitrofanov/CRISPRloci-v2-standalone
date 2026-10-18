<script setup>
import { computed } from "vue";
import { catalogNumber, catalogValue } from "../../features/catalog/index.js";

const props = defineProps({ summary: { type: Object, required: true } });
const coverage = computed(() => props.summary.coverage || {});
const fullCoverage = computed(
  () =>
    coverage.value.imported_assemblies != null &&
    coverage.value.imported_assemblies === coverage.value.catalog_assemblies,
);
const sourceLink = computed(() => {
  try {
    const url = new URL(props.summary.source_url);
    return url.protocol === "https:" && !url.username && !url.password ? url.href : null;
  } catch {
    return null;
  }
});
</script>

<template>
  <aside
    class="catalog-coverage"
    aria-label="Dataset coverage and provenance"
  >
    <div class="catalog-coverage-heading">
      <h2>{{ fullCoverage ? "Complete catalog imported" : "Imported dataset coverage" }}</h2>
      <span class="catalog-snapshot">{{ catalogValue(summary.snapshot_id) }}</span>
    </div>
    <p class="catalog-coverage-scope">
      <strong>{{ catalogNumber(coverage.imported_assemblies) }}</strong> genome records available
      here
      <span v-if="coverage.catalog_assemblies != null"
        >from a campaign catalog of
        {{ catalogNumber(coverage.catalog_assemblies) }} assemblies.</span
      >
    </p>
    <p
      v-if="!fullCoverage"
      class="catalog-coverage-note"
    >
      Searches and counts below cover this imported snapshot. An absent record does not establish
      absence in the full campaign.
    </p>
    <details class="catalog-provenance">
      <summary>Dataset source, scope and scientific limits</summary>
      <dl>
        <div>
          <dt>Source</dt>
          <dd>
            <a
              v-if="sourceLink"
              :href="sourceLink"
              target="_blank"
              rel="noopener noreferrer"
              >{{ summary.title || "CRISPR-Cas database" }}</a
            >
            <span v-else>{{ summary.title || "CRISPR-Cas database" }}</span>
          </dd>
        </div>
        <div>
          <dt>Snapshot created</dt>
          <dd>{{ catalogValue(summary.generated_at) }}</dd>
        </div>
        <div>
          <dt>Campaign successful analyses</dt>
          <dd>{{ catalogNumber(coverage.successful_assemblies) }}</dd>
        </div>
        <div>
          <dt>Campaign failed analyses</dt>
          <dd>{{ catalogNumber(coverage.failed_assemblies) }}</dd>
        </div>
        <div>
          <dt>Campaign unavailable assemblies</dt>
          <dd>{{ catalogNumber(coverage.unavailable_assemblies) }}</dd>
        </div>
        <div>
          <dt>Coordinates</dt>
          <dd>
            {{
              summary.coordinate_system ||
              "1-based inclusive, forward genomic source; deletion boundaries are 0-based interbase positions."
            }}
          </dd>
        </div>
        <div
          v-for="(value, key) in summary.provenance"
          :key="key"
        >
          <dt>{{ key.replaceAll("_", " ") }}</dt>
          <dd>{{ catalogValue(value) }}</dd>
        </div>
      </dl>
      <ul v-if="summary.limitations?.length">
        <li
          v-for="limitation in summary.limitations"
          :key="limitation"
        >
          {{ limitation }}
        </li>
      </ul>
    </details>
  </aside>
</template>
