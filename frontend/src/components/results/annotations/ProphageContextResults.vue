<script setup>
import { computed } from "vue";
const props = defineProps({
  available: Boolean,
  regions: { type: Array, default: () => [] },
  matches: { type: Array, default: () => [] },
  sources: { type: Array, default: () => [] },
  sourceArrays: { type: Array, default: () => [] },
  sourceSpacers: { type: Array, default: () => [] },
  sourceInventoryComplete: Boolean,
  searchPolicy: { type: Object, default: () => ({}) },
  counts: { type: Object, default: () => ({}) },
  regionsTruncated: Boolean,
  matchesTruncated: Boolean,
  compact: Boolean,
});
const shownRegions = computed(() => props.regions.slice(0, props.compact ? 5 : 100));
const shownMatches = computed(() => props.matches.slice(0, props.compact ? 5 : 100));
const number = (value) => (value == null ? "Not reported" : value.toLocaleString());
const nativeArray = (id) =>
  props.sourceArrays.find((row) => row.id === id)?.nativeId || "Not in preview";
const searchStatus = (status) =>
  ({
    eligible: "Searched",
    deletion: "Deletion occurrence",
    category_excluded: "Category outside search policy",
    unsupported_length: "Length outside 16–1000 nt",
    ambiguous_query: "Ambiguous query bases",
  })[status] || "Not reported";
</script>

<template>
  <section
    class="prophage-context"
    aria-label="Prophage and spacer context"
  >
    <h3>Candidate prophages and source-spacer matches</h3>
    <p>
      Candidate regions and exact sequence matches provide genomic context. Neither region overlap
      nor a match demonstrates active CRISPR targeting.
    </p>
    <p v-if="!available">Validated prophage-context evidence is unavailable.</p>
    <template v-else>
      <p>
        {{ number(counts.regions) }} candidate regions · {{ number(counts.genes) }} called genes ·
        {{ number(counts.spacer_matches) }} spacer matches
      </p>
      <p>
        {{ number(counts.assessed_sources) }} of {{ number(counts.sources) }} source records
        assessed for prophages. Records shorter than 20,000 nt are not assessed by this pinned
        Phigaro workflow. Spacer matching remains a separate same-record analysis.
      </p>
      <p v-if="searchPolicy.complete === true">
        Completed exact full-length search of unambiguous {{ searchPolicy.minLength }}–{{
          searchPolicy.maxLength
        }}
        nt spacers from {{ searchPolicy.queryCategories.join(" / ") }} arrays, on both strands of
        the same source record. Matches overlapping any reported array are excluded. PAM
        compatibility and targeting activity are not assessed.
      </p>
      <p v-else>Source-spacer search policy or completion is not reported.</p>
      <p v-if="counts.assessed_sources === 0">
        No source records were assessed for prophages; this is not a negative prediction.
      </p>
      <details v-if="sources.length && !compact">
        <summary>Inspect source assessment scope</summary>
        <ul>
          <li
            v-for="source in sources"
            :key="source.id"
          >
            {{ source.id }} · {{ number(source.length) }} nt ·
            {{
              source.status === "assessed"
                ? "Assessed"
                : source.status === "ineligible_short_sequence"
                  ? "Not assessed: below 20,000 nt"
                  : "Assessment status unavailable"
            }}
          </li>
        </ul>
      </details>
      <p v-if="!regions.length">
        {{
          regionsTruncated || counts.regions !== 0 || !(counts.assessed_sources > 0)
            ? "No candidate regions in this preview."
            : "No prophage candidates were reported with the pinned models. This does not establish absence."
        }}
      </p>
      <div
        v-else
        class="table-scroll"
      >
        <table>
          <caption>
            Predicted prophage regions
          </caption>
          <thead>
            <tr>
              <th>Source</th>
              <th>Region</th>
              <th>Interval</th>
              <th>Supporting genes</th>
              <th>Peak score</th>
              <th>Record edge</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="region in shownRegions"
              :key="region.id"
            >
              <td>{{ region.sourceId }}</td>
              <td>{{ region.id }}</td>
              <td>{{ region.start }}–{{ region.end }}</td>
              <td>{{ region.supportingGeneIds.length }}</td>
              <td>{{ number(region.peakScore) }}</td>
              <td>
                {{
                  region.leftEdge || region.rightEdge
                    ? "First or last called gene"
                    : region.leftEdge === false && region.rightEdge === false
                      ? "Not at an edge gene"
                      : "Not reported"
                }}
              </td>
            </tr>
          </tbody>
        </table>
      </div>
      <p v-if="!matches.length">
        {{
          matchesTruncated || counts.spacer_matches !== 0
            ? "No spacer matches in this preview."
            : "No qualifying source-spacer matches were reported."
        }}
      </p>
      <div
        v-else
        class="table-scroll"
      >
        <table>
          <caption>
            Spacer sequence matches
          </caption>
          <thead>
            <tr>
              <th>Source array / spacer</th>
              <th>Target source</th>
              <th>Match interval</th>
              <th>Strand</th>
              <th>Overlapping candidate regions</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="(match, index) in shownMatches"
              :key="`${match.queryId}:${index}`"
            >
              <td>{{ match.sourceId }} · {{ match.arrayId }} · {{ number(match.spacerIndex) }}</td>
              <td>{{ match.targetSourceId }}</td>
              <td>{{ match.start }}–{{ match.end }}</td>
              <td>{{ match.strand || "Unknown" }}</td>
              <td>{{ match.regionIds.length ? match.regionIds.join(", ") : "None reported" }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p
        v-if="
          regionsTruncated ||
          matchesTruncated ||
          shownRegions.length < regions.length ||
          shownMatches.length < matches.length
        "
      >
        This preview is limited. Download the complete report for all regions and sequence matches.
      </p>
      <details v-if="sourceSpacers.length && !compact">
        <summary>Trace source spacer sequences and array identifiers</summary>
        <p v-if="!sourceInventoryComplete">
          The source inventory is limited in this preview. Download the complete report for all
          eligible, excluded and deleted spacer occurrences.
        </p>
        <div class="table-scroll">
          <table>
            <thead>
              <tr>
                <th>Query / source</th>
                <th>Array / native array</th>
                <th>Source interval</th>
                <th>Source-forward spacer</th>
                <th>Search status</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="spacer in sourceSpacers"
                :key="spacer.id"
              >
                <td>{{ spacer.id }} · {{ spacer.sourceId }}</td>
                <td>
                  {{ spacer.arrayId }} / {{ nativeArray(spacer.arrayId) }} · spacer
                  {{ number(spacer.ordinal) }}
                </td>
                <td>
                  {{
                    spacer.boundary !== null
                      ? `Deletion boundary ${spacer.boundary} (0-based interbase)`
                      : `${number(spacer.start)}–${number(spacer.end)}`
                  }}
                </td>
                <td>{{ spacer.sequence ?? "No observed sequence" }}</td>
                <td>{{ searchStatus(spacer.status) }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </details>
      <p>
        Coordinates are 1-based and inclusive. Region boundaries are computational predictions;
        proximity, overlap and biological activity are different evidence.
      </p>
    </template>
  </section>
</template>

<style scoped>
.prophage-context {
  display: grid;
  gap: 0.8rem;
}
.prophage-context p {
  max-width: 90ch;
  margin: 0;
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
  text-align: left;
  padding: 0.65rem;
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
