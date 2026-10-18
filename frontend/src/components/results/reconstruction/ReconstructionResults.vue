<script setup>
import { computed } from "vue";

import { asArray, formatDuration, formatNumber } from "../../../utils/formatting.js";
import AppIcon from "../../common/AppIcon.vue";
import ReconstructionStory from "./ReconstructionStory.vue";

const props = defineProps({ summary: { type: Object, required: true } });
const orientation = computed(() => props.summary.orientation);
const reconstruction = computed(() => props.summary.reconstruction);
const rows = computed(() =>
  asArray(orientation.value?.selected_reconstructions).length
    ? asArray(orientation.value.selected_reconstructions)
    : asArray(reconstruction.value?.results).length
      ? asArray(reconstruction.value.results)
      : reconstruction.value?.selected_model
        ? [reconstruction.value.selected_model]
        : [],
);
const treePolicy = computed(
  () =>
    orientation.value?.tree_policy ||
    reconstruction.value?.tree_policy ||
    props.summary?.pipeline?.stages?.spacerplacer?.tree_source ||
    (reconstruction.value ? "estimated_from_arrays_uncalibrated_scale" : "not_reported"),
);
const deletionCount = (row) => row.deletions;
const noDeletionGroups = computed(() =>
  rows.value
    .filter((row) => typeof deletionCount(row) === "number" && deletionCount(row) === 0)
    .map((row) => row.name),
);
</script>

<template>
  <section
    v-if="rows.length"
    class="result-section reconstruction-section"
    aria-labelledby="reconstruction-heading"
  >
    <div class="result-heading">
      <div>
        <h3 id="reconstruction-heading">Spacer history</h3>
      </div>
      <p>SpacerPlacer reconstruction</p>
    </div>
    <div
      v-if="orientation?.reconstructions_truncated"
      class="history-truncation"
      role="note"
    >
      <AppIcon
        name="warning"
        :size="18"
      />
      <p>
        <strong>Some structured histories were omitted from this summary.</strong> Use result
        artifacts for complete detail.
      </p>
    </div>
    <div class="reconstruction-story-list">
      <ReconstructionStory
        v-for="(row, index) in rows"
        :key="row.name"
        :row="row"
        :index="index"
        :summary="summary"
      />
    </div>
    <details class="reconstruction-values">
      <summary>Exact SpacerPlacer estimates and runtime</summary>
      <div class="tree-policy">
        <span
          class="tree-glyph"
          aria-hidden="true"
          ><AppIcon name="tree"
        /></span>
        <div>
          <strong>Tree policy used: {{ String(treePolicy).replaceAll("_", " ") }}</strong>
          <p v-if="orientation?.comparisons?.length">
            Decisive groups report the supported hypothesis. Unresolved groups retain input order
            only as a reporting default.
          </p>
          <p v-else>
            Input orientation; not compared. The tree is inferred from the supplied arrays, with an
            uncalibrated branch scale.
          </p>
        </div>
      </div>

      <div class="table-wrap reconstruction-table">
        <table>
          <thead>
            <tr>
              <th>Group</th>
              <th>Preferred model</th>
              <th>BDM lnL</th>
              <th>Insertions</th>
              <th>Deletions</th>
              <th>BDM deletion rate</th>
              <th>Runtime</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="(row, index) in rows"
              :key="row.name || index"
            >
              <td>
                <strong>{{ row.name || `Group ${index + 1}` }}</strong>
              </td>
              <td>{{ row.preferred_model || "—" }}</td>
              <td>{{ formatNumber(row.bdm_log_likelihood, 3) }}</td>
              <td>{{ formatNumber(row.acquisitions) }}</td>
              <td>{{ formatNumber(deletionCount(row)) }}</td>
              <td>{{ formatNumber(row.deletion_rate_bdm, 4) }}</td>
              <td>{{ formatDuration(row.runtime_seconds) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
    </details>
    <div
      v-if="noDeletionGroups.length"
      class="warning-note"
    >
      <AppIcon name="warning" />
      <p>
        <strong
          >No deletion events were reconstructed for {{ noDeletionGroups.join(", ") }}.</strong
        >
        Deletion-rate estimates and orientation evidence may not be meaningful for those groups.
      </p>
    </div>
  </section>
  <p
    v-else
    class="empty-result"
  >
    No reconstruction was reported.
  </p>
</template>
