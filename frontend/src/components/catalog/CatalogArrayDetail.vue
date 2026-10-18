<script setup>
import { computed, nextTick, onMounted, ref } from "vue";
import {
  catalogCell,
  catalogInterval,
  catalogValue,
  sequenceSearchIssue,
} from "../../features/catalog/index.js";
import SequenceSearchPanel from "./SequenceSearchPanel.vue";

const props = defineProps({
  array: { type: Object, required: true },
  repeats: { type: Object, required: true },
  spacers: { type: Object, required: true },
  loading: { type: String, default: "" },
});
defineEmits(["close", "more"]);
const heading = ref(null);
const selected = ref(null);
let searchTrigger;
function search(section, row, event) {
  searchTrigger = event.currentTarget;
  selected.value = {
    kind: section.key,
    sequence: row.sequence,
    label: `${section.key === "repeats" ? "Repeat" : "Spacer"} ${row.ordinal} in array ${props.array.call_id}`,
  };
}
async function closeSearch() {
  selected.value = null;
  await nextTick();
  searchTrigger?.focus();
}
onMounted(() => {
  heading.value?.focus({ preventScroll: true });
  heading.value?.scrollIntoView({ block: "start" });
});
const units = computed(() =>
  [
    ...props.repeats.items.map((row) => ({ ...row, unit: "Repeat", position: row.start })),
    ...props.spacers.items.map((row) => ({
      ...row,
      unit: "Spacer",
      position: row.kind === "deletion" ? row.boundary_position : row.start,
    })),
  ]
    .filter((row) => row.position !== null && row.position !== undefined)
    .sort((left, right) => left.position - right.position)
    .slice(0, 40),
);
const sections = computed(() => [
  { key: "repeats", label: "Repeat occurrences", page: props.repeats },
  { key: "spacers", label: "Spacer occurrences", page: props.spacers },
]);
</script>

<template>
  <section
    class="catalog-array-detail"
    aria-label="Array details"
  >
    <div class="catalog-detail-heading">
      <div>
        <p>{{ array.accession }} / {{ catalogValue(array.sequence_id) }}</p>
        <h2
          ref="heading"
          tabindex="-1"
        >
          Array {{ array.call_id }}
        </h2>
      </div>
      <button
        class="catalog-button"
        type="button"
        @click="$emit('close')"
      >
        Close array
      </button>
    </div>
    <dl class="catalog-array-facts">
      <div>
        <dt>Coordinates</dt>
        <dd>{{ catalogInterval(array) }}</dd>
      </div>
      <div>
        <dt>Strand</dt>
        <dd>{{ catalogCell(array, "strand") }}</dd>
      </div>
      <div>
        <dt>Caller</dt>
        <dd>{{ catalogValue(array.tool) }} {{ array.tool_version }}</dd>
      </div>
      <div>
        <dt>Category</dt>
        <dd>{{ catalogValue(array.category) }}</dd>
      </div>
      <div>
        <dt>Caller score</dt>
        <dd>{{ catalogValue(array.certainty_score) }}</dd>
      </div>
      <div>
        <dt>Evidence level</dt>
        <dd>{{ catalogValue(array.evidence_level) }}</dd>
      </div>
    </dl>
    <p class="catalog-detail-note">
      Caller scores are not calibrated probabilities. Coordinates refer to the forward genomic
      source; unknown strand is not inferred.
    </p>
    <p class="catalog-consensus">
      <strong>Repeat consensus</strong><code>{{ catalogValue(array.repeat_consensus) }}</code>
    </p>
    <div
      v-if="units.length"
      class="catalog-unit-map"
      aria-label="Repeat and spacer occurrence map"
    >
      <div class="catalog-unit-strip">
        <span
          v-for="unit in units"
          :key="`${unit.unit}-${unit.id}`"
          :class="[
            'catalog-unit',
            unit.unit === 'Repeat' ? 'catalog-unit-repeat' : 'catalog-unit-spacer',
            { 'catalog-unit-deletion': unit.kind === 'deletion' },
          ]"
          :title="`${unit.unit} ${unit.ordinal}: ${catalogInterval(unit)}`"
          >{{ unit.kind === "deletion" ? "Δ" : unit.unit[0] }}{{ unit.ordinal }}</span
        >
      </div>
      <p>
        R = repeat; S = spacer; Δ = deletion. Up to 40 loaded occurrences in source-coordinate
        order. Widths are not to genomic scale.
      </p>
    </div>
    <div class="catalog-occurrence-sections">
      <section
        v-for="section in sections"
        :key="section.key"
        :aria-label="section.label"
      >
        <h3>{{ section.label }}</h3>
        <div
          v-if="section.page.items.length"
          class="catalog-table-scroll"
          tabindex="0"
          :aria-label="`${section.label} table`"
          role="region"
        >
          <table class="catalog-table catalog-occurrence-table sequence-responsive-table">
            <thead>
              <tr>
                <th scope="col">Position</th>
                <th scope="col">Coordinates</th>
                <th scope="col">Sequence</th>
                <th scope="col">Database</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="row in section.page.items"
                :key="row.id"
              >
                <td data-label="Position">{{ row.ordinal }}</td>
                <td data-label="Coordinates">{{ catalogInterval(row) }}</td>
                <td data-label="Sequence">
                  <code class="catalog-sequence">{{ catalogCell(row, "sequence") }}</code>
                </td>
                <td data-label="Database">
                  <button
                    type="button"
                    class="catalog-button"
                    :disabled="
                      Boolean(
                        sequenceSearchIssue(
                          row.sequence,
                          row.kind === 'deletion' ? 'deletion' : null,
                        ),
                      )
                    "
                    :aria-label="`Search in DB: ${section.key === 'repeats' ? 'repeat' : 'spacer'} ${row.ordinal}`"
                    @click="search(section, row, $event)"
                  >
                    Search in DB
                  </button>
                  <small
                    v-if="
                      sequenceSearchIssue(row.sequence, row.kind === 'deletion' ? 'deletion' : null)
                    "
                  >
                    {{
                      sequenceSearchIssue(row.sequence, row.kind === "deletion" ? "deletion" : null)
                    }}
                  </small>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
        <p v-else>No {{ section.key }} reported for this array.</p>
        <button
          v-if="section.page.next_cursor"
          type="button"
          class="catalog-button"
          :disabled="Boolean(loading)"
          @click="$emit('more', section.key)"
        >
          {{ loading === section.key ? "Loading…" : `Next ${section.key}` }}
        </button>
        <p class="catalog-detail-note">
          {{ section.page.items.length }} occurrences on this page{{
            section.page.next_cursor ? "; more available." : "."
          }}
        </p>
      </section>
    </div>
    <SequenceSearchPanel
      v-if="selected"
      :key="selected.label"
      v-bind="selected"
      @close="closeSearch"
    />
  </section>
</template>
