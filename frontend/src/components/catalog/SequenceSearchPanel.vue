<script setup>
import { computed, onMounted, ref, useId } from "vue";
import { catalogInterval, useSequenceSearch } from "../../features/catalog/index.js";

const props = defineProps({
  kind: { type: String, required: true },
  sequence: { type: String, required: true },
  label: { type: String, required: true },
});
defineEmits(["close"]);
const heading = ref(null);
const headingId = useId();
const { items, loading, error, nextCursor, page, retry, previous, next } = useSequenceSearch(
  computed(() => ({ kind: props.kind, sequence: props.sequence })),
);
onMounted(() => heading.value?.focus());
</script>

<template>
  <section
    class="sequence-search-panel"
    :aria-labelledby="headingId"
  >
    <div class="sequence-search-heading">
      <h4
        :id="headingId"
        ref="heading"
        tabindex="-1"
      >
        Similar candidates: {{ label }}
      </h4>
      <button
        type="button"
        class="catalog-button"
        @click="$emit('close')"
      >
        Close search
      </button>
    </div>
    <code class="sequence-search-query">{{ sequence }}</code>
    <p class="sequence-search-note">
      Full-length matches with up to one nucleotide substitution on either strand. No gaps or
      partial matches; IUPAC ambiguity symbols are compared literally. Identity describes sequence
      agreement, not biological function. Candidates are in database order.
    </p>
    <p
      v-if="loading"
      role="status"
    >
      Searching the database…
    </p>
    <div
      v-if="error"
      role="alert"
    >
      <p>{{ error }}</p>
      <button
        type="button"
        class="catalog-button"
        :disabled="loading"
        @click="retry"
      >
        Retry search
      </button>
    </div>
    <template v-if="!loading && !error">
      <div
        v-if="items.length"
        class="table-wrap"
        tabindex="0"
        role="region"
        aria-label="Similar sequence candidates"
      >
        <table class="sequence-responsive-table">
          <thead>
            <tr>
              <th>Genome / contig</th>
              <th>Array / coordinates</th>
              <th>Matched sequence</th>
              <th>Sequence comparison</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="row in items"
              :key="row.id"
            >
              <td data-label="Genome / contig">
                <strong>{{ row.accession }}</strong
                ><small>{{ row.sequence_id }}</small>
              </td>
              <td data-label="Array / coordinates">
                {{ row.array_id }}<small>{{ catalogInterval(row) }}</small>
              </td>
              <td data-label="Matched sequence">
                <code class="sequence-search-candidate">{{ row.sequence }}</code>
              </td>
              <td data-label="Sequence comparison">
                {{ row.comparison.identity }}% identity<small
                  >{{ row.comparison.substitutions }}
                  {{ row.comparison.substitutions === 1 ? "substitution" : "substitutions" }} ·
                  {{ row.comparison.strand }}</small
                >
              </td>
            </tr>
          </tbody>
        </table>
      </div>
      <p
        v-else
        role="status"
      >
        {{
          nextCursor
            ? "No candidates on this page. Continue to check the next database records."
            : "No candidates on this page within one substitution on either strand."
        }}
      </p>
      <nav
        class="sequence-search-pagination"
        aria-label="Sequence candidate pages"
      >
        <button
          type="button"
          class="catalog-button"
          :disabled="page === 0"
          @click="previous"
        >
          Previous candidates
        </button>
        <span role="status"
          >Page {{ page + 1 }} · {{ items.length }}
          {{ items.length === 1 ? "candidate" : "candidates" }}</span
        >
        <button
          type="button"
          class="catalog-button"
          :disabled="!nextCursor"
          @click="next"
        >
          Next candidates
        </button>
      </nav>
    </template>
  </section>
</template>
