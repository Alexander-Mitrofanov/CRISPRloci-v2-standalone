<script setup>
import { computed, ref } from "vue";
import { useCatalog } from "../../composables/useCatalog.js";
import { CATALOG_ENTITIES, catalogCsv, catalogFasta } from "../../features/catalog/index.js";
import { saveBlob } from "../../utils/download.js";
import CatalogCoverage from "./CatalogCoverage.vue";
import CatalogTable from "./CatalogTable.vue";
import CatalogArrayDetail from "./CatalogArrayDetail.vue";

const {
  summary,
  summaryLoading,
  summaryError,
  entity,
  activeEntity,
  searchType,
  submitted,
  searchError,
  query,
  appliedQuery,
  assembly,
  records,
  nextCursor,
  pageIndex,
  loading,
  error,
  detail,
  detailLoading,
  detailError,
  unitsLoading,
  initialize,
  chooseEntity,
  search,
  chooseAssembly,
  newSearch,
  nextPage,
  previousPage,
  retryPage,
  openArray: loadArray,
  closeDetail: clearDetail,
  moreUnits,
} = useCatalog();
const queryInput = ref(null);
const searchFields = {
  assemblies: {
    label: "Accession number",
    placeholder: "GCA_… or GCF_…",
    hint: "Enter the complete genome accession, including its version.",
  },
  repeats: {
    label: "Repeat sequence",
    placeholder: "Enter a nucleotide sequence",
    hint: "Exact sequence matches. IUPAC ambiguity codes are preserved.",
  },
  spacers: {
    label: "Spacer sequence",
    placeholder: "Enter a nucleotide sequence",
    hint: "Exact sequence matches. IUPAC ambiguity codes are preserved.",
  },
  cas_genes: {
    label: "Cas gene or profile",
    placeholder: "Enter a Cas gene or profile",
    hint: "Search a Cas gene name or the original annotation or profile.",
  },
};
const searchField = computed(() => searchFields[searchType.value]);
const fasta = computed(() => catalogFasta(records.value, entity.value));
let arrayTrigger;
function openArray(id) {
  arrayTrigger = document.activeElement;
  void loadArray(id);
}
function closeDetail() {
  clearDetail();
  if (arrayTrigger?.isConnected) {
    arrayTrigger.focus({ preventScroll: true });
    arrayTrigger.scrollIntoView({ block: "nearest" });
  }
}
function resetSearch() {
  newSearch();
  queryInput.value?.focus();
}
function exportPage(format) {
  const columns = [...new Set(records.value.flatMap((row) => Object.keys(row)))].map((key) => [
    key,
    key,
  ]);
  const text = format === "csv" ? catalogCsv(records.value, columns) : fasta.value;
  saveBlob(
    new Blob([text], {
      type: format === "csv" ? "text/csv;charset=utf-8" : "text/plain;charset=utf-8",
    }),
    `crisprloci-${entity.value}-page-${pageIndex.value + 1}.${format === "csv" ? "csv" : "fasta"}`,
  );
}
function tabKey(event, index) {
  let next = index;
  if (event.key === "ArrowRight") next = (index + 1) % CATALOG_ENTITIES.length;
  else if (event.key === "ArrowLeft")
    next = (index + CATALOG_ENTITIES.length - 1) % CATALOG_ENTITIES.length;
  else if (event.key === "Home") next = 0;
  else if (event.key === "End") next = CATALOG_ENTITIES.length - 1;
  else return;
  event.preventDefault();
  chooseEntity(CATALOG_ENTITIES[next].id);
  document.getElementById(`catalog-tab-${CATALOG_ENTITIES[next].id}`)?.focus();
}
</script>

<template>
  <section
    class="catalog-workspace"
    :class="{ 'catalog-workspace-empty': !submitted }"
    aria-labelledby="catalog-heading"
  >
    <div class="catalog-heading">
      <h1
        id="catalog-heading"
        tabindex="-1"
      >
        CRISPR–Cas database
      </h1>
      <p>Find arrays and Cas annotations by accession, repeat, spacer or Cas gene.</p>
    </div>
    <p
      v-if="summaryLoading && !summary"
      role="status"
      class="catalog-search-hint"
    >
      Connecting to the database…
    </p>
    <div
      v-else-if="summaryError"
      class="catalog-state catalog-error"
      role="alert"
    >
      <p>{{ summaryError }}</p>
      <button
        class="catalog-button"
        type="button"
        @click="initialize"
      >
        Retry database
      </button>
    </div>
    <div
      v-else-if="summary && !summary.available"
      class="catalog-state"
    >
      <h2>Database is not available yet</h2>
      <p>No dataset is currently available for searching.</p>
      <button
        class="catalog-button"
        type="button"
        @click="initialize"
      >
        Check availability
      </button>
    </div>
    <template v-else-if="summary?.available">
      <form
        class="catalog-search"
        aria-label="Database search"
        @submit.prevent="search"
      >
        <label class="catalog-search-type"
          >Search by
          <select
            v-model="searchType"
            name="search-type"
          >
            <option value="assemblies">Accession number</option>
            <option value="repeats">Repeat</option>
            <option value="spacers">Spacer</option>
            <option value="cas_genes">Cas gene</option>
          </select>
        </label>
        <label class="catalog-search-query"
          >{{ searchField.label }}
          <input
            ref="queryInput"
            v-model="query"
            name="query"
            type="search"
            autocomplete="off"
            spellcheck="false"
            :placeholder="searchField.placeholder"
            aria-describedby="catalog-search-hint"
            :aria-invalid="searchError ? 'true' : undefined"
          />
        </label>
        <button
          class="catalog-button catalog-button-primary"
          type="submit"
          :disabled="loading || !query.trim()"
        >
          Search
        </button>
      </form>
      <p
        id="catalog-search-hint"
        class="catalog-search-hint"
      >
        {{ searchField.hint }}
      </p>
      <p
        v-if="searchError"
        class="catalog-search-validation"
        role="alert"
      >
        {{ searchError }}
      </p>
      <section
        v-if="submitted"
        class="catalog-results"
        aria-label="Search results"
        :aria-busy="loading"
      >
        <div class="catalog-result-heading">
          <div>
            <h2>{{ assembly ? "Genome results" : activeEntity.label }}</h2>
            <p>{{ assembly || appliedQuery }}</p>
          </div>
          <button
            class="catalog-text-button"
            type="button"
            @click="resetSearch"
          >
            New search
          </button>
        </div>
        <div
          v-if="assembly"
          class="catalog-tabs"
          role="tablist"
          aria-label="Genome record type"
        >
          <button
            v-for="(item, index) in CATALOG_ENTITIES"
            :id="`catalog-tab-${item.id}`"
            :key="item.id"
            role="tab"
            type="button"
            :aria-selected="entity === item.id"
            aria-controls="catalog-records"
            :tabindex="entity === item.id ? 0 : -1"
            @click="chooseEntity(item.id)"
            @keydown="tabKey($event, index)"
          >
            {{ item.label }}
          </button>
        </div>
        <div
          id="catalog-records"
          :role="assembly ? 'tabpanel' : undefined"
          :aria-labelledby="assembly ? `catalog-tab-${entity}` : undefined"
        >
          <p
            v-if="loading"
            role="status"
            class="catalog-search-hint"
          >
            Searching…
          </p>
          <div
            v-else-if="error"
            class="catalog-state catalog-error"
            role="alert"
          >
            <p>{{ error }}</p>
            <button
              class="catalog-button"
              type="button"
              @click="retryPage"
            >
              Retry search
            </button>
          </div>
          <template v-else>
            <template v-if="records.length">
              <div class="catalog-results-toolbar">
                <p>
                  {{ records.length }} {{ records.length === 1 ? "record" : "records" }} on this
                  page
                </p>
                <div class="catalog-exports">
                  <button
                    type="button"
                    @click="exportPage('csv')"
                  >
                    Download this page as CSV
                  </button>
                  <button
                    v-if="fasta"
                    type="button"
                    @click="exportPage('fasta')"
                  >
                    {{
                      entity === "arrays"
                        ? "Download consensus FASTA"
                        : "Download this page as FASTA"
                    }}
                  </button>
                </div>
              </div>
              <CatalogTable
                :rows="records"
                :columns="activeEntity.columns"
                :label="`${activeEntity.label} records`"
                @array="openArray"
                @genome="chooseAssembly"
              />
              <p class="catalog-entity-description">{{ activeEntity.description }}</p>
            </template>
            <div
              v-else
              class="catalog-empty-result"
              role="status"
            >
              <p>
                {{
                  nextCursor
                    ? "No exact matches on this page. Continue to the next page."
                    : "No matching records found."
                }}
              </p>
            </div>
            <nav
              v-if="pageIndex > 0 || nextCursor"
              class="catalog-pagination"
              aria-label="Database pagination"
            >
              <button
                class="catalog-button"
                type="button"
                :disabled="pageIndex === 0"
                @click="previousPage"
              >
                Previous page
              </button>
              <span>Page {{ pageIndex + 1 }}</span>
              <button
                class="catalog-button"
                type="button"
                :disabled="!nextCursor"
                @click="nextPage"
              >
                Next page
              </button>
            </nav>
          </template>
        </div>
        <p
          v-if="detailLoading"
          role="status"
          class="catalog-search-hint"
        >
          Loading array and occurrences…
        </p>
        <p
          v-if="detailError"
          role="alert"
          class="catalog-state catalog-error"
        >
          {{ detailError }}
        </p>
        <CatalogArrayDetail
          v-if="detail"
          :key="detail.array.id"
          :array="detail.array"
          :repeats="detail.repeats"
          :spacers="detail.spacers"
          :loading="unitsLoading"
          @close="closeDetail"
          @more="moreUnits"
        />
        <details class="catalog-about">
          <summary>About this dataset</summary>
          <CatalogCoverage :summary="summary" />
        </details>
      </section>
    </template>
  </section>
</template>
