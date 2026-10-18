import { computed, onBeforeUnmount, onMounted, ref, watch } from "vue";
import { api } from "../api.js";
import { CATALOG_ENTITIES } from "../features/catalog/index.js";

export function useCatalog(client = api) {
  const summary = ref(null);
  const summaryLoading = ref(true);
  const summaryError = ref("");
  const entity = ref("assemblies");
  const searchType = ref("assemblies");
  const submitted = ref(false);
  const searchError = ref("");
  const query = ref("");
  const appliedQuery = ref("");
  const assembly = ref("");
  const records = ref([]);
  const nextCursor = ref(null);
  const cursors = ref([null]);
  const pageIndex = ref(0);
  const loading = ref(false);
  const error = ref("");
  const detail = ref(null);
  const detailLoading = ref(false);
  const detailError = ref("");
  const unitsLoading = ref("");
  let listController;
  let detailController;
  let summaryController;
  let listGeneration = 0;
  let requestedPage = 0;
  let detailGeneration = 0;
  const activeEntity = computed(() => CATALOG_ENTITIES.find((item) => item.id === entity.value));

  async function loadPage(index = 0) {
    if (!submitted.value) return;
    requestedPage = index;
    listController?.abort();
    listController = new AbortController();
    const generation = ++listGeneration;
    loading.value = true;
    error.value = "";
    try {
      const result = await client.catalogPage(
        entity.value,
        {
          q: appliedQuery.value,
          assembly_accession: assembly.value,
          limit: 25,
          cursor: cursors.value[index],
        },
        { signal: listController.signal },
      );
      if (generation !== listGeneration) return;
      records.value = result.items;
      nextCursor.value = result.next_cursor;
      pageIndex.value = index;
    } catch (failure) {
      if (generation === listGeneration && failure.name !== "AbortError") {
        error.value = failure.message || "Database records could not be loaded. Try again.";
      }
    } finally {
      if (generation === listGeneration) loading.value = false;
    }
  }

  function closeDetail() {
    detailController?.abort();
    detailGeneration += 1;
    detail.value = null;
    detailError.value = "";
    detailLoading.value = false;
    unitsLoading.value = "";
  }

  function resetPage() {
    listController?.abort();
    listGeneration += 1;
    loading.value = false;
    error.value = "";
    records.value = [];
    cursors.value = [null];
    nextCursor.value = null;
    pageIndex.value = 0;
    requestedPage = 0;
    closeDetail();
  }

  async function initialize() {
    summaryController?.abort();
    const controller = new AbortController();
    summaryController = controller;
    summaryLoading.value = true;
    summaryError.value = "";
    try {
      const result = await client.catalogSummary({ signal: controller.signal });
      if (controller.signal.aborted) return;
      summary.value = result;
    } catch (failure) {
      if (!controller.signal.aborted)
        summaryError.value = failure.message || "The database could not be reached.";
    } finally {
      if (!controller.signal.aborted) summaryLoading.value = false;
    }
  }

  function clearSearch() {
    submitted.value = false;
    appliedQuery.value = "";
    assembly.value = "";
    searchError.value = "";
    resetPage();
  }

  // Changing the input cancels previous requests; results always belong to a
  // deliberate submission, including clients that ignore AbortSignal.
  watch(query, clearSearch, { flush: "sync" });
  watch(
    searchType,
    () => {
      query.value = "";
      clearSearch();
    },
    { flush: "sync" },
  );

  function chooseEntity(value) {
    if (
      !submitted.value ||
      !assembly.value ||
      !CATALOG_ENTITIES.some((item) => item.id === value) ||
      value === entity.value
    )
      return;
    entity.value = value;
    resetPage();
    void loadPage();
  }

  function search() {
    clearSearch();
    let value = query.value.trim();
    if (!value) {
      searchError.value = "Enter an accession, sequence or Cas gene to search.";
      return;
    }
    if (searchType.value === "assemblies") {
      value = value.toUpperCase();
      if (!/^GC[AF]_\d+\.\d+$/u.test(value)) {
        searchError.value =
          "Enter a complete genome accession, including its version (GCA_… or GCF_…).";
        return;
      }
      assembly.value = value;
    } else if (["repeats", "spacers"].includes(searchType.value)) {
      value = value.replace(/\s/gu, "").toUpperCase();
      if (!/^[ACGTRYSWKMBDHVN]+$/u.test(value)) {
        searchError.value =
          "Enter a nucleotide sequence using A, C, G, T or IUPAC ambiguity codes.";
        return;
      }
    }
    if (value.length > 1000) {
      searchError.value = "Enter a search value of 1,000 characters or fewer.";
      return;
    }
    entity.value = searchType.value;
    appliedQuery.value = assembly.value ? "" : value;
    submitted.value = true;
    void loadPage();
  }

  function chooseAssembly(accession) {
    searchType.value = "assemblies";
    query.value = accession;
    clearSearch();
    assembly.value = accession;
    entity.value = "arrays";
    submitted.value = true;
    void loadPage();
  }

  function newSearch() {
    query.value = "";
    clearSearch();
  }

  function nextPage() {
    if (!nextCursor.value || loading.value) return;
    cursors.value[pageIndex.value + 1] = nextCursor.value;
    void loadPage(pageIndex.value + 1);
  }

  async function openArray(id) {
    closeDetail();
    detailController = new AbortController();
    const generation = detailGeneration;
    const options = { signal: detailController.signal };
    detailLoading.value = true;
    try {
      const [array, repeats, spacers] = await Promise.all([
        client.catalogArray(id, options),
        client.catalogPage("repeats", { array_id: id, limit: 100 }, options),
        client.catalogPage("spacers", { array_id: id, limit: 100 }, options),
      ]);
      if (generation === detailGeneration) detail.value = { array, repeats, spacers };
    } catch (failure) {
      if (generation === detailGeneration && failure.name !== "AbortError") {
        detailError.value =
          failure.message || "The array could not be loaded. Open it again to retry.";
      }
    } finally {
      if (generation === detailGeneration) detailLoading.value = false;
    }
  }

  async function moreUnits(kind) {
    if (!detail.value?.[kind]?.next_cursor || unitsLoading.value) return;
    const generation = detailGeneration;
    unitsLoading.value = kind;
    detailError.value = "";
    try {
      const result = await client.catalogPage(
        kind,
        {
          array_id: detail.value.array.id,
          limit: 100,
          cursor: detail.value[kind].next_cursor,
        },
        { signal: detailController.signal },
      );
      if (generation === detailGeneration) detail.value[kind] = result;
    } catch (failure) {
      if (generation === detailGeneration && failure.name !== "AbortError") {
        detailError.value = failure.message || "The next occurrences could not be loaded.";
      }
    } finally {
      if (generation === detailGeneration) unitsLoading.value = "";
    }
  }

  onMounted(initialize);
  onBeforeUnmount(() => {
    summaryController?.abort();
    listController?.abort();
    detailController?.abort();
    listGeneration += 1;
    detailGeneration += 1;
  });
  return {
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
    previousPage: () => loadPage(Math.max(0, pageIndex.value - 1)),
    retryPage: () => loadPage(requestedPage),
    openArray,
    closeDetail,
    moreUnits,
  };
}
