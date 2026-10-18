import { onBeforeUnmount, ref, watch } from "vue";
import { api } from "../../api.js";
import { sequenceComparison, sequenceSearchIssue } from "./sequenceSearch.js";

export function useSequenceSearch(source, client = api) {
  const items = ref([]);
  const loading = ref(false);
  const error = ref("");
  const nextCursor = ref(null);
  const page = ref(0);
  let cursors = [null];
  let requestedPage = 0;
  let controller;
  let generation = 0;

  async function load(index = 0) {
    controller?.abort();
    const current = ++generation;
    const selection = source.value;
    const issue = sequenceSearchIssue(selection.sequence);
    if (issue || !["repeats", "spacers"].includes(selection.kind)) {
      error.value = issue || "Choose a repeat or spacer.";
      loading.value = false;
      return;
    }
    controller = new AbortController();
    requestedPage = index;
    loading.value = true;
    error.value = "";
    try {
      const result = await client.catalogPage(
        selection.kind,
        {
          q: selection.sequence.toUpperCase(),
          match: "similar",
          limit: 25,
          cursor: cursors[index],
        },
        { signal: controller.signal },
      );
      if (current !== generation) return;
      items.value = result.items.map((row) => {
        const comparison = sequenceComparison(selection.sequence, row.sequence);
        if (!comparison) throw new Error("The database returned an inconsistent sequence match.");
        return { ...row, comparison };
      });
      nextCursor.value = result.next_cursor;
      page.value = index;
    } catch (failure) {
      if (current === generation && failure.name !== "AbortError")
        error.value = failure.message || "Database search failed. Try again.";
    } finally {
      if (current === generation) loading.value = false;
    }
  }

  watch(
    () => [source.value.kind, source.value.sequence],
    () => {
      items.value = [];
      nextCursor.value = null;
      cursors = [null];
      page.value = 0;
      void load();
    },
    { immediate: true },
  );
  onBeforeUnmount(() => {
    generation += 1;
    controller?.abort();
  });
  return {
    items,
    loading,
    error,
    nextCursor,
    page,
    retry: () => load(requestedPage),
    previous: () => !loading.value && page.value > 0 && load(page.value - 1),
    next: () => {
      if (loading.value || !nextCursor.value) return;
      cursors[page.value + 1] = nextCursor.value;
      void load(page.value + 1);
    },
  };
}
