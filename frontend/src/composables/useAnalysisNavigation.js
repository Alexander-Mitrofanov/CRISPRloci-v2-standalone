import { nextTick, onBeforeUnmount, ref, watch } from "vue";

import { ANALYSIS_MODES } from "../science.js";

const validMethod = (value) => ANALYSIS_MODES.some((item) => item.id === value);

// Query-based pages work on static hosts without competing with #job recovery.
export function useAnalysisNavigation(hasSession, browser = window) {
  const readMethod = () => new URL(browser.location.href).searchParams.get("method");
  const initialMethod = readMethod();
  const isDatabase = () => new URL(browser.location.href).searchParams.get("view") === "database";
  const method = ref(validMethod(initialMethod) ? initialMethod : null);
  const page = ref(
    isDatabase() ? "database" : hasSession.value ? "results" : method.value ? "input" : "methods",
  );

  async function focusPage() {
    await nextTick();
    const heading = browser.document
      .getElementById("main-content")
      ?.querySelector("h1, #results-heading, #job-heading");
    heading?.focus({ preventScroll: true });
    heading?.scrollIntoView({ block: "start" });
  }

  function navigate(nextPage, nextMethod = method.value, replace = false) {
    if (nextPage === "input" && !validMethod(nextMethod)) return;
    const url = new URL(browser.location.href);
    if (nextPage === "database") url.searchParams.set("view", "database");
    else url.searchParams.delete("view");
    if (nextPage === "methods") url.searchParams.delete("method");
    else if (validMethod(nextMethod)) url.searchParams.set("method", nextMethod);
    browser.history[replace ? "replaceState" : "pushState"](
      { ...browser.history.state, analysisPage: nextPage },
      "",
      url,
    );
    method.value = nextMethod;
    page.value = nextPage;
    void focusPage();
  }

  function restorePage(event) {
    const requested = readMethod();
    if (validMethod(requested)) method.value = requested;
    const savedPage = event.state?.analysisPage;
    page.value = isDatabase()
      ? "database"
      : savedPage === "methods"
        ? "methods"
        : savedPage === "input" && validMethod(requested)
          ? "input"
          : hasSession.value
            ? "results"
            : validMethod(requested)
              ? "input"
              : "methods";
    void focusPage();
  }

  browser.addEventListener("popstate", restorePage);
  onBeforeUnmount(() => browser.removeEventListener("popstate", restorePage));
  watch(hasSession, (active) => {
    if (active && page.value !== "database") navigate("results", method.value, true);
    else if (page.value === "results") navigate("methods", method.value, true);
  });

  return {
    method,
    page,
    chooseMethod: (id) => navigate("input", id),
    showMethods: () => navigate("methods"),
    showResults: (replace = false) => navigate("results", method.value, replace === true),
    showDatabase: () => navigate("database"),
  };
}
