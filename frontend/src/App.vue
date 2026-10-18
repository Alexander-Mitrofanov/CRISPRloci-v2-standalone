<script setup>
import { computed, ref } from "vue";
import ModeSelector from "./components/submission/ModeSelector.vue";
import { useAnalysisNavigation } from "./composables/useAnalysisNavigation.js";
import AnalysisForm from "./components/submission/AnalysisForm.vue";
import JobProgress from "./components/jobs/JobProgress.vue";
import ResultsView from "./components/results/ResultsView.vue";
import HeroHeader from "./components/shell/HeroHeader.vue";
import ReferencesSection from "./components/shell/ReferencesSection.vue";
import { useJobSession } from "./composables/useJobSession.js";
import { useServiceConfig } from "./composables/useServiceConfig.js";
import CatalogView from "./components/catalog/CatalogView.vue";
import { exampleAssetsFor } from "./example.js";

const { service, limits } = useServiceConfig();
const {
  credential,
  job,
  exampleSnapshot,
  pollError,
  cancelling,
  onSubmitted,
  onExampleLoaded,
  cancel,
  forget,
} = useJobSession();
const hasSession = computed(() => Boolean(exampleSnapshot.value || credential.value));
const publicBase = import.meta.env.BASE_URL;
const exampleAssets = computed(
  () =>
    exampleSnapshot.value &&
    exampleAssetsFor(
      exampleSnapshot.value.job.mode,
      exampleSnapshot.value.example.submission_options,
    ),
);
const references = ref(null);
const handoff = ref(null);
let handoffRevision = 0;
const { method, page, chooseMethod, showMethods, showResults, showDatabase } =
  useAnalysisNavigation(hasSession);

function submitted(nextCredential, initialJob) {
  onSubmitted(nextCredential, initialJob);
  if (page.value !== "database") showResults(true);
}

function exampleLoaded(snapshot) {
  onExampleLoaded(snapshot);
  if (snapshot && page.value !== "database") showResults(true);
}

function leaveJob() {
  handoff.value = null;
  forget();
  showMethods();
}

function continueAnalysis(input) {
  handoff.value = { ...input, revision: ++handoffRevision };
  chooseMethod(input.mode);
}
const activeHandoff = computed(() =>
  handoff.value?.mode === method.value && handoff.value.sourceJobId === job.value?.job_id
    ? handoff.value
    : null,
);
const blocksSubmission = computed(
  () =>
    Boolean(credential.value) &&
    !(
      activeHandoff.value &&
      ["completed", "completed_no_eligible_groups"].includes(job.value?.status)
    ),
);
</script>

<template>
  <div class="site-shell">
    <HeroHeader
      :service="service"
      :database-active="page === 'database'"
      @analyze="showMethods"
      @database="showDatabase"
      @references="references?.open()"
    />
    <main
      id="main-content"
      tabindex="-1"
    >
      <button
        v-if="hasSession && page !== 'results'"
        class="session-return"
        type="button"
        @click="showResults"
      >
        Return to current results
      </button>
      <ModeSelector
        v-if="page === 'methods'"
        :service="service"
        @select="chooseMethod"
        @database="showDatabase"
      />
      <KeepAlive>
        <CatalogView v-if="page === 'database'" />
      </KeepAlive>
      <KeepAlive :max="16">
        <AnalysisForm
          v-if="page === 'input' && method"
          :key="activeHandoff ? `${method}:${activeHandoff.revision}` : method"
          :initial-mode="method"
          :initial-input="activeHandoff"
          :service="service"
          :limits="limits"
          :has-active-job="blocksSubmission"
          @back="showMethods"
          @change-mode="chooseMethod"
          @submitted="submitted"
          @example-loaded="exampleLoaded"
        />
      </KeepAlive>
      <template v-if="exampleSnapshot && page === 'results'"
        ><p
          class="example-context"
          role="status"
        >
          Precomputed example result ready.
          {{ exampleSnapshot.example.source_note }}
          <a
            v-if="exampleAssets"
            :href="`${publicBase}${exampleAssets.input}`"
            download
            >Download example input</a
          >
          ·
          <a
            v-if="exampleAssets"
            :href="`${publicBase}${exampleAssets.result}`"
            download
            >Download precomputed output</a
          >
        </p>
        <div
          id="example-result"
          class="example-anchor"
        >
          <ResultsView
            :job="exampleSnapshot.job"
            :example-snapshot="exampleSnapshot"
          /></div
      ></template>
      <div
        v-if="credential && page === 'results'"
        id="job-status"
        class="job-anchor"
      >
        <JobProgress
          :job="job || { status: 'queued' }"
          :credential="credential"
          :cancelling="cancelling"
          @cancel="cancel"
          @forget="leaveJob"
        /><ResultsView
          :job="job"
          :credential="credential"
          :max-archive-bytes="limits.maxArchiveBytes"
          :service="service"
          :limits="limits"
          @continue-analysis="continueAnalysis"
        />
      </div>
      <p
        v-if="pollError"
        class="poll-error"
        role="alert"
      >
        {{ pollError }}
      </p>
      <ReferencesSection ref="references" />
    </main>
  </div>
</template>
