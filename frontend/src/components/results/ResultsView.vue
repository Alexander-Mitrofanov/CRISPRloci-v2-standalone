<script setup>
import { computed, nextTick, ref, useId, watch } from "vue";

import { normalizePublicResult, useAdapterMembership } from "../../features/results/index.js";
import { genomeAvailable } from "../../features/genome/index.js";
import GenomeExplorer from "./genome/GenomeExplorer.vue";
import { asArray } from "../../utils/formatting.js";
import { mergeAdapterMembership } from "../../utils/results.js";
import AppIcon from "../common/AppIcon.vue";
import LeaderResults from "./annotations/LeaderResults.vue";
import RepeatResults from "./repeats/RepeatResults.vue";
import RepeatTypeResults from "./repeats/RepeatTypeResults.vue";
import RepeatMapResults from "./repeats/RepeatMapResults.vue";
import ArrayCompareResults from "./comparison/ArrayCompareResults.vue";
import LocusCompareResults from "./comparison/LocusCompareResults.vue";
import ProphageContextResults from "./annotations/ProphageContextResults.vue";
import InsertionContextResults from "./InsertionContextResults.vue";
import DefenseContext from "./annotations/DefenseContext.vue";
import RnaFamilyResults from "./RnaFamilyResults.vue";
import CasContextRelations from "./annotations/CasContextRelations.vue";
import AssociationResults from "./protospacer/AssociationResults.vue";
import ViralResults from "./protospacer/ViralResults.vue";
import ProtospacerResults from "./protospacer/ProtospacerResults.vue";
import AnnotationResults from "./annotations/AnnotationResults.vue";
import DetectionSummary from "./overview/DetectionSummary.vue";
import EvidenceMap from "./overview/EvidenceMap.vue";
import PreflightSummary from "./overview/PreflightSummary.vue";
import ResultSynopsis from "./overview/ResultSynopsis.vue";
import OrientationResults from "./orientation/OrientationResults.vue";
import DownloadsPanel from "./publication/DownloadsPanel.vue";
import SequenceHandoff from "./publication/SequenceHandoff.vue";
import ProvenancePanel from "./publication/ProvenancePanel.vue";
import ReconstructionResults from "./reconstruction/ReconstructionResults.vue";

const props = defineProps({
  job: { type: Object, default: null },
  credential: { type: Object, default: null },
  maxArchiveBytes: { type: Number, default: 0 },
  exampleSnapshot: { type: Object, default: null },
  service: { type: Object, default: () => ({}) },
  limits: { type: Object, default: () => ({}) },
});
defineEmits(["continue-analysis"]);
const sourceSummary = computed(() =>
  normalizePublicResult(props.job?.summary || props.job?.result || {}, props.job || {}),
);
const currentJob = computed(() => props.job);
const currentCredential = computed(() => props.credential);
const { artifactGroups, membershipStatus } = useAdapterMembership(
  sourceSummary,
  currentJob,
  currentCredential,
);

const summary = computed(() => mergeAdapterMembership(sourceSummary.value, artifactGroups.value));
const completed = computed(() =>
  ["completed", "completed_no_eligible_groups"].includes(props.job?.status),
);
const noEligible = computed(() => props.job?.status === "completed_no_eligible_groups");
const detection = computed(() => summary.value.detection);
const arrays = computed(() => asArray(detection.value.arrays));
const annotationMode = computed(() =>
  ["cas", "tracrrna", "loci", "leader"].includes(props.job?.mode),
);
const repeatMode = computed(() => ["repeats", "repeat_context"].includes(props.job?.mode));
const comparisonMode = computed(() => props.job?.mode === "array_compare");
const contextModes = [
  "locus_compare",
  "prophage_context",
  "insertion_context",
  "defense_context",
  "rna_family",
];
const contextMode = computed(() => contextModes.includes(props.job?.mode));
const associationMode = computed(() => props.job?.mode === "spacer_association");
const typeMode = computed(() => props.job?.mode === "repeat_type");
const mapMode = computed(() => props.job?.mode === "repeat_map");
const viralMode = computed(() => props.job?.mode === "viral_search");
const protospacerMode = computed(() => props.job?.mode === "protospacer");
const hasDetection = computed(
  () =>
    ![
      "cas",
      "tracrrna",
      "repeats",
      "protospacer",
      "viral_search",
      "repeat_map",
      "repeat_type",
      "spacer_association",
      "array_compare",
      ...contextModes,
    ].includes(props.job?.mode),
);
const hasEvolution = computed(() => ["reconstruction", "orientation"].includes(props.job?.mode));
const cas = computed(() => summary.value.casandra || { status: "not_requested" });
const tracr = computed(() => summary.value.tracrrna || { status: "not_requested" });
const tabId = useId();
const selectedTab = ref("overview");
const genomeOpened = ref(false);
const genomeReady = ref(false);
const hasGenome = computed(
  () =>
    props.credential?.jobId === props.job?.job_id &&
    Boolean(props.credential?.accessToken) &&
    genomeAvailable(props.job?.mode, props.job?.artifacts),
);
const tabs = computed(() => [
  { id: "overview", label: "Overview" },
  ...(hasGenome.value ? [{ id: "genome", label: "Genome & locus" }] : []),
  ...(annotationMode.value && ["loci", "cas"].includes(props.job?.mode)
    ? [{ id: "cas", label: "Cas systems" }]
    : []),
  ...(annotationMode.value && ["loci", "tracrrna"].includes(props.job?.mode)
    ? [{ id: "tracrrna", label: "tracrRNA" }]
    : []),
  ...(["loci", "leader"].includes(props.job?.mode)
    ? [{ id: "leader", label: "Leader context" }]
    : []),
  ...(hasDetection.value ? [{ id: "arrays", label: "Arrays" }] : []),
  ...(repeatMode.value ? [{ id: "repeats", label: "Repeat evidence" }] : []),
  ...(comparisonMode.value ? [{ id: "array_compare", label: "Array comparison" }] : []),
  ...(contextMode.value
    ? [
        {
          id: props.job.mode,
          label: {
            locus_compare: "Gene neighborhoods",
            prophage_context: "Prophage context",
            insertion_context: "Insertion sequences",
            defense_context: "Defense repertoire",
            rna_family: "RNA family",
          }[props.job.mode],
        },
      ]
    : []),
  ...(associationMode.value ? [{ id: "associations", label: "Associations" }] : []),
  ...(typeMode.value ? [{ id: "repeat_type", label: "Repeat subtypes" }] : []),
  ...(mapMode.value ? [{ id: "repeat_map", label: "Repeat matches" }] : []),
  ...(viralMode.value ? [{ id: "viral", label: "Viral matches" }] : []),
  ...(protospacerMode.value ? [{ id: "protospacer", label: "Protospacer matches" }] : []),
  ...(props.job?.mode === "orientation" && !noEligible.value
    ? [{ id: "orientation", label: "Orientation" }]
    : []),
  ...(hasEvolution.value && !noEligible.value ? [{ id: "history", label: "History" }] : []),
  { id: "files", label: "Files & methods" },
]);
const activeTab = computed(() =>
  tabs.value.some((tab) => tab.id === selectedTab.value) ? selectedTab.value : "overview",
);
const warnings = computed(() => asArray(summary.value.warnings));

watch([() => props.job?.job_id, () => props.job?.mode], () => {
  selectedTab.value = "overview";
  genomeOpened.value = false;
  genomeReady.value = false;
});
watch(selectedTab, (value) => {
  if (value === "genome") genomeOpened.value = true;
});

function openGenomeEvidence() {
  selectedTab.value =
    { loci: "arrays", detection: "arrays", protospacer: "protospacer" }[props.job.mode] ||
    props.job.mode;
}

async function navigateTabs(event, index) {
  const last = tabs.value.length - 1;
  const next = {
    ArrowRight: (index + 1) % tabs.value.length,
    ArrowLeft: (index + last) % tabs.value.length,
    Home: 0,
    End: last,
  }[event.key];
  if (next == null) return;
  event.preventDefault();
  selectedTab.value = tabs.value[next].id;
  await nextTick();
  document.getElementById(`${tabId}-tab-${selectedTab.value}`)?.focus();
}
</script>

<template>
  <section
    v-if="completed"
    class="results"
    aria-labelledby="results-heading"
  >
    <div class="results-title">
      <div>
        <h2
          id="results-heading"
          tabindex="-1"
        >
          Analysis results
        </h2>
      </div>
      <span class="complete-stamp"><AppIcon name="check" /> Completed</span>
    </div>
    <div
      class="result-tabs"
      role="tablist"
      aria-label="Result sections"
    >
      <button
        v-for="(tab, index) in tabs"
        :id="`${tabId}-tab-${tab.id}`"
        :key="tab.id"
        type="button"
        role="tab"
        :aria-selected="activeTab === tab.id"
        :aria-controls="`${tabId}-panel-${tab.id}`"
        :tabindex="activeTab === tab.id ? 0 : -1"
        @click="selectedTab = tab.id"
        @keydown="navigateTabs($event, index)"
      >
        {{ tab.label }}
      </button>
    </div>
    <div
      v-if="noEligible"
      class="no-eligible"
      role="status"
    >
      <AppIcon name="info" />
      <div>
        <strong>No eligible evolutionary groups</strong>
        <p>
          Detection completed. No group met the selected filters; evolution and orientation were not
          evaluated.
        </p>
      </div>
    </div>
    <ul
      v-if="warnings.length"
      class="warning-list"
      aria-label="Analysis warnings"
    >
      <li
        v-for="(warning, index) in warnings"
        :key="index"
      >
        <AppIcon
          name="warning"
          :size="18"
        />
        <span
          ><strong>{{ warning.title || warning.code || "Analysis warning" }}</strong
          >{{ warning.message || String(warning) }}</span
        >
      </li>
    </ul>
    <section
      v-for="tab in tabs"
      :id="`${tabId}-panel-${tab.id}`"
      :key="tab.id"
      class="result-tab-panel"
      role="tabpanel"
      :aria-labelledby="`${tabId}-tab-${tab.id}`"
      :hidden="activeTab !== tab.id"
      tabindex="0"
    >
      <button
        v-if="tab.id === 'overview' && hasGenome"
        type="button"
        @click="selectedTab = 'genome'"
      >
        Explore genome context
      </button>
      <GenomeExplorer
        v-if="tab.id === 'genome' && genomeOpened"
        :key="job.job_id"
        :job-id="job.job_id"
        :mode="job.mode"
        :access-token="credential.accessToken"
        :artifacts="job.artifacts"
        :active="activeTab === 'genome'"
        @open-evidence="openGenomeEvidence"
        @ready="genomeReady = $event"
      />
      <ResultSynopsis
        v-if="
          tab.id === 'overview' &&
          !annotationMode &&
          !repeatMode &&
          !protospacerMode &&
          !viralMode &&
          !mapMode &&
          !typeMode &&
          !associationMode &&
          !comparisonMode &&
          !contextMode
        "
        :summary="summary"
        :no-eligible="noEligible"
        :has-evolution="hasEvolution"
        :has-orientation="job.mode === 'orientation'"
      />
      <AnnotationResults
        v-if="
          annotationMode &&
          job.mode !== 'leader' &&
          ['overview', 'cas', 'tracrrna'].includes(tab.id)
        "
        :cas="cas"
        :tracr="tracr"
        :arrays="arrays"
        :view="tab.id"
        :show-map="!hasGenome || !genomeReady"
      />
      <LeaderResults
        v-if="['leader', 'loci'].includes(job.mode) && ['overview', 'leader'].includes(tab.id)"
        :evidence="summary.crisprleader"
        :compact="tab.id === 'overview'"
      />
      <RepeatResults
        v-if="repeatMode && ['overview', 'repeats'].includes(tab.id)"
        :evidence="summary.repeats"
        :molecule="summary.repeats.molecule || job.options?.molecule || 'DNA'"
        :with-context="job.mode === 'repeat_context'"
        :compact="tab.id === 'overview'"
      />
      <DetectionSummary
        v-if="tab.id === 'arrays'"
        :key="job.job_id"
        :summary="detection"
        :arrays="arrays"
      />
      <ArrayCompareResults
        v-if="comparisonMode && ['overview', 'array_compare'].includes(tab.id)"
        v-bind="summary.array_compare"
        :compact="tab.id === 'overview'"
      />
      <LocusCompareResults
        v-if="job.mode === 'locus_compare' && ['overview', 'locus_compare'].includes(tab.id)"
        :available="summary.locus_compare.available"
        :loci="summary.locus_compare.loci"
        :links="summary.locus_compare.links"
        :counts="summary.locus_compare.counts"
        :identity-cutoff="summary.locus_compare.identityCutoff"
        :truncated="summary.locus_compare.truncated"
        :compact="tab.id === 'overview'"
      />
      <ProphageContextResults
        v-if="job.mode === 'prophage_context' && ['overview', 'prophage_context'].includes(tab.id)"
        v-bind="summary.prophage_context"
        :compact="tab.id === 'overview'"
      />
      <template
        v-if="
          job.mode === 'insertion_context' && ['overview', 'insertion_context'].includes(tab.id)
        "
      >
        <InsertionContextResults
          v-if="summary.insertion_context.available"
          :elements="summary.insertion_context.elements"
          :counts="summary.insertion_context.counts"
          :truncated="summary.insertion_context.truncated"
        />
        <p v-else>Validated insertion-sequence evidence is unavailable.</p>
      </template>
      <DefenseContext
        v-if="job.mode === 'defense_context' && ['overview', 'defense_context'].includes(tab.id)"
        :available="summary.defense_context.available"
        :systems="summary.defense_context.systems"
        :components="summary.defense_context.components"
        :hits="summary.defense_context.hits"
        :genes="summary.defense_context.genes"
        :counts="summary.defense_context.counts"
        :systems-truncated="summary.defense_context.systemsTruncated"
        :hits-truncated="summary.defense_context.hitsTruncated"
      />
      <template v-if="job.mode === 'rna_family' && ['overview', 'rna_family'].includes(tab.id)">
        <RnaFamilyResults
          v-if="summary.rna_family.available"
          :counts="summary.rna_family.counts"
          :policy="summary.rna_family.policy"
          :fold="summary.rna_family.fold"
          :alignment="summary.rna_family.alignment"
          :pairs="summary.rna_family.pairs"
          :warnings="summary.rna_family.warnings"
        />
        <p v-else>Validated RNA-family evidence is unavailable.</p>
      </template>
      <CasContextRelations
        v-if="['locus_compare', 'insertion_context', 'defense_context'].includes(tab.id)"
        :relations="summary[tab.id].casContext"
      />
      <AssociationResults
        v-if="associationMode && ['overview', 'associations'].includes(tab.id)"
        v-bind="summary.spacer_association"
        :compact="tab.id === 'overview'"
      />
      <RepeatTypeResults
        v-if="typeMode && ['overview', 'repeat_type'].includes(tab.id)"
        v-bind="summary.repeat_type"
        :compact="tab.id === 'overview'"
      />
      <RepeatMapResults
        v-if="mapMode && ['overview', 'repeat_map'].includes(tab.id)"
        :available="summary.repeat_map.available"
        :namespace="summary.repeat_map.namespace"
        :reference-sha256="summary.repeat_map.referenceSha256"
        :max-distance="summary.repeat_map.maxDistance"
        :top-k="summary.repeat_map.topK"
        :counts="summary.repeat_map.counts"
        :results="summary.repeat_map.results"
        :truncated="summary.repeat_map.truncated"
        :compact="tab.id === 'overview'"
      />
      <ViralResults
        v-if="viralMode && ['overview', 'viral'].includes(tab.id)"
        :compact="tab.id === 'overview'"
        :available="summary.viral_search.available"
        :completed="summary.viral_search.execution_completed"
        :outcome="summary.viral_search.outcome"
        :reference="summary.viral_search.reference"
        :counts="summary.viral_search.counts"
        :matches="summary.viral_search.matches"
        :candidates="summary.viral_search.candidate_viruses"
        :queries="summary.viral_search.queries"
        :substitutions="summary.viral_search.max_substitutions"
        :truncated="summary.viral_search.truncated"
      />
      <ProtospacerResults
        v-if="protospacerMode && ['overview', 'protospacer'].includes(tab.id)"
        :available="summary.protospacer.available"
        :search-complete="summary.protospacer.search_complete"
        :outcome="summary.protospacer.outcome"
        :reference="summary.protospacer.reference"
        :counts="summary.protospacer.counts"
        :matches="summary.protospacer.matches"
        :truncated="summary.protospacer.matches_truncated"
        :compact="tab.id === 'overview'"
      />
      <OrientationResults
        v-if="tab.id === 'orientation'"
        :summary="summary"
      />
      <ReconstructionResults
        v-if="tab.id === 'history'"
        :summary="summary"
      />
      <SequenceHandoff
        v-if="tab.id === 'files' && credential && !exampleSnapshot"
        :job="job"
        :credential="credential"
        :service="service"
        :limits="limits"
        @continue-analysis="$emit('continue-analysis', $event)"
      />
      <DownloadsPanel
        v-if="tab.id === 'files' && credential && !exampleSnapshot"
        :job="job"
        :credential="credential"
        :max-archive-bytes="maxArchiveBytes"
      />
      <details
        v-if="tab.id === 'files' && summary.publication_notes?.length"
        class="result-details publication-notes"
      >
        <summary>Publication notes ({{ summary.publication_notes.length }})</summary>
        <ul>
          <li
            v-for="(note, index) in summary.publication_notes"
            :key="index"
          >
            {{ note }}
          </li>
        </ul>
      </details>
      <details
        v-if="tab.id === 'files' && hasEvolution"
        class="result-details"
      >
        <summary>Filtering &amp; evidence</summary>
        <PreflightSummary :summary="summary" />
        <EvidenceMap
          v-if="!noEligible"
          :summary="summary"
          :membership-status="membershipStatus"
        />
      </details>
      <details
        v-if="tab.id === 'files'"
        class="result-details"
      >
        <summary>Methods &amp; provenance</summary>
        <ProvenancePanel
          :summary="summary"
          :show-warnings="false"
        />
      </details>
    </section>
  </section>
</template>
