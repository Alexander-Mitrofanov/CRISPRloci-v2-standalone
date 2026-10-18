<script setup>
import { computed, ref, useId, watch } from "vue";

import { asArray, finiteMetric, formatNumber } from "../../../utils/formatting.js";
import {
  entryRootGains,
  entryTreeHeight,
  historyLossCount,
  MAX_HISTORY_LEAVES,
  MAX_HISTORY_NODES,
  MAX_HISTORY_SPACER_COLUMNS,
  publicNodeName,
} from "../../../utils/history.js";
import { layoutNewick } from "../../../utils/newick.js";
import { comparisonDecisionFor, groupIdentity } from "../../../utils/results.js";
import AppIcon from "../../common/AppIcon.vue";
import HistoryCanvas from "./HistoryCanvas.vue";
import HistoryExactTable from "./HistoryExactTable.vue";
import HistoryNodeBrowser from "./HistoryNodeBrowser.vue";
import TreeDiagram from "./TreeDiagram.vue";

const props = defineProps({
  summary: { type: Object, required: true },
  group: { type: String, required: true },
  fallbackTree: { type: String, default: "" },
  reportedByDefault: Boolean,
});
const orientation = computed(() => props.summary.orientation || {});
const reconstructions = computed(() =>
  asArray(orientation.value.reconstructions).filter(
    (entry) =>
      String(entry?.group) === props.group && ["input", "reverse"].includes(entry?.hypothesis),
  ),
);
const comparison = computed(() =>
  asArray(orientation.value.comparisons).find(
    (item, index) => groupIdentity(item, index) === props.group,
  ),
);
const decision = computed(() =>
  comparison.value
    ? comparisonDecisionFor(comparison.value, orientation.value).label
    : "Unresolved",
);
const supportedHypothesis = computed(() =>
  decision.value === "Input order supported"
    ? "input"
    : decision.value === "Reverse input order supported"
      ? "reverse"
      : null,
);
const requiredHypothesis = computed(() => supportedHypothesis.value || "input");
const requiredHistoryAvailable = computed(() =>
  reconstructions.value.some((row) => row.hypothesis === requiredHypothesis.value),
);
const hypothesis = ref("input");
const scaleMode = ref("topology");
const fitCanvas = ref(
  Boolean(typeof window !== "undefined" && window.matchMedia?.("(max-width: 980px)")?.matches),
);
const selectedNodeName = ref("");
const nodeInspector = ref(null);
const descriptionId = useId();

watch(
  [reconstructions, requiredHypothesis],
  ([rows, required]) => {
    hypothesis.value = rows.some((row) => row.hypothesis === required)
      ? required
      : rows[0]?.hypothesis || required;
    selectedNodeName.value = "";
  },
  { immediate: true },
);

const entry = computed(
  () =>
    reconstructions.value.find((item) => item.hypothesis === hypothesis.value) ||
    reconstructions.value[0],
);
const completeOrder = computed(() =>
  asArray(entry.value?.spacer_order)
    .map(Number)
    .filter((value) => Number.isInteger(value) && value > 0),
);
const spacerOrder = computed(() => completeOrder.value.slice(0, MAX_HISTORY_SPACER_COLUMNS));
const omittedSpacerColumns = computed(() => completeOrder.value.length - spacerOrder.value.length);
const treeHeights = computed(() =>
  reconstructions.value.map(entryTreeHeight).filter((value) => value != null),
);
const sharedDistance = computed(() => Math.max(0, ...treeHeights.value));
const provisional = computed(() => layoutNewick(entry.value?.newick, 100));
const exceedsSizeLimit = computed(() =>
  Boolean(
    provisional.value &&
    (provisional.value.leaves.length > MAX_HISTORY_LEAVES ||
      provisional.value.nodes.length > MAX_HISTORY_NODES),
  ),
);
const height = computed(() => Math.max(310, (provisional.value?.leaves.length || 1) * 74 + 78));
const layout = computed(() =>
  provisional.value && !exceedsSizeLimit.value
    ? layoutNewick(entry.value?.newick, height.value, scaleMode.value, sharedDistance.value, 280)
    : null,
);
const nodeData = computed(() => asArray(entry.value?.nodes));
const dataByName = computed(
  () => new Map(nodeData.value.map((node) => [String(node?.name), node])),
);
const leafNames = computed(
  () => new Set(layout.value?.leaves.map((leaf) => String(leaf.name)) || []),
);
const internalData = computed(() =>
  nodeData.value.filter((node) => !leafNames.value.has(String(node?.name))),
);
const rootData = computed(() =>
  layout.value
    ? dataByName.value.get(String(layout.value.tree.name)) ||
      internalData.value[0] ||
      nodeData.value[0]
    : nodeData.value[0],
);
const selectedNode = computed(
  () =>
    nodeData.value.find((node) => String(node?.name) === selectedNodeName.value) || rootData.value,
);
const effectiveSelectedName = computed(() => String(selectedNode.value?.name || ""));
const likelihood = (kind) =>
  finiteMetric(
    comparison.value?.[
      kind === "input" ? "forward_ln_likelihood_bdm" : "reverse_ln_likelihood_bdm"
    ],
  );
const hypothesisLabel = (kind) =>
  kind === "input" ? "Input spacer order" : "Reversed spacer order";
const candidateFor = (kind) => reconstructions.value.find((item) => item.hypothesis === kind);
const inputEntry = computed(() => candidateFor("input"));
const reverseEntry = computed(() => candidateFor("reverse"));
const historyStatus = computed(() => {
  if (!requiredHistoryAvailable.value) return "Inspection only";
  if (supportedHypothesis.value)
    return hypothesis.value === supportedHypothesis.value
      ? "Supported history"
      : "Comparison history";
  return hypothesis.value === "input" ? "Reported default" : "Comparison history";
});
const unavailableTitle = computed(() =>
  supportedHypothesis.value
    ? `${hypothesisLabel(requiredHypothesis.value)} is supported, but its structured reconstruction is unavailable.`
    : "The input-order default reconstruction is unavailable for this unresolved group.",
);
const branchSummary = computed(
  () =>
    layout.value?.nodes
      .filter((node) => node.parent)
      .map((node) => {
        const data = dataByName.value.get(String(node.name));
        return [
          `${publicNodeName(node.parent.name)} to ${publicNodeName(node.name)}:`,
          `${asArray(data?.gains).length} gains,`,
          `${historyLossCount(data)} losses,`,
          `branch length ${formatNumber(node.length, 6)}`,
        ].join(" ");
      })
      .join("; ") || "No branch data available.",
);
const chooseNode = (name) => {
  selectedNodeName.value = String(name || "");
  if (nodeInspector.value) nodeInspector.value.open = true;
};
const chooseHypothesis = (kind) => {
  if (candidateFor(kind)) {
    hypothesis.value = kind;
    selectedNodeName.value = "";
  }
};
</script>

<template>
  <div
    v-if="reconstructions.length"
    class="history-explorer"
  >
    <p
      v-if="layout"
      :id="descriptionId"
      class="sr-only"
    >
      Ancestral reconstruction for {{ group }} under {{ hypothesisLabel(entry.hypothesis) }}. Leaves
      in display order: {{ layout.leaves.map((leaf) => publicNodeName(leaf.name)).join(", ") }}.
      Branches: {{ branchSummary }}. Use the node browser and exact node data table following the
      visual to inspect ordered spacer identities and events.
    </p>
    <div class="history-heading">
      <div>
        <h5>Ancestral arrays</h5>
      </div>
      <div class="history-legend">
        <span><i class="legend-gain" /> acquisition</span
        ><span><i class="legend-loss" /> deletion</span
        ><span><i class="legend-absence" /> absent</span>
      </div>
    </div>
    <div
      v-if="!requiredHistoryAvailable"
      class="history-availability"
      role="note"
    >
      <AppIcon
        name="warning"
        :size="18"
      />
      <p>
        <strong>{{ unavailableTitle }}</strong> The available hypothesis is shown for inspection
        only; it is not substituted for the missing reported history. Review the workflow warning
        and detailed artifacts.
      </p>
    </div>
    <div class="history-summary">
      <span
        ><small>History status</small><strong>{{ historyStatus }}</strong></span
      ><span
        ><small>Acquisitions</small
        ><strong>{{ formatNumber(entry?.acquisition_count) }}</strong></span
      ><span
        ><small>Deletions</small><strong>{{ formatNumber(entry?.deletion_count) }}</strong></span
      ><span
        ><small>Hypothesis shown</small><strong>{{ hypothesisLabel(hypothesis) }}</strong></span
      >
    </div>
    <details class="history-settings">
      <summary>History controls</summary>
      <div
        class="history-controls"
        :aria-label="`Reconstruction hypotheses for ${group}`"
      >
        <button
          v-for="kind in ['input', 'reverse']"
          :key="kind"
          type="button"
          :class="['history-hypothesis', { active: hypothesis === kind }]"
          :aria-pressed="hypothesis === kind"
          :disabled="!candidateFor(kind)"
          @click="chooseHypothesis(kind)"
        >
          <span
            >{{ hypothesisLabel(kind) }}<b v-if="kind === supportedHypothesis">supported</b
            ><b v-else-if="!supportedHypothesis && kind === 'input'">reported default</b
            ><b v-else-if="!candidateFor(kind)">history unavailable</b></span
          ><strong
            >{{ formatNumber(candidateFor(kind)?.acquisition_count) }} gains ·
            {{ formatNumber(candidateFor(kind)?.deletion_count) }} losses</strong
          ><small
            >BDM lnL {{ formatNumber(likelihood(kind), 3) }} · max root-to-tip
            {{ formatNumber(entryTreeHeight(candidateFor(kind)), 6) }}</small
          >
        </button>
      </div>
      <details
        v-if="inputEntry && reverseEntry"
        class="history-explanation"
      >
        <summary>Compare histories</summary>
        <div class="history-contrast">
          <AppIcon
            name="info"
            :size="18"
          />
          <p>
            <strong>Why the histories differ:</strong> input order needs
            {{ formatNumber(inputEntry.deletion_count) }} inferred deletions and places
            {{ formatNumber(entryRootGains(inputEntry)) }} acquisition{{
              entryRootGains(inputEntry) === 1 ? "" : "s"
            }}
            at the root; reversed order needs
            {{ formatNumber(reverseEntry.deletion_count) }} deletions and places
            {{ formatNumber(entryRootGains(reverseEntry)) }} at the root. CRISPR-evOr compares the
            full model likelihoods, not counts alone.
          </p>
        </div>
      </details>
      <div class="history-scale">
        <span>Tree layout</span
        ><button
          type="button"
          :class="{ active: scaleMode === 'topology' }"
          :aria-pressed="scaleMode === 'topology'"
          @click="scaleMode = 'topology'"
        >
          Readable topology</button
        ><button
          type="button"
          :class="{ active: scaleMode === 'branch' }"
          :aria-pressed="scaleMode === 'branch'"
          @click="scaleMode = 'branch'"
        >
          Shared branch scale</button
        ><small
          >Root-to-tip: {{ formatNumber(sharedDistance, 6) }}.<span v-if="omittedSpacerColumns > 0">
            Showing the first {{ MAX_HISTORY_SPACER_COLUMNS }} of {{ completeOrder.length }} spacer
            columns.</span
          ></small
        >
      </div>
      <div
        class="history-view-controls"
        role="group"
        aria-label="Ancestral history canvas view"
      >
        <span>Canvas view</span
        ><button
          type="button"
          :class="{ active: fitCanvas }"
          :aria-pressed="fitCanvas"
          @click="fitCanvas = true"
        >
          Fit overview</button
        ><button
          type="button"
          :class="{ active: !fitCanvas }"
          :aria-pressed="!fitCanvas"
          @click="fitCanvas = false"
        >
          Readable detail</button
        ><small>{{
          fitCanvas
            ? "Full tree; switch to detail for labels."
            : "Scroll horizontally to inspect every spacer."
        }}</small>
      </div>
    </details>
    <div
      v-if="exceedsSizeLimit"
      class="history-size-limit"
    >
      <AppIcon name="info" />
      <div>
        <strong>Structured history available in the result bundle</strong>
        <p>
          This group contains {{ formatNumber(provisional.leaves.length) }} leaves and
          {{ formatNumber(provisional.nodes.length) }} tree nodes. The inline browser is capped at
          {{ MAX_HISTORY_LEAVES }} leaves and {{ MAX_HISTORY_NODES }} nodes to keep this page
          responsive; use the detailed JSON/Newick artifacts for the complete reconstruction.
        </p>
      </div>
    </div>
    <template v-else-if="layout">
      <HistoryCanvas
        :layout="layout"
        :node-data="nodeData"
        :spacer-order="spacerOrder"
        :selected-node-name="effectiveSelectedName"
        :fit-canvas="fitCanvas"
        :group="group"
        :hypothesis-label="hypothesisLabel(entry.hypothesis)"
        :description-id="descriptionId"
        @select-node="chooseNode"
      />
      <details
        ref="nodeInspector"
        class="history-inspector"
      >
        <summary>Inspect nodes</summary>
        <HistoryNodeBrowser
          :layout="layout"
          :node-data="nodeData"
          :spacer-order="spacerOrder"
          :complete-order="completeOrder"
          :selected-node-name="effectiveSelectedName"
          @select-node="chooseNode"
        />
      </details>
      <HistoryExactTable
        :layout="layout"
        :node-data="nodeData"
        :group="group"
      />
      <p class="history-caveat">
        Array-derived estimates, not an independent organismal phylogeny. Orientation support does
        not establish transcription direction or a leader sequence.
      </p>
    </template>
    <TreeDiagram
      v-else-if="fallbackTree"
      :newick="fallbackTree"
      :group="group"
      :reported-by-default="reportedByDefault"
      :status-label="`${historyStatus} model tree`"
    />
  </div>
  <TreeDiagram
    v-else-if="fallbackTree"
    :newick="fallbackTree"
    :group="group"
    :reported-by-default="reportedByDefault"
    :status-label="
      comparison ? 'Displayed hypothesis model tree' : 'Input orientation; not compared'
    "
  />
</template>
