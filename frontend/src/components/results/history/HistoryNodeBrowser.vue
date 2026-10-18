<script setup>
import { computed } from "vue";

import { asArray } from "../../../utils/formatting.js";
import { historyLossCount, publicNodeName, spacerColorClass } from "../../../utils/history.js";

const props = defineProps({
  layout: { type: Object, required: true },
  nodeData: { type: Array, required: true },
  spacerOrder: { type: Array, required: true },
  completeOrder: { type: Array, required: true },
  selectedNodeName: { type: String, default: "" },
});
const emit = defineEmits(["select-node"]);
const leafNames = computed(() => new Set(props.layout.leaves.map((leaf) => String(leaf.name))));
const dataByName = computed(
  () => new Map(props.nodeData.map((node) => [String(node?.name), node])),
);
const internalData = computed(() =>
  props.nodeData.filter((node) => !leafNames.value.has(String(node?.name))),
);
const rootData = computed(
  () =>
    dataByName.value.get(String(props.layout.tree.name)) ||
    internalData.value[0] ||
    props.nodeData[0],
);
const selectedNode = computed(
  () =>
    props.nodeData.find((node) => String(node?.name) === props.selectedNodeName) || rootData.value,
);
const selectedSpacers = computed(() => new Set(asArray(selectedNode.value?.spacers).map(Number)));
const selectedSpacerIds = computed(() =>
  asArray(selectedNode.value?.spacers)
    .map(Number)
    .filter((value) => Number.isInteger(value) && value > 0),
);
const shownSelectedSpacers = computed(
  () => props.spacerOrder.filter((spacer) => selectedSpacers.value.has(spacer)).length,
);
const nodeKind = computed(() =>
  leafNames.value.has(String(selectedNode.value?.name))
    ? "Observed leaf"
    : String(selectedNode.value?.name) === String(props.layout.tree.name)
      ? "Inferred root"
      : "Inferred ancestor",
);
const spacerStripLabel = computed(() =>
  [
    `${nodeKind.value} ${publicNodeName(selectedNode.value?.name)} contains`,
    `${selectedSpacerIds.value.length} reconstructed spacers,`,
    `spacer IDs ${selectedSpacerIds.value.length ? selectedSpacerIds.value.join(", ") : "none"};`,
    `${shownSelectedSpacers.value} are visible in`,
    `${props.spacerOrder.length} of ${props.completeOrder.length} available columns`,
  ].join(" "),
);
</script>

<template>
  <div class="ancestor-browser">
    <div
      class="ancestor-tabs"
      role="group"
      aria-label="Observed and reconstructed nodes"
    >
      <button
        v-for="node in nodeData"
        :key="String(node.name)"
        type="button"
        :class="{ active: String(selectedNode?.name) === String(node.name) }"
        :aria-pressed="String(selectedNode?.name) === String(node.name)"
        @click="emit('select-node', node.name)"
      >
        {{ publicNodeName(node.name)
        }}<small
          >{{ leafNames.has(String(node.name)) ? "observed" : "inferred" }} ·
          {{ asArray(node.spacers).length }} spacers</small
        >
      </button>
    </div>
    <div
      v-if="selectedNode"
      class="ancestor-state"
    >
      <div>
        <span>{{ nodeKind }}</span
        ><strong>{{ publicNodeName(selectedNode.name) }}</strong
        ><small
          >{{ asArray(selectedNode.spacers).length }} spacers · +{{
            asArray(selectedNode.gains).length
          }}
          gains · −{{ historyLossCount(selectedNode) }} losses on the incoming branch<span
            v-if="completeOrder.length > spacerOrder.length"
          >
            · {{ shownSelectedSpacers }} present spacers shown in the first
            {{ spacerOrder.length }} columns</span
          ></small
        >
      </div>
      <div
        class="history-spacer-strip"
        role="img"
        tabindex="0"
        :aria-label="spacerStripLabel"
      >
        <i
          v-for="spacer in spacerOrder"
          :key="spacer"
          :class="selectedSpacers.has(spacer) ? ['present', spacerColorClass(spacer)] : ['absent']"
          ><span>{{ spacer }}</span></i
        >
      </div>
    </div>
  </div>
</template>
