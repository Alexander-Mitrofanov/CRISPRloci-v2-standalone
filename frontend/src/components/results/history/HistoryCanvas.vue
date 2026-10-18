<script setup>
import { computed } from "vue";

import { asArray } from "../../../utils/formatting.js";
import { historyLossCount, publicNodeName, spacerColorClass } from "../../../utils/history.js";

const props = defineProps({
  layout: { type: Object, required: true },
  nodeData: { type: Array, required: true },
  spacerOrder: { type: Array, required: true },
  selectedNodeName: { type: String, default: "" },
  fitCanvas: Boolean,
  group: { type: String, required: true },
  hypothesisLabel: { type: String, required: true },
  descriptionId: { type: String, required: true },
});
const emit = defineEmits(["select-node"]);
const dataByName = computed(
  () => new Map(props.nodeData.map((node) => [String(node?.name), node])),
);
const cellSize = computed(() => (props.spacerOrder.length > 50 ? 10 : 13));
const arrayStart = 485;
const width = computed(() =>
  Math.max(930, arrayStart + props.spacerOrder.length * cellSize.value + 30),
);
const dataFor = (node) => dataByName.value.get(String(node?.name)) || {};
const presentSpacers = (node) => new Set(asArray(dataFor(node).spacers).map(Number));
const chooseNode = (node) => emit("select-node", String(node?.name || ""));
</script>

<template>
  <div
    :class="['history-canvas', { 'is-fit': fitCanvas }]"
    role="region"
    tabindex="0"
    :aria-label="`Scrollable ancestral reconstruction canvas for ${group}`"
  >
    <svg
      :width="width"
      :height="layout.height"
      :viewBox="`0 0 ${width} ${layout.height}`"
      role="img"
      :aria-describedby="descriptionId"
      :aria-label="`Ancestral reconstruction for ${group} under ${hypothesisLabel} with ${layout.leaves.length} observed leaves`"
    >
      <text
        :x="arrayStart"
        y="24"
        class="history-axis-label"
      >
        aligned spacer identity →
      </text>
      <text
        v-for="(spacer, index) in spacerOrder"
        :key="`column-${spacer}`"
        :x="arrayStart + index * cellSize + (cellSize - 2) / 2"
        y="41"
        text-anchor="middle"
        class="history-column-label"
      >
        {{ spacer }}
      </text>
      <g
        v-for="node in layout.nodes.filter((item) => item.children.length)"
        :key="`history-edges-${node.id}`"
      >
        <line
          :x1="node.x"
          :x2="node.x"
          :y1="Math.min(...node.children.map((child) => child.y))"
          :y2="Math.max(...node.children.map((child) => child.y))"
          class="history-tree-line"
        />
        <line
          v-for="child in node.children"
          :key="child.id"
          :x1="node.x"
          :x2="child.x"
          :y1="child.y"
          :y2="child.y"
          class="history-tree-line"
        />
      </g>
      <rect
        v-for="leaf in layout.leaves"
        :key="`row-${leaf.id}`"
        x="320"
        :y="leaf.y - 18"
        :width="width - 338"
        height="36"
        rx="4"
        class="history-leaf-row"
      />
      <template
        v-for="node in layout.nodes.filter((item) => item.parent)"
        :key="`events-${node.id}`"
      >
        <g
          v-if="asArray(dataFor(node).gains).length || historyLossCount(dataFor(node))"
          class="history-event-badge"
          :transform="`translate(${Math.max(node.parent.x + 8, (node.parent.x + node.x) / 2 - 14)} ${node.y - 18})`"
        >
          <rect
            width="44"
            height="16"
            rx="8"
          />
          <text
            x="22"
            y="11"
            text-anchor="middle"
          >
            <tspan class="history-gain-text">+{{ asArray(dataFor(node).gains).length }}</tspan>
            <tspan class="history-loss-text">−{{ historyLossCount(dataFor(node)) }}</tspan>
          </text>
        </g>
      </template>
      <g
        v-if="asArray(dataFor(layout.tree).gains).length || historyLossCount(dataFor(layout.tree))"
        class="history-event-badge root-event"
        :transform="`translate(${Math.max(4, layout.tree.x - 8)} ${layout.tree.y - 25})`"
      >
        <rect
          width="44"
          height="16"
          rx="8"
        />
        <text
          x="22"
          y="11"
          text-anchor="middle"
        >
          <tspan class="history-gain-text">+{{ asArray(dataFor(layout.tree).gains).length }}</tspan>
          <tspan class="history-loss-text">−{{ historyLossCount(dataFor(layout.tree)) }}</tspan>
        </text>
      </g>
      <g
        v-for="node in layout.nodes"
        :key="`history-node-${node.id}`"
        class="history-node-target"
        role="button"
        tabindex="0"
        :aria-label="`Inspect node ${publicNodeName(node.name)}`"
        @click="chooseNode(node)"
        @keydown.enter.prevent="chooseNode(node)"
        @keydown.space.prevent="chooseNode(node)"
      >
        <circle
          :cx="node.x"
          :cy="node.y"
          :r="String(selectedNodeName) === String(node.name) ? 7 : node.children.length ? 5 : 4"
          :class="[
            node.children.length ? 'history-node' : 'history-node leaf',
            { selected: String(selectedNodeName) === String(node.name) },
          ]"
        >
          <title>
            {{ publicNodeName(node.name) }} · {{ asArray(dataFor(node).spacers).length }} spacers
          </title>
        </circle>
      </g>
      <g
        v-for="leaf in layout.leaves"
        :key="`leaf-state-${leaf.id}`"
      >
        <text
          x="330"
          :y="leaf.y + 4"
          class="history-leaf-label"
        >
          <title>{{ leaf.name }}</title>
          {{ publicNodeName(leaf.name).slice(0, 24) }}
        </text>
        <rect
          v-for="(spacer, index) in spacerOrder"
          :key="spacer"
          :x="arrayStart + index * cellSize"
          :y="leaf.y - 9"
          :width="cellSize - 2"
          height="18"
          rx="2"
          :class="
            presentSpacers(leaf).has(spacer)
              ? ['history-spacer', 'present', spacerColorClass(spacer)]
              : ['history-spacer', 'absent']
          "
        >
          <title>
            Spacer {{ spacer }}: {{ presentSpacers(leaf).has(spacer) ? "present" : "absent" }}
          </title>
        </rect>
      </g>
    </svg>
  </div>
</template>
