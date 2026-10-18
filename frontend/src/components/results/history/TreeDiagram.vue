<script setup>
import { computed } from "vue";

import { layoutNewick } from "../../../utils/newick.js";

const props = defineProps({
  newick: { type: String, default: "" },
  group: { type: String, default: "group" },
  selectedNode: { type: String, default: "" },
  reportedByDefault: Boolean,
  statusLabel: { type: String, default: "" },
});
const emit = defineEmits(["select-node"]);
const provisional = computed(() => layoutNewick(props.newick, 100));
const height = computed(() => Math.max(190, (provisional.value?.leaves.length || 1) * 58 + 58));
const layout = computed(() => layoutNewick(props.newick, height.value));
const internalEdges = computed(
  () => layout.value?.nodes.filter((node) => node.children.length) || [],
);
const displayLabel = computed(
  () =>
    props.statusLabel ||
    (props.reportedByDefault
      ? "Reported input-order model tree"
      : "Displayed hypothesis model tree"),
);
</script>

<template>
  <div
    v-if="layout"
    class="tree-graphic"
  >
    <div class="graphic-label">
      <span>{{ displayLabel }}</span
      ><small>Branch lengths scaled when available</small>
    </div>
    <svg
      :viewBox="`0 0 720 ${height}`"
      role="img"
      :aria-label="`${displayLabel} for ${group} with ${layout.leaves.length} leaves`"
    >
      <desc>Rooted model tree for {{ group }}. Exact Newick: {{ newick }}</desc>
      <g
        v-for="node in internalEdges"
        :key="`edges-${node.id}`"
      >
        <line
          :x1="node.x"
          :x2="node.x"
          :y1="Math.min(...node.children.map((child) => child.y))"
          :y2="Math.max(...node.children.map((child) => child.y))"
          class="tree-line"
        />
        <line
          v-for="child in node.children"
          :key="`edge-${child.id}`"
          :x1="node.x"
          :x2="child.x"
          :y1="child.y"
          :y2="child.y"
          class="tree-line"
        />
      </g>
      <g
        v-for="node in layout.nodes"
        :key="`node-${node.id}`"
        class="history-node-target"
        role="button"
        tabindex="0"
        :aria-label="`Inspect node ${node.name || `node ${node.id}`}`"
        @click="emit('select-node', node.name)"
        @keydown.enter.prevent="emit('select-node', node.name)"
      >
        <circle
          :cx="node.x"
          :cy="node.y"
          :r="String(node.name) === selectedNode ? 7 : node.children.length ? 4 : 5"
          :class="[
            node.children.length ? 'tree-node' : 'tree-leaf-node',
            { selected: String(node.name) === selectedNode },
          ]"
        />
        <text
          v-if="!node.children.length"
          :x="node.x + 11"
          :y="node.y + 4"
          class="tree-leaf-label"
        >
          {{ node.name.length > 30 ? `${node.name.slice(0, 28)}…` : node.name }}
        </text>
      </g>
    </svg>
  </div>
</template>
