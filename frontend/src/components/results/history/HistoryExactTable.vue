<script setup>
import { computed } from "vue";

import { historyValueList, publicNodeName, specialEventSummary } from "../../../utils/history.js";

const props = defineProps({
  layout: { type: Object, required: true },
  nodeData: { type: Array, required: true },
  group: { type: String, required: true },
});
const dataByName = computed(
  () => new Map(props.nodeData.map((node) => [String(node?.name), node])),
);
const leafNames = computed(() => new Set(props.layout.leaves.map((leaf) => String(leaf.name))));
const dataFor = (node) => dataByName.value.get(String(node?.name)) || {};
const nodeType = (node) =>
  leafNames.value.has(String(node.name))
    ? "Observed leaf"
    : node === props.layout.tree
      ? "Inferred root"
      : "Inferred ancestor";
</script>

<template>
  <details class="history-data-table">
    <summary>Exact node and branch data</summary>
    <div
      class="table-wrap"
      role="region"
      tabindex="0"
      :aria-label="`Scrollable exact ancestral reconstruction data for ${group}`"
    >
      <table>
        <thead>
          <tr>
            <th>Node</th>
            <th>Parent</th>
            <th>Type</th>
            <th>Branch length</th>
            <th>Ordered spacer IDs</th>
            <th>Gains</th>
            <th>Loss blocks</th>
            <th>Special-event candidates</th>
          </tr>
        </thead>
        <tbody>
          <tr
            v-for="node in layout.nodes"
            :key="`exact-${node.id}`"
          >
            <td>
              <strong>{{ publicNodeName(node.name) }}</strong>
            </td>
            <td>{{ node.parent ? publicNodeName(node.parent.name) : "—" }}</td>
            <td>{{ nodeType(node) }}</td>
            <td>
              {{
                node.parent
                  ? Number(node.length).toLocaleString(undefined, { maximumFractionDigits: 6 })
                  : "root"
              }}
            </td>
            <td>{{ historyValueList(dataFor(node).spacers) }}</td>
            <td>{{ historyValueList(dataFor(node).gains) }}</td>
            <td>{{ historyValueList(dataFor(node).loss_blocks) }}</td>
            <td>{{ specialEventSummary(dataFor(node)) }}</td>
          </tr>
        </tbody>
      </table>
    </div>
  </details>
</template>
