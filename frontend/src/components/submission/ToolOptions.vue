<script setup>
import { computed } from "vue";
import { toolOptionErrors, toolsForMode } from "../../toolOptions.js";
import ToolOptionField from "./ToolOptionField.vue";

const props = defineProps({
  mode: { type: String, required: true },
  modelValue: { type: Object, default: () => ({}) },
});
const emit = defineEmits(["update:modelValue"]);
const tools = computed(() => toolsForMode(props.mode, props.modelValue));
const errors = computed(() => toolOptionErrors(props.mode, props.modelValue));

function changedCount(tool) {
  return tool.flags.filter((field) => {
    const value = props.modelValue[tool.id]?.[field.flag];
    return value !== undefined && value !== field.default;
  }).length;
}

function update(tool, field, value) {
  const toolValues = { ...props.modelValue[tool.id] };
  if (value === undefined || value === field.default) delete toolValues[field.flag];
  else toolValues[field.flag] = value;
  const next = { ...props.modelValue };
  if (Object.keys(toolValues).length) next[tool.id] = toolValues;
  else delete next[tool.id];
  emit("update:modelValue", next);
}

function reset(tool) {
  const next = { ...props.modelValue };
  delete next[tool.id];
  emit("update:modelValue", next);
}
</script>

<template>
  <div class="tool-options">
    <p class="tool-options-intro">
      Tune the tools used by this analysis. Each flag shows its effect and default.
    </p>
    <details
      v-for="tool in tools"
      :key="tool.id"
      class="tool-options-tool"
    >
      <summary>
        <span
          ><strong>{{ tool.name }}</strong
          ><small>{{
            tool.flags.length
              ? `${tool.flags.length} adjustable ${tool.flags.length === 1 ? "flag" : "flags"}`
              : "Workflow-managed flags"
          }}</small></span
        >
        <span
          v-if="changedCount(tool)"
          class="tool-options-count"
          >{{ changedCount(tool) }} changed</span
        >
      </summary>
      <div class="tool-options-body">
        <div class="tool-options-about">
          <p>{{ tool.description }}</p>
          <a
            :href="tool.documentation"
            target="_blank"
            rel="noopener noreferrer"
            >CLI documentation<span class="sr-only"> for {{ tool.name }}</span></a
          >
        </div>
        <dl
          v-if="tool.cli_modes?.length"
          class="tool-cli-modes"
        >
          <template
            v-for="cliMode in tool.cli_modes"
            :key="cliMode.name"
          >
            <dt>{{ cliMode.name }}</dt>
            <dd>{{ cliMode.description }}</dd>
          </template>
        </dl>
        <button
          v-if="changedCount(tool)"
          type="button"
          class="tool-reset-all"
          @click="reset(tool)"
        >
          Reset {{ tool.name }} flags
        </button>
        <ToolOptionField
          v-for="field in tool.flags"
          :key="field.flag"
          :tool-id="tool.id"
          :field="field"
          :model-value="modelValue[tool.id]?.[field.flag]"
          :validation-error="errors[tool.id]?.[field.flag]"
          @update:model-value="update(tool, field, $event)"
        />
        <details
          v-if="tool.managed?.length"
          class="tool-managed"
        >
          <summary>Managed CLI arguments</summary>
          <p>These arguments are set by the selected workflow and service.</p>
          <dl>
            <template
              v-for="entry in tool.managed"
              :key="String(entry.flags)"
            >
              <dt>
                <code>{{ Array.isArray(entry.flags) ? entry.flags.join(", ") : entry.flags }}</code>
              </dt>
              <dd>{{ entry.description }}</dd>
            </template>
          </dl>
        </details>
      </div>
    </details>
  </div>
</template>
