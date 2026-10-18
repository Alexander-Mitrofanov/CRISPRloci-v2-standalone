<script setup>
import { computed } from "vue";
import { displayToolDefault, toolOptionError } from "../../toolOptions.js";

const props = defineProps({
  toolId: { type: String, required: true },
  field: { type: Object, required: true },
  modelValue: { type: [String, Number, Boolean], default: undefined },
  validationError: { type: String, default: "" },
});
const emit = defineEmits(["update:modelValue"]);
const id = computed(() => `tool-${props.toolId}-${props.field.flag.replace(/^-+/, "")}`);
const value = computed(() => props.modelValue ?? props.field.default);
const changed = computed(
  () => props.modelValue !== undefined && props.modelValue !== props.field.default,
);
const error = computed(
  () => props.validationError || toolOptionError(props.field, props.modelValue),
);

function updateNumber(event) {
  const raw = event.target.value;
  emit(
    "update:modelValue",
    raw === "" ? (props.field.default === null ? undefined : "") : Number(raw),
  );
}

function updateChoice(event) {
  emit(
    "update:modelValue",
    props.field.choices.find((choice) => String(choice) === event.target.value),
  );
}
</script>

<template>
  <div
    class="tool-option-field"
    :class="{ 'tool-option-changed': changed }"
  >
    <div class="tool-option-description">
      <label :for="id"
        ><strong>{{ field.label }}</strong> <code>{{ field.flag }}</code></label
      >
      <p :id="`${id}-help`">{{ field.description }}</p>
      <small :id="`${id}-default`">
        Default: {{ displayToolDefault(field.default) }}
        <template v-if="field.min !== undefined && field.max !== undefined">
          · Service range: {{ field.min }}–{{ field.max }}
        </template>
      </small>
    </div>
    <div class="tool-option-control">
      <select
        v-if="field.type === 'boolean'"
        :id="id"
        :value="String(value)"
        :aria-describedby="`${id}-help ${id}-default`"
        @change="emit('update:modelValue', $event.target.value === 'true')"
      >
        <option value="true">On</option>
        <option value="false">Off</option>
      </select>
      <select
        v-else-if="field.type === 'enum'"
        :id="id"
        :value="value"
        :aria-describedby="`${id}-help ${id}-default`"
        @change="updateChoice"
      >
        <option
          v-for="choice in field.choices"
          :key="choice"
          :value="choice"
        >
          {{ choice }}
        </option>
      </select>
      <input
        v-else
        :id="id"
        type="number"
        :min="field.min"
        :max="field.max"
        :step="field.type === 'integer' ? 1 : 'any'"
        :placeholder="field.default === null ? 'Automatic' : undefined"
        :value="value"
        :aria-describedby="`${id}-help ${id}-default${error ? ` ${id}-error` : ''}`"
        :aria-invalid="error ? 'true' : undefined"
        @input="updateNumber"
      />
      <button
        v-if="changed"
        type="button"
        class="tool-option-reset"
        :aria-label="`Reset ${field.flag} to default`"
        @click="emit('update:modelValue', undefined)"
      >
        Reset
      </button>
      <small
        v-if="error"
        :id="`${id}-error`"
        class="tool-option-error"
        role="alert"
        >{{ error }}</small
      >
    </div>
  </div>
</template>
