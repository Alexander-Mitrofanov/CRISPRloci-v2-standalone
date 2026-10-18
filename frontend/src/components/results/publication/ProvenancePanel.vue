<script setup>
import { computed } from "vue";

import { asArray } from "../../../utils/formatting.js";
import AppIcon from "../../common/AppIcon.vue";

const props = defineProps({
  summary: { type: Object, required: true },
  showWarnings: { type: Boolean, default: true },
});
const warnings = computed(() => asArray(props.summary.warnings));
const provenance = computed(() => props.summary.provenance);
const versions = computed(() => provenance.value.tool_versions);
const parameters = computed(() => provenance.value.parameters);
const displayValue = (value) =>
  typeof value === "boolean"
    ? value
      ? "enabled"
      : "disabled"
    : Array.isArray(value)
      ? value.join(", ")
      : value && typeof value === "object"
        ? JSON.stringify(value)
        : String(value);
</script>

<template>
  <section
    class="result-section provenance-section"
    aria-labelledby="provenance-heading"
  >
    <div class="result-heading">
      <div>
        <h3 id="provenance-heading">Methods & provenance</h3>
      </div>
    </div>
    <ul
      v-if="showWarnings && warnings.length"
      class="warning-list"
    >
      <li
        v-for="(warning, index) in warnings"
        :key="`${warning?.code || 'warning'}:${index}`"
      >
        <AppIcon
          name="warning"
          :size="18"
        /><span
          ><strong>{{ warning.title || warning.code || "Analysis warning" }}</strong
          >{{ warning.message || String(warning) }}</span
        >
      </li>
    </ul>
    <p
      v-else-if="showWarnings"
      class="no-warnings"
    >
      <AppIcon
        name="check"
        :size="17"
      />
      No workflow warnings were reported.
    </p>
    <div
      v-if="provenance.repositories?.length"
      class="repository-provenance"
    >
      <h4>Exact tool revisions</h4>
      <ul>
        <li
          v-for="entry in provenance.repositories"
          :key="entry.display_name"
        >
          {{ entry.display_name }}:
          <a
            :href="`${entry.repository.replace(/\.git$/, '')}/commit/${entry.commit}`"
            :aria-label="`${entry.display_name} revision ${entry.commit.slice(0, 12)}`"
            target="_blank"
            rel="noopener noreferrer"
            >{{ entry.commit.slice(0, 12) }}</a
          >
        </li>
      </ul>
    </div>
    <div class="provenance-grid">
      <div>
        <h4>Tool versions</h4>
        <dl v-if="Object.keys(versions).length">
          <div
            v-for="(value, key) in versions"
            :key="key"
          >
            <dt>{{ key.replaceAll("_", " ") }}</dt>
            <dd>{{ displayValue(value) }}</dd>
          </div>
        </dl>
        <p
          v-else
          class="muted"
        >
          See the provenance manifest in the result bundle.
        </p>
      </div>
      <div>
        <h4>Recorded policy</h4>
        <dl v-if="Object.keys(parameters).length">
          <div
            v-for="(value, key) in parameters"
            :key="key"
          >
            <dt>{{ key.replaceAll("_", " ") }}</dt>
            <dd>{{ displayValue(value) }}</dd>
          </div>
        </dl>
        <p
          v-else
          class="muted"
        >
          See the provenance manifest in the result bundle.
        </p>
      </div>
    </div>
  </section>
</template>
