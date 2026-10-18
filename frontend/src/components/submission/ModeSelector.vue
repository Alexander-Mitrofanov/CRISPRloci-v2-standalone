<script setup>
import { computed } from "vue";
import { ANALYSIS_MODES, analysisModeAvailable } from "../../science.js";
import { exampleAssetsFor } from "../../example.js";

const props = defineProps({ service: { type: Object, default: () => ({}) } });
const modes = computed(() =>
  ANALYSIS_MODES.filter(
    (item) =>
      !["protospacer", "viral_search"].includes(item.id) &&
      (analysisModeAvailable(item.id, props.service) || exampleAssetsFor(item.id)),
  ),
);

const spacerMode = computed(() =>
  ["protospacer", "viral_search"].find(
    (mode) => analysisModeAvailable(mode, props.service) || exampleAssetsFor(mode),
  ),
);

const emit = defineEmits(["select", "database"]);

function openDatabase(event) {
  if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey || event.button !== 0)
    return;
  event.preventDefault();
  emit("database");
}

function choose(event, id) {
  if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey || event.button !== 0)
    return;
  event.preventDefault();
  emit("select", id);
}
</script>

<template>
  <section
    class="method-picker"
    aria-labelledby="methods-heading"
  >
    <h1
      id="methods-heading"
      tabindex="-1"
    >
      Choose a method
    </h1>
    <div class="method-grid">
      <a
        class="method-link database-entry"
        href="?view=database"
        aria-label="Search the database"
        aria-describedby="database-entry-description"
        @click="openDatabase"
      >
        <strong>Search the database</strong>
        <span id="database-entry-description">Genomes, repeats, spacers &amp; Cas genes</span>
      </a>
      <a
        v-if="spacerMode"
        aria-label="Spacer searches"
        class="method-link"
        :href="`?method=${spacerMode}`"
        @click="choose($event, spacerMode)"
        ><strong>Spacer searches</strong
        ><span>Find spacers in DNA or candidate viruses for your spacers</span></a
      >
      <a
        v-for="item in modes"
        :key="item.id"
        class="method-link"
        :href="`?method=${item.id}`"
        :aria-label="item.title"
        :aria-describedby="`method-${item.id}-description`"
        @click="choose($event, item.id)"
      >
        <strong>{{ item.title }}</strong>
        <span :id="`method-${item.id}-description`">{{ item.short }}</span>
      </a>
    </div>
  </section>
</template>
