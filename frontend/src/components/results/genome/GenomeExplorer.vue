<script setup>
import { ref, watch } from "vue";
import { useGenomeData } from "../../../features/genome/index.js";
import GenomeBrowser from "./GenomeBrowser.vue";

const props = defineProps({
  jobId: { type: String, required: true },
  mode: { type: String, required: true },
  accessToken: { type: String, required: true },
  artifacts: { type: Array, required: true },
  active: { type: Boolean, required: true },
});
const emit = defineEmits(["open-evidence", "ready"]);
const { model, loading, error, load } = useGenomeData(props);
const rendererOpened = ref(false);
watch([model, () => props.active], () => {
  if (!model.value) {
    rendererOpened.value = false;
    emit("ready", false);
  } else if (props.active) rendererOpened.value = true;
});
</script>

<template>
  <section aria-label="Genome explorer">
    <h3>Genome &amp; locus</h3>
    <p>Explore reported features together on the submitted DNA, then zoom to individual bases.</p>
    <p
      v-if="loading"
      role="status"
    >
      Loading complete genome evidence…
    </p>
    <div
      v-if="error"
      role="alert"
    >
      <p>{{ error }}</p>
      <button
        type="button"
        :disabled="loading"
        @click="load"
      >
        Retry genome loading
      </button>
    </div>
    <GenomeBrowser
      v-if="model && rendererOpened"
      :key="jobId"
      :sources="model.sources"
      :tracks="model.tracks"
      :features="model.features"
      :notes="model.notes"
      :active="active"
      @open-evidence="emit('open-evidence')"
      @ready="emit('ready', $event)"
    />
  </section>
</template>
