<script setup>
import { toRefs } from "vue";

import { useArtifactDownloads } from "../../../features/results/index.js";
import { readableBytes } from "../../../utils/formatting.js";
import AppIcon from "../../common/AppIcon.vue";

const props = defineProps({
  job: { type: Object, required: true },
  credential: { type: Object, required: true },
  maxArchiveBytes: { type: Number, default: 0 },
});
const { job, credential } = toRefs(props);
const { bundle, individual, downloading, error, download } = useArtifactDownloads(job, credential);
</script>

<template>
  <section
    class="result-section downloads"
    aria-labelledby="downloads-heading"
  >
    <div class="result-heading">
      <div>
        <h3 id="downloads-heading">Downloads</h3>
      </div>
    </div>
    <p v-if="individual.some((artifact) => (artifact.filename || artifact.name) === 'result.json')">
      result.json is the compact scientific result for this analysis. Detailed evidence and
      provenance are available in the complete bundle.
    </p>
    <p
      v-if="
        individual.some((artifact) =>
          ['repeats.json', 'spacers.json', 'input-sequences.json'].includes(artifact.name),
        )
      "
    >
      Continue an analysis: upload repeats.json to repeat structure, reference mapping or subtype
      analysis; spacers.json to viral search or spacer-group association; arrays.json to observed
      array comparison; or input-sequences.json to a mode accepting the same input type. Full tool
      details are in the bundle.
    </p>
    <p class="download-memory-note">
      <AppIcon
        name="info"
        :size="16"
      />
      Authenticated files are buffered in this browser tab before saving<span v-if="maxArchiveBytes"
        >; the archive cap is {{ readableBytes(maxArchiveBytes) }}</span
      >.
    </p>
    <button
      v-if="bundle"
      class="bundle-button"
      type="button"
      :disabled="Boolean(downloading)"
      @click="download()"
    >
      <span><AppIcon name="download" /><i>ZIP</i></span
      ><span
        ><strong>{{
          downloading === "bundle" ? "Preparing download…" : "Download complete result bundle"
        }}</strong
        ><small>Tables · trees · alignments · warnings · provenance</small></span
      >
    </button>
    <div
      v-else
      class="archive-unavailable"
      role="note"
    >
      <AppIcon
        name="warning"
        :size="17"
      />
      <p><strong>Complete ZIP not available.</strong> Use any individual outputs listed below.</p>
    </div>
    <div
      v-if="individual.length"
      class="artifact-grid"
    >
      <button
        v-for="(artifact, index) in individual"
        :key="artifact.artifact_id || artifact.id || index"
        type="button"
        :disabled="Boolean(downloading)"
        @click="download(artifact)"
      >
        <AppIcon name="file" /><span
          ><strong>{{
            artifact.label || artifact.filename || artifact.name || `Artifact ${index + 1}`
          }}</strong
          ><small>{{
            (artifact.filename || artifact.name) === "result.json"
              ? "Compact scientific result · JSON"
              : artifact.media_type || artifact.kind || "Result file"
          }}</small></span
        ><AppIcon
          name="download"
          :size="17"
        />
      </button>
    </div>
    <p
      v-if="error"
      class="download-error"
      role="alert"
    >
      {{ error }}
    </p>
  </section>
</template>
