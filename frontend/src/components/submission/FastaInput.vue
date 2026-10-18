<script setup>
import { ref } from "vue";

import { readableBases } from "../../fasta.js";
import { readableBytes } from "../../utils/formatting.js";
import AppIcon from "../common/AppIcon.vue";

const props = defineProps({
  sequence: { type: String, required: true },
  filename: { type: String, required: true },
  inspection: { type: Object, required: true },
  molecule: { type: String, default: "DNA" },
  repeatInput: Boolean,
  spacerInput: Boolean,
  arrayInput: Boolean,
  familyInput: Boolean,
  locusInput: Boolean,
  loadingExample: Boolean,
  showExample: { type: Boolean, default: true },
  exampleDisabled: Boolean,
  maxRequestBytes: { type: Number, default: 0 },
});
const emit = defineEmits(["update:sequence", "update:filename", "load-example"]);
const fileInput = ref(null);
const fileError = ref("");

function updateSequence(value) {
  emit("update:sequence", value);
  emit("update:filename", "pasted-input.fasta");
  fileError.value = "";
}

async function onFile(event) {
  const [file] = event.target.files || [];
  event.target.value = "";
  if (!file) return;
  const margin = props.maxRequestBytes
    ? Math.min(65_536, Math.floor(props.maxRequestBytes * 0.05))
    : 0;
  const limit = props.maxRequestBytes ? Math.max(1, props.maxRequestBytes - margin) : 0;
  if (limit && file.size > limit) {
    fileError.value = `The selected file is ${readableBytes(file.size)}. It must be smaller than ${readableBytes(limit)} so the JSON request fits the service limit.`;
    return;
  }
  try {
    emit("update:sequence", await file.text());
    emit("update:filename", file.name);
    fileError.value = "";
  } catch {
    fileError.value = "The selected file could not be read as plain text.";
  }
}
</script>

<template>
  <div class="input-panel">
    <div
      id="example-entry"
      class="input-heading"
    >
      <div>
        <label for="fasta-input">{{
          arrayInput
            ? "Observed array JSON"
            : spacerInput
              ? "Spacer sequences"
              : repeatInput
                ? "Repeat sequences"
                : familyInput
                  ? "Homologous RNA-family sequences"
                  : locusInput
                    ? "Genomic neighborhoods"
                    : "Contigs or small genomes"
        }}</label>
        <p v-if="arrayInput">
          Upload arrays.json from detection. Each observed spacer counts toward the service record
          limit. Unknown orientation and uncertain boundaries are accepted.
        </p>
        <p v-else>{{ molecule }} FASTA or a compatible sequence JSON from a previous analysis.</p>
      </div>
      <button
        v-if="showExample"
        class="text-button example-button"
        type="button"
        :disabled="loadingExample || exampleDisabled"
        :title="
          exampleDisabled
            ? 'Finish or leave the current job before opening the example.'
            : undefined
        "
        @click="$emit('load-example')"
      >
        <AppIcon
          name="file"
          :size="16"
        />{{ loadingExample ? "Loading example…" : "Load example" }}
      </button>
    </div>
    <div class="upload-strip">
      <button
        class="upload-button"
        type="button"
        @click="fileInput?.click()"
      >
        <AppIcon
          name="upload"
          :size="18"
        />
        {{ arrayInput ? "Upload arrays.json" : "Upload FASTA or JSON" }}
      </button>
      <input
        ref="fileInput"
        type="file"
        tabindex="-1"
        accept=".fa,.fasta,.fna,.ffn,.fas,.txt,.json,text/plain,application/json"
        :aria-label="arrayInput ? 'Upload array JSON file' : 'Upload FASTA file'"
        @change="onFile"
      />
      <span class="filename">{{ sequence ? filename : "No file selected" }}</span>
      <span class="input-stats"
        ><b>{{ inspection.recordCount }}</b> records <i />
        <b>{{
          repeatInput
            ? `${inspection.baseCount.toLocaleString()} nt`
            : readableBases(inspection.baseCount)
        }}</b></span
      >
    </div>
    <textarea
      id="fasta-input"
      spellcheck="false"
      :value="sequence"
      :placeholder="
        arrayInput
          ? 'Paste detection arrays.json here'
          : repeatInput
            ? `>repeat_A\n${molecule === 'RNA' ? 'ACGU…' : 'ACGT…'}`
            : '>isolate_A\nACGT…\n>isolate_B\nACGT…'
      "
      aria-describedby="fasta-help fasta-errors"
      @input="updateSequence($event.target.value)"
    />
    <div
      id="fasta-help"
      class="input-foot"
    >
      <span v-if="arrayInput"
        >Observed spacers require unambiguous A/C/G/T; explicit deletions retain their slots. Input
        stays local until submission.</span
      >
      <span v-else
        >{{ familyInput ? "Unambiguous" : "IUPAC" }} {{ molecule }} accepted. Input stays local
        until submission.<template v-if="repeatInput">
          Ambiguous bases remain visible; affected analyses may be unsupported.</template
        ></span
      >
    </div>
    <div
      id="fasta-errors"
      class="field-errors"
      role="alert"
    >
      <p v-if="fileError">{{ fileError }}</p>
      <p
        v-for="item in sequence ? inspection.errors.slice(0, 3) : []"
        :key="item"
      >
        {{ item }}
      </p>
    </div>
  </div>
</template>
