<script setup>
import { computed, ref, useId } from "vue";
import { repeatStageText } from "../../../features/results/index.js";
import RepeatInstance from "./RepeatInstance.vue";

const props = defineProps({
  evidence: { type: Object, required: true },
  molecule: { type: String, default: "DNA" },
  withContext: Boolean,
  compact: Boolean,
});
const heading = useId();
const selectedId = ref(null);
const selected = computed(
  () =>
    props.evidence.instances.find((row) => row.id === selectedId.value) ||
    props.evidence.instances[0],
);
const orientation = (value) =>
  props.molecule === "RNA"
    ? "As supplied (RNA)"
    : value === "reverse"
      ? "Reverse-complement hypothesis"
      : value === "forward"
        ? "Forward hypothesis"
        : "Not reported";
const count = (value) => (value == null ? "Not available" : value.toLocaleString());
const metric = (value) => (value == null ? "Not available" : value.toFixed(4));
const interval = (row) =>
  row.start != null && row.end != null ? `[${row.start}, ${row.end})` : "No genomic interval";
</script>

<template>
  <section
    class="repeat-results"
    :aria-labelledby="heading"
  >
    <h3 :id="heading">Repeat evidence</h3>
    <p class="annotation-lead">
      {{
        withContext
          ? "Observed repeat instances, folded alone and in local array context."
          : "Isolated folding of supplied repeats. Array context was not supplied."
      }}
      Pair support is descriptive model evidence; it does not predict processing or immune activity.
    </p>
    <p
      v-if="evidence.status !== 'completed'"
      class="annotation-empty"
    >
      {{ repeatStageText(evidence) }}
    </p>
    <div class="annotation-counts">
      <div>
        <strong>{{ count(evidence.counts.sources) }}</strong
        ><span>Source records</span>
      </div>
      <div v-if="withContext">
        <strong>{{ count(evidence.counts.arrays) }}</strong
        ><span>Observed arrays</span>
      </div>
      <div>
        <strong>{{ count(evidence.counts.instances) }}</strong
        ><span>Repeat / orientation records</span>
      </div>
    </div>
    <p class="repeat-note">
      {{
        molecule === "RNA"
          ? "RNA is analyzed as the supplied transcribed 5′→3′ sequence."
          : "Forward and reverse-complement DNA hypotheses are separate alternatives, not a predicted transcription direction or independent repeats."
      }}
    </p>
    <p
      v-if="evidence.model.engine"
      class="repeat-note"
    >
      {{ evidence.model.engine }} {{ evidence.model.version }} · {{ evidence.model.parameters }} ·
      {{ evidence.model.temperature ?? "temperature not reported" }} °C. Temperature is a
      computational setting.
    </p>
    <dl class="repeat-classifications">
      <div
        v-for="(stage, name) in evidence.classification"
        :key="name"
      >
        <dt>
          {{
            name === "strand"
              ? "Strand inference"
              : name === "subtype"
                ? "Subtype association"
                : "Repeat family"
          }}
        </dt>
        <dd>{{ repeatStageText(stage) }}</dd>
      </div>
    </dl>
    <template v-if="!compact">
      <p
        v-if="evidence.interpretation"
        class="repeat-note"
      >
        {{ evidence.interpretation }}
      </p>
      <p
        v-if="evidence.truncated"
        class="annotation-notice"
      >
        This preview is limited to {{ evidence.instances.length }} repeat / orientation records.
        Download complete evidence from Files &amp; methods.
      </p>
      <p
        v-if="evidence.counts.instances === 0"
        class="annotation-empty"
      >
        No repeat instances were available for this analysis.
      </p>
      <p
        v-else-if="!evidence.instances.length"
        class="annotation-empty"
      >
        No instance preview was reported. Check the stage status and Files &amp; methods.
      </p>
      <template v-else>
        <p class="repeat-note">
          Choose a repeat and orientation. Source intervals use zero-based, half-open coordinates on
          the submitted forward sequence.
        </p>
        <div
          class="annotation-scroll"
          tabindex="0"
          role="region"
          aria-label="Repeat instances table"
        >
          <table class="annotation-table repeat-instance-table">
            <thead>
              <tr>
                <th scope="col">Inspect</th>
                <th scope="col">Source / array / repeat</th>
                <th scope="col">Orientation hypothesis</th>
                <th scope="col">Source interval</th>
                <th scope="col">Isolated support</th>
                <th scope="col">Context support</th>
                <th scope="col">Change</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="row in evidence.instances"
                :key="row.id"
                :class="{ 'repeat-selected': selected?.id === row.id }"
              >
                <td>
                  <input
                    type="radio"
                    :name="`${heading}-instance`"
                    :checked="selected?.id === row.id"
                    :aria-label="`Inspect ${row.source_id}, ${row.array_id || 'supplied repeat'}, ${row.repeat_index ?? ''}, ${orientation(row.orientation)}`"
                    @change="selectedId = row.id"
                  />
                </td>
                <td>
                  {{ row.source_id ?? "Not reported"
                  }}<small v-if="row.array_id"
                    >{{ row.array_id }} / {{ row.repeat_index ?? "Not reported" }}</small
                  >
                </td>
                <td>{{ orientation(row.orientation) }}</td>
                <td>{{ interval(row) }}</td>
                <td>{{ metric(row.comparison.isolated_pair_support) }}</td>
                <td>{{ metric(row.comparison.context_pair_support) }}</td>
                <td>{{ metric(row.comparison.support_change) }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <RepeatInstance
          v-if="selected"
          :key="selected.id"
          :instance="selected"
        />
      </template>
    </template>
  </section>
</template>
