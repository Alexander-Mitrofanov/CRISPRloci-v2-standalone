<script setup>
defineProps({ modelValue: { type: String, default: "" } });
defineEmits(["update:modelValue"]);
</script>
<template>
  <fieldset class="association-grouping">
    <legend>Choose how spacers are grouped</legend>
    <p>
      Grouping defines the association hypothesis. A sequence record or contig is not necessarily a
      whole host. No grouping is selected automatically.
    </p>
    <label for="association-grouping">Spacer grouping</label>
    <select
      id="association-grouping"
      :value="modelValue"
      @change="$emit('update:modelValue', $event.target.value)"
    >
      <option
        value=""
        disabled
      >
        Select a grouping
      </option>
      <option value="single_set">One group — I confirm all spacers belong together</option>
      <option value="per_record">Each spacer record separately (up to 100 groups)</option>
      <option value="per_source">By source record in a spacers.json export</option>
    </select>
    <p v-if="modelValue === 'per_source'">
      Paste or upload a complete spacers.json sequence-set export. Each occurrence must have source
      metadata. Sources remain separate; submitted origins are annotations, not independently
      verified host identities.
    </p>
    <p>
      Up to 1,000 occurrences, 200 bases each, and 100 groups. Duplicate occurrences are retained.
      Unambiguous spacers shorter than 27 bases and ambiguous bases are unsupported.
    </p>
  </fieldset>
</template>
<style scoped>
.association-grouping {
  margin: 0 0 1.5rem;
  padding: 1.2rem 1.4rem;
  min-width: 0;
  border: 1px solid var(--line-dark);
  border-radius: 0.6rem;
  background: #f0f8f7;
}
.association-grouping legend {
  padding: 0 0.4rem;
  color: var(--ink);
  font-weight: 600;
  font-size: 1rem;
}
.association-grouping p {
  max-width: 80ch;
  color: var(--ink-soft);
  font-size: 0.9rem;
  line-height: 1.65;
  margin: 0.5rem 0 1rem;
}
.association-grouping p:last-child {
  margin-bottom: 0;
}
.association-grouping label {
  display: block;
  margin-bottom: 0.5rem;
  font-weight: 600;
}
.association-grouping select {
  width: 100%;
  max-width: 100%;
  white-space: normal;
  padding: 0.7rem 0.8rem;
  min-height: 44px;
  border: 1px solid var(--line-dark);
  border-radius: 0.35rem;
  background: white;
  color: var(--ink);
  font: inherit;
}
</style>
