<script setup>
import { catalogCell, catalogEvidence } from "../../features/catalog/index.js";

defineProps({
  rows: { type: Array, required: true },
  columns: { type: Array, required: true },
  label: { type: String, required: true },
});
defineEmits(["array", "genome"]);
</script>

<template>
  <div
    class="catalog-table-scroll"
    role="region"
    :aria-label="label"
    tabindex="0"
  >
    <table class="catalog-table">
      <thead>
        <tr>
          <th
            v-for="column in columns"
            :key="column[0]"
            scope="col"
          >
            {{ column[1] }}
          </th>
        </tr>
      </thead>
      <tbody>
        <tr
          v-for="row in rows"
          :key="row.id ?? row.accession"
        >
          <td
            v-for="[key, , kind] in columns"
            :key="key"
          >
            <button
              v-if="kind === 'array'"
              class="catalog-record-link"
              type="button"
              :aria-label="`Open array ${row.call_id || row.array_id}`"
              @click="$emit('array', row.array_id ?? row.id)"
            >
              {{ catalogCell(row, key) }}
            </button>
            <button
              v-else-if="kind === 'genome'"
              class="catalog-record-link"
              type="button"
              :aria-label="`Browse genome ${row.accession}`"
              @click="$emit('genome', row.accession)"
            >
              {{ catalogCell(row, key) }}
            </button>
            <details
              v-else-if="kind === 'evidence'"
              class="catalog-call-evidence"
            >
              <summary :aria-label="`Evidence for ${row.call_id}`">View evidence</summary>
              <dl>
                <div
                  v-for="field in catalogEvidence(row, key)"
                  :key="field.label"
                >
                  <dt>{{ field.label }}</dt>
                  <dd>{{ field.value }}</dd>
                </div>
              </dl>
            </details>
            <code
              v-else-if="kind === 'sequence'"
              class="catalog-sequence"
              >{{ catalogCell(row, key) }}</code
            >
            <span v-else>{{ catalogCell(row, key) }}</span>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</template>
