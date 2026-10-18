<script setup>
import { computed, ref, watch } from "vue";

import { buildRecoveryUrl } from "../../features/recovery/index.js";
import AppIcon from "../common/AppIcon.vue";

const props = defineProps({ credential: { type: Object, required: true } });
const state = ref("idle");
const recoveryUrl = computed(() => buildRecoveryUrl(props.credential));
watch(
  () => `${props.credential.jobId}.${props.credential.accessToken}`,
  () => {
    state.value = "idle";
  },
);

async function copyLink() {
  const copiedUrl = recoveryUrl.value;
  try {
    await navigator.clipboard.writeText(copiedUrl);
    if (recoveryUrl.value !== copiedUrl) return;
    state.value = "copied";
  } catch {
    if (recoveryUrl.value !== copiedUrl) return;
    state.value = "failed";
  }
}
</script>

<template>
  <aside
    class="recovery-link"
    aria-label="Private job recovery link"
  >
    <AppIcon
      name="link"
      :size="20"
    />
    <div>
      <strong>Recovery link ready</strong>
      <p>Reopen this job until it expires. Anyone with the full link can access the job.</p>
    </div>
    <button
      type="button"
      @click="copyLink"
    >
      <AppIcon
        name="link"
        :size="16"
      />{{ state === "copied" ? "Link copied" : "Copy recovery link" }}
    </button>
    <small
      v-if="state !== 'idle'"
      class="recovery-link-status"
      role="status"
    >
      {{
        state === "copied"
          ? "Keep the link private; it is the key to this job."
          : "Clipboard access was unavailable. Copy the current address from the browser bar."
      }}
    </small>
  </aside>
</template>
