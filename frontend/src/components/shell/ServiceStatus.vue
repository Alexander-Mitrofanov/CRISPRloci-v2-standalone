<script setup>
defineProps({
  service: { type: Object, required: true },
  analysisContext: { type: Boolean, default: false },
});
</script>

<template>
  <div
    :class="['service-status', `service-${service.state}`]"
    role="status"
    :title="analysisContext ? `Analysis API: ${service.message}` : service.message"
  >
    <span
      class="status-pulse"
      aria-hidden="true"
    />
    <span>{{
      {
        online: analysisContext ? "Analysis ready" : "Service ready",
        offline: analysisContext ? "Analysis unavailable" : "Service unavailable",
        checking: analysisContext ? "Analysis connecting…" : "Connecting…",
      }[service.state] || "Connecting…"
    }}</span>
  </div>
</template>
