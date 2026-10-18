import { firstRecord } from "./values.js";

function asList(value) {
  return Array.isArray(value) ? value : value == null ? [] : [value];
}

export function normalizeProvenance(summary, context = {}) {
  const jobOptions = firstRecord(context.options);
  const requestOptions = firstRecord(context.request?.options);
  const source = firstRecord(summary?.provenance, context.provenance) || {};
  return {
    provenance: {
      ...source,
      parameters: firstRecord(source.parameters, jobOptions, requestOptions) || {},
      tool_versions: firstRecord(source.tool_versions, source.versions) || {},
    },
    warnings: [
      ...asList(summary?.warnings),
      ...asList(jobOptions?.warnings),
      ...asList(context.warnings),
    ].filter(Boolean),
  };
}
