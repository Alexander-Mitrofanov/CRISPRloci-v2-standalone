export function asRecord(value) {
  return value && typeof value === "object" && !Array.isArray(value) ? value : null;
}

export function asList(value) {
  return Array.isArray(value) ? value : [];
}

export function firstDefined(source, ...keys) {
  const record = asRecord(source);
  if (!record) return undefined;
  for (const key of keys) {
    if (record[key] !== undefined) return record[key];
  }
  return undefined;
}

export function firstRecord(...values) {
  return values.map(asRecord).find(Boolean) || null;
}
