export function formatNumber(value, digits = 0) {
  if (value == null || value === "") return "—";
  const number = Number(value);
  if (!Number.isFinite(number)) return "—";
  return number.toLocaleString(undefined, { maximumFractionDigits: digits });
}

export function formatDate(value) {
  const parsed = value ? new Date(value) : null;
  return parsed && !Number.isNaN(parsed.getTime()) ? parsed.toLocaleString() : "—";
}

export function formatDuration(value) {
  const seconds = Number(value);
  if (!Number.isFinite(seconds)) return "—";
  if (seconds < 60) return `${seconds.toFixed(seconds < 10 ? 1 : 0)} s`;
  return `${Math.floor(seconds / 60)} min ${Math.round(seconds % 60)} s`;
}

export function readableBytes(value) {
  const bytes = Number(value);
  if (!Number.isFinite(bytes) || bytes <= 0) return "the configured service limit";
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KiB`;
  return `${(bytes / 1024 ** 2).toFixed(bytes < 10 * 1024 ** 2 ? 1 : 0)} MiB`;
}

export function asArray(value) {
  if (Array.isArray(value)) return value;
  return value == null ? [] : [value];
}

export function getValue(source, ...keys) {
  for (const key of keys) {
    if (source && source[key] != null) return source[key];
  }
  return null;
}

export function finiteMetric(value) {
  if (value == null || value === "") return null;
  const number = Number(value);
  return Number.isFinite(number) ? number : null;
}

export function signedNumber(value, digits = 2) {
  const number = finiteMetric(value);
  if (number == null) return "not reported";
  return `${number > 0 ? "+" : ""}${number.toFixed(digits)}`;
}

export function downloadName(value, fallback) {
  const candidate = String(value || fallback)
    .split(/[\\/]/)
    .pop();
  return candidate.replace(/[^A-Za-z0-9._-]+/g, "_") || fallback;
}
