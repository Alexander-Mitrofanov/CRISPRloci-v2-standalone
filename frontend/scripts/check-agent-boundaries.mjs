import { readFile, readdir } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const sourceRoot = path.join(projectRoot, "src");
const componentsRoot = path.join(sourceRoot, "components");
const resultsRoot = path.join(sourceRoot, "components", "results");
const maxAgentSourceLines = 500;

async function filesBelow(root, extensions) {
  const entries = await readdir(root, { withFileTypes: true });
  const nested = await Promise.all(
    entries.map(async (entry) => {
      const target = path.join(root, entry.name);
      if (entry.isDirectory()) return filesBelow(target, extensions);
      return extensions.has(path.extname(entry.name)) ? [target] : [];
    }),
  );
  return nested.flat();
}

const failures = [];
const resultFiles = await filesBelow(resultsRoot, new Set([".vue"]));
const forbiddenResultPatterns = [
  ["getValue(", "presentation code must consume canonical result fields"],
  ["orientation_evidence", "orientation aliases belong in the result normalizer"],
  ["detected_arrays", "detection aliases belong in the result normalizer"],
  ["delta_lnL", "likelihood aliases belong in the result normalizer"],
  ["delta_ln_likelihood", "likelihood aliases belong in the result normalizer"],
  ["delta_log_likelihood", "likelihood aliases belong in the result normalizer"],
  ["summary?.spacerplacer", "reconstruction aliases belong in the result normalizer"],
  ["summary.spacerplacer", "reconstruction aliases belong in the result normalizer"],
  ["group_id", "group identity aliases belong in the result normalizer"],
  ["row.name || row.group", "reconstruction group aliases belong in the result normalizer"],
  ["provenance.value.versions", "provenance aliases belong in the result normalizer"],
  ["request?.options", "request policy aliases belong in the result normalizer"],
];

const resultPresentationFiles = [path.join(sourceRoot, "utils", "results.js"), ...resultFiles];
for (const file of resultPresentationFiles) {
  const source = await readFile(file, "utf8");
  for (const [pattern, message] of forbiddenResultPatterns) {
    if (source.includes(pattern)) failures.push(`${path.relative(projectRoot, file)}: ${message}`);
  }
  if (source.includes("features/results/model/")) {
    failures.push(
      `${path.relative(projectRoot, file)}: import the stable features/results façade, not model internals`,
    );
  }
}

const sourceFiles = await filesBelow(sourceRoot, new Set([".css", ".js", ".vue"]));
const componentFiles = await filesBelow(componentsRoot, new Set([".vue"]));
for (const file of componentFiles) {
  const source = await readFile(file, "utf8");
  if (/from\s+["'][^"']*api\.js["']/u.test(source)) {
    failures.push(
      `${path.relative(projectRoot, file)}: presentation components must use an orchestration composable instead of importing the API client`,
    );
  }
  if (/\bfetch\s*\(/u.test(source)) {
    failures.push(
      `${path.relative(projectRoot, file)}: presentation components must use an orchestration composable instead of calling fetch`,
    );
  }
}

for (const file of sourceFiles) {
  const lines = (await readFile(file, "utf8")).split(/\r?\n/u);
  if (lines.length > maxAgentSourceLines) {
    failures.push(
      `${path.relative(projectRoot, file)}: ${lines.length}-line source exceeds the ${maxAgentSourceLines}-line agent-context limit; split it by ownership`,
    );
  }
  lines.forEach((line, index) => {
    if (line.length > 180) {
      failures.push(
        `${path.relative(projectRoot, file)}:${index + 1}: ${line.length}-character line exceeds the 180-character agent-diff limit`,
      );
    }
  });
}

if (failures.length) {
  console.error(failures.join("\n"));
  process.exitCode = 1;
} else {
  console.log(
    JSON.stringify({
      canonical_result_boundary: "enforced",
      presentation_transport_boundary: "enforced",
      checked_components: componentFiles.length,
      checked_result_components: resultFiles.length,
      checked_result_presentation_files: resultPresentationFiles.length,
      checked_source_files: sourceFiles.length,
      max_line_length: 180,
      max_source_lines: maxAgentSourceLines,
    }),
  );
}
