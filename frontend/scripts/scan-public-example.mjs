#!/usr/bin/env node
import { existsSync, readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative, resolve } from "node:path";

import {
  EXAMPLE_FASTA_PATH,
  EXAMPLE_RESULT_PATH,
  EXAMPLE_ASSETS,
  validateExampleInput,
} from "../src/example.js";
import { hasForbiddenPublicIdentity } from "./public-identity-policy.mjs";
import { verifiedIgvBundle } from "./igv-vendor.mjs";

const roots = ["public", "dist"].filter((root) => existsSync(root));
const sequenceExtension = /\.(?:fa|fasta|fna|ffn|fas)$/i;

function fail(message) {
  throw new Error(message);
}

function walk(root) {
  const files = [];
  for (const name of readdirSync(root)) {
    const path = join(root, name);
    const stat = statSync(path);
    if (stat.isDirectory()) files.push(...walk(path));
    else if (stat.isFile()) files.push(path);
  }
  return files;
}

for (const root of roots) {
  const vendorBundle = root === "dist" ? verifiedIgvBundle(root) : null;
  const fastaPath = resolve(root, EXAMPLE_FASTA_PATH);
  const resultPath = resolve(root, EXAMPLE_RESULT_PATH);
  if (!existsSync(fastaPath) || !existsSync(resultPath)) {
    fail(`bound example assets are missing from ${root}`);
  }
  const sequenceAssets = walk(root).filter((path) => sequenceExtension.test(path));
  const expectedSequences = new Set([
    fastaPath,
    ...Object.values(EXAMPLE_ASSETS)
      .filter((asset) => sequenceExtension.test(asset.input))
      .map((asset) => resolve(root, asset.input)),
  ]);
  if (
    sequenceAssets.length !== expectedSequences.size ||
    sequenceAssets.some((path) => !expectedSequences.has(resolve(path)))
  ) {
    fail(
      `unexpected public sequence asset set in ${root}: ${sequenceAssets.map((path) => relative(root, path)).join(", ")}`,
    );
  }
  const fasta = readFileSync(fastaPath, "utf8");
  const snapshot = JSON.parse(readFileSync(resultPath, "utf8"));
  await validateExampleInput(snapshot, fasta);
  const headers = fasta.split(/\r?\n/).filter((line) => line.startsWith(">"));
  const expectedHeaders = snapshot.example.records.map((record) => `>${record.record_id}`);
  if (
    headers.length !== expectedHeaders.length ||
    headers.some((header, index) => header !== expectedHeaders[index])
  ) {
    fail(`example FASTA headers are not fully masked in ${root}`);
  }
  if (hasForbiddenPublicIdentity(fasta) || hasForbiddenPublicIdentity(JSON.stringify(snapshot))) {
    fail(`source identity metadata is present in ${root}`);
  }
  const validatedModules = new Set();
  for (const asset of Object.values(EXAMPLE_ASSETS)) {
    const input = resolve(root, asset.input);
    const result = resolve(root, asset.result);
    const packet = JSON.parse(readFileSync(result, "utf8"));
    if (packet.job.mode !== asset.key.split("-")[0]) fail(`wrong example method: ${asset.key}`);
    await validateExampleInput(packet, readFileSync(input, "utf8"));
    validatedModules.add(input);
    validatedModules.add(result);
  }
  for (const file of walk(root)) {
    if (!/\.(?:html|js|json|txt|css)$/i.test(file)) continue;
    const text = readFileSync(file, "utf8");
    const name = relative(root, file);
    const compiledBundle =
      (root === "dist" && /^assets\/[^/]+\.js$/i.test(name)) || validatedModules.has(resolve(file));
    if (
      (resolve(file) !== vendorBundle && hasForbiddenPublicIdentity(text, { compiledBundle })) ||
      hasForbiddenPublicIdentity(name)
    ) {
      fail(`source identity metadata is present: ${relative(root, file)}`);
    }
  }
}

console.log(
  JSON.stringify({
    scanned_roots: roots,
    stored_input: "sha256-bound",
    precomputed_result: "validated",
  }),
);
