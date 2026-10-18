/**
 * Stable public facade for the bundled-example feature.
 *
 * Components, tests, and build-time scanners import this file so internal
 * validators can evolve without creating cross-feature edit collisions.
 */
export {
  EXAMPLE_FASTA_PATH,
  EXAMPLE_RESULT_PATH,
  EXAMPLE_SCHEMA_VERSION,
} from "./features/example/contract.js";
export { validateExampleSnapshot } from "./features/example/snapshot.js";
export { validateExampleInput } from "./features/example/input.js";
export {
  EXAMPLE_ASSETS,
  exampleAssetsFor,
  exampleSubmissionMatches,
} from "./features/example/catalogue.js";
export { readExampleResponse } from "./features/example/response.js";
