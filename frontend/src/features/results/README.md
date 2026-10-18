# Public result model

`index.js` is the stable boundary between backend payloads and Vue result
components. `normalizePublicResult` accepts the supported canonical and legacy
field aliases once, then supplies one canonical object graph to every result
slice.

- `model/detection.js` owns detector arrays and category counts.
- `model/adapter.js` owns preflight and evolutionary-group membership.
- `model/orientation.js` owns hypothesis comparisons and structured histories.
- `model/provenance.js` owns warnings, tool versions, and recorded request policy.
- `model/reconstruction.js` owns canonical SpacerPlacer metrics and trees.
- `model/values.js` contains null-preserving object/list primitives.
- `useAdapterMembership.js` owns the optional legacy-manifest request lifecycle.
- `useArtifactDownloads.js` owns authenticated download state and transport calls.

Do not add backend alias lookup to a `.vue` file. Extend the owning normalizer,
add a focused model test, and keep components limited to canonical fields.
Result components are presentational: add request side effects to a feature
composable and export it through `index.js`.
