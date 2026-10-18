# Standalone package provenance

This repository contains runtime `0.3.0.dev2` and the selected cap15 bundle
`casandra-full361-v3-cap15-2026-09-20-v1`, extracted from CasAndraProject on
September 26, 2026. It supersedes the older `tool/` snapshot for this standalone
distribution. It does not update the original workspace or a running service.

The model manifest SHA-256 is
`966c4522d38585db6c3dbca14566aa25d520b7c165f487abe5f35845441eb1b3`.
The runtime release-manifest SHA-256 is
`a7efdacc98a36be684c8ca9fc70ae397b9d3a1f9782e2f5f01a9c60b39bf6e8e`.

[source.json](source.json) records the original workspace-relative paths and
every runtime module's original hash. These paths are historical references;
they are not required for inference or tests. Runtime-module hashes were
checked against the frozen release manifest before copying. All 19 model
artifacts were verified against the cap15 manifest.

Packaging changes:

- Replace the runtime's older default model with the byte-identical selected
  cap15 bundle. Preserve every Python runtime module and model artifact.
- Add standalone installation instructions, repository URLs, citation version,
  source-distribution contents, CI, and installed-wheel verification.
- Correct the minimum Python version to 3.11, matching the required
  scikit-learn 1.8.0 package metadata. Runtime Python modules are unchanged.
- Port study-relative tests to the packaged model and local fixtures. Preserve
  the original specificity evidence and V3/margin10 policy comparisons.
  Cap15 may retain calls withheld by V3, so the cap15 test checks cumulative
  rejection and unchanged evidence instead of requiring the old V3 labels.
  A separate test verifies all selected cap15 floor-map entries.
- Update model-specific assertions to the cap15 manifest and the September
  full361 cassette classifier. The canonical cassette still classifies as II-A;
  its classifier confidence is `0.5682433452220029` rather than the older
  model's `0.5449020585645451`.

[fixtures.json](fixtures.json) records the source and checksum of each small
fixture extracted from earlier studies. Tests exercise engineering contracts;
they do not establish biological accuracy or independent sensitivity.

Cap15 applies the selected score-floor cap to the inherited and later
specificity policies while retaining coverage and competitor checks. It can
withhold a baseline assignment; callers must preserve the accompanying raw
evidence and abstention semantics. No model was fitted or tuned during this
repository extraction. Existing studies remain in the research workspace.
