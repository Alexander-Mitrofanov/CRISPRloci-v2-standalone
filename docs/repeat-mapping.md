# Repeat reference mapping

The `repeat_map` workflow exposes the pinned CRISPRmap v2 `classify` command as
an isolated queued analysis. Its public label is **Map repeats to reference
families**. It maps supplied DNA repeats to an operator-selected immutable
annotated JSON reference. This is an independent historical-reference mapper in the [CRISPRmap-v2 fork](https://github.com/Alexander-Mitrofanov/CRISPRmap-v2), pinned by [the source lock](../tools/tools.lock.json). It does not recreate the original paper's classification pipeline or perform modern catalogue/MAP2 group search.

## Operator setup

Build the new toolset and application release using the normal operator workflow
in the root README. The lock's `runtime_profile` distinguishes this runner recipe
from earlier toolsets with the same scientific commits. Existing environments
and the active release are never modified. Provision the reference explicitly;
a web submission cannot supply a path, download a database or update a tool.

Set both variables in the local service environment before startup:

```text
CRISPR_WEB_REPEAT_REFERENCE=/absolute/operator/reference.json
CRISPR_WEB_REPEAT_REFERENCE_SHA256=<sha256-of-the-reference-file-bytes>
```

The reference must be a regular absolute file without symlinks and must follow
`crisprmap.reference/1.0.0`. The scientific release must contain the `crisprmap`
identity and runner. The API advertises this mode only when configured. The
worker verifies its file digest before execution, and the independent backend
validator checks it again before publishing results.

Use the standalone tool's explicit reference-import/build workflow. Its bundled
`examples/reference.json` is synthetic and is suitable for engineering checks
only. Historical reference terms and provenance are described in the tool's
README. The portal does not distribute reference datasets in the static frontend.

## Input, limits and interpretation

POST ordinary DNA FASTA using `mode: "repeat_map"`. At most 100 records of 200
bases each are accepted, further bounded by the server's general limits. RNA
input is not accepted in this mode. Ambiguous DNA remains visible as unsupported.

The effective policy permits 0–3 unit-cost global edits and 1–100 reported best ties per query (defaults: 3 edits and 20 ties), both reference orientations, and a 100,000-comparison native work budget. Results retain the requested distance/top-k, total best ties, reported tie count and native truncation. The browser displays up to 100 retained ties per query and separately marks preview truncation. Large diverse batches can exceed this work budget; submit smaller
batches. The native stage deadline is at most 300 seconds within the overall
worker deadline. A budget failure never masquerades as a complete no-hit result.

- **Exact:** a sequence or reverse complement exists in the reference. Assign a
  historical label only if all closest occurrences have the same nonmissing label.
- **Near:** closest references within the effective edit bound; family/motif values remain
  unassigned. Labels describe candidate/reference evidence only.
- **No hit:** no eligible reference within the retrieval bound; novelty is unknown.
- **Unsupported:** ambiguous query bases prevent unambiguous matching.

Conflicting and incomplete annotations are separate outcomes. Every tied
reference occurrence contributes to labels before the 20-hit display cap.
Reference-relative direction is not transcription strand. No calibrated
classification confidence, Cas subtype, host or structural motif is inferred.

Downloads are `repeat-mapping.json` and `repeat-mapping.tsv`, also in the protected
result ZIP. JSON includes all submitted query results, exact sequences, counts,
labels, displayed alignments, reference identity and truncation flags. TSV gives
one query-level row and neutralizes spreadsheet formulas. Native HTML, command
logs with private paths and intermediate FASTA are not published.

## Local verification

```bash
bash scripts/check.sh
backend/.venv/bin/python -m pytest -q --import-mode=importlib -o addopts='' backend/tests/test_repeat_map.py
```

The release gate adds an exact, no-hit and ambiguous-query control against the
pinned tool's synthetic reference. The browser script is
`frontend/scripts/browser_repeat_map.py`; serve the local frontend build with
`VITE_API_BASE_URL=https://analysis.example.org`, then run the script's `--help`
for its origin/output flags. It intercepts the API with a fixture produced by the
real mapper; the separate backend test exercises the real HTTP/worker/download
chain. Neither test contacts production.
