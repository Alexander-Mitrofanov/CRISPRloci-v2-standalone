# Precomputed web examples

Every public analysis mode has its own input and real precomputed result under
`frontend/public/examples/`. RNA repeats and Type V tracrRNA also have separate
examples, selected by the form's molecule/model option. No service connection
is needed to view the stored results.

Load example fills the current method's input and recorded options. View
precomputed result appears only while the exact input bytes and every submitted
scientific option match. An edit restores the usual live submission path.
The result uses the same Vue components, headings, tables, plots and uncertainty
semantics as a completed worker job. An example notice and input/output download
links sit outside the scientific result view. Authenticated job downloads and
continuation actions remain specific to live jobs.

Generate the examples with an existing verified local release:

```sh
backend/.venv/bin/python scripts/generate_examples.py --release .runtime/current
```

The command prepares small reference panels from pinned tool examples, starts
an isolated loopback API and worker, submits each exact input twice as separate
jobs, and compares their scientific summaries before saving public assets. It
changes neither production services nor the selected runtime. Raw job identities,
logs and comparison evidence stay under ignored `.runtime/example-generation/`.
Only explicitly constructed public input/result contracts are exported.

Stored input metadata includes raw and normalized SHA-256 hashes, record/base
counts, request options and the immutable release identity. Exact tool versions
and commits remain in the normal methods/provenance view. Execution durations
are preserved from the original run but excluded from independent-run equality;
they naturally vary between jobs.

The sequence/reference identifiers are masked before execution. The reference
search examples use three masked spacer records, three synthetic viral records,
three masked phage genomes, and the CRISPRmap engineering reference. To reproduce
them, use these exact panels and tool revisions: full production databases may
yield additional or different matches. The miniature phage panel cannot establish
calibration; its candidates retain null native FDRs. The CRISPRmap family labels
are explicitly synthetic. These are demonstrations of real execution and
presentation, not new biological validation claims.

`scripts/browser_examples.py` compares every scientific result tab with the
second independent job in headless Chromium at desktop and mobile widths. Its
screenshots and report remain local. `scripts/check.sh frontend` enforces the
coverage, input binding, option matching, normal rendering, public-data scan and
frontend build checks.

Verification on 11 October 2026: all 17 examples passed independent-run
scientific equality and desktop/mobile pixel comparisons of every scientific
tab, including expanded evidence tables. The frontend gate passed 350 tests,
formatting, lint, source boundaries, production build and public-data scans.
