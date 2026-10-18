# Repeat structure analysis

The scientific implementation is the independent [CRISPRrepeat repository](https://github.com/Alexander-Mitrofanov/CRISPRrepeat), checked out under `tools/CRISPRrepeat` and pinned in `tools/tools.lock.json`. Integration 0.5.0 invokes its CLI through an operator-owned runner. The API never imports ViennaRNA or executes folding.

Two queued web modes are available in this source release:

- `repeats`: supplied repeat FASTA, at most 200 nt/record and 1000 records within the service limits. Declare `molecule=DNA` (default, both forward and reverse-complement hypotheses) or `molecule=RNA` (as-transcribed forward sequence). Mixed T/U is rejected. Genomic coordinates and context scores remain missing.
- `repeat_context`: genomic DNA goes through CRISPRidentify, then each accepted observed repeat is folded both in isolation and within the reconstructed array. This mode does not require related arrays or evolutionary reconstruction. The category policy selects accepted candidates.

The web model is pinned to ViennaRNA 2.7.0, Turner 2004, 37 C, dangles 2, salt 1.021 M, local window 150 nt and maximum span 80 nt including endpoints, with lonely pairs allowed. The standalone tool supports explicit parameter sensitivity analyses. Context boundaries are the detected array boundaries; each local window may include spacers and other repeats. Short contexts use an explicitly reported shorter window/span. Limits of 5000 nt/array, 20000 total context nt, 1000 observed repeats, and 20 MB JSON cause explicit failures rather than silent trimming.

The isolated MFE pair set is an inspectable hypothesis. Reported values are mean support for that fixed set under the isolated/global partition function and the local-window calculation, plus their difference. An empty set gives null scores with `no_reference_pairs`. Local pair values are window averages and must never be summed as a partner distribution; accessibility is ViennaRNA's native single-base calculation. A context dot-bracket structure is not invented.

The UI shows up to 20 repeat/orientation records within a 700 KB preview budget, marks truncation, and shows at most 20 individual competing pairs with the original count. The complete independently validated JSON, TSV, RNA FASTA and isolated dot-bracket text are downloadable. Both Integration and the backend independently check source binding; backend checks include hashes, coordinate/strand transformations, pairing bounds, MFE pair identity and comparison arithmetic. Raw tool sidecars are excluded from publication.

No calibrated strand, subtype or family classifier is claimed. The standalone tool optionally searches an operator-supplied frozen provenance-bearing reference; the initial public web workflow reports `reference_not_configured`. Sequence similarity is not a subtype assignment. Nothing is automatically submitted to an external analysis service.

## Scientific verification

Standalone tests run real ViennaRNA calculations and cover single-window/global equivalence, context competition, source tampering, reverse mapping and missing values. An attributed offline benchmark includes the published Reimann et al. R6/R7 contrast and a cleaved repeat with very low support. The window sensitivity and counterexample are documented in [the scientific validation report](../tools/CRISPRrepeat/docs/REIMANN_REPRODUCTION.md). These selected cases are not a population accuracy benchmark or a cleavage classifier validation.

## Operator preparation

```sh
python3 scripts/sync_tools.py --check
python3 scripts/provision_tools.py --lock tools/tools.lock.json
python3 scripts/manage_release.py build --lock tools/tools.lock.json
python3 scripts/manage_release.py verify /absolute/path/to/staged/release
backend/.venv/bin/python scripts/browser_repeats.py /absolute/path/to/staged/release
```

The toolset provisions a separate repeat venv and checks the pinned Python binding/distribution version and native extension hash. Scientific manifest schema 4 adds the repeat identity while preserving prior release schemas. The browser check starts temporary local API/worker/frontend processes and tests DNA, RNA, real detected-array context, downloads, source provenance, reload, small-screen layout and compiled CSP. It does not activate the release or operate the deployed services.

## Verified candidate

CRISPRrepeat commit `4a98cc465656c4db3ed0faccaafe476707d20ddf` is pushed and GitHub CI passed. The application release `940c298927d5d1019cea39c0b604c2d6a02da9ccdc6a92e27952de7a65d98c22` is staged, verified and gated. Production was not activated.

Validation passed: 44 standalone tests; 468 backend/Integration/runtime tests; 145 frontend tests; source checks; production UI build/public scan; 13 real HTTP cases; three real browser workflows with mobile layout, CSP, source provenance and authenticated complete downloads. See [machine-readable evidence](repeat-verification.json) for exact identities and local evidence paths.
