# RNAalifold-v2

Bounded, offline comparative folding of a user-selected homologous RNA family.
The wrapper preserves the original sequence identities, aligns distinct sequences
with MAFFT core, and calls the unmodified ViennaRNA 2.7.0 comparative-folding API
used by RNAalifold. It does not infer homology, reverse-complement the inputs,
predict transcript boundaries, or establish RNA activity.

## Scientific purpose

CRISPRloci predicts tracrRNA candidates and describes individual repeats. A
researcher can then curate homologous, consistently oriented RNA candidates and
ask whether they retain a common structural core. MAFFT and RNAalifold were used
in published comparative tracrRNA analyses, including
[Gasiunas et al. (2020)](https://www.nature.com/articles/s41467-020-19344-1)
and [Liu et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906591/).
Subtype identity alone is not evidence that arbitrary candidates belong in one
family. A consensus single-family fold is not a repeat–anti-repeat duplex model.

## Installation and execution

Use an isolated Python 3.11+ environment, install `requirements.lock`, then
`pip install --no-deps .`. Provision MAFFT core separately and pin its exact
binary/version in the operator runtime inventory. The original local test used
MAFFT 7.505; the managed CRISPRloci runtime pins MAFFT 7.526. Inherited `MAFFT_BINARIES` is removed from both version and alignment child environments without changing the caller environment or temporary-directory settings. MAFFT RNA extensions are neither needed nor invoked.

```sh
rnaalifold-v2 --input family.fasta --output new-result-directory \
  --molecule RNA --mafft /operator/runtime/bin/mafft --threads 1
```

Inputs contain 2–50 FASTA records, each 30–500 unambiguous nucleotides, oriented
5′ to 3′. RNA accepts ACGU; DNA accepts ACGT and is converted to RNA for folding.
Identifiers are unique and restricted to safe ASCII tokens. At least two distinct
RNA sequences are required. Exact duplicates receive one equal representative
weight; original identifiers and duplicate membership are retained. Nonidentical
related sequences are not phylogenetically reweighted.

The alignment uses fixed FFT-NS-2 (`--nuc --inputorder --retree 2 --maxiterate 0`),
with a 1,000-column maximum. The comparative model uses RNA Turner 2004 parameters,
37 °C, covariance factor 1 and noncanonical factor 1. MFE scores include covariance
contributions; they should not be interpreted as the thermodynamic energy of a
single supplied RNA. The partition function supplies base-pair probabilities.
All predicted MFE pairs and other pairs with probability >=0.01 are retained.

`report.json` includes the source hashes, original-to-representative crosswalk,
RNA alignment, per-column counts, consensus and dot-bracket structure, comparative
MFE score, ensemble free energy, pair probabilities and descriptive pair counts.
A compensatory sequence-pair count means two canonical observed pairs differ at
both positions; it is not a phylogenetically corrected event count or a test of
statistical significance. Empty predicted pairing is a valid completed result.

The private `work/` directory contains native MAFFT alignment and folding evidence
for the independent server validation boundary. Reports are atomically installed
only after successful execution. Existing output directories are refused. There
are no job-time downloads, shell commands assembled from user data, or models
changed by submitted input. Process containment and the outer job deadline are
owned by the isolated worker.

## Refactoring and verification scope

This repository separates input/alignment/evidence transformations (`core.py`)
from orchestration (`cli.py`). It adds strict source reconstruction, deterministic
duplicate handling, explicit missing/zero-pair semantics, native subprocess error
propagation, and resource bounds. It does not fork or change ViennaRNA's numerical
engine and makes no claim to have eliminated all upstream bugs.

`core.py` owns FASTA normalization, representative/source crosswalks, alignment
reconstruction, column/pair evidence and the native fold. `cli.py:align` owns the
fixed MAFFT invocation; `cli.py:run` installs the validated native evidence and
report. Policy and bounds live in `core.py`. There is no persistent fold cache.
Native model defaults are explicit, and an evaluation-only comparative compound
captures family-scaled energy/Boltzmann tables before the folding compound uses
copies. Initialization runs inside the isolated worker process.

Run `python -m pytest -q` in the isolated environment with MAFFT on `PATH`.
The tests include real MAFFT/ViennaRNA execution, a compensatory-pair toy family,
exact duplicate handling, unpaired families, malformed records, mismatched
structure notation and failed native commands. Fresh-process native oracles
verify exact default energies and pair probabilities under changed process-global
settings; malformed native output and source/alignment boundaries fail closed. The synthetic fixture is not a
biological tracrRNA. The server independently checks native alignment reconstruction,
source hashes, dot-bracket pairing, probability bounds and mass, and recomputes
column and covariation counts without importing this package.

## Licenses and credit

Original wrapper code is MIT licensed. **That does not relicense ViennaRNA.**
ViennaRNA 2.7.0 has its own [upstream license](https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.0/COPYING),
preserved in `licenses/ViennaRNA-2.7.0-COPYING.txt`. It permits research,
educational and commercial use and modification with attribution and a restriction
on fee-based redistribution; inclusion in a commercial product asks users to
contact the authors. The intended free research web-server use preserves credit
to the authors and the Institute for Theoretical Chemistry, University of Vienna.
A commercial product or fee-based distribution needs separate review against
those conditions. This repository does not redistribute the ViennaRNA engine.

MAFFT core has a BSD-style license, preserved in `licenses/MAFFT-core-license.txt`.
Credit Kazutaka Katoh and the MAFFT developers. No MAFFT extensions or associated
third-party engines are redistributed here.

Cite [Bernhart et al., RNAalifold (2008)](https://doi.org/10.1186/1471-2105-9-474),
[Lorenz et al., ViennaRNA Package 2.0 (2011)](https://doi.org/10.1186/1748-7188-6-26),
and [Katoh and Standley, MAFFT v7 (2013)](https://doi.org/10.1093/molbev/mst010),
as appropriate. For exact license text, the preserved notices take precedence
over this operational summary.
