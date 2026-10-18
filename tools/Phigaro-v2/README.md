# Phigaro-v2

A bounded offline research-worker distribution of Phigaro 2.4.0, retaining its
MIT-licensed source and basic scoring method. Upstream:
https://github.com/bobeobibo/phigaro commit
`1ff5f85cee31e418bce24e4cd51c7528c43bc968`. This distribution is version 2.4.1.

It predicts candidate prophage regions in submitted sequences. In CRISPRloci it
joins these intervals to source-bound exact spacer matches outside reported
arrays. Sequence matches do not establish active targeting, inducibility,
infection, immunity, host range or a complete genome's absence of prophages.

## Reproducible worker runtime

Use Python 3.12, `requirements.lock`, Prodigal 2.6.3 and HMMER 3.4. Native programs
remain separate executables. Install this package, then run:

```
phigaro-v2 predict --input genomes.fasta --output new-result \
  --models-hmm /operator/models/pvogs.hmm --models-sha256 SHA256 \
  --prodigal /operator/bin/prodigal --hmmsearch /operator/bin/hmmsearch
```

Only operators acquire references. Jobs never download or update models. The
official pVOG archive is
https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
with SHA-256 `1c5bd9ec035fcae466760ab58c3e2dbd088e12cf1de608e8a1074b8b98bda5a9`.
Concatenate regular `.hmm` members in lexical full-member-name order (reject
unsafe names, links and duplicates). The 9,518-profile result has 767,855,419
bytes and SHA-256
`0ddeafff975f1cfd86686e95f52c7a9980f9ad10fa082797313097ca42be97ed`.
Model data are not redistributed in this Git repository.

Limits: 20 records, 2 Mb total, 10,000 predicted genes, 8 CPUs, 50 MB native HMM
rows, per-command timeout. Records shorter than 20 kb are explicitly ineligible
for the single-genome gene-calling model; they are not labelled negative.
Operators should bound the entire worker process group, output and deadline.
The JSON records source checksums, intervals, translated CDS evidence, partial
gene flags, retained native hits, profiles, scores and exact scoring parameters.
All coordinates refer to the 1-based inclusive forward input sequence.

## Changes and validation

Fixed a NumPy/list sentinel error that dropped terminal prophage intervals and
mutated list inputs, and the final-gene slicing error in the upstream batch
adapter. Added shell-free native execution, checked return codes, bounded
FASTA/model/output validation, deterministic HMM tie handling, source/protein
reconstruction and explicit missing/eligibility states. Only the bounded CLI is
installed; legacy web/batch entry points remain source history.

The scoring `gc` is Prodigal's reported gene GC covariate, preserving the
threshold calibration. Prodigal 2.6.3 reverse-strand GC can differ slightly from
the exact CDS fraction; `source_gc` reports the latter independently. These
values are intentionally distinguished. A freshly compiled Prodigal with modern
compiler optimization was observed to corrupt annotation metadata through
upstream overlapping `sprintf`; the pinned Conda build passed native metadata
checks. Never substitute an unvalidated native executable.

Run `python -m pytest tests`. The native lambda J02459.1 pilot
using all 9,518 profiles produced 61 CDSs and one candidate region. It checks
execution and evidence integrity, not sensitivity/specificity or applicability
to every organism. No claim is made that every upstream defect has been found.

## Code navigation

The installed `phigaro-v2` command enters `phigaro/bounded.py`: bounded FASTA
reading, native execution, source-bound CDS/HMM validation and JSON/TSV reports.
`phigaro/finder/v2.py` owns the original triangular window and two-threshold
interval rule; its vectorized row sums preserve the original addition order.
`phigaro/const.py` retains published defaults and pVOG penalty lists.
`tests/test_bounded.py` checks scientific evidence; `tests/test_refinement.py`
checks exact numeric parity, record isolation, deterministic ties and malformed
native outputs. The legacy `batch`, `cli`, `data` and `to_html` modules document
the original workflow and are not worker entry points.

Group genes once per source, and keep thresholds in `PARAMETERS`. Preserve
Prodigal's native `gc`, source-forward coordinates, partial flags and strict
ambiguous-codon translation. Short records remain unassessed. Reject protein
text before its native header and enforce gene-count/duplicate bounds while
reading, before accumulating an invalid report.

## Licenses and attribution

Phigaro: MIT, original [LICENSE.md](LICENSE.md) retained. Cite the original
Phigaro publication (Starikova et al., Bioinformatics 2020, doi:10.1093/bioinformatics/btaa250)
and pVOGs (Grazziotin et al., Nucleic Acids Research 2017,
doi:10.1093/nar/gkw975). The Phigaro MIT license does not relicense dependencies
or reference data.

Prodigal: GPL-3.0-or-later (University of Tennessee/UT-Battelle); HMMER: BSD-3-Clause
(Howard Hughes Medical Institute). They execute separately. NumPy and six are
BSD-3-Clause and MIT, respectively; Biopython uses the Biopython License and/or
BSD-3-Clause. Installed distributions retain their own notices. Full copied
native notices are under `licenses/`; provide corresponding source if binaries
are distributed. Hosting the free research service does not distribute these
worker binaries to browser users.

pVOGs is publicly available NCBI molecular research data. NCBI's reuse policy
https://www.ncbi.nlm.nih.gov/home/about/policies/ places no additional NCBI
restriction on molecular data and retains a third-party-rights caveat. Treat
that as a data-use policy, not a blanket MIT grant; preserve provenance and
citations. ViennaRNA and DefenseFinder have separate terms in their repositories.
