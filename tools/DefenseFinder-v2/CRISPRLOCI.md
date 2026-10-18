# DefenseFinder-v2

This is a modified GPL-3.0 derivative of mdmparis/defense-finder at
513ec1d4843592f8e4f3ca04a24e3675d89fba43. Modifications dated2026-10-12 add the
bounded offline `defense-context` command. Original authorship and LICENSE are
retained. The derivative does not claim to fix every possible upstream bug.

## License decision

The Python source is GPL-3.0. Model release1.2.3 at
5bdb41257fc7234fb8c3b70e3d362c52541860dd has an unambiguous GPL-3.0 LICENSE and
metadata. The exact model LICENSE is retained in third_party; the source asset
hash manifest is defense_context/models.lock.json. Releases1.2.4 through3.1.0
contain conflicting GPL-3.0 and CC BY-NC-SA4.0 license declarations in metadata.
They are not cleared or used here. The current upstream code can support those
newer models, but this adapter deliberately selects only the six1.2.3 families.

Primary sources:
- https://github.com/mdmparis/defense-finder/blob/513ec1d4843592f8e4f3ca04a24e3675d89fba43/LICENSE
- https://github.com/mdmparis/defense-finder-models/blob/1.2.3/metadata.yml
- https://github.com/mdmparis/defense-finder-models/blob/1.2.3/LICENSE
- https://github.com/mdmparis/defense-finder-models/blob/3.1.0/metadata.yml

MacSyFinder2.1.4 and Pyrodigal3.7.1 (including Prodigal) use GPL-3.0. HMMER3.4
uses BSD-3-Clause. Their notices are retained in third_party. Other minimal
runtime dependencies use permissive licenses except certifi's MPL-2.0 data/code.
Keep license texts, corresponding source and modification notices when conveying
a combined runtime. A private GitHub mirror does not remove recipient source
rights if binaries are distributed. Ordinary server execution does not itself
invoke GPL network-source terms. No optional pretrained embedding assets are
fetched or incorporated by this workflow.

## Offline provisioning and execution

Create an isolated Python3.12 environment and install requirements-context.lock.
Install this repository without dependency resolution (`pip install --no-deps .`),
or run `python -m defense_context` from the source tree. Provision HMMER3.4
separately and retain its source/notices. The classic upstream CLI additionally
requires requirements.txt; the new command does not import it.

An operator must fetch the model source at the exact1.2.3 commit into
`MODELS/defense-finder-models`. The CLI verifies all832 model/definition/license
files against its source-controlled manifest and refuses changed or extra model
files. Do not fetch or update models during a web job.

```
defense-context --input genome.fasta --output new-result-directory \
  --models-dir MODELS --hmmer /operator/runtime/bin/hmmsearch --workers 1
python scripts/native_context_smoke.py --models-dir MODELS \
  --hmmer /operator/runtime/bin/hmmsearch
```

## Scope and contract

Input: up to20 source DNA sequences and2Mb total. Pyrodigal3.7.1 meta mode calls
complete ordered proteins separately per contig. Identifiers are generated and
replicon boundaries are preserved with MacSyFinder gembase input; topology is
explicitly linear. Calls never connect the ends of an uploaded fragment.

Output: report.json schema defensefinder.context1.0.0 with source hashes,
one-based inclusive CDS coordinates, native partial flags, genetic codes,
translation evidence, model-rule systems, component assignments and independent
profile hits. Null/missing biological activity is not replaced by a functional
claim. System completeness is the native model-wholeness score. Only submitted
sequences and historical model1.2.3 are assessed. CasFinder and newer model
families are not assessed. Source-bound CasAndra evidence can be compared without
silently merging its classifications with these results.

The result is a repertoire annotation, not a prediction of activity, synergy,
virulence, or absence of all other defenses. Full genomes and richer current
model coverage require separate capacity and licensing review.

## Fixes and regression evidence

- Refactored upstream system summarization to avoid a many-to-many merge when
  one protein participates in multiple system calls; each system now appears once.
- Gene counts represent distinct component proteins, not duplicate profile rows.
- The new adapter consumes the actual `.res_hmm_extract` suffix, preserves isolated
  hits, rejects unbound source IDs, pins all model bytes and avoids network checks.
- Per-contig gene order, source translations, partial calls and fail-closed identity
  are covered by focused tests and the independent webserver boundary validators.
- A real40kb upstream nucleotide fragment produced42 proteins,1RM_Type_II system,
  2components and5profile hits with the full six-family native CPU scan.

See upstream README for method citations. Primary method:
https://doi.org/10.1038/s41467-022-30269-9

The CDS SHA-256 is calculated from the exact oriented submitted DNA, retaining
IUPAC ambiguity. Pyrodigal normalizes non-ACGT bases to N when calling genes;
its protein translation therefore emits X for a codon containing ambiguity.
The wrapper preserves that native protein and hashes its exact sequence. Native
regressions cover ambiguous CDSs on both input orientations. The server checks
complete codon frames, strand-aware partial boundaries and complete stop codons.

Refinement dated 2026-10-13 keeps scientific versions, model bytes and schemas
unchanged. FASTA sequence chunks are joined once, and input provenance hashes the
same bounded byte snapshot used for gene calling. Native tables stream under the
same size/symlink guard; HMM scores must be finite, protein intervals and global
FASTA ranks must bind their called proteins, and repeated system metadata must
agree. Overlapping systems, every component assignment and original profile tie
order remain intact. Native log handlers close after failed or interrupted scans.
The context boundary regressions also cover blank headers, source limits, shared
systems, conflicting metadata and the global-rank/local-position distinction.

## Configuration isolation (2.0.1)

The context entrypoint rejects undeclared configuration consumed by MacSyFinder
2.1.4 before execution: the selected VIRTUAL_ENV/MACSY_CONF/global path, home and
working-directory files, and selected family model configuration files. This
prevents ambient weights from changing system scores outside the declared policy.
The preflight does not mutate the caller's environment or scientific model assets.
