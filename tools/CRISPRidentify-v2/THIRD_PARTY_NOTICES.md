# Third-party notices

## CasAndra

CRISPRidentify v2 vendors CasAndra, an offline CPU-only Cas protein and
CRISPR-Cas cassette annotation tool. The supplied source snapshot identified
itself as version `0.1.0.dev0`; this integration identifies itself as
`0.1.0.dev0+crispridentify.3`.

The source tree was supplied in the development workspace without VCS
metadata, an immutable revision, or a published artifact locator. This
distribution therefore does not assert an upstream repository, commit, or
release. That provenance limitation is recorded explicitly in
`vendor/casandra/VENDORING.json`.

CasAndra copyright (C) 2026 CasAndra contributors. CasAndra is distributed
under the Apache License, Version 2.0. The complete license is included at
`vendor/casandra/LICENSE`, and the upstream copyright/license notice is
preserved at `vendor/casandra/UPSTREAM_LICENSE_NOTICE`.

The vendored snapshot has downstream changes by the CRISPRidentify
maintainers. They add model-manifest and artifact verification at every loading
boundary, expose deterministic profile-family evidence, normalize independent
protein and final cassette class/type/subtype winners while preserving their
raw evidence, and make the genome gate type-aware for subtype-null Type VII
candidates using `type_vote_fraction` and
`direct_profile_type_aggregation`. The final cassette gate consumes the
normalized claim, and a validated per-run count records masked cassette
conflicts. The changes also assign a downstream version suffix. Packaging
metadata was updated to require a PEP 639-aware setuptools release and to
include the license files. These changes are not an
official upstream CasAndra release. The known origin, its limitations,
downstream changes, and checksums are recorded in
`vendor/casandra/VENDORING.json` and `vendor/casandra/SOURCE_MANIFEST.sha256`.

## CasAndra model and reference-data attribution

The bundled deployment refit is
`casandra-cas-only-cpu-bundle-v4-fast-final`. It contains 341 Cas/hard-negative
HMM profiles and a cassette architecture ExtraTrees model. Its own manifests
mark it as a deployment refit, not an evaluation model; it must not be used to
re-estimate accuracy on records used to construct or fit it.

Type VII has limited protein-level test support in the bundled model card, and
the model card identifies an independent whole-genome benchmark as remaining
work. CRISPRidentify therefore labels subtype-null, type-only `CAS-VII`
cassette calls experimental and preserves their supporting evidence for review.

The classification and profile metadata draw on:

> K. S. Makarova, S. Shmakov, Y. I. Wolf, et al.  
> An updated evolutionary classification of CRISPR-Cas systems including rare
> variants. Nature Microbiology (2025).  
> DOI: 10.1038/s41564-025-02180-8

The associated reference-data release is archived at Zenodo under
DOI 10.5281/zenodo.17388109 and is attributed under CC BY 4.0. CasAndra
transforms the reference material through profile selection, HMM construction,
negative-profile selection, metadata normalization, and a frozen deployment
refit. The software citation is preserved in
`vendor/casandra/CITATION.cff`.

## Pyrodigal

CasAndra gene calling uses Pyrodigal 3.7.1, copyright Martin Larralde,
distributed under GNU GPL version 3 or later. Pyrodigal includes and adapts
Prodigal code and provides the in-process SIMD-compatible gene-calling API.
Source: https://github.com/althonos/pyrodigal.

## PyHMMER

CasAndra profile searches use PyHMMER 0.12.1, copyright Martin Larralde,
distributed under the MIT license. PyHMMER provides in-process bindings to the
HMMER3 search pipeline. Source: https://github.com/althonos/pyhmmer.

## Other Python runtime dependencies

Biopython, Joblib, NumPy, scikit-learn, SciPy, and their transitive
dependencies are installed from the pinned Conda environment rather than
copied into this source tree. Their package metadata and license files remain
authoritative. Review the explicit environment lock and those licenses before
redistributing a complete environment.
