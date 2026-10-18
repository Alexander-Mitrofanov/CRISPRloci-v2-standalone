# Reverse viral search council — 2026-09-28

The research, scientific-methods and integration reviewers discussed the goal
before implementation: accept user-supplied spacer occurrences and return
candidate viral sequence matches from a reproducible local BLAST database.

## Decisions

- Use all genomic nucleotide FASTA members of one NCBI viral RefSeq release,
  including phages. This is a defined reference snapshot, not all viral diversity
  or the complete GenBank collection. Keep accession versions and original titles;
  unprovided taxonomy/host information remains unknown.
- Freeze an operator download plan using the NCBI release number, full matching
  member list and official MD5 catalog. Add our own SHA-256 digests during download.
  Abort on changed/missing source bytes; do not silently follow a newer release.
- Separate metadata planning, network download, local BLAST indexing and offline
  search. Build/validate an immutable reference before publishing it. No API job
  downloads or updates a reference. Full viral data are provisioned on the VM only.
- Retain a normalized disk-backed source inventory so each reported alignment
  can be checked against actual reference bases, independently of BLAST's output.
  Bind all native index volumes and catalog bytes to a reference manifest digest.
- Use `blastn-short` with explicit settings. Accept one full-query ungapped HSP
  with at most two substitutions by default, or exact reported HSPs with zero.
  This is an exploratory computational tolerance, not validated biological
  targeting confidence. Do not combine partial HSPs into apparent full coverage.
- BLAST is heuristic: successful execution is not proof of exhaustive retrieval.
  Local clipping can omit terminal-mismatch matches. Record retention limits,
  effective parameters, database ordering and engine identity; do not report a
  unique/best virus or infer infection, host range, PAM compatibility or immunity.
- Bound queries, raw output, accepted hits, runtime and native process groups.
  Disable inherited NCBI configuration, remote data loaders and BLAST telemetry.
  Native/resource failures never become successful no-hit results.
- Preserve current standalone exact searches and the existing v2 web runtime pin.
  This release provisions a callable reverse search and VM reference; it does not
  activate a new web service or update a running worker.

Official sources and reproducible commands are listed in
[viral reference methods](viral-reference-methods.md). Tests use generated DNA
fixtures and installed BLAST binaries, never a downloaded viral corpus.
