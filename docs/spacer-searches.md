# Spacer searches in the web server

The **Spacer searches** workspace provides two queued workflows. Both accept
FASTA upload or pasted text, preserve separate drafts when changing direction,
and report coordinates on the forward target/reference sequence as 1-based
inclusive intervals.

| Direction | Input | Reference | Match rule |
| --- | --- | --- | --- |
| Find spacers in my sequence (`protospacer`) | Linear DNA contigs or small genomes | NCBI-hosted 2017 research spacer collection | Full-length exact spacer matches, both strands |
| Find viruses for my spacers (`viral_search`) | DNA spacers, one per FASTA record | NCBI viral RefSeq | BLAST short-query retrieval, full-query ungapped alignments, 0–2 substitutions |

Viral spacer eligibility is 18–80 nt and unambiguous A/C/G/T. Short or ambiguous
records are explicitly skipped; records over 80 nt are rejected at submission.
The ordinary configured record, request and job limits also apply. Results show
candidate accession/title, coordinates, strand, substitution count and distinct
spacer-sequence support. Repeated copies of one spacer do not add independent
sequence support. Viral records can include segments and are not species counts.

Viral retrieval is heuristic: clipped terminal mismatches may fail the
full-length filter, and a no-hit result does not demonstrate biological absence.
Neither workflow evaluates PAM compatibility, array overlap, functional
targeting, infection or host range. Missing and unfinished evidence remain
separate from completed searches with no accepted matches.

Complete validated results are available as JSON and TSV, alongside the
reproducibility archive. Public exports omit private native BLAST evidence,
executable paths and operator configuration. Result recovery and downloads retain
the existing per-job bearer capability protections.

## Operator setup

Use CRISPRspacer 0.3.0, Integration 0.7.0 and backend 2.3.0. The private standalone
repository remains [CRISPRspacer](https://github.com/Alexander-Mitrofanov/CRISPRspacer).
Its `docs/viral-reference-methods.md` documents NCBI download, integrity checking,
normalization, native index creation, deep validation and reproducible searches.
The existing `docs/viral-refseq-237-snapshot.json` records the verified full snapshot.
No full reference download is needed on a development workstation: local and
release tests build the bundled three-record fixtures.

Configure both path/digest pairs before starting the API and worker:

```text
CRISPR_WEB_SPACER_DATABASE=/absolute/spacer-bundle
CRISPR_WEB_SPACER_DATABASE_SHA256=<sha256 of spacers.sqlite3>
CRISPR_WEB_VIRAL_REFERENCE=/absolute/viral-reference
CRISPR_WEB_VIRAL_REFERENCE_SHA256=<sha256 of reference.json>
```

Keep these directories immutable and readable by the service account. Anonymous
requests cannot select filesystem paths, native commands or arbitrary reference
identities. The API advertises each direction only when its configuration is
present. Queue schema 6 stores the selected reference digest and substitution
policy per job. Only workers invoke the tools. A reference change invalidates a
job's execution instead of silently substituting different data.

The pinned runtime's `runners.json` supplies CRISPRspacer and the inventoried
`blastn` executable. The Integration adapter validates native output, then the
backend independently verifies source bindings, reference digest, raw HSPs,
subject sequence/coordinates, accepted matches, query states and aggregate counts.
Budget violations and invalid evidence fail the job instead of producing a
successful truncated result.

Run `scripts/check.sh backend` with native BLAST available and
`scripts/check.sh frontend`. Browser acceptance is documented in the frontend
README. `scripts/gate_release.py` builds tiny reference fixtures and exercises
both modes through the real queued API before allowing activation. Production
rollout additionally validates the full references under service isolation,
drains work, and backs up the full job/queue directory before schema migration.
