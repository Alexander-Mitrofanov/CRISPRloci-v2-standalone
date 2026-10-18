# CRISPRleader v2 in CRISPRloci

The server integrates the published `CRISPRleader-v2` commit
[`b5310163835649cdd89b0e03afb2a90b5421d675`](https://github.com/Alexander-Mitrofanov/CRISPRleader-v2/tree/b5310163835649cdd89b0e03afb2a90b5421d675) (package `2.0.0a1`), as recorded in [the authoritative source lock](../tools/tools.lock.json). It is a
**context extractor, not an available leader predictor**. The missing inference
policy, biological reference bundle and independent biological evaluation are
not supplied or simulated by this integration.

## Workflows

- `leader` / **Extract leader context**: CRISPRidentify detects arrays, then
  CRISPRleader adapts the source-bound canonical reports, extracts both flanks,
  and validates the result.
- `loci` / **Annotate a CRISPR locus**: the same stage consumes the retained
  Identify output from the complete tracrRNA workflow. Identify is not rerun.

`leader_flank_length` is an integer from 1 through 5,000, default 500. It is
persisted with the queued job and recorded in downloads and provenance. Both
Bona-fide and Possible arrays enter the leader context stage. No side is selected,
no confidence is invented, and a successful empty result does not establish
leader absence. The server's FASTA inputs do not provide circular-topology
assertions: windows clip at contig boundaries with topology explicitly unknown.

Left sequences are input-forward; right sequences are reverse complemented so
both point toward their adjacent terminal repeat. This is a sequence convention,
not an inference of transcription direction. Empty boundary windows remain in
the JSON/TSV evidence. FASTA and interval formats omit empty features.

## Contracts and publication

The API queues jobs; only the contained worker calls the immutable launcher.
Integration uses shell-free `adapt identify`, `extract`, and `validate` commands
in a dedicated CRISPRleader environment. It never invokes the unavailable
`predict` command. The backend independently checks the native source hashes,
accepted arrays, terminal-repeat anchors, reconstructed sequences, sides,
coordinates, counts, effective window length and absence of predictions.

The `crisprleader` summary is a bounded preview. Complete, capability-protected
artifacts are `leader-contexts.json`, `.tsv`, `.fasta`, `.gff3` and `.bed`.
JSON/BED are zero-based half-open; GFF3 is one-based inclusive. These are regenerated
from validated evidence, with public source IDs; private native directories and
process commands are excluded. TSV formula-like cells are escaped.

Scientific manifest version 3 adds the sixth Git source identity while retaining
read compatibility with older manifest versions. The runtime builder exports the
locked Git commit, never uncommitted tool changes. Existing operator update,
verification, staging and gate commands also manage CRISPRleader. Neither public
requests nor active jobs can update a tool or environment.

## Queue migration and verification

The additive queue migration retains jobs, capabilities, artifact records and
rate salts. Tests cover schema 2 and 3 upgrades. Back up the stopped SQLite queue
before activating a new schema; older application releases cannot open a newer
queue schema. A code rollback across that boundary requires its matching database
backup or a separate data root. Never delete existing job data to make a rollback
start.

The positive contract fixture was captured from the pinned tools; tests then
mutate sequences, coordinates, hashes, anchors, settings and prediction claims.
Checks also cover clipping, full exports beyond preview limits, API bounds,
queue persistence and result tabs. The real release gate includes positive and
empty leader inputs plus the complete locus combination. These are software
integration checks, not a biological validation of leader prediction.

## Recorded verification — 27 September 2026

[Verification record](leader-verification.json): 468 backend/Integration tests,
145 frontend tests, nine real HTTP/scientific cases, and the compiled Chromium
submission/download/recovery/mobile/CSP checks passed. Browser evidence is in
`.runtime/browser-leader/leader-final/`; the application candidate is
`.runtime/releases/11060f5a58f3f12dbef18a42ca349672a7b2f236f070087e34a1cccac5b2d774`.

The six-tool candidate is staged, not selected over the concurrent seven-tool
update. The current project lock additionally includes CRISPRrepeat. The next
combined release should be built and gated from that complete lock, preserving
this leader integration. This task did not change the running service selection
or its job data.
