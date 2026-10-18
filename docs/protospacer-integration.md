# Protospacer search and reference provisioning

CRISPRspacer 0.2.0 provides the scientific `database-search` command; Integration
0.6.0 adapts it into the queued `protospacer` workflow, with independent backend
2.2.0 validation and a frontend 2.2.0 results view. No scientific package is
imported into API request handlers. The new workflow searches uploaded target
DNA against an operator-selected, immutable spacer collection.

The NCBI source is the **NCBI-hosted Shmakov et al. 2017 spacerome snapshot**.
It is neither the complete current NCBI collection nor a virus reference
database. The [article data availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605939/)
documents the [NCBI source directory](https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/).
Download/install instructions live with the
[standalone tool](https://github.com/Alexander-Mitrofanov/CRISPRspacer/blob/main/docs/ncbi-spacerome.md).

## Reference staged on the de.NBI VM

On 2026-09-28, the full source was downloaded on `crispr-evor-web-server` only.
Local development and the tool repository retain three original sample records.

| Item | Value |
| --- | --- |
| Source directory | `/srv/crispr/catalogs/ncbi-spacerome-2017` |
| Final immutable index | `/srv/crispr/catalogs/ncbi-spacerome-2017/index-v0.2.0` |
| Source release | 2017-04-25 |
| Source file SHA-256 | `f92ae812231acb1b59035456379245bc824e402ef66f448e26f0932c0676d83b` |
| SQLite SHA-256 | `cedbad8bad4cda1cd4202e2a1708d00af7ea1ff7f5fa5b91cfdb10a3217f9998` |
| SQLite bytes | 389,877,760 |
| Source occurrences | 720,391 |
| Eligible occurrences | 718,792 |
| Excluded by length | 342 |
| Excluded by ambiguity after length filtering | 1,257 |

Import and semantic validation took about 27 seconds and 9.5 seconds respectively
on this VM. A synthetic 1,000,000-base target with the three sample spacers
inserted in both orientations recovered all six exact matches in 3.87 seconds;
deterministic result replay passed. This checks execution, not biological accuracy
or worst-case throughput. Repetitive inputs can hit explicit budgets.

Reference receipts and measurements are retained alongside the VM data.
Installing this reference does not alter the running API or worker release.

## Activate in an explicitly requested web release

1. Restore the eight pinned tool checkouts with `scripts/sync_tools.py`; provision
   isolated runtimes with `scripts/provision_tools.py`. CRISPRspacer needs only a
   separate Python environment, without external alignment software.
2. Build backend/Integration wheels using `scripts/manage_release.py build`.
   Scientific manifest schema `5.0.0` includes the `crisprspacer` source identity.
3. Set both operator environment values in the private deployment configuration:

   ```dotenv
   CRISPR_WEB_SPACER_DATABASE=/srv/crispr/catalogs/ncbi-spacerome-2017/index-v0.2.0
   CRISPR_WEB_SPACER_DATABASE_SHA256=cedbad8bad4cda1cd4202e2a1708d00af7ea1ff7f5fa5b91cfdb10a3217f9998
   ```

4. Run the staged release gate. It builds a separate index from exactly three
   bundled sample records, invokes the real queued pipeline, and checks positive,
   reverse-complement and no-match results and authenticated downloads.
5. Use the existing explicit activation/deployment procedure. Preserve the
   worker's no-network and read-only scientific/runtime/reference filesystem
   boundaries. Publish only `frontend/` to the public site.

The API advertises protospacer mode only when the operator reference is configured
and its scientific release includes CRISPRspacer. Each submitted job stores the
expected database digest in the durable queue; config drift fails before launch.
The scientific process verifies actual database bytes against the same digest.
Requests never supply reference paths, commands or download URLs. Reference
updates require a new directory and digest, with queued jobs drained or cancelled
before changing the configuration. No anonymous or running job updates data.

## Search contract

Database mode uses full-length exact comparison, both strands, ACGT spacers of
16–1,000 bases, linear targets and 1-based inclusive coordinates. Unknown PAM,
array overlap, confidence and E-value remain explicit. Uploaded source arrays
can therefore match; the result is candidate sequence evidence, not a functional
targeting or host assignment. The standalone small-input `search` and `self`
commands retain their existing substitution and circular-target options.

Integration verifies the tool report and input bindings. Backend validation
independently reconstructs every reported target interval, strand and sequence,
checks digests/counts/limits, and regenerates public JSON/TSV exports. Public
summaries cap rows/bytes; full bounded matches remain downloadable. Native
database paths, CLI commands, reference files and internal input bundles are
excluded from publication. TSV export escapes spreadsheet-executable prefixes.

The frontend renders distinct missing, incomplete, unassessable and zero-match
states, a bounded preview, snapshot provenance and limitations. The form accepts
target DNA; it does not let visitors provision or change references.

## Verification of this change

Standalone checks: 92 tests; GitHub CI passed on Python 3.10, 3.12 and 3.13.
`scripts/check.sh backend`: 530 passed; its separately opt-in installed-wheel
acceptance also passed with the final Integration and CRISPRspacer wheels.
`scripts/check.sh frontend`: 180 tests and build/source/public-export checks
passed. Thirteen compiled-browser checks passed using intercepted synthetic API
responses; actual scientific execution was checked separately through the real
queued worker and installed commands. No production service was activated.

The standalone VM executable is
`/srv/crispr/software/CRISPRspacer-0.2.0/environment/bin/crisprspacer`.
The reference and installation are root-owned and read-only for `crisprsvc`.
