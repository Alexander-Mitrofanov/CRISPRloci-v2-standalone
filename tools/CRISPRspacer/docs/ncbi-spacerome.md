# NCBI-hosted spacer collection

Use the **Shmakov et al. 2017 spacerome snapshot**, hosted by NCBI. It is a
published research collection, not all NCBI spacers or a continuously updated
NCBI database. The [2017 article's data availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605939/)
documents [NCBI's spacerome directory](https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/).
The [2020 article](https://www.nature.com/articles/s42003-020-1014-1) references
the earlier collection and an additional supplement. Its `Data.tar.gz` contains
experimental spacer sets and analysis tables; it is not a replacement general
spacer FASTA. NCBI's directory timestamps the selected file 2017-04-25.

## Full download: run on the data VM only

These are operator commands. Do not put them in API handlers or job workers.
Reserve at least 20 GiB free disk for the existing service. The source is about
13 MiB compressed; the SQLite index is larger. Use a new immutable directory.

```bash
# Run on the VM, with CRISPRspacer 0.2.0 installed in its own environment.
set -euo pipefail
mkdir -p /srv/crispr/catalogs/ncbi-spacerome-2017
cd /srv/crispr/catalogs/ncbi-spacerome-2017
test ! -e Supplementary_material_spacers.fna.gz
curl --fail --location --proto '=https' --tlsv1.2 --max-time 300 \
  'https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz' \
  -o Supplementary_material_spacers.fna.gz.partial
printf '%s  %s\n' \
  f92ae812231acb1b59035456379245bc824e402ef66f448e26f0932c0676d83b \
  Supplementary_material_spacers.fna.gz.partial | sha256sum --check --strict
gzip --test Supplementary_material_spacers.fna.gz.partial
mv Supplementary_material_spacers.fna.gz.partial Supplementary_material_spacers.fna.gz

crisprspacer database-index \
  --spacers Supplementary_material_spacers.fna.gz \
  --name ncbi-spacerome-2017 \
  --source-url 'https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz' \
  --source-release 2017-04-25 --output index
crisprspacer validate-database index
date -u '+%Y-%m-%dT%H:%M:%SZ' > retrieved-at.txt
```

Use `set -euo pipefail` when placing these commands in a script. The SHA-256
above was measured from NCBI HTTPS bytes on 2026-09-28; it is our snapshot pin,
not a checksum independently published by NCBI. If NCBI changes the bytes,
stop and review the new snapshot rather than removing checksum verification.
Record retrieval time, source URL, source checksum, importer version and command
in the deployment receipt. The source and index stay on the VM. Do not commit
them to GitHub, copy them to development machines or include them in job output.
Make the completed directory read-only to the worker account and to untrusted
users. Never replace its files while jobs are using it; build another versioned
directory for an update.

The full source contains **720,391 occurrence records** (measured at import,
not a count of unique sequences or validated arrays). Eligibility is ACGT only,
length 16–1,000 inclusive. The index retains all source records and exclusive
skip counts: length filtering precedes ambiguity filtering. No source IDs are
reinterpreted as validated taxonomy, source coordinates or array annotations.

## Development: exactly three downloaded sample records

Only `examples/ncbi-spacerome-three.fasta` is included here: the first three
records of the pinned source, in original order with original IDs. Its provenance
is in `examples/ncbi-spacerome-three.provenance.json`. It is a tiny mechanics
fixture, not a representative biological benchmark. Existing synthetic examples
remain separate. Do not run the full download above on the development machine.

```bash
crisprspacer database-index \
  --spacers examples/ncbi-spacerome-three.fasta \
  --name ncbi-spacerome-three-sample \
  --source-url 'https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz' \
  --source-release 2017-04-25-three-record-sample --output results/sample-index
crisprspacer database-search --database results/sample-index \
  --targets examples/ncbi-spacerome-three.fasta --output results/sample-search
crisprspacer validate-database-result results/sample-search --database results/sample-index
```

The sample index hashes the sample file. Its digest is intentionally different
from the full downloaded source; the provenance file identifies the parent.

## Standalone and web-server use

```bash
crisprspacer database-search \
  --database /srv/crispr/catalogs/ncbi-spacerome-2017/index \
  --expected-database-sha256 '<database_sha256 from index/reference.json>' \
  --targets submitted-targets.fasta --output results/protospacers
```

The same command is the scientific stage in CRISPRloci's queued `protospacer`
mode. Operators set `CRISPR_WEB_SPACER_DATABASE` to the immutable index directory
and `CRISPR_WEB_SPACER_DATABASE_SHA256` to the SQLite digest; the API stores that
digest with each job. Requests cannot choose reference paths or trigger downloads.
Workers use a fixed launcher with no network access. The standalone package has
no dependency on the web server. Deploying a new web runtime is a separate
operator release step; installing data alone does not activate a web feature.

Database search is **exact only**, full-length, both DNA orientations, and linear
target traversal. A 16-base prefix index nominates candidates; complete sequences
are compared before reporting them. Palindromes appear once as `both`; duplicate
spacer sequences with different occurrence IDs remain separate matches.
The existing small-input `search`/`self` commands still allow 0–2 substitutions
and explicit circular targets. Database search does not silently reuse those
policies and cannot detect circular-origin matches or substitutions.

Results are candidate sequence matches. PAMs and reported-array overlap are not
evaluated; self-array matches can remain. No biological targeting, host assignment,
infection, immunity, confidence or E-value is inferred. A no-match report is
relative to the selected snapshot and exact-search policy. Ambiguous target
windows never match. `seed_windows` counts all 16-base starts and
`ambiguous_seed_windows` counts starts whose prefix contains ambiguity; these
are search diagnostics, not full-spacer biological coverage measurements.

Index construction and `validate-database` check every spacer and its complete
forward/reverse seed inventory before operators pin the snapshot. The
read-only SQLite database is checked against the expected SHA-256 and its
internal metadata on every search. Output contains reference identity, source
target sequences, matched spacer evidence, JSON, TSV and a checksummed manifest;
it never copies the full database. `validate-database-result` replays against
the same separately retained database snapshot. Default limits: 5 million target
bases, 1,000 target records, 10,000 matches, 2 million nominated occurrence
candidates, 16 MiB of match evidence and 300 seconds. Exceeding a limit fails
without publishing a completed result. Shorter deadlines/caps can be requested.
