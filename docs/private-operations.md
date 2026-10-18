# Private connection records and publication safeguards

Record SSH aliases, VM names, management addresses, usernames, key **paths**,
jump-host chains, verified host-key fingerprints and dated connection steps only
in `.ops-private/`. The local `README.md` indexes the recovered records. Never
copy private key contents, passwords, bearer capabilities or recovery URLs into
notes. Keep keys in the operator's SSH agent/key store.

The directory is ignored at every nesting level. Its local permissions are 0700;
records are 0600. Back it up using an approved private encrypted backup, separately
from Git. Ignoring a file does not remove an existing tracked copy or Git history.

Before and after each connection, append an entry to the local `connections.md`:
UTC date, purpose, selected VM and SSH alias, jump route, exact command without
credentials, host-key verification source, result, services touched and evidence
paths. Distinguish recovered historical evidence from commands actually executed.
Do not weaken host-key checking. Add any new literal VM name, alias or management
address to `.ops-private/forbidden-values.json` (a JSON array of strings).

## Checks before publication

```bash
python3 scripts/install_git_guards.py
python3 scripts/check_private_data.py workspace
python3 scripts/check_private_data.py index
```

Setup installs `.githooks/pre-commit` and `pre-push` through local `core.hooksPath`.
Source archives can still be set up, with an explicit warning that Git guards
were not installed; install them when moving the sources into a real checkout.
The installer refuses to overwrite a different hooks policy. Pre-commit checks
the actual index, including force-added ignored files. Pre-push checks every
outgoing commit tree, including private files deleted in a later commit; new
references scan their entire reachable history. Both fail closed on inspection
errors. Checks report filenames and categories, never matched secrets.

`scripts/check.sh` checks the publishable working tree; when Git metadata is
available it also includes already tracked ignored files. Without Git metadata
it checks the source archive using its ignore rules. That is not a history audit.
The frontend exporter checks source and destination using the same policy before
the user commits or pushes. Only `frontend/` may be exported publicly.

The local signature register detects known connection details copied into any
file, including ordinary Markdown. Portable checks additionally reject private
paths, SSH configuration, key material and GitHub tokens. New arbitrary VM names
cannot be recognized without registration. A public HTTPS API origin necessarily
appears in the deployed browser bundle: use an approved public DNS name if it
must conceal the VM name. Keep the real origin out of source examples and supply
it through deployment configuration.

Hooks are local safeguards; `--no-verify`, a changed hooks configuration, another
clone or a direct API upload can bypass them. Install them in every authoring
checkout and do not bypass failures. CI provides a secondary source/index check;
it runs after upload and cannot undo disclosure. Before the first push from an
existing repository, use the pre-push check with the branch tip and an all-zero
remote object to audit all reachable history. If it fails, remove the data from
outgoing history before publishing; rotate any exposed credentials. Do not
rewrite shared history without coordinating with its owners.

This workspace's Git metadata was inaccessible during the 2026-09-28 audit.
The guards were tested in disposable Git repositories. Installation in the real
checkout and inspection of its previous history remain operator steps; neither
is claimed as complete by a successful source scan.
