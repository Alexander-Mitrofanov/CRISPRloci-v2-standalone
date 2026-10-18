"""Offline command line; jobs cannot download or alter scientific references."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

from .report import FIELDS, document, fasta, read_regular, sha

ROOT = Path(__file__).resolve().parent / "native"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Source-bound ISEScan insertion-sequence context"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    parser.add_argument("--version", action="version", version="ISEScan-v2 2.0.0")
    args = parser.parse_args(argv)
    if not 1 <= args.threads <= 8 or not 1 <= args.timeout_seconds <= 86400:
        parser.error("threads must be 1..8; timeout must be 1..86400 seconds")
    payload = read_regular(args.input, 10_000_000)
    records = fasta(payload)
    output = args.output.absolute()
    if output.exists():
        raise ValueError("Output directory must not exist")
    output.mkdir(parents=True)
    work = output / "work"
    work.mkdir()
    source = work / "source.fasta"
    source.write_text(
        "".join(f">seq{i:06d}\n{seq}\n" for i, seq in enumerate(records.values(), 1)),
        encoding="ascii",
    )
    models = {
        name: sha(read_regular(ROOT / "pHMMs" / name, 50_000_000))
        for name in ("clusters.faa.hmm", "clusters.single.faa")
    }
    # Integration's _run_process already makes this process a session leader.
    # A direct shell invocation gets its own process group so timeout containment
    # never signals the caller's shell. Native descendants remain in this group.
    if os.getpgrp() != os.getpid():
        os.setpgid(0, 0)
    with (work / "native.log").open("wb") as log:
        child = subprocess.Popen(
            [
                sys.executable,
                str(ROOT / "isescan.py"),
                "--seqfile",
                str(source),
                "--output",
                str(work / "native"),
                "--nthread",
                str(args.threads),
            ],
            cwd=work,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            code = child.wait(timeout=args.timeout_seconds)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgrp(), signal.SIGTERM)
            raise RuntimeError("Native ISEScan timed out")
    if code:
        raise RuntimeError(
            f"Native ISEScan failed (exit {code}); see private native log"
        )
    native = work / "native/work/source.fasta.tsv"
    if native.exists():
        raw = read_regular(native)
    else:
        # Upstream does not write a table on legitimate no-hit/no-ORF runs.
        # Require its explicit successful completion marker before normalizing.
        log = read_regular(work / "native.log")
        if b"ISEScan ends at" not in log or not any(
            marker in log
            for marker in (
                b"No hit was returned",
                b"no significant hit",
                b"No IS element was found",
            )
        ):
            raise ValueError(
                "Native output missing without a verified no-hit completion"
            )
        stream = io.StringIO()
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerow(FIELDS)
        raw = stream.getvalue().encode("ascii")
    result = document(
        payload, raw, {"threads": args.threads, "include_partial": True}, models
    )
    # The source table is preserved under a fixed bounded name; paths and raw
    # native diagnostic artifacts remain private under work/.
    (output / "native.tsv").write_bytes(raw)
    temporary = output / ".isescan.json.tmp"
    temporary.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf8"
    )
    temporary.replace(output / "isescan.json")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"ISEScan-v2: {exc}", file=sys.stderr)
        raise SystemExit(1)
