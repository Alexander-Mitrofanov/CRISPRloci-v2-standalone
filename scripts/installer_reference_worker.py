#!/usr/bin/env python3
"""Installer-only progress adapters; run inside one scientific environment."""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys

from installer.download import download
from runtime.progress import emit


def spacer_index(arguments):
    from crisprspacer import database
    from crisprspacer.cli import main

    path = Path(arguments[arguments.index("--spacers") + 1])
    original = database._records
    emit("Counting and validating source spacer records")
    total = sum(1 for _ in original(path))

    def records(source):
        for position, record in enumerate(original(source), 1):
            if position == 1 or position % 1000 == 0 or position == total:
                emit("Building the spacer database", position, total, "records")
            yield record
        emit("Finalizing and verifying the spacer database")

    database._records = records
    return main(["database-index", *arguments])


def viral(arguments, operation):
    from crisprspacer.cli import main
    if operation == "viral-download":
        from crisprspacer import ncbi_viral
        original = ncbi_viral._open

        class Response:
            def __init__(self, response, url):
                self.response = response
                self.label = url.rsplit("/", 1)[-1]
                raw = response.headers.get("Content-Length", "")
                self.total = int(raw) if raw.isdecimal() else None
                self.completed = 0

            def __getattr__(self, name):
                return getattr(self.response, name)

            def __enter__(self):
                return self

            def __exit__(self, *args):
                self.response.close()

            def read(self, size=-1):
                block = self.response.read(size)
                self.completed += len(block)
                emit("Downloading " + self.label, self.completed, self.total, "bytes")
                return block

        ncbi_viral._open = lambda url: Response(original(url), url)
    else:
        emit("Building the viral BLAST index; native indexing has no byte total")
    return main([operation, *arguments])


def historical(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(arguments)
    spec = importlib.util.spec_from_file_location("legacy_import", args.source / "scripts/fetch_legacy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cache = args.output.parent / "historical-downloads"
    for name, asset in module.ASSETS.items():
        download(asset["url"], cache / name, asset["sha256"], emit, limit=1024**2)
    return module.main(["--output-dir", str(args.output), "--cache-dir", str(cache)])


if __name__ == "__main__":
    operation, arguments = sys.argv[1], sys.argv[2:]
    if operation == "spacer-index":
        raise SystemExit(spacer_index(arguments))
    if operation in {"viral-download", "viral-index"}:
        raise SystemExit(viral(arguments, operation))
    if operation == "historical":
        raise SystemExit(historical(arguments))
    raise SystemExit("Unsupported installer reference operation")
