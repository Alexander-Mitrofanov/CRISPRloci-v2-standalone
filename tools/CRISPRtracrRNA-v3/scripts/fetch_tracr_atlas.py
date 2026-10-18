"""Stream the public 5.27 GB atlas; retain only original RNA records and provenance.

Requires ijson. Data remain CC BY-NC 4.0, separate from the software license.
"""

import hashlib
import json
import time
import urllib.request
from pathlib import Path

import ijson

URL = "https://storage.googleapis.com/crispr-cas-atlas-xy7q13lmk9/crispr-cas-atlas-v1.0.json"
EXPECTED_SHA256 = "5b4ba2fb99638d279e0c126100e19a4b77aba487b37b7df118e4bf4acd494720"
OUT = Path(__file__).resolve().parents[1] / "benchmarks/literature/raw"


class Reader:
    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()
        self.size = 0

    def read(self, size=-1):
        block = self.stream.read(size)
        self.digest.update(block)
        self.size += len(block)
        return block


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with urllib.request.urlopen(URL, timeout=90) as response:
        headers = dict(response.headers)
        reader = Reader(response)
        count = 0
        retained = 0
        temporary = OUT / "atlas_tracr.jsonl.partial"
        with temporary.open("w") as out:
            for entry in ijson.items(reader, "item", use_float=True):
                count += 1
                tracr = entry.get("tracr")
                if tracr and tracr.get("tracr"):
                    record = {
                        k: entry.get(k) for k in ("operon_id", "summary", "metadata", "tracr")
                    }
                    record["repeats"] = [c.get("crispr_repeat") for c in entry.get("crispr", [])]
                    out.write(json.dumps(record) + "\n")
                    retained += 1
                if count % 25000 == 0:
                    print(
                        f"{count} operons; {retained} tracrRNAs; "
                        f"{reader.size / 1e9:.2f} GB; {time.monotonic() - start:.0f}s",
                        flush=True,
                    )
        if reader.size != int(headers["Content-Length"]):
            raise ValueError("Truncated atlas response")
        if reader.digest.hexdigest() != EXPECTED_SHA256:
            raise ValueError("Atlas source has changed; refusing to replace frozen extraction")
        temporary.replace(OUT / "atlas_tracr.jsonl")
        manifest = {
            "url": URL,
            "retrieved": "2026-09-26",
            "headers": headers,
            "download_sha256": reader.digest.hexdigest(),
            "download_bytes": reader.size,
            "operons": count,
            "retained_tracr_records": retained,
            "extraction": (
                "Keep original operon_id, summary, metadata, tracr, "
                "and CRISPR repeat sequences only"
            ),
            "license": "CC BY-NC 4.0",
            "elapsed_seconds": time.monotonic() - start,
        }
        (OUT / "atlas_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
