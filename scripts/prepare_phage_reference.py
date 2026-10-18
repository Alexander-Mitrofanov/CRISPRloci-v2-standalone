#!/usr/bin/env python3
"""Operator-only, resumable NCBI snapshot acquisition and offline genome selection.

This command never changes an active reference or downloads data for a web job.
The snapshot includes all search results, then explicitly accounts for every
record as selected, duplicate or excluded. It is not all phages in nature.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
COMMON = 'txid10239[Organism:exp] AND "complete genome"[Title]'
HOST_QUERY = COMMON + ' AND ("vhost bacteria"[Filter] OR "vhost archaea"[Filter])'
QUERY = COMMON + ' AND (phage[Title] OR bacteriophage[Title] OR "vhost bacteria"[Filter] OR "vhost archaea"[Filter])'
ACCESSION = re.compile(r"[A-Z][A-Z0-9_]{1,24}\.[0-9]{1,5}")
MAX_CANDIDATES = 100_000
BATCH_SIZE = 250
MAX_METADATA_BYTES = 8 * 1024 * 1024
MAX_BATCH_BYTES = 512 * 1024 * 1024
MIN_FREE_BYTES = 20 * 1024**3
REQUEST_INTERVAL = 0.55
_request_lock = threading.Lock()
_next_request_at = 0.0
POLICY = (
    "NCBI Nucleotide complete-genome-title viral records with an explicit phage/bacteriophage title "
    "or NCBI bacterial/archaeal vhost annotation; explicit viral taxonomy, genomic DNA, one unsegmented record, 1,000–5,000,000 "
    "bases, ACGTN with at least one called base; exclude satellite/vector/partial annotations. "
    "Exact sequence duplicates collapse with RefSeq preferred, retaining all accession aliases. "
    "Isolation-host annotations are retained, not used to assert infection. No sampling."
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_json(path, value):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    os.replace(temp, path)


def read_json(path):
    with Path(path).open("rb") as stream:
        payload = stream.read(MAX_METADATA_BYTES + 1)
    if len(payload) > MAX_METADATA_BYTES:
        raise ValueError("Snapshot metadata exceeds bounded size")
    return json.loads(payload)


def require_space(destination, additional_bytes=0):
    if shutil.disk_usage(destination).free < MIN_FREE_BYTES + additional_bytes:
        raise ValueError("Reference preparation would consume the 20 GiB deployment reserve")


def search_ids(result):
    ids = result["idlist"]
    if result.get("errorlist") or not 0 < int(result["count"]) <= MAX_CANDIDATES:
        raise ValueError("Search is empty, invalid or exceeds supported capacity")
    if (
        not isinstance(ids, list)
        or len(ids) != int(result["count"])
        or any(not isinstance(uid, str) or not re.fullmatch(r"[0-9]+", uid) for uid in ids)
        or len(set(ids)) != len(ids)
    ):
        raise ValueError("Search did not enumerate every unique candidate")
    return ids


def frozen_searches(destination, snapshot):
    names = {"candidates.json", "host-qualified.json"}
    if (
        snapshot.get("schema_version") != "crisprloci.ncbi-phage-snapshot/1.0.0"
        or snapshot.get("query") != QUERY
        or snapshot.get("host_query") != HOST_QUERY
        or snapshot.get("batch_size") != BATCH_SIZE
        or set(snapshot.get("search_sha256", {})) != names
    ):
        raise ValueError("Existing snapshot uses another acquisition contract")
    sources = destination / "sources"
    for name, digest in snapshot["search_sha256"].items():
        if sha256(sources / name) != digest:
            raise ValueError("Frozen search receipt changed")
    candidates = search_ids(read_json(sources / "candidates.json"))
    host = search_ids(read_json(sources / "host-qualified.json"))
    if snapshot.get("candidate_count") != len(candidates) or not set(host) <= set(candidates):
        raise ValueError("Frozen search coverage changed")
    return candidates, set(host)


def wait_for_request_slot():
    """Share one request-start budget across search and all download threads."""
    global _next_request_at
    with _request_lock:
        delay = _next_request_at - time.monotonic()
        if delay > 0:
            time.sleep(delay)
        _next_request_at = time.monotonic() + REQUEST_INTERVAL


def request(endpoint, parameters):
    """At most two requests/sec across workers, with bounded transient retries."""
    payload = urllib.parse.urlencode({"tool": "crisprloci-reference-builder", **parameters}).encode()
    for attempt in range(6):
        if attempt:
            time.sleep(min(60, 2 ** attempt))
        wait_for_request_slot()
        try:
            response = urllib.request.urlopen(
                urllib.request.Request(BASE + endpoint, data=payload), timeout=180
            )
            host = urllib.parse.urlparse(response.url).hostname or ""
            if not (host == "ncbi.nlm.nih.gov" or host.endswith(".ncbi.nlm.nih.gov")):
                response.close()
                raise ValueError("Unexpected NCBI redirect")
            return response
        except (urllib.error.URLError, TimeoutError):
            if attempt == 5:
                raise
    raise RuntimeError("Unreachable request state")


def search(query):
    with request("esearch.fcgi", {"db": "nuccore", "term": query, "retmode": "json", "retmax": MAX_CANDIDATES}) as response:
        payload = response.read(MAX_METADATA_BYTES + 1)
    if len(payload) > MAX_METADATA_BYTES:
        raise ValueError("Search response exceeds metadata bound")
    result = json.loads(payload)["esearchresult"]
    search_ids(result)
    return result


def acquire_batch(destination, sources, selected_ids, start, download_workers):
    """Accept one distinct batch only after hashes and exact UID coverage pass."""
    path = sources / f"batch-{start:06d}.xml.gz"
    receipt_path = path.with_suffix(".json")
    if receipt_path.exists():
        receipt = read_json(receipt_path)
        if receipt["ids"] != selected_ids or sha256(path) != receipt["sha256"]:
            raise ValueError("Completed acquisition batch changed")
        return len(selected_ids)
    # Reserve room for the bounded in-flight batches as well as deployment data.
    require_space(destination, download_workers * (MAX_BATCH_BYTES + 1024 * 1024))
    temporary = path.with_suffix(".part")
    size = 0
    content_digest = hashlib.sha256()
    with request("efetch.fcgi", {"db": "nuccore", "id": ",".join(selected_ids), "rettype": "gbwithparts", "retmode": "xml"}) as response:
        with temporary.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            for block in iter(lambda: response.read(1024 * 1024), b""):
                size += len(block)
                if size > MAX_BATCH_BYTES:
                    raise ValueError("NCBI batch exceeds bounded response size")
                content_digest.update(block)
                compressed.write(block)
    os.replace(temporary, path)
    # Parse before accepting a batch; an HTML/error/truncated response is not a receipt.
    uids = []
    for record in iter_records(path):
        uids.append(record["uid"])
        require_viral_sequence(record)
    if len(uids) != len(selected_ids) or set(uids) != set(selected_ids):
        raise ValueError("NCBI batch omitted or substituted requested records")
    save_json(receipt_path, {"ids": selected_ids, "sha256": sha256(path), "content_sha256": content_digest.hexdigest(), "content_bytes": size,
                             "rettype": "gbwithparts", "retmode": "xml"})
    return len(selected_ids)


def acquire(destination, download_workers=1):
    if type(download_workers) is not int or not 1 <= download_workers <= 3:
        raise ValueError("Download workers must be between one and three")
    destination.mkdir(parents=True, exist_ok=True)
    sources = destination / "sources"
    sources.mkdir(exist_ok=True)
    snapshot_path = destination / "snapshot.json"
    if snapshot_path.exists():
        snapshot = read_json(snapshot_path)
    else:
        candidates, host = search(QUERY), search(HOST_QUERY)
        if not set(host["idlist"]) <= set(candidates["idlist"]):
            raise ValueError("NCBI changed during query capture; retry in a fresh snapshot")
        save_json(sources / "candidates.json", candidates)
        save_json(sources / "host-qualified.json", host)
        snapshot = {
            "schema_version": "crisprloci.ncbi-phage-snapshot/1.0.0",
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "query": QUERY, "host_query": HOST_QUERY, "endpoint": BASE,
            "batch_size": BATCH_SIZE,
            "candidate_count": len(candidates["idlist"]),
            "search_sha256": {name: sha256(sources / name) for name in ("candidates.json", "host-qualified.json")},
        }
        save_json(snapshot_path, snapshot)
    # A failed verification/resume must not leave an old completion marker behind.
    if snapshot.pop("acquisition_completed_at", None) is not None:
        save_json(snapshot_path, snapshot)
    candidates, _ = frozen_searches(destination, snapshot)
    completed = 0
    with ThreadPoolExecutor(max_workers=download_workers, thread_name_prefix="ncbi-download") as executor:
        futures = [executor.submit(acquire_batch, destination, sources, candidates[start:start + BATCH_SIZE], start, download_workers)
                   for start in range(0, len(candidates), BATCH_SIZE)]
        try:
            for future in as_completed(futures):
                completed += future.result()
                print(f"Acquired or verified {completed}/{len(candidates)} records", flush=True)
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    snapshot["acquisition_completed_at"] = datetime.now(timezone.utc).isoformat()
    save_json(snapshot_path, snapshot)
    return snapshot


def parse_record(node):
    def value(name):
        return node.findtext("GBSeq_" + name, default="")
    accession = value("accession-version")
    if not ACCESSION.fullmatch(accession):
        raise ValueError("Missing or invalid versioned nucleotide accession")
    identifiers = [child.text or "" for child in node.findall("./GBSeq_other-seqids/GBSeqid")]
    uids = [ident[3:] for ident in identifiers if re.fullmatch(r"gi\|[0-9]+", ident)]
    if len(uids) != 1:
        raise ValueError("NCBI record is not bound to exactly one requested nucleotide UID")
    source = [feature for feature in node.findall("./GBSeq_feature-table/GBFeature") if feature.findtext("GBFeature_key") == "source"]
    qualifiers = {}
    for feature in source:
        for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
            qualifiers.setdefault(qualifier.findtext("GBQualifier_name"), []).append(qualifier.findtext("GBQualifier_value", default=""))
    sequence = value("sequence").upper()
    declared_length = int(value("length"))
    contig = value("contig")
    if declared_length < 1 or (declared_length != len(sequence) and not (not sequence and contig)):
        raise ValueError("NCBI sequence length does not match its record")
    return {
        "uid": uids[0], "accession": accession, "title": value("definition"),
        "organism": value("organism"), "taxonomy": value("taxonomy"),
        "taxids": [item[6:] for item in qualifiers.get("db_xref", []) if item.startswith("taxon:")],
        "moltype": value("moltype"), "mol_type": qualifiers.get("mol_type", []),
        "segment": qualifiers.get("segment", []), "source_feature_count": len(source),
        "isolation_host_annotation": qualifiers.get("host", []),
        "isolation_source": qualifiers.get("isolation_source", []),
        "source_database": "RefSeq" if "_" in accession else "GenBank",
        "update_date": value("update-date"), "length": declared_length,
        "sequence_available": bool(sequence), "contig": contig or None,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest() if sequence else None,
        "sequence": sequence,
    }


def viral_taxonomy(record):
    return record["taxonomy"].split(";", 1)[0].strip() == "Viruses"


def require_viral_sequence(record):
    if viral_taxonomy(record) and not record["sequence_available"]:
        raise ValueError("Viral CONTIG record requires expanded NCBI sequence retrieval: " + record["accession"])


def iter_records(path):
    class BoundedXML:
        def __init__(self, stream):
            self.stream = stream
            self.size = 0

        def read(self, size=-1):
            remaining = MAX_BATCH_BYTES - self.size
            data = self.stream.read(min(size, remaining + 1) if size >= 0 else remaining + 1)
            self.size += len(data)
            if self.size > MAX_BATCH_BYTES:
                raise ValueError("NCBI XML batch exceeds bounded decompressed size")
            return data

    with gzip.open(path, "rb") as stream:
        root = None
        count = 0
        for event, element in ET.iterparse(BoundedXML(stream), events=("start", "end")):
            if root is None:
                root = element
                if root.tag != "GBSet":
                    raise ValueError("Expected an NCBI GenBank XML document")
            if event == "end" and element.tag == "GBSeq":
                count += 1
                if count > BATCH_SIZE:
                    raise ValueError("NCBI batch contains too many records")
                yield parse_record(element)
                root.clear()
        if count == 0:
            raise ValueError("NCBI batch contains no records")


def exclusion(record, host_ids):
    if not viral_taxonomy(record):
        return "not_viral_taxonomy"
    require_viral_sequence(record)
    title = record["title"].lower()
    if "complete genome" not in title:
        return "not_complete_genome_title"
    if record["uid"] not in host_ids and not re.search(r"\b(?:bacterio)?phage\b", title):
        return "no_phage_or_prokaryotic_vhost_evidence"
    if record["moltype"] != "DNA" or record["mol_type"] != ["genomic DNA"]:
        return "not_single_genomic_DNA_record"
    if record["source_feature_count"] != 1 or record["segment"] or re.search(r"\bsegment\b", title):
        return "segmented_or_multiple_source_features"
    if re.search(r"\b(?:satellite|vector|partial|synthetic)\b", title):
        return "satellite_vector_partial_or_synthetic"
    if not 1000 <= record["length"] <= 5_000_000:
        return "length_outside_supported_range"
    if set(record["sequence"]) - set("ACGTN"):
        return "unsupported_sequence_alphabet"
    if not set(record["sequence"]) & set("ACGT"):
        return "no_called_bases"
    return None


def select(snapshot_root):
    snapshot = read_json(snapshot_root / "snapshot.json")
    if not snapshot.get("acquisition_completed_at"):
        raise ValueError("Cannot select an incomplete NCBI acquisition")
    output = snapshot_root / "selected"
    if output.exists():
        raise ValueError("Selected snapshot already exists; never overwrite an immutable selection")
    sources = snapshot_root / "sources"
    ids, host_ids = frozen_searches(snapshot_root, snapshot)
    temporary = Path(tempfile.mkdtemp(prefix=".selection-", dir=snapshot_root))
    (temporary / "genomes").mkdir()
    db = sqlite3.connect(temporary / "selection.sqlite3")
    db.execute("CREATE TABLE records(accession TEXT PRIMARY KEY,uid TEXT UNIQUE,status TEXT,reason TEXT,representative TEXT,sequence_sha256 TEXT,metadata_json TEXT)")
    db.execute("CREATE TABLE candidates(accession TEXT PRIMARY KEY,priority INTEGER,sha256 TEXT,length INTEGER,file TEXT)")
    counts = Counter()
    source_receipts = []
    try:
        for start in range(0, len(ids), BATCH_SIZE):
            path = sources / f"batch-{start:06d}.xml.gz"
            receipt = read_json(path.with_suffix(".json"))
            if receipt["ids"] != ids[start:start + BATCH_SIZE] or sha256(path) != receipt["sha256"]:
                raise ValueError("Frozen sequence acquisition batch changed")
            source_receipts.append({"file": path.name, "sha256": receipt["sha256"]})
            seen = set()
            for record in iter_records(path):
                if record["uid"] not in receipt["ids"] or record["uid"] in seen:
                    raise ValueError("Acquisition identity mismatch")
                seen.add(record["uid"])
                reason = exclusion(record, host_ids)
                sequence = record.pop("sequence")
                record["inclusion_evidence"] = (["ncbi_prokaryotic_vhost_filter"] if record["uid"] in host_ids else []) + (["explicit_phage_title"] if re.search(r"\b(?:bacterio)?phage\b", record["title"], re.I) else [])
                metadata_json = json.dumps(record, sort_keys=True)
                fasta_bytes = len(sequence) + (len(sequence) + 79) // 80 + len(record["accession"]) + 2
                require_space(snapshot_root, 2 * len(metadata_json.encode()) + 8192 + (fasta_bytes if not reason else 0))
                db.execute("INSERT INTO records VALUES(?,?,?,?,?,?,?)", (record["accession"], record["uid"], "excluded" if reason else "eligible", reason, None, record["sequence_sha256"], metadata_json))
                if reason:
                    counts[reason] += 1
                    continue
                file = "genomes/" + record["accession"] + ".fasta"
                (temporary / file).write_text(">" + record["accession"] + "\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
                db.execute("INSERT INTO candidates VALUES(?,?,?,?,?)", (record["accession"], 0 if record["source_database"] == "RefSeq" else 1, record["sequence_sha256"], record["length"], file))
            if seen != set(receipt["ids"]):
                raise ValueError("Acquisition coverage mismatch")
            db.commit()
        total = db.execute("SELECT count(*) FROM records").fetchone()[0]
        if total != len(ids):
            raise ValueError("Selection failed to account for every NCBI record")
        representative = {}
        files = []
        bases = 0
        for accession, _, digest, length, file in db.execute("SELECT * FROM candidates ORDER BY priority,accession").fetchall():
            if digest in representative:
                db.execute("UPDATE records SET status='duplicate',representative=? WHERE accession=?", (representative[digest], accession))
                (temporary / file).unlink()
                continue
            representative[digest] = accession
            files.append(file)
            bases += length
            db.execute("UPDATE records SET status='selected',representative=? WHERE accession=?", (accession, accession))
        if not 1 <= len(files) <= MAX_CANDIDATES or bases > 10_000_000_000:
            raise ValueError("Selected reference exceeds the supported native database capacity")
        db.execute("DROP TABLE candidates")
        db.commit()
        statuses = dict(db.execute("SELECT status,count(*) FROM records GROUP BY status"))
        with (temporary / "accessions.tsv").open("w") as stream:
            stream.write("accession\tstatus\trepresentative\tsequence_sha256\treason\n")
            for row in db.execute("SELECT accession,status,representative,sequence_sha256,reason FROM records ORDER BY accession"):
                stream.write("\t".join(value or "" for value in row) + "\n")
        db.close()
        receipt = {
            "schema_version": "crisprloci.phage-selection/1.0.0", "policy": POLICY,
            "snapshot": snapshot, "source_files": source_receipts,
            "records": total, "statuses": statuses, "exclusions": dict(counts),
            "genome_count": len(files), "total_bases": bases,
            "selection_database_sha256": sha256(temporary / "selection.sqlite3"),
            "accession_inventory_sha256": sha256(temporary / "accessions.tsv"),
            "limitations": [
                "Complete-genome title and phage/vhost annotations define coverage; records lacking these annotations are absent.",
                "Title-based inclusion is not independently verified host taxonomy; isolation host does not establish viral host range.",
                "RNA, segmented, partial, satellite/vector/synthetic and unsupported records are excluded with reasons.",
                "Exact sequence deduplication is not species clustering or circular-rotation equivalence.",
                "A missing association is not absence of phages or infection; empirical calibration may be unavailable.",
            ],
        }
        save_json(temporary / "selection.json", receipt)
        (temporary / "genomes.txt").write_text("\n".join(str(output / file) for file in files) + "\n")
        save_json(temporary / "source.json", {
            "collection": "NCBI Nucleotide GenBank and RefSeq complete DNA phage/prokaryotic-virus records",
            "release": "snapshot-" + snapshot["retrieved_at"][:10], "retrieved_at": snapshot["retrieved_at"],
            "selection_policy": POLICY + " Selection receipt SHA256: " + sha256(temporary / "selection.json"),
        })
        os.replace(temporary, output)
        print(json.dumps({key: receipt[key] for key in ("genome_count", "total_bases", "statuses", "exclusions")}, indent=2), flush=True)
        return receipt
    finally:
        db.close()
        if temporary.exists():
            shutil.rmtree(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("acquire", "select"))
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--download-workers", type=int, choices=(1, 2, 3), default=1,
                        help="Concurrent acquisition batches, sharing at most two requests/sec (default: 1)")
    args = parser.parse_args()
    root = args.snapshot.expanduser().resolve()
    if args.command == "acquire":
        acquire(root, args.download_workers)
    else:
        select(root)


if __name__ == "__main__":
    main()
