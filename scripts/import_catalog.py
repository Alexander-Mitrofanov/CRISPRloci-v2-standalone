#!/usr/bin/env python3
"""Build a private, immutable catalogue snapshot from pinned source bundles.

Local mode accepts a bundle and SHA-256. Fleet mode consumes inventories produced
by catalog_source_inventory.py, downloads compact shards over operator SSH, and
publishes only after every selected source has been checksum-verified.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import threading

from catalog_snapshot import COLUMNS, TABLES, read_frame, transform_bundle

CATALOG_SIZE = 3_331_847
SUCCESSFUL = 3_331_684
SOURCE_FEATURE_COUNTS = {"arrays": 9_583_665, "cas_systems": 9_237_875, "cas_genes": 220_117_236}
PROFILE = "5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540"
NORMALIZER_VERSION = "2-feature-context"
LIMITATIONS = [
    "Predicted annotations from CRISPRidentify 2.3.2 region-screen and standalone CasAndra; not experimentally confirmed genes.",
    "Arrays include accepted Bona-fide and Possible calls. Low-score/rejected candidates and strand predictions were not retained by this campaign.",
    "Repeat and spacer counts describe occurrences, not unique sequences. Coordinates are 1-based inclusive on the forward source sequence; missing values remain unknown.",
    "Protein sequences and bulky native reports are not included in this compact catalogue; source checksums identify the preserved upstream exports.",
    "Sequence searches match exact nucleotide strings. No sequence alignment or similarity search is performed.",
    "Contigs without arrays or Cas annotations are omitted from feature shards; assembly sequence counts retain the complete upstream totals.",
]


class SnapshotWriter:
    def __init__(self, directory, snapshot_id):
        self.directory = Path(directory).absolute()
        self.directory.mkdir(parents=True, exist_ok=True)
        for path in [*self.directory.parents, self.directory]:
            if path.is_symlink():
                raise ValueError("snapshot directory cannot contain a symlink")
        if (self.directory / "catalog.sqlite3").exists():
            raise ValueError("snapshot is already published; choose a new directory")
        (self.directory / "shards").mkdir(exist_ok=True)
        (self.directory / "receipts").mkdir(exist_ok=True)
        for name in ("shards", "receipts", "catalog.building.sqlite3"):
            if (self.directory / name).is_symlink():
                raise ValueError("snapshot members cannot be symlinks")
        self.snapshot_id = snapshot_id
        self.lock = threading.Lock()
        self.db = sqlite3.connect(self.directory / "catalog.building.sqlite3", check_same_thread=False)
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.executescript(f"""
            CREATE TABLE IF NOT EXISTS export_metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS assemblies({TABLES['assemblies']});
            CREATE TABLE IF NOT EXISTS shards(id INTEGER PRIMARY KEY,relative_path TEXT NOT NULL,sha256 TEXT NOT NULL,uncompressed_bytes INTEGER NOT NULL,counts_json TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS assembly_shards(assembly_pk INTEGER PRIMARY KEY REFERENCES assemblies(assembly_pk),shard_id INTEGER NOT NULL REFERENCES shards(id));
        """)
        self.db.commit()
        identity = self.db.execute("SELECT value FROM export_metadata WHERE key='snapshot_id'").fetchone()
        if identity is not None and identity[0] != snapshot_id:
            raise ValueError("resume requires the same snapshot identifier")
        self.db.execute("INSERT OR IGNORE INTO export_metadata VALUES ('snapshot_id',?)", (snapshot_id,))
        version = self.db.execute("SELECT value FROM export_metadata WHERE key='normalizer_version'").fetchone()
        if version is not None and version[0] != NORMALIZER_VERSION:
            raise ValueError("resume requires the same normalization version")
        if version is None and self.db.execute("SELECT 1 FROM shards LIMIT 1").fetchone() is not None:
            raise ValueError("unversioned existing shards require an explicit offline migration")
        self.db.execute("INSERT OR IGNORE INTO export_metadata VALUES ('normalizer_version',?)", (NORMALIZER_VERSION,))
        self.db.commit()

    def completed(self):
        return {row[0] for row in self.db.execute("SELECT id FROM shards")}

    def verify_files(self):
        for shard_id, relative, expected in self.db.execute("SELECT id,relative_path,sha256 FROM shards ORDER BY id"):
            path = self.directory / relative
            if relative != f"shards/{shard_id:08d}.sqlite3.gz" or path.is_symlink() or not path.is_file():
                raise ValueError("committed shard path is unsafe or missing")
            with path.open("rb") as handle:
                actual = hashlib.file_digest(handle, "sha256").hexdigest()
            receipt_path = self.directory / "receipts" / f"{shard_id:08d}.json"
            if receipt_path.is_symlink() or not receipt_path.is_file() or receipt_path.stat().st_size > 1024**2:
                raise ValueError("committed shard receipt is unsafe or missing")
            receipt = json.loads(receipt_path.read_text())
            if actual != expected or receipt.get("sha256") != actual:
                raise ValueError("committed shard checksum mismatch")

    def accept(self, receipt, data):
        with self.lock:
            if shutil.disk_usage(self.directory).free - len(data) < 5 * 1024**3:
                raise ValueError("catalogue transfer stopped at the 5 GiB free-space reserve; resume on a larger volume")
            shard_id = receipt["id"]
            if self.db.execute("SELECT 1 FROM shards WHERE id=?", (shard_id,)).fetchone():
                raise ValueError("refusing to replace a committed shard")
            if hashlib.sha256(data).hexdigest() != receipt["sha256"]:
                raise ValueError("shard checksum mismatch")
            path = f"shards/{shard_id:08d}.sqlite3.gz"
            temporary = self.directory / (path + ".tmp")
            if temporary.is_symlink() or (self.directory / path).is_symlink():
                raise ValueError("shard output cannot be a symlink")
            with temporary.open("wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, self.directory / path)
            public_receipt = {key: value for key, value in receipt.items() if key != "assemblies"}
            receipt_path = self.directory / "receipts" / f"{shard_id:08d}.json"
            if receipt_path.is_symlink():
                raise ValueError("receipt output cannot be a symlink")
            with receipt_path.open("w") as handle:
                handle.write(json.dumps(public_receipt, sort_keys=True) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            for folder in ("shards", "receipts"):
                descriptor = os.open(self.directory / folder, os.O_RDONLY | os.O_DIRECTORY)
                try:
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
            with self.db:
                self.db.execute("INSERT INTO shards VALUES (?,?,?,?,?)", (shard_id, path, receipt["sha256"], receipt["uncompressed_bytes"], json.dumps(receipt["counts"])))
                self.db.executemany("INSERT INTO assemblies VALUES (" + ",".join("?" for _ in COLUMNS["assemblies"]) + ")", receipt["assemblies"])
                self.db.executemany("INSERT INTO assembly_shards VALUES (?,?)", ((row[0], shard_id) for row in receipt["assemblies"]))

    def publish(self, expected=None):
        self.verify_files()
        counts = {name: 0 for name in TABLES}
        identity = hashlib.sha256()
        for shard_id, digest, raw in self.db.execute("SELECT id,sha256,counts_json FROM shards ORDER BY id"):
            identity.update(f"{shard_id}:{digest}\n".encode())
            for key, value in json.loads(raw).items():
                counts[key] += value
        if expected is not None and counts["assemblies"] != expected:
            raise ValueError(f"catalogue coverage mismatch: {counts['assemblies']} != {expected}")
        statuses = dict(self.db.execute("SELECT status,count(*) FROM assemblies GROUP BY status"))
        if sum(statuses.values()) != counts["assemblies"] or self.db.execute("SELECT count(*) FROM assembly_shards").fetchone()[0] != counts["assemblies"]:
            raise ValueError("registry and shard counts disagree")
        if expected == CATALOG_SIZE and statuses != {"complete": SUCCESSFUL, "failed": 151, "unavailable": 12}:
            raise ValueError("final outcome counts differ from the completed campaign")
        if expected == CATALOG_SIZE and any(counts[key] != value for key, value in SOURCE_FEATURE_COUNTS.items()):
            raise ValueError("feature counts differ from the independent selected source summaries")
        if self.db.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("registry relations are invalid")
        summary = {
            "available": True, "snapshot_id": self.snapshot_id,
            "title": "CRISPRcasDb · bacterial and archaeal reference catalogue",
            "source_url": None, "generated_at": datetime.now(timezone.utc).isoformat(),
            "coordinate_system": "1-based inclusive; forward-source coordinates",
            "counts": {key: counts[key] for key in ("assemblies", "arrays", "repeats", "spacers", "cas_genes", "cas_systems")},
            "coverage": {"catalog_assemblies": CATALOG_SIZE, "successful_assemblies": SUCCESSFUL, "failed_assemblies": 151, "unavailable_assemblies": 12, "imported_assemblies": counts["assemblies"]},
            "limitations": LIMITATIONS + (["This is a partial imported snapshot; counts above describe only the imported records."] if counts["assemblies"] != CATALOG_SIZE else []),
            "provenance": {"source_project": "ClusterCRISPRCasDb / CRISPR Ensemble 0.6.2", "catalog_id": "c26a75991f66d0c8", "catalog_snapshot": "2026-09-02", "analysis_profile_sha256": PROFILE, "campaign_completed_at": "2026-09-22T05:41:03.397420Z", "campaign_verified_at": "2026-09-22T09:51:18.075175Z", "shards_sha256": identity.hexdigest(), "storage": "compressed-sqlite-shards-v1", "imported_status_counts": json.dumps(statuses, sort_keys=True)},
        }
        self.db.execute("INSERT OR REPLACE INTO export_metadata VALUES ('summary',?)", (json.dumps(summary),))
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_assembly_shards_shard ON assembly_shards(shard_id)")
        self.db.commit()
        if self.db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("catalogue integrity check failed")
        self.db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.close()
        manifest_path = self.directory / "manifest.json"
        if manifest_path.is_symlink():
            raise ValueError("manifest cannot be a symlink")
        with manifest_path.open("w") as handle:
            handle.write(json.dumps(summary, indent=2) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        with (self.directory / "catalog.building.sqlite3").open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(self.directory / "catalog.building.sqlite3", self.directory / "catalog.sqlite3")
        descriptor = os.open(self.directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return summary


def ssh_command(config, node):
    return ["ssh", "-F", str(Path(config).absolute()), "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "-o", "ControlMaster=no", "-o", "ServerAliveInterval=15", "-o", "ServerAliveCountMax=3", f"crispr-node-{node:02d}", "python3 -B -"]


def verify_receipt(job, receipt):
    if receipt["id"] != job["id"]:
        raise ValueError("received shard is not the requested source")
    if job.get("kind") == "raw_group":
        documents = [{"sha256": item["sha256"], "kind": item["kind"], "accessions": [entry["accession"] for entry in item["entries"]]} for item in job["items"]]
        expected_digest = hashlib.sha256(json.dumps(documents, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        expected_accessions = {accession for document in documents for accession in document["accessions"]}
    else:
        expected_digest = job["sha256"]
        expected_accessions = set(job["selected"])
    actual_accessions = [row[1] for row in receipt["assemblies"]]
    if (receipt["source_sha256"] != expected_digest or set(actual_accessions) != expected_accessions
            or len(actual_accessions) != len(expected_accessions)):
        raise ValueError("received shard differs from the selected source records")


def run_remote_transfer(config, node, jobs, writer, lane=0):
    source = Path(__file__).with_name("catalog_snapshot.py").read_text()
    if any(job.get("kind") == "raw_group" for job in jobs):
        raw_source = Path(__file__).with_name("catalog_raw.py").read_text()
        # The source module is executed in a namespace registered for its ordinary
        # import, keeping the exact same tested adapter on the source machines.
        source += "\nimport types,sys\n_adapter=types.ModuleType('catalog_snapshot')\n_adapter.__dict__.update(globals())\nsys.modules['catalog_snapshot']=_adapter\n"
        source += "\nexec(" + repr(raw_source) + ")\n"
    code = source + "\nimport sys\nJOBS=" + repr(jobs) + "\nfor job in JOBS:\n if job.get('kind')=='raw_group':\n  receipt,data=transform_raw(job['items'],job['id'])\n else:\n  receipt,data=transform_bundle(job['path'],job['sha256'],job['id'],set(job['selected']) if job.get('selected') is not None else None)\n emit_frame(sys.stdout.buffer,receipt,data)\n"
    error_path = writer.directory / f"transfer-vm{node:02d}-lane{lane}.stderr"
    with error_path.open("wb") as errors:
        process = subprocess.Popen(ssh_command(config, node), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=errors)
        def send():
            try:
                process.stdin.write(code.encode())
                process.stdin.close()
            except BrokenPipeError:
                pass
        thread = threading.Thread(target=send, daemon=True)
        thread.start()
        finished = 0
        try:
            while (frame := read_frame(process.stdout)) is not None:
                if finished >= len(jobs):
                    raise ValueError("received an unrequested shard")
                verify_receipt(jobs[finished], frame[0])
                writer.accept(*frame)
                finished += 1
                if finished % 50 == 0:
                    print(json.dumps({"node": node, "lane": lane, "shards": finished, "total": len(jobs)}), flush=True)
            if process.wait() != 0 or finished != len(jobs):
                raise ValueError(f"VM{node} transfer incomplete; inspect {error_path.name}")
        except BaseException:
            process.terminate()
            process.wait(timeout=30)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--snapshot-id", required=True)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--ssh-config", type=Path)
    parser.add_argument("--expected-assemblies", type=int)
    parser.add_argument("--workers-per-node", type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    writer = SnapshotWriter(args.output, args.snapshot_id)
    if args.bundle:
        if not args.sha256:
            parser.error("--bundle requires --sha256")
        prior = writer.db.execute("SELECT value FROM export_metadata WHERE key='local_bundle_sha256'").fetchone()
        if prior is not None and prior[0] != args.sha256:
            raise ValueError("resume requires the same local source bundle")
        writer.db.execute("INSERT OR IGNORE INTO export_metadata VALUES ('local_bundle_sha256',?)", (args.sha256,))
        writer.db.commit()
        if not writer.completed():
            writer.accept(*transform_bundle(args.bundle, args.sha256, 1))
    elif args.plan and args.ssh_config:
        opener = gzip.open if args.plan.suffix == ".gz" else open
        with opener(args.plan, "rt") as handle:
            plan = json.load(handle)
        if plan.get("schema") != "crisprloci.transfer-plan.v1" or plan.get("profile_sha256") != PROFILE:
            raise ValueError("unsupported transfer plan or scientific profile")
        with args.plan.open("rb") as handle:
            plan_digest = hashlib.file_digest(handle, "sha256").hexdigest()
        prior = writer.db.execute("SELECT value FROM export_metadata WHERE key='transfer_plan_sha256'").fetchone()
        if prior is not None and prior[0] != plan_digest:
            raise ValueError("resume requires the exact same transfer plan")
        writer.db.execute("INSERT OR IGNORE INTO export_metadata VALUES ('transfer_plan_sha256',?)", (plan_digest,))
        writer.db.commit()
        done = writer.completed()
        writer.verify_files()
        expected_ids = {job["id"] for jobs in plan["nodes"].values() for job in jobs}
        if len(expected_ids) != sum(len(jobs) for jobs in plan["nodes"].values()) or not done <= expected_ids:
            raise ValueError("snapshot shard membership does not match the transfer plan")
        with ThreadPoolExecutor(max_workers=8 * args.workers_per_node) as pool:
            futures = []
            for node, jobs in plan["nodes"].items():
                pending = [job for job in jobs if job["id"] not in done]
                for lane in range(args.workers_per_node):
                    partition = pending[lane::args.workers_per_node]
                    if partition:
                        futures.append(pool.submit(run_remote_transfer, args.ssh_config, int(node), partition, writer, lane))
            for future in as_completed(futures):
                future.result()
        if writer.completed() != expected_ids:
            raise ValueError("snapshot lacks requested source shards")
    else:
        parser.error("provide --bundle/--sha256 or --plan/--ssh-config")
    print(json.dumps(writer.publish(args.expected_assemblies), indent=2))


if __name__ == "__main__":
    main()
