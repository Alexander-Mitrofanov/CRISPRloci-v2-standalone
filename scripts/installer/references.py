"""Build user-selected immutable references using the pinned native tools."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

from .download import download, file_digest
from .engine import atomic_json
from runtime.progress import emit

SPACER_URL = "https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz"
SPACER_SHA256 = "f92ae812231acb1b59035456379245bc824e402ef66f448e26f0932c0676d83b"


def run(arguments):
    print("+ " + " ".join(str(part) for part in arguments), flush=True)
    subprocess.run([str(part) for part in arguments], check=True)


def sample_catalog(directory):
    """Generate a synthetic teaching catalogue with truthful local provenance."""
    from catalog_snapshot import create_tables, normalize_database
    from import_catalog import SnapshotWriter

    if (directory / "catalog.sqlite3").exists():
        return
    emit("Building the synthetic catalogue", 0, 1, "shards")
    connection = sqlite3.connect(":memory:")
    create_tables(connection)
    connection.execute("INSERT INTO assemblies VALUES(1,'GCA_000000000.0','complete','synthetic',1,1,0,0)")
    connection.execute("INSERT INTO sequences(sequence_pk,assembly_pk,sequence_id,length) VALUES(1,1,'synthetic-contig',120)")
    repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA"
    spacer = "ACGTACGTGATTACAGATTACAACGTACGT"
    repeat_length, spacer_length = len(repeat), len(spacer)
    end = 2 * repeat_length + spacer_length
    connection.execute("INSERT INTO arrays(array_pk,sequence_pk,tool,call_id,start,end,repeat_consensus,repeat_count,spacer_count,category) VALUES(1,1,'synthetic-example','example-array',1,?,?,2,1,'Bona-fide')", (end, repeat))
    connection.execute("INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,length,sequence) VALUES(1,1,1,1,?,?,?)", (repeat_length, repeat_length, repeat))
    connection.execute("INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,length,sequence) VALUES(2,1,2,?,?,?,?)", (repeat_length + spacer_length + 1, end, repeat_length, repeat))
    connection.execute("INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,start,end,length,sequence) VALUES(1,1,1,'observed',?,?,?,?)", (repeat_length + 1, repeat_length + spacer_length, spacer_length, spacer))
    connection.commit()
    compressed, registry, counts, size = normalize_database(connection.serialize(), 1)
    connection.close()
    writer = SnapshotWriter(directory, "standalone-synthetic-example-v1")
    writer.accept({"id": 1, "sha256": hashlib.sha256(compressed).hexdigest(),
                   "uncompressed_bytes": size, "counts": counts, "assemblies": registry}, compressed)
    published = writer.publish(expected=1)
    summary = {"available": True, "snapshot_id": "standalone-synthetic-example-v1",
               "title": "Synthetic installation example — not the research catalogue",
               "counts": published["counts"], "coverage": {"catalog_assemblies": 1, "imported_assemblies": 1},
               "limitations": ["Synthetic sequences and annotations demonstrate catalogue searches only; they are not biological findings or the hosted CRISPRcasDb collection."],
               "provenance": {"source_project": "CRISPRloci standalone synthetic installation example",
                              "storage": "compressed-sqlite-shards-v1",
                              "shards_sha256": published["provenance"]["shards_sha256"]}}
    with sqlite3.connect(directory / "catalog.sqlite3") as registry_db:
        registry_db.execute("UPDATE export_metadata SET value=? WHERE key='summary'", (json.dumps(summary),))
    atomic_json(directory / "manifest.json", summary)
    emit("Building the synthetic catalogue", 1, 1, "shards")


def build_reference(root, work, key, choice, release):
    document = json.loads((release / "release.json").read_text())
    toolset = Path(document["toolset"])
    tools = json.loads((toolset / "toolset.json").read_text())["tools"]
    runners = json.loads((toolset / "runners.json").read_text())
    sources = toolset / "sources"
    destination = work / "references" / key
    destination.parent.mkdir(parents=True, exist_ok=True)
    mode = choice["mode"]
    worker = root / "scripts/installer_reference_worker.py"
    expected = {}
    if mode in {"local", "import"}:
        source = Path(choice["path"])
        if source.is_file():
            expected["input_sha256"] = file_digest(source)
        else:
            marker = source / ("catalog.sqlite3" if key == "catalog" else "")
            if key == "catalog":
                expected["input_sha256"] = file_digest(marker)
            else:
                expected["input_sha256"] = hashlib.sha256(json.dumps(
                    {p.name: file_digest(p) for p in sorted(source.iterdir()) if p.is_file()}, sort_keys=True).encode()).hexdigest()
    receipt = work / (key + ".reference.json")
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved.get("input", {}) != expected:
            raise ValueError("Local reference input changed since this installation started; choose a new input path or clear its installer work directory")
    spacer_python = Path(tools["CRISPRspacer"]["environment"]) / "bin/python"
    blast_bin = Path(runners["blastn"][0]).parent
    env = {}
    if key == "spacer":
        source = sources / "CRISPRspacer/examples/ncbi-spacerome-three.fasta"
        if mode == "download":
            source = root / ".installer/cache/spacerome-2017.fna.gz"
            download(SPACER_URL, source, SPACER_SHA256, emit, limit=64 * 1024**2)
        elif mode == "local":
            source = Path(choice["path"])
        if not destination.exists():
            run([spacer_python, worker, "spacer-index", "--spacers", source, "--name", "standalone-" + mode,
                 "--source-url", SPACER_URL if mode != "local" else "https://example.org/user-supplied",
                 "--source-release", "2017-04-25" if mode == "download" else mode, "--output", destination])
        run([*runners["crisprspacer"], "validate-database", destination])
        reference = json.loads((destination / "reference.json").read_text())
        env = {"CRISPR_WEB_SPACER_DATABASE": str(destination),
               "CRISPR_WEB_SPACER_DATABASE_SHA256": reference["database_sha256"]}
    elif key == "viral":
        files = [sources / "CRISPRspacer/examples/synthetic-viral-targets.fasta"]
        release_name, source_url, manifest = mode, "https://example.org/standalone-synthetic", None
        if mode == "download":
            plan = work / "viral-plan.json"
            if not plan.exists():
                run([*runners["crisprspacer"], "viral-plan", "--output", plan])
            downloaded = work / "viral-sources"
            if not downloaded.exists():
                run([spacer_python, worker, "viral-download", "--plan", plan, "--output", downloaded])
            metadata = json.loads((downloaded / "download.json").read_text())
            files = [downloaded / item["name"] for item in metadata["files"]]
            release_name, source_url = metadata["source_release"], metadata["source_url"]
            manifest = downloaded / "download.json"
        elif mode == "local":
            files, source_url = [Path(choice["path"])], "https://example.org/user-supplied"
        if not destination.exists():
            command = [spacer_python, worker, "viral-index", "--sources", *files, "--name", "standalone-" + mode,
                       "--source-release", release_name, "--source-url", source_url, "--output", destination,
                       "--makeblastdb", blast_bin / "makeblastdb", "--blastdbcmd", blast_bin / "blastdbcmd"]
            if manifest:
                command += ["--source-manifest", manifest]
            run(command)
        run([*runners["crisprspacer"], "validate-viral-reference", destination, "--blastdbcmd", blast_bin / "blastdbcmd"])
        env = {"CRISPR_WEB_VIRAL_REFERENCE": str(destination),
               "CRISPR_WEB_VIRAL_REFERENCE_SHA256": file_digest(destination / "reference.json")}
    elif key == "association":
        args = []
        if mode == "sample":
            genomes = work / "sample-phages"
            genomes.mkdir(exist_ok=True)
            for source in sorted((sources / "SpacePHARER/examples").glob("*.gz")):
                target = genomes / source.stem
                if not target.exists():
                    with gzip.open(source, "rb") as input_stream, target.open("xb") as output:
                        shutil.copyfileobj(input_stream, output, 1024 * 1024)
                args += ["--genome", target]
        elif mode == "download":
            acquisition = work / "phage-acquisition"
            run([sys.executable, root / "scripts/prepare_phage_reference.py", "acquire", "--snapshot", acquisition])
            if not (acquisition / "selected/selection.json").exists():
                run([sys.executable, root / "scripts/prepare_phage_reference.py", "select", "--snapshot", acquisition])
            args += ["--genome-list", acquisition / "selected/genomes.txt", "--source-metadata", acquisition / "selected/source.json",
                     "--scope", "broad_phage_reference"]
        else:
            files = [p for p in sorted(Path(choice["path"]).iterdir()) if p.is_file() and p.name.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"))]
            if not files:
                raise ValueError("No genome FASTA files found in the chosen directory")
            listing = work / "local-genomes.txt"
            listing.write_text("\n".join(str(p) for p in files) + "\n")
            args += ["--genome-list", listing]
        emit("Building and validating the association reference")
        if not destination.exists():
            run([*runners["spacepharer"], "prepare-reference", *args, "--name", "standalone-" + mode, "--output", destination])
        run([*runners["spacepharer"], "verify-reference", "--reference", destination])
        env = {"CRISPR_WEB_ASSOCIATION_REFERENCE": str(destination),
               "CRISPR_WEB_ASSOCIATION_REFERENCE_SHA256": file_digest(destination / "reference.json")}
    elif key == "repeat":
        reference = sources / "CRISPRmap-v2/examples/reference.json"
        if mode == "historical":
            python = Path(tools["CRISPRmap-v2"]["environment"]) / "bin/python"
            run([python, worker, "historical", "--source", sources / "CRISPRmap-v2", "--output", destination])
            reference = destination / "reference.json"
        elif mode == "local":
            reference = Path(choice["path"])
        run([*runners["crisprmap"], "inspect-reference", "--reference", reference])
        env = {"CRISPR_WEB_REPEAT_REFERENCE": str(reference), "CRISPR_WEB_REPEAT_REFERENCE_SHA256": file_digest(reference)}
    elif key == "catalog":
        from index_catalog import build_index, read_registry

        if mode == "sample":
            sample_catalog(destination)
        elif mode == "local":
            destination = Path(choice["path"])
        else:
            if not (destination / "catalog.sqlite3").exists():
                run([sys.executable, root / "scripts/import_catalog.py", "--output", destination,
                     "--snapshot-id", "standalone-import", "--bundle", choice["path"], "--sha256", choice["sha256"]])
        summary, shards, _ = read_registry(destination)
        emit("Building the catalogue search index", 0, len(shards), "shards")
        build_index(destination, workers=2)
        emit("Building the catalogue search index", len(shards), len(shards), "shards")
        env = {"CRISPR_WEB_CATALOG_DATABASE": str(destination / "catalog.sqlite3")}
        print("Catalogue: " + summary["title"], flush=True)
    else:
        raise ValueError("Unknown reference task")
    atomic_json(receipt, {"environment": env, "input": expected, "choice": choice})
    return env
