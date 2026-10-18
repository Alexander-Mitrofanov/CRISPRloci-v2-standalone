#!/usr/bin/env python3
"""Generate public examples from two independent jobs in an isolated local release.

No production service or selected release is changed. Raw jobs and comparison
evidence stay in .runtime; only explicit public input/result contracts are saved.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import time
from urllib.error import URLError

from gate_release import request
from manage_release import environment, verify
from runtime.common import ROOT, RUNTIME, atomic_json, digest


def masked_fasta(text):
    rows, sequence = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if sequence:
                rows.append("".join(sequence))
            sequence = []
        else:
            sequence.append(line.strip())
    if sequence:
        rows.append("".join(sequence))
    return "".join(f">example_record_{i:02d}\n{seq}\n" for i, seq in enumerate(rows, 1))


def scientific_job(job):
    # Execution time is inherently different for independent jobs. Retain its
    # actual value in the published result, but exclude it from reproducibility.
    def stable(value):
        if isinstance(value, list):
            return [stable(row) for row in value]
        if isinstance(value, dict):
            return {key: stable(row) for key, row in value.items()
                    if key not in {"duration_seconds", "run_time", "runtime_seconds"}}
        return value
    return stable({key: job[key] for key in ("status", "mode", "input", "options", "summary")})


def prepare(release, work, document):
    env = environment(release, work / "data", verified=document)
    source = Path(document["toolset"]) / "sources"
    runners = json.loads((Path(document["toolset"]) / "runners.json").read_text())
    references = work / "reference-inputs"
    references.mkdir(exist_ok=True)

    def run(command):
        subprocess.run(command, env=env, check=True, timeout=180,
                       stdout=subprocess.DEVNULL)

    spacer = references / "spacers.fasta"
    spacer.write_text(masked_fasta((source / "CRISPRspacer/examples/ncbi-spacerome-three.fasta").read_text()))
    database = work / "spacer-database"
    if not database.exists():
        run([*runners["crisprspacer"], "database-index", "--spacers", str(spacer),
             "--name", "example-three-masked-spacers", "--source-url", "https://example.org/masked-demonstration",
             "--source-release", "example-v1", "--output", str(database)])
    env.update(CRISPR_WEB_SPACER_DATABASE=str(database),
               CRISPR_WEB_SPACER_DATABASE_SHA256=json.loads((database / "reference.json").read_text())["database_sha256"])
    viral_source = references / "viral.fasta"
    viral_source.write_text((source / "CRISPRspacer/examples/synthetic-viral-targets.fasta").read_text())
    viral = work / "viral-database"
    blast_bin = Path(runners["blastn"][0]).parent
    if not viral.exists():
        run([*runners["crisprspacer"], "viral-index", "--sources", str(viral_source),
             "--name", "example-synthetic-viruses", "--source-release", "example-v1",
             "--source-url", "https://example.org/synthetic", "--output", str(viral),
             "--makeblastdb", str(blast_bin / "makeblastdb"), "--blastdbcmd", str(blast_bin / "blastdbcmd")])
    env.update(CRISPR_WEB_VIRAL_REFERENCE=str(viral), CRISPR_WEB_VIRAL_REFERENCE_SHA256=digest(viral / "reference.json"))
    arguments = []
    for i, archive in enumerate(sorted((source / "SpacePHARER/examples").glob("*.gz")), 1):
        target = references / f"phage-{i}.fasta"
        target.write_text(masked_fasta(gzip.decompress(archive.read_bytes()).decode()).replace("example_record_01", f"example_phage_{i:02d}"))
        arguments.extend(("--genome", str(target)))
    panel = work / "association-database"
    if not panel.exists():
        run([*runners["spacepharer"], "prepare-reference", *arguments,
             "--name", "example-three-masked-phages", "--output", str(panel)])
    env.update(CRISPR_WEB_ASSOCIATION_REFERENCE=str(panel),
               CRISPR_WEB_ASSOCIATION_REFERENCE_SHA256=digest(panel / "reference.json"))
    reference = source / "CRISPRmap-v2/examples/reference.json"
    env.update(CRISPR_WEB_REPEAT_REFERENCE=str(reference), CRISPR_WEB_REPEAT_REFERENCE_SHA256=digest(reference))

    evolution = (ROOT / "frontend/public/example-input.fasta").read_text()
    single_array = evolution.split(">example_record_02")[0]
    complete = masked_fasta((source / "CRISPRtracrRNA-v3/examples/complete.fasta").read_text())
    repeat = json.loads(reference.read_text())["entries"][0]["sequence"]
    seq = spacer.read_text().splitlines()[1]
    reverse = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    lines = (source / "SpacePHARER/tests/fixtures/spacers.fasta").read_text().splitlines()
    records = {lines[i][1:]: lines[i+1] for i in range(0, len(lines), 2)}
    groups = json.loads((source / "SpacePHARER/tests/fixtures/groups.json").read_text())["groups"]
    rows = []
    for group_index, group in enumerate(groups, 1):
        source_id = f"example_source_{group_index:02d}"
        for i, key in enumerate(group["members"], 1):
            sequence = records[key]
            rows.append({"id": f"example_record_{len(rows)+1:02d}", "sequence": sequence, "origin": {
                "source_id": source_id, "source_sha256": hashlib.sha256(source_id.encode()).hexdigest(),
                "array_id": "example-array", "ordinal": i, "start": i*200+1, "end": i*200+len(sequence)}})
    association = json.dumps({"schema": {"name": "crisprloci.sequence-set", "version": "1.0.0"},
                              "kind": "spacer", "molecule": "DNA", "complete": True, "records": rows}, indent=2) + "\n"
    cases = {
        "array_compare": ((ROOT / "deploy/fixtures/verification-array-sharing.json").read_text(), {}),
        "spacer_association": (association, {"association_grouping": "per_source"}),
        "repeat_type": (masked_fasta((source / "RepeatTyper/examples/repeats.fasta").read_text()), {}),
        "repeat_map": (f">example_record_01\n{repeat}\n>example_record_02\n{'A'*24}\n>example_record_03\nNNNN\n", {}),
        "viral_search": (masked_fasta((source / "CRISPRspacer/examples/synthetic-viral-spacers.fasta").read_text()), {"viral_max_mismatches": 2}),
        "loci": (complete, {"tracr_model_type": "II", "leader_flank_length": 500}),
        "leader": (single_array, {"leader_flank_length": 500}),
        "cas": (complete, {}),
        "tracrrna": (masked_fasta((source / "CRISPRtracrRNA-v3/examples/model-II.fasta").read_text()), {"tracr_model_type": "II"}),
        "tracrrna-v": (masked_fasta((source / "CRISPRtracrRNA-v3/examples/model-V.fasta").read_text()), {"tracr_model_type": "V"}),
        "detection": (single_array, {}),
        "repeat_context": (single_array, {}),
        "repeats": (masked_fasta((ROOT / "backend/tests/fixtures/repeats/repeats.fasta").read_text()), {"molecule": "DNA"}),
        "repeats-rna": (masked_fasta((ROOT / "backend/tests/fixtures/repeats/rna.fasta").read_text()), {"molecule": "RNA"}),
        "protospacer": (f">example_record_01\nNN{seq}NN\n>example_record_02\n{reverse}\n", {}),
        "reconstruction": (evolution, {}),
        "orientation": (evolution, {}),
    }
    return env, cases


def generate(release, work, modes=None):
    release = release.resolve(strict=True)
    document = verify(release)
    work.mkdir(parents=True, exist_ok=True)
    env, cases = prepare(release, work, document)
    env.update(CRISPR_WEB_WORKER_CPU="2", CRISPR_WEB_MAX_RECORDS="100",
               CRISPR_WEB_RATE_LIMIT_SUBMISSIONS="100", CRISPR_WEB_STAGE_TIMEOUT_SECONDS="1800",
               CRISPR_WEB_JOB_TIMEOUT_SECONDS="5400", CRISPR_WEB_HSTS="false")
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    origin = f"http://127.0.0.1:{port}"
    python = str(release / "environment/bin/python")
    children, streams = [], []
    try:
        for name, command in [
            ("api", [python, "-m", "uvicorn", "crispr_web.api:create_app", "--factory",
                     "--host", "127.0.0.1", "--port", str(port), "--no-access-log", "--no-proxy-headers"]),
            ("worker", [python, "-m", "crispr_web.worker"]),
        ]:
            stream = (work / f"{name}.log").open("w")
            streams.append(stream)
            children.append(subprocess.Popen(command, cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT))
        deadline = time.monotonic() + 60
        while True:
            try:
                if request(origin, "/api/v1/health/fresh")["status"] == "ok":
                    break
            except (URLError, RuntimeError):
                if time.monotonic() > deadline or any(child.poll() is not None for child in children):
                    raise
            time.sleep(0.5)
        for key, (sequence, extra) in cases.items():
            if modes and key not in modes:
                continue
            mode = key.split("-")[0]
            filename = "input.json" if sequence.lstrip().startswith("{") else "input.fasta"
            submission = {"sequence": sequence, "filename": filename, "mode": mode,
                          "category_policy": "bona_fide_possible", "spacer_distance": 1,
                          "bias_corrections": True, **extra}
            jobs = []
            for run in range(2):
                saved = work / f"{key}-run-{run+1}.json"
                saved_request = work / f"{key}-request-{run+1}.json"
                if (saved.exists() and saved_request.exists() and json.loads(saved_request.read_text()) == submission
                        and json.loads(saved.read_text())["summary"]["provenance"]["release_id"] == document["release_id"]):
                    jobs.append(json.loads(saved.read_text()))
                    continue
                created = request(origin, "/api/v1/jobs", body=submission, expected=202)
                route = "/api/v1/jobs/" + created["job"]["job_id"]
                deadline = time.monotonic() + 1800
                while True:
                    job = request(origin, route, token=created["access_token"])
                    if job["status"] not in {"queued", "running"}:
                        break
                    if time.monotonic() > deadline:
                        raise RuntimeError(f"Example timed out: {key}; inspect {work}")
                    time.sleep(1)
                if job["status"] != "completed":
                    raise RuntimeError(f"Example failed: {key}: {job.get('error')}; inspect {work}")
                atomic_json(saved, job)
                atomic_json(saved_request, submission)
                jobs.append(job)
                print(f"Completed {key}: independent run {run+1}/2", flush=True)
            if jobs[0]["job_id"] == jobs[1]["job_id"] or scientific_job(jobs[0]) != scientific_job(jobs[1]):
                atomic_json(work / f"{key}-stable-1.json", scientific_job(jobs[0]))
                atomic_json(work / f"{key}-stable-2.json", scientific_job(jobs[1]))
                raise RuntimeError(f"Independent scientific results differ: {key}; inspect {work}")
            job = {field: copy.deepcopy(jobs[0][field]) for field in ("status", "mode", "input", "options", "summary")}
            job["artifacts"] = []
            options = {field: value for field, value in submission.items() if field not in {"sequence", "filename", "mode"}}
            snapshot = {"schema": {"name": "crisprloci.web-example", "version": "1.0.0"},
                        "example": {"input": {"filename": filename, "file_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                              "normalized_sha256": job["input"]["sha256"],
                                              "record_count": job["input"]["record_count"], "base_count": job["input"]["base_count"]},
                                    "submission_options": options,
                                    "source_note": "Pinned tool example with masked identifiers. Reference searches use bundled demonstration references, not the full production databases.",
                                    "verification": {"independent_runs": 2, "scientific_results_identical": True,
                                                     "release_id": document["release_id"], "variable_fields": ["execution durations"]}},
                        "job": job}
            target = ROOT / "frontend/public/examples" / key
            target.mkdir(parents=True, exist_ok=True)
            (target / filename).write_text(sequence)
            atomic_json(target / "result.json", snapshot)
            print(f"Verified and saved {key}", flush=True)
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
        for child in children:
            try:
                child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
        for stream in streams:
            stream.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, default=RUNTIME / "current")
    parser.add_argument("--work", type=Path, default=RUNTIME / "example-generation")
    parser.add_argument("--modes", nargs="+", help="Regenerate only these example keys")
    args = parser.parse_args()
    generate(args.release, args.work.resolve(), args.modes)
