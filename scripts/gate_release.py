#!/usr/bin/env python3
"""Exercise an immutable release through HTTP and real scientific subprocesses."""

from __future__ import annotations

import argparse
import io
import gzip
import json
import os
from pathlib import Path
import socket
import subprocess
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import zipfile

from manage_release import environment, verify
from runtime.common import ROOT, RUNTIME, atomic_json, digest
from runtime.context_gate import fixtures as context_fixtures, check_summary as check_context_summary, check_export as check_context_export, rejected_submissions


def request(origin, route, *, body=None, token=None, expected=200):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = "Bearer " + token
    call = Request(
        origin + route,
        data=None if body is None else json.dumps(body).encode(),
        headers=headers,
    )
    try:
        response = urlopen(call, timeout=30)
    except HTTPError as error:
        response = error
    with response:
        payload = response.read()
        if response.status != expected:
            raise RuntimeError(
                f"{route}: expected {expected}, got {response.status}: {payload[:500]!r}"
            )
        return (
            json.loads(payload)
            if "json" in response.headers.get("content-type", "")
            else payload
        )


def gate(release):
    os.umask(0o077)
    release = release.resolve(strict=True)
    document = verify(release)
    run_root = RUNTIME / "gates" / f"{release.name[:16]}-{time.time_ns()}"
    run_root.mkdir(parents=True)
    env = environment(release, run_root / "data", verified=document)
    spacer_fixtures = {}
    runners = json.loads((Path(document["toolset"]) / "runners.json").read_text())
    if "crisprspacer" in runners:
        # Release gates use exactly three bundled records; never download a database.
        sample = Path(document["toolset"]) / "sources/CRISPRspacer/examples/ncbi-spacerome-three.fasta"
        database = run_root / "spacer-sample-index"
        subprocess.run([
            *runners["crisprspacer"], "database-index", "--spacers", str(sample),
            "--name", "release-gate-three-record-sample", "--source-url",
            "https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz",
            "--source-release", "2017-04-25-three-record-sample", "--output", str(database),
        ], check=True, timeout=60, env=env)
        reference = json.loads((database / "reference.json").read_text())
        if reference["record_count"] != 3:
            raise RuntimeError("Spacer release gate requires exactly three sample records")
        env.update(CRISPR_WEB_SPACER_DATABASE=str(database),
                   CRISPR_WEB_SPACER_DATABASE_SHA256=reference["database_sha256"])
        seq = sample.read_text().splitlines()[1]
        reverse = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for name, text in {
            "spacer-positive": ">forward\nNN" + seq + "NN\n>reverse\n" + reverse + "\n",
            "spacer-negative": ">negative\n" + "A" * 100 + "\n",
        }.items():
            path = run_root / (name + ".fasta")
            path.write_text(text)
            spacer_fixtures[name] = path
    if "blastn" in runners and "crisprspacer" in runners:
        sample_root = Path(document["toolset"]) / "sources/CRISPRspacer/examples"
        viral = run_root / "viral-synthetic-reference"
        blast_bin = Path(runners["blastn"][0]).parent
        subprocess.run([
            *runners["crisprspacer"], "viral-index", "--sources", str(sample_root / "synthetic-viral-targets.fasta"),
            "--name", "synthetic-release-gate", "--source-release", "synthetic-3",
            "--source-url", "https://example.org/synthetic", "--output", str(viral),
            "--makeblastdb", str(blast_bin / "makeblastdb"), "--blastdbcmd", str(blast_bin / "blastdbcmd"),
        ], check=True, timeout=60, env=env)
        env.update(CRISPR_WEB_VIRAL_REFERENCE=str(viral), CRISPR_WEB_VIRAL_REFERENCE_SHA256=digest(viral / "reference.json"))
        spacer_fixtures["viral-positive"] = sample_root / "synthetic-viral-spacers.fasta"
        negative = run_root / "viral-negative.fasta"
        negative.write_text(">short\nACGT\n>ambiguous\n" + "N" * 25 + "\n")
        spacer_fixtures["viral-negative"] = negative
    if "spacepharer" in runners:
        sample = Path(document["toolset"]) / "sources/SpacePHARER"
        genomes = run_root / "association-genomes"
        genomes.mkdir()
        arguments = []
        for archive in sorted((sample / "examples").glob("*.gz")):
            target = genomes / archive.stem
            target.write_bytes(gzip.decompress(archive.read_bytes()))
            arguments.extend(("--genome", str(target)))
        panel = run_root / "association-panel"
        subprocess.run([*runners["spacepharer"], "prepare-reference", *arguments,
                        "--name", "release-gate-three-phages", "--output", str(panel)],
                       check=True, timeout=120, env=env)
        env.update(CRISPR_WEB_ASSOCIATION_REFERENCE=str(panel),
                   CRISPR_WEB_ASSOCIATION_REFERENCE_SHA256=digest(panel / "reference.json"))
        lines = (sample / "tests/fixtures/spacers.fasta").read_text().splitlines()
        records = {lines[i][1:]: lines[i+1] for i in range(0, len(lines), 2)}
        groups = json.loads((sample / "tests/fixtures/groups.json").read_text())["groups"]
        rows = []
        import hashlib
        for group in groups:
            for i, key in enumerate(group["members"], 1):
                seq = records[key]
                rows.append({"id": f"s{len(rows)+1}", "sequence": seq, "origin": {
                    "source_id": group["id"], "source_sha256": hashlib.sha256(group["id"].encode()).hexdigest(),
                    "array_id": "frozen-upstream", "ordinal": i, "start": i*200+1, "end": i*200+len(seq)}})
        packet = {"schema": {"name": "crisprloci.sequence-set", "version": "1.0.0"},
                  "kind": "spacer", "molecule": "DNA", "complete": True, "records": rows}
        positive = run_root / "association-positive.json"
        positive.write_text(json.dumps(packet))
        negative = run_root / "association-unsupported.fasta"
        negative.write_text(">unsupported\nNNNN\n")
        spacer_fixtures.update({"association-positive": positive, "association-unsupported": negative})
    if "repeattyper" in runners:
        spacer_fixtures["type-controls"] = Path(document["toolset"]) / "sources/RepeatTyper/examples/repeats.fasta"
    if "crisprmap" in runners:
        reference = Path(document["toolset"]) / "sources/CRISPRmap-v2/examples/reference.json"
        env.update(CRISPR_WEB_REPEAT_REFERENCE=str(reference), CRISPR_WEB_REPEAT_REFERENCE_SHA256=digest(reference))
        repeat = json.loads(reference.read_text())["entries"][0]["sequence"]
        query = run_root / "map-controls.fasta"
        query.write_text(f">exact\n{repeat}\n>no_hit\n{'A' * 24}\n>unsupported\nNNNN\n")
        spacer_fixtures["map-controls"] = query
    additional_fixtures, additional_cases = context_fixtures(Path(document["toolset"]), run_root, runners)
    spacer_fixtures.update(additional_fixtures)
    env.update(
        CRISPR_WEB_WORKER_CPU="2",
        CRISPR_WEB_MAX_RECORDS="100",
        CRISPR_WEB_RATE_LIMIT_SUBMISSIONS="50",
        CRISPR_WEB_STAGE_TIMEOUT_SECONDS="1800",
        CRISPR_WEB_JOB_TIMEOUT_SECONDS="5400",
        CRISPR_WEB_HSTS="false",
    )
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    origin = f"http://127.0.0.1:{port}"
    python = release / "environment/bin/python"
    children, streams, observations = [], [], []
    try:
        for name, command in [
            (
                "api",
                [
                    python,
                    "-m",
                    "uvicorn",
                    "crispr_web.api:create_app",
                    "--factory",
                    "--host",
                    "127.0.0.1",
                    "--port",
                    str(port),
                    "--no-access-log",
                    "--no-proxy-headers",
                ],
            ),
            ("worker", [python, "-m", "crispr_web.worker"]),
        ]:
            stream = (run_root / f"{name}.log").open("w")
            streams.append(stream)
            children.append(
                subprocess.Popen(
                    command,
                    cwd=run_root,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                )
            )
        deadline = time.monotonic() + 45
        while True:
            if any(process.poll() is not None for process in children):
                raise RuntimeError(f"Release process stopped; inspect {run_root}")
            try:
                health = request(origin, "/api/v1/health/fresh")
                if health["status"] == "ok":
                    break
            except (URLError, RuntimeError):
                if time.monotonic() >= deadline:
                    raise
            time.sleep(0.5)
        version = request(origin, "/api/v1/version")
        if (
            version["release_id"] != release.name
            or len(version["scientific_artifacts"])
            != len(
                json.loads((Path(document["toolset"]) / "toolset.json").read_text())[
                    "tools"
                ]
            )
            + 1
        ):
            raise RuntimeError("API release identity mismatch")
        manifest = json.loads((release / "scientific-manifest.json").read_text())
        for key, identity in manifest["scientific_artifacts"].items():
            observed = version["scientific_artifacts"].get(key, {})
            if any(observed.get(field) != identity.get(field)
                   for field in ("version", "commit", "repository", "sha256")):
                raise RuntimeError("API scientific source identity mismatch: " + key)
        rejected = rejected_submissions(runners)
        for body in rejected:
            request(origin, "/api/v1/jobs", body=body, expected=422)
        source = Path(document["toolset"]) / "sources/CRISPRtracrRNA-v3/examples"
        fixtures = {
            **spacer_fixtures,
            "no-arrays": ROOT / "deploy/fixtures/verification-no-arrays.fasta",
            "evolution": ROOT / "deploy/fixtures/verification-related-isolates.fasta",
            "evolution-extended": ROOT / "deploy/fixtures/verification-related-isolates.fasta",
            "complete": source / "complete.fasta",
            "model-II": source / "model-II.fasta",
            "model-V": source / "model-V.fasta",
            "repeat-DNA": ROOT / "backend/tests/fixtures/repeats/repeats.fasta",
            "repeat-RNA": ROOT / "backend/tests/fixtures/repeats/rna.fasta",
            "array-singleton": ROOT / "docs/examples/module-exchange/arrays.json",
            "array-sharing": ROOT / "deploy/fixtures/verification-array-sharing.json",
            "array-sharing-rc": ROOT / "deploy/fixtures/verification-array-sharing.json",
        }
        cases = [
            ("leader", "evolution", "II"),
            ("leader", "no-arrays", "II"),
            ("detection", "evolution", "II"),
            ("reconstruction", "evolution", "II"),
            ("orientation", "evolution", "II"),
            ("reconstruction", "evolution-extended", "II"),
            ("orientation", "evolution-extended", "II"),
            ("cas", "complete", "II"),
            ("tracrrna", "model-II", "II"),
            ("tracrrna", "model-V", "V"),
            ("loci", "complete", "II"),
        ]
        if "crisprrepeat" in version["scientific_artifacts"]:
            cases += [
                ("repeats", "repeat-DNA", "II"),
                ("repeats", "repeat-RNA", "II"),
                ("repeat_context", "evolution", "II"),
                ("repeat_context", "no-arrays", "II"),
            ]
        if "crisprspacer" in version["scientific_artifacts"]:
            cases += [("protospacer", "spacer-positive", "II"), ("protospacer", "spacer-negative", "II")]
        if "blastn" in runners and "crisprspacer" in version["scientific_artifacts"]:
            cases += [("viral_search", "viral-positive", "II"), ("viral_search", "viral-negative", "II")]
        if "spacepharer" in version["scientific_artifacts"]:
            cases += [("spacer_association", "association-positive", "II"), ("spacer_association", "association-unsupported", "II")]
        if "cctk" in version["scientific_artifacts"]:
            cases += [("array_compare", "array-singleton", "II"), ("array_compare", "array-sharing", "II"), ("array_compare", "array-sharing-rc", "II")]
        if "repeattyper" in version["scientific_artifacts"]:
            cases += [("repeat_type", "type-controls", "II")]
        if "crisprmap" in version["scientific_artifacts"]:
            cases += [("repeat_map", "map-controls", "II")]
        cases += additional_cases
        for mode, fixture, model in cases:
            created = request(
                origin,
                "/api/v1/jobs",
                body={
                    "sequence": fixtures[fixture].read_text(),
                    "filename": fixtures[fixture].name,
                    "mode": mode,
                    **(
                        {"molecule": "RNA" if fixture == "repeat-RNA" or mode == "rna_family" else "DNA"}
                        if mode in {"repeats", "repeat_context", "rna_family"}
                        else {}
                    ),
                    **({"association_grouping": "per_source" if fixture == "association-positive" else "per_record"}
                       if mode == "spacer_association" else {}),
                    **({"tool_options": {"cctk": {"--identity": "reverse-complement"}}}
                       if fixture == "array-sharing-rc" else {}),
                    **({"tool_options": {
                        "spacerplacer" if mode == "reconstruction" else "crispr_evor": {"--extend_branches": 0.01}
                    }} if fixture == "evolution-extended" else {}),
                    "tracr_model_type": model,
                    "leader_flank_length": 37 if mode == "leader" else 500,
                },
                expected=202,
            )
            identity, token = created["job"]["job_id"], created["access_token"]
            route = "/api/v1/jobs/" + identity
            request(origin, route, expected=401)
            deadline = time.monotonic() + 1800
            while True:
                job = request(origin, route, token=token)
                if job["status"] not in {"queued", "running"}:
                    break
                if time.monotonic() > deadline:
                    raise RuntimeError("Scientific release gate timed out: " + mode)
                time.sleep(1)
            if job["status"] != "completed":
                raise RuntimeError(
                    f"{mode}/{model} did not complete: {job.get('error')}; inspect {run_root}"
                )
            check_context_summary(mode, fixture, job["summary"])
            if mode == "spacer_association":
                evidence = job["summary"]["spacer_association"]
                expected = (2, 6, "unavailable") if fixture == "association-positive" else (0, 0, "not_run")
                if (evidence["counts"]["candidates"], evidence["counts"]["support_hits"], evidence["calibration_status"]) != expected:
                    raise RuntimeError("SpacePHARER real grouped calibration control differs")
                if evidence["counts"]["associations"] != 0 or any(row["native_fdr"] is not None for row in evidence["associations"]):
                    raise RuntimeError("Uncalibrated candidates were misrepresented")
            if mode == "array_compare":
                comparison = job["summary"]["array_compare"]
                counts = dict(zip(
                    ("arrays", "spacer_occurrences", "deletions", "edges", "components", "isolates"),
                    {"array-singleton": (1, 2, 0, 0, 1, 1), "array-sharing": (4, 6, 2, 1, 3, 2),
                     "array-sharing-rc": (4, 6, 2, 2, 2, 1)}[fixture],
                ))
                expected_policy = "reverse-complement" if fixture == "array-sharing-rc" else "exact"
                if comparison["counts"] != counts or comparison["options"]["identity"] != expected_policy:
                    raise RuntimeError("CCTK observed-sharing control differs")
                if fixture != "array-singleton":
                    arrays = comparison["arrays"]
                    edge = comparison["edges"][0]
                    if (edge["nshared"] != 1 or edge["jaccard"] != 1 / 3
                            or arrays[0]["strand"] is not None or arrays[1]["strand"] != "-"
                            or arrays[0]["array_id"] != arrays[1]["array_id"]
                            or arrays[0]["spacers"][2]["symbol"] is not None
                            or arrays[3]["evidence_status"] != "no_observed_spacers"):
                        raise RuntimeError("CCTK occurrence, gap or unknown-strand semantics changed")
            if mode == "repeat_type":
                typing = job["summary"]["repeat_type"]
                if typing["counts"] != {"records": 5, "distinct_classified_sequences": 3, "predicted": 3, "uncertain": 1, "unsupported": 1}:
                    raise RuntimeError("RepeatTyper control counts differ")
                if len(typing["model"]["labels"]) != 37 or typing["results"][3]["native_score"] is not None:
                    raise RuntimeError("RepeatTyper unknown/model semantics differ")
            if mode == "repeat_map":
                mapping = job["summary"]["repeat_map"]
                if mapping["reference"]["file_sha256"] != env["CRISPR_WEB_REPEAT_REFERENCE_SHA256"]:
                    raise RuntimeError("Repeat mapping reference identity mismatch")
                if [row["status"] for row in mapping["results"]] != ["exact", "no_hit", "unsupported"]:
                    raise RuntimeError("Repeat mapping controls failed")
            if mode == "viral_search":
                viral = job["summary"]["viral_search"]
                if viral["reference"]["reference_sha256"] != env["CRISPR_WEB_VIRAL_REFERENCE_SHA256"]:
                    raise RuntimeError("Viral reference identity changed")
                count = 3 if fixture == "viral-positive" else 0
                if viral["counts"]["accepted_hsps"] != count or viral["execution_completed"] is not True:
                    raise RuntimeError("Viral accepted evidence differs from synthetic fixture")
                if not count and viral["outcome"] != "no_eligible_queries":
                    raise RuntimeError("Ineligible input was presented as a negative search")
            if mode == "protospacer":
                spacer = job["summary"]["protospacer"]
                if spacer["reference"]["database_sha256"] != env["CRISPR_WEB_SPACER_DATABASE_SHA256"]:
                    raise RuntimeError("Spacer database identity changed")
                expected_count = 2 if fixture == "spacer-positive" else 0
                if spacer["counts"]["matches"] != expected_count:
                    raise RuntimeError("Spacer match count differs from release fixture")
                if expected_count and {(r["target_id"], r["target_start"], r["relative_strand"]) for r in spacer["matches"]} != {("forward", 3, "+"), ("reverse", 1, "-")}:
                    raise RuntimeError("Spacer target coordinate/orientation mismatch")
            if mode in {"repeats", "repeat_context"}:
                repeat = job["summary"]["repeats"]
                if repeat["status"] != "completed":
                    raise RuntimeError("Repeat evidence did not complete")
                if fixture == "no-arrays":
                    if repeat["counts"]["instances"] != 0:
                        raise RuntimeError(
                            "Negative repeat context fixture invented instances"
                        )
                elif repeat["counts"]["instances"] < 1:
                    raise RuntimeError("Positive repeat fixture lost instances")
                if repeat["molecule"] != ("RNA" if fixture == "repeat-RNA" else "DNA"):
                    raise RuntimeError("Repeat molecular identity changed")
            if (
                mode in {"cas", "loci"}
                and job["summary"]["casandra"]["cas_gene_count"] < 1
            ):
                raise RuntimeError("Positive Cas fixture lost its calls")
            if (
                mode in {"tracrrna", "loci"}
                and job["summary"]["tracrrna"]["prediction_count"] < 1
            ):
                raise RuntimeError("Positive tracrRNA fixture lost its candidates")
            if mode in {"leader", "loci"}:
                leader = job["summary"]["crisprleader"]
                if (
                    leader["status"] != "completed"
                    or leader["prediction_performed"] is not False
                    or leader["context_count"] != 2 * leader["array_count"]
                ):
                    raise RuntimeError("Leader context contract mismatch")
                if fixture == "evolution" and leader["array_count"] < 1:
                    raise RuntimeError("Positive leader fixture lost its arrays")
                if fixture == "no-arrays" and leader["array_count"] != 0:
                    raise RuntimeError("Empty leader fixture has unexpected arrays")
            bundle_bytes = request(origin, route + "/result", token=token)
            with zipfile.ZipFile(io.BytesIO(bundle_bytes)) as bundle:
                if (
                    bundle.testzip() is not None
                    or "provenance.json" not in bundle.namelist()
                ):
                    raise RuntimeError("Invalid result archive")
                provenance = json.loads(bundle.read("provenance.json"))
                if provenance["software"]["backend_version"] != version["version"]:
                    raise RuntimeError("Result provenance mismatch")
                if provenance["software"]["scientific_artifacts"] != version["scientific_artifacts"]:
                    raise RuntimeError("Result scientific source identity mismatch")
                if fixture == "evolution-extended" and provenance["options"]["tool_options"] != {
                    "spacerplacer" if mode == "reconstruction" else "crispr_evor": {"--extend_branches": 0.01}
                }:
                    raise RuntimeError("Branch-extension option provenance mismatch")
                compact = json.loads(bundle.read("results/result.json"))
                if (
                    compact.get("schema") != {"name": "crisprloci.module-result", "version": "1.0.0"}
                    or compact.get("module") != mode
                    or compact.get("complete") is not True
                    or set(compact) != {"schema", "module", "status", "complete", "data"}
                ):
                    raise RuntimeError("Compact scientific result is incomplete: " + mode)
                check_context_export(mode, bundle, job["summary"])
                if mode == "spacer_association":
                    evidence = json.loads(bundle.read("results/spacer-associations.json"))
                    if evidence["counts"] != job["summary"]["spacer_association"]["counts"]:
                        raise RuntimeError("SpacePHARER download differs from validated summary")
                    if any(name.startswith(("results/association/", "results/association-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private SpacePHARER working files leaked")
                if mode == "array_compare":
                    comparison = json.loads(bundle.read("results/array-comparison.json"))
                    if (comparison["counts"] != job["summary"]["array_compare"]["counts"]
                            or len(comparison["arrays"]) != comparison["counts"]["arrays"]
                            or len(comparison["edges"]) != comparison["counts"]["edges"]):
                        raise RuntimeError("CCTK export is incomplete")
                    if json.loads(bundle.read("results/arrays.json")) != json.loads(fixtures[fixture].read_text()):
                        raise RuntimeError("CCTK array handoff changed")
                    if any(name.startswith(("results/array-compare/", "results/array-compare-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private CCTK working files leaked")
                if mode == "repeat_type":
                    typing = json.loads(bundle.read("results/repeat-typing.json"))
                    if typing["counts"] != job["summary"]["repeat_type"]["counts"]:
                        raise RuntimeError("RepeatTyper export mismatch")
                    if any(name.startswith(("results/repeat-type/", "results/repeat-type-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private RepeatTyper worker output leaked")
                if mode == "repeat_map":
                    mapping = json.loads(bundle.read("results/repeat-mapping.json"))
                    if mapping["counts"] != job["summary"]["repeat_map"]["counts"]:
                        raise RuntimeError("Repeat mapping export mismatch")
                    if any(name.startswith(("results/repeat-map/", "results/repeat-map-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private repeat mapping outputs leaked")
                if mode == "viral_search":
                    viral_report = json.loads(bundle.read("results/viral-matches.json"))
                    if viral_report["counts"] != job["summary"]["viral_search"]["counts"]:
                        raise RuntimeError("Viral download differs from summary")
                    if any(name.startswith(("results/viral/", "results/viral-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private native viral files leaked")
                    if "blast" in viral_report or "/srv/" in json.dumps(viral_report):
                        raise RuntimeError("Private viral runtime identity leaked")
                if mode == "protospacer":
                    spacer_report = json.loads(bundle.read("results/protospacer-analysis.json"))
                    if spacer_report["counts"] != job["summary"]["protospacer"]["counts"]:
                        raise RuntimeError("Spacer download differs from summary")
                    if any(name.startswith(("results/protospacer/", "results/protospacer-input/")) for name in bundle.namelist()):
                        raise RuntimeError("Private native spacer files leaked")
                if mode in {"repeats", "repeat_context"}:
                    repeat_report = json.loads(
                        bundle.read("results/repeat-analysis.json")
                    )
                    if repeat_report["counts"] != job["summary"]["repeats"][
                        "counts"
                    ] or repeat_report["counts"]["instances"] != len(
                        repeat_report["instances"]
                    ):
                        raise RuntimeError("Repeat download is incomplete")
                    if any(
                        name.startswith("results/repeat/") for name in bundle.namelist()
                    ):
                        raise RuntimeError("Private native repeat files leaked")
                if mode in {"leader", "loci"}:
                    contexts = json.loads(bundle.read("results/leader-contexts.json"))
                    if (
                        contexts["complete"] is not True
                        or len(contexts["contexts"]) != contexts["context_count"]
                        or contexts["prediction_performed"] is not False
                    ):
                        raise RuntimeError(
                            "Leader download is incomplete or claims prediction"
                        )
                if mode in {"cas", "loci", "tracrrna"}:
                    annotations = json.loads(bundle.read("results/annotations.json"))
                    if annotations["complete"] is not True:
                        raise RuntimeError("Annotation download is incomplete")
                    if mode in {"tracrrna", "loci"}:
                        rna = annotations["tracrrna"]
                        if len(rna["candidates"]) != rna["prediction_count"]:
                            raise RuntimeError("RNA export omitted predictions")
            snapshot = run_root / f"{mode}-{fixture}-{model}.json"
            atomic_json(snapshot, job)
            observations.append(
                {
                    "mode": mode,
                    "fixture": fixture,
                    "model_type": model,
                    "status": job["status"],
                    "job_id": identity,
                    "fixture_sha256": digest(fixtures[fixture]),
                    "result_bytes": len(bundle_bytes),
                }
            )
            print(
                f"PASS {mode}/{model}: authenticated result, validated annotations and archive",
                flush=True,
            )
        verify(release)
        result = {
            "status": "passed",
            "release_sha256": digest(release / "release.json"),
            "completed_at": time.time(),
            "gate_script_sha256": digest(Path(__file__)),
            "context_gate_sha256": digest(ROOT / "scripts/runtime/context_gate.py"),
            "context_predispatch_rejections": len(rejected),
            "run_root": str(run_root),
            "cases": observations,
        }
        atomic_json(release / "gate.json", result)
        print("Release gate passed:", release)
    finally:
        for process in children:
            if process.poll() is None:
                process.terminate()
        for process in children:
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        for stream in streams:
            stream.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path)
    gate(parser.parse_args().release)
