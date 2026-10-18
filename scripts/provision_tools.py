#!/usr/bin/env python3
"""Build an isolated, commit-pinned scientific toolset without activating it."""

from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from runtime.common import (
    ROOT,
    RUNTIME,
    atomic_json,
    digest,
    environment_inventory,
    export_repository,
    load_lock,
    lock_id,
    run,
    source_inventory,
)

from runtime.context_tools import (
    install_defense_models,
    install_pvogs,
    native_identity,
    verify_context_native_receipt,
)
from runtime.progress import emit

CONTEXT_TOOLS = {
    "clinker-v2",
    "Phigaro-v2",
    "ISEScan-v2",
    "DefenseFinder-v2",
    "RNAalifold-v2",
}
CONTEXT_CONDA = {
    "Phigaro-v2": ["python=3.12", "pip", "prodigal=2.6.3", "hmmer=3.4"],
    "DefenseFinder-v2": ["python=3.12", "pip", "hmmer=3.4"],
    "ISEScan-v2": [
        "python=3.11",
        "pip",
        "numpy=1.26.4",
        "scipy=1.13.1",
        "biopython=1.85",
        "fastcluster=1.2.6",
        "fraggenescan=1.32",
        "hmmer=3.4",
        "blast=2.16.0",
    ],
    "RNAalifold-v2": ["python=3.12", "pip", "mafft=7.526"],
}
CONTEXT_REQUIREMENTS = {
    "clinker-v2": "requirements-runtime.lock",
    "Phigaro-v2": "requirements.lock",
    "DefenseFinder-v2": "requirements-context.lock",
    "RNAalifold-v2": "requirements.lock",
}


def context_command(name, package, prefix, command, env):
    """Operator-owned native references are baked into immutable launchers."""
    native = {}
    if name in {"ISEScan-v2", "DefenseFinder-v2", "Phigaro-v2"}:
        native["hmmsearch"] = native_identity(
            prefix / "bin/hmmsearch", "3.4", ["-h"], "HMMER 3.4 ", env
        )
    if name == "ISEScan-v2":
        native["phmmer"] = native_identity(
            prefix / "bin/phmmer", "3.4", ["-h"], "HMMER 3.4 ", env
        )
        metadata = list((prefix / "conda-meta").glob("fraggenescan-1.32-*.json"))
        if (
            len(metadata) != 1
            or json.loads(metadata[0].read_text()).get("version") != "1.32"
        ):
            raise ValueError("ISEScan requires a pinned FragGeneScan 1.32 package")
        executable = (prefix / "bin/FragGeneScan").resolve(strict=True)
        executable.relative_to(prefix.resolve(strict=True))
        native["FragGeneScan"] = {
            "version": "1.32",
            "path": str(executable),
            "sha256": digest(executable),
        }
        for binary in ("blastn", "blastp", "makeblastdb"):
            native[binary] = native_identity(
                prefix / "bin" / binary, "2.16.0", ["-version"], "2.16.0+", env
            )
    if name == "DefenseFinder-v2":
        models, native["defense-finder-models"] = install_defense_models(
            package, prefix
        )
        command += ["--models-dir", models, "--hmmer", prefix / "bin/hmmsearch"]
    if name == "Phigaro-v2":
        models, native["pVOGs"] = install_pvogs(prefix)
        native["prodigal"] = native_identity(
            prefix / "bin/prodigal", "2.6.3", ["-v"], "Prodigal V2.6.3", env
        )
        command += [
            "predict",
            "--models-hmm",
            models,
            "--models-sha256",
            native["pVOGs"]["sha256"],
            "--prodigal",
            prefix / "bin/prodigal",
            "--hmmsearch",
            prefix / "bin/hmmsearch",
        ]
    if name == "RNAalifold-v2":
        native.update(verify_repeat_runtime(prefix, env))
        native["mafft"] = native_identity(
            prefix / "bin/mafft", "7.526", ["--version"], "v7.526", env
        )
        command += ["--mafft", prefix / "bin/mafft"]
    verify_context_native_receipt({"name": name, "native_tools": native}, prefix)
    return command, native


def wrapper(path, command, environment_bin, extra_bins=()):
    # Paths are local operator-owned values, never request data.
    prefix = shlex.quote(":".join(str(item) for item in (environment_bin, *extra_bins)))
    path.write_text(
        "#!/bin/sh\nset -eu\nexport PYTHONNOUSERSITE=1\nunset PYTHONPATH\nunset MAFFT_BINARIES\n"
        f'export PATH={prefix}:"$PATH"\n'
        "export QT_QPA_PLATFORM=offscreen\n"
        + "exec "
        + " ".join(shlex.quote(str(arg)) for arg in command)
        + ' "$@"\n'
    )
    path.chmod(0o755)


def verify_repeat_runtime(prefix, env):
    """Bind the pinned Python distribution to its imported native extension."""
    probe = (
        "import importlib.metadata, json, RNA, RNA._RNA; "
        "print(json.dumps({'version': importlib.metadata.version('ViennaRNA'), "
        "'binding_version': RNA.__version__, 'path': RNA._RNA.__file__}))"
    )
    identity = json.loads(
        subprocess.check_output(
            [prefix / "bin/python", "-I", "-c", probe], text=True, env=env, timeout=30
        )
    )
    if (
        not isinstance(identity, dict)
        or set(identity) != {"version", "binding_version", "path"}
        or identity["version"] != "2.7.0"
        or identity["binding_version"] != "2.7.0"
    ):
        raise ValueError("CRISPRrepeat requires the ViennaRNA 2.7.0 Python binding")
    extension = Path(identity["path"]).resolve(strict=True)
    extension.relative_to(prefix.resolve(strict=True))
    return {
        "ViennaRNA": {**identity, "path": str(extension), "sha256": digest(extension)}
    }


def build_spacepharer_runtime(source, prefix, build_directory, version, env):
    """Compile only the pinned first-party/vendored source into its isolated env."""
    cmake = prefix / "bin/cmake"
    run(
        [
            cmake,
            "-S",
            source,
            "-B",
            build_directory,
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DHAVE_SSE4_1=1",
            "-DHAVE_AVX2=0",
            f"-DVERSION_OVERRIDE={version}",
            f"-DCMAKE_MAKE_PROGRAM={prefix / 'bin/ninja'}",
        ],
        env=env,
    )
    run([cmake, "--build", build_directory, "--parallel", "2"], env=env)
    run(
        [prefix / "bin/ctest", "--test-dir", build_directory, "--output-on-failure"],
        env=env,
    )
    executable = prefix / "bin/spacepharer-native"
    shutil.copy2(build_directory / "src/spacepharer-native", executable)
    actual = subprocess.check_output(
        [executable, "version"], text=True, env=env, timeout=30
    ).strip()
    if actual != version:
        raise ValueError("SpacePHARER native version differs from its source pin")
    return {
        "spacepharer-native": {
            "version": version,
            "path": str(executable),
            "sha256": digest(executable),
        }
    }


def build(lock_path, resume=False):
    lock = load_lock(lock_path)
    identity = lock_id(lock)
    target = RUNTIME / "toolsets" / identity
    if (target / "toolset.json").exists():
        from manage_release import verify_toolset

        verify_toolset(target)
        print("Toolset already built:", target)
        return target
    if target.exists() and not resume:
        raise ValueError(
            "Incomplete toolset exists; inspect its logs then use --resume"
        )
    target.mkdir(parents=True, exist_ok=True)
    atomic_json(target / "tools.lock.json", lock)
    conda = shutil.which("conda") or shutil.which("micromamba") or shutil.which("mamba")
    if not conda:
        raise ValueError(
            "Conda, Micromamba or Mamba is required for scientific native dependencies"
        )
    env = dict(os.environ, PIP_DISABLE_PIP_VERSION_CHECK="1", PYTHONNOUSERSITE="1")
    env.pop("PYTHONPATH", None)
    source = target / "sources"
    envs = target / "envs"
    bin_dir = target / "bin"
    for directory in (source, envs, bin_dir):
        directory.mkdir(exist_ok=True)
    for tool in lock["tools"]:
        destination = source / tool["name"]
        if not destination.exists():
            export_repository(tool, destination)
    # Reject incompatible combined-workflow pins before provisioning anything.
    upstream = json.loads((source / "CRISPRtracrRNA-v3/upstream-lock.json").read_text())
    for tool in lock["tools"]:
        if (
            tool["name"] in upstream
            and upstream[tool["name"]]["commit"] != tool["commit"]
        ):
            raise ValueError(
                f"{tool['name']} differs from CRISPRtracrRNA's validated dependency lock; update the compatible set together"
            )
    identities = {}
    for position, tool in enumerate(lock["tools"], 1):
        name = tool["name"]
        emit("Scientific environments: " + name, position - 1, len(lock["tools"]), "tools")
        package = source / name
        marker = target / f".{name}.installed.json"
        if marker.exists() and resume:
            record = json.loads(marker.read_text())
            receipt = target / f".{name}.environment.json"
            expected_prefix = (
                package / ".envs/crispridentify-core"
                if name == "CRISPRidentify-v2"
                else envs / name
            )
            launcher = bin_dir / name
            if (
                any(record.get(key) != value for key, value in tool.items())
                or record.get("source") != str(package)
                or record.get("environment") != str(expected_prefix)
                or record.get("command") != [str(launcher)]
                or record.get("source_files") != source_inventory(package)
                or record.get("launcher_sha256") != digest(launcher)
                or not receipt.is_file()
                or record.get("environment_sha256") != digest(receipt)
                or json.loads(receipt.read_text())
                != environment_inventory(expected_prefix)
            ):
                raise ValueError(
                    f"{name} source or installed environment receipt changed before resume"
                )
            if name == "CRISPRrepeat" and record.get(
                "native_tools"
            ) != verify_repeat_runtime(Path(record["environment"]), env):
                raise ValueError(
                    "CRISPRrepeat fold-engine receipt changed before resume"
                )
            if name in CONTEXT_TOOLS:
                verify_context_native_receipt(record, expected_prefix)
            identities[name] = record
            emit("Scientific environments: " + name, position, len(lock["tools"]), "tools")
            continue
        if name == "CRISPRidentify-v2":
            run(["bash", package / "install.sh"], cwd=package, env=env)
            prefix = package / ".envs/crispridentify-core"
            executable = package / "crispridentify"
            command = [executable]
        else:
            prefix = envs / name
            if not (prefix / "bin/python").exists():
                if name in {
                    "CasAndra",
                    "CRISPRleader-v2",
                    "CRISPRrepeat",
                    "CRISPRspacer",
                    "CRISPRmap-v2",
                    "RepeatTyper",
                    "SpacePHARER",
                    "CCTK",
                    "clinker-v2",
                }:
                    run([sys.executable, "-m", "venv", prefix], env=env)
                else:
                    dependencies = CONTEXT_CONDA.get(name, ["python=3.12", "pip"])
                    if name in CONTEXT_CONDA:
                        pass
                    elif name == "CRISPRtracrRNA-v3":
                        dependencies += [
                            "infernal=1.1.4",
                            "fasta3=36.3.8",
                            "blast>=2.12",
                        ]
                    else:
                        dependencies += ["mafft=7.526", "graphviz"]
                    run(
                        [
                            conda,
                            "create",
                            "--yes",
                            "--prefix",
                            prefix,
                            "--override-channels",
                            "-c",
                            "conda-forge",
                            "-c",
                            "bioconda",
                            *dependencies,
                        ],
                        env=env,
                    )
            python = prefix / "bin/python"
            install = [python, "-m", "pip", "install"]
            constraints = package / (
                "requirements-tested.txt"
                if name == "CasAndra"
                else "constraints-tested.txt"
            )
            if constraints.exists():
                install += ["-c", constraints]
            if name in {"RepeatTyper", "SpacePHARER", "CCTK"}:
                install += ["-r", package / "requirements.lock"]
            if name in CONTEXT_REQUIREMENTS:
                install += ["-r", package / CONTEXT_REQUIREMENTS[name]]
            if name == "CRISPRrepeat":
                install += ["ViennaRNA==2.7.0"]
            run([*install, package], cwd=package, env=env)
            run([python, "-m", "pip", "check"], env=env)
            command = [prefix / "bin" / tool["entrypoint"]]
        native_tools = (
            verify_repeat_runtime(prefix, env) if name == "CRISPRrepeat" else None
        )
        if name == "SpacePHARER":
            native_tools = build_spacepharer_runtime(
                package, prefix, target / "build" / name, tool["version"], env
            )
        if name in CONTEXT_TOOLS:
            command, native_tools = context_command(name, package, prefix, command, env)
        extra_bins = []
        if name == "CRISPRtracrRNA-v3":
            # Bioconda's IntaRNA 3.2.0 needs ViennaRNA tied to Python <=3.9.
            # Keep its native executable in a separate prefix; the v3 Python
            # orchestration package remains on supported Python 3.12.
            native_prefix = envs / "IntaRNA"
            if not (native_prefix / "bin/IntaRNA").exists():
                run(
                    [
                        conda,
                        "create",
                        "--yes",
                        "--prefix",
                        native_prefix,
                        "--override-channels",
                        "-c",
                        "conda-forge",
                        "-c",
                        "bioconda",
                        "intarna=3.2.0",
                    ],
                    env=env,
                )
            extra_bins.append(native_prefix / "bin")
            (target / "IntaRNA.conda-explicit.txt").write_text(
                subprocess.check_output(
                    [conda, "list", "--prefix", str(native_prefix), "--explicit"],
                    text=True,
                    env=env,
                )
            )
        launcher = bin_dir / name
        wrapper(launcher, command, prefix / "bin", extra_bins)
        run([launcher, "--help"], env=env)
        pip_freeze = subprocess.check_output(
            [prefix / "bin/python", "-m", "pip", "freeze"], text=True, env=env
        )
        (target / f"{name}.requirements.txt").write_text(pip_freeze)
        if (prefix / "conda-meta").is_dir():
            explicit = subprocess.check_output(
                [conda, "list", "--prefix", str(prefix), "--explicit"],
                text=True,
                env=env,
            )
            (target / f"{name}.conda-explicit.txt").write_text(explicit)
        record = {
            **tool,
            "source": str(package),
            "environment": str(prefix),
            "command": [str(launcher)],
            "launcher_sha256": digest(launcher),
            "source_files": source_inventory(package),
        }
        receipt = target / f".{name}.environment.json"
        atomic_json(receipt, environment_inventory(prefix))
        record["environment_sha256"] = digest(receipt)
        if native_tools is not None:
            record["native_tools"] = native_tools
        atomic_json(marker, record)
        identities[name] = record
        emit("Scientific environments: " + name, position, len(lock["tools"]), "tools")
    atomic_json(
        target / "runners.json",
        {
            "blastn": [str(envs / "CRISPRtracrRNA-v3/bin/blastn")],
            **{
                key: identities[name]["command"]
                for key, name in {
                    "crispridentify": "CRISPRidentify-v2",
                    "casandra": "CasAndra",
                    "crispr_tracrrna": "CRISPRtracrRNA-v3",
                    "spacerplacer": "SpacerPlacer",
                    "crispr_evor": "evOr",
                    "crisprleader": "CRISPRleader-v2",
                    "crisprrepeat": "CRISPRrepeat",
                    "crisprspacer": "CRISPRspacer",
                    "crisprmap": "CRISPRmap-v2",
                    "repeattyper": "RepeatTyper",
                    "spacepharer": "SpacePHARER",
                    "cctk": "CCTK",
                    "clinker": "clinker-v2",
                    "phigaro": "Phigaro-v2",
                    "isescan": "ISEScan-v2",
                    "defensefinder": "DefenseFinder-v2",
                    "rnaalifold": "RNAalifold-v2",
                }.items()
                if name in identities
            },
        },
    )
    environments = {
        str(Path(item["environment"]).relative_to(target)): environment_inventory(
            item["environment"]
        )
        for item in identities.values()
    }
    environments["envs/IntaRNA"] = environment_inventory(envs / "IntaRNA")
    atomic_json(target / "environments.json", environments)
    atomic_json(
        target / "toolset.json",
        {
            "schema_version": 1,
            "toolset_id": identity,
            "tools": identities,
            "environments_sha256": digest(target / "environments.json"),
            "runners_sha256": digest(target / "runners.json"),
        },
    )
    print("Scientific toolset staged:", target)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=ROOT / "tools/tools.lock.json")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    build(args.lock, args.resume)


if __name__ == "__main__":
    main()
