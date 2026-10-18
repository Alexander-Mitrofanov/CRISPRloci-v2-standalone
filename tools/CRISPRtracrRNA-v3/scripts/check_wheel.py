#!/usr/bin/env python3
"""Check one built wheel in a fresh environment without accessing a package index."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import tempfile
import venv
from pathlib import Path


def verify_resources(package: Path) -> int:
    """Verify distribution assets without invoking scientific executables."""
    if not (package / "py.typed").is_file():
        raise RuntimeError("Installed package is missing py.typed")
    data = package / "data"
    manifest = json.loads((data / "manifest.json").read_text())
    required = {"erpin/erpin", "erpin/rho-indep.epn"}
    required.update(
        f"models/{model}.cm{suffix}"
        for model in ("II", "V-1", "V-2", "V-3", "V-4")
        for suffix in ("", ".i1f", ".i1i", ".i1m", ".i1p")
    )
    if missing := required.difference(manifest):
        raise RuntimeError(f"Installed manifest omits required resources: {sorted(missing)}")
    for relative, expected in manifest.items():
        resource = data / relative
        if not resource.resolve().is_relative_to(data.resolve()):
            raise RuntimeError(f"Manifest resource is outside package data: {relative}")
        if not resource.is_file():
            raise RuntimeError(f"Installed resource is missing: {relative}")
        if hashlib.sha256(resource.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Installed resource checksum mismatch: {relative}")
    if sys.platform == "linux" and not os.access(data / "erpin/erpin", os.X_OK):
        raise RuntimeError("Installed ERPIN executable has lost executable permission")
    return len(manifest)


def verify_installation() -> None:
    import crisprtracrrna

    package = Path(crisprtracrrna.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()) or sys.prefix == sys.base_prefix:
        raise RuntimeError(f"Package did not import from the isolated environment: {package}")
    version = importlib.metadata.version("CRISPRtracrRNA")
    if crisprtracrrna.__version__ != version:
        raise RuntimeError("Installed package and distribution versions disagree")
    resources = verify_resources(package)
    print(json.dumps({"version": version, "resources": resources, "package": str(package)}))


def check_wheel(wheel: Path) -> None:
    wheel = wheel.resolve(strict=True)
    if wheel.suffix != ".whl":
        raise ValueError(f"Expected a wheel file: {wheel}")
    with tempfile.TemporaryDirectory(prefix="crisprtracrrna-wheel-") as temporary:
        work = Path(temporary)
        environment = work / "venv"
        venv.EnvBuilder(with_pip=True).create(environment)
        binaries = environment / ("Scripts" if os.name == "nt" else "bin")
        python = binaries / ("python.exe" if os.name == "nt" else "python")
        # The console launcher cannot receive -I. Remove inherited Python path
        # overrides and exclude the working directory from its import search.
        child_env = {
            key: value for key, value in os.environ.items() if not key.startswith("PYTHON")
        }
        child_env["PYTHONSAFEPATH"] = "1"
        subprocess.run(
            [
                str(python),
                "-I",
                "-m",
                "pip",
                "install",
                "--no-index",
                "--no-deps",
                "--disable-pip-version-check",
                str(wheel),
            ],
            cwd=work,
            env=child_env,
            check=True,
        )
        inspected = subprocess.run(
            [str(python), "-I", str(Path(__file__).resolve()), "--check-installed"],
            cwd=work,
            env=child_env,
            check=True,
            capture_output=True,
            text=True,
        )
        result = json.loads(inspected.stdout)
        console = binaries / ("crisprtracrrna.exe" if os.name == "nt" else "crisprtracrrna")
        for command in ([str(console)], [str(python), "-I", "-m", "crisprtracrrna"]):
            for option in ("--version", "--help"):
                completed = subprocess.run(
                    [*command, option],
                    cwd=work,
                    env=child_env,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                if not completed.stdout.strip():
                    raise RuntimeError(f"Installed command produced no output: {command} {option}")
                if option == "--version" and result["version"] not in completed.stdout:
                    raise RuntimeError(f"Installed command reported a different version: {command}")
        print(
            f"Wheel check passed: {wheel.name}, version {result['version']}, "
            f"{result['resources']} verified resources, console and module entry points."
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path, nargs="?")
    parser.add_argument("--check-installed", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.check_installed:
        verify_installation()
    elif args.wheel is not None:
        check_wheel(args.wheel)
    else:
        parser.error("provide the path to one built wheel")


if __name__ == "__main__":
    main()
