#!/usr/bin/env python3
"""Validate the small, supported CRISPRidentify runtime environment."""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys

try:
    from importlib import metadata as importlib_metadata
except ImportError:  # Python 3.7 audit/runtime compatibility.
    import importlib_metadata


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


CORE_PACKAGES = {
    "CasAndra": "0.1.0.dev0+crispridentify.3",
    "biopython": "1.86",
    "joblib": "1.5.3",
    "numpy": "1.26.4",
    "pyhmmer": "0.12.1",
    "pyrodigal": "3.7.1",
    "regex": "2024.11.6",
    "scikit-learn": "1.8.0",
    "scipy": "1.17.1",
}
CORE_EXECUTABLES = (
    "casandra",
    "mkvtree",
    "vmatch",
    "clustalo",
    "prodigal",
    "RNAfold",
    "blastn",
    "x86_64-conda-linux-gnu-cc",
)
EXECUTABLE_VERSION_PROBES = {
    "mkvtree": (("-version",), "mkvtree (Vmatch) 2.3.0"),
    "vmatch": (("-version",), "vmatch (Vmatch) 2.3.0"),
    "clustalo": (("--version",), "1.2.4"),
    "prodigal": (("-v",), "Prodigal V2.6.3"),
    "RNAfold": (("--version",), "RNAfold 2.7.0"),
    "hmmsearch": (("-h",), "HMMER 3.4"),
    "blastn": (("-version",), "blastn: 2.5.0+"),
    "x86_64-conda-linux-gnu-cc": (
        ("--version",), "(conda-forge gcc 14.4.0-1) 14.4.0"
    ),
}
RESOURCE_SHA256 = {
    "vendor/casandra/src/casandra/models/protein/profiles.hmm": (
        "034841fc34bc11e2fac675858a530bbe03ac02134cbd8d018947e67c5a5906b7"
    ),
    "vendor/casandra/src/casandra/models/cassette_architecture/model.joblib": (
        "afaf3cd7014b7221abb2b36c389c8f8e562bd7ebcff5728ac0c0ec5d18f74e9a"
    ),
    "tools/hmm_search/models_is_element.hmm": (
        "dd2bb9bb113d44646c1a84ae42210cfc8c144ae6fb45867d05fc62ad09602e5b"
    ),
}
CASANDRA_INSTALLED_RESOURCE_SHA256 = {
    "models/protein/profiles.hmm": RESOURCE_SHA256[
        "vendor/casandra/src/casandra/models/protein/profiles.hmm"
    ],
    "models/cassette_architecture/model.joblib": RESOURCE_SHA256[
        "vendor/casandra/src/casandra/models/cassette_architecture/model.joblib"
    ],
}
BLAST_DATABASES = (
    "tools/blasting/Verified_repeats_dataset1.fa",
)
BLAST_SUFFIXES = (
    "", ".nhr", ".nin", ".nog", ".nsd", ".nsi", ".nsq",
)


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _environment_executable(name):
    sibling = Path(sys.executable).resolve().parent / name
    if sibling.is_file() and os.access(str(sibling), os.X_OK):
        return sibling
    discovered = shutil.which(name)
    return Path(discovered).resolve() if discovered else None


def _executable_version_text(path, arguments):
    """Return combined version output without invoking a shell."""

    completed = subprocess.run(
        [str(path)] + list(arguments),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        timeout=15,
        check=False,
    )
    return completed.stdout


def collect_environment_errors(profile="core"):
    """Return actionable failures for one supported runtime profile."""

    if profile not in ("core", "is-elements"):
        raise ValueError("unknown environment profile: {}".format(profile))

    errors = []
    if sys.version_info[:2] != (3, 12):
        errors.append(
            "Python 3.12 is required; found {}".format(
                ".".join(str(part) for part in sys.version_info[:3])
            )
        )

    for distribution, expected_version in sorted(CORE_PACKAGES.items()):
        try:
            observed_version = importlib_metadata.version(distribution)
        except importlib_metadata.PackageNotFoundError:
            errors.append("missing Python package: {}".format(distribution))
            continue
        if observed_version != expected_version:
            errors.append(
                "{} version {} is untested; expected {}".format(
                    distribution, observed_version, expected_version
                )
            )

    executables = list(CORE_EXECUTABLES)
    if profile == "is-elements":
        executables.extend(("hmmsearch", "fasta36"))
    for executable in executables:
        path = _environment_executable(executable)
        if path is None:
            errors.append("missing executable: {}".format(executable))
            continue
        probe = EXECUTABLE_VERSION_PROBES.get(executable)
        if probe is None:
            continue
        arguments, expected_marker = probe
        try:
            version_text = _executable_version_text(path, arguments)
        except (OSError, subprocess.SubprocessError) as error:
            errors.append(
                "could not check {} version: {}".format(executable, error)
            )
            continue
        if expected_marker not in version_text:
            errors.append(
                "{} version is untested; expected output containing {!r}"
                .format(executable, expected_marker)
            )

    required_hashes = [
        "vendor/casandra/src/casandra/models/protein/profiles.hmm",
        "vendor/casandra/src/casandra/models/cassette_architecture/model.joblib",
    ]
    if profile == "is-elements":
        required_hashes.append("tools/hmm_search/models_is_element.hmm")
    for relative_path in required_hashes:
        path = PROJECT_ROOT / relative_path
        if not path.is_file():
            errors.append("missing resource: {}".format(relative_path))
        elif _sha256(path) != RESOURCE_SHA256[relative_path]:
            errors.append("resource checksum mismatch: {}".format(relative_path))

    for database in BLAST_DATABASES:
        for suffix in BLAST_SUFFIXES:
            relative_path = database + suffix
            if not (PROJECT_ROOT / relative_path).is_file():
                errors.append("missing BLAST database file: {}".format(relative_path))

    try:
        from components import literal_fuzzy
        if not literal_fuzzy.AVAILABLE:
            errors.append(
                "native fuzzy matcher unavailable: {}".format(
                    literal_fuzzy.LOAD_ERROR or "unknown loader error"
                )
            )
    except Exception as error:  # pragma: no cover - diagnostic boundary
        errors.append("native fuzzy matcher check failed: {}".format(error))

    try:
        from components.sequence_path_vectorizer import (
            NATIVE_ASCII_PATH_FEATURES_AVAILABLE,
            NATIVE_PACKED_PATH_FEATURES_AVAILABLE,
            NATIVE_PATH_FEATURE_COUNTS_AVAILABLE,
            NATIVE_PATH_FEATURES_AVAILABLE,
        )
        if not NATIVE_PATH_FEATURES_AVAILABLE:
            errors.append("native path-feature accelerator unavailable")
        elif not NATIVE_PATH_FEATURE_COUNTS_AVAILABLE:
            errors.append(
                "native path-feature count accelerator unavailable"
            )
        elif not NATIVE_PACKED_PATH_FEATURES_AVAILABLE:
            errors.append(
                "native packed path-feature accelerator unavailable"
            )
        elif not NATIVE_ASCII_PATH_FEATURES_AVAILABLE:
            errors.append(
                "native ASCII path-feature accelerator unavailable"
            )
    except Exception as error:  # pragma: no cover - diagnostic boundary
        errors.append(
            "native path-feature accelerator check failed: {}".format(error)
        )

    try:
        import numpy as np
        from components.feature_schema import COMPACT_FEATURE_SCHEMA
        from components.compact_model import predict_compact_scores

        scores = predict_compact_scores(
            np.zeros((1, len(COMPACT_FEATURE_SCHEMA)), dtype=np.float64)
        )
        if (
                scores.shape != (1,)
                or not np.all(np.isfinite(scores))
                or not np.all((scores >= 0.0) & (scores <= 1.0))):
            raise ValueError("compact model returned an invalid score")
    except Exception as error:  # pragma: no cover - diagnostic boundary
        errors.append("compact model check failed: {}".format(error))

    try:
        import json
        import joblib
        import casandra
        from casandra.production_v2 import load_model

        if casandra.__version__ != "0.1.0.dev0+crispridentify.3":
            raise ValueError("unexpected CasAndra runtime version")
        model_root = Path(casandra.__file__).resolve().parent / "models"
        manifest_path = model_root / "manifest.json"
        bundle_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if bundle_manifest.get("schema_version") != 4:
            raise ValueError("CasAndra bundle schema is not version 4")
        if bundle_manifest.get("bundle_id") != "casandra-cas-only-cpu-bundle-v4-fast-final":
            raise ValueError("unexpected CasAndra model bundle")
        package_root = Path(casandra.__file__).resolve().parent
        for relative_path, expected_hash in CASANDRA_INSTALLED_RESOURCE_SHA256.items():
            installed_path = package_root / relative_path
            if not installed_path.is_file():
                raise ValueError(
                    "installed CasAndra resource is missing: {}".format(
                        relative_path
                    )
                )
            if _sha256(installed_path) != expected_hash:
                raise ValueError("installed CasAndra checksum mismatch: {}".format(relative_path))
        _, _, metadata, hmms = load_model(model_root / "protein")
        if len(metadata) != 341 or len(hmms) != 341:
            raise ValueError("CasAndra protein profile bundle is incomplete")
        architecture_model = joblib.load(
            model_root / "cassette_architecture" / "model.joblib"
        )
        classifier = architecture_model.named_steps.get("extratreesclassifier")
        if len(getattr(classifier, "estimators_", ())) != 500:
            raise ValueError("CasAndra cassette model canary changed")
    except Exception as error:  # pragma: no cover - diagnostic boundary
        errors.append("CasAndra check failed: {}".format(error))

    try:
        from components.portable_strand import load_default_portable_strand_model

        model = load_default_portable_strand_model()
        labels = model.predict_strands(
            ["GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"]
        )
        if labels != ["Forward"]:
            raise ValueError("portable strand canary changed: {!r}".format(labels))
    except Exception as error:  # pragma: no cover - diagnostic boundary
        errors.append("portable strand check failed: {}".format(error))

    return errors


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="validate CRISPRidentify dependencies and model artifacts"
    )
    parser.add_argument(
        "--profile", choices=("core", "is-elements"), default="core"
    )
    args = parser.parse_args(argv)

    errors = collect_environment_errors(args.profile)
    if errors:
        print("CRISPRidentify environment check failed:", file=sys.stderr)
        for error in errors:
            print("  - {}".format(error), file=sys.stderr)
        return 1

    print(
        "CRISPRidentify {} environment: OK (Python {})".format(
            args.profile,
            ".".join(str(part) for part in sys.version_info[:3]),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
