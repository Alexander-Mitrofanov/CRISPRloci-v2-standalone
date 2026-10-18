import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT
VENDOR = V2 / "vendor" / "casandra"
MODEL_ROOT = VENDOR / "src" / "casandra" / "models"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_artifact_manifest(test, directory):
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    artifacts = manifest.get("artifacts")
    test.assertIsInstance(artifacts, dict)
    test.assertTrue(artifacts)
    for relative_name, record in artifacts.items():
        relative = PurePosixPath(relative_name)
        test.assertFalse(relative.is_absolute())
        test.assertNotIn("..", relative.parts)
        path = directory.joinpath(*relative.parts)
        test.assertTrue(path.is_file(), relative_name)
        test.assertEqual(record["size"], path.stat().st_size, relative_name)
        test.assertEqual(record["sha256"], sha256(path), relative_name)
    return manifest


class V2DistributionTests(unittest.TestCase):
    def test_native_matcher_build_normalizes_shared_library_permissions(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "release"
            scripts = project / "scripts"
            native = project / "components" / "native"
            scripts.mkdir(parents=True)
            native.mkdir(parents=True)

            build = scripts / "build_native_matcher.sh"
            build.write_text(
                (V2 / "scripts" / "build_native_matcher.sh").read_text(
                    encoding="utf-8"
                ),
                encoding="utf-8",
            )
            build.chmod(0o755)
            (native / "literal_fuzzy.c").write_text("fixture\n")
            (native / "path_features.c").write_text("fixture\n")

            compiler = project / "fake-cc"
            compiler.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "output=''\n"
                "while (($#)); do\n"
                "  if [[ \"$1\" == '-o' ]]; then\n"
                "    output=\"$2\"\n"
                "    shift 2\n"
                "    continue\n"
                "  fi\n"
                "  shift\n"
                "done\n"
                "printf 'fake shared library\\n' > \"${output}\"\n"
                "chmod 0600 \"${output}\"\n",
                encoding="utf-8",
            )
            compiler.chmod(0o755)
            environment = dict(os.environ)
            environment["CC"] = str(compiler)
            environment.pop("CONDA_PREFIX", None)

            subprocess.run(
                [str(build)],
                check=True,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            library = native / "libliteral_fuzzy.so"
            self.assertTrue(library.is_file())
            self.assertEqual(0o755, stat.S_IMODE(library.stat().st_mode))
            self.assertEqual([], list(native.glob(".libliteral_fuzzy.*.so")))

    def test_v2_is_a_complete_release_tree_and_legacy_cas_is_absent(self):
        self.assertEqual("2.3.2", (V2 / "VERSION").read_text().strip())
        ignore = (V2 / ".gitignore").read_text(encoding="utf-8")
        self.assertNotIn("vendor/", ignore)
        self.assertNotIn("trained_models/", ignore)

        required = (
            V2 / "CRISPRidentify.py",
            V2 / "environment.yml",
            V2 / "locks" / "conda-linux-64-core.lock",
            V2 / "scripts" / "setup_local.sh",
            V2 / "scripts" / "run_crispridentify.sh",
            V2 / "scripts" / "run_crispridentify_batch.py",
            V2 / "install.sh",
            V2 / "crispridentify",
            V2 / "crispridentify-batch",
            V2 / "docs" / "CLUSTER.md",
            V2 / ".github" / "workflows" / "ci.yml",
            V2 / "tests" / "__init__.py",
            V2 / "components" / "casandra_adapter.py",
            VENDOR / "VENDORING.json",
            VENDOR / "SOURCE_MANIFEST.sha256",
            MODEL_ROOT / "protein" / "profiles.hmm",
            MODEL_ROOT / "cassette_architecture" / "model.joblib",
        )
        for path in required:
            self.assertTrue(path.is_file(), str(path))

        retired = (
            V2 / "components" / "cas_prediction.py",
            V2 / "components" / "cas_model.py",
            V2 / "trained_models" / "cas",
            V2 / "licenses" / "CRISPRcasIdentifier-GPL-3.0.txt",
        )
        for path in retired:
            self.assertFalse(path.exists(), str(path))

        generated_suffixes = (".orig", ".rej", ".egg-info")
        generated_names = {"build", "dist", ".pytest_cache", "CasAndra.egg-info"}
        for path in V2.rglob("*"):
            relative_parts = path.relative_to(V2).parts
            if (
                    relative_parts[0] in {".envs", ".git"}
                    or "__pycache__" in path.parts
                    or path.name.endswith((".pyc", ".pyo"))):
                continue
            if path.is_file():
                self.assertLess(
                    path.stat().st_size,
                    100 * 1024 * 1024,
                    "file exceeds GitHub's ordinary per-file limit: {}".format(path),
                )
            self.assertNotIn(path.name, generated_names)
            self.assertFalse(path.name.endswith(generated_suffixes), str(path))


    def test_release_modes_and_large_models_are_plain_git_files(self):
        launchers = {
            V2 / "install.sh",
            V2 / "crispridentify",
            V2 / "crispridentify-batch",
            V2 / "scripts" / "setup_local.sh",
            V2 / "scripts" / "run_crispridentify.sh",
            V2 / "scripts" / "build_native_matcher.sh",
        }
        for path in V2.rglob("*"):
            if not path.is_file():
                continue
            relative_parts = path.relative_to(V2).parts
            if (
                    relative_parts[0] in {".envs", ".git"}
                    or "__pycache__" in relative_parts
                    or path.name.endswith((".pyc", ".pyo", ".so"))):
                continue
            executable = bool(stat.S_IMODE(path.stat().st_mode) & 0o111)
            self.assertEqual(path in launchers, executable, str(path))

        models = {
            V2 / "trained_models" / "compact_model.npz": 50_000,
            V2 / "trained_models" / "strand" / "crisprstrand_model_r_v1.npz": 100_000,
            MODEL_ROOT / "protein" / "profiles.hmm": 50_000_000,
            MODEL_ROOT / "cassette_architecture" / "model.joblib": 4_000_000,
        }
        lfs_prefix = b"version https://git-lfs.github.com/spec/v1"
        for path, minimum_size in models.items():
            self.assertGreater(path.stat().st_size, minimum_size, str(path))
            self.assertFalse(path.read_bytes()[:128].startswith(lfs_prefix), str(path))

    def test_convenience_launcher_preserves_caller_relative_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            caller = root / "caller with spaces"
            env_dir = root / "fake environment"
            caller.mkdir()
            (env_dir / "bin").mkdir(parents=True)
            fake_python = env_dir / "bin" / "python"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'CWD=%s\\n' \"$PWD\"\n"
                "for argument in \"$@\"; do "
                "printf 'ARG=%s\\n' \"$argument\"; done\n",
                encoding="utf-8",
            )
            fake_python.chmod(0o755)
            environment = os.environ.copy()
            environment["CRISPRIDENTIFY_ENV"] = str(env_dir)
            absolute_manifest = root / "absolute topology.tsv"
            command = [
                str(V2 / "crispridentify"),
                "--file", "input file.fasta",
                "--input_folder=input directory",
                "--input_folder_multifasta", "multifasta directory",
                "--result_folder=output directory",
                "--pickle_report", "pickle directory",
                "--json_report=json directory",
                "--sequence_topology_manifest", str(absolute_manifest),
                "--cas", "True",
            ]
            completed = subprocess.run(
                command,
                cwd=caller,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                [
                    "CWD={}".format(V2),
                    "ARG={}".format(V2 / "CRISPRidentify.py"),
                    "ARG=--file", "ARG={}".format(caller / "input file.fasta"),
                    "ARG=--input_folder={}".format(caller / "input directory"),
                    "ARG=--input_folder_multifasta",
                    "ARG={}".format(caller / "multifasta directory"),
                    "ARG=--result_folder={}".format(caller / "output directory"),
                    "ARG=--pickle_report",
                    "ARG={}".format(caller / "pickle directory"),
                    "ARG=--json_report={}".format(caller / "json directory"),
                    "ARG=--sequence_topology_manifest",
                    "ARG={}".format(absolute_manifest),
                    "ARG=--cas", "ARG=True",
                ],
                completed.stdout.splitlines(),
            )

            default_result = subprocess.run(
                [str(V2 / "crispridentify"), "--file", "input.fasta"],
                cwd=caller,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.splitlines()
            self.assertEqual(
                ["ARG=--result_folder", "ARG={}".format(caller / "Results")],
                default_result[-2:],
            )

    def test_casandra_runtime_and_nested_model_manifests_are_exact(self):
        top = assert_artifact_manifest(self, MODEL_ROOT)
        self.assertEqual(4, top["schema_version"])
        self.assertEqual(
            "casandra-cas-only-cpu-bundle-v4-fast-final",
            top["bundle_id"],
        )
        protein = assert_artifact_manifest(self, MODEL_ROOT / "protein")
        self.assertEqual(
            "casandra-production-profile-hmm-type-v-v4-fast",
            protein["model_id"],
        )
        architecture = assert_artifact_manifest(
            self, MODEL_ROOT / "cassette_architecture"
        )
        self.assertNotIn(
            "production_constructed_cassettes.jsonl",
            architecture["artifacts"],
        )
        self.assertIn(
            "production_constructed_cassettes.jsonl",
            architecture["omitted_training_artifacts"],
        )

    def test_vendored_source_manifest_covers_every_release_file(self):
        manifest_path = VENDOR / "SOURCE_MANIFEST.sha256"
        records = {}
        for line_number, line in enumerate(
                manifest_path.read_text(encoding="utf-8").splitlines(), 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, separator, relative_name = line.partition("  ")
            self.assertEqual("  ", separator, line_number)
            self.assertEqual(64, len(digest), line_number)
            self.assertTrue(
                all(character in "0123456789abcdef" for character in digest),
                line_number,
            )
            relative = PurePosixPath(relative_name)
            self.assertFalse(relative.is_absolute(), line_number)
            self.assertNotIn("..", relative.parts, line_number)
            self.assertNotIn(relative_name, records)
            records[relative_name] = digest
        self.assertTrue(records)

        actual = {
            path.relative_to(VENDOR).as_posix()
            for path in VENDOR.rglob("*")
            if (
                path.is_file()
                and path != manifest_path
                and "__pycache__" not in path.parts
                and not path.name.endswith((".pyc", ".pyo"))
            )
        }
        self.assertEqual(actual, set(records))
        for relative_name, digest in records.items():
            self.assertEqual(digest, sha256(VENDOR / relative_name), relative_name)

    def test_setup_and_lock_pin_the_validated_offline_stack(self):
        environment = (V2 / "environment.yml").read_text(encoding="utf-8")
        is_environment = (V2 / "environments" / "is-elements.yml").read_text(encoding="utf-8")
        lock = (V2 / "locks" / "conda-linux-64-core.lock").read_text(
            encoding="utf-8"
        )
        pyproject = (VENDOR / "pyproject.toml").read_text(encoding="utf-8")
        exact = {
            "biopython": "1.86",
            "joblib": "1.5.3",
            "numpy": "1.26.4",
            "pyhmmer": "0.12.1",
            "pyrodigal": "3.7.1",
            "scikit-learn": "1.8.0",
            "scipy": "1.17.1",
        }
        for package, version in exact.items():
            for profile in (environment, is_environment):
                self.assertIn("- {}={}".format(package, version), profile)
            self.assertIn("{}=={}".format(package, version), pyproject)
            self.assertIn("/{}-{}-".format(package, version), lock)

        setup = (V2 / "scripts" / "setup_local.sh").read_text(encoding="utf-8")
        for token in (
                "micromamba", "env update", "--prune", "--no-index",
                "--no-deps", "--no-build-isolation", "pip check"):
            self.assertIn(token, setup)
        for profile in (environment, is_environment):
            self.assertIn("- c-compiler=1.11.0", profile)
            self.assertIn("- gcc_linux-64=14.4.0", profile)
        self.assertIn("/c-compiler-1.11.0-", lock)
        self.assertIn("/gcc_linux-64-14.4.0-", lock)
        lock_records = [line for line in lock.splitlines() if line.startswith("https://")]
        self.assertEqual(len(lock_records), 80)
        self.assertEqual(len(lock_records), len(set(lock_records)))
        for record in lock_records:
            url, separator, digest = record.partition("#")
            self.assertEqual(separator, "#")
            self.assertEqual(len(digest), 64)
            self.assertTrue(all(character in "0123456789abcdef" for character in digest))
            self.assertTrue(
                url.startswith((
                    "https://conda.anaconda.org/conda-forge/",
                    "https://conda.anaconda.org/bioconda/",
                ))
            )
        self.assertIn('host_os="$(uname -s)"', setup)
        self.assertIn('host_arch="$(uname -m)"', setup)
        self.assertLess(
            setup.index('host_os="$(uname -s)"'),
            setup.index('solver="$(command -v mamba)"'),
        )
        self.assertLess(
            setup.index('for hook in "${env_dir}"/etc/conda/activate.d/*.sh'),
            setup.index('"${project_dir}/scripts/build_native_matcher.sh"'),
        )
        native_build = (V2 / "scripts" / "build_native_matcher.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("x86_64-conda-linux-gnu-cc", native_build)
        self.assertIn("compiler resolved outside", native_build)

        doctor = (V2 / "scripts" / "check_environment.py").read_text(encoding="utf-8")
        self.assertIn("x86_64-conda-linux-gnu-cc", doctor)
        self.assertIn("(conda-forge gcc 14.4.0-1) 14.4.0", doctor)

if __name__ == "__main__":
    unittest.main()
