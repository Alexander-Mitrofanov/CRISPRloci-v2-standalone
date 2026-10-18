from __future__ import annotations

import ast
import importlib
import importlib.util
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1] / "src" / "crispr_web"


def _module_name(path: Path) -> tuple[str, bool]:
    relative = path.relative_to(PACKAGE_ROOT).with_suffix("")
    parts = ["crispr_web", *relative.parts]
    is_package = parts[-1] == "__init__"
    if is_package:
        parts.pop()
    return ".".join(parts), is_package


def _runtime_modules() -> dict[str, tuple[Path, bool]]:
    return {
        module: (path, is_package)
        for path in PACKAGE_ROOT.rglob("*.py")
        for module, is_package in [_module_name(path)]
    }


def _import_targets(
    module: str,
    path: Path,
    is_package: bool,
    known_modules: set[str],
) -> set[str]:
    package = module if is_package else module.rpartition(".")[0]
    targets: set[str] = set()
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            targets.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                dotted = "." * node.level + (node.module or "")
                target = importlib.util.resolve_name(dotted, package)
            else:
                target = node.module or ""
            if target:
                targets.add(target)
            if node.module is None:
                targets.update(
                    candidate
                    for alias in node.names
                    if (candidate := f"{target}.{alias.name}") in known_modules
                )
    return targets


def _dependency_graph() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    modules = _runtime_modules()
    known = set(modules)
    all_imports: dict[str, set[str]] = {}
    internal: dict[str, set[str]] = {}
    for module, (path, is_package) in modules.items():
        targets = _import_targets(module, path, is_package, known)
        all_imports[module] = targets
        internal[module] = {target for target in targets if target in known}
    return all_imports, internal


def _assert_acyclic(graph: dict[str, set[str]]) -> None:
    visiting: list[str] = []
    visited: set[str] = set()

    def visit(module: str) -> None:
        if module in visited:
            return
        if module in visiting:
            start = visiting.index(module)
            cycle = " -> ".join([*visiting[start:], module])
            raise AssertionError(f"internal import cycle: {cycle}")
        visiting.append(module)
        for dependency in sorted(graph[module]):
            visit(dependency)
        visiting.pop()
        visited.add(module)

    for module in sorted(graph):
        visit(module)


def test_runtime_import_boundaries_and_graph() -> None:
    all_imports, internal = _dependency_graph()

    for module, targets in all_imports.items():
        assert not any(
            target == "crispr_integration" or target.startswith("crispr_integration.")
            for target in targets
        ), module
        if module.startswith("crispr_web.http"):
            assert not any(
                target.startswith(("crispr_web.worker", "crispr_web.execution"))
                for target in targets
            ), module
        if module.startswith("crispr_web.execution"):
            assert not any(
                target.startswith("crispr_web.http") for target in targets
            ), module

    _assert_acyclic(internal)


def test_public_compatibility_facades() -> None:
    expected = {
        "crispr_web.artifacts": {
            "_ArchiveLimitExceeded",
            "_BUNDLE_SUFFIXES",
            "_CappedArchiveWriter",
            "_DENIED_PARTS",
            "_EMBEDDED_ABSOLUTE_PATH",
            "_IDENTIFY_FILENAMES",
            "_INDIVIDUAL_EXACT",
            "_INTERNAL_JSON",
            "_MAX_PUBLIC_JSON_BYTES",
            "_MEDIA_TYPES",
            "_OMIT",
            "_PRIVATE_MANIFEST_KEYS",
            "_STATIC_PRIVATE_ROOTS",
            "_allowed_individual",
            "_bundle_outputs",
            "_checkpoint",
            "_contains_private_path",
            "_fsync_published",
            "_private_markers",
            "_public_json_safe",
            "_public_text_safe",
            "_read_json_bounded",
            "_record",
            "_safe_file",
            "_sanitize_manifest_value",
            "_write_sanitized_manifest",
            "_write_zip",
            "collect_artifacts",
            "file_sha256",
            "resolve_registered_file",
        },
        "crispr_web.api": {
            "API_VERSION",
            "JsonErrorMiddleware",
            "RequestLimitMiddleware",
            "SecurityHeadersMiddleware",
            "create_app",
        },
        "crispr_web.summary": {
            "_MAX_PUBLIC_SUMMARY_BYTES",
            "build_summary",
            "sanitize_public_summary",
        },
        "crispr_web.worker": {
            "ProcessContainmentError",
            "StorageCleanupError",
            "Worker",
            "_accounted_storage_bytes",
            "_adapter_group_names",
            "_job_tree_within_budget",
            "_supported_schema_major",
            "main",
        },
    }

    for module_name, names in expected.items():
        module = importlib.import_module(module_name)
        assert set(module.__all__) == names, f"{module_name} changed its stable __all__"
        missing = sorted(name for name in names if not hasattr(module, name))
        assert not missing, f"{module_name} lost compatibility exports: {missing}"
