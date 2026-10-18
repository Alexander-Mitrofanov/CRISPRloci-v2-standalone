from __future__ import annotations

import importlib
import json
from pathlib import Path

import pytest

from crispr_web.worker import _adapter_group_names, _supported_schema_major


@pytest.mark.parametrize("fixture_name", ["current_report.json", "future_report.json"])
def test_real_integration_adapter_matches_backend_major_contract(
    tmp_path, monkeypatch, fixture_name
):
    repository_root = Path(__file__).resolve().parents[2]
    integration_root = repository_root / "Integration"
    monkeypatch.syspath_prepend(str(integration_root))
    adapt_reports = importlib.import_module("crispr_integration.adapter").adapt_reports
    fixture = integration_root / "tests" / "fixtures" / fixture_name
    report = json.loads(fixture.read_text(encoding="utf-8"))
    source = report["source"]
    expected_fasta = tmp_path / "expected.fasta"
    expected_fasta.write_text(
        f">{source['id']} contract-test\nACGTACGTACGT\n", encoding="ascii"
    )

    manifest = adapt_reports(
        fixture,
        tmp_path / "adapter",
        expected_fasta_path=expected_fasta,
    )

    assert _supported_schema_major(report, "CRISPRidentify-report", 1)
    assert _supported_schema_major(manifest, "crispr-tool-integration-adapter", 1)
    assert manifest["source_binding"]["status"] == "verified"
    assert _adapter_group_names(manifest, manifest["emitted_group_count"]) is not None
    serialized = json.dumps(manifest)
    assert str(tmp_path) not in serialized
    assert source["sha256_scope"] == "uppercase_sequence"
