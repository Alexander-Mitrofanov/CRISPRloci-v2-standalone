"""Exercise installation choices, bounded downloads and failure/HTTP boundaries."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys
import threading
from types import SimpleNamespace
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from installer import download as downloads
from installer import engine as engines
from installer import plan
from installer import server as servers


@pytest.mark.parametrize("value", [
    {"command": "do something"}, {"cpu": True}, {"cpu": 0}, {"cpu": 33},
    {"developer": "yes"}, {"references": {"unknown": {"mode": "sample"}}},
    {"references": {"spacer": {"mode": "shell"}}},
    {"references": {"spacer": {"mode": "local", "path": "/missing-installer-input.fna"}}},
    {"references": {"spacer": {"mode": "sample", "command": "custom"}}},
])
def test_unknown_or_invalid_choices_never_become_commands(value):
    with pytest.raises((ValueError, TypeError)):
        plan.choices(value)


def test_local_references_require_regular_paths_and_import_digest(tmp_path):
    source = tmp_path / "source.fna"
    source.write_text(">a\nACGT\n")
    alias = tmp_path / "link.fna"
    alias.symlink_to(source)
    with pytest.raises(ValueError, match="symlink"):
        plan.choices({"references": {"spacer": {"mode": "local", "path": str(alias)}}})
    with pytest.raises(ValueError, match="SHA-256"):
        plan.choices({"references": {"catalog": {"mode": "import", "path": str(source)}}})
    assert plan.choices({"references": {"spacer": {"mode": "local", "path": str(source)}}})["references"]["spacer"]["path"] == str(source)
    assert plan.stages(plan.defaults())[-1]["id"] == "activate"


class Response(io.BytesIO):
    def __init__(self, value, length=True):
        super().__init__(value)
        self.headers = {"Content-Length": str(len(value))} if length else {}

    def geturl(self):
        return "https://example.org/download"


def test_download_reports_real_bytes_and_reuses_only_verified_cache(tmp_path, monkeypatch):
    data = b"reference" * 300_000
    monkeypatch.setattr(downloads, "urlopen", lambda *a, **k: Response(data))
    monkeypatch.setattr(downloads.shutil, "disk_usage", lambda _: SimpleNamespace(free=10 * 1024**3))
    events = []
    destination = tmp_path / "data.gz"
    expected = hashlib.sha256(data).hexdigest()
    downloads.download("https://example.org/file", destination, expected, lambda *event: events.append(event))
    assert destination.read_bytes() == data
    assert events[0][1:3] == (0, len(data))
    assert events[-1][1:3] == (len(data), len(data))
    assert any(0 < event[1] < len(data) for event in events)
    destination.write_bytes(b"changed")
    with pytest.raises(ValueError, match="Cached download"):
        downloads.download("https://example.org/file", destination, expected, lambda *a: None)
    assert destination.read_bytes() == b"changed"


@pytest.mark.parametrize("failure", ["checksum", "cancel", "limit"])
def test_failed_download_never_publishes_partial_bytes(tmp_path, monkeypatch, failure):
    data = b"reference data"
    monkeypatch.setattr(downloads, "urlopen", lambda *a, **k: Response(data))
    monkeypatch.setattr(downloads.shutil, "disk_usage", lambda _: SimpleNamespace(free=10 * 1024**3))
    destination = tmp_path / "data.gz"
    with pytest.raises((ValueError, InterruptedError)):
        downloads.download("https://example.org/file", destination,
                           "0" * 64 if failure == "checksum" else hashlib.sha256(data).hexdigest(),
                           lambda *a: None, lambda: failure == "cancel",
                           limit=1 if failure == "limit" else 100)
    assert not destination.exists() and not destination.with_suffix(".gz.partial").exists()


def test_failed_gate_cannot_activate_or_mark_installation_complete(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, "system", lambda *args: {"errors": []})
    engine = engines.Engine(tmp_path)
    seen = []

    def command(arguments, **kwargs):
        for stage in plan.stages(plan.defaults()):
            if stage["id"] in arguments:
                seen.append(stage["id"])
                if stage["id"] == "gate":
                    raise ValueError("Real-science gate failed")

    monkeypatch.setattr(engine, "command", command)
    monkeypatch.setattr(engine, "bootstrap", lambda: None)
    engine.start({})
    engine.thread.join(timeout=5)
    assert engine.snapshot()["status"] == "failed"
    assert "gate" in seen and "activate" not in seen
    assert not (tmp_path / ".installer/installation.json").exists()


def test_cancellation_stops_before_later_stages(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, "system", lambda *args: {"errors": []})
    engine = engines.Engine(tmp_path)
    seen = []

    def command(arguments, **kwargs):
        seen.append(arguments)
        engine.cancel()

    monkeypatch.setattr(engine, "command", command)
    engine.start({})
    engine.thread.join(timeout=5)
    assert engine.snapshot()["status"] == "cancelled" and len(seen) == 1


def test_subprocess_progress_is_structured_and_logs_are_bounded(tmp_path):
    engine = engines.Engine(tmp_path)
    code = "import json; print('CRISPR_PROGRESS '+json.dumps({'label':'Indexing','completed':2,'total':5,'unit':'records'})); [print('log '+str(i)) for i in range(180)]"
    engine.command([sys.executable, "-c", code])
    state = engine.snapshot()
    assert state["progress"] == {"label": "Indexing", "completed": 2, "total": 5, "unit": "records"}
    assert len(state["logs"]) == 150 and state["logs"][-1] == "log 179"


@pytest.fixture
def installer_server(tmp_path, monkeypatch):
    monkeypatch.setattr(servers, "system", lambda *args: {"errors": [], "free_bytes": 100 * 1024**3})
    server = servers.InstallerServer(root=tmp_path)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield server
    server.shutdown()
    server.server_close()
    thread.join(timeout=5)


def request(server, route, value=None, **headers):
    call = Request(server.origin + route, data=None if value is None else json.dumps(value).encode(),
                   headers={"Authorization": "Bearer " + server.token, "Content-Type": "application/json", **headers})
    try:
        return urlopen(call, timeout=5)
    except HTTPError as error:
        return error


def test_installer_api_rejects_cross_origin_bad_token_and_host(installer_server):
    server = installer_server
    with request(server, "/api/state", Authorization="Bearer wrong") as response:
        assert response.status == 401
    with request(server, "/api/install", {}, Origin="https://outside.example") as response:
        assert response.status == 403
    with request(server, "/api/state", Host="rebound.example") as response:
        assert response.status == 403
    with request(server, "/api/plan", {"command": "custom"}) as response:
        assert response.status == 400
    assert server.engine.snapshot()["status"] == "idle"


def test_wizard_plan_is_read_only_and_assets_cannot_escape(installer_server):
    server = installer_server
    with request(server, "/api/plan", {}) as response:
        body = json.load(response)
        assert response.status == 200
        assert body["stages"][-1]["id"] == "activate"
    assert not (server.root / ".installer").exists()
    with request(server, "/") as response:
        assert response.status == 200 and b"Reference data" in response.read()
        assert "frame-ancestors 'none'" in response.headers["Content-Security-Policy"]
    with request(server, "/../../.env") as response:
        assert response.status == 404


def test_local_launcher_does_not_inherit_deployment_configuration(tmp_path, monkeypatch):
    monkeypatch.setenv("CRISPR_WEB_DATA_ROOT", "/unrelated")
    monkeypatch.setenv("CRISPR_WEB_SPACER_DATABASE", "/unrelated/reference")
    monkeypatch.setenv("VITE_API_BASE_URL", "https://other.example")
    environment = engines.local_environment(tmp_path)
    assert not any(key.startswith(("CRISPR_WEB_", "VITE_")) for key in environment)
    assert environment["MAMBA_ROOT_PREFIX"] == str(tmp_path / ".installer/mamba")


def test_resume_space_check_is_tied_to_choices_and_exact_source_lock(tmp_path, monkeypatch):
    (tmp_path / "tools").mkdir()
    lock = tmp_path / "tools/tools.lock.json"
    lock.write_text('{"tools": []}')
    monkeypatch.setattr(plan.shutil, "disk_usage", lambda _: SimpleNamespace(free=10 * plan.GIB))
    value = plan.defaults()
    assert plan.system(tmp_path, value)["errors"]
    work = plan.work_directory(tmp_path, value)
    work.mkdir(parents=True)
    (work / "completed.json").write_text('["preflight","bootstrap"]')
    (tmp_path / ".installer/bootstrap.json").write_text('{}')
    assert plan.system(tmp_path, value)["resuming"]
    assert not plan.system(tmp_path, value)["errors"]
    lock.write_text('{"tools": ["changed"]}')
    assert not plan.system(tmp_path, value)["resuming"]


@pytest.mark.parametrize("invalid", ["modes", "catalog", None])
def test_candidate_api_checks_all_modes_and_catalogue_before_activation(tmp_path, monkeypatch, invalid):
    from installer import verify

    root = tmp_path
    (root / "frontend/src/contracts").mkdir(parents=True)
    (root / "frontend/src/contracts/module-ports-v1.json").write_text('{"detection":"target"}')
    release = root / "candidate"
    process = SimpleNamespace(poll=lambda: None, send_signal=lambda *_: None, wait=lambda **_: 0)
    monkeypatch.setattr(verify.subprocess, "Popen", lambda *a, **k: process)

    def response(url, **kwargs):
        if url.endswith("health"):
            payload = {"release_id": "candidate", "worker": "ok"}
        elif url.endswith("config"):
            payload = {"modes": [] if invalid == "modes" else ["detection"]}
        else:
            payload = {"available": invalid != "catalog"}
        return io.BytesIO(json.dumps(payload).encode())

    monkeypatch.setattr(verify, "urlopen", response)
    if invalid:
        with pytest.raises(ValueError):
            verify.verify_local(root, root, release, {"cpu": 2, "references": {}})
    else:
        verify.verify_local(root, root, release, {"cpu": 2, "references": {}})


def test_sample_catalogue_has_truthful_provenance_and_searchable_rows(tmp_path):
    import sqlite3
    from installer.references import sample_catalog
    from index_catalog import build_index, read_registry
    from crispr_web.catalog import Catalog

    sample_catalog(tmp_path)
    summary, shards, _ = read_registry(tmp_path)
    assert summary["counts"]["assemblies"] == 1
    assert summary["coverage"]["catalog_assemblies"] == 1
    assert "Synthetic" in summary["title"]
    assert "ClusterCRISPRCasDb" not in json.dumps(summary)
    assert Catalog(tmp_path / "catalog.sqlite3").summary().available
    assert len(shards) == 1
    result = build_index(tmp_path, workers=1, min_free=0)
    assert result["status"] == "complete"
    with sqlite3.connect(tmp_path / "exact-search.sqlite3") as connection:
        assert connection.execute("SELECT count(*) FROM hits").fetchone()[0] > 0
