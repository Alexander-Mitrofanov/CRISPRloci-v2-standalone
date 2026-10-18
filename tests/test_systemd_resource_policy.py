from __future__ import annotations

from pathlib import Path


DEPLOY = Path(__file__).resolve().parents[1] / "deploy"


def source(relative: str) -> str:
    return (DEPLOY / relative).read_text(encoding="utf-8")


def unit_section(relative: str, name: str) -> str:
    document = source(relative)
    marker = f"[{name}]\n"
    assert document.count(marker) == 1
    return document.split(marker, 1)[1].split("\n[", 1)[0]


def assert_directives(section: str, directives: tuple[str, ...]) -> None:
    lines = section.splitlines()
    for directive in directives:
        assert lines.count(directive) == 1, directive


def test_services_bound_restart_log_swap_and_io_pressure() -> None:
    services = (
        "systemd/crisprloci-api.service",
        "systemd/crisprloci-worker.service",
        "systemd/crisprloci-cleanup.service",
    )
    for relative in services:
        assert_directives(
            unit_section(relative, "Unit"),
            ("StartLimitIntervalSec=5min", "StartLimitBurst=3"),
        )
        assert_directives(
            unit_section(relative, "Service"),
            (
                "LogRateLimitIntervalSec=1min",
                "LogRateLimitBurst=60",
                "MemoryAccounting=yes",
                "MemorySwapMax=0",
                "MemoryZSwapMax=0",
                "IOAccounting=yes",
            ),
        )


def test_api_and_worker_have_explicit_shared_storage_budgets() -> None:
    api = unit_section("systemd/crisprloci-api.service", "Service")
    assert_directives(
        api,
        (
            "IOWeight=50",
            "IOReadBandwidthMax=/srv/crispr 64M",
            "IOWriteBandwidthMax=/srv/crispr 16M",
            "IOReadIOPSMax=/srv/crispr 1000",
            "IOWriteIOPSMax=/srv/crispr 500",
        ),
    )
    assert "Bandwidth suffixes are decimal SI units as defined by systemd." in api
    assert "1000/500 read/write IOPS ceilings" in api

    worker = unit_section("systemd/crisprloci-worker.service", "Service")
    assert_directives(
        worker,
        (
            "IOWeight=25",
            "IOReadBandwidthMax=/srv/crispr 128M",
            "IOWriteBandwidthMax=/srv/crispr 64M",
            "IOReadIOPSMax=/srv/crispr 2000",
            "IOWriteIOPSMax=/srv/crispr 1000",
            "TasksMax=1024",
        ),
    )
    assert "Bound both bulk throughput and small-file pressure" in worker
    assert "Bandwidth suffixes are decimal SI units as defined by systemd." in worker
    assert "2000/1000 read/write IOPS ceilings" in worker


def test_cleanup_is_bounded_and_runs_at_low_resource_priority() -> None:
    cleanup = unit_section("systemd/crisprloci-cleanup.service", "Service")
    assert_directives(
        cleanup,
        (
            "TimeoutStartSec=20min",
            "Nice=10",
            "CPUAccounting=yes",
            "CPUWeight=10",
            "IOWeight=10",
            "MemoryMax=512M",
            "TasksMax=64",
        ),
    )

