"""Synthetic network fixtures only: never download viral sequence data in tests."""
import gzip
import hashlib
import io
import json
from types import SimpleNamespace

import pytest

from crisprspacer import ncbi_viral as ncbi
from crisprspacer.io import InputError, LimitError, digest


class Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url
        self.status = 200

    def geturl(self):
        return self.url


@pytest.fixture
def remote(monkeypatch):
    sources = {"viral.1.1.genomic.fna.gz": gzip.compress(b">NC_TEST1.1 synthetic record\nACGTACGT\n", mtime=0),
               "viral.2.1.genomic.fna.gz": gzip.compress(b">NC_TEST2.1 synthetic record\nTGCATGCA\n", mtime=0)}
    catalog_url = ncbi.BASE + "release-catalog/release237.files.installed"
    catalog = b"0" * 32 + b"\tunrelated.protein.faa.gz\n"
    catalog += "".join(hashlib.md5(data, usedforsecurity=False).hexdigest() + "\t" + name + "\n"
                        for name, data in sources.items()).encode()
    listing = "<html>" + "".join(f'<a href="{name}">{name}</a>' for name in sources)
    listing += '<a href="viral.1.genomic.gbff.gz">unused annotation</a></html>'
    payloads = {ncbi.RELEASE_URL: b"237\n", ncbi.VIRAL_URL: listing.encode(), catalog_url: catalog,
                **{ncbi.VIRAL_URL + name: data for name, data in sources.items()}}
    calls = []

    def open_fixture(url):
        calls.append(url)
        return Response(payloads[url], url)

    monkeypatch.setattr(ncbi, "_open", open_fixture)
    return payloads, calls, sources, catalog_url


def make_plan(tmp_path):
    path = tmp_path / "plan.json"
    return path, ncbi.plan_ncbi_viral(path)


def test_plan_is_metadata_only_complete_and_pinned(tmp_path, remote):
    payloads, calls, sources, catalog_url = remote
    path, plan = make_plan(tmp_path)
    assert json.loads(path.read_bytes()) == plan
    assert [row["name"] for row in plan["files"]] == sorted(sources)
    assert all(row["bytes"] is None for row in plan["files"])
    assert set(calls) == {ncbi.RELEASE_URL, ncbi.VIRAL_URL, catalog_url}
    assert len(calls) == 6
    assert plan["catalog"]["sha256"] == digest(payloads[catalog_url])
    assert plan["source_release"] == "237"
    assert plan["schema"] == ncbi.PLAN_SCHEMA


@pytest.mark.parametrize("changed", ["release", "listing", "catalog"])
def test_plan_fails_when_metadata_changes(tmp_path, remote, monkeypatch, changed):
    payloads, _, _, catalog_url = remote
    target = {"release": ncbi.RELEASE_URL, "listing": ncbi.VIRAL_URL, "catalog": catalog_url}[changed]
    visits = 0

    def unstable(url):
        nonlocal visits
        if url == target:
            visits += 1
        data = payloads[url]
        if url == target and visits > 1:
            data += b"\n"
        return Response(data, url)

    monkeypatch.setattr(ncbi, "_open", unstable)
    with pytest.raises(InputError, match="changed"):
        make_plan(tmp_path)
    assert not (tmp_path / "plan.json").exists()


def test_plan_rejects_catalog_listing_disagreement(tmp_path, remote):
    payloads, _, _, _ = remote
    payloads[ncbi.VIRAL_URL] = b'<a href="viral.1.1.genomic.fna.gz">one only</a>'
    with pytest.raises(InputError, match="disagree"):
        make_plan(tmp_path)


def test_download_receipt_binds_all_sources_and_original_plan(tmp_path, remote):
    payloads, calls, sources, catalog_url = remote
    plan_path, plan = make_plan(tmp_path)
    calls.clear()
    output = tmp_path / "sources"
    receipt = ncbi.download_ncbi_viral(plan_path, output, min_free_bytes=0)
    assert receipt["schema"] == ncbi.DOWNLOAD_SCHEMA
    assert (output / "plan.json").read_bytes() == plan_path.read_bytes()
    assert (output / "release.files.installed").read_bytes() == payloads[catalog_url]
    assert receipt["plan_sha256"] == digest(plan_path.read_bytes())
    assert json.loads((output / "download.json").read_bytes()) == receipt
    for row in receipt["files"]:
        data = sources[row["name"]]
        assert (output / row["name"]).read_bytes() == data
        assert row["bytes"] == len(data)
        assert row["sha256"] == digest(data)
        assert row["md5"] == hashlib.md5(data, usedforsecurity=False).hexdigest()
    assert calls.count(catalog_url) == 2
    assert ncbi.RELEASE_URL not in calls  # Historical plans need their pinned catalog, not today's release.
    assert receipt["source_release"] == plan["source_release"]


@pytest.mark.parametrize("field,value", [("url", "http://127.0.0.1/private"),
                                         ("name", "../../outside"), ("md5", "not-a-hash")])
def test_untrusted_plan_cannot_select_url_or_path(tmp_path, remote, field, value):
    _, calls, _, _ = remote
    plan_path, plan = make_plan(tmp_path)
    plan["files"][0][field] = value
    plan_path.write_text(json.dumps(plan))
    calls.clear()
    with pytest.raises(InputError):
        ncbi.download_ncbi_viral(plan_path, tmp_path / "sources", min_free_bytes=0)
    assert calls == []
    assert not (tmp_path / "sources").exists()


def test_plan_cannot_omit_catalog_sources(tmp_path, remote):
    _, calls, _, catalog_url = remote
    plan_path, plan = make_plan(tmp_path)
    plan["files"] = plan["files"][:1]
    plan_path.write_text(json.dumps(plan))
    calls.clear()
    with pytest.raises(InputError, match="inventory changed"):
        ncbi.download_ncbi_viral(plan_path, tmp_path / "sources", min_free_bytes=0)
    assert calls == [catalog_url]


@pytest.mark.parametrize("corruption", ["catalog", "source"])
def test_changed_source_never_publishes_partial_directory(tmp_path, remote, corruption):
    payloads, _, sources, catalog_url = remote
    plan_path, _ = make_plan(tmp_path)
    key = catalog_url if corruption == "catalog" else ncbi.VIRAL_URL + next(iter(sources))
    payloads[key] += b"tampered"
    output = tmp_path / "sources"
    with pytest.raises(InputError, match="catalog|MD5"):
        ncbi.download_ncbi_viral(plan_path, output, min_free_bytes=0)
    assert not output.exists()
    assert not (tmp_path / "sources.partial").exists()
    assert not list(tmp_path.glob(".crisprspacer-*"))


def test_valid_official_checksum_does_not_replace_gzip_validation(tmp_path, remote):
    payloads, _, sources, catalog_url = remote
    name = next(iter(sources))
    corrupt = b"not gzip"
    old_md5 = hashlib.md5(sources[name], usedforsecurity=False).hexdigest().encode()
    new_md5 = hashlib.md5(corrupt, usedforsecurity=False).hexdigest().encode()
    payloads[catalog_url] = payloads[catalog_url].replace(old_md5, new_md5)
    payloads[ncbi.VIRAL_URL + name] = corrupt
    plan_path, _ = make_plan(tmp_path)
    with pytest.raises(InputError, match="gzip"):
        ncbi.download_ncbi_viral(plan_path, tmp_path / "sources", min_free_bytes=0)
    assert not (tmp_path / "sources").exists()


@pytest.mark.parametrize("budget", ["compressed", "decompressed", "disk"])
def test_download_resource_bounds(tmp_path, remote, monkeypatch, budget):
    plan_path, _ = make_plan(tmp_path)
    if budget == "compressed":
        monkeypatch.setattr(ncbi, "MAX_SOURCE_BYTES", 1)
    elif budget == "decompressed":
        monkeypatch.setattr(ncbi, "MAX_UNCOMPRESSED_BYTES", 1)
    else:
        monkeypatch.setattr(ncbi.shutil, "disk_usage", lambda _: SimpleNamespace(free=0))
    with pytest.raises(LimitError):
        ncbi.download_ncbi_viral(plan_path, tmp_path / "sources", min_free_bytes=1)
    assert not (tmp_path / "sources").exists()


def test_metadata_budget_is_enforced_before_parsing(remote):
    with pytest.raises(LimitError, match="metadata"):
        ncbi._fetch(ncbi.VIRAL_URL, 10)


def test_official_catalog_format_and_duplicate_rejection():
    line = b"104124827b73be6bfd87af42596bf4be\tviral.1.1.genomic.fna.gz\n"
    assert ncbi._catalog(line)[0]["md5"] == "104124827b73be6bfd87af42596bf4be"
    for data in (line + line, b"x\tviral.1.1.genomic.fna.gz", b"extra " + line):
        with pytest.raises(InputError, match="checksum"):
            ncbi._catalog(data)


def test_no_arbitrary_network_urls_or_redirects():
    assert ncbi._allowed_url(ncbi.VIRAL_URL + "viral.1.1.genomic.fna.gz")
    for url in ("http://ftp.ncbi.nlm.nih.gov/refseq/release/viral/", ncbi.VIRAL_URL + "../private",
                ncbi.VIRAL_URL + "viral.1.1.genomic.fna.gz?x=1", "https://127.0.0.1/"):
        assert not ncbi._allowed_url(url)
    with pytest.raises(InputError, match="redirects"):
        ncbi._NoRedirect().redirect_request(None, None, 302, "", {}, "https://127.0.0.1/")


def test_existing_plan_and_destination_preserved(tmp_path, remote):
    plan_path, _ = make_plan(tmp_path)
    saved = plan_path.read_bytes()
    with pytest.raises(InputError, match="already exists"):
        ncbi.plan_ncbi_viral(plan_path)
    assert plan_path.read_bytes() == saved
    output = tmp_path / "sources"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("keep")
    with pytest.raises(InputError, match="already exists"):
        ncbi.download_ncbi_viral(plan_path, output, min_free_bytes=0)
    assert marker.read_text() == "keep"
