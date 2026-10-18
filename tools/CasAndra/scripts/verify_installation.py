"""Verify the installed wheel and exercise its CLI without workspace imports."""

from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

import casandra
from casandra.api import Analyzer, DEFAULT_MODEL_BUNDLE
from casandra.model import verify_artifacts
from casandra.runtime import parse_fasta


ROOT = Path(__file__).resolve().parents[1]


def run(*args, cwd):
    completed = subprocess.run(
        [sys.executable, "-m", "casandra", *map(str, args)],
        cwd=cwd, check=True, text=True, capture_output=True,
    )
    return json.loads(completed.stdout)


def main():
    package = Path(casandra.__file__).resolve().parent
    assert not package.is_relative_to(ROOT / "src"), "Install the wheel before verification"
    provenance = json.loads((ROOT / "provenance/source.json").read_text())
    assert casandra.__version__ == provenance["runtime_version"]
    for relative, digest in provenance["runtime_modules"].items():
        installed = package / Path(relative).relative_to("src/casandra")
        assert hashlib.sha256(installed.read_bytes()).hexdigest() == digest, relative
    assert hashlib.sha256((DEFAULT_MODEL_BUNDLE / "manifest.json").read_bytes()).hexdigest() == (
        provenance["bundle_manifest_sha256"]
    )
    summary = Analyzer.from_bundle().inspect()
    assert summary["bundle_id"] == provenance["bundle_id"]
    assert summary["profile_hmms"] == 361
    assert summary["specificity_overlay"]["search_spaces"] == {
        "core": 361, "inherited": 370, "new": 387,
    }
    fixture = ROOT / "tests/fixtures/canonical_type_iia_cassette.faa"
    with tempfile.TemporaryDirectory(prefix="casandra-installed-") as temporary:
        work = Path(temporary)
        assert run("inspect-model", cwd=work) == summary
        run("predict-proteins", "--input", fixture, "--output", work / "proteins.jsonl",
            "--threads", 1, cwd=work)
        records = [json.loads(line) for line in (work / "proteins.jsonl").read_text().splitlines()]
        assert len(records) == 4
        assert [r["result"] for r in records] == ["Cas9", "Cas1", "Cas2", "no cas"]
        for command in ("annotate-proteins", "classify-cassette"):
            run(command, "--input", fixture, "--output", work / command, "--threads", 1, cwd=work)

        # Synthetic reverse translation checks input plumbing, not biological accuracy.
        from Bio.Data import CodonTable

        codons = CodonTable.unambiguous_dna_by_id[11].back_table
        genome = ""
        gff = ["##gff-version 3"]
        for name, sequence in parse_fasta(fixture):
            start = len(genome) + 1
            coding = "".join(codons[amino] for amino in sequence) + "TAA"
            genome += coding
            gff.append(f"synthetic\tsmoke\tCDS\t{start}\t{len(genome)}\t.\t+\t0\tID={name}")
            genome += "ATGC" * 25
        (work / "genome.fna").write_text(">synthetic\n" + genome + "\n")
        (work / "genes.gff3").write_text("\n".join(gff) + "\n")
        run("predict-genome", "--genome", work / "genome.fna", "--gff3", work / "genes.gff3",
            "--output", work / "genome", "--threads", 1, cwd=work)
        run("predict-genome", "--genome", work / "genome.fna", "--gene-mode", "meta",
            "--output", work / "metagenome", "--threads", 1, cwd=work)
        for name in ("annotate-proteins", "classify-cassette", "genome", "metagenome"):
            output = work / name
            manifest = json.loads((output / "manifest.json").read_text())
            verify_artifacts(output, manifest, field="files")
            metadata = json.loads((output / "run.json").read_text())
            assert metadata["program_version"] == provenance["runtime_version"]
            assert metadata["bundle_id"] == provenance["bundle_id"]
    print(json.dumps({
        "status": "passed", "runtime": casandra.__version__, "bundle_id": summary["bundle_id"],
        "runtime_modules_verified": len(provenance["runtime_modules"]),
        "cli_checks": ["inspect-model", "predict-proteins", "annotate-proteins", "classify-cassette",
                       "predict-genome (supplied CDS)", "predict-genome (metagenomic CDS calling)"],
    }, indent=2))


if __name__ == "__main__":
    main()
