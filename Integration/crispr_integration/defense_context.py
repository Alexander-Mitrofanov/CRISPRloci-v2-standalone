"""Offline DefenseFinder execution and source-bound output gate."""
import hashlib
import math
from pathlib import Path
import re

from .protospacer import _records, _payload, _json

MODEL = {"id": "defense-finder-models-1.2.3", "version": "1.2.3",
         "commit": "5bdb41257fc7234fb8c3b70e3d362c52541860dd", "license": "GPL-3.0",
         "scope": "DefenseFinder_1-5_and_RM;CasFinder_and_newer_models_not_assessed",
         "manifest_sha256": "bb4dc6de8dfce9179e5d14b1649d1af1e27ca5b0e14840787e883aa46a8d7348"}
POLICY = {"coordinates": "one_based_inclusive", "topology": "linear",
          "gene_caller": "pyrodigal-3.7.1-meta", "coverage_profile": 0.4,
          "scope": "submitted_sequences_only", "functional_activity": "not_assessed"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
CODONS = dict(zip((a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                 "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", strict=True))


def sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def validate_report(doc, records, input_sha256):
    if (doc["schema"] != {"name": "defensefinder.context", "version": "1.0.0"}
            or doc["tool"] not in ({"name": "DefenseFinder-v2", "version": "2.0.0"}, {"name": "DefenseFinder-v2", "version": "2.0.1"})
            or doc["execution_completed"] is not True or doc["model"] != MODEL or doc["policy"] != POLICY
            or doc["input"] != {"sha256": input_sha256, "record_count": len(records)}
            or doc["sources"] != [{"id": k, "length": len(v), "sha256": sha(v)} for k, v in records.items()]):
        raise ValueError("DefenseFinder identity or source mismatch")
    genes, system_ids, positions = {}, set(), {}
    if len(doc["genes"]) > 10000 or len(doc["hits"]) > 100000 or len(doc["components"]) > 100000:
        raise ValueError("DefenseFinder evidence count exceeded")
    source_order = {key: n for n, key in enumerate(records, 1)}
    for gene in doc["genes"]:
        source, start, end = gene["source_id"], gene["start"], gene["end"]
        if (source not in records or type(start) is not int or type(end) is not int
                or not 1 <= start <= end <= len(records[source]) or type(gene["strand"]) is not int
                or gene["strand"] not in (-1, 1) or gene["translation_table"] not in (4, 11)
                or any(type(gene[key]) is not bool for key in ("partial_begin", "partial_end", "initiator_methionine"))):
            raise ValueError("DefenseFinder CDS bounds/strand mismatch")
        position, previous_start = positions.get(source, (0, 0))
        position += 1
        replicon = f"r{source_order[source]:06d}"
        if (gene["id"] != f"{replicon}_{position:06d}" or gene["replicon"] != replicon
                or type(gene["position"]) is not int or gene["position"] != position
                or start < previous_start or gene["id"] in genes):
            raise ValueError("DefenseFinder CDS order mismatch")
        positions[source] = position, start
        cds = records[source][start-1:end]
        if gene["strand"] == -1:
            cds = cds.translate(COMPLEMENT)[::-1]
        partial_5prime = gene["partial_begin"] if gene["strand"] == 1 else gene["partial_end"]
        partial_3prime = gene["partial_end"] if gene["strand"] == 1 else gene["partial_begin"]
        stops = {"TAA", "TAG"} if gene["translation_table"] == 4 else {"TAA", "TAG", "TGA"}
        if (len(cds) % 3 or gene["initiator_methionine"] != (not partial_5prime)
                or not partial_3prime and cds[-3:] not in stops):
            raise ValueError("DefenseFinder CDS frame or complete boundary mismatch")
        protein = "".join("W" if gene["translation_table"] == 4 and cds[i:i+3] == "TGA"
                          else CODONS.get(cds[i:i+3], "X") for i in range(0, len(cds)-2, 3)).removesuffix("*")
        if gene["initiator_methionine"]:
            if cds[:3] not in ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"):
                raise ValueError("DefenseFinder invalid translation start")
            protein = "M" + protein[1:]
        if (gene["cds_sha256"] != sha(cds) or gene["protein"] != protein
                or gene["protein_sha256"] != sha(protein) or not protein or "*" in protein):
            raise ValueError("DefenseFinder translation mismatch")
        genes[gene["id"]] = gene
    for system in doc["systems"]:
        ids = system["gene_ids"]
        if (system["id"] in system_ids or not ids or len(set(ids)) != len(ids)
                or any(key not in genes or genes[key]["source_id"] != system["source_id"] for key in ids)
                or ids != sorted(ids, key=lambda key: genes[key]["position"])
                or system["start"] != min(genes[key]["start"] for key in ids)
                or system["end"] != max(genes[key]["end"] for key in ids)
                or not re.fullmatch(r"[A-Za-z0-9_.:/-]{1,300}", system["model"])
                or system["type"] != system["model"].split("/")[-2]
                or system["subtype"] != system["model"].split("/")[-1]):
            raise ValueError("DefenseFinder system mapping mismatch")
        system_ids.add(system["id"])
        for key in ("wholeness", "score"):
            if type(system[key]) not in (int, float) or not math.isfinite(system[key]):
                raise ValueError("Invalid system score")
        if not 0 <= system["wholeness"] <= 1:
            raise ValueError("Invalid system completeness")
    membership = {(s["id"], key) for s in doc["systems"] for key in s["gene_ids"]}
    if membership != {(c["system_id"], c["gene_id"]) for c in doc["components"]}:
        raise ValueError("System component membership mismatch")
    for row in doc["components"] + doc["hits"]:
        if row["gene_id"] not in genes or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,200}", row["profile"]):
            raise ValueError("Unbound profile hit")
        if any(type(row[key]) not in (int, float) or not math.isfinite(row[key]) for key in ("score", "i_evalue", "profile_coverage")):
            raise ValueError("Invalid profile evidence")
        if row["i_evalue"] < 0 or not 0 <= row["profile_coverage"] <= 1:
            raise ValueError("Invalid profile range")
    for hit in doc["hits"]:
        if (type(hit["protein_start"]) is not int or type(hit["protein_end"]) is not int
                or not 1 <= hit["protein_start"] <= hit["protein_end"] <= len(genes[hit["gene_id"]]["protein"])
                or type(hit["sequence_coverage"]) not in (int, float) or not 0 <= hit["sequence_coverage"] <= 1):
            raise ValueError("Invalid HMM alignment interval")
    counts = {key: len(doc[key]) for key in ("sources", "genes", "systems", "components", "hits")}
    if doc["counts"] != counts or any(type(v) is not int for v in doc["counts"].values()):
        raise ValueError("DefenseFinder count mismatch")


def run_defense_context_pipeline(input_path, output_path, *, runner_config, cpu=1, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    records = _records(Path(input_path).absolute())
    if not 1 <= len(records) <= 20 or sum(map(len, records.values())) > 2000000:
        raise ValueError("Supply at most 20 genomic sequences and 2 Mb")
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "defense-input").mkdir()
    fasta = job / "defense-input/genomes.fasta"
    fasta.write_text("".join(f">{key}\n{seq}\n" for key, seq in records.items()), encoding="ascii")
    result = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
              "workflow": "defense_context", "status": "running", "stages": {}, "tool_options": runner_config.tool_options}
    _atomic_json(job / "integration_result.json", result)
    try:
        from .context_annotations import run_cas_context
        result["stages"]["casandra"] = run_cas_context(job, records, runner_config, timeout_seconds)
        stage = _run_process("defensefinder", [*runner_config.defensefinder,
            "--input", str(fasta), "--output", str(job / "defense-context"), "--workers", str(min(cpu, 8))],
            job / "logs", timeout_seconds=timeout_seconds)
        doc = _json(_payload(job, job / "defense-context/report.json"))
        validate_report(doc, records, sha(fasta.read_bytes()))
        stage.update(counts=doc["counts"], model_id=MODEL["id"], model_sha256=MODEL["manifest_sha256"])
        result["stages"]["defensefinder"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("defensefinder", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
