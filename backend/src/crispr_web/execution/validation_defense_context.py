"""Independent source and native-table validation for defense context.

No imports from Integration or scientific tool packages.
"""
import csv
import hashlib
import io
import itertools
import math
import re

from .validation_annotations import _payload, _strict_json, _fasta

MODEL = {"id": "defense-finder-models-1.2.3", "version": "1.2.3",
         "commit": "5bdb41257fc7234fb8c3b70e3d362c52541860dd", "license": "GPL-3.0",
         "scope": "DefenseFinder_1-5_and_RM;CasFinder_and_newer_models_not_assessed",
         "manifest_sha256": "bb4dc6de8dfce9179e5d14b1649d1af1e27ca5b0e14840787e883aa46a8d7348"}
POLICY = {"coordinates": "one_based_inclusive", "topology": "linear", "gene_caller": "pyrodigal-3.7.1-meta",
          "coverage_profile": 0.4, "scope": "submitted_sequences_only", "functional_activity": "not_assessed"}
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASE = "TCAG"
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_TOKEN = re.compile(r"[A-Za-z0-9_.:-]{1,200}\Z")


def _sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def _finite(value, lower=None, upper=None):
    return (type(value) in (int, float) and math.isfinite(value)
            and (lower is None or value >= lower) and (upper is None or value <= upper))


def _translation(cds, code, start):
    amino = []
    for i in range(0, len(cds)-2, 3):
        codon = cds[i:i+3]
        if code == 4 and codon == "TGA":
            amino.append("W")
        elif all(base in _BASE for base in codon):
            amino.append(_AA[_BASE.index(codon[0])*16 + _BASE.index(codon[1])*4 + _BASE.index(codon[2])])
        else:
            amino.append("X")
    if amino and amino[-1] == "*":
        amino.pop()
    if start:
        if cds[:3] not in ("ATG", "TTG", "GTG", "CTG", "ATA", "ATC", "ATT") or not amino:
            raise ValueError("Invalid CDS start")
        amino[0] = "M"
    return "".join(amino)


def validate_defense_context_pipeline(claimed, root, expected, integration, input_root, checkpoint=None):
    try:
        return _validate(claimed, root, expected, integration, input_root, checkpoint)
    except (OSError, ValueError, TypeError, KeyError, IndexError, AttributeError, UnicodeError, RecursionError):
        return False


def _validate(claimed, root, expected, integration, input_root, checkpoint):
    if input_root is None or claimed.molecule != "DNA" or integration.get("workflow") != "defense_context":
        return False
    records = {}
    for path in sorted(input_root.iterdir()):
        key, seq = _fasta(_payload(input_root, path, checkpoint))
        if key in records or expected.get(key) != (len(seq), _sha(seq)):
            return False
        records[key] = seq
    if set(records) != set(expected) or not 1 <= len(records) <= 20 or sum(map(len, records.values())) > 2000000:
        return False
    doc = _strict_json(_payload(root, root / "defense-context/report.json", checkpoint))
    if (doc["schema"] != {"name": "defensefinder.context", "version": "1.0.0"}
            or doc["tool"] not in ({"name": "DefenseFinder-v2", "version": "2.0.0"}, {"name": "DefenseFinder-v2", "version": "2.0.1"})
            or doc["execution_completed"] is not True or doc["model"] != MODEL or doc["policy"] != POLICY
            or doc["input"] != {"sha256": _sha("".join(f">{k}\n{s}\n" for k, s in records.items())), "record_count": len(records)}
            or doc["sources"] != [{"id": k, "length": len(s), "sha256": _sha(s)} for k, s in records.items()]):
        return False
    for key, limit in (("genes", 10000), ("systems", 20000), ("components", 100000), ("hits", 100000)):
        if not isinstance(doc[key], list) or len(doc[key]) > limit:
            return False
    genes, previous = {}, {}
    source_order = {key: n for n, key in enumerate(records, 1)}
    for row in doc["genes"]:
        if checkpoint:
            checkpoint()
        key = row["source_id"]
        left, right = row["start"], row["end"]
        if (key not in records or type(left) is not int or type(right) is not int
                or not 1 <= left <= right <= len(records[key]) or type(row["strand"]) is not int
                or row["strand"] not in (-1, 1) or type(row["translation_table"]) is not int
                or row["translation_table"] not in (4, 11)
                or any(type(row[p]) is not bool for p in ("partial_begin", "partial_end", "initiator_methionine"))):
            return False
        previous_position, previous_start = previous.get(key, (0, 0))
        position = previous_position + 1
        replicon = f"r{source_order[key]:06d}"
        if (row["id"] != f"{replicon}_{position:06d}" or row["replicon"] != replicon
                or type(row["position"]) is not int or row["position"] != position
                or left < previous_start or row["id"] in genes):
            return False
        previous[key] = position, left
        cds = records[key][left-1:right]
        if row["strand"] == -1:
            cds = cds.translate(_COMPLEMENT)[::-1]
        # Pyrodigal flags identify genomic left/right edges, regardless of strand.
        incomplete_start = row["partial_begin"] if row["strand"] == 1 else row["partial_end"]
        incomplete_stop = row["partial_end"] if row["strand"] == 1 else row["partial_begin"]
        stop_codons = {"TAA", "TAG"} if row["translation_table"] == 4 else {"TAA", "TAG", "TGA"}
        if (len(cds) % 3 or row["initiator_methionine"] != (not incomplete_start)
                or not incomplete_stop and cds[-3:] not in stop_codons):
            return False
        protein = _translation(cds, row["translation_table"], row["initiator_methionine"])
        if (not protein or "*" in protein or row["protein"] != protein
                or row["cds_sha256"] != _sha(cds) or row["protein_sha256"] != _sha(protein)):
            return False
        genes[row["id"]] = row
    systems = {}
    for row in doc["systems"]:
        members = row["gene_ids"]
        if (not _TOKEN.fullmatch(row["id"]) or row["id"] in systems or not members
                or len(set(members)) != len(members) or any(g not in genes for g in members)
                or any(genes[g]["source_id"] != row["source_id"] for g in members)
                or members != sorted(members, key=lambda g: genes[g]["position"])
                or type(row["start"]) is not int or type(row["end"]) is not int
                or row["start"] != min(genes[g]["start"] for g in members)
                or row["end"] != max(genes[g]["end"] for g in members)
                or not _finite(row["wholeness"], 0, 1) or not _finite(row["score"])
                or not _TOKEN.fullmatch(row["type"]) or not _TOKEN.fullmatch(row["subtype"])
                or not re.fullmatch(r"defense-finder-models/(DefenseFinder_[1-5]|RM)/[A-Za-z0-9_./-]{1,200}", row["model"])
                or row["type"] != row["model"].split("/")[-2] or row["subtype"] != row["model"].split("/")[-1]):
            return False
        systems[row["id"]] = row
    pairs = {(s["id"], g) for s in systems.values() for g in s["gene_ids"]}
    if pairs != {(c["system_id"], c["gene_id"]) for c in doc["components"]}:
        return False
    for row in doc["components"] + doc["hits"]:
        if (row["gene_id"] not in genes or not _TOKEN.fullmatch(row["profile"])
                or not _finite(row["score"]) or not _finite(row["i_evalue"], 0)
                or not _finite(row["profile_coverage"], 0, 1)):
            return False
    for row in doc["hits"]:
        if (type(row["protein_start"]) is not int or type(row["protein_end"]) is not int
                or not 1 <= row["protein_start"] <= row["protein_end"] <= len(genes[row["gene_id"]]["protein"])
                or not _finite(row["sequence_coverage"], 0, 1)):
            return False
    if not _native_matches(root, doc, genes, systems, checkpoint):
        return False
    from .validation_context_annotations import validate_cas_context
    validate_cas_context(root, records, integration["stages"]["casandra"], checkpoint)
    counts = {key: len(doc[key]) for key in ("sources", "genes", "systems", "components", "hits")}
    stage = integration["stages"]["defensefinder"]
    return (all(type(value) is int for value in doc["counts"].values()) and doc["counts"] == counts
            and stage["counts"] == counts and stage["status"] == "completed"
            and stage["model_id"] == MODEL["id"] and stage["model_sha256"] == MODEL["manifest_sha256"])


def _native_matches(root, doc, genes, systems, checkpoint):
    if not genes:
        return not doc["systems"] and not doc["components"] and not doc["hits"]
    native_proteins = _payload(root, root / "defense-context/work/proteins.faa", checkpoint).decode("ascii")
    expected_proteins = "".join(f'>{gene["id"]}\n{gene["protein"]}\n' for gene in doc["genes"])
    if native_proteins != expected_proteins:
        return False
    ranks = {gene["id"]: (rank, len(gene["protein"])) for rank, gene in enumerate(doc["genes"], 1)}
    components, hits = [], {}
    native_bytes = 0
    found_systems = set()
    for family in [f"DefenseFinder_{n}" for n in range(1, 6)] + ["RM"]:
        folder = root / "defense-context/work" / family
        payload = _payload(root, folder / "best_solution.tsv", checkpoint).decode("utf-8")
        native_bytes += len(payload.encode("utf-8"))
        for row in csv.DictReader((line for line in io.StringIO(payload) if not line.startswith("#")), delimiter="\t"):
            system = systems[row["sys_id"]]
            if (row["model_fqn"] != system["model"] or row["replicon"] != genes[row["hit_id"]]["replicon"]
                    or (int(row["hit_pos"]), int(row["hit_seq_len"])) != ranks[row["hit_id"]]
                    or float(row["sys_wholeness"]) != system["wholeness"] or float(row["sys_score"]) != system["score"]):
                return False
            found_systems.add(row["sys_id"])
            components.append({"system_id": row["sys_id"], "gene_id": row["hit_id"], "profile": row["gene_name"],
                               "status": row["hit_status"], "score": float(row["hit_score"]),
                               "i_evalue": float(row["hit_i_eval"]), "profile_coverage": float(row["hit_profile_cov"])})
        paths = sorted(itertools.islice((folder / "hmmer_results").glob("*.res_hmm_extract"), 1001))
        if len(paths) > 1000:
            return False
        for path in paths:
            payload = _payload(root, path, checkpoint).decode("utf-8")
            native_bytes += len(payload.encode("utf-8"))
            if native_bytes > 20_000_000:
                return False
            for row in csv.reader((line for line in io.StringIO(payload) if not line.startswith("#")), delimiter="\t"):
                if not row:
                    continue
                if (len(row) != 11 or row[0] not in genes or row[1] != genes[row[0]]["replicon"]
                        or (int(row[2]), int(row[3])) != ranks[row[0]]):
                    return False
                hit = {"gene_id": row[0], "profile": row[4], "i_evalue": float(row[5]), "score": float(row[6]),
                       "profile_coverage": float(row[7]), "sequence_coverage": float(row[8]),
                       "protein_start": int(row[9]), "protein_end": int(row[10])}
                key = row[0], row[4]
                if key not in hits or hit["score"] > hits[key]["score"]:
                    hits[key] = hit
    components.sort(key=lambda row: (row["system_id"], row["gene_id"], row["profile"]))
    return (found_systems == set(systems) and components == doc["components"]
            and [hits[key] for key in sorted(hits)] == doc["hits"])
