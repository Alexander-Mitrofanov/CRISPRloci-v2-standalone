"""Integration source/TSV gate; does not import the scientific tool."""

from pathlib import Path
import csv
import hashlib
import io
import math
import re

SCHEMA = {"name": "isescan.context", "version": "1.0.0"}
INTERPRETATION = "Predicted insertion sequences and terminal inverted repeats; overlap does not establish gene disruption or loss of function. Target-site duplications were not assessed."
FIELDS = [
    "seqID",
    "family",
    "cluster",
    "isBegin",
    "isEnd",
    "isLen",
    "ncopy4is",
    "start1",
    "end1",
    "start2",
    "end2",
    "score",
    "irId",
    "irLen",
    "nGaps",
    "orfBegin",
    "orfEnd",
    "strand",
    "orfLen",
    "E-value",
    "E-value4copy",
    "type",
    "ov",
    "tir",
]
DNA = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
MAX_BYTES = 20_000_000


def sha(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def read_regular(path, limit=MAX_BYTES):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size > limit:
        raise ValueError("Expected bounded regular file")
    with path.open("rb") as handle:
        payload = handle.read(limit + 1)
    if len(payload) > limit:
        raise ValueError("File exceeds limit")
    return payload


def fasta(payload):
    if len(payload) > 10_000_000:
        raise ValueError("FASTA exceeds byte limit")
    records = {}
    current = None
    for line in payload.decode("ascii").splitlines():
        if line.startswith(">"):
            words = line[1:].split()
            if not words or len(words[0]) > 200 or words[0] in records:
                raise ValueError("Invalid or duplicate sequence identifier")
            current = words[0]
            if not re.fullmatch(r"[A-Za-z0-9_.:|+\-]+", current):
                raise ValueError("Invalid sequence identifier")
            records[current] = []
        elif line.strip():
            if current is None:
                raise ValueError("Sequence before FASTA header")
            records[current].append(line.strip().upper())
    records = {key: "".join(parts) for key, parts in records.items()}
    if not 1 <= len(records) <= 200 or sum(map(len, records.values())) > 2_000_000:
        raise ValueError("Input count or total length exceeds limit")
    if any(not seq or set(seq) - DNA for seq in records.values()):
        raise ValueError("Invalid DNA sequence")
    return records


def integer(value, minimum=0):
    if not re.fullmatch(r"-?\d+", value):
        raise ValueError("Invalid integer")
    number = int(value)
    if number < minimum:
        raise ValueError("Integer out of range")
    return number


def metric(value, missing=False):
    number = float(value)
    if missing and number == -1:
        return None
    if not math.isfinite(number) or number < 0:
        raise ValueError("Invalid numeric evidence")
    return number


def interval(start, end, length):
    start, end = integer(start, 1), integer(end, 1)
    if not start <= end <= length:
        raise ValueError("Interval outside source")
    return start, end


def normalize(tsv, records, aliases):
    text = tsv.decode("ascii")
    table = csv.DictReader(io.StringIO(text), delimiter="\t")
    if table.fieldnames != FIELDS:
        raise ValueError("Unsupported native TSV columns")
    items = []
    seen = set()
    for ordinal, row in enumerate(table, 1):
        if ordinal > 10000 or None in row or any(v is None for v in row.values()):
            raise ValueError("Invalid native TSV row or count")
        ident = aliases[row["seqID"]]
        seq = records[ident]
        begin, end = interval(row["isBegin"], row["isEnd"], len(seq))
        if integer(row["isLen"], 1) != end - begin + 1:
            raise ValueError("IS length mismatch")
        native_type = row["type"]
        if native_type not in {"c", "p"}:
            raise ValueError("Unknown native IS type")
        family, cluster = row["family"], row["cluster"]
        if any(
            not value or len(value) > 200 or not re.fullmatch(r"[A-Za-z0-9_./|+\-]+", value)
            for value in (family, cluster)
        ):
            raise ValueError("Invalid model identifier")
        scope = (ident, begin, end, cluster)
        if scope in seen:
            raise ValueError("Duplicate native prediction")
        seen.add(scope)
        strand = row["strand"]
        orf = None
        if row["orfBegin"] == "0" and row["orfEnd"] == "0":
            if strand not in {"", "."} or row["orfLen"] not in {"0", "1"} or row["E-value4copy"] not in {"-1", "-1.0"}:
                raise ValueError("Inconsistent missing transposase evidence")
        else:
            left, right = interval(row["orfBegin"], row["orfEnd"], len(seq))
            if strand not in {"+", "-"} or integer(row["orfLen"], 1) != right - left + 1:
                raise ValueError("Inconsistent transposase interval")
            orf = {"start": left, "end": right, "strand": strand}
        tir = None
        coords = [row[k] for k in ("start1", "end1", "start2", "end2")]
        if coords == ["0"] * 4:
            if any(integer(row[k]) != 0 for k in ("score", "irId", "irLen", "nGaps")) or row["tir"] != "-:-":
                raise ValueError("Inconsistent missing inverted repeats")
        else:
            left1, right1 = interval(coords[0], coords[1], len(seq))
            left2, right2 = interval(coords[2], coords[3], len(seq))
            if not begin <= left1 <= right1 < left2 <= right2 <= end:
                raise ValueError("Inverted repeats outside IS")
            pair = row["tir"].upper().split(":")
            if len(pair) != 2 or len(pair[0]) != len(pair[1]) or len(pair[0]) > 10000:
                raise ValueError("Invalid inverted repeat alignment")
            a, b = pair
            if any(x == y == "-" for x, y in zip(a, b, strict=True)):
                raise ValueError("Double-gap inverted-repeat column")
            if (
                set(a + b) - DNA - {"-"}
                or a.replace("-", "") != seq[left1 - 1 : right1]
                or b.replace("-", "") != seq[left2 - 1 : right2].translate(COMPLEMENT)[::-1]
            ):
                raise ValueError("Inverted repeats differ from source DNA")
            matches = sum(x == y and x != "-" for x, y in zip(a, b, strict=True))
            gaps = sum(x == "-" or y == "-" for x, y in zip(a, b, strict=True))
            if integer(row["irLen"], 1) != len(a) or integer(row["irId"]) != matches or integer(row["nGaps"]) != gaps:
                raise ValueError("Inverted repeat alignment metrics mismatch")
            tir = {
                "left": {"start": left1, "end": right1},
                "right": {"start": left2, "end": right2},
                "left_alignment": a,
                "right_reverse_complement_alignment": b,
                "score": metric(row["score"]),
                "identical_positions": matches,
                "alignment_length": len(a),
                "gap_positions": gaps,
            }
        items.append(
            {
                "id": f"is{ordinal:06d}",
                "source_id": ident,
                "start": begin,
                "end": end,
                "length": end - begin + 1,
                "sequence_sha256": sha(seq[begin - 1 : end]),
                "family": family,
                "cluster": cluster,
                "native_type": native_type,
                "contig_edge": begin == 1 or end == len(seq),
                "copy_count": integer(row["ncopy4is"], 1),
                "transposase_orf": orf,
                "transposase_missing_reason": "not_identified_by_native_tool" if orf is None else None,
                "evalue_best_copy": metric(row["E-value"]),
                "evalue_this_copy": metric(row["E-value4copy"], missing=True),
                "terminal_inverted_repeats": tir,
                "target_site_duplication": None,
                "target_site_duplication_status": "not_assessed",
            }
        )
    return items


def document(payload, tsv, options, model_hashes):
    records = fasta(payload)
    aliases = {f"seq{i:06d}": ident for i, ident in enumerate(records, 1)}
    return {
        "schema": SCHEMA,
        "tool": {"name": "ISEScan-v2", "version": "2.0.0", "upstream_version": "1.7.3"},
        "coordinate_system": "1-based-inclusive-source-forward",
        "input_sha256": sha(payload),
        "sources": [{"id": key, "length": len(seq), "sha256": sha(seq)} for key, seq in records.items()],
        "options": options,
        "models": model_hashes,
        "elements": normalize(tsv, records, aliases),
        "interpretation": INTERPRETATION,
    }
