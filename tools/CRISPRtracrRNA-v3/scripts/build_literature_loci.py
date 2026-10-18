"""Map published oriented reference RNAs exactly into frozen natural loci.

No fuzzy coordinate assignment: unmatched or multiply matched records stay out
of the localization denominator. Source downloads are indexed separately.
"""

import difflib
import hashlib
import json
from pathlib import Path

from crisprtracrrna.sequence import read_fasta, reverse_complement

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks/literature"


def main():
    rows = []
    for line in (DATA / "catalog.jsonl").open():
        reference = json.loads(line)
        if reference["source"] != "gasiunas2020":
            continue
        path = DATA / "raw/dooley_loci" / (reference["genome_label"] + ".fasta")
        contigs = list(read_fasta(path))
        matches = []
        sequence = reference["sequence"]
        for contig in contigs:
            for strand, query in (("+", sequence), ("-", reverse_complement(sequence))):
                offset = 0
                while (start := contig.bases.find(query, offset)) != -1:
                    matches.append(
                        {
                            "contig": contig.id,
                            "start": start + 1,
                            "end": start + len(query),
                            "strand": strand,
                        }
                    )
                    offset = start + 1
        best = None
        if not matches:
            best = max(
                difflib.SequenceMatcher(None, seq, contig.bases, autojunk=False)
                .find_longest_match()
                .size
                for contig in contigs
                for seq in (sequence, reverse_complement(sequence))
            )
        rows.append(
            {
                "id": reference["id"],
                "label": reference["label"],
                "source": reference["source"],
                "fasta": str(path.relative_to(ROOT)),
                "sequence": sequence,
                "repeat": reference["repeat"],
                "status": "exact_unique"
                if len(matches) == 1
                else "ambiguous"
                if matches
                else "unmatched",
                "matches": matches,
                "best_exact_block": best,
                "fasta_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "reference_sequence_sha256": reference["sequence_sha256"],
            }
        )
    (DATA / "loci.json").write_text(json.dumps(rows, indent=2) + "\n")
    mapped = sum(r["status"] == "exact_unique" for r in rows)
    print(f"{len(rows)} references; {mapped} uniquely mapped")


if __name__ == "__main__":
    main()
