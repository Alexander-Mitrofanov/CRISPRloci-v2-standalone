"""Recover the 41 natural tracrRNA predictions from Briner Figure S2B.

The publisher PDF has an incorrect text map: its displayed alignment hyphen
decodes as capital N. The nucleotide glyphs A/C/G/U decode correctly. Page 3
was visually checked before applying this source-specific conversion. Do not
reuse it on arbitrary RNA data or interpret N as an unknown nucleotide here.

The last twenty rows have a separately positioned continuation with hidden
repeated species labels. Only sequence characters after x=466 in that second
block belong in the visible alignment. Keep both decoded and raw aligned forms
for audit. No sequence is inferred from another source.
"""

import csv
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path

import pdfplumber
from pypdf import PdfReader, PdfWriter
from pypdf.generic import IndirectObject

ROOT = Path(__file__).resolve().parents[1] / "benchmarks/literature/raw/alternate_retrieval"
PDF = ROOT / "briner_mmc1.pdf"
EXPECTED_PDF_SHA256 = "ddac48d05e4806ce9a140ac9d42d9b60d32cc32fef237d8fce2a78b4dd17c836"
LABELS = [
    "S. pyogenes",
    "S. dysagalactiae",
    "S. equi",
    "S. thermophilus CR3",
    "S. salivarius",
    "S. gallolyticus",
    "S. lutetiensis",
    "S. anginosis B",
    "S. mitis",
    "S. sanguinis",
    "S. oralis",
    "S. mutans",
    "S. intermedius",
    "S. anginosus A",
    "S. thermophilus CR1",
    "S. vestibularis",
    "S. gordonii",
    "S. parasanguinis",
    "S. orisratti",
    "S. henryi",
    "S. infantarius",
    "L. brevis",
    "L. buchneri A",
    "L. mucosae",
    "L. curvatus",
    "L. fermentum",
    "L. reuteri",
    "L. plantarum",
    "L. pentosus",
    "L. coryniformis",
    "L. buchneri B",
    "L. otakiensis",
    "L. paracasei",
    "L. casei",
    "L. rhamnosus",
    "L. gasseri",
    "L. hominis",
    "L. johnsonii",
    "L. jensenii",
    "L. delbrueckii",
    "L. rossiae",
]


def extract_alignments(chars):
    start = next(
        i for i, c in enumerate(chars) if c["text"] == "S" and c["x0"] < 15 and 441 < c["top"] < 443
    )
    continuation = next(
        i
        for i, c in enumerate(chars)
        if c["text"] == "L" and 403 < c["x0"] < 404 and 604 < c["top"] < 606
    )
    primary = chars[start:continuation]
    ys = [
        c["top"]
        for c in primary
        if c["text"] in {"S", "L"} and c["x0"] < 15 and 440 < c["top"] < 710
    ]
    assert len(ys) == 41, len(ys)
    alignments = []
    for label, y in zip(LABELS, ys, strict=True):
        first = [c for c in primary if abs(c["top"] - y) < 0.05 and c["x0"] >= 80]
        second = (
            [
                c
                for c in chars[continuation:]
                if abs(c["top"] - y) < 0.3 and c["x0"] >= 466 and c["top"] < 710
            ]
            if label.startswith("L.")
            else []
        )
        raw = "".join(c["text"] for c in sorted(first, key=lambda c: c["x0"]))
        raw += "".join(c["text"] for c in sorted(second, key=lambda c: c["x0"]))
        alignments.append((label, y, raw))
    return alignments


def main():
    assert hashlib.sha256(PDF.read_bytes()).hexdigest() == EXPECTED_PDF_SHA256
    with pdfplumber.open(PDF) as document:
        original = extract_alignments(document.pages[2].chars)
    # Independently verify gap decoding using the font's actual WinAnsi bytes.
    # The erroneous Unicode map converts byte 0x2d (hyphen) to U+004e (N).
    # Remove only that font's map from an in-memory page clone, then compare
    # every aligned nucleotide and gap against the positional extraction.
    reader = PdfReader(PDF)
    font = IndirectObject(940, 0, reader).get_object()
    assert font["/BaseFont"] == "/RKACQV+ArialMT"
    assert font["/Encoding"] == "/WinAnsiEncoding"
    unicode_map = font["/ToUnicode"].get_object().get_data().decode()
    assert "<2d> <004e>" in unicode_map
    del font["/ToUnicode"]
    writer = PdfWriter()
    writer.add_page(reader.pages[2])
    memory = io.BytesIO()
    writer.write(memory)
    memory.seek(0)
    with pdfplumber.open(memory) as verified:
        corrected = extract_alignments(verified.pages[0].chars)
    rows = []
    for i, ((label, y, raw), (label2, _, aligned)) in enumerate(
        zip(original, corrected, strict=True), 1
    ):
        assert label == label2
        assert set(raw) <= set("ACGUN"), (label, raw)
        assert set(aligned) <= set("ACGU-"), (label, aligned)
        assert raw.replace("N", "-") == aligned, label
        rna = aligned.replace("-", "")
        sequence = rna.replace("U", "T")
        rows.append(
            {
                "label": label,
                "source_locator": f"Supplemental Figure S2B, PDF page 3, RNA row {i}",
                "raw_pdf_text_alignment": raw,
                "aligned_rna": aligned,
                "rna_sequence": rna,
                "sequence": sequence,
                "length": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "subtype": "II-A",
                "evidence": "computational_prediction_in_primary_functional_study",
                "training_overlap": "source_collection_used_to_build_bundled_Type_II_models",
                "source_doi": "10.1016/j.molcel.2014.09.019",
                "source_url": "https://ars.els-cdn.com/content/image/1-s2.0-S1097276514007515-mmc1.pdf",
                "normalization": (
                    "Visually verified PDF gap glyph decodes as N; map to alignment gap, "
                    "remove gaps, convert U to T"
                ),
                "figure_y": y,
            }
        )
    (ROOT / "briner2014_figure_sequences.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (ROOT / "briner2014_figure_sequences.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (ROOT / "briner2014_figure_sequences.fasta").open("w") as f:
        for i, row in enumerate(rows, 1):
            f.write(f">Briner2014_{i:02d} {row['label']}\n{row['sequence']}\n")
    (ROOT / "briner2014_extraction_qc.json").write_text(
        json.dumps(
            {
                "source_pdf_sha256": EXPECTED_PDF_SHA256,
                "source_url": rows[0]["source_url"],
                "extractor": "scripts/extract_briner_figure.py",
                "dependencies": {
                    name: importlib.metadata.version(name) for name in ("pdfplumber", "pypdf")
                },
                "visual_review": "Rendered PDF page 3, Figure S2B; 41 natural RNA alignment rows",
                "font_validation": (
                    "WinAnsi byte 0x2d is hyphen; erroneous ToUnicode maps it to N. "
                    "In-memory removal of that font map independently reconstructs "
                    "identical aligned RNA/gap strings for all 41 rows."
                ),
                "sequence_characters": "ACGU only after removing printed alignment gaps",
                "ambiguous_residues": 0,
                "preserved": (
                    "Raw PDF-decoded aligned strings and corrected aligned RNA retained per row; "
                    "source PDF is unchanged"
                ),
                "boundaries": (
                    "Published predicted native RNA regions; "
                    "do not infer experimental boundaries for all 41 systems"
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "records": len(rows),
                "unique": len({r["sequence"] for r in rows}),
                "length_range": [min(r["length"] for r in rows), max(r["length"] for r in rows)],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
