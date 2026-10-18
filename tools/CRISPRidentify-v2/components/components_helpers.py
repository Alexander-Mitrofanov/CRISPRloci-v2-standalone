"""Strict, non-mutating FASTA input preparation helpers."""

import os
import re
import tempfile
from os import listdir
from os.path import isfile, join


DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")


def process_string_from_header(input_string):
    """Retain the historical NCBI version-stripping filename convention."""

    def replace_match(match):
        if match.group(2).isdigit():
            return match.group(1)
        return match.group(1) + "-" + match.group(2)

    return re.sub(r"(\w+)\.(\w+)", replace_match, input_string)


def _parse_fasta_records(file_path):
    """Return validated ``(header_without_marker, sequence)`` records."""

    records = []
    header = None
    sequence_parts = []

    def finish_record(line_number):
        if header is None:
            return
        if not sequence_parts:
            raise ValueError(
                "{} has an empty FASTA record {!r} before line {}".format(
                    file_path, header, line_number
                )
            )
        records.append((header, "".join(sequence_parts)))

    with open(file_path, encoding="utf-8") as fasta_file:
        for line_number, raw_line in enumerate(fasta_file, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record(line_number)
                header = line[1:].strip()
                sequence_parts = []
                if not header:
                    raise ValueError(
                        "{} line {} has an empty FASTA header".format(
                            file_path, line_number
                        )
                    )
                continue
            if header is None:
                raise ValueError(
                    "{} line {} contains sequence before a FASTA header".format(
                        file_path, line_number
                    )
                )
            invalid = sorted(set(line).difference(DNA_ALPHABET))
            if invalid:
                raise ValueError(
                    "{} line {} contains non-IUPAC DNA symbols: {}".format(
                        file_path,
                        line_number,
                        ", ".join(repr(symbol) for symbol in invalid),
                    )
                )
            sequence_parts.append(line)

    finish_record(line_number + 1 if "line_number" in locals() else 1)
    if not records:
        raise ValueError("{} contains no FASTA records".format(file_path))
    return records


def _safe_record_name(header, flag_ncbi_formatting):
    name_source = header.split()[0] if flag_ncbi_formatting else header
    if flag_ncbi_formatting:
        name_source = process_string_from_header(name_source)
    safe_name = re.sub(r"[^A-Za-z0-9_-]+", "_", name_source).strip("_")
    if not safe_name:
        raise ValueError("FASTA header {!r} has no safe identifier".format(header))
    return safe_name


def _write_records(records, destination, flag_ncbi_formatting):
    seen_headers = set()
    seen_names = set()
    for header, sequence in records:
        identifier = header.split()[0]
        if identifier in seen_headers:
            raise ValueError(
                "duplicate FASTA record identifier: {!r}".format(identifier)
            )
        seen_headers.add(identifier)

        safe_name = _safe_record_name(header, flag_ncbi_formatting)
        if safe_name in seen_names:
            raise ValueError(
                "FASTA record names collide after sanitization: {!r}".format(
                    safe_name
                )
            )
        seen_names.add(safe_name)

        # Filename normalization is a presentation concern.  The staged
        # record remains the same biological source, so preserve its complete
        # header (including an accession version) for source-bound provenance.
        output_header = header
        output_path = os.path.join(destination, safe_name + ".fa")
        with open(output_path, "x", encoding="utf-8") as output_file:
            output_file.write(">{}\n{}\n".format(output_header, sequence))


def _temporary_record_directory(label):
    safe_label = re.sub(r"[^A-Za-z0-9_-]+", "_", label).strip("_")
    return tempfile.mkdtemp(prefix=(safe_label or "fasta") + "_")


def multiline_fasta_handle_python(file, flag_ncbi_formatting=False):
    records = _parse_fasta_records(file)
    base_name = os.path.basename(file).split(".")[0]
    destination = _temporary_record_directory(base_name)
    try:
        _write_records(records, destination, flag_ncbi_formatting)
    except Exception:
        # The directory contains only newly generated staging files.  Keep a
        # failed split from leaking partial inputs into a later run.
        import shutil
        shutil.rmtree(destination, ignore_errors=True)
        raise
    return destination


def folder_of_multifasta_handle(folder_multifasta):
    """Split all FASTA records in a folder without concatenating raw files."""

    input_files = sorted(
        join(folder_multifasta, name)
        for name in listdir(folder_multifasta)
        if isfile(join(folder_multifasta, name))
    )
    if not input_files:
        raise ValueError(
            "{} contains no input files".format(folder_multifasta)
        )
    records = []
    for input_file in input_files:
        records.extend(_parse_fasta_records(input_file))

    destination = _temporary_record_directory("multifasta_folder")
    try:
        _write_records(records, destination, flag_ncbi_formatting=False)
    except Exception:
        import shutil
        shutil.rmtree(destination, ignore_errors=True)
        raise
    return destination
