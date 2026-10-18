"""Shared, fail-closed wrappers for CRISPRidentify external programs."""

import math
import os
import re
import shutil
import subprocess
import sys


_CANDIDATE_IDENTIFIER = re.compile(r"^candidate_(\d+)$")


class ExternalToolError(RuntimeError):
    """A required external program failed or returned malformed output."""


def resolve_executable(name):
    environment_binary = os.path.join(
        os.path.dirname(os.path.abspath(sys.executable)), name
    )
    if (
            os.path.isfile(environment_binary)
            and os.access(environment_binary, os.X_OK)):
        return environment_binary
    executable = shutil.which(name)
    if executable is None:
        raise FileNotFoundError(
            "required executable {!r} was not found; run "
            "scripts/setup_local.sh and use scripts/run_crispridentify.sh"
            .format(name)
        )
    return os.path.abspath(executable)


def require_nonempty_file(path, description):
    absolute_path = os.path.abspath(path)
    if not os.path.isfile(absolute_path):
        raise FileNotFoundError(
            "required {} is missing: {}".format(description, absolute_path)
        )
    if os.path.getsize(absolute_path) == 0:
        raise ExternalToolError(
            "required {} is empty: {}".format(description, absolute_path)
        )
    return absolute_path


def run_external_tool(name, arguments, stdin=None, stdout=None):
    command = [resolve_executable(name)] + [
        os.fspath(argument) for argument in arguments
    ]
    try:
        completed = subprocess.run(
            command,
            stdin=stdin,
            stdout=stdout if stdout is not None else subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError as error:
        raise ExternalToolError(
            "could not execute {}: {}".format(name, error)
        ) from error
    if completed.returncode != 0:
        stderr = completed.stderr or b""
        if not isinstance(stderr, str):
            stderr = stderr.decode("utf-8", errors="replace")
        detail = stderr.strip()
        if len(detail) > 1000:
            detail = detail[-1000:]
        raise ExternalToolError(
            "{} failed with exit status {}{}".format(
                name,
                completed.returncode,
                ": {}".format(detail) if detail else "",
            )
        )
    return completed


def parse_prodigal_scores(output_path, sequence_lengths):
    """Parse Prodigal GenBank output into local one-based CDS intervals."""

    expected_identifiers = {
        "candidate_{}".format(index): index
        for index in range(len(sequence_lengths))
    }
    seen_identifiers = set()
    scores = [dict() for _ in sequence_lengths]
    current_identifier = None
    pending_interval = None

    with open(output_path) as output_file:
        for line_number, raw_line in enumerate(output_file, 1):
            line = raw_line.rstrip("\n")
            if line.startswith("DEFINITION"):
                if pending_interval is not None:
                    raise ExternalToolError(
                        "Prodigal CDS is missing its confidence qualifier "
                        "before line {}".format(line_number)
                    )
                header_match = re.search(r'seqhdr="([^"]+)"', line)
                if header_match is None:
                    raise ExternalToolError(
                        "Prodigal DEFINITION has no seqhdr at line {}"
                        .format(line_number)
                    )
                current_identifier = header_match.group(1)
                if current_identifier not in expected_identifiers:
                    raise ExternalToolError(
                        "Prodigal returned an unknown candidate: {}".format(
                            current_identifier
                        )
                    )
                if current_identifier in seen_identifiers:
                    raise ExternalToolError(
                        "Prodigal returned a duplicate candidate: {}".format(
                            current_identifier
                        )
                    )
                seen_identifiers.add(current_identifier)
                pending_interval = None
                continue

            if line.strip() == "//" and pending_interval is not None:
                raise ExternalToolError(
                    "Prodigal CDS is missing its confidence qualifier "
                    "before line {}".format(line_number)
                )

            if line.startswith("     CDS"):
                if current_identifier is None:
                    raise ExternalToolError(
                        "Prodigal CDS precedes its DEFINITION at line {}"
                        .format(line_number)
                    )
                if pending_interval is not None:
                    raise ExternalToolError(
                        "Prodigal CDS is missing its confidence qualifier "
                        "before line {}".format(line_number)
                    )
                location = line[len("     CDS"):].strip()
                location_match = re.fullmatch(
                    r"<?([0-9]+)\.\.>?([0-9]+)", location
                )
                is_complement = False
                if location_match is None:
                    location_match = re.fullmatch(
                        r"complement\(<?([0-9]+)\.\.>?([0-9]+)\)",
                        location,
                    )
                    is_complement = location_match is not None
                if location_match is None:
                    raise ExternalToolError(
                        "unsupported Prodigal CDS location at line {}: {}"
                        .format(line_number, location)
                    )
                start = int(location_match.group(1))
                end = int(location_match.group(2))
                candidate_index = expected_identifiers[current_identifier]
                if (
                        start < 1
                        or end < start
                        or end > sequence_lengths[candidate_index]):
                    raise ExternalToolError(
                        "Prodigal CDS lies outside candidate {}: {}..{} of {}"
                        .format(
                            current_identifier,
                            start,
                            end,
                            sequence_lengths[candidate_index],
                        )
                    )
                pending_interval = (start, end, is_complement)
                continue

            if pending_interval is not None and "/note=" in line:
                confidence_match = re.search(
                    r"(?:^|;)conf=([0-9]+(?:\.[0-9]+)?);", line
                )
                if confidence_match is None:
                    raise ExternalToolError(
                        "Prodigal CDS note has no confidence at line {}"
                        .format(line_number)
                    )
                candidate_index = expected_identifiers[current_identifier]
                if pending_interval in scores[candidate_index]:
                    raise ExternalToolError(
                        "Prodigal returned a duplicate CDS interval for {}: {}"
                        .format(current_identifier, pending_interval)
                    )
                confidence = float(confidence_match.group(1))
                if not math.isfinite(confidence):
                    raise ExternalToolError(
                        "Prodigal returned a non-finite confidence"
                    )
                scores[candidate_index][pending_interval] = confidence
                pending_interval = None

    if pending_interval is not None:
        raise ExternalToolError(
            "Prodigal output ended before the final CDS confidence qualifier"
        )
    missing = sorted(set(expected_identifiers) - seen_identifiers)
    if missing:
        raise ExternalToolError(
            "Prodigal did not return every candidate: {}".format(missing)
        )
    return scores
