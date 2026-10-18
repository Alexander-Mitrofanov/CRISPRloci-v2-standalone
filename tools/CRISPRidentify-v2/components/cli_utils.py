import argparse
import os
from pathlib import Path
import shutil


def parse_boolean(value):
    """Parse an explicit, case-insensitive True/False CLI value."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == "true":
            return True
        if normalized == "false":
            return False
    raise argparse.ArgumentTypeError(
        "expected a boolean value: True or False"
    )


def positive_integer(value):
    """Parse a strictly positive integer before pipeline work begins."""

    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise argparse.ArgumentTypeError(
            "expected a positive integer"
        )
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            "expected a positive integer"
        )
    if parsed < 1:
        raise argparse.ArgumentTypeError(
            "expected a positive integer"
        )
    return parsed


def readable_input_file(value):
    """Return an absolute readable regular-file path for ``argparse``."""

    path = Path(value).expanduser().absolute()
    if not path.exists():
        raise argparse.ArgumentTypeError(
            "input file does not exist: {}".format(path)
        )
    if not path.is_file():
        raise argparse.ArgumentTypeError(
            "input file is not a regular file: {}".format(path)
        )
    if not os.access(str(path), os.R_OK):
        raise argparse.ArgumentTypeError(
            "input file is not readable: {}".format(path)
        )
    return path


def readable_input_directory(value):
    """Return an absolute readable/searchable directory for ``argparse``."""

    path = Path(value).expanduser().absolute()
    if not path.exists():
        raise argparse.ArgumentTypeError(
            "input folder does not exist: {}".format(path)
        )
    if not path.is_dir():
        raise argparse.ArgumentTypeError(
            "input folder is not a directory: {}".format(path)
        )
    if not os.access(str(path), os.R_OK | os.X_OK):
        raise argparse.ArgumentTypeError(
            "input folder is not readable/searchable: {}".format(path)
        )
    return path


def _executable_path(name, project_candidates, executable_finder):
    for candidate in project_candidates:
        candidate = Path(candidate)
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return executable_finder(name)


def optional_analysis_errors(
        project_root,
        is_element=False,
        cas=False,
        executable_finder=None):
    """Return actionable dependency failures for requested optional analyses."""

    project_root = Path(project_root)
    if executable_finder is None:
        executable_finder = shutil.which

    errors = []
    if is_element:
        fasta36 = _executable_path(
            "fasta36",
            (project_root / "tools" / "fasta" / "fasta36",),
            executable_finder,
        )
        if fasta36 is None:
            errors.append(
                "--is_element True requires fasta36, but it was not found; "
                "run scripts/setup_local.sh --profile is-elements and launch "
                "with CRISPRIDENTIFY_ENV=.envs/crispridentify-is-elements"
            )

    if cas:
        from components.casandra_adapter import casandra_installation_error

        error = casandra_installation_error()
        if error is not None:
            errors.append(error)

    return errors
