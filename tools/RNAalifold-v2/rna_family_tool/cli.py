"""Offline bounded worker entry point; MAFFT is operator supplied."""
import argparse
import contextlib
import json
import os
import re
import signal
import subprocess
import sys
import threading
from pathlib import Path

from . import __version__
from .core import (
    LIMITS,
    POLICY,
    alignment_rows,
    columns,
    fasta,
    fold,
    pair_evidence,
    representatives,
    sha,
    structure_pairs,
)


def _stop_native(signum, frame):
    # Raising lets the command's finally block reap its isolated group.
    raise SystemExit(128 + signum)


def _run_native(command, *, timeout, check=True, capture_output=False, **kwargs):
    """Contain the native command and every descendant on all exit paths."""
    if capture_output:
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    child = subprocess.Popen(command, start_new_session=True, **kwargs)
    previous = None
    if threading.current_thread() is threading.main_thread():
        previous = signal.signal(signal.SIGTERM, _stop_native)
    try:
        stdout, stderr = child.communicate(timeout=timeout)
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(child.pid, signal.SIGKILL)
        child.wait()
        for stream in (child.stdin, child.stdout, child.stderr):
            if stream is not None:
                stream.close()
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)
    if check and child.returncode:
        raise subprocess.CalledProcessError(child.returncode, command, stdout, stderr)
    return subprocess.CompletedProcess(command, child.returncode, stdout, stderr)


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def align(groups, work, *, mafft, threads, timeout_seconds):
    """Run fixed FFT-NS-2 and validate source reconstruction before folding."""
    native_input = "".join(f'>{row["id"]}\n{row["sequence"].replace("U", "T")}\n' for row in groups)
    (work / "representatives.fasta").write_text(native_input, encoding="ascii")
    child_environment = os.environ.copy()
    child_environment.pop("MAFFT_BINARIES", None)
    version_run = _run_native([mafft, "--version"], capture_output=True, text=True, check=True, timeout=10, env=child_environment)
    match = re.search(r"v(7\.[0-9]+)\b", version_run.stdout + version_run.stderr)
    if not match:
        raise ValueError("Unsupported MAFFT version")
    with (work / "aligned.fasta").open("wb") as out, (work / "mafft.stderr").open("wb") as err:
        _run_native([mafft, "--nuc", "--inputorder", "--thread", str(threads), "--retree", "2", "--maxiterate", "0",
                        str(work / "representatives.fasta")], stdout=out, stderr=err, check=True, timeout=timeout_seconds, env=child_environment)
    aligned_path = work / "aligned.fasta"
    if aligned_path.stat().st_size > LIMITS["max_input_bytes"]:
        raise ValueError("Alignment exceeds byte limit")
    with aligned_path.open("rb") as stream:
        aligned = fasta(stream.read(LIMITS["max_input_bytes"] + 1), molecule="DNA", aligned=True)
    return alignment_rows(groups, aligned), match[1]


def run(input_path, output, *, molecule="RNA", mafft="mafft", threads=1, timeout_seconds=300):
    if type(threads) is not int or not 1 <= threads <= 8:
        raise ValueError("Threads must be 1–8")
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 14400:
        raise ValueError("Timeout must be 1–14400 seconds")
    source = Path(input_path)
    if source.is_symlink() or not source.is_file() or source.stat().st_size > LIMITS["max_input_bytes"]:
        raise ValueError("Invalid RNA family input file")
    with source.open("rb") as stream:
        payload = stream.read(LIMITS["max_input_bytes"] + 1)
    records = fasta(payload, molecule=molecule)
    sources, groups = representatives(records)
    destination = Path(output)
    if destination.exists() or destination.is_symlink():
        raise ValueError("Output directory must not already exist")
    destination.mkdir(parents=True)
    work = destination / "work"
    work.mkdir()
    alignment, mafft_version = align(groups, work, mafft=mafft, threads=threads, timeout_seconds=timeout_seconds)
    sequences = [row["sequence"] for row in alignment]
    native = fold(sequences)
    native["alignment_sha256"] = sha("".join(f'>{row["id"]}\n{row["sequence"]}\n' for row in alignment))
    write_json(work / "fold.json", native)
    consensus, column_evidence = columns(sequences)
    probability = {(row["left"], row["right"]): row["probability"] for row in native["probabilities"]}
    pairs = [pair_evidence(sequences, left, right, probability[left, right]) for left, right in structure_pairs(native["structure"])]
    report = {"schema": {"name": "rnaalifold.family", "version": "1.0.0"},
              "tool": {"name": "RNAalifold-v2", "version": __version__}, "status": "completed",
              "input": {"sha256": sha(payload), "molecule": molecule}, "policy": POLICY, "limits": LIMITS,
              "engine": {"name": "ViennaRNA", "version": "2.7.0"},
              "aligner": {"name": "MAFFT", "version": mafft_version, "strategy": "FFT-NS-2", "threads": threads},
              "sources": sources, "alignment": alignment, "columns": column_evidence,
              "fold": {"consensus": consensus, "structure": native["structure"], "temperature_celsius": 37,
                       "mfe_score_kcal_mol": native["mfe_score_kcal_mol"],
                       "ensemble_free_energy_kcal_mol": native["ensemble_free_energy_kcal_mol"]},
              "pairs": pairs, "probabilities": native["probabilities"],
              "counts": {"sources": len(sources), "unique_sequences": len(groups), "alignment_columns": len(consensus),
                         "mfe_pairs": len(pairs), "retained_probabilities": len(native["probabilities"])},
              "warnings": ["Homology, biological orientation and RNA boundaries are supplied by the user, not inferred.",
                           "Exact duplicates contribute one sequence weight; related nonidentical sequences are not phylogenetically reweighted.",
                           "Comparative folding scores include covariance contributions. Descriptive covariation is not a significance test or activity evidence."]}
    write_json(destination / "report.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--molecule", choices=("DNA", "RNA"), default="RNA")
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--threads", type=int, choices=range(1, 9), default=1)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    args = parser.parse_args(argv)
    if not 1 <= args.timeout_seconds <= 14400:
        parser.error("Timeout must be 1–14400 seconds")
    try:
        run(args.input, args.output, molecule=args.molecule, mafft=args.mafft, threads=args.threads,
            timeout_seconds=args.timeout_seconds)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"RNA family analysis failed: {type(exc).__name__}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
