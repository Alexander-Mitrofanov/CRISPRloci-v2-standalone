"""Optional source-forward FASTA views of accepted contextual Possible arrays."""
from pathlib import Path

from components.contextual_spacer_rescue import report_contextual_spacer_summary


FASTA_FILENAMES = (
    "Contextual_Arrays.fasta", "Contextual_Repeats.fasta", "Contextual_Spacers.fasta",
)


def stage_contextual_spacer_outputs(publication, result_path, report, *, enabled):
    """Write actual displayed members; never concatenate alternative endpoints."""
    selected = []
    for array in report["arrays"]:
        summary = report_contextual_spacer_summary(array)
        if summary is not None:
            selected.append((array, summary))
    for filename, kind in zip(FASTA_FILENAMES, ("array", "repeat", "spacer")):
        path = Path(result_path) / filename
        if not enabled or not selected:
            publication.delete(path)
            continue

        def write(stream, kind=kind):
            for array, summary in selected:
                suffix = (
                    " category=Possible confidence_ceiling=Possible"
                    " view=actual_checked_member strand=input_forward"
                    " member={display_geometry_id} admissible_geometries={member_count}"
                    " own_score={display_own_score} evidence=report.json"
                ).format(**summary)
                suffix += " boundary_uncertain={}".format(
                    "true" if summary["member_count"] > 1 else "false")
                parts = ([{"source_interval": array["source_interval"],
                           "sequence": array["source_sequence"], "index": 1}]
                         if kind == "array" else array[kind + "s"])
                for part in parts:
                    interval = part["source_interval"]
                    stream.write(
                        ">Contextual_{}_{}_{}_{}_{}{}\n{}\n".format(
                            array["id"], kind, part["index"], interval["start"],
                            interval["end"], suffix, part["sequence"]))

        publication.stage_writer(path, write, encoding="ascii")
