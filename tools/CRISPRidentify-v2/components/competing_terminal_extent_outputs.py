"""Presentation of unaccepted terminal observations, never accepted arrays."""

import csv
import html
import json
from pathlib import Path
import shutil


CSV_FILENAME = "Competing_terminal_extents.csv"
CSV_FIELDS = (
    "Accepted array ID", "Accepted repeat count", "Accepted own score",
    "Observation ID", "Geometry ID", "Candidate ID", "Observed category",
    "Observed repeat count", "Observed own score", "Disposition",
    "Block reasons", "Recorded negative assessment", "Repeat intervals (0-based half-open)",
    "Report evidence",
)


def _json(value):
    return json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":"))


def _score(value):
    return format(value, ".6g")


def _escape(value):
    return html.escape(str(value), quote=True)


def render_competing_extent_details(observations):
    """Render already validated summaries under their accepted array card."""
    if not observations:
        return ""
    sections = []
    for item in observations:
        intervals = "; ".join("[{}, {})".format(*interval)
                              for interval in item["source_intervals"])
        paths = []
        for index, (interval, repeat) in enumerate(zip(
                item["source_intervals"], item["repeats"])):
            paths.append("Repeat {} [{}, {}): {}".format(index + 1, *interval, repeat))
            if index < len(item["spacers"]):
                paths.append("Spacer {}: {}".format(index + 1, item["spacers"][index]))
        reasons = _json(item["block_reasons"]) if item["block_reasons"] else "None recorded"
        sections.append("""\
              <section class="competing-terminal-extent">
                <h4>Unaccepted competing terminal extent</h4>
                <p>Accepted call: {accepted_count} repeats, Possible; own score {accepted_score}.<br>
                   Observed path: {count} repeats, {category}; own score {score}.</p>
                <p>Disposition: <strong>{status}</strong>. Block reasons: {reasons}.</p>
                <p>Recorded negative assessment: <code>{negatives}</code>.
                   <code>not_recorded</code> means unknown, not a passing check.
                   Scores are retained from existing final rows; they were not recomputed for this view.</p>
                <p>This observation is excluded from accepted counts.
                   Aliases are observations, not independent support votes.</p>
                <details class="nested-details">
                  <summary>Complete observed source path ({count} repeats)</summary>
                  <p>Coordinates: 0-based, half-open, input forward strand.<br>{intervals}</p>
                  <pre class="sequence-block">{path}</pre>
                  <p>Observation: <code>{observation_id}</code><br>
                     Geometry: <code>{geometry_id}</code><br>
                     Candidate: <code>{candidate_id}</code><br>
                     Full evidence: <code>{pointer}</code></p>
                </details>
              </section>""".format(
            accepted_count=item["accepted_repeat_count"],
            accepted_score=_escape(_score(item["accepted_own_score"])),
            count=item["repeat_count"], category=_escape(item["category"]),
            score=_escape(_score(item["own_score"])), status=_escape(item["status"]),
            reasons=_escape(reasons), negatives=_escape(_json(item["negative_review"])),
            intervals=_escape(intervals), path=_escape("\n".join(paths)),
            observation_id=_escape(item["observation_id"]),
            geometry_id=_escape(item["geometry_id"]),
            candidate_id=_escape(item["candidate_id"]),
            pointer=_escape(item["report_pointer"])))
    return "\n".join(sections)


def _text_notes(summaries):
    for array_id, observations in summaries.items():
        for item in observations:
            yield (
                "\nUnaccepted competing terminal extent for {}\n"
                "Accepted call: {} repeats, Possible; own score {}.\n"
                "Observed path: {} repeats, {}; own score {}.\n"
                "Disposition: {}. Block reasons: {}.\n"
                "Recorded negative assessment: {}. not_recorded means unknown.\n"
                "Stored scores were not recomputed. This observation does not change accepted counts.\n"
                "Full source, alignment and provenance: {}\n"
            ).format(array_id, item["accepted_repeat_count"], item["accepted_own_score"],
                     item["repeat_count"], item["category"], item["own_score"],
                     item["status"], _json(item["block_reasons"]),
                     _json(item["negative_review"]), item["report_pointer"])


def stage_competing_extent_outputs(publication, result_path, summaries, *, machine_profile=False):
    """Stage the companion and text notes after the current legacy writer.

    The writer replaces Possible_Candidates.txt each run. Copy that new text
    in bounded blocks, then append the disclosure within the report transaction.
    Machine mode's normal cleanup owns deletion of the legacy text file.
    """
    root = Path(result_path)
    csv_path = root / CSV_FILENAME
    if machine_profile or not any(summaries.values()):
        publication.delete(csv_path)
        return

    def write_csv(stream):
        writer = csv.writer(stream)
        writer.writerow(CSV_FIELDS)
        for array_id, observations in summaries.items():
            for item in observations:
                writer.writerow((
                    array_id, item["accepted_repeat_count"], item["accepted_own_score"],
                    item["observation_id"], item["geometry_id"], item["candidate_id"],
                    item["category"], item["repeat_count"], item["own_score"],
                    item["status"], _json(item["block_reasons"]),
                    _json(item["negative_review"]), _json(item["source_intervals"]),
                    item["report_pointer"],
                ))

    publication.stage_writer(csv_path, write_csv, newline="")
    text_path = root / "Possible_Candidates.txt"

    def write_text(stream):
        with text_path.open(encoding="utf-8") as original:
            shutil.copyfileobj(original, stream, length=64 * 1024)
        stream.writelines(_text_notes(summaries))

    publication.stage_writer(text_path, write_text)
