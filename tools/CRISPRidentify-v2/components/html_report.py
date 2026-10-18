"""Dependency-free HTML output with exact source-reconstruction checks."""

import hashlib
import html
import json
import math
import os

from components.biological_loci import EXTENSION_SCHEMA_NAME
from components.cas_html import CAS_STYLES
from components.cas_html import render_cas_section
from components.biological_loci import EXTENSION_SCHEMA_VERSION
from components.feature_schema import feature_display_name
from components.report_model import CATEGORY_KEYS
from components.report_model import CATEGORY_NAMES
from components.report_model import ReportModelBuilder
from components.report_model import SCHEMA_NAME
from components.output_transaction import OutputTransaction
from components.output_transaction import remove_completion_marker
from components.report_extensions import CLASSIFICATION_SCHEMA_NAME
from components.report_extensions import CLASSIFICATION_SCHEMA_VERSION
from components.report_extensions import LIKELY_TANDEM_SCHEMA_NAME
from components.report_extensions import LIKELY_TANDEM_SCHEMA_VERSION
from components.report_extensions import SOURCE_GEOMETRY_SCHEMA_NAME
from components.report_extensions import SOURCE_GEOMETRY_SCHEMA_VERSION
from components.report_extensions import _likely_tandem_core_architecture_matches
from components.report_extensions import attach_report_extensions
from components.reconstruction_validator import read_single_fasta
from components.reconstruction_boundary_family import report_boundary_family_summary
from components.contextual_spacer_rescue import report_contextual_spacer_summary
from components.possible_reconstruction import (
    EVIDENCE_FILENAME as POSSIBLE_RECONSTRUCTION_EVIDENCE_FILENAME,
    report_possible_reconstruction_summary,
)
from components.competing_terminal_extents import competing_terminal_extent_summaries
from components.competing_terminal_extent_outputs import (
    CSV_FILENAME as COMPETING_EXTENT_CSV_FILENAME,
    render_competing_extent_details,
)
from components.likely_tandem import LIKELY_TANDEM_POLICY
from components.likely_tandem import evaluate_likely_tandem
from components.likely_tandem import RULE_NAME as LIKELY_TANDEM_RULE_NAME
from components.likely_tandem import RULE_VERSION as LIKELY_TANDEM_RULE_VERSION
from components.likely_tandem_transaction import (
    CLASSIFICATION_BASIS_METHOD,
    TANDEM_TAG,
    TRANSACTION_RULE_NAME,
    TRANSACTION_RULE_VERSION,
    resolve_likely_tandem_audit_inputs,
)


def canonical_report_json(model):
    """Serialize a report model as strict, deterministic JSON."""

    return json.dumps(
        model,
        allow_nan=False,
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
    ) + "\n"


def _write_canonical_report_json(stream, model):
    """Write byte-identical canonical JSON with bounded encoder buffering."""
    encoder = json.JSONEncoder(
        allow_nan=False, ensure_ascii=False, indent=2, sort_keys=True,
    )
    chunks = []
    size = 0
    for chunk in encoder.iterencode(model):
        chunks.append(chunk)
        size += len(chunk)
        if size >= 64 * 1024:
            stream.write("".join(chunks))
            chunks.clear()
            size = 0
    if chunks:
        stream.write("".join(chunks))
    stream.write("\n")


def _publish_report_model(model, result_path, write_html, publication=None):
    """Publish the detached model built by either validated writer path.

    Kept private: callers must finish ReportModelBuilder and extension/audit
    validation before running legacy writers or entering this publication step.
    The standalone CanonicalReportMaker still performs its own validation.
    """
    json_path = os.path.join(result_path, "report.json")
    html_path = os.path.join(result_path, "report.html")
    remove_completion_marker(json_path)
    html_content = render_html_report(model) if write_html else None

    def stage_and_publish(transaction):
        companion = os.path.abspath(os.path.join(result_path, COMPETING_EXTENT_CSV_FILENAME))
        if companion not in {os.fspath(path) for path in transaction.targets}:
            # Standalone JSON/HTML writers do not generate the CSV companion.
            # Retire a previous one instead of pairing it with a fresh report.
            transaction.delete(companion)
        possible_audit = os.path.abspath(os.path.join(
            result_path, POSSIBLE_RECONSTRUCTION_EVIDENCE_FILENAME,
        ))
        if possible_audit not in {os.fspath(path) for path in transaction.targets}:
            # A standalone writer has per-array proofs but no run-wide audit.
            transaction.delete(possible_audit)
        if html_content is None:
            transaction.delete(html_path)
        else:
            transaction.stage_text(html_path, html_content)
        transaction.stage_writer(
            json_path, lambda stream: _write_canonical_report_json(stream, model)
        )
        transaction.publish(completion_target=json_path)

    if publication is None:
        with OutputTransaction() as transaction:
            stage_and_publish(transaction)
    else:
        stage_and_publish(publication)


def _script_safe_json(model):
    """Serialize JSON without allowing an embedded ``</script>`` terminator."""

    serialized = json.dumps(
        model,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    return (
        serialized
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def _escape(value):
    return html.escape(str(value), quote=True)


def _format_score(score):
    return "Not scored" if score is None else "{:.4f}".format(score)


def _array_tags(array):
    tags = array.get("tags", [])
    if (
        not isinstance(tags, list)
        or any(
            not isinstance(tag, str)
            or not tag
            or len(tag) > 64
            or tag[0] not in "abcdefghijklmnopqrstuvwxyz"
            or any(
                character not in "abcdefghijklmnopqrstuvwxyz0123456789_-"
                for character in tag
            )
            for tag in tags
        )
        or tags != sorted(set(tags))
    ):
        raise ValueError(
            "array tags must be unique sorted lowercase identifiers"
        )
    return tuple(tags)


def _tag_badges(array):
    tags = _array_tags(array)
    badges = [
        '<span class="tag-badge">{}</span>'.format(_escape(tag))
        for tag in tags
    ]
    if report_boundary_family_summary(array) is not None:
        badges.append('<span class="tag-badge">Uncertain repeat ends</span>')
    contextual = report_contextual_spacer_summary(array)
    if contextual is not None:
        badges.append('<span class="tag-badge">Localized identical-spacer block</span>')
        if contextual["member_count"] > 1:
            badges.append('<span class="tag-badge">Uncertain repeat ends</span>')
    if not badges:
        return '<span class="muted">None</span>'
    return " ".join(badges)


def _boundary_family_details(array):
    family = report_boundary_family_summary(array)
    if family is None:
        return ""
    members = "".join('<li><code>{}</code></li>'.format(_escape(member))
                      for member in family["admissible_geometry_ids"])
    return """\
              <section class="interpretation-card">
                <h3>Shared extent, uncertain repeat ends</h3>
                <p>{count:,} supported repeat positions span {start:,}&ndash;{end:,}.
                   {members_count} complete boundary alternatives remain admissible.
                   The sequences, coordinates and score shown here belong to one
                   actual member; they do not identify a uniquely resolved segmentation.</p>
                <p>Original accepted core: {core_count} repeats.
                   Member own scores: {minimum_score}&ndash;{maximum_score}.</p>
                <details class="nested-details">
                  <summary>Displayed member and complete alternatives</summary>
                  <p>Family: <code>{family_id}</code><br>
                     Displayed member: <code>{display_id}</code></p>
                  <ul>{members}</ul>
                  <p>Complete coupled geometries and their evidence are in report.json.
                     Observed repeat-end choices cannot be mixed into new paths.</p>
                </details>
              </section>""".format(
        count=family["supported_repeat_count"], start=family["union_span"][0] + 1,
        end=family["union_span"][1], members_count=family["member_count"],
        core_count=family["core_repeat_count"],
        minimum_score=_escape(_format_score(family["member_score_min"])),
        maximum_score=_escape(_format_score(family["member_score_max"])),
        family_id=_escape(family["family_id"]),
        display_id=_escape(family["display_geometry_id"]), members=members,
    )


def _orientation_text(orientation):
    if orientation["status"] == "not_computed":
        return "Not computed (?)"
    if orientation["strand"] == "+":
        return "Forward (+)"
    return "Reverse (-)"


def _possible_reconstruction_details(array):
    summary = report_possible_reconstruction_summary(array)
    if summary is None:
        return ""
    return """\
              <section class="interpretation-card">
                <h3>Checked extent of a Possible core</h3>
                <p>The original {core_count}-repeat Possible core supports this
                   checked {count}-repeat path. The result remains Possible.</p>
                <p>Original core score: {core_score}; freshly checked core score:
                   {fresh_core_score}; full path's own score: {own_score}.
                   These scores are not calibrated probabilities.</p>
                <p>Source, structure and applicable negative checks passed.
                   Their scope does not establish biological authenticity.
                   The complete evidence is recorded in report.json.</p>
              </section>""".format(
        core_count=summary["core_repeat_count"], count=summary["repeat_count"],
        core_score=_escape(_format_score(summary["core_original_own_score"])),
        fresh_core_score=_escape(_format_score(summary["core_fresh_own_score"])),
        own_score=_escape(_format_score(summary["own_score"])),
    )


def _contextual_spacer_details(array):
    summary = report_contextual_spacer_summary(array)
    if summary is None:
        return ""
    return """\
              <section class="interpretation-card">
                <h3>Localized identical-spacer block</h3>
                <p>The complete {count}-repeat path and both maximal flanks
                   passed the contextual checks. The result remains Possible.</p>
                <p>The flank checks are correlated diagnostic context, not
                   independent votes. The displayed path's own score is {score}.
                   Scores are not calibrated probabilities.</p>
                <p>{members} complete boundary interpretations are retained.
                   Coordinates and sequences show actual member <code>{display}</code>.
                   Alternative boundaries cannot be mixed into a new path.</p>
                <p>The source, structure and applicable negative checks support
                   this limited interpretation; they do not establish biological
                   authenticity. Complete evidence is in report.json.</p>
              </section>""".format(
        count=summary["repeat_count"], score=_escape(_format_score(summary["display_own_score"])),
        members=summary["member_count"], display=_escape(summary["display_geometry_id"]),
    )


def _wrap_sequence(sequence, width=100):
    return "\n".join(
        sequence[index:index + width]
        for index in range(0, len(sequence), width)
    )


def _summary_table(model):
    rows = []
    for array in model["arrays"]:
        interval = array["source_interval"]
        rows.append(
            """\
          <tr data-category="{category_key}">
            <th scope="row"><a href="#{array_id}">{array_id}</a></th>
            <td><span class="badge badge-{category_key}">{category}</span></td>
            <td>{tags}</td>
            <td class="numeric">{start:,}&ndash;{end:,}</td>
            <td>{orientation}</td>
            <td class="numeric">{repeat_count}</td>
            <td class="numeric">{spacer_count}</td>
            <td><code>{consensus}</code></td>
            <td class="numeric">{score}</td>
            <td><span class="validation validation-valid">Passed</span></td>
          </tr>""".format(
                category_key=_escape(array["category_key"]),
                array_id=_escape(array["id"]),
                category=_escape(array["category"]),
                tags=_tag_badges(array),
                start=interval["start"],
                end=interval["end"],
                orientation=_escape(_orientation_text(array["orientation"])),
                repeat_count=array["repeat_count"],
                spacer_count=array["spacer_count"],
                consensus=_escape(array["consensus_repeat"]),
                score=_escape(_format_score(array["certainty_score"])),
            )
        )

    empty_row = ""
    if not rows:
        empty_row = """\
          <tr class="empty-table-row">
            <td colspan="10">No CRISPR arrays were reported for this sequence.</td>
          </tr>"""

    return """\
      <div class="table-scroll" tabindex="0" aria-label="Scrollable array summary">
        <table id="array-summary">
          <caption>Reported CRISPR arrays in source-coordinate order</caption>
          <thead>
            <tr>
              <th scope="col">Array ID</th>
              <th scope="col">Category</th>
              <th scope="col">Tags</th>
              <th scope="col">Source coordinates</th>
              <th scope="col">Orientation</th>
              <th scope="col">Repeats</th>
              <th scope="col">Spacers</th>
              <th scope="col">Consensus repeat</th>
              <th scope="col">Confidence score</th>
              <th scope="col">Exact source reconstruction</th>
            </tr>
          </thead>
          <tbody>
{rows}
{empty_row}
          </tbody>
        </table>
      </div>""".format(
        rows="\n".join(rows),
        empty_row=empty_row,
    )


def _filters(model):
    counts = model["summary"]["category_counts"]
    show_tagged_low_score = any(
        array.get("category") == "Low score"
        and "tandem" in _array_tags(array)
        for array in model.get("arrays", ())
    )
    checkboxes = []
    for name, key in zip(CATEGORY_NAMES, CATEGORY_KEYS):
        checked = (
            ""
            if key == "low-score" and not show_tagged_low_score
            else " checked"
        )
        checkboxes.append(
            """\
          <label class="filter-option">
            <input type="checkbox" data-category-filter="{key}"{checked}>
            <span class="badge badge-{key}">{name}</span>
            <span aria-label="{count} arrays">({count})</span>
          </label>""".format(
                key=_escape(key),
                checked=checked,
                name=_escape(name),
                count=counts[name],
            )
        )
    return """\
      <fieldset class="filters">
        <legend>Filter arrays by category</legend>
{checkboxes}
        <span class="filter-actions">
          <button type="button" data-filter-action="all">Show all</button>
          <button type="button" data-filter-action="none">Hide all</button>
        </span>
      </fieldset>
      <p id="visible-status" class="status-line" aria-live="polite"></p>""".format(
        checkboxes="\n".join(checkboxes)
    )


def _architecture(array):
    segments = []
    for element in array["architecture"]:
        if element["type"] == "spacer-deletion":
            boundary = element["source_boundary"]
            label = (
                "Complete spacer deletion {}: 0-based interbase boundary {}, "
                "between repeats {} and {}"
            ).format(
                element["index"],
                boundary["position"],
                element["flanking_repeat_indices"][0],
                element["flanking_repeat_indices"][1],
            )
            segments.append(
                """\
              <span class="segment segment-spacer-deletion"
                    title="{label}"
                    aria-hidden="true"></span>""".format(
                    label=_escape(label),
                )
            )
            continue

        interval = element["source_interval"]
        label = "{} {}: {}-{}, {} bases, {}".format(
            element["type"].capitalize(),
            element["index"],
            interval["start"],
            interval["end"],
            interval["length"],
            element["sequence"],
        )
        segments.append(
            """\
              <span class="segment segment-{element_type}"
                    style="flex-grow:{length}"
                    title="{label}"
                    aria-hidden="true"></span>""".format(
                element_type=_escape(element["type"]),
                length=interval["length"],
                label=_escape(label),
            )
        )

    component_count = len(array["architecture"])
    return """\
          <div class="architecture-legend" aria-hidden="true">
            <span><i class="legend-repeat"></i> R: repeat</span>
            <span><i class="legend-spacer"></i> S: spacer</span>
            <span><i class="legend-deletion"></i> D: complete spacer deletion</span>
          </div>
          <div class="architecture-viewport" tabindex="0"
               role="img"
               aria-label="{aria_label}">
            <div class="architecture" style="--component-count:{component_count}">
{segments}
            </div>
          </div>
          <p class="scroll-hint">
            The architecture is proportional to source length. Large arrays
            remain horizontally scrollable.
          </p>""".format(
        aria_label=_escape(
            "{} has {} repeats, {} non-deleted spacers, and {} complete "
            "spacer deletions from {} to {}".format(
                array["id"],
                array["repeat_count"],
                array.get("non_deleted_spacer_count", array["spacer_count"]),
                array.get("spacer_deletion_count", 0),
                array["source_interval"]["start"],
                array["source_interval"]["end"],
            )
        ),
        component_count=component_count,
        segments="\n".join(segments),
    )


def _component_table(array):
    rows = []
    for element in array["architecture"]:
        if element["type"] == "spacer-deletion":
            boundary = element["source_boundary"]
            rows.append(
                """\
                  <tr class="spacer-deletion-row">
                    <th scope="row">Spacer deletion {index}</th>
                    <td class="numeric">interbase {position:,} (0-based)</td>
                    <td class="numeric">0</td>
                    <td><span class="deletion-event-label">&empty; Complete spacer deletion</span></td>
                  </tr>""".format(
                    index=element["index"],
                    position=boundary["position"],
                )
            )
            continue

        interval = element["source_interval"]
        rows.append(
            """\
                  <tr>
                    <th scope="row">{element_type} {index}</th>
                    <td class="numeric">{start:,}&ndash;{end:,}</td>
                    <td class="numeric">{length}</td>
                    <td><code class="sequence-inline">{sequence}</code></td>
                  </tr>""".format(
                element_type=_escape(element["type"].capitalize()),
                index=element["index"],
                start=interval["start"],
                end=interval["end"],
                length=interval["length"],
                sequence=_escape(element["sequence"]),
            )
        )
    return """\
          <details class="nested-details">
            <summary>Repeat and spacer source intervals ({count} components)</summary>
            <div class="component-table-wrap" tabindex="0"
                 aria-label="Scrollable repeat and spacer table">
              <table class="component-table">
                <caption>Independently source-sliced array components</caption>
                <thead>
                  <tr>
                    <th scope="col">Component</th>
                    <th scope="col">Source coordinates</th>
                    <th scope="col">Length</th>
                    <th scope="col">Source DNA</th>
                  </tr>
                </thead>
                <tbody>
{rows}
                </tbody>
              </table>
            </div>
          </details>""".format(
        count=len(array["architecture"]),
        rows="\n".join(rows),
    )


def _feature_tables(array):
    if not array["feature_sets"]:
        return """\
          <p class="muted">Evaluated feature vectors were not available.</p>"""

    tables = []
    for feature_set in array["feature_sets"]:
        rows = []
        for feature in feature_set["features"]:
            value = json.dumps(
                feature["value"],
                allow_nan=False,
                ensure_ascii=False,
                sort_keys=True,
            )
            rows.append(
                """\
                    <tr>
                      <th scope="row" title="Internal feature key: {key}">{name}</th>
                      <td><code>{value}</code></td>
                    </tr>""".format(
                    key=_escape(feature["name"]),
                    name=_escape(feature_display_name(feature["name"])),
                    value=_escape(value),
                )
            )
        tables.append(
            """\
              <table class="feature-table">
                <caption>Input value table</caption>
                <thead>
                  <tr>
                    <th scope="col">Feature</th>
                    <th scope="col">Value</th>
                  </tr>
                </thead>
                <tbody>
{rows}
                </tbody>
              </table>""".format(
                rows="\n".join(rows),
            )
        )
    return """\
          <details class="nested-details">
            <summary>Classifier input values</summary>
            <div class="feature-grid">
{tables}
            </div>
          </details>""".format(tables="\n".join(tables))


_EXTENSION_CONTRACTS = {
    "classification_semantics": (
        CLASSIFICATION_SCHEMA_NAME,
        CLASSIFICATION_SCHEMA_VERSION,
        "declarative",
    ),
    "source_geometry": (
        SOURCE_GEOMETRY_SCHEMA_NAME,
        SOURCE_GEOMETRY_SCHEMA_VERSION,
        "computed",
    ),
    "biological_loci": (
        EXTENSION_SCHEMA_NAME,
        EXTENSION_SCHEMA_VERSION,
        "computed",
    ),
    "likely_tandem": (
        LIKELY_TANDEM_SCHEMA_NAME,
        LIKELY_TANDEM_SCHEMA_VERSION,
        "applied",
    ),
}

_LIKELY_TANDEM_STATS_FIELDS = frozenset((
    "total_rows_before",
    "total_rows_after",
    "accepted_rows_reviewed",
    "accepted_rows_before",
    "accepted_rows_after",
    "likely_tandem_count",
    "retained_count",
    "bona_fide_quarantined",
    "possible_quarantined",
    "diagnostic_rows_before",
    "diagnostic_rows_after",
))


def _is_plain_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _strict_json_equal(left, right):
    """Compare JSON values without Python's bool/int numeric coercion."""

    try:
        options = {
            "allow_nan": False,
            "ensure_ascii": True,
            "separators": (",", ":"),
            "sort_keys": True,
        }
        return json.dumps(left, **options) == json.dumps(right, **options)
    except (TypeError, ValueError):
        return False


def _validate_source_geometry_extension(model, extension):
    source_sha256 = model.get("source", {}).get("sha256")
    if extension.get("source_sha256") != source_sha256:
        raise ValueError(
            "source_geometry extension does not match report source"
        )
    arrays = extension.get("arrays")
    scope = extension.get("scope")
    if not isinstance(arrays, list) or not isinstance(scope, dict):
        raise ValueError(
            "source_geometry extension arrays/scope are malformed"
        )
    accepted_ids = {
        array["id"]
        for array in model.get("arrays", [])
        if array.get("category") in ("Bona-fide", "Possible")
    }
    geometry_ids = []
    for position, geometry in enumerate(arrays):
        if not isinstance(geometry, dict):
            raise ValueError(
                "source_geometry array {} must be a mapping".format(
                    position
                )
            )
        array_id = geometry.get("array_id")
        if not isinstance(array_id, str) or not array_id:
            raise ValueError(
                "source_geometry array {} has no id".format(position)
            )
        if geometry.get("source_sha256") != source_sha256:
            raise ValueError(
                "source_geometry array source hash mismatch"
            )
        statistics = geometry.get("statistics")
        advisories = geometry.get("advisories")
        if not isinstance(statistics, dict) or not isinstance(
            advisories, list
        ):
            raise ValueError(
                "source_geometry array statistics/advisories are malformed"
            )
        for advisory in advisories:
            if (
                not isinstance(advisory, dict)
                or not isinstance(advisory.get("code"), str)
                or not isinstance(advisory.get("message"), str)
                or advisory.get("decision_effect") != "none"
            ):
                raise ValueError(
                    "source_geometry advisory violates its non-decision "
                    "contract"
                )
        geometry_ids.append(array_id)
    if (
        len(geometry_ids) != len(set(geometry_ids))
        or set(geometry_ids) != accepted_ids
        or scope.get("array_count") != len(arrays)
    ):
        raise ValueError(
            "source_geometry extension coverage disagrees with accepted arrays"
        )


def _validate_biological_loci_extension(model, extension):
    source = model.get("source", {})
    source_sha256 = source.get("sha256")
    source_length = source.get("length")
    if extension.get("source_sha256") != source_sha256:
        raise ValueError(
            "biological_loci extension does not match report source"
        )
    counts = extension.get("counts")
    loci = extension.get("loci")
    if not isinstance(counts, dict) or not isinstance(loci, list):
        raise ValueError(
            "biological_loci extension counts/loci are malformed"
        )

    accepted_ids = {
        array["id"]
        for array in model.get("arrays", [])
        if array.get("category") in ("Bona-fide", "Possible")
    }
    seen_members = []
    singleton_count = 0
    for position, locus in enumerate(loci):
        if not isinstance(locus, dict):
            raise ValueError(
                "biological_loci locus {} must be a mapping".format(
                    position
                )
            )
        locus_id = locus.get("id")
        interval = locus.get("source_interval")
        member_ids = locus.get("member_ids")
        if (
            not isinstance(locus_id, str)
            or not locus_id
            or not isinstance(interval, dict)
            or not isinstance(member_ids, list)
            or not member_ids
            or not all(
                isinstance(member_id, str) and member_id
                for member_id in member_ids
            )
            or len(member_ids) != len(set(member_ids))
            or locus.get("member_count") != len(member_ids)
        ):
            raise ValueError(
                "biological_loci locus {} is malformed".format(position)
            )
        start = interval.get("start")
        end = interval.get("end")
        length = interval.get("length")
        if (
            not _is_plain_int(start)
            or not _is_plain_int(end)
            or not _is_plain_int(length)
            or start < 1
            or end < start
            or length != end - start + 1
            or not _is_plain_int(source_length)
            or end > source_length
        ):
            raise ValueError(
                "biological_loci locus interval is malformed"
            )
        seen_members.extend(member_ids)
        singleton_count += len(member_ids) == 1

    if (
        len(seen_members) != len(set(seen_members))
        or set(seen_members) != accepted_ids
    ):
        raise ValueError(
            "biological_loci membership disagrees with accepted arrays"
        )
    expected_counts = {
        "raw_accepted_array_count": len(accepted_ids),
        "biological_locus_count": len(loci),
        "singleton_locus_count": singleton_count,
        "multi_member_locus_count": len(loci) - singleton_count,
    }
    if counts != expected_counts:
        raise ValueError(
            "biological_loci counts disagree with memberships"
        )


def _validate_likely_tandem_ledger_architecture(record, source_length):
    interval = record.get("source_interval")
    architecture = record.get("source_architecture")
    if not isinstance(interval, dict) or not isinstance(architecture, dict):
        raise ValueError("likely_tandem ledger geometry is malformed")
    start = interval.get("start")
    end = interval.get("end")
    length = interval.get("length")
    if (
            not _is_plain_int(start)
            or not _is_plain_int(end)
            or not _is_plain_int(length)
            or start < 1
            or end < start
            or length != end - start + 1
            or not _is_plain_int(source_length)
            or end > source_length
            or architecture.get("coordinate_system")
            != "zero_based_half_open"
            or architecture.get("reference") != "input_forward_strand"):
        raise ValueError("likely_tandem ledger interval is malformed")
    repeats = architecture.get("repeats")
    spacers = architecture.get("spacers")
    spacer_copy = record.get("spacers_source_forward")
    if (
            not isinstance(repeats, list)
            or not isinstance(spacers, list)
            or not isinstance(spacer_copy, list)
            or len(repeats) != len(spacers) + 1
            or len(spacers) != len(spacer_copy)
            or record.get("repeat_count") != len(repeats)):
        raise ValueError("likely_tandem ledger cardinality is malformed")

    def validate_elements(elements, expected, label):
        result = []
        for position, element in enumerate(elements, 1):
            source_interval = (
                element.get("source_interval")
                if isinstance(element, dict) else None
            )
            sequence = element.get("sequence") if isinstance(
                element, dict
            ) else None
            if not isinstance(source_interval, dict):
                raise ValueError(
                    "likely_tandem ledger {} geometry is malformed".format(
                        label
                    )
                )
            element_start = source_interval.get("start")
            element_end = source_interval.get("end")
            if (
                    not _is_plain_int(element.get("index"))
                    or element.get("index") != position
                    or not isinstance(sequence, str)
                    or not sequence
                    or not _is_plain_int(element_start)
                    or not _is_plain_int(element_end)
                    or element_start < start - 1
                    or element_end <= element_start
                    or element_end > end
                    or element_end - element_start != len(sequence)
                    or element.get("sequence_sha256")
                    != hashlib.sha256(sequence.encode("utf-8")).hexdigest()
                    or (
                        expected is not None
                        and sequence != expected[position - 1]
                    )):
                raise ValueError(
                    "likely_tandem ledger {} geometry is malformed".format(
                        label
                    )
                )
            result.append((element_start, element_end, sequence))
        return result

    repeat_values = validate_elements(repeats, None, "repeat")
    spacer_values = validate_elements(spacers, spacer_copy, "spacer")
    reconstructed = []
    for position, repeat in enumerate(repeat_values):
        reconstructed.append(repeat[2])
        if position < len(spacer_values):
            spacer = spacer_values[position]
            if repeat[1] != spacer[0] or spacer[1] != repeat_values[
                    position + 1][0]:
                raise ValueError(
                    "likely_tandem ledger architecture is not contiguous"
                )
            reconstructed.append(spacer[2])
    reconstructed = "".join(reconstructed)
    if (
            repeat_values[0][0] != start - 1
            or repeat_values[-1][1] != end
            or len(reconstructed) != length
            or architecture.get("reconstructed_source_span_sha256")
            != hashlib.sha256(reconstructed.encode("utf-8")).hexdigest()
            or not _strict_json_equal(
                record.get("evidence"),
                evaluate_likely_tandem(spacer_copy).as_dict(),
            )):
        raise ValueError(
            "likely_tandem ledger source/evidence replay failed"
        )


def _validate_no_likely_tandem_accepted_arrays(model):
    """Independently reject qualifying rows left in accepted categories."""

    for array in model.get("arrays", []):
        if array.get("category") not in ("Bona-fide", "Possible"):
            continue
        architecture = array.get("architecture")
        if not isinstance(architecture, list):
            raise ValueError(
                "accepted array architecture is malformed for tandem replay"
            )
        spacer_slots = []
        derived_repeats = []
        derived_spacers = []
        derived_deletions = []
        for position, element in enumerate(architecture):
            if not isinstance(element, dict):
                raise ValueError(
                    "accepted array architecture is malformed for tandem replay"
                )
            kind = element.get("type")
            if position % 2 == 0:
                if kind != "repeat":
                    raise ValueError(
                        "accepted array architecture is malformed for tandem replay"
                    )
                derived_repeats.append({
                    key: value for key, value in element.items()
                    if key != "type"
                })
                continue
            if kind not in ("spacer", "spacer-deletion"):
                raise ValueError(
                    "accepted array architecture is malformed for tandem replay"
                )
            sequence = element.get("sequence")
            if not isinstance(sequence, str):
                raise ValueError(
                    "accepted array spacer is malformed for tandem replay"
                )
            spacer_slots.append(sequence)
            target = (
                derived_spacers if kind == "spacer" else derived_deletions
            )
            target.append({
                key: value for key, value in element.items()
                if key != "type"
            })
        if (
                not architecture
                or len(derived_repeats) != len(spacer_slots) + 1
                or array.get("repeats") != derived_repeats
                or array.get("spacers") != derived_spacers
                or array.get("spacer_deletions") != derived_deletions
                or array.get("repeat_count") != len(derived_repeats)
                or array.get("spacer_count") != len(spacer_slots)
                or array.get("non_deleted_spacer_count")
                != len(derived_spacers)
                or array.get("spacer_deletion_count")
                != len(derived_deletions)):
            raise ValueError(
                "accepted array architecture/tables differ for tandem replay"
            )
        if evaluate_likely_tandem(spacer_slots).likely_tandem:
            raise ValueError(
                "accepted array {!r} still qualifies as likely_tandem"
                .format(array.get("id"))
            )


def _validate_likely_tandem_extension(model, extension):
    source_sha256 = model.get("source", {}).get("sha256")
    expected_rule = {
        "name": LIKELY_TANDEM_RULE_NAME,
        "version": LIKELY_TANDEM_RULE_VERSION,
        "linkage": "complete_link_all_spacer_pairs",
        "identity": "1 - global_levenshtein/max(sequence_lengths)",
        "policy": LIKELY_TANDEM_POLICY.as_dict(),
    }
    if (
            extension.get("source_sha256") != source_sha256
            or extension.get("disposition") != "likely_tandem"
            or extension.get("tag") != TANDEM_TAG
            or extension.get("reported_category") != "Low score"
            or not _strict_json_equal(extension.get("rule"), expected_rule)
            or extension.get("external_tool_support_used_for_decision")
            is not False):
        raise ValueError("likely_tandem extension envelope is malformed")
    records = extension.get("arrays")
    stats = extension.get("transaction_stats")
    if (
            not isinstance(records, list)
            or not isinstance(stats, dict)
            or frozenset(stats) != _LIKELY_TANDEM_STATS_FIELDS
            or any(
                not _is_plain_int(value) or value < 0
                for value in stats.values()
            )
            or not _is_plain_int(extension.get("array_count"))
            or extension.get("array_count") != len(records)
            or stats.get("likely_tandem_count") != len(records)
            or stats.get("accepted_rows_before")
            != stats.get("accepted_rows_after") + len(records)
            or stats.get("accepted_rows_reviewed")
            != stats.get("accepted_rows_before")
            or stats.get("retained_count") != stats.get("accepted_rows_after")
            or stats.get("total_rows_before") != stats.get("total_rows_after")
            or stats.get("diagnostic_rows_after")
            != stats.get("diagnostic_rows_before") + len(records)
            or stats.get("accepted_rows_after")
            != model.get("summary", {}).get("accepted_array_count")):
        raise ValueError("likely_tandem extension counts are malformed")

    source_counts = {"Bona-fide": 0, "Possible": 0}
    seen_record_keys = set()
    for position, record in enumerate(records, 1):
        if not isinstance(record, dict):
            raise ValueError("likely_tandem ledger row is malformed")
        source_index = record.get("source_category_index")
        source_name = record.get("source_category_name")
        expected_source_name = {0: "Bona-fide", 2: "Possible"}.get(
            source_index
        )
        evidence = record.get("evidence")
        score = record.get("score")
        if (
                record.get("id") != "LT{}".format(position)
                or record.get("action") != "quarantined"
                or record.get("disposition") != "likely_tandem"
                or record.get("tag") != TANDEM_TAG
                or not _is_plain_int(source_index)
                or source_name != expected_source_name
                or record.get("target_category_index") != 4
                or record.get("target_category_name") != "Low score"
                or not _is_plain_int(record.get("row_ordinal"))
                or record.get("row_ordinal") < 0
                or isinstance(score, bool)
                or not isinstance(score, (int, float))
                or not math.isfinite(float(score))
                or not isinstance(evidence, dict)
                or evidence.get("decision") != "likely_tandem"
                or evidence.get("complete_link_threshold_satisfied") is not True):
            raise ValueError("likely_tandem ledger row is malformed")
        key = (
            json.dumps(
                record.get("source_interval"),
                sort_keys=True,
                separators=(",", ":"),
            ),
            source_index,
            json.dumps(
                record.get("locus_key"),
                sort_keys=True,
                separators=(",", ":"),
            ),
            record.get("row_ordinal"),
        )
        if key in seen_record_keys:
            raise ValueError("likely_tandem ledger contains duplicate rows")
        seen_record_keys.add(key)
        source_counts[source_name] += 1
        _validate_likely_tandem_ledger_architecture(
            record, model.get("source", {}).get("length")
        )
    if (
            stats.get("bona_fide_quarantined")
            != source_counts["Bona-fide"]
            or stats.get("possible_quarantined")
            != source_counts["Possible"]):
        raise ValueError(
            "likely_tandem source-category counts are malformed"
        )

    accepted_only = model.get("run", {}).get("parameters", {}).get(
        "param_accepted_only"
    ) is True
    tagged_arrays = [
        array for array in model.get("arrays", [])
        if TANDEM_TAG in _array_tags(array)
    ]
    _validate_no_likely_tandem_accepted_arrays(model)
    low_score_count = model.get("summary", {}).get(
        "category_counts", {}
    ).get("Low score")
    if accepted_only:
        if (
                tagged_arrays
                or low_score_count != 0
                or any(
                    record.get("reported_array_id") is not None
                    or record.get("core_projection")
                    != "omitted_by_accepted_only"
                    for record in records
                )):
            raise ValueError(
                "accepted-only likely_tandem projection is malformed"
            )
        return
    if low_score_count != stats.get("diagnostic_rows_after"):
        raise ValueError("likely_tandem Low score count is malformed")
    core_by_id = {}
    for array in tagged_arrays:
        core_by_id.setdefault(array.get("id"), []).append(array)
    consumed = set()
    for record in records:
        reported_id = record.get("reported_array_id")
        matches = core_by_id.get(reported_id, ())
        if (
                record.get("core_projection") != "reported_as_low_score"
                or not isinstance(reported_id, str)
                or len(matches) != 1
                or reported_id in consumed):
            raise ValueError(
                "likely_tandem core projection is not one-to-one"
            )
        array = matches[0]
        basis = array.get("classification_basis")
        transition = basis.get("transition") if isinstance(
            basis, dict
        ) else None
        if (
                array.get("category") != "Low score"
                or not _strict_json_equal(
                    array.get("source_interval"),
                    record.get("source_interval"),
                )
                or not _strict_json_equal(
                    array.get("legacy_interval_key"),
                    record.get("locus_key"),
                )
                or isinstance(array.get("certainty_score"), bool)
                or not isinstance(array.get("certainty_score"), (int, float))
                or not math.isfinite(float(array["certainty_score"]))
                or not _strict_json_equal(
                    array.get("certainty_score"), record.get("score")
                )
                or not isinstance(basis, dict)
                or basis.get("method") != CLASSIFICATION_BASIS_METHOD
                or basis.get("disposition") != "likely_tandem"
                or basis.get("tag") != TANDEM_TAG
                or not _strict_json_equal(
                    basis.get("likely_tandem"), record.get("evidence")
                )
                or not isinstance(transition, dict)
                or transition.get("rule_name") != TRANSACTION_RULE_NAME
                or transition.get("rule_version") != TRANSACTION_RULE_VERSION
                or not _strict_json_equal(
                    transition.get("source_category_index"),
                    record.get("source_category_index"),
                )
                or transition.get("source_category_name")
                != record.get("source_category_name")
                or transition.get("target_category_index") != 4
                or transition.get("target_category_name") != "Low score"
                or not _strict_json_equal(
                    transition.get("row_ordinal"),
                    record.get("row_ordinal"),
                )
                or not _strict_json_equal(
                    transition.get("locus_key"), record.get("locus_key")
                )
                or transition.get("reported_disposition")
                != "likely_tandem"
                or transition.get("reported_tag") != TANDEM_TAG
                or not _likely_tandem_core_architecture_matches(
                    array, record
                )):
            raise ValueError(
                "likely_tandem core row differs from its ledger"
            )
        consumed.add(reported_id)
    if consumed != set(core_by_id):
        raise ValueError(
            "tandem-tagged core arrays are missing ledger rows"
        )


def _extension(model, name):
    extensions = model.get("extensions")
    if extensions is None:
        return None
    if not isinstance(extensions, dict):
        raise ValueError("report extensions must be a mapping")
    value = extensions.get(name)
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("{} extension must be a mapping".format(name))

    expected_name, expected_version, expected_status = (
        _EXTENSION_CONTRACTS[name]
    )
    schema = value.get("schema")
    # Unknown or future self-versioned extensions are deliberately ignored;
    # an exact supported schema is validated before any field is rendered.
    if (
        not isinstance(schema, dict)
        or schema.get("name") != expected_name
        or schema.get("version") != expected_version
    ):
        return None
    if value.get("status") != expected_status:
        raise ValueError(
            "{} extension has invalid status".format(name)
        )
    if name == "source_geometry":
        _validate_source_geometry_extension(model, value)
    elif name == "biological_loci":
        _validate_biological_loci_extension(model, value)
    elif name == "likely_tandem":
        _validate_likely_tandem_extension(model, value)
    return value


def _classification_interpretation(model):
    semantics = _extension(model, "classification_semantics") or {}
    score = semantics.get("certainty_score")
    reconstruction = semantics.get("validation")

    reconstruction_text = (
        "Exact source reconstruction checks that reported coordinates and DNA "
        "agree with the original input sequence. It does not establish "
        "biological truth or classification correctness."
    )
    if isinstance(reconstruction, dict) and isinstance(
        reconstruction.get("interpretation"), str
    ):
        reconstruction_text = reconstruction["interpretation"]

    score_bands = ""
    if isinstance(score, dict):
        conditions = []
        for band in score.get("initial_score_bands", []):
            if isinstance(band, dict) and isinstance(
                band.get("condition"), str
            ):
                conditions.append(
                    "<code>{}</code>".format(_escape(band["condition"]))
                )
        assignment = score.get("category_assignment")
        assignment_text = None
        if isinstance(assignment, dict) and isinstance(
            assignment.get("interpretation"), str
        ):
            assignment_text = assignment["interpretation"]
        if conditions:
            score_bands = """\
        <p>
          <strong>Initial score bands:</strong> {conditions}.{assignment}
        </p>""".format(
                conditions=", ".join(conditions),
                assignment=(
                    " " + _escape(assignment_text)
                    if assignment_text else ""
                ),
            )

    return """\
    <section aria-labelledby="interpretation-heading">
      <h2 id="interpretation-heading">How to interpret this report</h2>
      <div class="interpretation-card">
        <ul>
          <li>{reconstruction_text}</li>
        </ul>
{score_bands}
      </div>
    </section>""".format(
        reconstruction_text=_escape(reconstruction_text),
        score_bands=score_bands,
    )


def _source_geometry_by_id(model):
    geometry = _extension(model, "source_geometry")
    if geometry is None:
        return {}
    result = {}
    for item in geometry.get("arrays", []):
        if not isinstance(item, dict):
            continue
        array_id = item.get("array_id")
        if isinstance(array_id, str):
            result[array_id] = item
    return result


def _format_metric(value):
    if value is None:
        return "Not available"
    if isinstance(value, float):
        return "{:.4g}".format(value)
    return str(value)


def _format_number_list(values):
    if not isinstance(values, list):
        return "Not available"
    return json.dumps(values, allow_nan=False, separators=(",", ":"))


def _source_geometry_details(geometry):
    if not isinstance(geometry, dict):
        return ""

    statistics = geometry.get("statistics")
    if not isinstance(statistics, dict):
        statistics = {}
    gap_stats = statistics.get("spacer_gap_length_bp")
    if not isinstance(gap_stats, dict):
        gap_stats = {}
    non_deleted_stats = statistics.get("non_deleted_spacer_gap_length_bp")
    if not isinstance(non_deleted_stats, dict):
        non_deleted_stats = {}
    period_stats = statistics.get("repeat_start_to_start_period_bp")
    if not isinstance(period_stats, dict):
        period_stats = {}

    advisories = []
    for advisory in geometry.get("advisories", []):
        if not isinstance(advisory, dict):
            continue
        code = advisory.get("code", "source_geometry_advisory")
        message = advisory.get("message", "")
        advisories.append(
            """\
                <li>
                  <strong>{code}</strong>: {message}
                  <span class="muted">Decision effect: {effect}.</span>
                </li>""".format(
                code=_escape(code),
                message=_escape(message),
                effect=_escape(advisory["decision_effect"]),
            )
        )
    if advisories:
        advisory_output = """\
              <ul class="advisory-list">
{items}
              </ul>""".format(items="\n".join(advisories))
    else:
        advisory_output = """\
              <p class="muted">
                No non-decision source-geometry advisories were emitted.
              </p>"""

    return """\
              <details class="nested-details geometry-details">
                <summary>Accepted-array source geometry and diagnostics</summary>
                <p class="muted">
                  Derived only from exactly reconstructed source intervals.
                  These diagnostics do not change the category or score.
                </p>
                <dl class="array-metadata">
                  <div><dt>Repeat starts (1-based)</dt>
                    <dd><code>{repeat_starts}</code></dd></div>
                  <div><dt>Repeat inclusive ends (1-based)</dt>
                    <dd><code>{repeat_ends}</code></dd></div>
                  <div><dt>Repeat start-to-start periods (bp)</dt>
                    <dd><code>{periods}</code></dd></div>
                  <div><dt>All spacer/gap slots (bp; deletions are 0)</dt>
                    <dd><code>{all_gaps}</code></dd></div>
                  <div><dt>Non-deleted spacer/gap slots (bp)</dt>
                    <dd><code>{non_deleted_gaps}</code></dd></div>
                  <div><dt>Spacer-deletion slot indices</dt>
                    <dd><code>{deletion_slots}</code></dd></div>
                  <div><dt>All-gap mean / median / CV</dt>
                    <dd>{gap_mean} / {gap_median} / {gap_cv}</dd></div>
                  <div><dt>Non-deleted-gap mean / median</dt>
                    <dd>{non_deleted_mean} / {non_deleted_median}</dd></div>
                  <div><dt>Repeat-period median (bp)</dt>
                    <dd>{period_median}</dd></div>
                </dl>
                <h3>Non-decision advisories</h3>
{advisory_output}
              </details>""".format(
        repeat_starts=_escape(_format_number_list(
            geometry.get("repeat_starts_1_based")
        )),
        repeat_ends=_escape(_format_number_list(
            geometry.get("repeat_inclusive_ends_1_based")
        )),
        periods=_escape(_format_number_list(
            geometry.get("repeat_start_to_start_periods_bp")
        )),
        all_gaps=_escape(_format_number_list(
            geometry.get("spacer_gap_lengths_bp")
        )),
        non_deleted_gaps=_escape(_format_number_list(
            geometry.get("non_deleted_spacer_gap_lengths_bp")
        )),
        deletion_slots=_escape(_format_number_list(
            geometry.get("spacer_deletion_slot_indices")
        )),
        gap_mean=_escape(_format_metric(gap_stats.get("mean"))),
        gap_median=_escape(_format_metric(gap_stats.get("median"))),
        gap_cv=_escape(_format_metric(
            gap_stats.get("coefficient_of_variation")
        )),
        non_deleted_mean=_escape(_format_metric(
            non_deleted_stats.get("mean")
        )),
        non_deleted_median=_escape(_format_metric(
            non_deleted_stats.get("median")
        )),
        period_median=_escape(_format_metric(period_stats.get("median"))),
        advisory_output=advisory_output,
    )


def _biological_loci_section(model):
    extension = _extension(model, "biological_loci")
    if extension is None:
        return ""

    counts = extension.get("counts")
    if not isinstance(counts, dict):
        counts = {}
    rows = []
    for locus in extension.get("loci", []):
        if not isinstance(locus, dict):
            continue
        interval = locus.get("source_interval")
        if not isinstance(interval, dict):
            continue
        member_links = []
        for member_id in locus.get("member_ids", []):
            if isinstance(member_id, str):
                member_links.append(
                    '<a href="#{target}"><code>{label}</code></a>'.format(
                        target=_escape(member_id),
                        label=_escape(member_id),
                    )
                )
        rows.append(
            """\
              <tr>
                <th scope="row"><code>{locus_id}</code></th>
                <td class="numeric">{start:,}&ndash;{end:,}</td>
                <td class="numeric">{member_count}</td>
                <td class="member-links">{members}</td>
              </tr>""".format(
                locus_id=_escape(locus.get("id", "")),
                start=interval.get("start", 0),
                end=interval.get("end", 0),
                member_count=_escape(locus.get("member_count", "")),
                members=", ".join(member_links),
            )
        )

    return """\
    <section aria-labelledby="loci-heading">
      <h2 id="loci-heading">Heuristic biological-locus presentation groups</h2>
      <div class="loci-card">
        <p>
          This optional view is a <strong>heuristic presentation grouping</strong>,
          not a detector decision. It does not merge, suppress, reclassify, or
          select one member over another. All members and the raw accepted-array
          count remain visible.
        </p>
        <dl class="source-grid">
          <div><dt>Raw accepted arrays</dt><dd>{raw_count}</dd></div>
          <div><dt>Heuristic locus groups</dt><dd>{locus_count}</dd></div>
          <div><dt>Singleton groups</dt><dd>{singleton_count}</dd></div>
          <div><dt>Multi-member groups</dt><dd>{multi_count}</dd></div>
        </dl>
        <div class="table-scroll" tabindex="0"
             aria-label="Scrollable heuristic biological-locus groups">
          <table id="biological-loci">
            <caption>All accepted-array members in every heuristic group</caption>
            <thead>
              <tr>
                <th scope="col">Stable group ID</th>
                <th scope="col">Source envelope</th>
                <th scope="col">Members</th>
                <th scope="col">All member array IDs</th>
              </tr>
            </thead>
            <tbody>
{rows}
            </tbody>
          </table>
        </div>
      </div>
    </section>""".format(
        raw_count=_escape(counts.get("raw_accepted_array_count", "")),
        locus_count=_escape(counts.get("biological_locus_count", "")),
        singleton_count=_escape(counts.get("singleton_locus_count", "")),
        multi_count=_escape(counts.get("multi_member_locus_count", "")),
        rows="\n".join(rows),
    )


def _likely_tandem_section(model):
    extension = _extension(model, "likely_tandem")
    accepted_only = model.get("run", {}).get("parameters", {}).get(
        "param_accepted_only"
    ) is True
    tagged = any(
        TANDEM_TAG in _array_tags(array)
        for array in model.get("arrays", [])
    )
    _validate_no_likely_tandem_accepted_arrays(model)
    semantics_schema = model.get("extensions", {}).get(
        "classification_semantics", {}
    ).get("schema", {}) if isinstance(
        model.get("extensions"), dict
    ) else {}
    audit_contract = model.get("extensions", {}).get(
        "classification_semantics", {}
    ).get("likely_tandem_audit", {}) if isinstance(
        model.get("extensions"), dict
    ) else {}
    tandem_aware_producer = (
        semantics_schema.get("name") == CLASSIFICATION_SCHEMA_NAME
        and semantics_schema.get("version") == CLASSIFICATION_SCHEMA_VERSION
        and isinstance(audit_contract, dict)
        and audit_contract.get("required") is True
        and audit_contract.get("extension") == "likely_tandem"
    )
    if extension is None:
        if accepted_only or tagged or tandem_aware_producer:
            raise ValueError(
                "report requires a supported likely_tandem audit ledger"
            )
        return ""

    rows = []
    for record in extension.get("arrays", []):
        interval = record["source_interval"]
        evidence = record.get("evidence", {})
        reported_id = record.get("reported_array_id")
        if reported_id is None:
            reported = "Omitted by accepted-only projection"
        else:
            reported = '<a href="#{target}"><code>{label}</code></a>'.format(
                target=_escape(reported_id),
                label=_escape(reported_id),
            )
        rows.append(
            """\
              <tr>
                <th scope="row"><code>{ledger_id}</code></th>
                <td class="numeric">{start:,}&ndash;{end:,}</td>
                <td>{source_category}</td>
                <td class="numeric">{spacer_count}</td>
                <td class="numeric">{minimum_identity}</td>
                <td>{projection}</td>
                <td>{reported}</td>
              </tr>""".format(
                ledger_id=_escape(record.get("id", "")),
                start=interval["start"],
                end=interval["end"],
                source_category=_escape(record.get(
                    "source_category_name", ""
                )),
                spacer_count=_escape(evidence.get("spacer_count", "")),
                minimum_identity=_escape(_format_metric(
                    evidence.get("minimum_pairwise_edit_identity")
                )),
                projection=_escape(record.get("core_projection", "")),
                reported=reported,
            )
        )
    projection_note = (
        "This accepted-only report omits demoted Low score rows from the "
        "core array list; the source-bound audit ledger remains visible."
        if accepted_only else
        "Every ledger row maps one-to-one to a tagged Low score core array."
    )
    return """\
    <section aria-labelledby="tandem-audit-heading">
      <h2 id="tandem-audit-heading">Tandem demotion audit</h2>
      <div class="loci-card">
        <p>
          <strong>{count}</strong> accepted candidate(s) met the bounded
          complete-link low-spacer-diversity rule. Their score is preserved,
          their public category is <strong>Low score</strong>, and their tag is
          <span class="tag-badge">tandem</span>. External-tool agreement was
          not used for this decision.
        </p>
        <p class="muted">{projection_note}</p>
        <div class="table-scroll" tabindex="0"
             aria-label="Scrollable tandem demotion audit ledger">
          <table id="likely-tandem-ledger">
            <caption>Source-bound likely-tandem transaction ledger</caption>
            <thead>
              <tr>
                <th scope="col">Ledger ID</th>
                <th scope="col">Source interval</th>
                <th scope="col">Prior category</th>
                <th scope="col">Spacers</th>
                <th scope="col">Minimum identity</th>
                <th scope="col">Core projection</th>
                <th scope="col">Reported array</th>
              </tr>
            </thead>
            <tbody>
{rows}
            </tbody>
          </table>
        </div>
      </div>
    </section>""".format(
        count=_escape(extension.get("array_count", 0)),
        projection_note=_escape(projection_note),
        rows="\n".join(rows),
    )


def _array_cards(model):
    geometry_by_id = _source_geometry_by_id(model)
    competing_by_id = competing_terminal_extent_summaries(model)
    cards = []
    for array in model["arrays"]:
        interval = array["source_interval"]
        orientation = _orientation_text(array["orientation"])
        cards.append(
            """\
        <article id="{array_id}" class="array-card"
                 data-category="{category_key}">
          <details>
            <summary>
              <span class="badge badge-{category_key}">{category}</span>
              {tags}
              <strong>{array_id}</strong>
              <span>{start:,}&ndash;{end:,}</span>
              <span>{orientation}</span>
            </summary>
            <div class="array-body">
              <dl class="array-metadata">
                <div><dt>Coordinate contract</dt>
                  <dd>1-based inclusive intervals; 0-based interbase deletion boundaries</dd></div>
                <div><dt>Confidence score</dt>
                  <dd>{score}</dd></div>
                <div><dt>Tags</dt>
                  <dd>{tags}</dd></div>
                <div><dt>Consensus repeat</dt>
                  <dd><code>{consensus}</code></dd></div>
                <div><dt>Complete spacer deletions</dt>
                  <dd>{spacer_deletion_count}</dd></div>
                <div><dt>Exact source reconstruction</dt>
                  <dd><span class="validation validation-valid">
                    Reconstructed exactly from source DNA
                  </span></dd></div>
              </dl>
{boundary_family}
{possible_reconstruction}
{contextual_spacer}
{competing_extents}
{source_geometry}
{architecture}
{components}
{features}
              <details class="nested-details">
                <summary>Complete source DNA ({length:,} bases)</summary>
                <pre class="sequence-block">{source_sequence}</pre>
              </details>
            </div>
          </details>
        </article>""".format(
                array_id=_escape(array["id"]),
                category_key=_escape(array["category_key"]),
                category=_escape(array["category"]),
                tags=_tag_badges(array),
                start=interval["start"],
                end=interval["end"],
                orientation=_escape(orientation),
                score=_escape(_format_score(array["certainty_score"])),
                consensus=_escape(array["consensus_repeat"]),
                spacer_deletion_count=array.get("spacer_deletion_count", 0),
                boundary_family=_boundary_family_details(array),
                possible_reconstruction=_possible_reconstruction_details(array),
                contextual_spacer=_contextual_spacer_details(array),
                competing_extents=render_competing_extent_details(
                    competing_by_id.get(array["id"], [])
                ),
                source_geometry=_source_geometry_details(
                    geometry_by_id.get(array["id"])
                ),
                architecture=_architecture(array),
                components=_component_table(array),
                features=_feature_tables(array),
                length=interval["length"],
                source_sequence=_escape(
                    _wrap_sequence(array["source_sequence"])
                ),
            )
        )
    if not cards:
        return """\
        <p class="empty-state">
          No CRISPR arrays were reported. The source record and zero-array
          source-reconstruction result are still captured in this report.
        </p>"""
    return "\n".join(cards)


_STYLES = r"""
    :root {
      color-scheme: light dark;
      --bg: #f7f8fa;
      --surface: #ffffff;
      --text: #18212b;
      --muted: #536273;
      --border: #c8d0d9;
      --accent: #2457a6;
      --focus: #d97706;
      --repeat: #1d70b8;
      --spacer: #f3c13a;
      --deletion: #c2415d;
      --valid: #176b3a;
      --bf-bg: #d8f3dc;
      --bf-fg: #14532d;
      --possible-bg: #fff3bf;
      --possible-fg: #6b4f00;
      --low-bg: #e5e7eb;
      --low-fg: #374151;
    }
    @media (prefers-color-scheme: dark) {
      :root {
        --bg: #111827;
        --surface: #1f2937;
        --text: #f3f4f6;
        --muted: #cbd5e1;
        --border: #64748b;
        --accent: #93c5fd;
        --focus: #fbbf24;
        --valid: #86efac;
        --bf-bg: #14532d;
        --bf-fg: #dcfce7;
        --possible-bg: #713f12;
        --possible-fg: #fef3c7;
        --low-bg: #374151;
        --low-fg: #f3f4f6;
      }
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI",
        sans-serif;
      line-height: 1.5;
    }
    a { color: var(--accent); }
    a:focus-visible, button:focus-visible, input:focus-visible,
    [tabindex="0"]:focus-visible, summary:focus-visible {
      outline: 3px solid var(--focus);
      outline-offset: 3px;
    }
    .skip-link {
      position: absolute;
      left: 0.5rem;
      top: -5rem;
      background: var(--surface);
      color: var(--text);
      padding: 0.75rem;
      z-index: 10;
    }
    .skip-link:focus { top: 0.5rem; }
    header, main, footer {
      width: min(96rem, calc(100% - 2rem));
      margin-inline: auto;
    }
    header { padding-block: 2rem 1rem; }
    h1, h2 { line-height: 1.2; }
    .source-card, .interpretation-card, .loci-card,
    .filters, .array-card, .empty-state {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 0.5rem;
      padding: 1rem;
      margin-block: 1rem;
    }
    .source-grid, .array-metadata {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(15rem, 1fr));
      gap: 0.75rem 1.5rem;
      margin: 0;
    }
    .source-grid div, .array-metadata div { min-width: 0; }
    dt { color: var(--muted); font-weight: 700; }
    dd { margin: 0.2rem 0 0; overflow-wrap: anywhere; }
    .coordinate-contract {
      border-left: 0.35rem solid var(--accent);
      padding-left: 0.8rem;
    }
    .filters { display: flex; flex-wrap: wrap; gap: 0.75rem 1rem; }
    .filters legend { font-weight: 700; padding-inline: 0.25rem; }
    .filter-option { display: inline-flex; gap: 0.35rem; align-items: center; }
    .filter-actions { display: inline-flex; gap: 0.5rem; margin-left: auto; }
    button {
      background: var(--surface);
      color: var(--text);
      border: 1px solid var(--border);
      border-radius: 0.25rem;
      padding: 0.35rem 0.65rem;
      cursor: pointer;
    }
    .badge, .validation, .tag-badge {
      display: inline-block;
      border-radius: 999px;
      font-size: 0.875rem;
      font-weight: 700;
      padding: 0.15rem 0.55rem;
      white-space: nowrap;
    }
    .badge-bona-fide { background: var(--bf-bg); color: var(--bf-fg); }
    .badge-possible { background: var(--possible-bg); color: var(--possible-fg); }
    .badge-low-score { background: var(--low-bg); color: var(--low-fg); }
    .tag-badge { background: #e6d8ff; color: #45256f; }
    .validation-valid {
      border: 1px solid var(--valid);
      color: var(--valid);
    }
    .status-line, .muted, .scroll-hint { color: var(--muted); }
    .table-scroll, .component-table-wrap {
      overflow: auto;
      max-width: 100%;
      border: 1px solid var(--border);
      border-radius: 0.35rem;
      background: var(--surface);
    }
    table { width: 100%; border-collapse: collapse; }
    caption {
      text-align: left;
      font-weight: 700;
      padding: 0.75rem;
    }
    th, td {
      border-top: 1px solid var(--border);
      padding: 0.6rem 0.75rem;
      text-align: left;
      vertical-align: top;
    }
    thead th {
      border-top: 0;
      background: var(--surface);
      position: sticky;
      top: 0;
      z-index: 1;
    }
    .numeric { font-variant-numeric: tabular-nums; white-space: nowrap; }
    code, pre {
      font-family: ui-monospace, SFMono-Regular, Consolas, monospace;
      font-size: 0.9em;
    }
    code { overflow-wrap: anywhere; }
    .array-card { padding: 0; scroll-margin-top: 1rem; }
    .array-card > details > summary {
      cursor: pointer;
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 0.6rem 1rem;
      padding: 1rem;
    }
    .array-body { border-top: 1px solid var(--border); padding: 1rem; }
    .architecture-legend {
      display: flex;
      gap: 1rem;
      margin-top: 1.25rem;
      font-size: 0.9rem;
    }
    .architecture-legend span { display: inline-flex; align-items: center; gap: 0.3rem; }
    .architecture-legend i {
      display: inline-block;
      width: 1.5rem;
      height: 0.8rem;
      border: 1px solid var(--text);
    }
    .legend-repeat {
      background: repeating-linear-gradient(
        135deg, var(--repeat), var(--repeat) 3px,
        transparent 3px, transparent 6px
      );
    }
    .legend-spacer { background: var(--spacer); }
    .legend-deletion {
      width: 0.3rem !important;
      background: var(--deletion);
    }
    .architecture-viewport {
      overflow-x: auto;
      padding: 0.75rem 0.15rem;
      margin-top: 0.3rem;
    }
    .architecture {
      display: flex;
      align-items: center;
      gap: 2px;
      height: 3rem;
      min-width: max(100%, calc(var(--component-count) * 10px));
    }
    .segment {
      display: block;
      height: 2rem;
      border: 1px solid var(--text);
    }
    .segment-repeat {
      min-width: 6px;
      max-width: 3rem;
      background: repeating-linear-gradient(
        135deg, var(--repeat), var(--repeat) 3px,
        transparent 3px, transparent 6px
      );
      clip-path: polygon(50% 0, 100% 50%, 50% 100%, 0 50%);
    }
    .segment-spacer {
      min-width: 10px;
      max-width: 8rem;
      background: var(--spacer);
    }
    .segment-spacer-deletion {
      flex: 0 0 4px;
      min-width: 4px;
      height: 2.75rem;
      border: 0;
      border-radius: 2px;
      background: var(--deletion);
      box-shadow: 0 0 0 1px var(--surface), 0 0 0 2px var(--deletion);
    }
    .spacer-deletion-row th { color: var(--deletion); }
    .deletion-event-label { font-weight: 700; }
    .nested-details { margin-top: 1rem; }
    .nested-details > summary { cursor: pointer; font-weight: 700; }
    .geometry-details h3 { margin-bottom: 0.25rem; }
    .advisory-list { margin-top: 0.25rem; }
    .member-links { overflow-wrap: anywhere; }
    .component-table-wrap { max-height: 26rem; margin-top: 0.6rem; }
    .sequence-inline { word-break: break-all; }
    .sequence-block {
      background: var(--bg);
      border: 1px solid var(--border);
      border-radius: 0.35rem;
      max-height: 26rem;
      overflow: auto;
      padding: 0.75rem;
      white-space: pre-wrap;
      overflow-wrap: anywhere;
    }
    .feature-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(20rem, 1fr));
      gap: 1rem;
      margin-top: 0.6rem;
    }
    footer { color: var(--muted); padding-block: 2rem; }
    [hidden] { display: none !important; }
    @media print {
      :root { color-scheme: light; }
      body { background: #fff; color: #000; }
      .skip-link, .filters, .status-line, .scroll-hint { display: none; }
      header, main, footer { width: 100%; }
      .array-card, .source-card, .interpretation-card, .loci-card {
        break-inside: avoid;
      }
      .table-scroll, .component-table-wrap, .architecture-viewport,
      .sequence-block { overflow: visible; max-height: none; }
      details > * { display: block !important; }
    }
"""


_SCRIPT = r"""
    (() => {
      const reportElement = document.getElementById("report-data");
      window.CRISPRidentifyReport = JSON.parse(reportElement.textContent);

      const filters = Array.from(
        document.querySelectorAll("[data-category-filter]")
      );
      const filterable = Array.from(
        document.querySelectorAll("[data-category]")
      );
      const status = document.getElementById("visible-status");
      const totalArrays = window.CRISPRidentifyReport.summary.array_count;

      function applyFilters() {
        const visibleCategories = new Set(
          filters.filter(input => input.checked)
            .map(input => input.dataset.categoryFilter)
        );
        let visibleCards = 0;
        filterable.forEach(element => {
          const visible = visibleCategories.has(element.dataset.category);
          element.hidden = !visible;
          if (visible && element.matches("article.array-card")) {
            visibleCards += 1;
          }
        });
        status.textContent =
          `Showing ${visibleCards} of ${totalArrays} reported arrays with ` +
          `exact source reconstruction.`;
      }

      filters.forEach(input => input.addEventListener("change", applyFilters));
      document.querySelectorAll("[data-filter-action]").forEach(button => {
        button.addEventListener("click", () => {
          const checked = button.dataset.filterAction === "all";
          filters.forEach(input => { input.checked = checked; });
          applyFilters();
        });
      });
      applyFilters();
    })();
"""


def render_html_report(model):
    """Render a canonical report as one offline HTML document."""

    schema = model.get("schema", {})
    if schema.get("name") != SCHEMA_NAME:
        raise ValueError("unsupported report schema")
    if model.get("validation", {}).get("status") != "valid":
        raise ValueError(
            "refusing to render a report without exact source reconstruction"
        )
    for array in model.get("arrays", []):
        if array.get("validation", {}).get("status") != "valid":
            raise ValueError(
                "refusing to render an array without exact source reconstruction"
            )

    source = model["source"]
    summary = model["summary"]
    return """\
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="generator" content="CRISPRidentify">
  <title>CRISPRidentify report: {source_id}</title>
  <style>
{styles}
  </style>
</head>
<body>
  <a class="skip-link" href="#main-content">Skip to report content</a>
  <header>
    <h1>CRISPRidentify report with exact source reconstruction</h1>
    <p>
      <strong>{array_count}</strong> arrays reported for
      <strong>{source_id}</strong>.
    </p>
  </header>
  <main id="main-content">
    <section aria-labelledby="source-heading">
      <h2 id="source-heading">Source record</h2>
      <div class="source-card">
        <dl class="source-grid">
          <div><dt>FASTA ID</dt><dd>{source_id}</dd></div>
          <div><dt>Full header</dt><dd>{source_header}</dd></div>
          <div><dt>Length</dt><dd>{source_length:,} bases</dd></div>
          <div><dt>Sequence SHA-256</dt><dd><code>{source_hash}</code></dd></div>
          <div><dt>Schema</dt><dd>{schema_name} {schema_version}</dd></div>
          <div><dt>Source reconstruction check</dt>
            <dd><span class="validation validation-valid">
              {validated_count} of {array_count} array sequences reconstructed exactly
            </span></dd></div>
        </dl>
        <p class="coordinate-contract">
          Displayed intervals are <strong>1-based and inclusive</strong> on the
          <strong>input-forward source sequence</strong>. A complete spacer
          deletion has no interval and is displayed instead as an explicitly
          labelled <strong>0-based interbase boundary</strong>. Predicted
          biological orientation is separate and never reverses source
          coordinates.
        </p>
        <p>
          Exact reconstruction verifies reported coordinates and DNA against
          this input record only. It does <strong>not</strong> establish that an
          array is biologically real, correctly classified, or associated with
          a calibrated probability.
        </p>
      </div>
    </section>
{interpretation}
    <section aria-labelledby="summary-heading">
      <h2 id="summary-heading">Array summary</h2>
{filters}
{summary_table}
    </section>
{biological_loci}
{likely_tandem_audit}
{cas_section}
    <section aria-labelledby="details-heading">
      <h2 id="details-heading">Array details</h2>
{cards}
    </section>
  </main>
  <footer>
    Generated as a self-contained offline report. No network resources are
    loaded.
  </footer>
  <script id="report-data" type="application/json">{report_json}</script>
  <script>
{script}
  </script>
</body>
</html>
""".format(
        source_id=_escape(source["id"]),
        source_header=_escape(source["header"]),
        source_length=source["length"],
        source_hash=_escape(source["sha256"]),
        schema_name=_escape(model["schema"]["name"]),
        schema_version=_escape(model["schema"]["version"]),
        array_count=summary["array_count"],
        validated_count=model["validation"]["validated_array_count"],
        filters=_filters(model),
        summary_table=_summary_table(model),
        interpretation=_classification_interpretation(model),
        biological_loci=_biological_loci_section(model),
        likely_tandem_audit=_likely_tandem_section(model),
        cas_section=render_cas_section(model),
        cards=_array_cards(model),
        report_json=_script_safe_json(model),
        styles=_STYLES + CAS_STYLES,
        script=_SCRIPT,
    )




class CanonicalReportMaker:
    """Write the canonical JSON report and, optionally, its offline HTML view."""

    def __init__(
        self,
        file_path,
        result_path,
        categories,
        non_array_data,
        list_feature_names,
        parameters=None,
        flags=None,
        write_html=False,
        sequence_topology_evidence=None,
        circular_array_join_result=None,
        publication=None,
        likely_tandem_records=None,
        likely_tandem_stats=None,
        likely_tandem_source_categories=None,
    ):
        self.file_path = file_path
        self.result_path = os.fspath(result_path)
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_feature_names = list_feature_names
        self.parameters = parameters
        self.flags = flags
        self.write_html = bool(write_html)
        self.sequence_topology_evidence = sequence_topology_evidence
        self.circular_array_join_result = circular_array_join_result
        self.likely_tandem_records = likely_tandem_records
        self.likely_tandem_stats = likely_tandem_stats
        self.likely_tandem_source_categories = (
            likely_tandem_source_categories
        )
        self.publication = publication
        self.model = None
        self.json_path = os.path.join(self.result_path, "report.json")
        self.html_path = os.path.join(self.result_path, "report.html")
        self._write()

    def _write(self):
        remove_completion_marker(self.json_path)
        (
            self.likely_tandem_records,
            self.likely_tandem_stats,
            self.likely_tandem_source_categories,
        ) = resolve_likely_tandem_audit_inputs(
            self.categories,
            self.likely_tandem_records,
            self.likely_tandem_stats,
            self.likely_tandem_source_categories,
            accepted_only=bool(
                (self.parameters or {}).get("param_accepted_only", False)
            ),
        )
        _source_header, source_sequence = read_single_fasta(self.file_path)
        core_model = ReportModelBuilder(
            file_path=self.file_path,
            categories=self.categories,
            non_array_data=self.non_array_data,
            list_feature_names=self.list_feature_names,
            parameters=self.parameters,
            flags=self.flags,
            provenance_categories=self.likely_tandem_source_categories,
        ).build()
        include_biological_loci = bool(
            (self.flags or {}).get("flag_biological_loci", False)
        )
        self.model = attach_report_extensions(
            core_model,
            include_biological_loci=include_biological_loci,
            sequence_topology_evidence=self.sequence_topology_evidence,
            circular_array_join_result=self.circular_array_join_result,
            likely_tandem_records=self.likely_tandem_records,
            likely_tandem_stats=self.likely_tandem_stats,
            likely_tandem_source_categories=(
                self.likely_tandem_source_categories
            ),
            source_sequence=source_sequence,
            competing_terminal_extent_context={
                "accepted_categories": self.categories,
                "provenance_categories": self.likely_tandem_source_categories,
                "parameters": self.parameters or {},
                "feature_names": self.list_feature_names,
            },
        )
        _publish_report_model(
            self.model, self.result_path, self.write_html, self.publication
        )

    def output(self):
        return {
            "model": self.model,
            "json_path": self.json_path,
            "html_path": self.html_path if self.write_html else None,
        }
