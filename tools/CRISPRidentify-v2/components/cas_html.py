"""Accessible, provider-neutral visualization of integrated Cas calls."""

import html
import importlib.util
from pathlib import Path

try:
    from components.cas_output import normalize_cas_data
except ModuleNotFoundError as error:
    # The distribution tests also load this module directly by file path,
    # without first creating the standalone package namespace.
    if error.name != "components.cas_output":
        raise
    _SPEC = importlib.util.spec_from_file_location(
        "v2_cas_output_for_html",
        Path(__file__).with_name("cas_output.py"),
    )
    _CAS_OUTPUT = importlib.util.module_from_spec(_SPEC)
    _SPEC.loader.exec_module(_CAS_OUTPUT)
    normalize_cas_data = _CAS_OUTPUT.normalize_cas_data


def _escape(value):
    return html.escape(str(value), quote=True)


def _validate(model):
    source = model.get("source")
    if not isinstance(source, dict):
        raise ValueError("Cas report source must be an object")
    source_id = source.get("id")
    source_length = source.get("length")
    if not isinstance(source_id, str) or not source_id:
        raise ValueError("Cas report source ID must be a non-empty string")
    if (
            isinstance(source_length, bool)
            or not isinstance(source_length, int)
            or source_length < 1):
        raise ValueError("Cas report source length must be a positive integer")
    cas = normalize_cas_data(
        model.get("cas", {"status": "not_requested", "systems": []}),
        source_id=source_id,
        source_length=source_length,
    )
    return cas, cas["systems"]


def _interval_style(start, end, container_start, container_end):
    span = container_end - container_start + 1
    left = 100.0 * (start - container_start) / span
    width = 100.0 * (end - start + 1) / span
    return "left:{:.6f}%;width:{:.6f}%".format(left, width)


def _gene_lanes(genes):
    lane_ends = []
    placements = []
    for gene in genes:
        lane = next(
            (
                index for index, lane_end in enumerate(lane_ends)
                if gene["start"] > lane_end
            ),
            None,
        )
        if lane is None:
            lane = len(lane_ends)
            lane_ends.append(gene["end"])
        else:
            lane_ends[lane] = gene["end"]
        placements.append((gene, lane))
    return placements, max(1, len(lane_ends))


def _score_text(gene):
    if gene["score"] is None:
        return "not reported"
    kind = gene["score_kind"] or "provider score"
    return "{:.3f} ({})".format(gene["score"], kind)


def _gene_track(system):
    placements, lane_count = _gene_lanes(system["genes"])
    segments = []
    for gene, lane in placements:
        annotation = gene["annotation"]
        biological_role = "Cas-positive" if gene["is_cas"] else "non-Cas context"
        label = "{}; {}: {}-{} ({}); evidence {}".format(
            biological_role,
            annotation,
            gene["start"],
            gene["end"],
            gene["strand"],
            _score_text(gene),
        )
        segments.append(
            """\
                <span class="cas-gene cas-gene-{annotation_class} cas-strand-{strand_class}"
                      style="{interval_style};top:{top:.2f}rem"
                      title="{title}">
                  <span>{annotation}</span>
                </span>""".format(
                annotation_class="annotated" if gene["is_cas"] else "unknown",
                strand_class=(
                    "forward" if gene["strand"] == "+" else "reverse"
                ),
                interval_style=_interval_style(
                    gene["start"],
                    gene["end"],
                    system["start"],
                    system["end"],
                ),
                top=0.35 + lane * 2.35,
                title=_escape(label),
                annotation=_escape(annotation),
            )
        )
    return """\
              <div class="cas-track-viewport" tabindex="0">
                <div class="cas-gene-track"
                     style="height:{height:.2f}rem"
                     role="img"
                     aria-label="{aria_label}">
{segments}
                </div>
              </div>
              <p class="scroll-hint">
                Gene positions and widths are proportional to the
                1-based inclusive cassette interval. Arrowheads show coding
                strand; overlapping genes use separate lanes. Grey genes are
                non-Cas context retained by the provider.
              </p>""".format(
        height=0.7 + lane_count * 2.35,
        aria_label=_escape(
            "{} cassette from {} to {} with {} context genes".format(
                system["type"],
                system["start"],
                system["end"],
                len(system["genes"]),
            )
        ),
        segments="\n".join(segments),
    )


def _classification_evidence(system):
    classification = system["classification"]
    confidence = classification["confidence"]
    if confidence is None:
        confidence_text = "No classification confidence was reported."
    else:
        confidence_text = (
            "Provider-reported confidence: {:.3f} ({}) for method {}. "
            "Its meaning is method-specific; no cross-provider calibration "
            "is implied."
        ).format(
            confidence,
            classification["confidence_kind"],
            classification["method"],
        )

    legacy_rows = []
    for evidence in system.get("probabilities", [])[:5]:
        legacy_rows.append(
            """\
                  <div class="cas-evidence-row">
                    <span>{cas_type}</span>
                    <span class="cas-evidence-track" aria-hidden="true">
                      <i style="width:{width:.3f}%"></i>
                    </span>
                    <span class="numeric">{value:.3f}</span>
                  </div>""".format(
                cas_type=_escape(evidence["cas_type"]),
                width=evidence["probability"] * 100.0,
                value=evidence["probability"],
            )
        )
    legacy = ""
    if legacy_rows:
        legacy = """\
              <div class="cas-evidence"
                   aria-label="Legacy classifier evidence scores">
{rows}
              </div>
              <p class="muted">
                These legacy classifier scores are retained for compatibility;
                they are not calibrated biological probabilities.
              </p>""".format(rows="\n".join(legacy_rows))
    return """\
              <p class="cas-confidence">{confidence}</p>
{legacy}""".format(
        confidence=_escape(confidence_text),
        legacy=legacy,
    )


def _gene_table(system):
    rows = []
    for gene in system["genes"]:
        rows.append(
            """\
                    <tr>
                      <th scope="row"><code>{protein_id}</code></th>
                      <td>{is_cas}</td>
                      <td>{annotation}</td>
                      <td><code>{profile}</code></td>
                      <td class="numeric">{start:,}&ndash;{end:,}</td>
                      <td>{strand}</td>
                      <td>{score}</td>
                    </tr>""".format(
                protein_id=_escape(gene["protein_id"]),
                is_cas="Cas-positive" if gene["is_cas"] else "Context",
                annotation=_escape(gene["annotation"]),
                profile=(
                    "&mdash;" if gene["profile"] is None
                    else _escape(gene["profile"])
                ),
                start=gene["start"],
                end=gene["end"],
                strand=_escape(gene["strand"]),
                score=_escape(_score_text(gene)),
            )
        )
    return """\
              <details class="nested-details">
                <summary>Context genes ({count})</summary>
                <div class="component-table-wrap" tabindex="0">
                  <table class="cas-gene-table">
                    <thead>
                      <tr>
                        <th scope="col">Protein ID</th>
                        <th scope="col">Call</th>
                        <th scope="col">Operational annotation</th>
                        <th scope="col">Provider profile</th>
                        <th scope="col">Source coordinates</th>
                        <th scope="col">Strand</th>
                        <th scope="col">Evidence score</th>
                      </tr>
                    </thead>
                    <tbody>
{rows}
                    </tbody>
                  </table>
                </div>
              </details>""".format(
        count=len(system["genes"]),
        rows="\n".join(rows),
    )


def _genome_overview(model, systems):
    source_length = model["source"]["length"]
    array_markers = []
    for array in model.get("arrays", []):
        interval = array["source_interval"]
        array_markers.append(
            """\
              <a class="genome-marker genome-marker-array"
                 href="#{array_id}" style="{style}"
                 title="{title}" aria-label="{title}"></a>""".format(
                array_id=_escape(array["id"]),
                style=_interval_style(
                    interval["start"], interval["end"], 1, source_length
                ),
                title=_escape(
                    "CRISPR array {}: {}-{}".format(
                        array["id"], interval["start"], interval["end"]
                    )
                ),
            )
        )
    cas_markers = []
    for system in systems:
        system_id = "cas-cassette-{}".format(system["cassette_id"])
        cas_markers.append(
            """\
              <a class="genome-marker genome-marker-cas"
                 href="#{system_id}" style="{style}"
                 title="{title}" aria-label="{title}"></a>""".format(
                system_id=_escape(system_id),
                style=_interval_style(
                    system["start"], system["end"], 1, source_length
                ),
                title=_escape(
                    "{} cassette {}: {}-{}".format(
                        system["type"],
                        system["cassette_id"],
                        system["start"],
                        system["end"],
                    )
                ),
            )
        )
    return """\
        <div class="genome-context-legend" aria-hidden="true">
          <span><i class="legend-array-locus"></i> CRISPR array</span>
          <span><i class="legend-cas-locus"></i> Cas cassette</span>
        </div>
        <div class="genome-context-track" role="img"
             aria-label="Genome-wide positions of CRISPR arrays and Cas cassettes">
          <div class="genome-axis genome-axis-arrays">
            <span class="genome-axis-label">Arrays</span>
{array_markers}
          </div>
          <div class="genome-axis genome-axis-cas">
            <span class="genome-axis-label">Cas</span>
{cas_markers}
          </div>
        </div>
        <div class="genome-context-scale" aria-hidden="true">
          <span>1</span><span>{midpoint:,}</span><span>{length:,} bp</span>
        </div>""".format(
        array_markers="\n".join(array_markers),
        cas_markers="\n".join(cas_markers),
        midpoint=(source_length + 1) // 2,
        length=source_length,
    )


def _provider_label(method):
    provider = method.get("provider", "Cas annotation provider")
    version = method.get("provider_version")
    return "{} {}".format(provider, version) if version else provider


def _confidence_cell(system):
    classification = system["classification"]
    if classification["confidence"] is None:
        return "&mdash;"
    return "{:.3f}<br><small>{}</small>".format(
        classification["confidence"],
        _escape(classification["confidence_kind"]),
    )


def _experimental_badge(system):
    if not system["experimental"]:
        return ""
    return """\
<span class="badge badge-cas-experimental"
      title="Experimental provider classification">Experimental</span>"""


def _experimental_notice(system):
    if not system["experimental"]:
        return ""
    return """\
              <div class="cas-experimental-note" role="note">
                <strong>Experimental Type VII call.</strong> {reason}
              </div>""".format(
        reason=_escape(system["experimental_reason"]),
    )


def render_cas_section(model):
    cas, systems = _validate(model)
    if cas.get("status") == "not_requested":
        return """\
    <section aria-labelledby="cas-heading">
      <h2 id="cas-heading">CRISPR-Cas systems</h2>
      <p class="empty-state">
        Cas prediction was not requested. Run with <code>--cas True</code> to
        add cassette and Cas-protein annotations.
      </p>
    </section>"""

    method = cas["method"]
    provider_label = _provider_label(method)
    if not systems:
        return """\
    <section aria-labelledby="cas-heading">
      <h2 id="cas-heading">CRISPR-Cas systems</h2>
      <p class="empty-state">
        {provider} completed annotation, but no cassette passed its genome-level
        acceptance rules. {gene_count} Cas-positive proteins were retained in
        the canonical JSON and CSV outputs.
      </p>
    </section>""".format(
            provider=_escape(provider_label),
            gene_count=cas.get("gene_count", 0),
        )

    summary_rows = []
    cards = []
    accepted_positive_ids = set()
    for system_index, system in enumerate(systems):
        classification = system["classification"]
        cas_gene_count = sum(gene["is_cas"] for gene in system["genes"])
        accepted_positive_ids.update(
            gene["protein_id"] for gene in system["genes"] if gene["is_cas"]
        )
        non_cas_count = len(system["genes"]) - cas_gene_count
        system_id = "cas-cassette-{}".format(system["cassette_id"])
        summary_rows.append(
            """\
                <tr>
                  <th scope="row"><a href="#{system_id}">Cassette {cassette_id}</a></th>
                  <td>
                    <span class="badge badge-cas">{cas_type}</span>
                    {experimental_badge}
                  </td>
                  <td class="numeric">{start:,}&ndash;{end:,}</td>
                  <td class="numeric">{cas_gene_count}</td>
                  <td class="numeric">{context_count}</td>
                  <td class="numeric">{confidence}</td>
                </tr>""".format(
                system_id=_escape(system_id),
                cassette_id=system["cassette_id"],
                cas_type=_escape(system["type"]),
                experimental_badge=_experimental_badge(system),
                start=system["start"],
                end=system["end"],
                cas_gene_count=cas_gene_count,
                context_count=non_cas_count,
                confidence=_confidence_cell(system),
            )
        )
        cards.append(
            """\
        <article id="{system_id}" class="cas-card">
          <details{open_attr}>
            <summary>
              <span class="badge badge-cas">{cas_type}</span>
              {experimental_badge}
              <strong>Cassette {cassette_id}</strong>
              <span>{start:,}&ndash;{end:,}</span>
              <span>{cas_gene_count} Cas-positive genes</span>
            </summary>
            <div class="cas-body">
              <dl class="array-metadata">
                <div><dt>Classification</dt><dd>{cas_type}</dd></div>
                <div><dt>Source interval</dt>
                  <dd>{start:,}&ndash;{end:,} (1-based inclusive)</dd></div>
                <div><dt>Provider cassette ID</dt>
                  <dd><code>{provider_cassette_id}</code></dd></div>
                <div><dt>Provider</dt><dd>{provider}</dd></div>
                <div><dt>Classification method</dt>
                  <dd><code>{classification_method}</code></dd></div>
                <div><dt>Classification status</dt>
                  <dd>{classification_status}</dd></div>
                <div><dt>Non-Cas context genes</dt>
                  <dd>{non_cas_count}</dd></div>
              </dl>
{experimental_notice}
{gene_track}
{classification_evidence}
{gene_table}
            </div>
          </details>
        </article>""".format(
                system_id=_escape(system_id),
                open_attr=" open" if system_index == 0 else "",
                cas_type=_escape(system["type"]),
                experimental_badge=_experimental_badge(system),
                cassette_id=system["cassette_id"],
                start=system["start"],
                end=system["end"],
                cas_gene_count=cas_gene_count,
                provider_cassette_id=_escape(system["provider_cassette_id"]),
                provider=_escape(provider_label),
                classification_method=_escape(classification["method"]),
                classification_status=(
                    "Experimental" if system["experimental"]
                    else "Standard provider result"
                ),
                non_cas_count=non_cas_count,
                experimental_notice=_experimental_notice(system),
                gene_track=_gene_track(system),
                classification_evidence=_classification_evidence(system),
                gene_table=_gene_table(system),
            )
        )

    orphan_count = sum(
        gene["protein_id"] not in accepted_positive_ids for gene in cas["genes"]
    )
    return """\
    <section aria-labelledby="cas-heading">
      <h2 id="cas-heading">CRISPR-Cas systems</h2>
      <div class="cas-overview-card">
        <p>
          <strong>{system_count}</strong> Cas cassettes were classified by
          <strong>{provider}</strong>. The genome track relates their positions
          to reported CRISPR arrays. {orphan_count} additional Cas-positive
          proteins were not part of an accepted cassette.
        </p>
{overview}
        <div class="table-scroll" tabindex="0"
             aria-label="Scrollable Cas cassette summary">
          <table id="cas-summary">
            <caption>Predicted Cas cassettes in source-coordinate order</caption>
            <thead>
              <tr>
                <th scope="col">Cassette</th>
                <th scope="col">Classification</th>
                <th scope="col">Source coordinates</th>
                <th scope="col">Cas-positive genes</th>
                <th scope="col">Non-Cas context</th>
                <th scope="col">Provider confidence</th>
              </tr>
            </thead>
            <tbody>
{summary_rows}
            </tbody>
          </table>
        </div>
      </div>
{cards}
    </section>""".format(
        system_count=len(systems),
        provider=_escape(provider_label),
        orphan_count=orphan_count,
        overview=_genome_overview(model, systems),
        summary_rows="\n".join(summary_rows),
        cards="\n".join(cards),
    )


CAS_STYLES = r"""
    :root {
      --cas: #0f766e;
      --cas-soft: #ccfbf1;
      --cas-experimental: #9a3412;
      --cas-experimental-soft: #ffedd5;
      --cas-unknown: #64748b;
      --cas-track-bg: #e2e8f0;
    }
    @media (prefers-color-scheme: dark) {
      :root {
        --cas: #5eead4;
        --cas-soft: #134e4a;
        --cas-experimental: #fdba74;
        --cas-experimental-soft: #7c2d12;
        --cas-unknown: #64748b;
        --cas-track-bg: #475569;
      }
    }
    .cas-overview-card, .cas-card {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 0.5rem;
      margin-block: 1rem;
    }
    .cas-overview-card { padding: 1rem; }
    .cas-card { padding: 0; scroll-margin-top: 1rem; }
    .cas-card > details > summary {
      cursor: pointer;
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 0.6rem 1rem;
      padding: 1rem;
    }
    .cas-body { border-top: 1px solid var(--border); padding: 1rem; }
    .badge-cas { background: var(--cas-soft); color: var(--cas); }
    .badge-cas-experimental {
      background: var(--cas-experimental-soft);
      color: var(--cas-experimental);
      margin-inline-start: 0.3rem;
    }
    .cas-experimental-note {
      background: var(--cas-experimental-soft);
      border-inline-start: 0.3rem solid var(--cas-experimental);
      border-radius: 0.25rem;
      color: var(--text);
      margin-block: 1rem;
      padding: 0.75rem 1rem;
    }
    .genome-context-legend {
      display: flex;
      gap: 1.25rem;
      font-size: 0.9rem;
      margin-top: 1rem;
    }
    .genome-context-legend span {
      align-items: center; display: inline-flex; gap: 0.35rem;
    }
    .genome-context-legend i {
      border: 1px solid var(--text);
      display: inline-block;
      height: 0.75rem;
      width: 1.6rem;
    }
    .legend-array-locus { background: var(--repeat); }
    .legend-cas-locus { background: var(--cas); }
    .genome-context-track { margin: 0.45rem 0 0; padding-left: 4.5rem; }
    .genome-axis {
      background: var(--cas-track-bg);
      height: 1.35rem;
      margin: 0.65rem 0;
      position: relative;
    }
    .genome-axis-label {
      color: var(--muted);
      font-size: 0.82rem;
      font-weight: 700;
      left: -4.5rem;
      position: absolute;
      width: 4rem;
    }
    .genome-marker {
      display: block;
      height: 100%;
      min-width: 3px;
      position: absolute;
    }
    .genome-marker-array { background: var(--repeat); }
    .genome-marker-cas { background: var(--cas); }
    .genome-context-scale {
      display: flex;
      font-size: 0.78rem;
      justify-content: space-between;
      margin-left: 4.5rem;
    }
    .cas-track-viewport {
      border: 1px solid var(--border);
      border-radius: 0.35rem;
      margin-top: 1.1rem;
      min-width: 0;
      overflow-x: auto;
      padding: 0.35rem;
    }
    .cas-gene-track {
      background: linear-gradient(
        to bottom,
        transparent calc(50% - 1px),
        var(--border) calc(50% - 1px),
        var(--border) calc(50% + 1px),
        transparent calc(50% + 1px)
      );
      min-width: 42rem;
      position: relative;
    }
    .cas-gene {
      align-items: center;
      display: flex;
      height: 1.75rem;
      justify-content: center;
      min-width: 3px;
      overflow: hidden;
      padding-inline: 0.55rem;
      position: absolute;
      white-space: nowrap;
    }
    .cas-gene span {
      color: #fff;
      font-size: 0.72rem;
      font-weight: 700;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .cas-gene-annotated { background: var(--cas); }
    .cas-gene-unknown { background: var(--cas-unknown); }
    .cas-strand-forward {
      clip-path: polygon(0 0, calc(100% - 0.7rem) 0, 100% 50%,
        calc(100% - 0.7rem) 100%, 0 100%, 0.3rem 50%);
    }
    .cas-strand-reverse {
      clip-path: polygon(0 50%, 0.7rem 0, 100% 0,
        calc(100% - 0.3rem) 50%, 100% 100%, 0.7rem 100%);
    }
    .cas-confidence { margin-top: 1rem; }
    .cas-evidence {
      display: grid;
      gap: 0.35rem;
      margin-top: 1rem;
      max-width: 42rem;
    }
    .cas-evidence-row {
      align-items: center;
      display: grid;
      gap: 0.65rem;
      grid-template-columns: 7rem minmax(8rem, 1fr) 3.5rem;
    }
    .cas-evidence-track {
      background: var(--cas-track-bg);
      border-radius: 999px;
      height: 0.75rem;
      overflow: hidden;
    }
    .cas-evidence-track i {
      background: var(--cas);
      display: block;
      height: 100%;
    }
    @media print {
      .cas-card, .cas-overview-card { break-inside: avoid; }
      .cas-track-viewport { overflow: visible; }
    }
"""
