"""Portable, escaped evidence reports; report.json is the completion artifact."""

from __future__ import annotations

import csv
import html
import io
import json
from pathlib import Path
import shutil
import tempfile

from .io import atomic_json


def tsv_report(report: dict) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow([
        "id", "sequence", "status", "best_distance", "best_hit_count",
        "family_status", "historical_family", "family_candidates", "motif_status",
        "historical_motif", "motif_evidence", "reason",
    ])
    for row in report["results"]:
        family, motif = row["annotations"]["family"], row["annotations"]["motif"]
        values = [row["id"], row["sequence"], row["status"], row["best_distance"],
                  row["best_hit_count"], family["status"], family["value"],
                  ";".join(family["values"]), motif["status"], motif["value"],
                  ";".join(motif["values"]), row["reason"]]
        # Prevent spreadsheet formula execution when users open TSV exports.
        writer.writerow(["" if v is None else ("'" + v if isinstance(v, str)
                         and v.startswith(("=", "+", "-", "@")) else v) for v in values])
    return stream.getvalue()


def html_report(report: dict) -> str:
    def escape(value):
        return html.escape(str(value), quote=True)
    rows = []
    for row in report["results"]:
        family = row["annotations"]["family"]
        motif = row["annotations"]["motif"]
        hits = []
        for hit in row["hits"]:
            hits.append(
                f'<li><strong>{escape(hit["reference_id"])}</strong> · '
                f'{escape(hit["orientation"])} · {hit["distance"]} edits'
                f'<pre>{escape(hit["query_aligned"])}\n'
                f'{escape(hit["reference_aligned"])}</pre></li>'
            )
        family_text = ", ".join(family["values"]) or "Missing"
        motif_text = ", ".join(motif["values"]) or "Missing"
        rows.append(
            f'<tr><td><strong>{escape(row["id"])}</strong><code>{escape(row["sequence"])}</code></td>'
            f'<td><span class="badge {escape(row["status"])}">{escape(row["status"])}</span>'
            f'<small>{escape(row["reason"] or "")}</small></td>'
            f'<td>{escape(family_text)}<small>{escape(family["status"])}</small></td>'
            f'<td>{escape(motif_text)}<small>{escape(motif["status"])}</small></td>'
            f'<td>{escape(row["best_distance"] if row["best_distance"] is not None else "—")}</td>'
            f'<td><details><summary>{row["reported_hit_count"]} of '
            f'{row["best_hit_count"]} best hits</summary><ul>{"".join(hits)}</ul></details></td></tr>'
        )
    counts = {status: sum(r["status"] == status for r in report["results"])
              for status in ("exact", "near", "no_hit", "unsupported")}
    cards = "".join(f'<div><b>{count}</b><span>{status.replace("_", " ")}</span></div>'
                    for status, count in counts.items())
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'">
<title>CRISPRmap-v2 · Repeat evidence</title><style>
:root{{color-scheme:light}}*{{box-sizing:border-box}}body{{margin:0;background:#f6f4ee;color:#22332f;
font:16px/1.5 system-ui,sans-serif}}main{{max-width:1500px;margin:auto;padding:48px 32px}}
h1{{font:normal 42px/1.15 Georgia,serif;margin:10px 0}}.eyebrow{{text-transform:uppercase;
letter-spacing:.14em;font-size:12px;color:#53675e}}.intro{{max-width:900px}}.metrics{{display:flex;
gap:12px;margin:32px 0}}.metrics div{{background:white;border-top:3px solid #437e6c;padding:15px;
min-width:140px}}.metrics b{{font:34px Georgia,serif;display:block}}.metrics span{{font-size:13px}}
.table{{overflow:auto;background:white}}table{{border-collapse:collapse;width:100%;font-size:14px}}
td,th{{text-align:left;border-bottom:1px solid #dcded6;padding:15px;vertical-align:top}}
th{{font-size:11px;text-transform:uppercase;letter-spacing:.07em;color:#53675e}}code,small{{display:block}}
code{{font-size:11px;overflow-wrap:anywhere;margin-top:6px}}small{{font-size:12px;color:#596760;
max-width:220px;margin-top:7px}}.badge{{border-radius:4px;padding:3px 7px;background:#e9ece7}}
.exact{{background:#d6ebe0}}.near{{background:#faedc8}}summary{{cursor:pointer}}ul{{padding-left:16px}}
pre{{font-size:11px;line-height:1.4}}footer{{margin-top:25px;font-size:12px;overflow-wrap:anywhere}}
@media(max-width:650px){{main{{padding:25px 16px}}.metrics{{flex-wrap:wrap}}h1{{font-size:32px}}}}
</style></head><body><main><div class="eyebrow">Prokaryotic repeat conservation · 2026 implementation</div>
<h1>CRISPRmap-v2</h1><p class="intro">Historical annotation evidence and full-length sequence matches.
Near matches are uncalibrated candidates. Reference-relative orientation does not establish transcription
direction. Cas subtype and host are not predicted. Motif evidence on near matches describes reference
records only; no structural motif is assigned.</p>
<p>Reference: <strong>{escape(report["reference"]["namespace"])}</strong></p>
<section class="metrics" aria-label="Result counts">{cards}</section>
<div class="table"><table><thead><tr><th>Repeat</th><th>Evidence</th><th>Sequence family</th>
<th>Historical motif evidence</th><th>Edits</th><th>Alignments</th></tr></thead>
<tbody>{"".join(rows)}</tbody></table></div><footer>Reference content SHA-256:
{escape(report["reference"]["content_sha256"])}<br>Input SHA-256:
{escape(report.get("input", {}).get("sha256", "unavailable"))}<br>
JSON contains the full method, limits, provenance and missing-value states.</footer></main></body></html>'''


def write_report(report: dict, destination: Path) -> None:
    """Refuse existing paths. Stage all artifacts, then write completion JSON last."""
    destination = Path(destination)
    if any(p.is_symlink() for p in (destination.absolute(), *destination.absolute().parents)):
        raise ValueError("output paths must not traverse symlinks")
    if destination.exists() or destination.is_symlink():
        raise ValueError("output path already exists; choose a new directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".crisprmap-", dir=destination.parent) as tmp:
        stage = Path(tmp)
        artifacts = {
            "report.tsv": tsv_report(report), "report.html": html_report(report),
            "report.json": json.dumps(report, indent=2, allow_nan=False) + "\n",
        }
        if sum(len(value.encode()) for value in artifacts.values()) > 128 * 1024 * 1024:
            raise ValueError("report exceeds 128 MiB output limit")
        for name, value in artifacts.items():
            (stage / name).write_text(value, encoding="utf-8")
        destination.mkdir()  # Exclusive reservation; never overwrite a raced-in directory.
        created = []
        try:
            for name in ("report.tsv", "report.html"):
                with (destination / name).open("xb") as handle:
                    created.append(destination / name)
                    with (stage / name).open("rb") as source:
                        shutil.copyfileobj(source, handle)
            atomic_json(destination / "report.json", report)
        except BaseException:
            for path in created:
                path.unlink(missing_ok=True)
            destination.rmdir()
            raise
