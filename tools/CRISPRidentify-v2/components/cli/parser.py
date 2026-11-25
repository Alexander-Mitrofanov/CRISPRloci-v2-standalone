"""Build the detector CLI without parsing argv or initializing the pipeline."""

import argparse

from components.cli.presets import BALANCED_DEFAULTS, PRESET_CONFIGS
from components.cli.validation import (
    parse_boolean,
    positive_integer,
    readable_input_directory,
    readable_input_file,
)
from components.paths import read_release_version


def build_parser():
    """Return a fresh parser; callers own argument parsing and execution."""
    parser = argparse.ArgumentParser(
        description=(
            "Detect and classify CRISPR arrays, with optional CasAndra Cas annotation"
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version="CRISPRidentify {}".format(read_release_version()),
    )
    parser.add_argument(
        "--preset",
        choices=tuple(PRESET_CONFIGS),
        default="balanced",
        help=(
            "atomic analysis preset: balanced keeps the default high-recall fast "
            "analysis, region-screen is the validated high-throughput region-"
            "recall screen, and full selects historical full search "
            "(default: balanced)"
        ),
    )
    parser.add_argument(
        "--print-effective-config",
        action="store_true",
        help="print resolved scientific configuration as JSON and exit",
    )
    parser.add_argument("--repeat-template-fasta", "--repeat_template_fasta",
                        help="optional small FASTA panel of 15..63-base literal repeat templates")
    parser.add_argument("--boundary-consensuses", metavar="none|SEQ[,SEQ...]",
                        help="override bundled repeat boundaries on newly extended arrays, or use none for source-only resolution")
    parser.add_argument("--_repeat_template_sha256", help=argparse.SUPPRESS)
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "--quiet",
        action="store_true",
        help="suppress routine progress output",
    )
    verbosity_group.add_argument(
        "--verbose",
        action="store_true",
        help="show detailed pipeline stages and diagnostic telemetry",
    )
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument(
        "--input_folder",
        type=readable_input_directory,
        default=None,
        help="input folder (default: None)",
    )

    input_group.add_argument(
        "--file",
        type=readable_input_file,
        default=None,
        help="input file (default: None)",
    )

    input_group.add_argument(
        "--input_folder_multifasta",
        type=readable_input_directory,
        default=None,
        help="input folder of multifasta (default: None)",
    )

    parser.add_argument(
        "--result_folder",
        type=str,
        default="Results",
        help="folder with the result (default: Results)",
    )

    parser.add_argument(
        "--pickle_report",
        type=str,
        default="",
        help="pickled report file (default: None)",
    )

    parser.add_argument(
        "--json_report",
        type=str,
        default="",
        help=(
            "legacy JSON export folder; canonical report.json is "
            "always written in each result directory (default: None)"
        ),
    )

    parser.add_argument(
        "--fasta_report",
        type=parse_boolean,
        default=False,
        help="fasta report file (default: False)",
    )

    parser.add_argument(
        "--html_report",
        type=parse_boolean,
        default=False,
        help=(
            "write optional offline report.html; canonical "
            "source-reconstructed report.json is always written "
            "(default: False)"
        ),
    )

    parser.add_argument(
        "--output_profile",
        choices=("standard", "machine"),
        default="standard",
        help=(
            "standard writes legacy text/CSV plus canonical JSON; machine writes "
            "only the transactional canonical report for large-scale ingestion "
            "(default: standard)"
        ),
    )

    parser.add_argument(
        "--biological_loci",
        type=parse_boolean,
        default=False,
        help=(
            "add optional heuristic biological-locus presentation groups to the "
            "canonical JSON and optional HTML report (default: False)"
        ),
    )

    parser.add_argument(
        "--sequence_topology_manifest",
        type=readable_input_file,
        default=None,
        help=(
            "strict per-source topology manifest backed by frozen INSDC "
            "flatfiles; sources absent from the manifest remain unknown "
            "(default: None)"
        ),
    )

    parser.add_argument(
        "--strand",
        type=parse_boolean,
        default=True,
        help="CRISPR array orientation prediction (default: True)",
    )

    parser.add_argument(
        "--is_element",
        type=parse_boolean,
        default=False,
        help="IS-element annotation; requires the is-elements environment (default: False)",
    )

    parser.add_argument(
        "--cas",
        type=parse_boolean,
        default=False,
        help="CasAndra CRISPR-Cas system annotation (default: False)",
    )

    parser.add_argument(
        "--cas_gene_mode",
        choices=("auto", "single", "meta"),
        default="auto",
        help="CasAndra/Pyrodigal gene mode (default: auto)",
    )
    parser.add_argument(
        "--cas_timeout_seconds",
        type=positive_integer,
        default=3600,
        help=(
            "maximum runtime for each CasAndra genome annotation before its "
            "process group is terminated (default: 3600)"
        ),
    )
    parser.add_argument(
        "--cas_completeness",
        choices=("complete", "partial"),
        default=None,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--parallel",
        type=parse_boolean,
        default=True,
        help="parallel computations (default: True)",
    )

    parser.add_argument(
        "--cpu", type=str, default="ALL", help="global CPU budget (default: ALL)"
    )

    parser.add_argument(
        "--folder_jobs",
        type=int,
        default=1,
        help="persistent, affinity-isolated folder lanes (default: 1)",
    )

    parser.add_argument("--_folder_lane", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--_folder_queue", type=str, default=None, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--_folder_lane_id", type=int, default=None, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--_folder_cpu_ids", type=str, default=None, help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--fast_run",
        type=parse_boolean,
        default=True,
        help=(
            "high-recall fast path; pass --fast_run False for the "
            "historical full search (default: True)"
        ),
    )

    parser.add_argument(
        "--fast_run_seed_profile",
        choices=("balanced", "legacy", "sensitive"),
        default="balanced",
        help=(
            "Vmatch -l seed/pair-gap settings in fast mode: balanced 23 -10 100; "
            "legacy 23 25 60; sensitive 21 -10 100. Full mode always uses "
            "21 -10 100. These are not guaranteed biological spacer bounds "
            "(default: balanced)"
        ),
    )

    parser.add_argument(
        "--fast_run_rescue",
        type=parse_boolean,
        default=True,
        help=(
            "run the bounded classifier-routed start/end rescue and independent "
            "source-component selection in fast mode (default: True)"
        ),
    )

    parser.add_argument(
        "--fast_run_rescue_seed_scope",
        choices=("lineage", "locus", "family", "cluster"),
        default="cluster",
        help=(
            "seed domain for bounded fast rescue: lineage trims only the raw "
            "seeds that generated routed rows; locus adds Vmatch seeds observed "
            "inside the routed source interval; family adds the three nearest "
            "same-cluster seeds; cluster is the broader diagnostic control "
            "(default: cluster)"
        ),
    )

    parser.add_argument(
        "--fast_run_rescue_trim_profile",
        choices=("full", "total3", "one_sided"),
        default="full",
        help=(
            "terminal trim hypotheses per routed seed: full uses all 15 nonzero "
            "pairs in 0..3 x 0..3; total3 caps their summed trim; one_sided "
            "changes one edge only (default: full)"
        ),
    )

    parser.add_argument(
        "--fast_run_phase_repair",
        type=parse_boolean,
        default=True,
        help=(
            "batch-score the hard-capped source-backed phase grid and commit only "
            "Possible-to-Bona-fide promotions or new accepted calls in unresolved "
            "components (default: True)"
        ),
    )

    parser.add_argument(
        "--fast_two_repeat_terminal_repair",
        type=parse_boolean,
        default=False,
        help=(
            "opt in to a source-only, hard-capped terminal completion of narrow "
            "two-repeat accepted-only candidates; every proposal must pass the "
            "ordinary structural filter and model (default: False)"
        ),
    )

    parser.add_argument(
        "--fast_accepted_family_rescue",
        type=parse_boolean,
        default=False,
        help=(
            "opt in to two hard-capped source-only accepted-family rescues: "
            "strict family-guided phase/terminal completion and an exact-template "
            "three-repeat satellite; disabled by default pending an independent "
            "biological holdout (default: False)"
        ),
    )

    parser.add_argument(
        "--split_competitive_components",
        type=parse_boolean,
        default=True,
        help=(
            "select independent accepted winners for source-disjoint scored "
            "components in fast mode; required by the default fast-run rescue "
            "(default: True)"
        ),
    )

    parser.add_argument(
        "--degenerated",
        type=parse_boolean,
        default=True,
        help=(
            "enable degenerated/terminal-repeat enhancement; "
            "False also disables architecture rescue (default: True)"
        ),
    )

    parser.add_argument(
        "--min_len_rep",
        type=int,
        default=21,
        help="later structural minimum average length of the repeats (default: 21)",
    )

    parser.add_argument(
        "--max_len_rep",
        type=int,
        default=55,
        help="later structural maximum average length of the repeats (default: 55)",
    )

    parser.add_argument(
        "--min_len_spacer",
        type=int,
        default=18,
        help="later structural minimum average length of spacers (default: 18)",
    )

    parser.add_argument(
        "--max_len_spacer",
        type=int,
        default=78,
        help="later structural maximum average length of spacers (default: 78)",
    )

    parser.add_argument(
        "--min_repeats", type=int, default=3, help="min number of repeats (default: 3)"
    )

    parser.add_argument(
        "--enhancement_max_min",
        type=parse_boolean,
        default=True,
        help="enhancement with filter (default: True)",
    )

    parser.add_argument(
        "--enhancement_start_end",
        type=parse_boolean,
        default=True,
        help="enhancement with start end omitting (default: True)",
    )

    parser.add_argument(
        "--enhancement_intermediate",
        type=parse_boolean,
        default=True,
        help="intermediate substring enhancement (default: True)",
    )

    parser.add_argument(
        "--max_identical_spacers",
        type=int,
        default=4,
        help="maximum number of identical spacers in the array (default: 4)",
    )

    parser.add_argument(
        "--max_identical_cluster_spacers",
        type=int,
        default=3,
        help="maximum number of consecutive identical spacers in the array (default: 3)",
    )

    parser.add_argument(
        "--margin_degenerated",
        type=int,
        default=30,
        help="maximum length of the spacer margin for the degenerated search (default: 30)",
    )

    parser.add_argument(
        "--max_edit_distance_enhanced",
        type=int,
        default=6,
        help="ordinary enhancement edit limit (default: 6). With --degenerated, the anchored second pass allows up to four additional edits, capped at one third of repeat length",
    )

    parser.add_argument(
        "--accepted_only",
        type=parse_boolean,
        default=False,
        help="retain and report only Bona-fide/Possible arrays (default: False)",
    )

    parser.add_argument(
        "--retain_candidate_evidence",
        type=parse_boolean,
        default=False,
        help=(
            "write all five final candidate categories to a compressed "
            "research sidecar; does not add accepted arrays (default: False)"
        ),
    )

    parser.add_argument(
        "--reconstruction_policy",
        "--reconstruction-policy",
        choices=("legacy", "conservative"),
        default=BALANCED_DEFAULTS["reconstruction_policy"],
        help=(
            "separate Bona-fide core confidence from array "
            "reconstruction; conservative validates compatible Alternatives and "
            "preserves ambiguous phases (default: conservative for balanced; "
            "legacy for full, region-screen, or architecture rescue)"
        ),
    )

    parser.add_argument(
        "--architecture_rescue",
        type=parse_boolean,
        default=False,
        help=(
            "enable source-backed terminal/shared-repeat assembly, "
            "coding-context adjudication, and conservative array rescue; "
            "requires --degenerated True (default: False)"
        ),
    )

    parser.add_argument(
        "--fuzzy_search_mode",
        choices=("full", "family", "native"),
        default="native",
        help="fuzzy repeat search strategy (default: native)",
    )

    parser.add_argument(
        "--fuzzy_search_flank",
        type=int,
        default=100,
        help="bases added around each source discovery window (default: 100)",
    )

    parser.add_argument(
        "--fuzzy_family_padding",
        type=int,
        default=15,
        help="bases added around family-anchor search windows (default: 15)",
    )

    parser.add_argument(
        "--fuzzy_boundary_paths",
        type=parse_boolean,
        choices=(True, False),
        default=True,
        help=("build fuzzy forced-boundary path alternatives (default: True)"),
    )
    return parser
