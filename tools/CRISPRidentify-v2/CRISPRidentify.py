import argparse
import contextlib
import json
import math
import shutil
import signal
import subprocess
import sys
import tempfile
import warnings
import os


# Keep all remaining set/dictionary traversal deterministic.  The model-facing
# EDeN hashes themselves are reproduced explicitly by components.legacy_hash,
# so they no longer depend on the interpreter minor version.
_HASH_REEXEC_MARKER = "_CRISPRIDENTIFY_HASH_REEXEC"
if __name__ == "__main__" and os.environ.get("PYTHONHASHSEED") != "0":
    if os.environ.get(_HASH_REEXEC_MARKER) == "1":
        raise RuntimeError(
            "unable to start CRISPRidentify with PYTHONHASHSEED=0"
        )
    compatible_environment = os.environ.copy()
    compatible_environment["PYTHONHASHSEED"] = "0"
    compatible_environment[_HASH_REEXEC_MARKER] = "1"
    os.execve(
        sys.executable,
        [sys.executable] + sys.argv,
        compatible_environment,
    )

from pathlib import Path
from os import listdir
from os.path import isfile, join
from time import time

from components.cli_utils import optional_analysis_errors
from components.cli_utils import parse_boolean
from components.cli_utils import positive_integer
from components.cli_utils import readable_input_directory
from components.cli_utils import readable_input_file
from components.folder_scheduler import get_allowed_cpu_ids
from components.folder_scheduler import partition_cpu_lanes
from components.folder_scheduler import pipeline_result_key
from components.folder_scheduler import resolve_cpu_ids
from components.folder_scheduler import rewrite_lane_argv
from components.folder_scheduler import validate_unique_result_keys
from components.runtime_resources import configure_native_threads
from components.runtime_resources import native_thread_environment
from components.runtime_presets import BALANCED_DEFAULTS
from components.runtime_presets import PRESET_CONFIGS
from components.runtime_presets import explicit_option_destinations
from components.runtime_presets import preset_provenance
from components.runtime_presets import resolve_preset
from components.legacy_hash import py37_hash

warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

if py37_hash(("EDeN", 42)) != -1821860980875793120:
    raise RuntimeError("the portable EDeN hash compatibility check failed")

FLAG_DEVELOPER_MODE = False


class InputValidationError(ValueError):
    """An input folder cannot be consumed safely without user action."""

parser = argparse.ArgumentParser(
    description=(
        'Detect and classify CRISPR arrays, with optional CasAndra Cas '
        'annotation'
    )
)
try:
    _release_version = Path(__file__).with_name("VERSION").read_text(
        encoding="utf-8"
    ).strip()
except OSError:
    _release_version = "unknown"
parser.add_argument(
    "--version",
    action="version",
    version="CRISPRidentify {}".format(_release_version),
)
parser.add_argument(
    '--preset',
    choices=tuple(PRESET_CONFIGS),
    default='balanced',
    help=(
        'atomic analysis preset: balanced keeps the default high-recall fast '
        'analysis, region-screen is the validated high-throughput region-'
        'recall screen, and full selects historical full search '
        '(default: balanced)'
    ),
)
parser.add_argument(
    '--print-effective-config',
    action='store_true',
    help='print resolved scientific configuration as JSON and exit',
)
verbosity_group = parser.add_mutually_exclusive_group()
verbosity_group.add_argument(
    '--quiet',
    action='store_true',
    help='suppress routine progress output',
)
verbosity_group.add_argument(
    '--verbose',
    action='store_true',
    help='show detailed pipeline stages and diagnostic telemetry',
)
input_group = parser.add_mutually_exclusive_group(required=False)
input_group.add_argument('--input_folder', type=readable_input_directory,
                    default=None,
                    help='input folder (default: None)')

input_group.add_argument('--file', type=readable_input_file, default=None,
                    help='input file (default: None)')

input_group.add_argument('--input_folder_multifasta',
                    type=readable_input_directory, default=None,
                    help='input folder of multifasta (default: None)')

parser.add_argument('--result_folder', type=str, default="Results",
                    help='folder with the result (default: Results)')

parser.add_argument('--pickle_report', type=str, default='',
                    help='pickled report file (default: None)')

parser.add_argument('--json_report', type=str, default='',
                    help=(
                        'legacy JSON export folder; canonical report.json is '
                        'always written in each result directory (default: None)'
                    ))

parser.add_argument('--fasta_report', type=parse_boolean, default=False,
                    help='fasta report file (default: False)')

parser.add_argument('--html_report', type=parse_boolean, default=False,
                    help=(
                        'write optional offline report.html; canonical '
                        'source-reconstructed report.json is always written '
                        '(default: False)'
                    ))

parser.add_argument(
    '--output_profile',
    choices=('standard', 'machine'),
    default='standard',
    help=(
        'standard writes legacy text/CSV plus canonical JSON; machine writes '
        'only the transactional canonical report for large-scale ingestion '
        '(default: standard)'
    ),
)

parser.add_argument(
    '--biological_loci',
    type=parse_boolean,
    default=False,
    help=(
        'add optional heuristic biological-locus presentation groups to the '
        'canonical JSON and optional HTML report (default: False)'
    ),
)

parser.add_argument(
    '--sequence_topology_manifest',
    type=readable_input_file,
    default=None,
    help=(
        'strict per-source topology manifest backed by frozen INSDC '
        'flatfiles; sources absent from the manifest remain unknown '
        '(default: None)'
    ),
)

parser.add_argument('--strand', type=parse_boolean, default=True,
                    help='CRISPR array orientation prediction (default: True)')

parser.add_argument('--is_element', type=parse_boolean, default=False,
                    help='IS-element annotation; requires the is-elements environment (default: False)')

parser.add_argument('--cas', type=parse_boolean, default=False,
                    help='CasAndra CRISPR-Cas system annotation (default: False)')

parser.add_argument('--cas_gene_mode', choices=('auto', 'single', 'meta'),
                    default='auto',
                    help='CasAndra/Pyrodigal gene mode (default: auto)')
parser.add_argument(
    '--cas_timeout_seconds',
    type=positive_integer,
    default=3600,
    help=(
        'maximum runtime for each CasAndra genome annotation before its '
        'process group is terminated (default: 3600)'
    ),
)
parser.add_argument('--cas_completeness', choices=('complete', 'partial'),
                    default=None, help=argparse.SUPPRESS)

parser.add_argument('--parallel', type=parse_boolean, default=True,
                    help='parallel computations (default: True)')

parser.add_argument('--cpu', type=str, default="ALL",
                    help='global CPU budget (default: ALL)')

parser.add_argument('--folder_jobs', type=int, default=1,
                    help='persistent, affinity-isolated folder lanes (default: 1)')

parser.add_argument('--_folder_lane', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--_folder_queue', type=str, default=None, help=argparse.SUPPRESS)
parser.add_argument('--_folder_lane_id', type=int, default=None, help=argparse.SUPPRESS)
parser.add_argument('--_folder_cpu_ids', type=str, default=None, help=argparse.SUPPRESS)

parser.add_argument('--fast_run', type=parse_boolean, default=True,
                    help=(
                        'high-recall fast path; pass --fast_run False for the '
                        'historical full search (default: True)'
                    ))

parser.add_argument(
    '--fast_run_seed_profile',
    choices=('balanced', 'legacy', 'sensitive'),
    default='balanced',
    help=(
        'Vmatch seed/gap profile used by fast_run: balanced keeps the fast '
        '23-nt seed while restoring the accepted spacer domain; legacy uses '
        'the historical narrow fast search; sensitive uses normal-mode seeds '
        '(default: balanced)'
    ),
)

parser.add_argument(
    '--fast_run_rescue',
    type=parse_boolean,
    default=True,
    help=(
        'run the bounded classifier-routed start/end rescue and independent '
        'source-component selection in fast mode (default: True)'
    ),
)

parser.add_argument(
    '--fast_run_rescue_seed_scope',
    choices=('lineage', 'locus', 'family', 'cluster'),
    default='cluster',
    help=(
        'seed domain for bounded fast rescue: lineage trims only the raw '
        'seeds that generated routed rows; locus adds Vmatch seeds observed '
        'inside the routed source interval; family adds the three nearest '
        'same-cluster seeds; cluster is the broader diagnostic control '
        '(default: cluster)'
    ),
)

parser.add_argument(
    '--fast_run_rescue_trim_profile',
    choices=('full', 'total3', 'one_sided'),
    default='full',
    help=(
        'terminal trim hypotheses per routed seed: full uses all 15 nonzero '
        'pairs in 0..3 x 0..3; total3 caps their summed trim; one_sided '
        'changes one edge only (default: full)'
    ),
)

parser.add_argument(
    '--fast_run_phase_repair',
    type=parse_boolean,
    default=True,
    help=(
        'batch-score the hard-capped source-backed phase grid and commit only '
        'Possible-to-Bona-fide promotions or new accepted calls in unresolved '
        'components (default: True)'
    ),
)

parser.add_argument(
    '--fast_two_repeat_terminal_repair',
    type=parse_boolean,
    default=False,
    help=(
        'opt in to a source-only, hard-capped terminal completion of narrow '
        'two-repeat accepted-only candidates; every proposal must pass the '
        'ordinary structural filter and model (default: False)'
    ),
)

parser.add_argument(
    '--fast_accepted_family_rescue',
    type=parse_boolean,
    default=False,
    help=(
        'opt in to two hard-capped source-only accepted-family rescues: '
        'strict family-guided phase/terminal completion and an exact-template '
        'three-repeat satellite; disabled by default pending an independent '
        'biological holdout (default: False)'
    ),
)

parser.add_argument(
    '--split_competitive_components',
    type=parse_boolean,
    default=True,
    help=(
        'select independent accepted winners for source-disjoint scored '
        'components in fast mode; required by the default fast-run rescue '
        '(default: True)'
    ),
)

parser.add_argument('--degenerated', type=parse_boolean, default=True,
                    help=(
                        'enable degenerated/terminal-repeat enhancement; '
                        'False also disables architecture rescue (default: True)'
                    ))

parser.add_argument('--min_len_rep', type=int, default=21,
                    help='min avg. length of the repeats (default: 21)')

parser.add_argument('--max_len_rep', type=int, default=55,
                    help='max avg. length of the repeats (default: 55)')

parser.add_argument('--min_len_spacer', type=int, default=18,
                    help='min avg. length of spacers (default: 18)')

parser.add_argument('--max_len_spacer', type=int, default=78,
                    help='max avg. length of spacers (default: 78)')

parser.add_argument('--min_repeats', type=int, default=3,
                    help='min number of repeats (default: 3)')

parser.add_argument('--enhancement_max_min', type=parse_boolean, default=True,
                    help='enhancement with filter (default: True)')

parser.add_argument('--enhancement_start_end', type=parse_boolean, default=True,
                    help='enhancement with start end omitting (default: True)')

parser.add_argument(
    '--enhancement_intermediate',
    type=parse_boolean,
    default=True,
    help='intermediate substring enhancement (default: True)',
)

parser.add_argument('--max_identical_spacers', type=int, default=4,
                    help='maximum number of identical spacers in the array (default: 4)')

parser.add_argument('--max_identical_cluster_spacers', type=int, default=3,
                    help='maximum number of consecutive identical spacers in the array (default: 3)')

parser.add_argument('--margin_degenerated', type=int, default=30,
                    help='maximum length of the spacer margin for the degenerated search (default: 30)', )

parser.add_argument('--max_edit_distance_enhanced', type=int, default=6,
                    help='maximum edit distance for the evaluated array enhancement (default: 6)')

parser.add_argument('--accepted_only', type=parse_boolean, default=False,
                    help='retain and report only Bona-fide/Possible arrays (default: False)')

parser.add_argument('--retain_candidate_evidence', type=parse_boolean, default=False,
                    help=('write all five final candidate categories to a compressed '
                          'research sidecar; does not add accepted arrays (default: False)'))

parser.add_argument(
    '--reconstruction_policy', '--reconstruction-policy',
    choices=('legacy', 'conservative'),
    default=BALANCED_DEFAULTS['reconstruction_policy'],
    help=('separate Bona-fide core confidence from array '
          'reconstruction; conservative validates compatible Alternatives and '
          'preserves ambiguous phases (default: conservative for balanced; '
          'legacy for full, region-screen, or architecture rescue)'),
)

parser.add_argument(
    '--architecture_rescue',
    type=parse_boolean,
    default=False,
    help=(
        'enable source-backed terminal/shared-repeat assembly, '
        'coding-context adjudication, and conservative array rescue; '
        'requires --degenerated True (default: False)'
    ),
)

parser.add_argument('--fuzzy_search_mode', choices=('full', 'family', 'native'), default='native',
                    help='fuzzy repeat search strategy (default: native)')

parser.add_argument('--fuzzy_search_flank', type=int, default=100,
                    help='bases added around each Vmatch cluster (default: 100)')

parser.add_argument('--fuzzy_family_padding', type=int, default=15,
                    help='bases added around family-anchor search windows (default: 15)')

parser.add_argument(
    '--fuzzy_boundary_paths',
    type=parse_boolean,
    choices=(True, False),
    default=True,
    help=(
        'build fuzzy forced-boundary path alternatives '
        '(default: True)'
    ),
)


# Spawned matcher workers import this module without starting a second CLI
# or creating/changing into a second temporary working directory.
if __name__ == "__main__":
    script_absolute_path = os.path.dirname(os.path.abspath(__file__))
    work_directory = os.getcwd()
    pid = os.getpid()

    args = parser.parse_args()

    # Resolve lane affinity and native limits before Pipeline imports NumPy. This
    # covers the direct Python entry point as well as shell and batch launchers.
    if args._folder_cpu_ids:
        try:
            folder_cpu_ids = {
                int(cpu_id) for cpu_id in args._folder_cpu_ids.split(",")
            }
        except ValueError:
            parser.error("invalid internal folder CPU list")
        if not folder_cpu_ids or min(folder_cpu_ids) < 0:
            parser.error("invalid internal folder CPU list")
        if hasattr(os, "sched_setaffinity"):
            try:
                os.sched_setaffinity(0, folder_cpu_ids)
            except OSError as error:
                parser.error("unable to set folder lane CPU affinity: {}".format(error))
    try:
        runtime_cpu_ids = resolve_cpu_ids(
            args.cpu, allowed_cpus=get_allowed_cpu_ids()
        )
    except ValueError as error:
        parser.error(str(error))
    flag_cpu = str(len(runtime_cpu_ids))
    configure_native_threads(len(runtime_cpu_ids))

    from components.pipeline import Pipeline
    from components.module_evaluation import resolve_competitive_component_split
    from components.sequence_topology import (
        SequenceTopologyError,
        load_sequence_topology_manifest,
    )
    from components import literal_fuzzy
    from components.components_output_maker import CompleteFastaOutputMaker
    from components.components_output_maker import CompleteFolderSummaryMaker
    from components.components_output_maker import CompleteJsonOutputMaker
    from components.components_helpers import multiline_fasta_handle_python
    from components.components_helpers import folder_of_multifasta_handle

    explicit_destinations = explicit_option_destinations(parser, sys.argv[1:])
    runtime_values = {
        destination: getattr(args, destination)
        for destination in BALANCED_DEFAULTS
    }
    try:
        resolved_runtime_values = resolve_preset(
            args.preset,
            runtime_values,
            explicit_destinations=explicit_destinations,
        )
    except ValueError as error:
        parser.error(str(error))
    for destination, value in resolved_runtime_values.items():
        setattr(args, destination, value)

    try:
        split_competitive_components = resolve_competitive_component_split(
            fast_run=args.fast_run,
            requested=args.split_competitive_components,
            fast_run_rescue=args.fast_run_rescue,
        )
    except (TypeError, ValueError) as error:
        parser.error(str(error))

    if args.output_profile == 'machine' and (
            args.html_report
            or args.fasta_report
            or args.pickle_report
            or args.json_report):
        parser.error(
            '--output_profile machine is incompatible with HTML, FASTA, pickle, '
            'and legacy JSON exports'
        )

    cas_gene_mode = args.cas_gene_mode
    if args.cas_completeness is not None:
        legacy_gene_mode = {
            "complete": "single",
            "partial": "meta",
        }[args.cas_completeness]
        if cas_gene_mode != "auto" and cas_gene_mode != legacy_gene_mode:
            parser.error(
                "--cas_gene_mode conflicts with deprecated --cas_completeness"
            )
        cas_gene_mode = legacy_gene_mode

    if not args.print_effective_config and not (
            args.input_folder
            or args.file
            or args.input_folder_multifasta):
        parser.error(
            "one input source is required unless --print-effective-config is used"
        )

    complete_path_folder = (args.input_folder)
    if complete_path_folder:
        complete_path_folder = Path(complete_path_folder).absolute()

    complete_path_file = args.file
    if complete_path_file:
        complete_path_file = Path(complete_path_file).absolute()

    complete_folder_multifasta = args.input_folder_multifasta
    if complete_folder_multifasta:
        complete_folder_multifasta = Path(complete_folder_multifasta).absolute()

    sequence_topology_manifest = None
    if args.sequence_topology_manifest:
        try:
            sequence_topology_manifest = load_sequence_topology_manifest(
                Path(args.sequence_topology_manifest).absolute()
            )
        except SequenceTopologyError as error:
            parser.error("invalid sequence topology manifest: {}".format(error))

    optional_errors = optional_analysis_errors(
        project_root=script_absolute_path,
        is_element=args.is_element,
        cas=args.cas,
    )
    if optional_errors:
        parser.error("\n".join(optional_errors))

    if args.fuzzy_search_mode == "native" and not literal_fuzzy.AVAILABLE:
        parser.error(
            "the validated native fuzzy matcher is unavailable: {}; run "
            "scripts/setup_local.sh or scripts/build_native_matcher.sh"
            .format(literal_fuzzy.LOAD_ERROR or "unknown loader error")
        )

    folder_result = args.result_folder
    if folder_result:
        folder_result = Path(folder_result).absolute()

    pickle_folder = args.pickle_report
    if pickle_folder:
        pickle_folder = Path(pickle_folder).absolute()

    json_folder = args.json_report
    if json_folder:
        json_folder = Path(json_folder).absolute()


    flag_enhancement_max_min = args.enhancement_max_min
    flag_enhancement_start_end = args.enhancement_start_end
    flag_enhancement_intermediate = args.enhancement_intermediate

    flag_parallel = args.parallel
    flag_fast_run = args.fast_run

    strand_flag = args.strand
    is_flag = args.is_element
    cas_flag = args.cas
    degenerated_flag = args.degenerated
    fasta_report = args.fasta_report

    flags = {"flag_parallel": flag_parallel,
             "flag_cpu": flag_cpu,
             "flag_fast_run": flag_fast_run,
             "fast_run_seed_profile": args.fast_run_seed_profile,
             "flag_fast_run_rescue": args.fast_run_rescue,
             "fast_run_rescue_seed_scope": args.fast_run_rescue_seed_scope,
             "fast_run_rescue_trim_profile": (
                 args.fast_run_rescue_trim_profile
             ),
             "flag_fast_run_phase_repair": args.fast_run_phase_repair,
             "flag_fast_two_repeat_terminal_repair": (
                 args.fast_two_repeat_terminal_repair
             ),
             "flag_fast_accepted_family_rescue": (
                 args.fast_accepted_family_rescue
             ),
             "flag_strand": strand_flag,
             "flag_is": is_flag,
             "flag_cas": cas_flag,
             "cas_gene_mode": cas_gene_mode,
             "cas_timeout_seconds": args.cas_timeout_seconds,
             "flag_fasta_report": fasta_report,
             "flag_html_report": args.html_report,
             "flag_biological_loci": args.biological_loci,
             "output_profile": args.output_profile,
             "flag_degenerated": degenerated_flag,
             "flag_enhancement_min_max": flag_enhancement_max_min,
             "flag_enhancement_start_end": flag_enhancement_start_end,
             "flag_enhancement_intermediate": flag_enhancement_intermediate,
             "flag_fuzzy_boundary_paths": args.fuzzy_boundary_paths,
             }

    min_rep = args.min_len_rep
    max_rep = args.max_len_rep
    max_spacer = args.max_len_spacer
    min_spacer = args.min_len_spacer
    min_repeats = args.min_repeats
    max_identical_spacers = args.max_identical_spacers
    max_identical_cluster_spacers = args.max_identical_cluster_spacers
    margin_degenerated = args.margin_degenerated
    max_edit_distance_enhancement = args.max_edit_distance_enhanced
    if args.folder_jobs < 1:
        parser.error('--folder_jobs must be a positive integer')
    if min_rep < 1 or max_rep < min_rep:
        parser.error("repeat lengths must satisfy 1 <= min_len_rep <= max_len_rep")
    if min_spacer < 0 or max_spacer < min_spacer:
        parser.error(
            "spacer lengths must satisfy 0 <= min_len_spacer <= max_len_spacer"
        )
    if min_repeats < 2:
        parser.error("--min_repeats must be at least 2")
    if max_identical_spacers < 0 or max_identical_cluster_spacers < 0:
        parser.error("identical-spacer limits must be non-negative")
    if margin_degenerated < 0:
        parser.error("--margin_degenerated must be non-negative")
    if max_edit_distance_enhancement < 0:
        parser.error("--max_edit_distance_enhanced must be non-negative")
    if args.fuzzy_search_flank < 0:
        parser.error('--fuzzy_search_flank must be non-negative')
    if args.fuzzy_family_padding < 0:
        parser.error('--fuzzy_family_padding must be non-negative')

    effective_scientific_config = {
        destination: getattr(args, destination)
        for destination in BALANCED_DEFAULTS
    }
    effective_scientific_config.update({
        "biological_loci": args.biological_loci,
        "cas_gene_mode": cas_gene_mode,
        "cas_timeout_seconds": args.cas_timeout_seconds,
        "enhancement_intermediate": args.enhancement_intermediate,
        "enhancement_max_min": args.enhancement_max_min,
        "enhancement_start_end": args.enhancement_start_end,
        "fasta_report": args.fasta_report,
        "fast_run_rescue_seed_scope": args.fast_run_rescue_seed_scope,
        "fast_run_rescue_trim_profile": args.fast_run_rescue_trim_profile,
        "html_report": args.html_report,
        "legacy_json_report": bool(args.json_report),
        "margin_degenerated": margin_degenerated,
        "max_edit_distance_enhanced": max_edit_distance_enhancement,
        "max_identical_cluster_spacers": max_identical_cluster_spacers,
        "max_identical_spacers": max_identical_spacers,
        "max_len_rep": max_rep,
        "max_len_spacer": max_spacer,
        "min_len_rep": min_rep,
        "min_len_spacer": min_spacer,
        "min_repeats": min_repeats,
        "pickle_report": bool(args.pickle_report),
        "sequence_topology_manifest": bool(args.sequence_topology_manifest),
        "split_competitive_components": split_competitive_components,
    })
    runtime_preset_provenance = preset_provenance(
        args.preset,
        effective_scientific_config,
        explicit_destinations=explicit_destinations,
        preset_explicit="preset" in explicit_destinations,
    )
    effective_config_document = {
        "schema_version": 1,
        "release_version": _release_version,
        "preset": args.preset,
        "preset_provenance": runtime_preset_provenance,
        "effective_config": effective_scientific_config,
    }
    if args.retain_candidate_evidence:
        effective_config_document["research_outputs"] = {
            "retain_candidate_evidence": True,
        }

    if args.print_effective_config:
        print(json.dumps(
            effective_config_document,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        ))
        raise SystemExit(0)

    parameters = {
        "param_min_avg_repeat_length": min_rep,
        "param_max_avg_repeat_length": max_rep,
        "param_max_avg_spacer_length": max_spacer,
        "param_min_avg_spacer_length": min_spacer,
        "param_min_repeats": min_repeats,
        "param_max_identical_spacers": max_identical_spacers,
        "param_max_identical_cluster_spacers": max_identical_cluster_spacers,
        "param_spacer_margin_degenerated_search": margin_degenerated,
        "param_max_edit_distance": max_edit_distance_enhancement,
        "param_runtime_preset": args.preset,
        "param_release_version": _release_version,
        "param_runtime_preset_provenance": runtime_preset_provenance,
        "param_effective_config": effective_scientific_config,
        "param_fuzzy_boundary_paths": args.fuzzy_boundary_paths,
        "param_strand": args.strand,
        "param_is_element": args.is_element,
        "param_cas": args.cas,
        "param_biological_loci": args.biological_loci,
        "param_enhancement_max_min": args.enhancement_max_min,
        "param_enhancement_start_end": args.enhancement_start_end,
        "param_enhancement_intermediate": args.enhancement_intermediate,
        "param_source_backed_terminal_completion": (
            args.architecture_rescue and args.degenerated
        ),
        "param_source_backed_architecture_rescue": args.architecture_rescue,
        "param_source_terminal_completion_rule": (
            "one_step_substitution_only_rate_0.20_cap_8_anchor_10_"
            "spacer_median_3mad_10_110"
        ),
        "param_accepted_only": args.accepted_only,
        "param_fast_run": args.fast_run,
        "param_fast_run_seed_profile": args.fast_run_seed_profile,
        "param_fast_run_rescue": bool(
            args.fast_run and args.fast_run_rescue
        ),
        "param_fast_run_rescue_seed_scope": args.fast_run_rescue_seed_scope,
        "param_fast_run_rescue_trim_profile": (
            args.fast_run_rescue_trim_profile
        ),
        "param_fast_run_phase_repair": bool(
            args.fast_run and args.fast_run_phase_repair
        ),
        "param_fast_two_repeat_terminal_repair": bool(
            args.fast_run and args.fast_two_repeat_terminal_repair
        ),
        "param_fast_accepted_family_rescue": bool(
            args.fast_run and args.fast_accepted_family_rescue
        ),
        "param_degenerated": args.degenerated,
        "param_output_profile": args.output_profile,
        "param_split_competitive_components": split_competitive_components,
        "param_fuzzy_search_mode": args.fuzzy_search_mode,
        "param_fuzzy_search_flank": args.fuzzy_search_flank,
        "param_fuzzy_family_padding": args.fuzzy_family_padding,
        "param_cas_gene_mode": cas_gene_mode,
        "param_cas_timeout_seconds": args.cas_timeout_seconds,
        "param_reconstruction_policy": args.reconstruction_policy,
    }


    # Runtime intermediates belong in the operating system temporary area, not
    # beside the source tree.  All user-selected input/output paths were made
    # absolute above, and project assets are resolved from ``script_absolute_path``.
    pid_work_directory = tempfile.mkdtemp(
        prefix="crispridentify_{}_".format(pid)
    )
    os.chdir(pid_work_directory)


def _status(message):
    """Write concise progress unless routine output was disabled."""

    if not args.quiet:
        print(message, flush=True)


def _execute_pipeline(**kwargs):
    """Run one genome with detailed pipeline chatter only when requested."""

    kwargs["retain_candidate_evidence"] = args.retain_candidate_evidence
    if args.verbose:
        return Pipeline(**kwargs)
    with open(os.devnull, "w") as sink:
        with contextlib.redirect_stdout(sink):
            return Pipeline(**kwargs)


def _prepare_input_files(folder):
    """Return deterministic source names without modifying the input folder."""

    folder = os.fspath(folder)
    files = sorted(
        name for name in listdir(folder) if isfile(join(folder, name))
    )
    unsafe = [
        name for name in files
        if any(ord(character) < 32 or ord(character) == 127 for character in name)
    ]
    if unsafe:
        raise InputValidationError(
            "input filename contains a control character and cannot be "
            "processed safely: {}; rename it explicitly; CRISPRidentify did "
            "not modify the file".format(repr(unsafe[0]))
        )

    validate_unique_result_keys(join(folder, name) for name in files)
    return files


def _ensure_output_roots(result_folder, pickle_folder, json_folder):
    os.makedirs(os.fspath(result_folder), exist_ok=True)
    if pickle_folder:
        os.makedirs(os.fspath(pickle_folder), exist_ok=True)
    if json_folder:
        os.makedirs(os.fspath(json_folder), exist_ok=True)


def _finalize_folder_outputs(result_folder, json_folder, result_keys):
    """Build global artifacts once, after every genome has completed."""
    if args.output_profile == "machine":
        return
    result_keys = tuple(sorted(result_keys))
    CompleteFolderSummaryMaker(
        folder_result=result_folder, sub_folders=result_keys
    )
    if fasta_report:
        CompleteFastaOutputMaker(
            folder_result=result_folder, sub_folders=result_keys
        )
    if json_folder:
        CompleteJsonOutputMaker(
            folder_json_result=json_folder,
            folder_text_tesult=result_folder,
            result_keys=result_keys,
        )


def run_over_folder_of_files(
        folder, result_folder, pickle_folder, json_folder,
        chunk_number=None, number_of_chunks=None, finalize=True,
):
    files = _prepare_input_files(folder)
    _ensure_output_roots(result_folder, pickle_folder, json_folder)

    if number_of_chunks:
        chunk_size = math.ceil(len(files) / number_of_chunks)
        chunk_start = (chunk_number - 1) * chunk_size
        chunk_end = chunk_number * chunk_size
        chunk = files[chunk_start:chunk_end]
        if args.verbose:
            print(chunk_start)
            print(chunk_end)
    else:
        chunk = files

    for index, file_name in enumerate(chunk, 1):
        if args.verbose:
            print(
                "\n\n\n\t\t\t\tExecuting file {} out of {} ({})\n\n\n"
                .format(index, len(chunk), file_name)
            )
        else:
            _status("Processing {}/{}: {}".format(
                index, len(chunk), file_name
            ))
        _execute_pipeline(
            result_folder_path="{}/".format(result_folder),
            pickle_folder_path="{}".format(pickle_folder),
            json_folder_path="{}".format(json_folder),
            file_path=join(folder, file_name),
            parameters=parameters,
            flags=flags,
            flag_dev_mode=FLAG_DEVELOPER_MODE,
            sequence_topology_manifest=sequence_topology_manifest,
        )
        if not args.verbose:
            _status("Completed {}/{}: {}".format(
                index, len(chunk), file_name
            ))

    if finalize:
        _finalize_folder_outputs(
            result_folder,
            json_folder,
            (pipeline_result_key(name) for name in chunk),
        )


def _tail_text(path, line_count=80):
    try:
        with open(path, errors="replace") as stream:
            return "".join(stream.readlines()[-line_count:])
    except OSError as error:
        return "unable to read lane log: {}".format(error)


def _claim_next_queue_file(queue_root, lane_index):
    """Atomically claim the next largest-first queue token."""
    pending = os.path.join(queue_root, "pending")
    claimed = os.path.join(queue_root, "claimed", "lane_{:02d}".format(lane_index))
    while True:
        tokens = sorted(listdir(pending))
        if not tokens:
            return None
        for token in tokens:
            source = os.path.join(pending, token)
            destination = os.path.join(claimed, token)
            try:
                os.rename(source, destination)
            except FileNotFoundError:
                if not os.path.isdir(claimed):
                    raise RuntimeError(
                        "dynamic queue lane directory disappeared: {}".format(claimed)
                    )
                if os.path.lexists(source):
                    raise
                continue
            return destination


def run_dynamic_queue_lane(
        queue_root, lane_index, result_folder, pickle_folder, json_folder,
):
    """Consume a shared atomic queue without reloading runtime state."""
    _ensure_output_roots(result_folder, pickle_folder, json_folder)
    completed = 0
    while True:
        claimed_token = _claim_next_queue_file(queue_root, lane_index)
        if claimed_token is None:
            break
        source_path = os.readlink(claimed_token)
        completed += 1
        if args.verbose:
            print(
                "\n\n\n\t\t\t\tLane {} executing queue item {} ({})\n\n\n"
                .format(lane_index, completed, os.path.basename(source_path))
            )
        else:
            _status("Lane {} processing {}: {}".format(
                lane_index, completed, os.path.basename(source_path)
            ))
        _execute_pipeline(
            result_folder_path="{}/".format(result_folder),
            pickle_folder_path="{}".format(pickle_folder),
            json_folder_path="{}".format(json_folder),
            file_path=source_path,
            parameters=parameters,
            flags=flags,
            flag_dev_mode=FLAG_DEVELOPER_MODE,
            sequence_topology_manifest=sequence_topology_manifest,
        )
    _status("Lane {} processed {} genome(s)".format(lane_index, completed))


def _claimed_queue_inputs(queue_root, lane_index):
    claimed = os.path.join(queue_root, "claimed", "lane_{:02d}".format(lane_index))
    names = []
    for token in sorted(listdir(claimed)):
        token_path = os.path.join(claimed, token)
        try:
            names.append(os.path.basename(os.readlink(token_path)))
        except OSError:
            names.append(token)
    return names


def run_over_folder_in_parallel(folder, result_folder, pickle_folder, json_folder):
    """Run persistent subprocess lanes with a dynamic largest-first queue."""
    file_names = _prepare_input_files(folder)
    _ensure_output_roots(result_folder, pickle_folder, json_folder)
    if not file_names:
        _finalize_folder_outputs(result_folder, json_folder, ())
        return

    selected_cpus = resolve_cpu_ids(args.cpu, allowed_cpus=get_allowed_cpu_ids())
    lane_count = min(args.folder_jobs, len(file_names), len(selected_cpus))
    if lane_count <= 1:
        run_over_folder_of_files(
            folder, result_folder, pickle_folder, json_folder, finalize=True
        )
        return

    cpu_lanes = partition_cpu_lanes(selected_cpus, lane_count)
    input_paths = tuple(Path(folder) / name for name in file_names)
    ordered_paths = sorted(
        input_paths, key=lambda path: (-path.stat().st_size, str(path))
    )
    queue_root = tempfile.mkdtemp(prefix="folder_queue_", dir=pid_work_directory)
    pending = os.path.join(queue_root, "pending")
    claimed_root = os.path.join(queue_root, "claimed")
    os.makedirs(pending)
    os.makedirs(claimed_root)
    for lane_index in range(1, lane_count + 1):
        os.makedirs(os.path.join(
            claimed_root, "lane_{:02d}".format(lane_index)
        ))
    for rank, source in enumerate(ordered_paths):
        token = os.path.join(pending, "{:08d}".format(rank))
        # Point at the input path rather than its resolved target: Pipeline's
        # result key is intentionally derived from the benchmark input name.
        os.symlink(str(source.absolute()), token)

    processes = []
    try:
        for lane_index, cpu_ids in enumerate(cpu_lanes, 1):
            lane_folder = os.path.join(
                claimed_root, "lane_{:02d}".format(lane_index)
            )
            child_arguments = rewrite_lane_argv(
                sys.argv[1:],
                input_folder=lane_folder,
                result_folder=result_folder,
                cpu=len(cpu_ids),
                folder_jobs=1,
            ) + (
                "--_folder_queue", queue_root,
                "--_folder_lane_id", str(lane_index),
                "--_folder_cpu_ids", ",".join(str(cpu_id) for cpu_id in cpu_ids),
            )
            command = (
                sys.executable,
                os.path.join(script_absolute_path, "CRISPRidentify.py"),
            ) + child_arguments
            log_path = os.path.join(queue_root, "lane_{:02d}.log".format(lane_index))
            log_handle = open(log_path, "w")
            environment = os.environ.copy()
            environment.update(native_thread_environment(len(cpu_ids)))
            try:
                process = subprocess.Popen(
                    command,
                    cwd=work_directory,
                    env=environment,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except Exception:
                log_handle.close()
                raise
            processes.append((
                lane_index, process, log_handle, log_path, cpu_ids
            ))
            _status(
                "Folder lane {}: dynamic queue, CPUs {}, log {}".format(
                    lane_index,
                    ",".join(str(cpu_id) for cpu_id in cpu_ids),
                    log_path,
                )
            )

        failures = []
        for lane_index, process, log_handle, log_path, cpu_ids in processes:
            return_code = process.wait()
            log_handle.close()
            claimed_names = _claimed_queue_inputs(queue_root, lane_index)
            if return_code == 0:
                _status(
                    "Folder lane {} completed {} genome(s)".format(
                        lane_index, len(claimed_names)
                    )
                )
            else:
                failures.append((
                    lane_index, return_code, log_path, claimed_names
                ))

        if failures:
            details = []
            for lane_index, return_code, log_path, claimed_names in failures:
                details.append(
                    "lane {} failed with exit code {} after claiming {}\n{}".format(
                        lane_index,
                        return_code,
                        ", ".join(claimed_names),
                        _tail_text(log_path),
                    )
                )
            raise RuntimeError("\n\n".join(details))

        remaining = listdir(pending)
        processed_count = sum(
            len(_claimed_queue_inputs(queue_root, lane_index))
            for lane_index in range(1, lane_count + 1)
        )
        if remaining or processed_count != len(file_names):
            raise RuntimeError(
                "dynamic folder queue incomplete: {} pending, {} of {} claimed".format(
                    len(remaining), processed_count, len(file_names)
                )
            )
    except BaseException:
        for _lane_index, process, log_handle, _log_path, _cpu_ids in processes:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        for _lane_index, process, log_handle, _log_path, _cpu_ids in processes:
            if process.poll() is None:
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
            if not log_handle.closed:
                log_handle.close()
        raise

    _finalize_folder_outputs(
        result_folder,
        json_folder,
        (pipeline_result_key(name) for name in file_names),
    )




def _dispatch_folder(folder):
    if args._folder_queue:
        if args._folder_lane_id is None:
            raise ValueError("dynamic folder lane requires a lane id")
        run_dynamic_queue_lane(
            args._folder_queue,
            args._folder_lane_id,
            folder_result,
            pickle_folder,
            json_folder,
        )
    elif args.folder_jobs > 1 and not args._folder_lane:
        run_over_folder_in_parallel(
            folder, folder_result, pickle_folder, json_folder
        )
    else:
        run_over_folder_of_files(
            folder,
            folder_result,
            pickle_folder,
            json_folder,
            finalize=not args._folder_lane,
        )


def main():
    start_time = time()
    temporary_input_folder = None
    try:
        if complete_path_file:
            temporary_input_folder = multiline_fasta_handle_python(
                complete_path_file, flag_ncbi_formatting=True
            )
            _status("Prepared input sequence(s)")
            _dispatch_folder(temporary_input_folder)
        elif complete_path_folder:
            _dispatch_folder(complete_path_folder)
        elif complete_folder_multifasta:
            _status("Preparing multifasta inputs")
            temporary_input_folder = folder_of_multifasta_handle(
                complete_folder_multifasta
            )
            _dispatch_folder(temporary_input_folder)
        else:
            parser.error(
                "one of --file, --input_folder, or "
                "--input_folder_multifasta is required"
            )
    except InputValidationError as error:
        parser.error(str(error))
    finally:
        if temporary_input_folder and os.path.exists(temporary_input_folder):
            shutil.rmtree(temporary_input_folder)

    _status("Finished in {:.2f} s; results: {}".format(
        time() - start_time, folder_result
    ))


if __name__ == "__main__":
    try:
        main()
    finally:
        os.chdir(work_directory)
        shutil.rmtree(pid_work_directory, ignore_errors=True)
