from argparse import ArgumentParser, ArgumentTypeError
import math
import os
import json
from pathlib import Path

from sp_model.input_parser import InputParser
from sp_model.helpers.misc import create_logger
from sp_model.summary import compose_summary_dict, write_summary
from sp_model.run_experiments import run_pickled_data, run_multiple_groups
import sp_model.helpers.import_data as import_data
from sp_model.crisprfinder_to_spacerfasta import CRISPRFinderParser

__version__ = '1.0.4'


def _positive_float(value):
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0:
        raise ArgumentTypeError('value must be finite and greater than zero')
    return converted


def _non_negative_float(value):
    converted = float(value)
    if not math.isfinite(converted) or converted < 0:
        raise ArgumentTypeError('value must be finite and non-negative')
    return converted


def _non_negative_int(value):
    converted = int(value)
    if converted < 0:
        raise ArgumentTypeError('value must be non-negative')
    return converted


def _positive_int(value):
    converted = int(value)
    if converted <= 0:
        raise ArgumentTypeError('value must be a positive integer')
    return converted


def _probability(value):
    converted = float(value)
    if not 0 < converted < 1:
        raise ArgumentTypeError('value must be strictly between zero and one')
    return converted


def parse_args(argv=None, prog='spacerplacer', description=None):
    parser = ArgumentParser(prog=prog,
                            description=description or
                            'SpacerPlacer: An Ancestral Reconstruction Algorithm for CRISPR Arrays',
                            )
    ############################################################################################################
    # General:
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('--verbosity', type=int, default=1, choices=[0, 1, 2],
                        help='Verbosity level: 0: no output, 1: minimal output, 2: maximal output (default 1).')
    parser.add_argument('--seed', type=int, default=2357, help='Seed for the random number generator '
                                                               'only impact on MAFFT.')
    parser.add_argument('--alignment_workers', type=_positive_int, default=1,
                        help='Number of independent CRISPR groups to align concurrently. '
                             'Use the CPU allocation supplied to the web job (default 1).')
    parser.add_argument('--aligned_groups_path', type=str, default=None,
                        help='Validated aligned-groups JSON from a previous SpacerPlacer stage. '
                             'When supplied, MAFFT is skipped and these exact alignments are reused.')

    # parser.add_argument('--group_by_metadata', type=str, default=None,
    #                     help='Not supported.')
    ############################################################################################################
    # Input/Output locations + format:
    parser.add_argument('input_path', type=str,
                        help='Path to input file or folder containing input files. The required input format is described '
                             'in the readme.')
    parser.add_argument('output_path', type=str,
                        help='Path to output folder. The output folder will be created, if it does not exist. '
                             'The output folder will contain files as described in the readme.')
    parser.add_argument('--tree_path', type=str, default=None,
                        help='Path to tree json file or folder with newick files. '
                             'If none is provided, trees are estimated by SpacerPlacer. '
                             'The trees can be given in a newick format or as a dictionary in a json file '
                             '(such a json file is also returned by SpacerPlacer).')
    parser.add_argument('-it', '--input_type', type=str, choices=['spacer_fasta', 'sf', 'pickled', 'ccdb',
                                                                  'crisprcasdb',
                                                                  'ccf_json', 'crisprcasfinder_json'],
                        default='spacer_fasta',  # think about this name
                        help='Determines the input type, i.e. either already preprocessed fasta style spacer arrays, '
                             'a pickled file with CRISPR group(s) (our own data structure), '
                             'DNA data (extracted from CRISPRCasdb) in a structure described in the '
                             'github repository or directly from result.json files '
                             'returned by CRISPRCasFinder.')
    parser.add_argument('--cluster_spacers', action='store_true',
                        help='If given, the spacers of the input data are clustered (using levenshtein distance).'
                             'The clustering process is described in the paper.')
    parser.add_argument('--cluster_spacers_max_distance', type=_non_negative_int, default=1,
                        help='Determines the maximum Levenshtein distance for spacers to be clustered, i.e. identified '
                             'with each other (same spacer id) within repeat groups. ')
    # If the user provides CRISPRCasFinder result json as input format (requires --input type ccf_result_json):
    parser.add_argument('--min_evidence_level', type=int, choices=range(1, 5), default=4,
                        help='Determines the minimum evidence level for the CRISPRCasFinder results to be included in '
                             'the reconstruction. The evidence level is given in the CRISPRCasFinder output. '
                             'The default (4) is the highest possible value.')

    ############################################################################################################
    # Reconstruction:
    parser.add_argument('--rec_model', type=str, choices=['gtr', ],
                        default='gtr', help='Determines the reconstruction model. Currently redundant, '
                                            'as only "gtr" is implemented.')
    parser.add_argument('--insertion_rate', type=_positive_float, default=1e-5,
                        help='Insertion rate of the reconstruction model. Generally should be chosen to be small '
                             'compared to the deletion rate. Otherwise may lead to high number of independent '
                             'acquisitions and worse reconstructions. The ratio between insertion and deletion rate '
                             'can be tuned in case of reconstructions with excessive independent acquisitions or excessive '
                             'accumulations of insertions at the root. '
                             'Default insertion and deletion rates are provided and '
                             'work well for the tested datasets.')
    parser.add_argument('--deletion_rate', type=_positive_float, default=1e-1,
                        help='Deletion rate of the reconstruction model. Generally should be chosen to be large '
                             'compared to the insertion rate. Otherwise may lead to high number of independent '
                             'acquisitions and worse reconstructions. The ratio between insertion and deletion rate '
                             'can be tuned in case of reconstructions with excessive independent acquisitions or excessive '
                             'accumulations of insertions at the root. '
                             'Default insertion and deletion rates are provided and '
                             'work well for the tested datasets.')
    # Don't know if this is needed:
    parser.add_argument('--extend_branches', type=_non_negative_float, default=0.0,
                        help='Extends branches of the tree by the given value. This is useful, if the tree is not '
                             'well resolved to allow placement of events on '
                             'small (or length=0) branches.')
    # parser.add_argument('--stationary_dist_at_root', action='store_true', default=False,
    #                     help="If given, the likelihood at the root is computed based on the stationary distribution. "
    #                          "Otherwise the maximum of the child's branch lengths and the standard likelihood "
    #                          "computation is used. The stationary distribution might lead to undesired results, "
    #                          "that are not parsimonious. We do not recommend to use this flag.")
    # ########################################################################################################### Tree
    # tree estimation:
    parser.add_argument('--tree_distance_function', type=str, choices=['likelihood'],
                        default='likelihood', help='Determines the distance function used for the tree estimation '
                                                   'with UPGMA. Currently redundant, as only "likelihood" is implemented.')
    parser.add_argument('--tree_construction_method', type=str, choices=['upgma', 'nj'],
                        default='upgma', help='Determines the tree construction method used for the tree estimation. '
                                              'Currently UPGMA (upgma) and neighbor joining (nj) are implemented.')
    parser.add_argument('--tree_insertion_rate', type=_positive_float, default=None,
                        help='The user can provide their own insertion rate for the tree estimation based on the '
                             'block deletion likelihood function. '
                             'In the default case, well performing parameters are chosen '
                             'depending on the used likelihood function ("simple", "ode_based").')
    parser.add_argument('--tree_deletion_rate', type=_positive_float, default=None,
                        help='The user can provide their own deletion rate for the tree estimation based on the '
                             'block deletion likelihood function. '
                             'In the default case, well performing parameters are chosen '
                             'depending on the used likelihood function ("simple", "ode_based").'
                        )
    parser.add_argument('--tree_alpha', type=_positive_float, default=None,
                        help='The user can provide their own deletion rate for the tree estimation based on the '
                             'block deletion likelihood function. '
                             'In the default case, well performing parameters are chosen '
                             'depending on the used likelihood function ("simple", "ode_based").'
                        )
    parser.add_argument('--tree_lh_fct', type=str, choices=['simple', 'ode_based'], default='simple',
                        help='Determines the likelihood function used for the tree estimation with UPGMA. '
                             'Details about the likelihood functions can be found in the paper. The default "simple" has '
                             'been found to perform better in simulations.')
    parser.add_argument('--combine_non_unique_arrays', action='store_true',
                        help='If given, arrays with the exactly the same spacer content are combined before the '
                             'tree estimation and reconstruction. This might be helpful reduce clutter '
                             '(especially in visualization).'
                             'Arrays with the same spacer content can not be separated in the tree '
                             'estimation anyway. '
                             'The default is NOT to combine arrays with the same spacer content. '
                             'If arrays are combined, the combined array names are saved in detailed results csv under '
                             'column "combined array names".')

    ############################################################################################################
    # Visualization:
    parser.add_argument('--no_plot_reconstruction', action='store_true',
                        help='If given, the reconstruction is not plotted.')
    parser.add_argument('--no_plot_order_graph', action='store_true',
                        help='If given, the Partial Spacer Insertion Order is not plotted.')
    parser.add_argument('--dpi_rec', type=int, default=90,
                        help='DPI for the reconstruction plot pdf. Default is 90.')
    parser.add_argument('--figsize_rec', nargs=3, default=[None, None, 'px'], metavar=('WIDTH', 'HEIGHT',
                                                                                       '{px, mm, in}'),
                        help='Provide the width and height of the reconstruction plot in the provided unit "px" (pixel), "mm" '
                             '(millimeter) or "in" (inches). '
                             'By default the size is determined by the drawing function.')
    parser.add_argument('--show_spacer_name_row', action='store_false',
                        help='If given, the spacer original names are hidden in the reconstruction plot (above '
                             'the alignment and the respective spacer ids). '
                             'By default the spacer names are shown.')
    # Does this work?
    parser.add_argument('--do_show', action='store_true',
                        help='If given, the plots are shown directly. Might lead to crash due to some issue with '
                             'pyqt5.')

    ############################################################################################################
    # Evaluation (LRT parameters):
    parser.add_argument('--lh_fct', type=str, choices=['simple', 'ode_based'], default='simple',
                        help='Determines the likelihood function used for the Block deletion model estimates and '
                             'the likelihood ratio test. "simple" allows bias '
                             'corrections in rho and alpha and was found to perform better in simulations. '
                             '"ode_based" is the less performant, but, in theory, more precise likelihood function.')
    parser.add_argument('--no_alpha_bias_correction', action='store_true',
                        help='If given, omits the alpha bias correction. Best performance is achieved with '
                             'the bias correction. Bias correction is only possible, '
                             'if simple likelihood function is used.')
    parser.add_argument('--no_rho_bias_correction', action='store_true',
                        help='If given, omits the rho bias correction. Best performance is achieved with '
                             'the bias correction. Bias correction is only possible, '
                             'if simple likelihood function is used.')
    parser.add_argument(
        '--rho_bias_failure_policy',
        choices=['fail', 'fallback_uncorrected'],
        default='fail',
        help=(
            'Policy when the optional unobserved-spacer rho correction has no finite fit. '
            'fail preserves strict behavior; fallback_uncorrected retries only the affected group/orientation '
            'without rho correction and records the effective policy (default: fail). Rate estimates change; '
            'optimizer likelihoods and the IDM/BDM LRT are not recomputed by rho correction.'
        ),
    )
    parser.add_argument('--significance_level', type=_probability, default=0.05,
                        help='Determines the significance level for the likelihood ratio test between the Independent '
                             'Deletion Model and the Block Deletion Model. '
                             'The test statistic (-2*ln(lh_idm/ln_bdm) = 2*ln(lh_bdm - ln_idm) is chi-squared '
                             "distributed with one single degree of freedom by Wilks' Theorem.")

    ############################################################################################################
    # Orientation determination:
    parser.add_argument('--determine_orientation', action='store_true',
                        help='If given, both a forward and a reverse '
                             'reconstruction is made and the orientation of the arrays is determined by SpacerPlacer. '
                             'Otherwise, the orientation is '
                             'accepted as provided. Note, if the tree is estimated by SpacerPlacer, '
                             'a separate tree ist estimated for the reversed array, i.e. the trees between forward and '
                             'reverse orientation might differ.'
                        )
    parser.add_argument('--orientation_decision_boundary', type=_non_negative_float, default=5,
                        help='Determines the decision boundary for the orientation decision. If the modulus of the '
                             'difference between the log likelihoods of the two orientations is less than or equal to '
                             'the decision boundary, the orientation is set to "not determined" (ND). '
                             'Otherwise, the orientation '
                             'is set to the orientation with the higher likelihood. The decision boundary is given in '
                             'logscale.'
                             'The default value %(default)s was determined empirically and works well for the tested '
                             'datasets.')
    ############################################################################################################
    # Additional data:
    parser.add_argument('--save_reconstructed_events', action='store_true',
                        help='If given, the tree and reconstructed events along the tree are saved in '
                             '"reconstructed_events". On the basis of this data, the reconstruction is '
                             'visualized and the events can be analyzed in more detail. Event export also works with '
                             'reconstruction plotting disabled.')

    parser.add_argument('--allow_unsafe_pickle', action='store_true',
                        help='Allow loading user-supplied pickle files. Pickle can execute arbitrary code and is '
                             'disabled by default for web-server safety.')
    parser.add_argument('--allow_existing_output', action='store_true',
                        help='Allow writing into a non-empty output directory. Use only for trusted resume workflows.')

    return parser.parse_args(argv)


def validate_paths_and_security(args):
    """Validate user-controlled paths before creating or writing output files."""
    input_path = Path(args.input_path).expanduser().resolve()
    output_path = Path(args.output_path).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f'Input path does not exist: {input_path}')
    if input_path == output_path:
        raise ValueError('Input and output paths must be different')
    if input_path.is_dir() and input_path in output_path.parents:
        raise ValueError('Output path must not be located inside the input directory')
    if output_path.exists() and not output_path.is_dir():
        raise ValueError(f'Output path exists and is not a directory: {output_path}')
    if output_path.exists() and any(output_path.iterdir()) and not args.allow_existing_output:
        raise FileExistsError(
            f'Output directory is not empty: {output_path}. Use --allow_existing_output only for trusted resumes.'
        )
    if args.input_type == 'pickled' and not args.allow_unsafe_pickle:
        raise ValueError('Pickled input is disabled by default; pass --allow_unsafe_pickle only for trusted data')
    if args.dpi_rec <= 0:
        raise ValueError('dpi_rec must be greater than zero')
    if args.tree_path is not None:
        tree_path = Path(args.tree_path).expanduser().resolve()
        if not tree_path.exists():
            raise FileNotFoundError(f'Tree path does not exist: {tree_path}')
        if tree_path.suffix.lower() in {'.pkl', '.pickle'} and not args.allow_unsafe_pickle:
            raise ValueError('Pickled trees are disabled by default; pass --allow_unsafe_pickle only for trusted data')
    return input_path, output_path


def check_and_parse_input_data(input_type, input_path, output_path, tree_path, logger, cluster_spacers=False,
                               cluster_spacers_max_distance=1, min_evidence_level=4,):
    if input_type == 'spacer_fasta' or input_type == 'sf':
        if os.path.splitext(input_path)[-1] in {'.fa', '.fasta', '.fna'}:
            ls_path_to_spacer_fasta = [input_path]
        else:
            ls_path_to_spacer_fasta = [os.path.join(input_path, group) for group in sorted(os.listdir(input_path))
                                       if group.split('.')[-1].lower() in {'fa', 'fasta', 'fna'}
                                       and os.path.isfile(os.path.join(input_path, group))]

        if not ls_path_to_spacer_fasta:
            logger.error(f'No fasta file found in {input_path}.')
            raise ValueError(f'No fasta file found in {input_path}.')
    elif input_type in ['ccf_json', 'crisprcasfinder_json',]:
        if os.path.isfile(input_path):
            ls_result_json_paths = [input_path]
        else:
            direct_result = os.path.join(input_path, 'result.json')
            ls_result_json_paths = [direct_result] if os.path.isfile(direct_result) else []
            ls_result_json_paths.extend(
                os.path.join(input_path, child, 'result.json')
                for child in sorted(os.listdir(input_path))
                if os.path.isfile(os.path.join(input_path, child, 'result.json'))
            )
        if not ls_result_json_paths:
            raise ValueError(f'No CRISPRCasFinder result.json found in {input_path}.')
        output_folder = os.path.join(output_path, 'additional_data', 'ccf_json_clustered_groups')
        output_folder_too_small = os.path.join(output_path, 'additional_data', 'ccf_json_too_small_clustered_groups')
        additional_output = os.path.join(output_path, 'additional_data', 'ccf_json_clustering_info')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        parser = CRISPRFinderParser(ls_result_json_paths, output_folder=output_folder,
                                    output_folder_for_too_small_groups=output_folder_too_small,
                                    min_evidence_level=min_evidence_level,
                                    reorient_by_direction_pred=True, logger=logger,
                                    cluster_by_spacer_overlap=True,
                                    filter_more_than_x_in_group=2, # might want to adjust this
                                    max_levenshtein_distance_spacers=cluster_spacers_max_distance
                                    if cluster_spacers else 0)
        parser.write_spacer_cluster_to_file(additional_output)
        parser.write_spacer_seq_to_spacer_id(additional_output)
        ls_path_to_spacer_fasta = [os.path.join(output_folder, group) for group in sorted(os.listdir(output_folder))
                                   if os.path.isfile(os.path.join(output_folder, group))]
    else:
        ls_path_to_spacer_fasta = []
        path_to_spacer_fasta_folder = os.path.join(output_path, 'additional_data', 'spacer_fasta')
        if not os.path.exists(path_to_spacer_fasta_folder):
            os.makedirs(path_to_spacer_fasta_folder)
        for group in sorted(os.listdir(input_path)):
            if not os.path.isdir(os.path.join(input_path, group)):
                continue
            path_to_spacer_fasta = os.path.join(path_to_spacer_fasta_folder, group + '.fa')
            path_to_n_seq_file = os.path.join(path_to_spacer_fasta_folder,
                                              group + '_spacer_name_to_seq.fa')
            InputParser(os.path.join(input_path, group), path_to_spacer_fasta,
                        spacer_number_to_seq_file=path_to_n_seq_file,
                        cluster_spacers=cluster_spacers,
                        cluster_spacers_max_distance=cluster_spacers_max_distance,
                        logger=logger)
            ls_path_to_spacer_fasta.append(path_to_spacer_fasta)

    output_path_to_tree = os.path.join(output_path, 'additional_data')
    output_tree_path = os.path.join(output_path, 'additional_data', 'dict_nwk_trees.json')
    if tree_path is None:
        output_tree_path = None
    elif os.path.isfile(tree_path):
        logger.info(f'Using tree file {tree_path}.')
        tree_extension = os.path.splitext(tree_path)[1].lower()
        if tree_extension in ['.json', '.pkl', '.pickle']:
            output_tree_path = tree_path
        elif len(ls_path_to_spacer_fasta) == 1:
            dict_trees = {os.path.splitext(os.path.basename(ls_path_to_spacer_fasta[0]))[0]:
                              import_data.load_single_tree(tree_path).format('newick')}
            os.makedirs(output_path_to_tree, exist_ok=True)
            with open(output_tree_path, 'w', encoding='utf-8') as tree_file:
                json.dump(dict_trees, tree_file)
        else:
            raise ValueError('A single Newick tree file can only be used with one input group')
    else:
        logger.info(f'Using tree folder {tree_path}.')
        dict_trees = {}
        tree_files_by_stem = {}
        for filename in sorted(os.listdir(tree_path)):
            candidate = os.path.join(tree_path, filename)
            if os.path.isfile(candidate) and not filename.lower().endswith(('.fa', '.fasta', '.fna')):
                tree_files_by_stem.setdefault(os.path.splitext(filename)[0], []).append(candidate)
        for fasta_path in ls_path_to_spacer_fasta:
            group_no_ext = os.path.splitext(os.path.basename(fasta_path))[0]
            matches = tree_files_by_stem.get(group_no_ext, [])
            if len(matches) != 1:
                raise ValueError(
                    f'Expected exactly one tree named {group_no_ext} (+ extension) in {tree_path}; '
                    f'found {len(matches)}.'
                )
            dict_trees[group_no_ext] = import_data.load_single_tree(matches[0]).format('newick')
        os.makedirs(output_path_to_tree, exist_ok=True)
        with open(output_tree_path, 'w', encoding='utf-8') as tree_file:
            json.dump(dict_trees, tree_file)
    return ls_path_to_spacer_fasta, output_tree_path

def main(argv=None, prog='spacerplacer', description=None):
    args = parse_args(argv, prog=prog, description=description)

    input_path, output_path = validate_paths_and_security(args)
    args.input_path = str(input_path)
    args.output_path = str(output_path)
    if args.tree_path is not None:
        args.tree_path = str(Path(args.tree_path).expanduser().resolve())
    if args.aligned_groups_path is not None:
        aligned_groups_path = Path(args.aligned_groups_path).expanduser().resolve()
        if not aligned_groups_path.is_file():
            raise FileNotFoundError(f'Aligned-groups JSON does not exist: {aligned_groups_path}')
        args.aligned_groups_path = str(aligned_groups_path)

    if args.verbosity == 0:
        logging_level = 'WARNING'
    elif args.verbosity == 1:
        logging_level = 'INFO'
    else:
        logging_level = 'DEBUG'
    logfile_path = os.path.join(args.output_path, '0_logger.log')

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    logger = create_logger('spacerplacer', logging_level, outfile=logfile_path)

    logger.info(f'Starting SpacerPlacer on Input path: {args.input_path}')
    # start with function for single group -> multiple groups
    logger.debug(f'Input type: {args.input_type}')

    rec_parameter_dict = {'model': args.rec_model, 'gain_rate': args.insertion_rate, 'loss_rate': args.deletion_rate,
                          }

    logger.debug(vars(args))

    figsize_rec = [float(args.figsize_rec[0]) if args.figsize_rec[0] is not None else None,
                   float(args.figsize_rec[1]) if args.figsize_rec[1] is not None else None,
                   args.figsize_rec[2]]
    if args.figsize_rec[2] not in ['px', 'mm', 'in']:
        logger.error(f'Figure size unit {args.figsize_rec[2]} not recognized. '
                     f'Please use "px", "mm" or "in".')
        raise ValueError(f'Figure size unit {args.figsize_rec[2]} not recognized. '
                         f'Please use "px", "mm" or "in".')

    if args.input_type == 'pickled':
        run_pickled_data(rec_parameter_dict,
                                              lh_fct=args.lh_fct,
                                              data_path=args.input_path,
                                              save_path=args.output_path,
                                              plot_tree=not args.no_plot_reconstruction,
                                              logfile_path=os.path.join(args.output_path, '0_logger.log'),
                                              plot_order=not args.no_plot_order_graph,
                                              combine_non_unique_arrays=args.combine_non_unique_arrays,
                                              tree_path=args.tree_path,
                                              tree_save_path=args.output_path,
                                              hide_unobserved_spacers=False,
                                              selection_fct=None,
                                              significance_level=args.significance_level,
                                              finite_bl_at_root=True,
                                              alignment_options_dict=None,
                                              group_by=None,
                                              extend_branches=args.extend_branches,
                                              tree_distance_function=args.tree_distance_function,
                                              tree_lh_fct=args.tree_lh_fct,
                                              tree_construction_method=args.tree_construction_method,
                                              tree_gain_rate=args.tree_insertion_rate,
                                              tree_loss_rate=args.tree_deletion_rate,
                                              tree_alpha=args.tree_alpha,
                                              determine_orientation=args.determine_orientation,
                                              orient_boundary=args.orientation_decision_boundary,
                                              alpha_bias_correction=not args.no_alpha_bias_correction,
                                              rho_bias_correction=not args.no_rho_bias_correction,
                                              rho_bias_failure_policy=args.rho_bias_failure_policy,
                                              core_genome_trees=True if args.tree_path is not None else False,
                                              save_reconstructed_events=args.save_reconstructed_events,
                                              dpi=args.dpi_rec,
                                              figsize_rec=figsize_rec,
                                              show_spacer_name_row=args.show_spacer_name_row,
                                              )
    elif args.input_type in ['ccdb', 'crisprcasdb', 'spacer_fasta', 'sf', 'ccf_json', 'crisprcasfinder_json']:
        ls_path_to_spacer_fasta, tree_path = check_and_parse_input_data(args.input_type, args.input_path,
                                                                        args.output_path, args.tree_path, logger,
                                                                        cluster_spacers=args.cluster_spacers,
                                                                        cluster_spacers_max_distance=args.cluster_spacers_max_distance,
                                                                        min_evidence_level=args.min_evidence_level)
        df_results_wo_details = run_multiple_groups(ls_path_to_spacer_fasta, args.output_path,
                                                    rec_parameter_dict=rec_parameter_dict,
                                                        logger=logger,
                                                    lh_fct=args.lh_fct,
                                                    plot_tree=not args.no_plot_reconstruction,
                                                    do_show=args.do_show,
                                                    determine_orientation=args.determine_orientation,
                                                    orientation_decision_boundary=args.orientation_decision_boundary,
                                                    tree_path=tree_path,
                                                    plot_order=not args.no_plot_order_graph,
                                                    significance_level=args.significance_level,
                                                    extend_branches=args.extend_branches,
                                                    tree_distance_fct=args.tree_distance_function,
                                                    tree_construction_method=args.tree_construction_method,
                                                    tree_lh_fct=args.tree_lh_fct,
                                                    tree_insertion_rate=args.tree_insertion_rate,
                                                    tree_deletion_rate=args.tree_deletion_rate,
                                                    tree_alpha=args.tree_alpha,
                                                    alpha_bias_correction=not args.no_alpha_bias_correction,
                                                    rho_bias_correction=not args.no_rho_bias_correction,
                                                    rho_bias_failure_policy=args.rho_bias_failure_policy,
                                                    combine_non_unique_arrays=args.combine_non_unique_arrays,
                                                    seed=args.seed,
                                                    save_reconstructed_events=args.save_reconstructed_events,
                                                    dpi=args.dpi_rec,
                                                    figsize_rec=figsize_rec,
                                                    show_spacer_name_row=args.show_spacer_name_row,
                                                    alignment_workers=args.alignment_workers,
                                                    aligned_groups_path=args.aligned_groups_path,
                                                    )
        summary_dict = compose_summary_dict(df_results_wo_details, dict(vars(args)))
        write_summary(summary_dict, os.path.join(args.output_path, 'summary.txt'))

    else:
        message = f'Input type {args.input_type} not recognized.'
        logger.error(message)
        raise ValueError(message)

    return args

if __name__ == '__main__':
    main()
