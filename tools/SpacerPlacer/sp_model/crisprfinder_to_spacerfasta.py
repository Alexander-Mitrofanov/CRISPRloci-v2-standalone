import json
import os.path

from sp_model.helpers.misc import create_logger
from sp_model.helpers.sequence import cluster_by_edit_distance, cluster_by_overlap


def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def reverse_complement(dna_seq):
    """Backward-compatible wrapper around the shared IUPAC-aware helper."""
    from sp_model.helpers.sequence import reverse_complement as _reverse_complement

    return _reverse_complement(dna_seq)


class CRISPRFinderParser:
    def __init__(self, crisprfinder_result_json_file, output_folder=None, output_folder_for_too_small_groups=None,
                 min_evidence_level=4,
                 reorient_by_direction_pred=True, logger=None,
                 cluster_by_spacer_overlap=True, filter_more_than_x_in_group=None, max_levenshtein_distance_spacers=0):
        if not 1 <= min_evidence_level <= 4:
            raise ValueError("min_evidence_level must be between 1 and 4")
        if max_levenshtein_distance_spacers < 0:
            raise ValueError("max_levenshtein_distance_spacers must be non-negative")
        self.crisprfinder_result_json_file = crisprfinder_result_json_file
        self.min_evidence_level = min_evidence_level
        self.reorient_by_direction_pred = reorient_by_direction_pred
        self.max_levenshtein_distance_spacers = max_levenshtein_distance_spacers
        self.crispr_dictionary = {}
        self.clustered_by_repeat_arrays = {}
        self.clustered_groups = {}
        self.clustered_groups_spacer_ids = {}
        self.cluster_spacer_seq_to_spacer_id = {}
        self.spacer_id_to_clustered_spacer_seq_for_each_group = {}
        self.spacer_cluster = None
        self.sample_summary = {}
        self.logger = logger if logger is not None else create_logger('CRISPRFinderParser', 1)
        self.cluster_by_spacer_overlap = cluster_by_spacer_overlap
        self.filter_more_than_x_in_group = filter_more_than_x_in_group
        self.output_folder = output_folder
        self.output_folder_for_too_small_groups = output_folder_for_too_small_groups

        self._parse_ls_crisprfinder_files()
        self._cluster_array_groups_by_repeat()
        if max_levenshtein_distance_spacers > 0:
            self._cluster_spacers_by_levenshtein_distance()
        if self.cluster_by_spacer_overlap:
            self._cluster_groups_by_spacer_overlap()
        self._distribute_spacer_ids()
        if self.output_folder is not None:
            os.makedirs(self.output_folder, exist_ok=True)
            if self.output_folder_for_too_small_groups is not None:
                os.makedirs(self.output_folder_for_too_small_groups, exist_ok=True)
            self.write_spacer_fasta(
                self.output_folder, filter_more_than_x_in_group=self.filter_more_than_x_in_group
            )
            # self.write_spacer_seq_to_spacer_id(self.output_folder)


    def write_spacer_fasta(self, output_folder, filter_more_than_x_in_group=None):
        self.logger.info(f'Writing spacer sequences to FASTA files in {output_folder}')
        for group_name, ls_array_tuples in self.clustered_groups_spacer_ids.items():
            of = output_folder
            if filter_more_than_x_in_group is not None:
                if len(ls_array_tuples) < filter_more_than_x_in_group:
                    self.logger.info(f'Group {group_name} has fewer than {filter_more_than_x_in_group} and is saved to {self.output_folder_for_too_small_groups}.')
                    of = self.output_folder_for_too_small_groups
                    if of is None:
                        self.logger.warning('Skipping undersized group %s because no small-group output folder was set.',
                                            group_name)
                        continue
            os.makedirs(of, exist_ok=True)
            with open(f'{of}/{group_name}.fa', 'w') as f:
                for array_name, spacers, direction in ls_array_tuples:
                    if not array_name or any(character in array_name for character in '\r\n'):
                        raise ValueError(f'Invalid FASTA array name: {array_name!r}')
                    f.write(f'>{array_name} \n')
                    f.write(f'{", ".join([str(s) for s in spacers])} \n')

    def write_spacer_seq_to_spacer_id(self, output_folder):
        self.logger.info(f'Writing dictionary of spacer sequences -> spacer ids to {output_folder}')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for group_name, dict_spacer_seq_to_spacer_id in self.cluster_spacer_seq_to_spacer_id.items():
            with open(f'{output_folder}/{group_name}_spacer_seq_to_spacer_id.json', 'w') as f:
                json.dump(dict_spacer_seq_to_spacer_id, f)

    def write_spacer_cluster_to_file(self, output_folder):
        """Saves spacer clusters to file. The first spacer in each cluster is the representative spacer."""
        self.logger.info(f'Writing spacer clusters to {output_folder}')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if self.spacer_cluster is not None:
            with open(f'{output_folder}/spacer_clusters.json', 'w') as f:
                json.dump(self.spacer_cluster, f, default=set_default)

    def get_spacer_cluster(self):
        return self.spacer_cluster

    def get_clustered_groups_spacer_ids(self):
        return self.clustered_groups_spacer_ids

    def get_cluster_spacer_seq_to_spacer_id(self):
        return self.cluster_spacer_seq_to_spacer_id

    def get_cluster_spacer_id_to_spacer_seq(self):
        return {
            group_name: {spacer_id: sequence for sequence, spacer_id in mapping.items()}
            for group_name, mapping in self.cluster_spacer_seq_to_spacer_id.items()
        }

    def _parse_ls_crisprfinder_files(self):
        if isinstance(self.crisprfinder_result_json_file, list):
            for file_path in self.crisprfinder_result_json_file:
                self._parse_crisprfinder_file(file_path)
        else:
            self._parse_crisprfinder_file(self.crisprfinder_result_json_file)

    def _parse_crisprfinder_file(self, file_path):
        self.logger.info(f'Parsing CRISPRFinder JSON file: {file_path}')
        with open(file_path, 'r') as f:
            data = json.load(f)
            ls_results = data.get('Sequences')
            if not isinstance(ls_results, list):
                raise ValueError(f"CRISPRCasFinder file {file_path} must contain a 'Sequences' list")

            for result in ls_results:
                if not isinstance(result, dict) or 'Id' not in result:
                    raise ValueError(f'Each sequence entry in {file_path} must be an object with an Id')
                sequence_id = str(result['Id'])
                sequence_length = result.get('Length', 'unknown')
                self.logger.info(f'Parsing CRISPRFinder result for sequence {sequence_id} with length {sequence_length}, '
                                 f'skipping CRISPRs with evidence level < {self.min_evidence_level}:')

                sequence_crisprs = result.get('Crisprs', [])
                if len(sequence_crisprs) == 0:
                    self.logger.info(f'No CRISPRs found for sequence {sequence_id}')
                    continue
                else:
                    self.logger.info(f'Found {len(sequence_crisprs)} CRISPR(s) for sequence {sequence_id}')

                for crispr in sequence_crisprs:
                    required_fields = {'Evidence_Level', 'Name', 'Start', 'End', 'DR_Consensus',
                                       'CRISPRDirection', 'Regions'}
                    missing_fields = required_fields - crispr.keys()
                    if missing_fields:
                        raise ValueError(f'CRISPR entry in {file_path} is missing {sorted(missing_fields)}')
                    if crispr['Evidence_Level'] < self.min_evidence_level:
                        self.logger.debug(f'Skipping CRISPR with evidence level {crispr["Evidence_Level"]}')
                        continue
                    crispr_name = str(crispr['Name'])
                    if crispr_name in self.crispr_dictionary:
                        base_name = f"{sequence_id}::{crispr_name}"
                        crispr_name = base_name
                        suffix = 2
                        while crispr_name in self.crispr_dictionary:
                            crispr_name = f"{base_name}::{suffix}"
                            suffix += 1
                    crispr_start = crispr['Start']
                    crispr_end = crispr['End']

                    self.logger.info(f'Parsing CRISPR array {crispr_name} with {crispr.get("Spacers", "unknown")} spacers,'
                                     f'consensus repeat: {crispr["DR_Consensus"]}, '
                                     f'direction: {crispr["CRISPRDirection"]},'
                                     f'potential orientation: {crispr.get("Potential_Orientation")}.')
                    dr_consensus = str(crispr['DR_Consensus']).upper()
                    invalid_repeat = set(dr_consensus) - set('ACGTURYSWKMBDHVN')
                    if not dr_consensus or invalid_repeat:
                        raise ValueError(f'Invalid DR consensus for CRISPR {crispr_name}: {dr_consensus!r}')
                    crispr_direction = crispr['CRISPRDirection']
                    potential_orientation = crispr.get('Potential_Orientation')

                    repeats = []
                    spacers = []

                    for val in crispr['Regions']:
                        if not isinstance(val, dict) or 'Type' not in val or 'Sequence' not in val:
                            raise ValueError(f'Malformed region in CRISPR {crispr_name}')
                        if val['Type'] == 'DR':
                            repeats.append(val['Sequence'].upper())
                        elif val['Type'] == 'Spacer':
                            spacers.append(val['Sequence'].upper())

                    if not spacers:
                        raise ValueError(f'CRISPR {crispr_name} contains no spacer regions')
                    invalid_spacers = set().union(*(set(spacer) for spacer in spacers)) - set('ACGTURYSWKMBDHVN')
                    if invalid_spacers:
                        raise ValueError(f'Invalid spacer characters in CRISPR {crispr_name}: {sorted(invalid_spacers)}')

                    self.crispr_dictionary[crispr_name] = {'start': crispr_start,
                                                           'end': crispr_end,
                                                           'direction': crispr_direction,
                                                           'potential_orientation': potential_orientation,
                                                           'dr_consensus': dr_consensus,
                                                           'repeats': repeats,
                                                           'spacers': spacers}

    def _cluster_array_groups_by_repeat(self):
        """
        Clusters CRISPR arrays based on equal consensus repeat.
        :return:
        """
        self.logger.info('Clustering CRISPR arrays based on equal consensus repeat')
        for array_name, array_data in self.crispr_dictionary.items():
            direction = array_data['direction']
            dr_consensus = array_data['dr_consensus']
            spacers = array_data['spacers']
            # Not entirely sure how the output of CRISPRCasFinder is structured. For now, I assume that the spacers
            # are in the same orientation as in the input sequence (and thus repeats/arrays/spacers
            # need to be reversed for overlap search).
            if self.reorient_by_direction_pred and direction == '-':
                dr_consensus = reverse_complement(dr_consensus)
                spacers = [reverse_complement(s) for s in spacers[::-1]]
            if dr_consensus not in self.clustered_by_repeat_arrays:
                self.clustered_by_repeat_arrays[dr_consensus] = []
            self.clustered_by_repeat_arrays[dr_consensus].append((array_name, spacers,
                                                                  array_data['direction']))


    def _cluster_spacers_by_levenshtein_distance(self):
        """Cluster spacers into transitive edit-distance components."""
        self.logger.info(
            "Clustering spacers with maximum Levenshtein distance %s",
            self.max_levenshtein_distance_spacers,
        )
        self.spacer_cluster = {}
        clustered_arrays = {}
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            unique_spacers = [
                spacer
                for _, spacers, _ in array_tuples
                for spacer in spacers
            ]
            clusters = cluster_by_edit_distance(
                unique_spacers, self.max_levenshtein_distance_spacers
            )
            representative_by_spacer = {
                spacer: cluster[0]
                for cluster in clusters
                for spacer in cluster
            }
            self.spacer_cluster[repeat] = {
                cluster[0]: set(cluster) for cluster in clusters
            }
            clustered_arrays[repeat] = [
                (
                    array_name,
                    [representative_by_spacer[spacer] for spacer in spacers],
                    direction,
                )
                for array_name, spacers, direction in array_tuples
            ]
        self.clustered_by_repeat_arrays = clustered_arrays

    def _cluster_groups_by_spacer_overlap(self):
        """Cluster arrays into transitive spacer-overlap components."""
        self.logger.info("Clustering CRISPR arrays based on spacer overlap")
        counter = 0
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            if len(array_tuples) == 1:
                self.logger.info(
                    "Only one array with consensus repeat %s was found.", repeat
                )
            spacer_arrays = [spacers for _, spacers, _ in array_tuples]
            for component in cluster_by_overlap(spacer_arrays):
                group_name = f"g_{counter}"
                self.clustered_groups[group_name] = [
                    array_tuples[index] for index in component
                ]
                counter += 1

    def _distribute_spacer_ids(self):
        """
        Assigns an ID (number) to each unique spacer sequence in each group (IDs are not globally unique across groups).
        (If spacers are clustered by overlap the ID stands for the representative spacer (which in turn stands for all
        spacers in the cluster).) A dictionary is created mapping (representative) spacer sequences to spacer IDs.
        :return:
        """
        self.logger.info('Distributing spacer IDs')
        if self.cluster_by_spacer_overlap:
            groups = self.clustered_groups
        else:
            groups = self.clustered_by_repeat_arrays
        for group_name, ls_array_tuples in groups.items():
            dict_spacer_seq_to_spacer_id = {}
            spacer_counter = 0
            ls_tuples = []
            for array_name, spacers, direction in ls_array_tuples:
                spacers_as_id = []
                for spacer_seq in spacers:
                    if spacer_seq not in dict_spacer_seq_to_spacer_id:
                        spacer_id = spacer_counter
                        dict_spacer_seq_to_spacer_id[spacer_seq] = spacer_id
                        spacer_counter += 1
                    else:
                        spacer_id = dict_spacer_seq_to_spacer_id[spacer_seq]
                    spacers_as_id.append(spacer_id)
                ls_tuples.append((array_name, spacers_as_id, direction))
            self.clustered_groups_spacer_ids[group_name] = ls_tuples
            self.cluster_spacer_seq_to_spacer_id[group_name] = dict_spacer_seq_to_spacer_id

            # TODO implement dictionary spacer_id -> representative spacer seq -> all spacers in cluster
            # This is missing the repeat information, could use repeat as group name (+ additional counter for uniqueness)?
            # dict_spacer_id_to_spacer_clustered_seq = {idy: self.spacer_cluster[repeat][seq] for seq, idy in dict_spacer_seq_to_spacer_id.items()}
            # self.spacer_id_to_clustered_spacer_seq_for_each_group[group_name] = dict_spacer_id_to_spacer_clustered_seq
