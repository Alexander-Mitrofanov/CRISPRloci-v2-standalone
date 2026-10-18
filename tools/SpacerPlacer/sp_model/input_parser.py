import os

from sp_model.helpers.sequence import cluster_by_edit_distance, reverse_complement


def levenshtein_distance(seq1, seq2):
    """Backward-compatible wrapper around the optimized shared helper."""
    from sp_model.helpers.sequence import levenshtein_distance as _distance

    return _distance(seq1, seq2)

class InputParser:
    def __init__(self, input_folder, restructured_output_path, spacer_number_to_seq_file=None,
                 cluster_spacers=False, check_orientation_for_highest_spacer_overlap=False,
                 split_into_groups_by_spacer_overlap=False, cluster_spacers_max_distance=0, logger=None):
        self.input_folder = input_folder
        self.restructured_output_path = restructured_output_path
        self.spacer_number_to_seq_file = spacer_number_to_seq_file
        self.cluster_spacers_max_distance = cluster_spacers_max_distance
        self.spacer_seq_to_number = None
        # self.spacer_number_to_seq = None
        self.file_name_to_spacer_indices = None
        self.spacer_orientations = None
        self.logger = logger

        self.cluster_spacers = cluster_spacers
        self.check_orientation_for_highest_spacer_overlap = check_orientation_for_highest_spacer_overlap
        self.split_into_groups_by_spacer_overlap = split_into_groups_by_spacer_overlap

        self._enumerate_spacers()
        self.clustering_steps()
        self.write_output_to_file()

    def _enumerate_spacers(self):
        spacers = {}
        spacers_ori = {}
        spacer_count = 0
        result = {}
        records = {}

        forward_orientation = os.path.join(self.input_folder, "pos_strand")
        reverse_orientation = os.path.join(self.input_folder, "neg_strand")
        orientation_folders = [
            (forward_orientation, False, "pos_strand"),
            (reverse_orientation, True, "neg_strand"),
        ]
        existing_folders = [folder for folder, _, _ in orientation_folders if os.path.isdir(folder)]
        if not existing_folders:
            message = (
                "No pos_strand or neg_strand folder found in input group "
                f"'{self.input_folder}'."
            )
            if self.logger:
                self.logger.error(message)
            raise FileNotFoundError(message)

        for folder, should_reverse, orientation in orientation_folders:
            if not os.path.isdir(folder):
                if self.logger:
                    self.logger.info("Optional input folder '%s' is absent.", folder)
                continue
            for file_name in sorted(os.listdir(folder)):
                file_path = os.path.join(folder, file_name)
                if not os.path.isfile(file_path):
                    continue
                array_name = os.path.splitext(file_name)[0]
                if not array_name:
                    raise ValueError(f"Input filename has no usable array name: {file_path}")
                with open(file_path, "r", encoding="utf-8") as file:
                    raw_sequences = [
                        line.strip().upper()
                        for line in file
                        if line.strip() and not line.startswith(">")
                    ]
                if not raw_sequences:
                    raise ValueError(f"No spacer sequences found in {file_path}")
                invalid = set().union(*(set(spacer) for spacer in raw_sequences)) - set("ACGTURYSWKMBDHVN")
                if invalid:
                    raise ValueError(f"Invalid DNA character(s) {sorted(invalid)} in {file_path}")
                if array_name in records:
                    previous_sequences, _, previous_orientation, previous_path = records[array_name]
                    if raw_sequences != previous_sequences or previous_orientation == orientation:
                        raise ValueError(
                            f"Conflicting duplicate array name '{array_name}' in {previous_path} and {file_path}."
                        )
                    if self.logger:
                        self.logger.warning(
                            "Array '%s' occurs identically in both orientation folders; using %s from %s.",
                            array_name, orientation, file_path,
                        )
                records[array_name] = (raw_sequences, should_reverse, orientation, file_path)

        for array_name, (raw_sequences, should_reverse, orientation, _) in records.items():
            sequences = (
                [reverse_complement(sequence) for sequence in reversed(raw_sequences)]
                if should_reverse else raw_sequences
            )
            result[array_name] = []
            for spacer in sequences:
                if spacer not in spacers:
                    spacers[spacer] = spacer_count
                    spacers_ori[spacer] = orientation
                    spacer_count += 1
                result[array_name].append(spacers[spacer])
        self.spacer_orientations = spacers_ori
        self.spacer_seq_to_number = spacers
        # self.spacer_number_to_seq = {v: k for k, v in spacers.items()}
        self.file_name_to_spacer_indices = result
        return

    def write_output_to_file(self):
        with open(self.restructured_output_path, 'w') as file:
            for file_name, spacer_indices in self.file_name_to_spacer_indices.items():
                file.write(f">{file_name}\n")
                file.write(", ".join(f"{index}" for index in spacer_indices))
                file.write("\n")

        if self.spacer_number_to_seq_file is not None:
            with open(self.spacer_number_to_seq_file, 'w') as file:
                for spacer_seq, spacer_number in self.spacer_seq_to_number.items():
                    file.write(f">{spacer_number}, {self.spacer_orientations[spacer_seq]}\n")
                    file.write(f"{spacer_seq}")
                    file.write("\n")

    def clustering_steps(self):
        """
        :return:
        """
        if self.check_orientation_for_highest_spacer_overlap:
            # Accept the orientation as is, as a start. Then check for each array, if reversing would increase overlap.
            # If so, reverse the array. If not, keep as is. Run whole ls of arrays this way one after the other.
            # Could be problematic, if the arrays are split in orientation. But probably most efficient.
            # Can cluster spacers in each step. But probably sufficient to cluster at the end.
            raise NotImplementedError
        if self.cluster_spacers:
            # cluster spacers with levenstein distance. Start with single clusters. Then merge clusters with some
            # distance (x) in mutations iteratively, until no clusters are mergeable.
            new_clusters = self.cluster_spacers_with_levenshtein()
            self.rename_spacers(new_clusters)
        if self.split_into_groups_by_spacer_overlap:
            # Split arrays, if they have no overlap in spacers. Need extra cases for no overlap. Throw warning and then
            # main whole group (w/o split)? Or main only overlapping groups and have additional group without overlap?
            raise NotImplementedError
        return

    def cluster_spacers_with_levenshtein(self):
        return cluster_by_edit_distance(
            self.spacer_seq_to_number.keys(), self.cluster_spacers_max_distance
        )

    def rename_spacers(self, clusters):
        # Create new numbering by identifying each cluster with a number
        old_numbering = {spacer: number for number, spacer in self.spacer_seq_to_number.items()}

        new_numbering = {}
        for cluster_number, cluster in enumerate(clusters):
            for spacer in cluster:
                new_numbering[spacer] = cluster_number
        self.file_name_to_spacer_indices = {file_name: [new_numbering[old_numbering[spacer]] for spacer in spacers]
                                            for file_name, spacers in self.file_name_to_spacer_indices.items()}
        self.spacer_seq_to_number = new_numbering


    @staticmethod
    def _rev_com(sequence):
        return reverse_complement(sequence)
