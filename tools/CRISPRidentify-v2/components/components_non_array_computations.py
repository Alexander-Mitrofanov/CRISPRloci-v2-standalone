import os
import subprocess
import shutil
from os.path import isfile
import math

from components.coordinate_utils import crispr_candidate_span
from components.coordinate_utils import is_forward_strand
from components.reconstruction_validator import ReconstructionValidationError
from components.reconstruction_validator import validate_candidate_reconstruction
from components.portable_strand import load_default_portable_strand_model
from components.sequence_utils import reverse_complement


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve_executable(executable_name, *project_relative_candidates):
    """Resolve an executable or raise instead of returning an empty result."""
    candidates = [
        os.path.join(PROJECT_ROOT, candidate)
        for candidate in project_relative_candidates
    ]
    environment_candidate = shutil.which(executable_name)
    if environment_candidate:
        candidates.append(environment_candidate)
    for candidate in candidates:
        if isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError(
        "required executable was not found: {}".format(executable_name)
    )


def resolve_project_file(path):
    """Resolve a required resource against the job cwd and project root."""
    candidates = [path]
    if not os.path.isabs(path):
        candidates.append(os.path.join(PROJECT_ROOT, path))
    for candidate in candidates:
        if isfile(candidate):
            return os.path.abspath(candidate)
    raise FileNotFoundError("required resource was not found: {}".format(path))

#              Strand computation
#######################################################
#######################################################


class StrandComputation:
    def __init__(self, list_of_crisprs):
        self.list_of_crisprs = list_of_crisprs
        self.dict_strands = {}

        self._compute_all_strands()

    def _compute_all_strands(self):
        if not self.list_of_crisprs:
            return

        model = load_default_portable_strand_model()
        strands = model.predict_strands(
            [crispr.consensus for crispr in self.list_of_crisprs]
        )
        if len(strands) != len(self.list_of_crisprs):
            raise ValueError(
                "portable CRISPRstrand returned {} predictions for {} arrays"
                .format(len(strands), len(self.list_of_crisprs))
            )
        self.dict_strands.update(enumerate(strands))

    def output(self):
        return self.dict_strands



#      IS element computation
################################################
################################################


class FastaMatch:
    def __init__(self, first_id, second_id, similarity, coverage):
        self.first_id = int(first_id)
        self.second_id = int(second_id)
        self.similarity = similarity
        self.coverage = coverage


class FastaSimilarity():
    def __init__(self, list_sequences, similarity, coverage):
        self.list_sequences = list_sequences
        self.similarity = similarity
        self.coverage = coverage

        self.fasta_info = None

        self._run_fasta()
        self._parse_fasta()

        self._filter_by_similarity()
        self._filter_by_coverage()

    def _run_fasta(self):
        with open("sequences.fa", "w") as f:
            for index, repeat in enumerate(self.list_sequences):
                f.write(">{}\n".format(index))
                f.write("{}\n".format(repeat.strip()))
        try:
            executable = resolve_executable(
                "fasta36", "tools/fasta/fasta36"
            )
            with open("fasta_similarity.fastab", "w") as output_file:
                subprocess.run(
                    [executable, "sequences.fa", "sequences.fa", "-m", "8"],
                    stdout=output_file,
                    stderr=subprocess.PIPE,
                    check=True,
                    universal_newlines=True,
                )
        except Exception:
            try:
                os.remove("fasta_similarity.fastab")
            except OSError:
                pass
            raise
        finally:
            try:
                os.remove("sequences.fa")
            except OSError:
                pass

    def _parse_fasta(self):
        with open("fasta_similarity.fastab", "r") as f:
            lines = f.readlines()

        self.fasta_info = [
            line.strip().split("\t") for line in lines if line.strip()
        ]
        malformed = [line for line in self.fasta_info if len(line) < 10]
        if malformed:
            raise ValueError("malformed FASTA36 tabular output")
        self.fasta_info = [fasta_line_match for fasta_line_match in self.fasta_info
                           if int(fasta_line_match[0]) < int(fasta_line_match[1])]

        os.remove("fasta_similarity.fastab")

    def _filter_by_similarity(self):
        self.fasta_info = [fasta_line_match for fasta_line_match in self.fasta_info
                           if float(fasta_line_match[2]) / 100 >= self.similarity - 0.00005]

    def _filter_by_coverage(self):
        filtered_by_coverage = []
        for fasta_line_match in self.fasta_info:
            interval_first_element = abs(int(fasta_line_match[6]) - int(fasta_line_match[7])) + 1
            interval_second_element = abs(int(fasta_line_match[8]) - int(fasta_line_match[9])) + 1

            coverage_first = float(interval_first_element) / len(self.list_sequences[int(fasta_line_match[0])])
            coverage_second = float(interval_second_element) / len(self.list_sequences[int(fasta_line_match[1])])
            min_coverage = min(coverage_first, coverage_second)
            if min_coverage >= self.coverage - 0.00005:
                min_coverage = int(min_coverage * 100)
                fasta_line_match.append(min_coverage)
                filtered_by_coverage.append(fasta_line_match)

        self.fasta_info = filtered_by_coverage

    def output(self):
        list_fasta_info_objects = []
        for line in self.fasta_info:
            fasta_line_match_important = [line[i] for i in (0, 1, 2, -1)]
            fasta_match = FastaMatch(*fasta_line_match_important)
            list_fasta_info_objects.append(fasta_match)
        return list_fasta_info_objects


class HMMMatch:
    def __init__(self, target, query, e_value, score):
        self.target = target
        self.query = query
        self.e_value = e_value
        self.score = score

    def __repr__(self):
        return "{} {} {} {}".format(self.target, self.query, self.e_value, self.score)


class HMMMatchProteinCoordinates:
    def __init__(self, target, query,  e_value, score, protein_start, protein_end, protein_strand):
        self.target = target
        self.query = query
        self.e_value = e_value
        self.score = score
        self.protein_start = protein_start
        self.protein_end = protein_end
        self.protein_strand = protein_strand

    def __repr__(self):
        return "{} {} {} {} {} {} {}".format(self.target, self.query, self.e_value,
                                             self.score, self.protein_start,
                                             self.protein_end, self.protein_strand)


class HMMResultParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.list_matches = []
        self._input_parsing()

    def _input_parsing(self):
        with open(self.file_path, "r") as f:
            for line_number, line in enumerate(f, 1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                fields = stripped.split()
                if len(fields) < 6:
                    raise ValueError(
                        "malformed HMMER table row {} in {}: expected at "
                        "least 6 columns".format(line_number, self.file_path)
                    )
                try:
                    e_value = float(fields[4])
                    score = float(fields[5])
                except ValueError as error:
                    raise ValueError(
                        "malformed HMMER numeric value on row {} in {}"
                        .format(line_number, self.file_path)
                    ) from error
                if not math.isfinite(e_value) or not math.isfinite(score):
                    raise ValueError(
                        "non-finite HMMER numeric value on row {} in {}"
                        .format(line_number, self.file_path)
                    )

                self.list_matches.append(
                    HMMMatch(fields[0], fields[2], e_value, score)
                )

    def output(self):
        return self.list_matches


class ISElementSearch:
    def __init__(self, dna_seq, hmm_model):
        self.dna_seq = dna_seq
        self.hmm_model = hmm_model

        self.dict_protein_description = {}
        self.hmm_matches = []
        self.hmm_protein_match = []

        self._protein_extraction()
        self._hmm_extraction()
        self._hmm_match_protein_creation()

    def _protein_extraction(self):
        try:
            os.remove("protein.fa")
        except OSError:
            pass
        with open('dna.fa', 'w') as f:
            f.write('>seq_to_score\n')
            f.write(self.dna_seq)

        try:
            executable = resolve_executable(
                "prodigal", "tools/prodigal/prodigal"
            )
            subprocess.run(
                [
                    executable, "-i", "dna.fa", "-p", "meta",
                    "-a", "protein.fa",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                universal_newlines=True,
            )
        except Exception:
            try:
                os.remove("protein.fa")
            except OSError:
                pass
            raise
        finally:
            try:
                os.remove("dna.fa")
            except OSError:
                pass

        with open("protein.fa", "r") as f:
            lines = f.readlines()
        header_lines = [line[1:] for line in lines if ">" in line]
        keys = [line.split()[0] for line in header_lines]
        coordinates = [(int(line.split()[2]), int(line.split()[4])) for line in header_lines]
        strands = [int(line.split()[6]) for line in header_lines]
        strands = ["Forward" if x == 1 else "Reverse" for x in strands]
        values = [(*coordinate, strand) for coordinate, strand in zip(coordinates, strands)]

        self.dict_protein_description = {key: value for key, value in zip(keys, values)}

    def _hmm_extraction(self):
        try:
            if os.stat("protein.fa").st_size != 0:
                executable = resolve_executable(
                    "hmmsearch", "tools/hmm_search/hmmsearch"
                )
                hmm_model = resolve_project_file(self.hmm_model)
                subprocess.run(
                    [
                        executable, "--tblout", "result_hmm.out",
                        hmm_model, "protein.fa",
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True,
                    universal_newlines=True,
                )

                hmm_result_parser = HMMResultParser("result_hmm.out")
                self.hmm_matches = hmm_result_parser.output()
        finally:
            for path in ("result_hmm.out", "protein.fa"):
                try:
                    os.remove(path)
                except OSError:
                    pass

    def _hmm_match_protein_creation(self):
        for hmm_match in self.hmm_matches:
            hmm_target = hmm_match.target
            hmm_query = hmm_match.query
            hmm_e_value = hmm_match.e_value
            hmm_score = hmm_match.score
            protein_start, protein_end, protein_strand = self.dict_protein_description[hmm_target]
            self.hmm_protein_match.append(HMMMatchProteinCoordinates(hmm_target, hmm_query, hmm_e_value,
                                                                     hmm_score, protein_start, protein_end,
                                                                     protein_strand))

    def output(self):
        return self.hmm_protein_match


def project_local_prodigal_interval(
        segment_start_zero_based, segment_length, local_start, local_end):
    """Project a local one-based inclusive Prodigal interval to the genome."""
    if segment_start_zero_based < 0 or segment_length < 0:
        raise ValueError("invalid zero-based source segment")
    if local_start < 1 or local_end < local_start or local_end > segment_length:
        raise ValueError(
            "Prodigal interval {}..{} is outside a {}-base source segment"
            .format(local_start, local_end, segment_length)
        )
    return (
        segment_start_zero_based + local_start,
        segment_start_zero_based + local_end,
    )


class FullISElementSearch:
    def __init__(self, full_dna, list_of_crisprs, hmm_model, min_similarity, min_coverage):
        self.full_dna = full_dna
        self.list_of_crisprs = list_of_crisprs
        self.hmm_model = hmm_model
        self.min_similarity = min_similarity
        self.min_coverage = min_coverage

        self.crispr_pairs = []
        self.is_results = {}
        self.consensus_repeats = []
        self.source_ordered_crisprs = []

        self._compute_all_is_elements()

    def _compute_all_is_elements(self):
        if self.list_of_crisprs:
            self._prepare_source_order()
        if len(self.list_of_crisprs) < 2:
            return
        self._get_consensus_repeats()
        self._get_crispr_pairs()
        self._perform_hmm_search()

    def _prepare_source_order(self):
        source_records = []
        for input_index, candidate in enumerate(self.list_of_crisprs):
            reconstruction = validate_candidate_reconstruction(
                self.full_dna, candidate
            )
            if not reconstruction.valid:
                raise ReconstructionValidationError(reconstruction)
            source_records.append(
                (
                    reconstruction.array_start,
                    reconstruction.array_end,
                    input_index,
                    candidate,
                )
            )

        self.source_ordered_crisprs = sorted(
            source_records,
            key=lambda record: (record[0], record[1], record[2]),
        )

    def _get_consensus_repeats(self):
        self.consensus_repeats = [
            record[3].consensus for record in self.source_ordered_crisprs
        ]

    def _get_crispr_pairs(self):
        fasta_similarity = FastaSimilarity(self.consensus_repeats, self.min_similarity, self.min_coverage)
        fasta_results = fasta_similarity.output()
        for fasta_match in fasta_results:
            first_id, second_id = fasta_match.first_id, fasta_match.second_id
            if second_id - first_id == 1:
                cr_pair = (first_id, second_id)
                self.crispr_pairs.append(cr_pair)

    def _perform_hmm_search(self):
        for cr_pair in self.crispr_pairs:
            first_index, second_index = cr_pair[0], cr_pair[1]
            first_record = self.source_ordered_crisprs[first_index]
            second_record = self.source_ordered_crisprs[second_index]
            end_first = first_record[1]
            begin_second = second_record[0]

            if begin_second < end_first:
                raise ValueError(
                    "IS search requires source-ordered, non-overlapping arrays"
                )

            dna_interval = self.full_dna[end_first:begin_second]
            is_element_search = ISElementSearch(dna_interval, self.hmm_model)
            is_search_result = is_element_search.output()
            sorted_is_search_result = sorted(is_search_result, key=lambda x: x.e_value)
            if sorted_is_search_result:
                best_hit = sorted_is_search_result[0]
                best_interval_start = best_hit.protein_start
                best_interval_end = best_hit.protein_end
                best_interval_strand = best_hit.protein_strand
                best_interval_target = best_hit.target
                best_interval_query = best_hit.query
                genomic_start, genomic_end = project_local_prodigal_interval(
                    end_first,
                    len(dna_interval),
                    best_interval_start,
                    best_interval_end,
                )
                original_input_index = first_record[2]
                self.is_results[original_input_index] = (
                    genomic_start,
                    genomic_end,
                    best_interval_strand,
                    best_interval_target,
                    best_interval_query,
                )

    def output(self):
        return self.is_results


#        Leader seq search
##########################################################
##########################################################


def rev_compliment_seq(seq):
    """Compatibility wrapper for the historically misspelled public helper."""

    return reverse_complement(seq)


class FullLeaderSeqSearch:
    def __init__(self, list_crispr_candidates, list_corresponding_strands, full_dna):
        self.list_crispr_candidates = list_crispr_candidates
        self.list_strands = list_corresponding_strands
        self.full_dna = full_dna

        self.dict_leader_data = {}
        self.dict_downstream_data = {}

        self._compute_all_leaders()

    def _compute_all_leaders(self):
        if len(self.list_strands) != len(self.list_crispr_candidates):
            raise ValueError("every CRISPR candidate must have a strand annotation")

        for index, crispr_candidate in enumerate(self.list_crispr_candidates):
            strand = self.list_strands[index]
            leader, downstream = LeaderSeqSearch(crispr_candidate, strand, self.full_dna).output()
            self.dict_leader_data[index] = leader
            self.dict_downstream_data[index] = downstream

    def output(self):
        return self.dict_leader_data, self.dict_downstream_data


class LeaderSeqSearch:
    def __init__(self, crispr_candidate, strand, full_dna):
        self.crispr_candidate = crispr_candidate
        self.strand = strand
        self.full_dna = full_dna

        self.leader = None
        self.downstream_region = None

        self._compute_leader_seq()

    def _compute_leader_seq(self):
        array_start, array_end = crispr_candidate_span(self.crispr_candidate)
        sequence_length = len(self.full_dna)

        upstream = self.full_dna[max(0, array_start - 200):array_start]
        downstream = self.full_dna[array_end:min(sequence_length, array_end + 200)]

        if is_forward_strand(self.strand):
            self.leader = upstream
            self.downstream_region = downstream
        else:
            self.leader = rev_compliment_seq(downstream)
            self.downstream_region = rev_compliment_seq(upstream)

    def output(self):
        return self.leader, self.downstream_region
