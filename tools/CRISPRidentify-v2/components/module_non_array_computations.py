from components.components_non_array_computations import FullISElementSearch
from components.components_non_array_computations import FullLeaderSeqSearch
from components.components_non_array_computations import StrandComputation
from components.casandra_adapter import CasAndraPredictor
from components.reconstruction_validator import read_single_fasta


class NonArrayComputations:
    def __init__(self, file_path, categories, flags_non_arrays_computations):
        self.file_path = file_path
        self.categories = categories
        self.flags_non_arrays_computations = flags_non_arrays_computations

        self.list_of_crisprs_bona_fide = [
            row[1]
            for key in sorted(self.categories[0].keys())
            for row in self.categories[0][key]
        ]
        self.list_of_crisprs_alternative = [el[1] for key in sorted(self.categories[1].keys())
                                            for el in self.categories[1][key]]
        self.list_of_crisprs_possible = [el[1] for key in sorted(self.categories[2].keys())
                                         for el in self.categories[2][key]]

        self.hmm_model_is_elements = "tools/hmm_search/models_is_element.hmm"

        self.is_element_result = {}
        self.cas_result = {"status": "not_requested", "systems": []}
        self.strand_results = {}
        self.data_with_all_computations = {}

        self._get_complete_dna()
        self._calculate_all_non_array_values()

    def _get_complete_dna(self):
        self.source_header, self.dna = read_single_fasta(self.file_path)


    def _calculate_all_non_array_values(self):
        self._calculate_strand()
        self._calculate_leader()

        if self.flags_non_arrays_computations["flag_is"]:
            self._calculate_is_elements()

        if self.flags_non_arrays_computations.get("flag_cas", False):
            self._calculate_cas()

        self.data_with_all_computations = {"IS": self.is_element_result,
                                           "Cas": self.cas_result,
                                           "Strand": self.strand_results,
                                           "Leader": [self.leader_results_bona_fide, self.leader_results_alternative, self.leader_results_possible],
                                           "Downstream": [self.downstream_results_bona_fide, self.downstream_results_alternative, self.downstream_results_possible]}

    def _calculate_cas(self):
        predictor = CasAndraPredictor(
            input_fasta=self.file_path,
            source_id=self.source_header.split()[0],
            source_dna=self.dna,
            cpu_count=int(self.flags_non_arrays_computations.get("flag_cpu", 1)),
            gene_mode=self.flags_non_arrays_computations.get(
                "cas_gene_mode", "auto"
            ),
            timeout_seconds=int(self.flags_non_arrays_computations.get(
                "cas_timeout_seconds", 3600
            )),
        )
        self.cas_result = predictor.run()

    def _calculate_is_elements(self):
        fies = FullISElementSearch(full_dna=self.dna, list_of_crisprs=self.list_of_crisprs_bona_fide,
                                   hmm_model=self.hmm_model_is_elements, min_similarity=0.9, min_coverage=0.9)

        self.is_element_result = fies.output()

    def _calculate_strand(self):
        if self.flags_non_arrays_computations["flag_strand"]:
            categorized_crisprs = (
                ("Bona-fide", self.list_of_crisprs_bona_fide),
                ("Alternative", self.list_of_crisprs_alternative),
                ("Possible", self.list_of_crisprs_possible),
            )
            all_crisprs = [
                crispr
                for _, category_crisprs in categorized_crisprs
                for crispr in category_crisprs
            ]

            # Predict all categories in one batch, then restore the category-local
            # zero-based indexes expected by the reporting and leader components.
            strand_computation = StrandComputation(
                list_of_crisprs=all_crisprs,
            )
            all_strands = strand_computation.output()

            offset = 0
            for category_name, category_crisprs in categorized_crisprs:
                self.strand_results[category_name] = {
                    category_index: all_strands[offset + category_index]
                    for category_index in range(len(category_crisprs))
                }
                offset += len(category_crisprs)

        else:
            self.strand_results["Bona-fide"] = {index: "Forward (Orientation was not computed)"
                                                for index in range(len(self.list_of_crisprs_bona_fide))}
            self.strand_results["Alternative"] = {index: "Forward (Orientation was not computed)"
                                                  for index in range(len(self.list_of_crisprs_alternative))}
            self.strand_results["Possible"] = {index: "Forward (Orientation was not computed)"
                                               for index in range(len(self.list_of_crisprs_possible))}

    def _calculate_leader(self):
        flss_bona_fide = FullLeaderSeqSearch(self.list_of_crisprs_bona_fide, self.strand_results["Bona-fide"], self.dna)
        self.leader_results_bona_fide, self.downstream_results_bona_fide = flss_bona_fide.output()

        flss_alternative = FullLeaderSeqSearch(self.list_of_crisprs_alternative, self.strand_results["Alternative"],
                                               self.dna)
        self.leader_results_alternative, self.downstream_results_alternative = flss_alternative.output()

        flss_possible = FullLeaderSeqSearch(self.list_of_crisprs_possible, self.strand_results["Possible"], self.dna)
        self.leader_results_possible, self.downstream_results_possible = flss_possible.output()

    def output(self):
        return self.data_with_all_computations
