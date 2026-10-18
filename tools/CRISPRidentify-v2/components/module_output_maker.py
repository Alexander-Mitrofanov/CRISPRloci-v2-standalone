import json
import os
from os.path import join

from components.components_output_maker import SimpleOutputMaker
from components.components_output_maker import SummaryOutputMaker
from components.components_output_maker import SummaryMakerCSV
from components.components_output_maker import PickleOutputMaker
from components.components_output_maker import FastaOutputArrayMaker
from components.components_output_maker import JsonOutputMaker

from components.cas_output import normalize_cas_data
from components.cas_output import stage_cas_outputs
from components.html_report import _publish_report_model
from components.likely_tandem_transaction import (
    resolve_likely_tandem_audit_inputs,
)
from components.output_transaction import OutputTransaction
from components.output_transaction import remove_completion_marker
from components.report_extensions import attach_report_extensions
from components.report_model import ReportModelBuilder
from components.reconstruction_validator import read_single_fasta
from components.candidate_evidence import stage_candidate_evidence, EVIDENCE_FILENAME
from components.reconstruction_transaction import (
    EVIDENCE_FILENAME as RECONSTRUCTION_EVIDENCE_FILENAME,
)
from components.reconstruction_boundary_family import validate_boundary_family_audit
from components.possible_reconstruction import (
    EVIDENCE_FILENAME as POSSIBLE_RECONSTRUCTION_EVIDENCE_FILENAME,
    validate_possible_reconstruction_audit,
)
from components.contextual_spacer_rescue import (
    EVIDENCE_FILENAME as CONTEXTUAL_SPACER_EVIDENCE_FILENAME,
    validate_contextual_spacer_audit,
)
from components.competing_terminal_extents import competing_terminal_extent_summaries
from components.competing_terminal_extent_outputs import stage_competing_extent_outputs
from components.contextual_spacer_outputs import stage_contextual_spacer_outputs


_MACHINE_PROFILE_LEGACY_ARTIFACTS = (
    "Alternative_Candidates.txt",
    "Arrays.fasta",
    "Bona-Fide_Candidates.txt",
    "Cas_genes.csv",
    "Cas_systems.csv",
    "Low_Score_Candidates.txt",
    "Possible_Candidates.txt",
    "Possible_Discarded_Candidates.txt",
    "Repeats.fasta",
    "Spacers.fasta",
    "Summary.csv",
)


class OutputMaker:
    def __init__(self, file_path, parameters, flags, result_path, pickle_result_path,
                 json_result_path, categories, non_array_data, list_features, header,
                 sequence_topology_evidence=None,
                 circular_array_join_result=None,
                 likely_tandem_records=None,
                 likely_tandem_stats=None,
                 likely_tandem_source_categories=None,
                 candidate_evidence_categories=None,
                 reconstruction_evidence=None,
                 possible_reconstruction_evidence=None,
                 contextual_spacer_evidence=None):
        self.file_path = file_path
        self.parameters = parameters
        self.flags = flags
        self.result_path = result_path
        self.pickle_result_path = pickle_result_path
        self.json_result_path = json_result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features
        self.header = header
        self.candidate_evidence_categories = candidate_evidence_categories
        self.reconstruction_evidence = reconstruction_evidence
        self.possible_reconstruction_evidence = possible_reconstruction_evidence
        self.contextual_spacer_evidence = contextual_spacer_evidence
        self.sequence_topology_evidence = sequence_topology_evidence
        self.circular_array_join_result = circular_array_join_result
        self.likely_tandem_records = likely_tandem_records
        self.likely_tandem_stats = likely_tandem_stats
        self.likely_tandem_source_categories = (
            likely_tandem_source_categories
        )

        self._make_output()

    def _make_output(self):
        completion_path = join(self.result_path, "report.json")
        remove_completion_marker(completion_path)
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
        publication = OutputTransaction()
        try:
            self._make_output_transaction(publication)
        finally:
            publication.discard()

    def _make_output_transaction(self, publication):
        source_header, source_dna = read_single_fasta(self.file_path)
        source_id = source_header.split()[0]
        cas_data = normalize_cas_data(
            (self.non_array_data or {}).get(
                "Cas", {"status": "not_requested", "systems": []}
            ),
            source_id=source_id,
            source_length=len(source_dna),
        )
        self.non_array_data = dict(self.non_array_data or {})
        self.non_array_data["Cas"] = cas_data

        # Fail before direct legacy writers publish anything if category tags,
        # tandem provenance, source binding, or extension ledgers are invalid.
        core_model = ReportModelBuilder(
            file_path=self.file_path,
            categories=self.categories,
            non_array_data=self.non_array_data,
            list_feature_names=self.list_features,
            parameters=self.parameters,
            flags=self.flags,
            provenance_categories=self.likely_tandem_source_categories,
        ).build()
        validate_boundary_family_audit(core_model, self.reconstruction_evidence)
        report_model = attach_report_extensions(
            core_model,
            include_biological_loci=bool(
                self.flags.get("flag_biological_loci", False)
            ),
            sequence_topology_evidence=self.sequence_topology_evidence,
            circular_array_join_result=self.circular_array_join_result,
            likely_tandem_records=self.likely_tandem_records,
            likely_tandem_stats=self.likely_tandem_stats,
            likely_tandem_source_categories=(
                self.likely_tandem_source_categories
            ),
            source_sequence=source_dna,
            competing_terminal_extent_context={
                "accepted_categories": self.categories,
                "provenance_categories": self.likely_tandem_source_categories,
                "parameters": self.parameters,
                "feature_names": self.list_features,
            },
        )
        # The extension builder returns a detached validated snapshot. Reuse it
        # for publication after legacy writers instead of rebuilding everything.
        del core_model
        validate_possible_reconstruction_audit(
            report_model, self.possible_reconstruction_evidence,
            provenance_categories=self.likely_tandem_source_categories,
            likely_tandem_records=self.likely_tandem_records,
        )
        validate_contextual_spacer_audit(
            report_model, self.contextual_spacer_evidence,
            source_sequence=source_dna,
            provenance_categories=self.likely_tandem_source_categories,
            likely_tandem_records=self.likely_tandem_records,
        )
        competing_extents = competing_terminal_extent_summaries(report_model)

        reconstruction_path = join(self.result_path, RECONSTRUCTION_EVIDENCE_FILENAME)
        if self.reconstruction_evidence is not None:
            publication.stage_text(reconstruction_path, json.dumps(
                self.reconstruction_evidence, allow_nan=False, sort_keys=True,
                separators=(",", ":"),
            ) + "\n")
        else:
            publication.delete(reconstruction_path)

        possible_path = join(self.result_path, POSSIBLE_RECONSTRUCTION_EVIDENCE_FILENAME)
        if self.possible_reconstruction_evidence is not None:
            publication.stage_text(possible_path, json.dumps(
                self.possible_reconstruction_evidence, allow_nan=False,
                sort_keys=True, separators=(",", ":"),
            ) + "\n")
        else:
            publication.delete(possible_path)

        contextual_path = join(self.result_path, CONTEXTUAL_SPACER_EVIDENCE_FILENAME)
        if self.contextual_spacer_evidence is not None:
            publication.stage_text(contextual_path, json.dumps(
                self.contextual_spacer_evidence, allow_nan=False,
                sort_keys=True, separators=(",", ":"),
            ) + "\n")
        else:
            publication.delete(contextual_path)

        if self.candidate_evidence_categories is not None:
            stage_candidate_evidence(
                publication, join(self.result_path, EVIDENCE_FILENAME),
                source_header=source_header, source_sequence=source_dna,
                categories=self.candidate_evidence_categories,
                parameters=self.parameters, feature_names=self.list_features,
                canonical_report=report_model,
            )
        else:
            # Do not leave evidence from an earlier run beside a fresh report.
            publication.delete(join(self.result_path, EVIDENCE_FILENAME))

        machine_profile = self.flags.get("output_profile") == "machine"
        if machine_profile:
            # A reused result directory must not masquerade as machine-only
            # output.  Deletions and report.json are one transaction, with the
            # canonical report still published strictly last.
            for name in _MACHINE_PROFILE_LEGACY_ARTIFACTS:
                publication.delete(join(self.result_path, name))
        else:
            SimpleOutputMaker(categories=self.categories,
                              result_path=self.result_path,
                              non_array_data=self.non_array_data,
                              list_features=self.list_features)

            SummaryOutputMaker(result_path=self.result_path,
                               categories=self.categories,
                               non_array_data=self.non_array_data,
                               header=self.header,
                               list_feature_names=self.list_features)

            SummaryMakerCSV(result_path=self.result_path,
                            categories=self.categories,
                            non_array_data=self.non_array_data)

            stage_cas_outputs(
                publication,
                self.result_path,
                cas_data,
                source_id=source_id,
                source_length=len(source_dna),
            )

            if self.flags["flag_fasta_report"] is True:
                FastaOutputArrayMaker(folder_result=self.result_path,
                                      categories=self.categories,
                                      non_array_data=self.non_array_data)

            else:
                self._remove_generated_files(
                    "Repeats.fasta", "Spacers.fasta", "Arrays.fasta"
                )

        stage_competing_extent_outputs(
            publication, self.result_path, competing_extents,
            machine_profile=machine_profile,
        )
        stage_contextual_spacer_outputs(
            publication, self.result_path, report_model,
            enabled=not machine_profile and self.flags.get("flag_fasta_report") is True,
        )

        if self.pickle_result_path:
            PickleOutputMaker(
                file_path=self.file_path,
                pickle_result_folder=self.pickle_result_path,
                categories=self.categories,
            )

        if self.json_result_path:
            JsonOutputMaker(
                file_path=self.file_path,
                json_result_folder=self.json_result_path,
                categories=self.categories,
                non_array_data=self.non_array_data,
            )

        _publish_report_model(
            report_model,
            self.result_path,
            write_html=self.flags.get("flag_html_report", False),
            publication=publication,
        )

    def _remove_generated_files(self, *names):
        """Remove known optional artifacts disabled for the current run."""

        for name in names:
            path = join(self.result_path, name)
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
