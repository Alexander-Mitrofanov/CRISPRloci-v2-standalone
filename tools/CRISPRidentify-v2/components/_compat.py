"""Lazy aliases for the pre-package module names.

Keep legacy imports (including old pickles and external scripts) bound to the
same module object as canonical imports. Do not eagerly import the scientific
stack: CLI workers must be able to import components before CPU limits are set.
New code must use the canonical paths below.
"""

import sys
from importlib import import_module
from importlib.abc import Loader, MetaPathFinder
from importlib.util import spec_from_loader

LEGACY_MODULES = {
    "accepted_core_boundary_attestation": "rescue.accepted_core_boundary_attestation",
    "accepted_core_boundary_completion": "rescue.accepted_core_boundary_completion",
    "accepted_core_boundary_transaction": "rescue.accepted_core_boundary_transaction",
    "accepted_family_tandem_transaction": "rescue.accepted_family_tandem_transaction",
    "accepted_family_tandem_validator": "rescue.accepted_family_tandem_validator",
    "accepted_source_intact_conflict": "rescue.accepted_source_intact_conflict",
    "architecture_rescue_coding_context": "rescue.architecture_rescue_coding_context",
    "architecture_rescue_integration": "rescue.architecture_rescue_integration",
    "architecture_rescue_policy": "rescue.architecture_rescue_policy",
    "biological_loci": "evidence.biological_loci",
    "boundary_assessment": "evidence.boundary_assessment",
    "bounded_subprocess": "runtime.bounded_subprocess",
    "candidate_evidence": "reporting.candidate_evidence",
    "candidate_origins": "evidence.candidate_origins",
    "cas_html": "reporting.cas_html",
    "cas_output": "reporting.cas_output",
    "casandra_adapter": "annotation.casandra_adapter",
    "circular_array_join": "rescue.circular_array_join",
    "circular_array_join_integration": "rescue.circular_array_join_integration",
    "cli_utils": "cli.validation",
    "compact_feature_extraction": "scoring.compact_feature_extraction",
    "compact_ml_runtime": "scoring.compact_ml_runtime",
    "compact_model": "scoring.compact_model",
    "competing_terminal_extent_outputs": "reporting.competing_terminal_extent_outputs",
    "competing_terminal_extents": "reconstruction.competing_terminal_extents",
    "components_detection": "detection.search",
    "components_detection_refinement": "refinement.candidates",
    "components_evaluated_arrays_enhancement": "refinement.enhancement_search",
    "components_helpers": "sequence.fasta",
    "components_non_array_computations": "annotation.sequence_annotations",
    "components_output_maker": "reporting.legacy_exports",
    "contextual_spacer_outputs": "reporting.contextual_spacer_outputs",
    "contextual_spacer_rescue": "reconstruction.contextual_spacer_rescue",
    "coordinate_utils": "sequence.coordinates",
    "external_tools": "runtime.external_tools",
    "fast_accepted_family_rescue": "rescue.fast_accepted_family_rescue",
    "fast_run_rescue": "rescue.fast_run_rescue",
    "fast_two_repeat_terminal_repair": "rescue.fast_two_repeat_terminal_repair",
    "feature_schema": "scoring.feature_schema",
    "folder_scheduler": "runtime.folder_scheduler",
    "fuzzy_alignment": "detection.fuzzy_alignment",
    "fuzzy_alignment_merge": "detection.fuzzy_alignment_merge",
    "fuzzy_boundary_lattice": "detection.fuzzy_boundary_lattice",
    "fuzzy_lattice_storage": "detection.fuzzy_lattice_storage",
    "fuzzy_occurrence_paths": "detection.fuzzy_occurrence_paths",
    "fuzzy_search_utils": "detection.fuzzy_search_utils",
    "html_report": "reporting.html_report",
    "initial_competition": "evidence.initial_competition",
    "internal_completion_authority": "refinement.internal_completion_authority",
    "internal_local_core_bridge": "refinement.internal_local_core_bridge",
    "internal_repeat_completion": "refinement.internal_repeat_completion",
    "internal_repeat_transaction": "refinement.internal_repeat_transaction",
    "legacy_hash": "scoring.legacy_hash",
    "likely_tandem": "rescue.likely_tandem",
    "likely_tandem_transaction": "rescue.likely_tandem_transaction",
    "literal_fuzzy": "detection.literal_fuzzy",
    "localized_spacer_context": "reconstruction.localized_spacer_context",
    "module_detection": "detection.stage",
    "module_detection_refinement": "refinement.stage",
    "module_evaluated_arrays_enhancement": "refinement.enhancement",
    "module_evaluation": "scoring.evaluation",
    "module_non_array_computations": "annotation.stage",
    "module_output_maker": "reporting.stage",
    "occurrence_protection": "evidence.occurrence_protection",
    "output_transaction": "reporting.output_transaction",
    "portable_strand": "annotation.portable_strand",
    "possible_reconstruction": "reconstruction.possible_reconstruction",
    "reconstruction_boundary_family": "reconstruction.reconstruction_boundary_family",
    "reconstruction_reference": "reconstruction.reconstruction_reference",
    "reconstruction_selection": "reconstruction.reconstruction_selection",
    "reconstruction_transaction": "reconstruction.reconstruction_transaction",
    "reconstruction_validation": "reconstruction.reconstruction_validation",
    "reconstruction_validator": "sequence.validation",
    "repeat_boundary_repair": "refinement.repeat_boundary_repair",
    "report_extensions": "reporting.report_extensions",
    "report_model": "reporting.report_model",
    "report_provenance_validation": "reporting.report_provenance_validation",
    "runtime_presets": "cli.presets",
    "runtime_resources": "runtime.runtime_resources",
    "sequence_path_vectorizer": "scoring.sequence_path_vectorizer",
    "sequence_topology": "sequence.topology",
    "sequence_utils": "sequence.utils",
    "shared_repeat_augmentation": "rescue.shared_repeat_augmentation",
    "source_intact_existence_rescue": "rescue.source_intact_existence_rescue",
    "source_intact_recovery_transaction": "rescue.source_intact_recovery_transaction",
    "source_intact_terminal_authorization": "rescue.source_intact_terminal_authorization",
    "source_periodic_lattice_category_transaction": "rescue.source_periodic_lattice_category_transaction",
    "source_periodic_lattice_repair": "rescue.source_periodic_lattice_repair",
    "source_regular_subpath_recovery": "rescue.source_regular_subpath_recovery",
    "source_regular_subpath_recovery_transaction": "rescue.source_regular_subpath_recovery_transaction",
    "source_superlocus_augmentation": "rescue.source_superlocus_augmentation",
    "terminal_repair_completion": "refinement.terminal_repair_completion",
    "vmatch_occurrence_chain": "detection.vmatch_occurrence_chain",
}


class _LegacyLoader(Loader):
    def __init__(self, target):
        self.target = target

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        # Replace the temporary alias, keeping canonical import metadata intact.
        sys.modules[module.__name__] = import_module(self.target)


class _LegacyFinder(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        package, _, name = fullname.rpartition(".")
        if package != "components" or name not in LEGACY_MODULES:
            return None
        loader = _LegacyLoader("components." + LEGACY_MODULES[name])
        return spec_from_loader(fullname, loader)


def install_legacy_imports():
    if not any(isinstance(finder, _LegacyFinder) for finder in sys.meta_path):
        sys.meta_path.append(_LegacyFinder())
