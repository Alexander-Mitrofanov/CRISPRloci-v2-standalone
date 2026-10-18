import os
import pickle
import json
import csv
from os import listdir
from os.path import basename
from os.path import isfile, join
from urllib.parse import quote
from components.components_detection_refinement import CrisprCandidate
from components.coordinate_utils import crispr_candidate_output_coordinates
from components.coordinate_utils import is_forward_strand
from components.folder_scheduler import pipeline_result_key
from components.likely_tandem_transaction import validated_row_tags
from components.sequence_utils import reverse_complement
from components.reconstruction_boundary_family import boundary_family_summary
from components.possible_reconstruction import possible_reconstruction_summary
from components.contextual_spacer_rescue import contextual_spacer_summary


SUMMARY_FIELDS = (
    "ID", "Region index", "Start", "End", "Length", "Consensus repeat",
    "Repeat Length", "Average Spacer Length", "Number of spacers", "Strand",
    "Category", "Confidence score", "Geometry status", "Boundary family ID",
    "Display geometry ID", "Admissible geometry count", "Supported repeat count",
    "Original core repeat count", "Extent start", "Extent end",
)


def _boundary_family_note(record):
    family = boundary_family_summary(record)
    if family is None:
        return ""
    return (
        "Boundary uncertainty: {supported_repeat_count} supported repeat positions; "
        "{member_count} admissible full geometries. Coordinates, sequences and score "
        "below describe actual member {display_geometry_id}, not a uniquely resolved "
        "segmentation. Family: {family_id}. Complete coupled alternatives: report.json.\n\n"
    ).format(**family)


def _possible_reconstruction_note(record):
    summary = possible_reconstruction_summary(record)
    if summary is None:
        return ""
    return (
        "Checked extent of a Possible core: {core_repeat_count} original repeats "
        "to {repeat_count} checked repeats. The result remains Possible. "
        "Original core score: {core_original_own_score:.6f}; fresh core score: "
        "{core_fresh_own_score:.6f}; full path own score: {own_score:.6f}. "
        "Scores are not calibrated probabilities. Complete evidence: report.json.\n\n"
    ).format(**summary)


def _boundary_family_fasta_suffix(record):
    family = boundary_family_summary(record)
    if family is None:
        return ""
    return (
        " boundary_uncertain=true view=actual_admissible_member"
        " family={family_id} member={display_geometry_id}"
        " admissible_geometries={member_count} evidence=report.json"
    ).format(**family)


def _boundary_family_csv(record):
    contextual = contextual_spacer_summary(record)
    if contextual is not None:
        return ["context_checked_boundary_family" if contextual["member_count"] > 1
                else "context_checked_geometry", "", contextual["display_geometry_id"],
                str(contextual["member_count"]), str(contextual["repeat_count"]), "",
                str(contextual["union_span"][0] + 1), str(contextual["union_span"][1])]
    family = boundary_family_summary(record)
    if family is None:
        return ["reported_geometry", "", "", "", "", "", "", ""]
    return [family["status"], family["family_id"], family["display_geometry_id"],
            str(family["member_count"]), str(family["supported_repeat_count"]),
            str(family["core_repeat_count"]), str(family["union_span"][0] + 1),
            str(family["union_span"][1])]


def _contextual_spacer_note(record):
    summary = contextual_spacer_summary(record)
    if summary is None:
        return ""
    return (
        "Localized identical-spacer block: the complete {repeat_count}-repeat "
        "path and both maximal flanks were checked. The result remains Possible. "
        "Flank checks are correlated diagnostic context, not independent votes. "
        "Displayed path own score: {display_own_score:.6f}; "
        "{member_count} complete boundary interpretations retained. Coordinates "
        "and sequences describe actual member {display_geometry_id}. "
        "Scores are not calibrated probabilities. Complete evidence: report.json.\n\n"
    ).format(**summary)


def rev_compliment_seq(seq):
    """Compatibility wrapper for the historically misspelled public helper."""

    return reverse_complement(seq)


def iter_category_records(category):
    for interval in sorted(category):
        for record in category[interval]:
            yield interval, record


REGION_INDEX_CATEGORY_PRIORITY = (0, 2, 3, 4, 1)


def build_crispr_region_indexes(categories):
    """Return stable region IDs shared by all output representations."""

    region_indexes = {}
    next_index = 1
    for category_index in REGION_INDEX_CATEGORY_PRIORITY:
        for interval in sorted(categories[category_index]):
            if interval not in region_indexes:
                region_indexes[interval] = next_index
                next_index += 1
    return region_indexes



def _resolve_subfolders(folder_result, sub_folders=None):
    """Return deterministic, validated per-input result folders.

    An explicit list scopes aggregation to the current invocation.  This keeps
    successful old runs in a reused result root from leaking into a new global
    report.
    """

    if sub_folders is None:
        return tuple(sorted(
            folder for folder in listdir(folder_result)
            if not isfile(join(folder_result, folder))
        ))

    resolved = []
    seen = set()
    for folder in sub_folders:
        folder = os.fspath(folder)
        if (
            not folder
            or folder in (".", "..")
            or basename(os.path.normpath(folder)) != folder
        ):
            raise ValueError("invalid result subfolder name: {!r}".format(folder))
        if folder in seen:
            raise ValueError("duplicate result subfolder: {!r}".format(folder))
        path = join(folder_result, folder)
        if not os.path.isdir(path):
            raise ValueError("missing result subfolder: {}".format(path))
        seen.add(folder)
        resolved.append(folder)
    return tuple(sorted(resolved))


def _iter_fasta_records(path):
    """Yield strict FASTA records without line-pair truncation."""

    header = None
    sequence_lines = []
    with open(path, encoding="utf-8") as input_file:
        for line_number, raw_line in enumerate(input_file, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not sequence_lines:
                        raise ValueError(
                            "empty FASTA record {!r} in {}".format(header, path)
                        )
                    yield header, "".join(sequence_lines)
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        "empty FASTA header in {} line {}".format(
                            path, line_number
                        )
                    )
                sequence_lines = []
            else:
                if header is None:
                    raise ValueError(
                        "FASTA sequence before header in {} line {}".format(
                            path, line_number
                        )
                    )
                if ">" in line:
                    raise ValueError(
                        "embedded FASTA header marker in {} line {}".format(
                            path, line_number
                        )
                    )
                sequence_lines.append(line)

    if header is None:
        return
    if not sequence_lines:
        raise ValueError("empty FASTA record {!r} in {}".format(header, path))
    yield header, "".join(sequence_lines)


def _aggregate_csv_reports(
    folder_result,
    sub_folders,
    source_name,
    output_name,
    default_fields,
    require_sources=False,
):
    """Aggregate generated CSVs without string concatenation or truncation."""

    output_path = join(folder_result, output_name)
    expected_fields = None
    global_counter = 1
    with open(output_path, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.writer(output_file)
        for sub_folder in sub_folders:
            source_path = join(folder_result, sub_folder, source_name)
            if not isfile(source_path):
                if require_sources:
                    raise ValueError(
                        "missing current result report: {}".format(source_path)
                    )
                continue

            with open(
                source_path, newline="", encoding="utf-8"
            ) as source_file:
                reader = csv.reader(source_file)
                try:
                    fields = next(reader)
                except StopIteration:
                    raise ValueError("empty CSV report: {}".format(source_path))
                if not fields or any(not field for field in fields):
                    raise ValueError(
                        "malformed CSV header in {}".format(source_path)
                    )
                if expected_fields is None:
                    expected_fields = fields
                    writer.writerow(["Name", "Global ID"] + fields)
                elif fields != expected_fields:
                    raise ValueError(
                        "inconsistent CSV header in {}".format(source_path)
                    )

                for line_number, row in enumerate(reader, 2):
                    if not row or not any(row):
                        continue
                    if len(row) != len(fields):
                        raise ValueError(
                            "malformed CSV row in {} line {}".format(
                                source_path, line_number
                            )
                        )
                    writer.writerow(
                        [sub_folder, global_counter] + row
                    )
                    global_counter += 1

        if expected_fields is None:
            writer.writerow(["Name", "Global ID"] + list(default_fields))


class RevComComputation:
    def __init__(self, crispr_candidate):
        self.forward_crispr_candidate = crispr_candidate
        self.rev_com_candidate = None
        self._compute_rev_com_candidate()

    def _compute_rev_com_candidate(self):
        list_repeats_f = self.forward_crispr_candidate.list_repeats
        list_repeats_gaped_f = self.forward_crispr_candidate.list_repeats_gaped
        list_repeat_starts = self.forward_crispr_candidate.list_repeat_starts
        list_spacers = self.forward_crispr_candidate.list_spacers

        list_repeats_rev_com = [rev_compliment_seq(repeat) for repeat in list_repeats_f][::-1]

        list_repeats_gaped_rev_com = [rev_compliment_seq(repeat_gaped)
                                      for repeat_gaped in list_repeats_gaped_f][::-1]

        list_repeat_starts_rev_com = list_repeat_starts[::-1]

        list_spacers_rev_com = [rev_compliment_seq(spacer) for spacer in list_spacers][::-1]

        self.rev_com_candidate = CrisprCandidate(list_repeats=list_repeats_rev_com,
                                                 list_repeats_gaped=list_repeats_gaped_rev_com,
                                                 list_repeat_starts=list_repeat_starts_rev_com,
                                                 list_spacers=list_spacers_rev_com)

    def output(self):
        return self.rev_com_candidate


class OutputCrispr:
    def __init__(self, crispr_candidate, coordinate_candidate=None):
        coordinate_candidate = coordinate_candidate or crispr_candidate
        coordinates = crispr_candidate_output_coordinates(coordinate_candidate)
        self.start = coordinates["start"]
        self.end = coordinates["end"]
        self.number_of_repeats = len(crispr_candidate.list_repeats)
        self.avg_repeat_length = len(crispr_candidate.consensus)
        if crispr_candidate.list_spacers:
            self.avg_spacer_length = int(
                sum(len(spacer) for spacer in crispr_candidate.list_spacers)
                / len(crispr_candidate.list_spacers)
            )
        else:
            self.avg_spacer_length = 0
        self.consensus = crispr_candidate.consensus
        self.list_repeats = crispr_candidate.list_repeats
        self.list_spacers = crispr_candidate.list_spacers

        self.list_repeat_starts = [x + 1 for x in crispr_candidate.list_repeat_starts]
        represented_length = sum(len(repeat) for repeat in self.list_repeats)
        represented_length += sum(len(spacer) for spacer in self.list_spacers)
        if represented_length != coordinates["length"]:
            raise ValueError(
                "oriented CRISPR components do not reconstruct the source span"
            )

        source_starts = list(crispr_candidate.list_repeat_starts)
        self.list_spacer_starts = []
        if self.list_spacers:
            if all(
                    source_starts[index] < source_starts[index + 1]
                    for index in range(len(source_starts) - 1)):
                for index, spacer in enumerate(self.list_spacers):
                    spacer_start = source_starts[index] + len(
                        self.list_repeats[index]
                    )
                    if spacer_start + len(spacer) != source_starts[index + 1]:
                        raise ValueError(
                            "forward spacer indexes do not reconstruct source geometry"
                        )
                    self.list_spacer_starts.append(spacer_start + 1)
            elif all(
                    source_starts[index] > source_starts[index + 1]
                    for index in range(len(source_starts) - 1)):
                for index, spacer in enumerate(self.list_spacers):
                    spacer_start = source_starts[index + 1] + len(
                        self.list_repeats[index + 1]
                    )
                    if spacer_start + len(spacer) != source_starts[index]:
                        raise ValueError(
                            "reverse spacer indexes do not reconstruct source geometry"
                        )
                    self.list_spacer_starts.append(spacer_start + 1)
            else:
                raise ValueError(
                    "oriented repeat starts must be strictly monotonic"
                )
        self.dot_representation = crispr_candidate.dot_repr()
        self.dot_representation_web_server = crispr_candidate.dot_repr_web_server()


class SimpleOutputMaker:
    def __init__(self, categories, non_array_data, result_path, list_features):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_simple_txt_files()

    def _index_crispr_candidates(self):
        self.dict_crispr_indexes = build_crispr_region_indexes(self.categories)

    def _write_simple_txt_files(self):
        if not os.path.exists(self.result_path):
            os.makedirs(self.result_path)

        for obsolete_name in (
                "Alternative_Candidates.txt",
                "Possible_Discarded_Candidates.txt"):
            obsolete_path = os.path.join(self.result_path, obsolete_name)
            if os.path.exists(obsolete_path):
                os.remove(obsolete_path)

        f_possible = open(self.result_path + '/Possible_Candidates.txt', 'w')
        f_low_score = open(self.result_path + '/Low_Score_Candidates.txt', 'w')

        file_names = [f_possible, f_low_score]
        category_specs = ((2, "Possible"), (4, "Low score"))

        for (category_index, category_name), file_name in zip(
                category_specs, file_names):

            records = list(iter_category_records(self.categories[category_index]))

            if records:
                for index, (interval, record) in enumerate(records):
                    score, array, feature_info = record[:3]
                    array_index = self.dict_crispr_indexes[interval]
                    if category_name in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"][category_name][index]
                    else:
                        strand = "Forward (Orientation was not computed)"
                    if not is_forward_strand(strand):
                        crispr = RevComComputation(array).output()
                    else:
                        crispr = array

                    coordinates = crispr_candidate_output_coordinates(array)
                    crispr_stats = array.compute_stats()
                    file_name.write(
                        "{} CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                            .format(category_name, array_index, coordinates["start"],
                                    coordinates["end"], crispr_stats["number_repeats"],
                                    crispr_stats["avg_repeat"], crispr_stats["avg_spacer"]))

                    tags = validated_row_tags(record)
                    if tags:
                        file_name.write("Tags: {}\n\n".format(", ".join(tags)))

                    file_name.write(_boundary_family_note(record))
                    file_name.write(_possible_reconstruction_note(record))
                    file_name.write(_contextual_spacer_note(record))
                    file_name.write(crispr.dot_repr())

                    file_name.write(f"\nStrand: {strand}\n\n")

                    file_name.write("\n")
                    list_reported_features = []
                    for feature_index, feature_list in enumerate(self.list_features):
                        for feature, value in zip(feature_list, feature_info[feature_index][0]):
                            if feature not in list_reported_features:
                                file_name.write("{}: {}\n".format(feature, value))
                                list_reported_features.append(feature)

                    file_name.write("\n")
                    file_name.write("Confidence score: {}\n\n\n".format(score))

                    file_name.write('\n{}\n\n'.format('=' * 100))

            file_name.close()


class SummaryOutputMaker:
    def __init__(self, result_path, categories, non_array_data, header, list_feature_names):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.header = header
        self.list_feature_names = list_feature_names

        self._make_text_summary()

    def _make_text_summary(self):
        result_path = self.result_path + '/Bona-Fide_Candidates.txt'
        records = [
            record for _, record in iter_category_records(self.categories[0])
        ]
        list_crisprs = [record[1] for record in records]
        list_scores = [record[0] for record in records]
        list_feature_vectors = [record[2] for record in records]

        with open(result_path, "w") as f:
            f.write(self.header)
            f.write("\n")

            for index, array in enumerate(list_crisprs):
                strand = self.non_array_data["Strand"]["Bona-fide"][index]
                if is_forward_strand(strand):
                    output_crispr = OutputCrispr(array)
                else:
                    output_crispr = OutputCrispr(
                        RevComComputation(array).output(),
                        coordinate_candidate=array,
                    )

                crispr_index = index + 1
                start, end = output_crispr.start, output_crispr.end
                number_of_repeats = output_crispr.number_of_repeats
                avg_length_repeat, avg_length_spacer = output_crispr.avg_repeat_length, output_crispr.avg_spacer_length

                if index > 0:
                    f.write('\n{}\n\n'.format('=' * 100))
                f.write(
                    "CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                        .format(str(crispr_index), start, end, number_of_repeats,
                                avg_length_repeat, avg_length_spacer))

                f.write(_boundary_family_note(records[index]))
                f.write(_possible_reconstruction_note(records[index]))
                f.write(_contextual_spacer_note(records[index]))
                f.write(output_crispr.dot_representation)

                f.write("\n")

                f.write("Leader region\n")

                leader = self.non_array_data["Leader"][0][index]

                f.write(leader)

                f.write("\n\n")

                f.write("Downstream region\n")

                downstream = self.non_array_data["Downstream"][0][index]

                f.write(downstream)

                f.write("\n\nStrand: {}\n\n".format(strand))

                f.write("#   Array features:\n")
                list_reported_features = []
                for group_of_features, feature_vector in zip(self.list_feature_names,
                                                             list_feature_vectors[index]):
                    for feature_name, feature_value in zip(group_of_features, feature_vector[0]):
                        if feature_name not in list_reported_features:
                            f.write("#   {}: {}\n".format(feature_name, feature_value))
                            list_reported_features.append(feature_name)

                f.write("_" * 30)
                f.write("\n")
                f.write("#   Confidence score: {}\n\n".format(list_scores[index]))

                if index in self.non_array_data["IS"]:
                    f.write("IS Element: {} [{}-{}]\n\n".format(self.non_array_data["IS"][index][4],
                                                                self.non_array_data["IS"][index][0],
                                                                self.non_array_data["IS"][index][1]))


class SummaryMakerCSV:
    def __init__(self, result_path, categories, non_array_data):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_csv_summary()

    def _index_crispr_candidates(self):
        self.dict_crispr_indexes = build_crispr_region_indexes(self.categories)


    def _write_csv_summary(self):

        all_arrays = [
            record[1] for category_index in (0, 2)
            for _, record in iter_category_records(self.categories[category_index])
        ]

        if all_arrays:

            result_csv_path = self.result_path + '/Summary.csv'
            with open(result_csv_path, "w") as f:
                f.write(",".join(SUMMARY_FIELDS))
                f.write("\n")
                global_index = 1
                for category_index, category in zip([0, 2], ["Bona-fide", "Possible"]):
                    records = list(iter_category_records(self.categories[category_index]))
                    for index, (interval, record) in enumerate(records):
                        score, array = record[:2]
                        array_index = self.dict_crispr_indexes[interval]
                        strand = self.non_array_data["Strand"][category][index]
                        crispr = array
                        coordinates = crispr_candidate_output_coordinates(crispr)
                        crispr_stats = crispr.compute_stats()
                        crispr_index = str(array_index)
                        start = str(coordinates["start"])
                        end = str(coordinates["end"])
                        length = str(coordinates["length"])
                        if is_forward_strand(strand):
                            consensus_repeat = crispr.consensus
                        else:
                            consensus_repeat = rev_compliment_seq(crispr.consensus)

                        repeat_length = str(crispr_stats["avg_repeat"])
                        average_spacer_length = str(crispr_stats["avg_spacer"])
                        number_of_spacers = str(crispr_stats["number_repeats"] - 1)

                        string_to_write = ",".join([str(global_index), crispr_index, start, end, length,
                                                    consensus_repeat, repeat_length, average_spacer_length,
                                                    number_of_spacers, strand, category, str(score)]
                                                   + _boundary_family_csv(record))
                        global_index += 1
                        f.write(string_to_write)
                        f.write("\n")
        else:
            # Explicitly replace any prior positive result.  A header-only
            # file is an unambiguous current run with zero accepted arrays.
            result_csv_path = self.result_path + "/Summary.csv"
            with open(result_csv_path, "w") as f:
                f.write(",".join(SUMMARY_FIELDS))
                f.write("\n")


class PickleOutputMaker:
    def __init__(self, file_path, pickle_result_folder, categories):
        self.file_name = file_path
        self.pickle_result_folder = pickle_result_folder
        self.categories = categories

        self._write_pickle()

    def _write_pickle(self):
        try:
            os.mkdir(self.pickle_result_folder)
        except OSError:
            pass

        file_base = basename(self.file_name)
        acc_num = pipeline_result_key(file_base)
        reported_categories = list(self.categories)
        reported_categories[1] = {}
        reported_categories[3] = {}
        with open(self.pickle_result_folder + '/' + acc_num + '.pkl', "wb") as output:
            pickle.dump(reported_categories, output)


class JsonOutputMaker:
    def __init__(self, file_path, json_result_folder, categories, non_array_data):
        self.file_path = file_path
        self.json_result_folder = json_result_folder
        self.categories = categories
        self.non_array_data = non_array_data

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_json()

    def _index_crispr_candidates(self):
        self.dict_crispr_indexes = build_crispr_region_indexes(self.categories)


    def _write_json(self):
        try:
            os.mkdir(self.json_result_folder)
        except OSError:
            pass
        json_dictionary = {}
        global_index = 1
        dictionary_all_arrays = {}

        category_specs = (
            (0, "Bona-fide"),
            (2, "Possible"),
            (4, "Low score"),
        )

        for category_index, category_name in category_specs:
            dict_category = {}

            records = list(iter_category_records(self.categories[category_index]))

            if records:
                for index, (interval, record) in enumerate(records):
                    score, array = record[:2]
                    array_index = self.dict_crispr_indexes[interval]
                    if category_name in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"][category_name][index]
                    else:
                        strand = "Forward (Orientation was not computed)"
                    is_forward_strand(strand)

                    coordinates = crispr_candidate_output_coordinates(array)

                    dict_crispr = self.crispr_candidate_to_dictionary(array)
                    dict_crispr["strand"] = strand
                    dict_crispr.update(coordinates)
                    dict_crispr["array_index"] = array_index
                    dict_crispr["certainty_score"] = score
                    tags = validated_row_tags(record)
                    if tags:
                        dict_crispr["tags"] = list(tags)
                    family = boundary_family_summary(record)
                    if family is not None:
                        dict_crispr["boundary_family"] = family
                        dict_crispr["boundary_family_evidence"] = "report.json:classification_basis.boundary_family"
                    possible_summary = possible_reconstruction_summary(record)
                    if possible_summary is not None:
                        dict_crispr["possible_reconstruction"] = possible_summary
                    contextual_summary = contextual_spacer_summary(record)
                    if contextual_summary is not None:
                        dict_crispr["contextual_spacer_context"] = contextual_summary

                    dict_category[global_index] = dict_crispr
                    global_index += 1

            dictionary_all_arrays[category_name] = dict_category

        file_base = basename(self.file_path)
        acc_num = pipeline_result_key(file_base)

        json_dictionary["Arrays"] = dictionary_all_arrays
        json_dictionary["coordinate_system"] = {
            "indexing": 1,
            "end_inclusive": True,
            "reference": "input_forward_strand",
        }

        json_representation = json.dumps(json_dictionary)
        with open(self.json_result_folder + '/' + acc_num + '.json', "w") as f:
            f.write(json_representation)

    @staticmethod
    def crispr_candidate_to_dictionary(crispr_candidate):
        list_repeat_starts = [
            start + 1
            for start in crispr_candidate.list_repeat_starts
        ]
        list_repeat_ends = [
            start + len(repeat)
            for start, repeat in zip(
                crispr_candidate.list_repeat_starts,
                crispr_candidate.list_repeats,
            )
        ]
        list_repeats = crispr_candidate.list_repeats
        list_spacers = crispr_candidate.list_spacers
        consensus_repeat = crispr_candidate.consensus
        dot_representation = crispr_candidate.dot_repr()
        dot_representation_web_server = crispr_candidate.dot_repr_web_server()

        dict_crispr = {"list_repeats": list_repeats,
                       "list_spacers": list_spacers,
                       "list_repeat_starts": list_repeat_starts,
                       "list_repeat_ends": list_repeat_ends,
                       "consensus_repeat": consensus_repeat,
                       "dot_representation": dot_representation,
                       "dot_representation_web_server": dot_representation_web_server}

        return dict_crispr




class CompleteFolderSummaryMaker:
    def __init__(self, folder_result, sub_folders=None):
        self.folder_result = folder_result
        self.requested_sub_folders = sub_folders

        self._list_sub_folders_()
        self._make_complete_summary()

    def _list_sub_folders_(self):
        self.sub_folders = _resolve_subfolders(
            self.folder_result, self.requested_sub_folders
        )

    def _make_complete_summary(self):
        _aggregate_csv_reports(
            self.folder_result,
            self.sub_folders,
            "Summary.csv",
            "Complete_summary.csv",
            SUMMARY_FIELDS,
            require_sources=self.requested_sub_folders is not None,
        )


class FastaOutputArrayMaker:
    def __init__(self, folder_result, categories, non_array_data):
        self.folder_result = folder_result
        self.categories = categories
        self.non_array_data = non_array_data

        self._make_fasta_files()

    def _make_fasta_files(self):
        result_path_repeats = join(self.folder_result, 'Repeats.fasta')
        result_path_spacers = join(self.folder_result, 'Spacers.fasta')
        result_path_arrays = join(self.folder_result, 'Arrays.fasta')
        records = [record for _, record in iter_category_records(self.categories[0])]
        list_crisprs = [record[1] for record in records]

        # Always replace prior optional FASTA artifacts in a reused per-genome
        # result directory, including when the current run has no bona-fide
        # arrays.
        for output_path in (
            result_path_repeats, result_path_spacers, result_path_arrays
        ):
            with open(output_path, "w"):
                pass

        if list_crisprs:

            with open(result_path_repeats, "w") as fr:
                with open(result_path_spacers, "w") as fs:
                    with open(result_path_arrays, "w") as fa:
                        for index, array in enumerate(list_crisprs):
                            family_suffix = _boundary_family_fasta_suffix(records[index])

                            strand = self.non_array_data["Strand"]["Bona-fide"][index]
                            if is_forward_strand(strand):
                                output_crispr = OutputCrispr(array)
                            else:
                                output_crispr = OutputCrispr(
                                    RevComComputation(array).output(),
                                    coordinate_candidate=array,
                                )

                            crispr_index = index + 1
                            start, end = output_crispr.start, output_crispr.end

                            repeats = output_crispr.list_repeats
                            spacers = output_crispr.list_spacers
                            array_seq = "".join([r+s for r, s in zip(repeats, spacers)]) + repeats[-1]

                            for repeat_number, repeat in enumerate(repeats, 1):
                                header_line = f">CRISPR_{crispr_index}_{start}_{end}_repeat_{repeat_number}{family_suffix}\n"
                                repeat_line = repeat + "\n"
                                fr.write(header_line)
                                fr.write(repeat_line)

                            for spacer_number, spacer in enumerate(spacers, 1):
                                header_line = f">CRISPR_{crispr_index}_{start}_{end}_spacer_{spacer_number}{family_suffix}\n"
                                repeat_line = spacer + "\n"
                                fs.write(header_line)
                                fs.write(repeat_line)

                            fa.write(f">CRISPR_{crispr_index}_{start}_{end}{family_suffix}\n")
                            fa.write(array_seq + "\n")


class CompleteFastaOutputMaker:
    def __init__(self, folder_result, sub_folders=None):
        self.folder_result = folder_result
        self.requested_sub_folders = sub_folders

        self._list_sub_folders_()
        self._make_complete_summary()

    def _list_sub_folders_(self):
        self.sub_folders = _resolve_subfolders(
            self.folder_result, self.requested_sub_folders
        )

    def _make_complete_summary(self):
        summary_path_repeats = join(self.folder_result, "Complete_repeat_dataset.fasta")
        summary_path_spacers = join(self.folder_result, "Complete_spacer_dataset.fasta")
        summary_path_arrays = join(self.folder_result, "Complete_array_dataset.fasta")

        with open(summary_path_repeats, "w") as fr:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Repeats.fasta")
                if os.path.exists(complete_path):
                    for header, sequence in _iter_fasta_records(complete_path):
                        accession = quote(sub_folder.strip(), safe="")
                        new_header = ">" + accession + "_" + header
                        fr.write(new_header + "\n")
                        fr.write(sequence + "\n")

        with open(summary_path_spacers, "w") as fs:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Spacers.fasta")
                if os.path.exists(complete_path):
                    for header, sequence in _iter_fasta_records(complete_path):
                        accession = quote(sub_folder.strip(), safe="")
                        new_header = ">" + accession + "_-_" + header
                        fs.write(new_header + "\n")
                        fs.write(sequence + "\n")

        with open(summary_path_arrays, "w") as fa:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Arrays.fasta")
                if os.path.exists(complete_path):
                    for header, sequence in _iter_fasta_records(complete_path):
                        accession = quote(sub_folder.strip(), safe="")
                        new_header = ">" + accession + "_" + header
                        fa.write(new_header + "\n")
                        fa.write(sequence + "\n")

class CompleteJsonOutputMaker:
    def __init__(self, folder_json_result, folder_text_tesult, result_keys=None):
        self.folder_json_result = folder_json_result
        self.folder_text_result = folder_text_tesult
        self.result_keys = result_keys

        self._make_complete_json_summary()

    def _make_complete_json_summary(self):
        if self.result_keys is None:
            json_files_in_folder = sorted(
                join(self.folder_json_result, file)
                for file in listdir(self.folder_json_result)
                if (
                    isfile(join(self.folder_json_result, file))
                    and file.endswith(".json")
                )
            )
        else:
            json_files_in_folder = []
            for result_key in sorted(self.result_keys):
                json_file = join(
                    self.folder_json_result, result_key + ".json"
                )
                if not isfile(json_file):
                    raise ValueError(
                        "missing JSON result for current input: {}".format(
                            json_file
                        )
                    )
                json_files_in_folder.append(json_file)

        complete_json = {}
        for json_file in json_files_in_folder:
            with open(json_file) as f:
                json_data = json.load(f)
                file_name = pipeline_result_key(json_file)
                complete_json[file_name] = json_data

        with open(join(self.folder_text_result, "Complete_json_summary.json"), "w") as f:
            json.dump(complete_json, f, indent=4)
