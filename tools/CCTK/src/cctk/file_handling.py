from __future__ import annotations

from .io import safe_open, read_arrays, read_fasta
from .errors import CCTKError
from dataclasses import dataclass, field
import os
import sys
import json
import re
from collections import Counter, defaultdict

# cctkpkg imports
from . import sequence_operations


class AssemblyCRISPRs:
    """Class for reading and storing contents of minced and blast output.

    Args:
      CRISPR_types_dict (dict):
        dictionary of fasta sequence of CRISPR repeats to compare to
        minced-identified repeats
      minced_file (str):
      Path to the minced output file to parse into this class instance

    Attributes:
      accession (str):
        name of file (up to extension) that minced was run on
      array_count (int):
        number of distinct arrays identified
      array_ids (dict):
        ID numbers of each array {Array_number : ID number}
      array_locations (dict):
        Locations of each array {Array_number : [contig, start, stop]}
      array_orientations (dict):
        dict of array numbers and bool of whether they are in the
        reverse direction {Array_number : reverse?}
      array_sizes (dict):
        Number of spacers in each array {Array_number : Spacer_count}
      arrays (dict):
        Sequences of spacers in each array
        {Array_number : [Spacer_sequences]}
      arrays_encoded (dict):
        ID numbers of spacers in each array
        {Array_number : [Spacer_IDs]}
      has_CRISPR (bool):
        Whether any CRISPR arrays were identified
      repeat (dict):
        For each array, the seqeunce of the mose common repeat
        {Array_number : repeat_sequence}
      repeat_scores (dict):
        For each array, the hamming distance when compared to the most
        similar repeat in the provided repeats file
        {Array_number : hamming_dist}
      repeat_types (dict):
        For each array, the most similar repeat in the provided
        repeats file {Array_number : repeat_type}
      repeat_arrays (dict):
        For each array, the repeat array ID and sequences of repeats
        found. {Repeat_array_ID : repeat_sequences}
    """

    def __init__(self, crispr_types_dict=None, minced_file=None):
        self.accession = ""
        self.array_count = 0
        self.arrays = {}
        self.has_crispr = False
        if minced_file:
            self.read_minced(crispr_types_dict, minced_file)

    def read_minced(self, crispr_types_dict, minced_file):
        self.accession = minced_file.split("/")[-1].split("_minced")[0]
        if os.stat(minced_file).st_size == 0:
            self.has_crispr = False
            return
        else:
            self.has_crispr = True
        with safe_open(minced_file, "r") as f:
            for line in f.readlines():
                if "Sequence" in line and "bp)" in line:
                    contig = re.match(r"Sequence '(.*)' \([0-9]+ bp\)", line)[1]
                    continue
                if "CRISPR" in line and "Range:" in line:
                    self.array_count += 1
                    array_num = int(re.match(r"CRISPR ([0-9]+) ", line)[1])
                    array = FoundArray()
                    array.genome = self.accession
                    array.contig = contig
                    array.start, array.stop = [
                        i for i in re.findall(r"Range:\s+(\d+)\W+(\d+)\W+", line)[0]
                    ]
                    continue
                if "[" in line:
                    array.repeats.append(re.match(r"\d+\s+(\w+)\s+", line)[1])
                    array.spacers.append(re.match(r"\d+\s+\w+\s+(\w+)\s+", line)[1])
                    continue
                if "Repeats" in line:
                    repeat = Counter(array.repeats).most_common(1)[0][0]
                    (array.repeat_id, array.repeat_score, array.reverse) = (
                        sequence_operations.get_repeat_info(crispr_types_dict, repeat)
                    )

                    if array.reverse:
                        array.repeats = [
                            sequence_operations.rev_comp(repeat)
                            for repeat in reversed(array.repeats)
                        ]
                        array.spacers = [
                            sequence_operations.rev_comp(spacer)
                            for spacer in reversed(array.spacers)
                        ]
                        array.start, array.stop = array.stop, array.start

                    # If there are lots of mismatches between the
                    # repeat and its best match. It's probably wrong.
                    if array.repeat_score > 5:
                        array.repeat_id = "unknown_CRISPR_type({})".format(array.repeat_id)
                    self.arrays[array_num] = array


@dataclass
class FoundArray:
    """Typed observed array state shared by BLAST and MinCED adapters."""

    genome: str = ""
    contig: str = ""
    start: int = 0
    stop: int = 0
    repeats: list[str] = field(default_factory=list)
    repeat_array_id: str = ""
    repeat_id: str = ""
    repeat_score: int = 0
    reverse: bool = False
    spacers: list[str] = field(default_factory=list)
    spacer_ids: list[str] = field(default_factory=list)
    id: str = ""


class BlastResult:
    """
    A class to store column contents from a blast result in an easily
    retrieved form.

    Attributes:
        qseqid (str): query (e.g., unknown gene) sequence id
        sseqid (str): subject (e.g., reference genome) sequence id
        pident (float): percentage of identical matches
        length (int): alignment length (sequence overlap)
        mismatch (int): number of mismatches
        gapopen (int): number of gap openings
        qstart (int): start of alignment in query
        qend (int): end of alignment in query
        sstart (int): start of alignment in subject
        send (int): end of alignment in subject
        evalue (str): expect value
        bitscore (float): bit score
        qlen (int): length of query sequence
        slen (int): length of subject sequence
        sseq (str): sequence of subject
        strand (str): Whether the blast hit was on the top (plus) or
          bottom (minus) strand of the DNA
    """

    def __init__(self, blast_line):
        bits = blast_line.split("\t")
        if len(bits) not in {14, 15}:
            raise CCTKError("BLAST row must have 14 or 15 fields")
        self.qseqid = bits[0]
        self.sseqid = bits[1]
        self.pident = float(bits[2])
        self.length = int(bits[3])
        self.mismatch = int(bits[4])
        self.gapopen = bits[5]
        self.qstart = int(bits[6])
        self.qend = int(bits[7])
        self.sstart = int(bits[8])
        self.send = int(bits[9])
        self.evalue = bits[10]
        self.bitscore = bits[11]
        self.qlen = int(bits[12])
        self.slen = int(bits[13])
        self.sseq = bits[14] if len(bits) == 15 else ""
        self.strand = "plus" if self.sstart < self.send else "minus"


def read_array_file(file: str) -> dict[str, list[str]]:
    return read_arrays(file)


def fasta_to_dict(FASTA_file: str) -> dict[str, str]:
    return read_fasta(FASTA_file)


def read_colour_scheme(infile):
    with safe_open(infile, "r") as fin:
        spacer_colours = json.load(fin)
    return spacer_colours


def write_colour_scheme(outfile, colours):
    with safe_open(outfile, "w", encoding="utf-8") as fout:
        json.dump(colours, fout, ensure_ascii=False, indent=4)


def read_colour_file(colour_file):
    with safe_open(colour_file, "r") as fin:
        cf_list = [i.strip().replace('"', "").replace("'", "") for i in fin.readlines()]
    return cf_list


def write_array_file(array_dict, outfile):
    with safe_open(outfile, "w") as fout:
        fout.write(
            "\n".join(
                [
                    "{}\t{}".format(k, " ".join([str(s) for s in v.spacers]))
                    for k, v in array_dict.items()
                ]
            )
            + "\n"
        )


def write_spacer_fasta(spacer_dict, outfile):
    with safe_open(outfile, "w") as fout:
        for seq, sp_id in spacer_dict.items():
            fout.write(f">{sp_id}\n{seq}\n")


def write_array_id_seq_file(spacer_id_dict, array_id_dict, outdir):

    # seq file
    with safe_open(outdir + "Array_seqs.txt", "w") as fout:
        for sps, array_id in array_id_dict.items():
            fout.write(f"{array_id}\t{sps}\n")

    # ID file
    with safe_open(outdir + "Array_IDs.txt", "w") as fout:
        for sps, array_id in array_id_dict.items():
            fout.write(
                "{}\t{}\n".format(
                    array_id, " ".join([str(spacer_id_dict[sp]) for sp in sps.split()])
                )
            )


def write_cr_sum_tabs(all_assemblies, outfile, append=False):
    # Determine delimiters based on file extension
    if outfile[-3:] == "csv":
        maj_delim = ","
        min_delim = "\t"
        inc_count = True
        # additional nonsense to stop excel assuming dates
        prepend_stuff = '"=""'
        append_stuff = '"""'
    else:
        maj_delim = "\t"
        min_delim = "|"
        inc_count = False

    if append:
        # If appending, store existing contents to avoid redundant lines
        existing_lines = []
        if os.path.isfile(outfile):
            with safe_open(outfile) as fin:
                existing_lines = [line for line in fin]
        outcontents = "".join(existing_lines)  # still have newlines
    else:
        outcontents = (
            maj_delim.join(
                [
                    "Sequence_ID",
                    "Has_CRISPR",
                    "Array_count",
                    "Spacers",
                    "Spacer_IDs",
                    "Array_IDs",
                    "Array_locations",
                    "Repeat_sequences",
                    "Array_CRISPR_types",
                    "Array_repeats",
                    "Array_repeat_score",
                ]
            )
            + "\n"
        )

    for entry in all_assemblies:
        # Construct line components
        if inc_count:  # If this is for excel, add the array number to the info
            sequence_id = prepend_stuff + entry.accession + append_stuff
            has_crispr = str(entry.has_crispr)
            array_count = str(entry.array_count)
            spacers = min_delim.join(
                "{}: {}".format(array_num, " ".join(array.spacers))
                for array_num, array in entry.arrays.items()
            )
            spacer_ids = (
                prepend_stuff
                + min_delim.join(
                    "{}: {}".format(array_num, " ".join(array.spacer_ids))
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )
            array_ids = (
                prepend_stuff
                + min_delim.join(
                    "{}: {}".format(array_num, array.id)
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )
            array_locations = (
                prepend_stuff
                + min_delim.join(
                    "{}: {} {}-{}".format(array_num, array.contig, array.start, array.stop)
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )
            repeat_sequences = min_delim.join(
                "{}: {}".format(array_num, Counter(array.repeats).most_common(1)[0][0])
                for array_num, array in entry.arrays.items()
            )
            array_crispr_types = (
                prepend_stuff
                + min_delim.join(
                    "{}: {}".format(array_num, array.repeat_id)
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )
            array_repeats = (
                prepend_stuff
                + min_delim.join(
                    "{}: {}".format(array_num, array.repeat_array_id)
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )
            array_repeat_scores = (
                prepend_stuff
                + min_delim.join(
                    "{}: {}({})".format(array_num, array.repeat_id, array.repeat_score)
                    for array_num, array in entry.arrays.items()
                )
                + append_stuff
            )

        else:  # If txt, array number just complicates downstream work
            sequence_id = entry.accession
            has_crispr = str(entry.has_crispr)
            array_count = str(entry.array_count)
            spacers = min_delim.join(
                "{}".format(" ".join(array.spacers)) for array in entry.arrays.values()
            )
            spacer_ids = min_delim.join(
                "{}".format(" ".join(array.spacer_ids)) for array in entry.arrays.values()
            )
            array_ids = min_delim.join("{}".format(array.id) for array in entry.arrays.values())
            array_locations = min_delim.join(
                "{} {}-{}".format(array.contig, array.start, array.stop)
                for array in entry.arrays.values()
            )
            repeat_sequences = min_delim.join(
                "{}".format(Counter(array.repeats).most_common(1)[0][0])
                for array in entry.arrays.values()
            )
            array_crispr_types = min_delim.join(
                "{}".format(array.repeat_id) for array in entry.arrays.values()
            )
            array_repeats = min_delim.join(
                "{}".format(array.repeat_array_id) for array in entry.arrays.values()
            )
            array_repeat_scores = min_delim.join(
                "{}({})".format(array.repeat_id, array.repeat_score)
                for array in entry.arrays.values()
            )

        line = [
            sequence_id,
            has_crispr,
            array_count,
            spacers,
            spacer_ids,
            array_ids,
            array_locations,
            repeat_sequences,
            array_crispr_types,
            array_repeats,
            array_repeat_scores,
        ]

        # If appending (i.e. open_mode = 'a'), check if this id already
        # in outfile. Quicker to only check open_mode if not appending
        if append:
            if any([sequence_id in line for line in existing_lines]):
                continue

        outcontents += maj_delim.join(line) + "\n"

    with safe_open(outfile, "w") as fout:
        fout.write(outcontents)


def write_array_reps_file(all_assemblies, outdir):
    array_reps_dict = defaultdict(list)
    for assembly in all_assemblies:
        acc = assembly.accession
        for array in assembly.arrays.values():
            array_reps_dict[array.id].append(acc)

    outcontents = []
    for arr, accs in array_reps_dict.items():
        outcontents.append("{}\t{}".format(arr, " ".join(accs)))

    with safe_open(outdir + "Array_representatives.txt", "w") as fout:
        fout.write("\n".join(outcontents))


def write_repeat_id_seq_file(repeat_id_lookup, outdir):
    # repeat file
    with safe_open(outdir + "Array_repeats.txt", "w") as fout:
        for rep_array_id in sorted(
            repeat_id_lookup.keys(),
            key=lambda x: (int(x.split("_")[0]), x.split("_")[1]),
        ):
            reps = repeat_id_lookup[rep_array_id]
            fout.write(f"{rep_array_id}\t{reps}\n")


def write_CRISPR_files(
    all_assemblies,
    spacer_id_dict,
    array_id_dict,
    repeat_id_lookup,
    cluster_reps_dict,
    outdir,
    append=False,
    min_shared=1,
):
    """Write output files from CRISPR identification scripts"""
    # Spacers in fasta format
    write_spacer_fasta(spacer_id_dict, outdir + "CRISPR_spacers.fna")

    # Array IDs and spacers
    write_array_id_seq_file(spacer_id_dict, array_id_dict, outdir)

    # Repeat array IDs and repeat seqs
    write_repeat_id_seq_file(repeat_id_lookup, outdir)

    # Array representatives
    write_array_reps_file(all_assemblies, outdir)

    # Cluster representatives if snp_threshold used to group spacers
    if len(cluster_reps_dict) > 0:
        write_clus_reps(cluster_reps_dict, spacer_id_dict, outdir, append)

    # CSV summary table
    write_cr_sum_tabs(all_assemblies, outdir + "CRISPR_summary_table.csv", append)
    # TXT summary table
    write_cr_sum_tabs(all_assemblies, outdir + "CRISPR_summary_table.txt", append)

    # Network

    # Make list of FoundArray instances representing all unique arrays
    array_dict = {}
    for assembly in all_assemblies:
        for array in assembly.arrays.values():
            if array.id not in array_dict:
                array_dict[array.id] = array
    array_list = [a for a in array_dict.values()]

    # Write bed file of array locations
    write_array_loc_bed(all_assemblies, outdir)

    if len(array_list) > 2000:
        print(
            "A large number of arrays found. Building the array network will take a very long time and is being skipped."
            + " If you would like to build the array network later, use 'cctk network'."
        )
        return

    # Build network of array spacer sharing
    network = sequence_operations.build_network(array_list)

    # Write network file
    write_network_file(network, outdir + "Array_network.txt", min_shared)

    # Identify network clusters for cluster file
    network_filt = [(edge.a, edge.b) for edge in network if edge.nshared >= min_shared]

    clusters = sequence_operations.identify_network_clusters(network_filt)

    # need to unpack clustes from dict for file
    cluster_list = [[k for k in clus.keys()] for clus in clusters]

    write_array_cluster_file(cluster_list, outdir + "Array_clusters.txt")


def read_assembly_list_file(filename):
    assemblies = []
    with safe_open(filename, "r") as fin:
        for line in fin:
            assemblies.append(line.strip())

    return assemblies


def read_array_types_file(filename: str) -> dict[str, str]:
    """Read complete two-column annotations without silently replacing duplicates."""
    array_types: dict[str, str] = {}
    with safe_open(filename, "r") as fin:
        for line in fin:
            bits = line.split()
            if not bits:
                continue
            if len(bits) != 2 or bits[0] in array_types:
                raise CCTKError("array types require two columns and unique identifiers")
            array_types[bits[0]] = bits[1]

    return array_types


def write_network_file(network, filename, min_shared=1):
    if len(network) == 0:
        return
    # If no types provided, don't include those columns in the outfile
    if network[0].a_type != "":
        outcontents = [
            "\t".join(
                [
                    "Array_A",
                    "Array_B",
                    "Shared_spacers",
                    "Jaccard_similarity",
                    "Array_A_type",
                    "Array_B_type",
                ]
            )
        ]
        for edge in network:
            if edge.nshared < min_shared:
                continue

            outcontents.append(
                "\t".join(
                    [
                        str(i)
                        for i in [
                            edge.a,
                            edge.b,
                            edge.nshared,
                            edge.jaccard,
                            edge.a_type,
                            edge.b_type,
                        ]
                    ]
                )
            )
    else:
        outcontents = ["\t".join(["Array_A", "Array_B", "Shared_spacers", "Jaccard_similarity"])]
        for edge in network:
            if edge.nshared < min_shared:
                continue

            outcontents.append(
                "\t".join([str(i) for i in [edge.a, edge.b, edge.nshared, edge.jaccard]])
            )

    with safe_open(filename, "w") as fout:
        fout.write("\n".join(outcontents) + "\n")


def write_array_cluster_file(cluster_list, filename):
    with safe_open(filename, "w") as fout:
        for cluster in cluster_list:
            fout.write(" ".join([str(c) for c in cluster]) + "\n")


def write_array_loc_bed(all_assemblies, outdir):
    bed_contents = [
        "#" + "\t".join(["contig", "contigStart", "contigEnd", "name", "score", "strand"])
    ]
    for assembly in all_assemblies:
        for array in assembly.arrays.values():
            if array.reverse:
                bed_contents.append(
                    "\t".join(
                        [
                            array.contig,
                            str(int(array.stop) - 1),  # 0 base
                            str(array.start),
                            str(array.id),
                            "0",
                            "-",
                        ]
                    )
                )
            else:
                bed_contents.append(
                    "\t".join(
                        [
                            array.contig,
                            str(int(array.start) - 1),  # 0 base
                            str(array.stop),
                            str(array.id),
                            "0",
                            "+",
                        ]
                    )
                )

    with safe_open(outdir + "Array_locations.bed", "w") as fout:
        fout.write("\n".join(bed_contents))


def read_genome_reps_file(filename):
    """Read file with genome name and the array of interest for constrain"""
    genome_array_dict = defaultdict(list)
    outgroup_taxons = []
    with safe_open(filename, "r") as fin:
        for line in fin.readlines():
            bits = line.split()
            if bits[0].lower() == "outgroup":
                # replace underscores with spaces to relate to newick
                outgroup_taxons.append(bits[1].replace("_", " "))
                continue
            genome_array_dict[bits[1].replace("_", " ")].append(bits[0])

    return genome_array_dict, outgroup_taxons


def write_clus_reps(cluster_reps_dict, spacer_id_dict, outdir, append):

    # Fist check if we need to append
    if append:
        with safe_open(outdir + "Spacer_cluster_members.txt") as fin:
            # Store existing cluster information
            existing_dict = {}
            for line in fin:
                bits = line.split()
                existing_dict[bits[0]] = bits[1:]
        # Add new information to relevant clusters
        for reps_dict in cluster_reps_dict.values():
            # No need to include CR_type
            for rep, cluster in reps_dict.items():
                rep_id = spacer_id_dict[rep]
                if rep_id not in existing_dict:
                    # If this rep is new add it
                    existing_dict[rep_id] = cluster
                    continue

                # Else check if any new cluster members need appending
                for member in cluster:
                    if member not in existing_dict[rep_id]:
                        existing_dict[rep_id].append(member)
        # Prepare modified file
        outcontents = ""
        for rep, cluster in existing_dict.items():
            outcontents += "{}\t{}\n".format(rep, " ".join(cluster))

    else:
        # Otherwise, write it the easy way
        outcontents = ""
        for reps_dict in cluster_reps_dict.values():
            # No need to include CR_type
            for rep, cluster in reps_dict.items():
                outcontents += "{}\t{}\n".format(spacer_id_dict[rep], " ".join(cluster))

    with safe_open(outdir + "Spacer_cluster_members.txt", "w") as fout:
        fout.write(outcontents)


def check_append(outdir):
    """Check that the necessary files exist"""
    if outdir[-1] != "/":
        outdir += "/"

    if not os.path.isfile(outdir + "CRISPR_spacers.fna"):
        sys.stderr.write(
            "Append ERROR: CRISPR_spacers.fna not found in the "
            "provided outdir: '{}'. To append to a previous run, you must "
            "provide the CRISPR_spacers.fna from that run in the '{}'.\n"
            "".format(outdir, outdir)
        )
        raise CCTKError("analysis cannot continue; see diagnostic above")
