import os
import pickle
import Bio.Phylo as Phylo
import io
import math


def unpickle(f):
    """
    Reads complete pkl file from given file handle.
    :param f:
    :return:
    """
    ls = []
    while True:
        try:
            ls.append(pickle.load(f))
        except EOFError:  # noqa: PERF203 - stream length is not encoded; EOF terminates concatenated pickles
            break
    return ls


def load_crispr_arrays_from_pickle(path):
    with open(path, "rb") as f:
        dict_crispr_by_repeats = pickle.load(f)
    return dict_crispr_by_repeats


def load_real_data_tree(path, ls_file_names):
    ls_trees = []
    for file_name in ls_file_names:
        ls_trees += Phylo.parse(os.path.join(path, file_name), format='newick')
    return ls_trees


def load_single_tree(path):
    with open(path, 'r', encoding='utf-8') as tree_file:
        return _validate_tree(Phylo.read(tree_file, format='newick'))


def load_single_tree_from_string(string):
    return _validate_tree(Phylo.read(io.StringIO(string), 'newick'))


def _validate_tree(tree):
    """Reject incomplete or nonphysical supplied branch lengths and leaf IDs."""
    names = [leaf.name for leaf in tree.get_terminals()]
    if any(not isinstance(name, str) or not name.strip() for name in names) or len(names) != len(set(names)):
        raise ValueError('Supplied tree leaves must have unique non-empty names')
    for node in tree.find_clades():
        if node is tree.root and node.branch_length is None:
            # There is no incoming root edge in a rooted phylogeny.
            node.branch_length = 0.0
        elif node.branch_length is None:
            raise ValueError('Supplied tree must specify every non-root branch length')
        if not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Supplied tree branch lengths must be finite and non-negative')
    return tree


def read_old_data(path, ls_files):
    ls_metadata = []
    ls_arrays = []
    for file in ls_files:
        print(os.path.join(path, file))
        with open(os.path.join(path, file)) as f:
            lines = (line.rstrip() for line in f.readlines())
            lines = [line for line in lines if line]
            head = lines[0].split(')')[0]
            head = head.strip('(')
            head = [h.strip().strip("'") for h in head.split(',')]
            for line in lines[1:]:
                line = line.split(':')
                metadata = line[0]
                arrays = line[1]
                arrays = arrays.split(',')
                arrays = [a.strip() for a in arrays]
                ls_metadata.append(metadata)
                ls_arrays.append(arrays)
    return {'head': head, 'metadata': ls_metadata, 'arrays': ls_arrays}
