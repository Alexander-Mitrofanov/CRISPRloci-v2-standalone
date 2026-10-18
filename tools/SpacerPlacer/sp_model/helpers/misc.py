
import numpy as np
from timeit import default_timer as timer
import logging
import os
import sys

from Bio import Phylo

from types import ModuleType, FunctionType
from gc import get_referents

BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: ' + str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def create_logger(name, level, outfile=None):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False
    for existing_handler in list(logger.handlers):
        logger.removeHandler(existing_handler)
        existing_handler.close()

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    if outfile:
        outfile = os.fspath(outfile)
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
        fh = logging.FileHandler(outfile)

        fh.setLevel(level)
        fh.setFormatter(formatter)

        logger.addHandler(fh)
    return logger


class RunTimer:
    def __init__(self):
        self.start_time = timer()
        self.checkpoints = [self.start_time]

    def time_from_start(self):
        self.checkpoints.append(timer())
        return self.checkpoints[-1] - self.start_time

    def time_between_checkpoints(self, i, j):
        interval = self.checkpoints[j] - self.checkpoints[i]
        return interval

    def time_from_last_checkpoint(self):
        self.checkpoints.append(timer())
        interval = self.checkpoints[-1] - self.checkpoints[-2]
        return interval

    def __str__(self):
        time_from_start = [self.time_between_checkpoints(0, i) for i in range(len(self.checkpoints))]
        return f'Time since start of Timer {self.checkpoints[-1] - self.start_time};  ' \
               f'Checkpoint times (wrt. start): {time_from_start}'


def remove_completely_same_arrays(ls_arrays, ls_names):
    """Remove duplicate arrays in linear time while preserving input order."""
    if len(ls_arrays) != len(ls_names):
        raise ValueError("ls_arrays and ls_names must have the same length")
    cleaned_arrays = []
    cleaned_names = []
    combined_names = {}
    representative_by_array = {}
    for array, name in zip(ls_arrays, ls_names):
        key = tuple(array)
        representative = representative_by_array.get(key)
        if representative is None:
            representative_by_array[key] = name
            cleaned_arrays.append(array)
            cleaned_names.append(name)
        else:
            combined_names.setdefault(representative, []).append(name)
    return cleaned_arrays, cleaned_names, combined_names


class GraphContainsCycleException(Exception):
    pass


def topological_sort(s, adj_matrix, label_dict):
    """Consume a binary order graph, retaining the historical LIFO tie order.

    Cache incoming-edge counts once instead of summing a full matrix column
    after every removed edge. The supplied roots and matrix are still mutated.
    """
    incoming = adj_matrix.sum(axis=0)
    sorted_spacers = []
    while s:
        node = s.pop()
        sorted_spacers.append(node)
        for child in np.flatnonzero(adj_matrix[node] == 1).tolist():
            adj_matrix[node, child] = 0
            incoming[child] -= 1
            if incoming[child] == 0:
                s.append(child)
    if np.any(adj_matrix):
        raise GraphContainsCycleException('Spacer arrays contain (at least one) cycle!')
    return [label_dict[node] for node in sorted_spacers]


def estimate_gr_based_on_estimated_lr(ls_avg_array_len, ls_loss_rate, deletion_model, ls_alpha=None):
    ls_gain_rate = []
    for (i, length), lr in zip(enumerate(ls_avg_array_len), ls_loss_rate):
        if deletion_model == 'independent':
            ls_gain_rate.append(lr * length)
        elif deletion_model == 'block':
            ls_gain_rate.append(lr * length * ls_alpha[i])
        else:
            raise NotImplementedError('Estimating gr is not implemented for this model!')
    return ls_gain_rate


def write_tree_to_file(tree, path, fmt='newick'):
    with open(path, 'w') as file:
        Phylo.write(tree, file, format=fmt)
    return


def keep_len_constant(ls_gain_rate, ls_avg_array_length, ls_deletion_model, ls_alpha=None):
    ls_loss_rate = []
    for (i, gr), avg_len, del_m in zip(enumerate(ls_gain_rate), ls_avg_array_length, ls_deletion_model):
        if del_m == 'independent':
            ls_loss_rate.append(gr / avg_len)
        elif del_m == 'block':
            ls_loss_rate.append(gr / (avg_len * ls_alpha[i]))
        else:
            raise NotImplementedError('Keeping lr constant is not implemented for this model!')
    return ls_loss_rate


def uniqify_ls_order_preserving(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result
