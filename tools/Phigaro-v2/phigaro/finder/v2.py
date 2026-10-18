from __future__ import division
import numpy as np

from .base import AbstractFinder, Phage


class V2Finder(AbstractFinder):
    def __init__(self, window_len, threshold_min, threshold_max, mode):
        """
        :type window_len: int
        :type threshold_min: float
        :type threshold_max: float
        """
        self.window_len = window_len
        self.threshold_min = threshold_min
        self.threshold_max = threshold_max
        self.mode = mode

    def find_phages(self, bacteria_npn, bacteria_gc):
        scores = calc_scores(bacteria_npn, self.window_len)
        if self.mode != "without_gc":
            gc_scores = calc_scores(bacteria_gc, self.window_len)
            ranges = scan_phages(
                np.array(scores) * np.array(gc_scores),
                self.threshold_min,
                self.threshold_max,
            )
        else:
            ranges = scan_phages(
                np.array(scores), self.threshold_min, self.threshold_max
            )
        for begin, end in ranges:
            yield Phage(
                begin=begin,
                end=end,
                is_prophage=True,
            )


def calc_scores(npn, window_len):
    """Zero-padded triangular sums in the original per-window addition order."""
    if len(npn) == 0:
        return []
    half = window_len // 2
    padded = np.concatenate([np.zeros(half), npn, np.zeros(half)])
    kernel = tri_kernel(np.arange(half * 2 + 1) / (half * 2))
    windows = np.lib.stride_tricks.sliding_window_view(padded, len(kernel))
    # Summing rows preserves the native loop's floating-point reduction order.
    return np.sum(windows * kernel, axis=1).tolist()


def tri_kernel(pos):
    return 1 - abs(pos - 0.5) / 0.5


def scan_phages(scores, threshold_min, threshold_max):
    ranges = []

    begin = None
    end = None
    has_peak = False

    # A list sentinel closes a terminal interval without mutating caller data.
    for i, v in enumerate([*scores, 0.0]):
        if v >= threshold_max:
            has_peak = True
        if v >= threshold_min:
            if begin is None:
                begin = i
            end = i
        else:
            if begin is not None and has_peak:
                ranges.append((begin, end))
            begin = None
            end = None
            has_peak = False

    return ranges
