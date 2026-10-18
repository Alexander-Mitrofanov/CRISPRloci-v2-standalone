import math
from statistics import NormalDist

import numpy as np


def test_significance_ratio_chi2(value_to_test, significance_value, n=1):
    if not 0 < significance_value < 1:
        raise ValueError('significance_value must be between zero and one')
    if n <= 0:
        raise ValueError('chi-square degrees of freedom must be greater than zero')
    if n == 1 and 1 - significance_value / 2 < 1:
        # A one-degree chi-square variable is the square of a standard normal.
        quantile = NormalDist().inv_cdf(1 - significance_value / 2) ** 2
    elif n == 2:
        # A two-degree chi-square variable is exponential with scale two.
        quantile = -2 * math.log(significance_value)
    else:
        # Preserve uncommon/fractional degrees of freedom without paying the
        # SciPy special-function import cost on every normal web job.
        from scipy.special import gammainccinv

        quantile = 2 * gammainccinv(n / 2, significance_value)
    return True if value_to_test > quantile else False, quantile


def test_ls_significance_ratio_chi2(ls_values_to_test, significance_value, n=1):
    ls_quantiles = []
    ls_test_results = []
    for val in ls_values_to_test:
        r, q = test_significance_ratio_chi2(val, significance_value, n=n)
        ls_quantiles.append(q)
        ls_test_results.append(r)
    return ls_test_results, ls_quantiles


def np_rate_exponential(rate, *args, **kwargs):
    if rate == 0:
        return np.random.exponential(np.inf, *args, **kwargs)
    return np.random.exponential(1 / rate, *args, **kwargs)


def p_stationary_frag_length(n, rate):
    prod = 1
    p = np.ones(n + 1)
    for i in range(n + 1):
        val = (i + 1) * (i + 2)
        prod *= val / (2 * rate) + 1
        p[i] = val / (2 * rate * prod)
    return p


def sample_stationary_frag_length(sample_size, rate, dist_limit=1000):
    """
    Samples the stationary distribution of the fragment model. This distribution is created empirically, the dist_limit
    determines the precision (in theory dist_limit = infty). (see Kupczok et al.)
    :param sample_size:
    :param rate:
    :param dist_limit:
    :return:
    """
    p = p_stationary_frag_length(dist_limit, rate)
    return np.random.choice(range(dist_limit + 1), sample_size, p=p / sum(p), replace=True)


def ev_stationary_frag_length(rate, dist_limit=1000):
    p = p_stationary_frag_length(dist_limit, rate)
    return sum([n * p[n] for n in range(len(p))])


def var_stationary_frag_length(rate, dist_limit=1000):
    p = p_stationary_frag_length(dist_limit, rate)
    ev = sum([n * p[n] for n in range(len(p))])
    ev_2 = sum([n ** 2 * p[n] for n in range(len(p))])
    return ev_2 - ev ** 2
