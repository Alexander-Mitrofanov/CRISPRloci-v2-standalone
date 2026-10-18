"""Native thread policy shared by CLI entry points and detector workers.

Configure this before importing numerical libraries. Detector processes own
the parallelism; their small classifier matrices do not need BLAS pools.
Sequential external tools may use the resolved per-genome CPU allocation.
"""

import os


_EXTERNAL_THREAD_BUDGET = 1


def native_thread_environment(cpu_budget):
    if (
            isinstance(cpu_budget, bool)
            or not isinstance(cpu_budget, int)
            or cpu_budget < 1):
        raise ValueError("native CPU budget must be a positive integer")
    return {
        "OMP_NUM_THREADS": str(cpu_budget),
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    }


def configure_native_threads(cpu_budget):
    """Set native limits before imports, or before work in a pool child."""
    global _EXTERNAL_THREAD_BUDGET
    os.environ.update(native_thread_environment(cpu_budget))
    _EXTERNAL_THREAD_BUDGET = cpu_budget


def external_thread_budget():
    return _EXTERNAL_THREAD_BUDGET
