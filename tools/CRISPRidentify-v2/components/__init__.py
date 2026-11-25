"""CRISPRidentify analysis packages; start with docs/CODE_MAP.md.

Import concrete modules from their functional package. Historical flat module
names remain lazy aliases for external callers and serialized candidates.
"""

from ._compat import install_legacy_imports as _install_legacy_imports

_install_legacy_imports()
del _install_legacy_imports
