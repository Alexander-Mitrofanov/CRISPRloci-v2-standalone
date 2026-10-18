"""Shared, dependency-free exact terms for offline indexing and API verification.

Cas aliases are complete gene-name tokens explicitly present in source labels:
``C25_Cas5_cmr3_1`` supplies ``cas5`` and ``cmr3``. Variant tokens are preserved;
``Cas5f1`` does not supply ``cas5``. No family or subtype is inferred.
"""

from functools import lru_cache
import hashlib
import re


TERM_VERSION = "source-label-aliases-v1"
KINDS = {"arrays": 1, "repeats": 2, "spacers": 3, "cas-systems": 4, "cas-genes": 5}
_DNA = re.compile(r"[ACGTRYSWKMBDHVN]+", re.IGNORECASE)
_GENE_TOKEN = re.compile(r"(?<![a-z0-9])(?:cas|cmr|cs[a-z])[0-9]+[a-z0-9]*(?![a-z0-9])")


def normalize_term(kind: int, value: str) -> str:
    if kind == 5:
        return value.casefold()
    if kind in {1, 2, 3} and _DNA.fullmatch(value):
        return value.upper()
    return value


@lru_cache(maxsize=4096)
def _gene_terms(value: str) -> frozenset[str]:
    normalized = value.casefold()
    return frozenset({normalized, *_GENE_TOKEN.findall(normalized)})


def terms_for(kind: int, value: str | None) -> frozenset[str]:
    if not value:
        return frozenset()
    if kind == 5:
        # Source labels are normally short. Do not retain anomalously large
        # operator data in the shared cache.
        if len(value) > 1024:
            normalized = value.casefold()
            return frozenset({normalized, *_GENE_TOKEN.findall(normalized)})
        return _gene_terms(value)
    return frozenset({normalize_term(kind, value)})


def query_hash(kind: int, query: str) -> bytes:
    return hashlib.sha256(normalize_term(kind, query).encode("utf-8")).digest()[:16]


def gene_matches(annotation: str | None, profile: str | None, query: str) -> bool:
    normalized = normalize_term(5, query)
    return normalized in terms_for(5, annotation) or normalized in terms_for(5, profile)
