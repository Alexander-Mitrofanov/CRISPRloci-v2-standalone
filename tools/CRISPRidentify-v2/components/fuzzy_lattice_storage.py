"""Lossless, process-local storage for idle immutable fuzzy lattices.

References own their compressed records. A small shared cache retains decoded
records only; it cannot keep an unbounded registry of discarded references.
This changes storage, never the set of nodes, paths or completeness evidence.
Public candidate serialization continues to expose the original record types.
"""
from collections import OrderedDict
from dataclasses import dataclass, fields
import hashlib
import io
import pickle
import weakref
import zlib

from components.fuzzy_boundary_lattice import (
    FuzzyBoundaryDomain, FuzzyBoundaryLattice, FuzzyBoundaryNode,
    FuzzyBoundaryPath, FuzzyBoundarySearchPolicy, FuzzyBoundarySeed,
    FuzzyEditLimits,
)
from components.fuzzy_occurrence_paths import (
    OccurrencePathLattice, OccurrencePathNode, OccurrencePathRecord,
)

_DATACLASSES = frozenset((FuzzyBoundaryDomain, FuzzyBoundaryLattice,
                        FuzzyBoundaryNode, FuzzyBoundaryPath,
                        FuzzyBoundarySearchPolicy, FuzzyBoundarySeed,
                        FuzzyEditLimits))
_TUPLES = frozenset((tuple, OccurrencePathLattice,
                    OccurrencePathNode, OccurrencePathRecord))
_ROOTS = frozenset((FuzzyBoundaryLattice, OccurrencePathLattice))
_CLASSES = {(cls.__module__, cls.__name__): cls
            for cls in (_DATACLASSES | _TUPLES) if cls is not tuple}


def _require_immutable(record):
    seen, pending = set(), [record]
    while pending:
        value = pending.pop()
        if id(value) in seen:
            continue
        seen.add(id(value))
        kind = type(value)
        if kind in (str, int, bool, type(None)):
            continue
        if kind in _TUPLES:
            pending.extend(value)
        elif kind in _DATACLASSES:
            pending.extend(getattr(value, field.name) for field in fields(value))
        else:
            raise TypeError("fuzzy lattice storage requires known immutable records")


class _LatticeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            return _CLASSES[(module, name)]
        except KeyError:
            raise pickle.UnpicklingError("unsupported class in private lattice record") from None


@dataclass(frozen=True, eq=False, slots=True, weakref_slot=True)
class CompressedLattice:
    payload: bytes
    digest: str
    record_type: type
    cache: object

    def resolve(self):
        return self.cache.resolve(self)


class LatticeDecodeCache:
    """Keep at most four decoded records, excluding caller-held references.

    This is an entry bound, not a byte bound. Eviction loses only reuse; the
    complete compressed record stays with each live candidate reference.
    """

    def __init__(self, capacity=4):
        if type(capacity) is not int or capacity < 1:
            raise ValueError("lattice decode cache capacity must be positive")
        self.capacity = capacity
        self._decoded = OrderedDict()

    def compact(self, lattice):
        if type(lattice) not in _ROOTS:
            raise TypeError("unsupported fuzzy lattice record")
        _require_immutable(lattice)
        payload = zlib.compress(pickle.dumps(lattice, protocol=5), level=1)
        return CompressedLattice(payload, hashlib.sha256(payload).hexdigest(),
                                 type(lattice), self)

    def resolve(self, reference):
        identity = id(reference)
        entry = self._decoded.get(identity)
        if entry is not None and entry[0]() is reference and (
                entry[1] is reference.payload and entry[2] == reference.digest
                and entry[3] is reference.record_type):
            self._decoded.move_to_end(identity)
            return entry[4]
        if hashlib.sha256(reference.payload).hexdigest() != reference.digest:
            raise ValueError("compressed fuzzy lattice digest mismatch")
        if reference.record_type not in _ROOTS:
            raise TypeError("unsupported compressed fuzzy lattice type")
        self._decoded.pop(identity, None)
        while len(self._decoded) >= self.capacity:
            self._decoded.popitem(last=False)
        decoder = zlib.decompressobj()
        decoded = decoder.decompress(reference.payload) + decoder.flush()
        if not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
            raise ValueError("compressed fuzzy lattice framing mismatch")
        stream = io.BytesIO(decoded)
        value = _LatticeUnpickler(stream).load()
        if stream.read(1) or type(value) is not reference.record_type:
            raise ValueError("compressed fuzzy lattice type or framing mismatch")
        _require_immutable(value)

        def discard(weak):
            current = self._decoded.get(identity)
            if current is not None and current[0] is weak:
                self._decoded.pop(identity)

        self._decoded[identity] = (weakref.ref(reference, discard),
                                   reference.payload, reference.digest,
                                   reference.record_type, value)
        return value
