"""Bounded source reads and explicit field projection for compact exports."""

from ...exchange import MAX_EXCHANGE_BYTES, strict_json
from ...summaries.bounded_io import _read_bounded


class Reader:
    def __init__(self, root, checkpoint=None):
        self.root = root
        self.checkpoint = checkpoint
        self.remaining = 2 * MAX_EXCHANGE_BYTES
        self.documents = {}

    def bytes(self, relative):
        payload = _read_bounded(
            self.root,
            self.root / relative,
            max_bytes=min(self.remaining, MAX_EXCHANGE_BYTES),
            checkpoint=self.checkpoint,
        )
        self.remaining -= len(payload)
        return payload

    def json(self, relative):
        if relative not in self.documents:
            self.documents[relative] = strict_json(self.bytes(relative))
        self.step()
        return self.documents[relative]

    def step(self):
        if self.checkpoint:
            self.checkpoint()


def fields(row, names):
    """Keep named scientific values, including explicit null/false/zero."""
    if not isinstance(row, dict):
        raise ValueError("expected a scientific record")
    return {key: row[key] for key in names.split() if key in row}


def sources(rows):
    return [fields(row, "id length sha256 topology") for row in rows]


def register(table, identity, row):
    if identity in table and table[identity] != row:
        raise ValueError("inconsistent scientific identity")
    table[identity] = row


def array_data(read):
    doc = read.json("arrays.json")
    return fields(doc, "coordinate_system score_semantics categories sources arrays")
