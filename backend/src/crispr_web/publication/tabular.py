"""Spreadsheet shielding with narrow typed scientific-column exceptions."""
import re


def scientific_cell(value, *, strand=False, protein_alignment=False):
    if strand and value in {"+", "-", "both"}:
        return value
    if protein_alignment and isinstance(value, str) and re.fullmatch(r"[ABCDEFGHIKLMNPQRSTVWXYZJUO-]+", value):
        return value
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return value
