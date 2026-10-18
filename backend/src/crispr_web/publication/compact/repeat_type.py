"""Repeat subtype evidence joined to preserved occurrence origins."""


def repeat_type_data(read):
    doc = read.json("repeat-typing.json")
    source = read.json("input-sequences.json")
    origins = {row["id"]: row.get("origin") for row in source["records"]}
    return {"model": doc["model"], "policy": doc["policy"], "counts": doc["counts"],
            "results": [{**row, **({"origin": origins[row["id"]]} if origins.get(row["id"]) is not None else {})}
                        for row in doc["results"]]}
