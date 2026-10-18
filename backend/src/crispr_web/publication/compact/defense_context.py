"""Complete defense evidence retains source, protein and profile identities."""


def defense_context_data(read):
    doc = read.json("defense-context.json")
    return {key: doc[key] for key in ("cas_context", "model", "policy", "counts", "sources", "genes", "systems", "components", "hits")}
