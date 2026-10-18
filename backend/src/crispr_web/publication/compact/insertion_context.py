"""Complete insertion evidence for stable module result consumers."""


def insertion_context_data(read):
    return {
        key: value
        for key, value in read.json("insertion-context.json").items()
        if key not in {"schema", "execution_completed"}
    }
