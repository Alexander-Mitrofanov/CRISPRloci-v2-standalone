"""Complete validated prophage context for module-result consumers."""


def prophage_context_data(read):
    return {
        k: v
        for k, v in read.json("prophage-context.json").items()
        if k not in {"schema", "execution_completed"}
    }
