"""Complete observed comparison, preserving source and occurrence annotations."""


def array_compare_data(read):
    document = read.json("array-comparison.json")
    return {
        key: value
        for key, value in document.items()
        if key not in {"schema", "execution_completed"}
    }
