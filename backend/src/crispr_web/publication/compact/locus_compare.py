"""Complete source-bound comparison in the compact public contract."""


def locus_compare_data(read):
    return {
        key: value
        for key, value in read.json("locus-comparison.json").items()
        if key not in {"schema", "execution_completed"}
    }
