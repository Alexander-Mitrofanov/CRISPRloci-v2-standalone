"""Complete group/target support joined to submitted spacer origins."""


def association_data(read):
    doc = read.json("spacer-associations.json")
    source = read.json("input-sequences.json")
    return {
        key: doc[key]
        for key in (
            "reference",
            "requested_fdr",
            "calibration_status",
            "calibration",
            "groups",
            "associations",
            "counts",
            "interpretation",
        )
    } | {
        "occurrences": source["records"],
        "origin_semantics": "Submitted origin annotations are preserved, not independently verified host identities.",
    }
