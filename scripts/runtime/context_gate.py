"""Real-tool positive and negative controls for the five context workflows."""

import json

MODES = {
    "locus_compare",
    "insertion_context",
    "defense_context",
    "prophage_context",
    "rna_family",
}
EXPORTS = {
    "locus_compare": "locus-comparison.json",
    "insertion_context": "insertion-context.json",
    "defense_context": "defense-context.json",
    "prophage_context": "prophage-context.json",
    "rna_family": "rna-family.json",
}
CONTROL_CASES = {
    "clinker": ("locus_compare", "locus-positive", "locus-negative"),
    "isescan": ("insertion_context", "insertion-positive", "insertion-negative"),
    "defensefinder": ("defense_context", "defense-positive", "defense-negative"),
    "phigaro": ("prophage_context", "prophage-positive", "prophage-negative", "prophage-ineligible"),
    "rnaalifold": ("rna_family", "rna-positive", "rna-negative"),
}


def require_coverage(artifacts, gate):
    required = {
        (mode, fixture)
        for runner, (mode, *fixtures) in CONTROL_CASES.items()
        if runner in artifacts
        for fixture in fixtures
    }
    completed = {
        (item.get("mode"), item.get("fixture"))
        for item in gate.get("cases", [])
        if item.get("status") == "completed"
    }
    if required and (
        not required <= completed
        or gate.get("context_predispatch_rejections", 0) < len(set(artifacts) & set(CONTROL_CASES))
    ):
        raise ValueError("Context release gate is missing native controls or input rejections")


def fixtures(toolset, run_root, runners):
    paths = {}
    cases = []

    def write(name, text):
        path = run_root / (name + ".fasta")
        path.write_text(text, encoding="ascii")
        paths[name] = path

    # Empty-gene computational controls are deliberately shorter than a CDS.
    # Kilobase homopolymers are unsuitable: Pyrodigal calls partial ORFs in them.
    if "clinker" in runners:
        source = toolset / "sources/CRISPRtracrRNA-v3/examples/complete.fasta"
        lines = []
        started = False
        for line in source.read_text().splitlines():
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                lines.append(line.strip())
        seq = "".join(lines)[:10000]
        if len(seq) != 10000:
            raise ValueError("Clinker control needs 10kb from one pinned source record")
        write("locus-positive", f">neighborhood1\n{seq}\n>neighborhood2\n{seq}\n")
        write(
            "locus-negative",
            ">no_genes1\nACGTACGTACGT\n>no_genes2\nTGCATGCATGCA\n",
        )
        cases.extend(
            [
                ("locus_compare", "locus-positive", "II"),
                ("locus_compare", "locus-negative", "II"),
            ]
        )
    if "isescan" in runners:
        paths["insertion-positive"] = toolset / "sources/ISEScan-v2/tests/NC_012624.fna"
        write("insertion-negative", ">negative\n" + "ACGT" * 250 + "\n")
        cases.extend(
            [
                ("insertion_context", "insertion-positive", "II"),
                ("insertion_context", "insertion-negative", "II"),
            ]
        )
    if "defensefinder" in runners:
        paths["defense-positive"] = (
            toolset / "sources/DefenseFinder-v2/tests/data/context_fragment.fna"
        )
        write("defense-negative", ">negative\nACGTACGTACGT\n")
        cases.extend(
            [
                ("defense_context", "defense-positive", "II"),
                ("defense_context", "defense-negative", "II"),
            ]
        )
    if "phigaro" in runners:
        paths["prophage-positive"] = (
            toolset / "sources/Phigaro-v2/tests/fixtures/lambda.fasta"
        )
        # Confirmed with the pinned single-genome Prodigal policy: one partial
        # ORF, no pVOG-supported region. Short input is a separate unassessed case.
        write("prophage-negative", ">negative\n" + "A" * 20000 + "\n")
        write("prophage-ineligible", ">short\n" + "A" * 1000 + "\n")
        cases.extend(
            [
                ("prophage_context", "prophage-positive", "II"),
                ("prophage_context", "prophage-negative", "II"),
                ("prophage_context", "prophage-ineligible", "II"),
            ]
        )
    if "rnaalifold" in runners:
        # Synthetic computational control from the pinned tool repository.
        paths["rna-positive"] = (
            toolset / "sources/RNAalifold-v2/tests/data/family.fasta"
        )
        write(
            "rna-negative",
            ">unpaired1\n" + "A" * 40 + "\n>unpaired2\n" + "A" * 39 + "C\n",
        )
        cases.extend(
            [("rna_family", "rna-positive", "II"), ("rna_family", "rna-negative", "II")]
        )
    if any(not path.is_file() for path in paths.values()):
        raise ValueError("A pinned context release fixture is missing")
    return paths, cases


def rejected_submissions(runners):
    """Inputs that must fail before any scientific process is dispatched."""
    result = []
    if "clinker" in runners:
        result.append({"mode": "locus_compare", "sequence": ">single\n" + "ACGT" * 250})
    if "rnaalifold" in runners:
        result.append({"mode": "rna_family", "molecule": "RNA", "sequence": ">a\n" + "A" * 40 + "\n>b\n" + "A" * 40})
    for runner, mode in (("isescan", "insertion_context"), ("defensefinder", "defense_context"), ("phigaro", "prophage_context")):
        if runner in runners:
            result.append({"mode": mode, "molecule": "RNA", "sequence": ">invalid_molecule\n" + "ACGU" * 250})
    return result


def check_summary(mode, fixture, summary):
    if mode not in MODES:
        return
    doc = summary[mode]
    if doc.get("available") is not True:
        raise RuntimeError("Context result unavailable: " + mode)
    positive = fixture.endswith("positive")
    counts = doc["counts"]
    if mode == "insertion_context":
        if counts["elements"] != (1 if positive else 0):
            raise RuntimeError("ISEScan control count differs")
        if positive:
            row = doc["elements"][0]
            if (
                (row["family"], row["start"], row["end"], row["native_type"])
                != ("IS5", 30355, 31343, "p")
                or row["target_site_duplication"] is not None
                or row["terminal_inverted_repeats"] is None
            ):
                raise RuntimeError("ISEScan source-bound positive evidence differs")
    if mode == "defense_context":
        expected = (
            {"sources": 1, "genes": 42, "systems": 1, "components": 2, "hits": 5}
            if positive
            else {"sources": 1, "genes": 0, "systems": 0, "components": 0, "hits": 0}
        )
        if counts != expected or doc["policy"]["functional_activity"] != "not_assessed":
            raise RuntimeError("DefenseFinder control or uncertainty differs")
        if positive and doc["systems"][0]["subtype"] != "RM_Type_II":
            raise RuntimeError("DefenseFinder positive model differs")
    if mode == "locus_compare":
        if (
            counts["loci"] != 2
            or (positive and counts["links"] < 1)
            or (not positive and (counts["genes"] or counts["links"]))
        ):
            raise RuntimeError(
                "Clinker positive/negative neighborhood evidence differs"
            )
    if mode == "rna_family":
        if (
            counts["sources"] != (4 if positive else 2)
            or counts["unique_sequences"] != (3 if positive else 2)
        ):
            raise RuntimeError("RNA-family source/deduplication evidence differs")
        if (positive and counts["mfe_pairs"] < 1) or (
            not positive and counts["mfe_pairs"] != 0
        ):
            raise RuntimeError("RNA-family pairing control differs")
        if (
            doc["engine"]["version"] != "2.7.0"
            or doc["policy"]["activity"] != "not_assessed"
        ):
            raise RuntimeError("RNA-family engine or interpretation differs")
    if mode == "prophage_context":
        ineligible = fixture == "prophage-ineligible"
        expected = (0, 0, 0) if ineligible else ((1, 61, 1) if positive else (1, 1, 0))
        if (
            (counts["assessed_sources"], counts["genes"], counts["regions"]) != expected
            or doc["parameters"]["gene_mode"] != "single"
        ):
            raise RuntimeError("Phigaro assessment, gene policy or region control differs")


def check_export(mode, bundle, summary):
    if mode not in MODES:
        return
    report = json.loads(bundle.read("results/" + EXPORTS[mode]))
    if report["counts"] != summary[mode]["counts"]:
        raise RuntimeError(
            "Context protected export differs from validated summary: " + mode
        )
    forbidden = (
        "results/locus-compare/",
        "results/locus-compare-input/",
        "results/isescan/",
        "results/insertion-input/",
        "results/defense-context/",
        "results/defense-input/",
        "results/phigaro/",
        "results/phigaro-input/",
        "results/source-spacer/",
        "results/prophage-input/",
        "results/rna-family/",
        "results/rna-family-input/",
    )
    if any(name.startswith(forbidden) for name in bundle.namelist()):
        raise RuntimeError("Private context native files leaked")
