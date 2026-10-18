# Original CRISPRmap evidence and provenance

Audit date: 30 September 2026. The sources below concern the Freiburg
prokaryotic repeat-classification tool described by Lange et al. (2013), not the
unrelated optical-screening method with the same name.

## What was recovered

The official [REPEATS FASTA](https://rna.informatik.uni-freiburg.de/CRISPRmap/oridata/CRISPRmap_data.fa)
contains 3,527 records, 2,106 distinct oriented sequence strings and lengths
19–55 nucleotides. Its identifiers encode historical motif and family labels;
`X` means no annotation. This is the **2013** reference with 40 families and 33
motifs. It is not the 2014 reference with 4,719 records, 24 families and 18 motifs.
The [official documentation](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRmap)
explains those different partitions.

The original [16-query example](https://rna.informatik.uni-freiburg.de/CRISPRmap/Result.jsp?toolName=CRISPRmap&jobID=8888684)
retains its input FASTA and `Results.tab`. Its archived configuration identifies
`v1.0-2013` and orientation optimization enabled. The downloaded input and output
agree up to the original orientation reversal for query 13. Fourteen output
rows declare database membership; two are treated as new queries by the original
server. Our benchmark pairs original input and output by ordinal, preserving the
full input identifiers rather than the original UI's truncated identifiers.

| Frozen asset | Bytes | SHA256 |
| --- | ---: | --- |
| 2013 REPEATS FASTA | 305908 | `87d0c59372648ef1a2fbb8dd0bc06e887dbf1fbc5d4615e3e02eb488801c972e` |
| Original example input FASTA | 3781 | `c6c99cba8b451212db420604104af566f31d429d49a294482d8e8f1fc09e2cb8` |
| Original example Results.tab | 1407 | `26ae52c701eacc5c32a5ce123eba94951fdacfb3274d4517b3c75fcfa5efbe25` |

The historical [eight-query 2014 example](https://rna.informatik.uni-freiburg.de/CRISPRmap/Result.jsp?toolName=CRISPRmap&jobID=0788874)
was also inspected. Its configuration identifies `v2.1-2014`, with orientation
optimization disabled. It is excluded from this benchmark because we have not
recovered its complete reference. Substituting the 2013 reference would mix class
namespaces and invalidate a direct class comparison.

## Original execution is unavailable

On the audit date, the official
[input page](https://rna.informatik.uni-freiburg.de/CRISPRmap/Input.jsp) reported
that the service was offline because of incompatibility with its newer
infrastructure. The published article and official documentation provided a web
service and downloadable reference data. The repository search and official
software links examined in this audit did not yield a redistributable local
CRISPRmap executable, source release or commit. Saved original configuration
identifies a historical shell entrypoint but provides no executable content.

Consequently, the benchmark compares CRISPRmap-v2 with **authentic published
original outputs**. It does not rerun the original executable. Original runtime,
peak memory, throughput and speedup remain unavailable. No newly implemented
nearest-neighbor method is described as the original CRISPRmap algorithm.

## Acquisition and rights separation

`python scripts/fetch_legacy.py` is an explicit operator action. It downloads
only the three pinned assets, verifies their byte counts and SHA256 digests,
parses the label schema and creates a local reference manifest. Runtime search
is offline and never invokes this script. Modified upstream files are rejected,
not silently adopted as a new reference release.

The [2013 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3783184/) is distributed
under CC BY 3.0. The downloadable server assets do not carry a separately
verified data/software license in the files we recovered. The article's license
is not assumed to license unpublished server code. This repository therefore
contains the downloader, provenance and aggregate benchmark evidence, and does
not redistribute downloaded raw assets. They remain in ignored `data/legacy/`.
The MIT license applies to this project's new source code, not third-party
reference assets.

Historical labels are reference annotations, not independent experimental truth.
The original example is small and substantially overlaps the reference. Its
concordance cannot validate current biological coverage or establish accuracy on
modern metagenomes. See [SCIENCE.md](SCIENCE.md) for the interpretation of the
separate held-out and synthetic-control evaluations.
