# Five contextual analyses for CRISPRloci

CRISPRloci 2.12.0 adds five workflows to help researchers interpret a detected
CRISPR–Cas locus, its mobile-element context and related RNA sequences. They
reuse the existing detection, annotation and sequence contracts. The API queues
jobs; isolated workers execute pinned scientific environments. These additions
follow the [literature and community-benefit council](../improvenments/10-literature-workflow-council-2026-10-12.md).

## Scientific purpose and outputs

| Workflow | Why it belongs in CRISPRloci | Result and interpretation |
| --- | --- | --- |
| **clinker locus comparison** | Researchers compare gene organization after identifying a candidate locus. CasAndra provides the source-bound CDSs and Cas cassettes, avoiding manual reannotation and coordinate conversion. | Protein alignments and gene-order maps for 2–8 submitted neighborhoods, each ≤50 kb. Similarity links are not proof of orthology, function or shared ancestry. Full protein alignments and source hashes are downloadable. |
| **Phigaro prophage context** | A spacer match becomes easier to interpret when researchers can see whether its location overlaps a candidate prophage. The workflow connects CRISPRidentify → source-bound CRISPRspacer search → Phigaro intervals. | Full-length exact matches on both strands within the same submitted record; matches overlapping reported arrays are excluded. Candidate prophage overlap does not establish active targeting, inducibility, infection or host range. Records shorter than 20 kb are explicitly unassessed by the Phigaro gene-calling model. |
| **ISEScan insertion context** | Researchers investigating unusual Cas neighborhoods need evidence for nearby insertion sequences before designing follow-up experiments. Source-coordinate joins place predicted elements relative to CasAndra cassettes. | IS family, complete/partial native call, transposase interval, available terminal inverted repeats, and nearest/overlapping Cas cassette coordinates. Target-site duplications are not assessed; overlap alone does not establish gene disruption or loss of function. |
| **DefenseFinder defense context** | CRISPR–Cas is one component of a sequence's defensive repertoire. Reporting other model-supported systems alongside Cas cassettes helps researchers select loci and compare submitted regions. | Source-bound defense-system components, HMM evidence and Cas cassette proximity. Model set 1.2.3 is explicitly pinned; later models and CasFinder are not assessed. No population enrichment, defense-island classification or functional activity is inferred. |
| **RNAalifold RNA families** | After curating homologous, consistently oriented tracrRNA candidates, researchers can evaluate a possible shared structural core without manually transferring candidates between alignment and folding programs. | MAFFT alignment and ViennaRNA comparative fold for 2–50 sequences of 30–500 nt, at least two distinct. Exact duplicates share one representative weight. Pair probabilities and descriptive covariation are shown; the method does not establish homology, transcript boundaries, structural significance or repeat–anti-repeat duplex pairing. |

All genomic context workflows operate on the submitted sequence coverage, capped
at 2 Mb. A zero call is not evidence of absence from an organism's complete
genome. RNA-family homology and 5′→3′ orientation are user-supplied assumptions;
subtype labels alone are insufficient to form a family. Missing values remain
missing. Coordinate-neighbor ties have an explicit count and truncation state.

Primary methods: [clinker](https://doi.org/10.1093/bioinformatics/btab007),
[Phigaro](https://doi.org/10.1093/bioinformatics/btaa250),
[ISEScan](https://doi.org/10.1093/bioinformatics/btx433),
[DefenseFinder](https://doi.org/10.1038/s41467-022-30269-9), and
[RNAalifold](https://doi.org/10.1186/1471-2105-9-474).
Application-study evidence and the distinction between community practice and
this implementation's validation are recorded in the linked literature review.

## License decisions

Each tool has a separate private repository under Alexander-Mitrofanov. Exact
source commits are in [the tool lock](../tools/tools.lock.json); original
copyright notices and applicable license texts remain with their source.

| Repository | Relevant terms | Integration decision |
| --- | --- | --- |
| [clinker-v2](https://github.com/Alexander-Mitrofanov/clinker-v2) | Upstream MIT; dependencies retain their own terms, including separately installed LGPL `argh`. | Run in its own worker environment. Preserve notices; render validated data in the server's interface. |
| [Phigaro-v2](https://github.com/Alexander-Mitrofanov/Phigaro-v2) | MIT source; separate Prodigal GPL-3.0-or-later and HMMER BSD-3-Clause executables. | Hosted research use with attribution. The pVOG data are operator-provisioned and are not relicensed as MIT or committed to the repository. |
| [ISEScan-v2](https://github.com/Alexander-Mitrofanov/ISEScan-v2) | Apache-2.0 source/models as distributed upstream; SSW and kseq MIT. FragGeneScan 1.32's actual source notice says GNU GPL with no version specified, despite conflicting Bioconda BSD metadata. | Preserve the authoritative source notice and run native dependencies separately. No independently acquired ISfinder database is bundled. |
| [DefenseFinder-v2](https://github.com/Alexander-Mitrofanov/DefenseFinder-v2) | GPL-3.0 source and the selected 1.2.3 model release. Later model metadata contains conflicting GPL and noncommercial declarations. | Use the unambiguous older model release, preserve its exact identity and disclose its narrower coverage. Do not silently update to the disputed releases. |
| [RNAalifold-v2](https://github.com/Alexander-Mitrofanov/RNAalifold-v2) | New wrapper MIT; unmodified ViennaRNA 2.7 retains its custom license. It permits research/educational/commercial use with credit and restrictions on fee-based redistribution; commercial-product inclusion asks users to contact the authors. MAFFT core retains its BSD-style terms. | Suitable for the intended free research service with University of Vienna/authors' credit. The wrapper's MIT license does not relicense ViennaRNA. Commercial distribution requires separate consideration of the original conditions. |

Authoritative notices are retained in the repositories. For data,
[NCBI's policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) places no
additional NCBI restriction on molecular-data reuse while retaining a
third-party-rights caveat. Reference acquisition preserves source, release and
hash provenance. Any later distribution of worker binaries must also meet their
corresponding-source and notice obligations; the public web export contains only
frontend assets.

## Implementation and verification

The adapters invoke shell-free, operator-owned commands and validate source
binding before producing pipeline output. The backend repeats validation through
its own code and original submitted inputs, without importing Integration or
scientific-package validators. Worker cancellation, deadlines, per-job
capabilities, bounded reads and retention remain in force. Model acquisition and
updates are operator actions, never job actions.

Known defects fixed include clinker alignment-option mutation and clustering
bridges, Phigaro terminal-region handling and final-gene slicing, ISEScan native
command quoting/error propagation and source-relative data paths, and duplicate
DefenseFinder system summaries for shared proteins. The RNA wrapper adds explicit
deduplication, empty-input handling and DNA/RNA source reconstruction. Complete
CDSs require appropriate start/stop evidence. Protected exports are atomically
installed so tool-created symlinks cannot redirect postprocessing writes.

Phigaro retains the native Prodigal GC covariate used for its published threshold
calibration and separately reports exact source-CDS GC. The two can differ for
reverse-strand calls in Prodigal 2.6.3. Native executables are pinned and tested;
an independently compiled Prodigal that corrupted header metadata was rejected.

Local native pilots cover genuine executable paths: CasAndra→clinker protein
alignments; a positive and zero-call ISEScan control; DefenseFinder's RM Type II
system; Phigaro's lambda reference with 61 CDSs and one candidate region; and
MAFFT/ViennaRNA paired and unpaired families. A synthetic spacer fixture verifies
an array exclusion plus both-strand external matches and interval joins. These
are execution and evidence checks, not a new sensitivity/specificity benchmark
or a guarantee that every upstream bug has been eliminated.

The release gate additionally runs positive and zero-evidence controls through
the packaged HTTP API, protected downloads and independent validators for all
five workflows, plus invalid-input checks and all existing scientific modes.
Only the exact candidate that passes the full application, browser, native and
service-containment checks may be activated. Deployment outcomes and connection
records stay in ignored private operations storage.

CDS evidence retains ambiguity explicitly. Phigaro matches native Prodigal
translation: every codon containing an ambiguous base is reported as X, including
GCN and AAR, while source-DNA hashes retain the submitted symbols. DefenseFinder
CDS hashes bind the
exact oriented submitted DNA, including IUPAC symbols, while Pyrodigal protein
translation reports X for codons containing ambiguous bases. CasAndra's native
CDS nucleotide fields instead contain the caller's non-ACGT-to-N normalization;
its exact source FASTA hash remains separately bound and the independent backend
checks that normalization against the submitted coordinates. Neither workflow
reconstructs an amino acid by resolving an ambiguous codon.

## Effective context policies

Phigaro uses joint Prodigal single-mode training with genetic code 11. Submit a genome/assembly or closely related strains; combining unrelated genomes can alter predictions for unchanged records. The 20 kb minimum is an execution eligibility threshold, not an independently validated accuracy cutoff. Native `gc` and source-derived `source_gc` remain distinct.

Full and compact prophage-context evidence publishes `search_policy`: same-record, full-length exact ACGT spacers of 16–1000 bases from Bona-fide/Possible arrays, both strands, zero substitutions/gaps and overlap exclusion against all reported arrays with zero padding. Palindromes have strand `both`. Completed search has no result truncation; PAM is not evaluated and targeting activity is not assessed. Short sources are searched independently of prophage eligibility. Excluded raw work rows are not promised in downloads. Arrays/repeats/spacers reuse validated detection evidence through the existing typed ports.

DefenseFinder rejects ambient configuration paths consumed by pinned MacSyFinder 2.1.4 before scientific work, including selected global/environment, user, local and selected family model settings. Native global FASTA rank and protein length are checked independently of per-source gene position. RNAalifold strips inherited `MAFFT_BINARIES` only from child environments; managed MAFFT is 7.526. The earlier 7.505 example run is historical, not a managed-release identity.

Current derivative source identities after the 14 October engineering repairs (the [lock](../tools/tools.lock.json) is authoritative):

- [evOr 1.0.4](https://github.com/Alexander-Mitrofanov/evOr/tree/c46f6f9be2ef80f26168197d25a50b6fd645a24e), commit `c46f6f9be2ef80f26168197d25a50b6fd645a24e`.
- [clinker-v2 2.0.1](https://github.com/Alexander-Mitrofanov/clinker-v2/tree/249881b31c5a062ef47efd9c061b00a2831f8a63), commit `249881b31c5a062ef47efd9c061b00a2831f8a63`.
- [DefenseFinder-v2 2.0.1](https://github.com/Alexander-Mitrofanov/DefenseFinder-v2/tree/f579941ca35fb4bd26fcd4c35fcffc23df9d91a2), commit `f579941ca35fb4bd26fcd4c35fcffc23df9d91a2`.
- [RNAalifold-v2 1.0.1](https://github.com/Alexander-Mitrofanov/RNAalifold-v2/tree/0498c888203ba7bd6a8e6c6e35a27229adc61d1e), commit `0498c888203ba7bd6a8e6c6e35a27229adc61d1e`.

The independent readers explicitly accept the preceding and current patch report versions for clinker, DefenseFinder and RNAalifold. Retained historical fixtures keep their original tool identities; current executions and release manifests identify the new pinned source.
