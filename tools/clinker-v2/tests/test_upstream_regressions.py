import numpy as np
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clinker.align import Globaligner, Alignment, Link, assign_groups, get_pairs
from clinker.classes import Cluster, Locus, Gene, merge_locations, shift_origin


def cluster(name):
    return Cluster(name, [Locus(name, [Gene(label=f"{name}{i}", translation="MKKAA") for i in range(2)])])


def test_alignment_does_not_mutate_caller_config():
    ga = Globaligner({"mode": "global", "substitution_matrix": "PAM250", "open_gap_score": -10, "extend_gap_score": -.5})
    before = dict(ga.aligner_config)
    ga.align_clusters(cluster("a"), cluster("b"))
    assert ga.aligner_config == before


def test_synteny_preserves_weight_and_is_symmetric():
    a, b = cluster("a"), cluster("b")
    ga = Globaligner()
    alignment = Alignment(query=a, target=b)
    for x, y in zip(a.loci[0].genes, b.loci[0].genes):
        alignment.add_link(x, y, 1, 1)
    ga.add_alignment(alignment)
    assert ga.synteny(a.uid, b.uid, i=.5) == 2.5
    assert np.array_equal(ga.matrix(i=3), [[0, 5], [5, 0]])
    assert list(get_pairs(a)) == [(0, 1)]


def test_bridge_link_merges_entire_groups():
    a, b, c, d = [Gene() for _ in range(4)]
    assign_groups([Link(query=a,target=b,identity=1), Link(query=c,target=d,identity=1), Link(query=b,target=c,identity=1)])
    assert len({g._group for g in [a,b,c,d]}) == 1


def test_merge_keeps_nonadjacent_and_odd_location_parts():
    parts = [FeatureLocation(0,2,1), FeatureLocation(4,6,1), FeatureLocation(6,9,1)]
    assert [(int(p.start), int(p.end)) for p in merge_locations(parts)] == [(0,2),(4,9)]


def test_alignment_roundtrip_retains_identity():
    a,b=cluster("a"),cluster("b")
    original=Alignment(query=a,target=b)
    assert Alignment.from_dict(original.to_dict()).uid == original.uid


def test_circular_multipart_shift_preserves_extracted_sequence():
    record=SeqRecord(Seq("ACGT"*25),id="circle")
    feature=SeqFeature(CompoundLocation([FeatureLocation(80,90,1),FeatureLocation(5,10,1)]),type="CDS",qualifiers={"locus_tag":["g"]})
    record.features=[feature]
    before=feature.extract(record.seq)
    shifted=shift_origin(record)
    assert shifted.features[-1].extract(shifted.seq)==before


def test_feature_fallback_uses_code_and_explicit_zero_start():
    record=SeqRecord(Seq("AATGTGATAA"),id="x")
    feature=SeqFeature(FeatureLocation(1,10,1),type="CDS",qualifiers={"transl_table":["4"]})
    gene=Gene.from_seqfeature(feature,record,start=0)
    assert gene.translation=="MW"
    assert gene.start==0
