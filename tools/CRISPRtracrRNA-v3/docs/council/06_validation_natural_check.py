"""Independent read-only recount of the frozen council natural panels."""
import copy
import csv
import hashlib
import json
import math
import statistics
import zipfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ROUND = ROOT / 'benchmarks/council_round2'

def read(p): return json.loads(p.read_text())
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def fasta(p):
    records = {}
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]; records[key] = ''
        elif line.strip(): records[key] += line.strip().upper().replace('U','T')
    return records

def rc(s): return s.translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]
def oriented(s, strand):
    assert strand in ['+','-']
    return s if strand=='+' else rc(s)
def interval(iv,n):
    assert 1 <= iv['start'] <= iv['end'] <= n
    assert iv['length']==iv['end']-iv['start']+1
    return iv['start']-1,iv['end']
def matched(p,ref,contig):
    if contig!=ref['contig'] or p['strand']!=ref['strand']: return 0
    iv=p['interval']; overlap=max(0,min(iv['end'],ref['end'])-max(iv['start'],ref['start'])+1)
    shorter=min(iv['end']-iv['start']+1,ref['end']-ref['start']+1)
    return overlap if overlap>=20 and overlap>=shorter/2 else 0

def normalize(p,experimental=False):
    p=copy.deepcopy(p)
    for k in ['admission_branch','evidence_class']:p.pop(k,None)
    if experimental:p.pop('export_id',None)
    if p.get('interaction'):
        for k in ['target_paired_positions','repeat_paired_positions']:p['interaction'].pop(k,None)
    return p

def geometry(p,source):
    s,e=interval(p['interval'],len(source)); seq=oriented(source[s:e],p['strand'])
    assert p['sequence']==seq
    for hit in p.get('model_hits',[])+p.get('terminator_hits',[]):
        a,b=interval(hit['interval'],len(seq));x,y=interval(hit['genomic_interval'],len(source))
        assert (x,y)==((s+a,s+b) if p['strand']=='+' else (e-b,e-a))
        assert oriented(seq[a:b],hit['strand'])==oriented(source[x:y],hit['genomic_strand'])
    anti=p.get('anti_repeat')
    if anti:
        a,b=interval(anti['local_interval'],len(seq));x,y=interval(anti['interval'],len(source))
        assert (x,y)==((s+a,s+b) if p['strand']=='+' else (e-b,e-a))
        assert anti['sequence']==seq[a:b]==oriented(source[x:y],p['strand'])

entries=read(ROOT/'benchmarks/literature/loci.json')
assert len(entries)==79 and sum(e['status']=='exact_unique' for e in entries)==78
native=read(ROUND/'native_execution_manifest.json')
source_archives={
    'archived_detector_pinned_native': ROUND/'baseline/source.zip',
    'current_default_pinned_native': ROUND/'current_default_pinned_native_source.zip',
    'current_experimental_pinned_native': ROUND/'current_experimental_pinned_native_source.zip',
}
stages={}; all_predictions={}; archive_checks={}; evidence_checks={}
for stage in ['historical','archived_detector_pinned_native','current_default_pinned_native','current_experimental_pinned_native']:
    rows=[]; lengths=[]; classes=Counter(); admissions=[]; count_evidence=0; predmap={}
    m=None if stage=='historical' else read(ROUND/stage/'run_manifest.json')
    if m:
        assert len(m['inputs'])==len(m['results'])==79
        assert [e['label'] for e in m['inputs']]==[e['label'] for e in entries]
        assert all(e['status']=='complete' for e in m['results'])
        assert m['settings']['threads']==m['threads_per_locus']==1 and m['workers']==2
        assert m['settings']['fasta36']==native['binary']
        assert m['native_execution_manifest_sha256']==digest(ROUND/'native_execution_manifest.json')
        x=m.get('executables',m.get('native_sha256')); assert x['fasta36']['sha256']==native['binary_sha256']
        with zipfile.ZipFile(source_archives[stage]) as z:
            assert all(hashlib.sha256(z.read(path)).hexdigest()==sha for path,sha in m['source_sha256'].items())
        archive_checks[stage]={'source_files':len(m['source_sha256']),'archive_sha256':digest(source_archives[stage])}
    for entry in entries:
        path=ROUND/'baseline/reports'/f"{entry['label']}.json" if stage=='historical' else ROUND/stage/entry['label']/'report.json'
        d=read(path); assert d['status']=='complete'; ranked=[]; flat=[]
        source=fasta(ROOT/entry['fasta']); assert digest(ROOT/entry['fasta'])==entry['fasta_sha256']
        for rec in d['records']:
            assert rec['source']['file_sha256']==entry['fasta_sha256']
            seq=source[rec['source']['id']]
            assert hashlib.sha256(seq.encode()).hexdigest()==rec['source']['sha256']
            for p in rec['predictions']:
                geometry(p,seq);ranked.append((rec['source']['id'],p));flat.append(p);lengths.append(len(p['sequence']))
                classes[p.get('evidence_class','historical_unspecified')]+=1
            for p in rec.get('short_rescue_proposals',[]):
                geometry(p,seq)
                if p['short_rescue_admission']['accepted']:
                    assert any(q['id']==p['id'] for q in rec['predictions'])
                    admissions.append({'label':entry['label'],'id':p['id'],'class':p['evidence_class']})
        assert len(flat)==d['summary']['prediction_count']
        if m:
            assert d['configuration']==m['settings']
            manifest=read(path.parent/'manifest.json')['files']
            for f,sha in manifest.items():assert digest(path.parent/f)==sha
            count_evidence+=len(manifest)
            csvrows=list(csv.DictReader((path.parent/'summary.csv').open()))
            seqrows=fasta(path.parent/'predictions.fasta')
            assert len(csvrows)==len(seqrows)==len(flat)
            for row,p in zip(csvrows,flat):
                assert row['prediction_id']==p['export_id']
                assert row['sequence']==seqrows[p['export_id']]==p['sequence']
                assert row['strand']==p['strand'] and int(row['start'])==p['interval']['start'] and int(row['end'])==p['interval']['end']
        predmap[entry['label']]=flat
        ranked.sort(key=lambda cp:(-cp[1]['ranking']['score'],cp[1]['id']))
        truth=entry['matches'][0] if entry['status']=='exact_unique' else None
        matches=[(rank,p,matched(p,truth,c)) for rank,(c,p) in enumerate(ranked,1) if truth and matched(p,truth,c)]
        rows.append({'label':entry['label'],'detected':bool(matches) if truth else None,'rank':matches[0][0] if matches else None,'coverage':max((o/(truth['end']-truth['start']+1) for _,p,o in matches),default=0) if truth else None,'count':len(ranked)})
    mapped=[r for r in rows if r['detected'] is not None];counts=[r['count'] for r in rows]
    stages[stage]={'complete':79,'localized':sum(r['detected'] for r in mapped),'top_k':{k:sum(r['rank'] is not None and r['rank']<=k for r in mapped) for k in [1,3,5,10]},'mrr':sum(1/r['rank'] if r['rank'] else 0 for r in mapped)/78,'coverage90':sum(r['coverage']>=.9 for r in mapped),'count':sum(counts),'median':statistics.median(counts),'p95':sorted(counts)[math.ceil(.95*len(counts))-1],'max':max(counts),'length_median':statistics.median(lengths),'length_p95':sorted(lengths)[math.ceil(.95*len(lengths))-1],'length_max':max(lengths),'admitted':admissions,'classes':dict(classes),'rows':rows,'batch_wall_seconds':m['batch_wall_seconds'] if m else None}
    evidence_checks[stage]=count_evidence;all_predictions[stage]=predmap
old=all_predictions['archived_detector_pinned_native'];new=all_predictions['current_default_pinned_native'];exp=all_predictions['current_experimental_pinned_native']
for label in old:
    assert [normalize(p) for p in old[label]]==[normalize(p) for p in new[label]]
    assert [normalize(p,True) for p in new[label]]==[normalize(p,True) for p in exp[label] if p['admission_branch']!='short_rescue']
    assert not any(p.get('admission_branch')=='short_rescue' for p in new[label])
def bylabel(stage):return {r['label']:r for r in stages[stage]['rows']}
hist=bylabel('historical');arch=bylabel('archived_detector_pinned_native');before=bylabel('current_default_pinned_native');after=bylabel('current_experimental_pinned_native')
required={label for label,r in hist.items() if r['detected'] is True}
assert len(required)==68
retention={stage: sorted(label for label in required if bylabel(stage)[label]['detected'] is not True) for stage in stages if stage!='historical'}
assert not any(retention.values())
gains=[label for label in before if after[label]['detected'] is True and before[label]['detected'] is not True]
losses=[label for label in before if before[label]['detected'] is True and after[label]['detected'] is not True]
deltas=[{'label':label,'old_count':hist[label]['count'],'archived_count':arch[label]['count'],'old_rank':hist[label]['rank'],'archived_rank':arch[label]['rank']} for label in hist if hist[label]!=arch[label]]
result={'source_archive_linkage':archive_checks,'evidence_hashes_verified':evidence_checks,'same_dependency_default_exact':True,'experimental_preserves_all_standard_predictions':True,'historical68_losses':retention,'dependency_deltas':deltas,'experimental_gained_labels':gains,'experimental_lost_labels':losses,'added_max_per_locus':max(after[k]['count']-before[k]['count'] for k in before),'wall_ratio':stages['current_experimental_pinned_native']['batch_wall_seconds']/stages['current_default_pinned_native']['batch_wall_seconds'],'stages':stages}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v if k!='stages' else {s:{a:b for a,b in data.items() if a!='rows'} for s,data in v.items()} for k,v in result.items()},indent=2))
