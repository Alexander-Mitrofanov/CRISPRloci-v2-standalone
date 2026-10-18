"""Read-only validation of completed controls; incomplete panels remain not measured."""
import hashlib,json,statistics,zipfile
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
ROUND=ROOT/'benchmarks/council_round2'
load=lambda p:json.loads(p.read_text())
hashfile=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
control=load(ROUND/'controls/manifest.json');native=load(ROUND/'native_execution_manifest.json')
assert hashfile(ROUND/'controls/manifest.json')=='67efd2ede9053e2638de91c2646c08c157c347c13c833a8aa4a66b51725ebb55'
def fasta(p):
 out={}
 for l in p.read_text().splitlines():
  if l.startswith('>'):key=l[1:].split()[0];out[key]=''
  elif l.strip():out[key]+=l.strip().upper().replace('U','T')
 return out

def provenance(stage, m):
 assert m['native_dependency_profile']=='council_seed20260927_threads1'
 assert m['native_execution_manifest_sha256']==hashfile(ROUND/'native_execution_manifest.json')
 assert m['settings']['fasta36']==native['binary']
 assert m['native_sha256']['fasta36']['sha256']==native['binary_sha256']
 assert m['settings']['threads']==m['threads_per_locus']==1 and m['workers']==2
 assert m['settings']['short_anti_repeat_rescue']==('experimental' in stage or stage.startswith('conditional'))
 assert m['protocol_sha256']==hashfile(ROUND/'protocol.json')
 assert m['protocol_amendment_sha256']==hashfile(ROUND/'protocol_amendment_v2.json')
 archive=ROUND/(stage+'_source.zip')
 if not archive.exists(): archive=ROUND/'current_default_pinned_native_source.zip'
 if archive.exists():
  with zipfile.ZipFile(archive) as z:
   assert all(hashlib.sha256(z.read(p)).hexdigest()==h for p,h in m['source_sha256'].items())
  return {'source_archive':str(archive.relative_to(ROOT)),'source_archive_sha256':hashfile(archive)}
 # A per-run source digest is acceptable only if the current bytes still match;
 # final delivery should retain that version too.
 assert all(hashfile(ROOT/p)==h for p,h in m['source_sha256'].items())
 return {'source_current_bytes_match':True,'source_archive_not_yet_present':True}

result={'controls_manifest_sha256':hashfile(ROUND/'controls/manifest.json'),'panels':{}}
for stage in ['genomic_default','genomic_experimental','conditional_dinucleotide_1','conditional_dinucleotide_2','conditional_derangement']:
 path=ROUND/stage/'run_manifest.json'
 if not path.exists() or not load(path).get('completed_at'):
  result['panels'][stage]={'status':'not_measured','reason':'native panel incomplete'};continue
 m=load(path);artifact=provenance(stage,m);rows=[];nhashes=0
 if stage.startswith('genomic'):
  expected={r['label']:r for r in control['genomic_controls']}
  assert len(m['inputs'])==len(m['results'])==19 and all(r['status']=='complete' for r in m['results'])
  assert {r['label'] for r in m['inputs']}==set(expected)
  for entry in m['inputs']:
   label=entry['label'];f=ROOT/entry['fasta'];assert hashfile(f)==entry['fasta_sha256']==expected[label]['sha256']
   d=load(ROUND/stage/label/'report.json');assert d['status']=='complete' and d['configuration']==m['settings']
   predictions=[];proposals=[]
   for rec in d['records']:
    assert rec['source']['file_sha256']==entry['fasta_sha256'];predictions.extend(rec['predictions']);proposals.extend(rec.get('short_rescue_proposals',[]))
   for p,h in load(ROUND/stage/label/'manifest.json')['files'].items():assert hashfile(ROUND/stage/label/p)==h;nhashes+=1
   admitted=[p for p in proposals if p['short_rescue_admission']['accepted']]
   assert len(predictions)==d['summary']['prediction_count']
   assert sum(p['admission_branch']=='short_rescue' for p in predictions)==len(admitted)
   rows.append({'label':label,'predictions':len(predictions),'proposals':len(proposals),'admitted':len(admitted),'input_sha256':entry['fasta_sha256']})
  exported=load(ROUND/stage/'metrics.json')
  assert exported['complete']==19 and all(r['detected'] is None and r['coverage'] is None for r in exported['rows'])
 else:
  name=stage.removeprefix('conditional_');expected={}
  for r in control['repeat_controls']:
   if r['panel']==name:expected.setdefault(r['label'],{})[r['array_id']]=r['control_repeat']
  assert len(m['inputs'])==len(m['results'])==m['complete']==79
  assert all(r['status']=='complete' for r in m['results'])
  assert sum(len(v) for v in expected.values())==84
  assert m['control_manifest_sha256']==hashfile(ROUND/'controls/manifest.json')
  for entry in m['inputs']:
   label=entry['label'];work=ROUND/stage/label;d=load(work/'result.json')
   assert d['status']=='complete' and d['control_repeats']==expected[label] and d['arrays']==len(expected[label])
   assert list(fasta(ROOT/entry['fasta']).values()) == list(fasta(work/'input.fasta').values())
   for p,h in load(work/'manifest.json')['files'].items():assert hashfile(work/p)==h;nhashes+=1
   admitted=[p for p in d['proposals'] if p['short_rescue_admission']['accepted']]
   assert len(admitted)==d['short_admitted'] and len(d['proposals'])==d['short_proposals']
   rows.append({'label':label,'arrays':d['arrays'],'proposals':len(d['proposals']),'admitted':len(admitted),'admissions':[{'id':p['id'],'evidence_class':p['evidence_class'],'interval':p['interval'],'strand':p['strand']} for p in admitted]})
 counts=[r['admitted'] for r in rows]
 result['panels'][stage]={'status':'complete_verified','inputs':len(rows),'proposals':sum(r['proposals'] for r in rows),'admitted':sum(counts),'median_admitted':statistics.median(counts),'max_admitted':max(counts),'evidence_hashes_verified':nhashes,'rows':rows,**artifact}
 if stage.startswith('genomic'):result['panels'][stage]['predictions']=sum(r['predictions'] for r in rows)
 else:
  assert sum(counts)==m['short_admitted'] and statistics.median(counts)==m['short_median'] and max(counts)==m['short_max']
  result['panels'][stage]['burden_gate_passes']=statistics.median(counts)<=1 and max(counts)<=10
p=result['panels']
if all(p[n]['status']=='complete_verified' for n in ['genomic_default','genomic_experimental']):
 a={r['label']:r for r in p['genomic_default']['rows']};b={r['label']:r for r in p['genomic_experimental']['rows']}
 assert all(a[k]['input_sha256']==b[k]['input_sha256'] for k in a)
 result['paired_genomic']={'identical_inputs':True,'burden_increased_inputs':sum(b[k]['predictions']>a[k]['predictions'] for k in a),'default_total':p['genomic_default']['predictions'],'experimental_total':p['genomic_experimental']['predictions'],'new_admissions':p['genomic_experimental']['admitted']}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v if k!='panels' else {s:{a:b for a,b in d.items() if a!='rows'} for s,d in v.items()} for k,v in result.items()},indent=2))
