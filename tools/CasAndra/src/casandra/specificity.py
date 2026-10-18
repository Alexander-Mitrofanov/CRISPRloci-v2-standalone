"""Optional, hash-verified protein specificity overlay (local successor runtime).

Core and auxiliary HMM libraries remain separate, at their frozen search spaces.
No profile, hierarchy, score or coverage is retrained or silently repaired here.
"""
from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

from casandra._specificity_prior import apply as apply_prior_floor
from casandra._specificity_competition import apply_competition
from casandra._specificity_v3 import apply_cumulative, validate_card


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def domain_record(domain):
    a = domain.alignment
    return {'score': float(domain.score), 'independent_evalue': float(domain.i_evalue),
            'conditional_evalue': float(domain.c_evalue), 'reported': bool(domain.reported),
            'included': bool(domain.included), 'env_from': domain.env_from, 'env_to': domain.env_to,
            'hmm_from': a.hmm_from, 'hmm_to': a.hmm_to, 'hmm_length': a.hmm_length,
            'target_from': a.target_from, 'target_to': a.target_to, 'target_length': a.target_length}


def hit_record(hit, index, metadata=None, *, old=False):
    pid = _text(hit.name)
    r = {'reference_id' if old else 'profile': pid, 'ordinal': index,
         'score': float(hit.score), 'evalue': float(hit.evalue),
         'pre_score': float(hit.pre_score), 'sum_score': float(hit.sum_score), 'bias': float(hit.bias),
         'reported': bool(hit.reported), 'included': bool(hit.included),
         'domains': [domain_record(d) for d in hit.domains]}
    if metadata is not None:
        m = metadata[pid]
        r.update(intended_use=m['intended_use'], labels=copy.deepcopy(m['labels']),
                 families=copy.deepcopy(m.get('families', [])))
    if old:
        r.pop('ordinal')
        r.pop('included')
    return r


@dataclass(frozen=True)
class SpecificityOverlay:
    directory: Path
    manifest: dict
    floor: dict
    inherited_card: dict
    card: dict
    inherited_hmms: tuple
    new_hmms: tuple
    extensions: tuple = ()

    @classmethod
    def from_protein_directory(cls, directory, configuration, metadata, core_card, core_hmms):
        directory = Path(directory).resolve()
        if type(configuration.get('schema_version')) is not int or configuration['schema_version'] != 1:
            raise ValueError('Unsupported specificity overlay schema')
        name = configuration.get('manifest')
        if not isinstance(name, str):
            raise ValueError('Specificity manifest path is required')
        path = (directory/name).resolve()
        if not path.is_relative_to(directory) or not path.is_file():
            raise ValueError('Specificity manifest escapes the verified model')
        if _sha(path) != configuration.get('manifest_sha256'):
            raise ValueError('Specificity manifest checksum mismatch')
        manifest = json.loads(path.read_text())
        if type(manifest.get('schema_version')) is not int or manifest['schema_version'] != 1 or manifest.get('mode') not in {'v3', 'margin10'}:
            raise ValueError('Unknown specificity implementation or mode')
        if manifest.get('search_spaces') != {'core': 361, 'inherited': 370, 'new': 387}:
            raise ValueError('Frozen specificity search spaces mismatch')
        if len(core_hmms) != 361 or core_card['report_evalue'] != 1000.0 or core_card['detection_threshold'] != 0.0:
            raise ValueError('Specificity requires the unchanged full361 core settings')
        base = path.parent
        artifacts = manifest.get('artifacts', {})
        required = {'inherited_floor.json','inherited_competition.json','repair_policy.json',
                    'inherited_profiles.hmm','new_profiles.hmm','new_profile_metadata.json'}
        if not required.issubset(artifacts):
            raise ValueError('Specificity manifest is incomplete')
        for relative, pin in artifacts.items():
            target = (base/relative).resolve()
            if not target.is_relative_to(base) or not target.is_file() or target.stat().st_size != pin['size'] or _sha(target) != pin['sha256']:
                raise ValueError('Specificity artifact checksum mismatch: '+relative)
        for relative, digest in manifest['core_files_sha256'].items():
            target = (directory/relative).resolve()
            if not target.is_relative_to(directory) or _sha(target) != digest:
                raise ValueError('Core library differs from the frozen specificity parent')
        load = lambda name: json.loads((base/name).read_text())
        floor, inherited, card = load('inherited_floor.json'),load('inherited_competition.json'),load('repair_policy.json')
        validate_card(card)
        if inherited['Z'] != 370 or card['auxiliary_Z'] != 387 or card['inherited_auxiliary_Z'] != 370:
            raise ValueError('Specificity card search spaces mismatch')
        primary = {pid for pid,m in metadata.items() if m['intended_use']=='primary_search'}
        if set(card['positive_profile_registry']) != primary or set(card['profiles']) != primary:
            raise ValueError('Specificity positive profile registry mismatch')
        if _sha(base/'inherited_profiles.hmm') != inherited['hmm_sha256'] or _sha(base/'new_profiles.hmm') != card['negative_bank_hmm_sha256']:
            raise ValueError('Specificity HMM/card identity mismatch')
        if _sha(base/'new_profile_metadata.json') != card['negative_bank_metadata_sha256']:
            raise ValueError('Specificity metadata/card identity mismatch')
        with pyhmmer.plan7.HMMFile(base/'inherited_profiles.hmm') as handle:
            old = tuple(handle)
        with pyhmmer.plan7.HMMFile(base/'new_profiles.hmm') as handle:
            new = tuple(handle)
        old_ids,new_ids = [_text(h.name) for h in old],[_text(h.name) for h in new]
        if len(old)!=9 or len(set(old_ids))!=9 or set(inherited['allowed_references'])-set(old_ids):
            raise ValueError('Inherited reference registry mismatch')
        if len(new)!=26 or len(set(new_ids))!=26 or set(new_ids)!=set(card['negative_profiles']):
            raise ValueError('New reference registry mismatch')
        pointers=manifest.get('extensions',[])
        if not isinstance(pointers,list):
            raise ValueError('Additional specificity bank pointers must be a list')
        from casandra.specificity_extension import AdditionalSpecificityBank
        parent={'mode':manifest['mode'],'repair_policy_sha256':_sha(base/'repair_policy.json'),
            'inherited_floor_sha256':_sha(base/'inherited_floor.json'),
            'inherited_competition_sha256':_sha(base/'inherited_competition.json')}
        extensions=tuple(AdditionalSpecificityBank.load(base,p,metadata,parent) for p in pointers)
        if len({x.card['id'] for x in extensions}) != len(extensions):
            raise ValueError('Additional specificity policy IDs must be unique')
        return cls(base,manifest,floor,inherited,card,old,new,extensions)

    def _scan(self, queries, hmms, Z, threads, *, old=False):
        alphabet = pyhmmer.easel.Alphabet.amino()
        digital = [pyhmmer.easel.TextSequence(name=q['sequence_id'].encode(),sequence=q['sequence']).digitize(alphabet) for q in queries]
        stream = pyhmmer.hmmer.hmmscan(digital,list(hmms),cpus=threads,Z=Z,E=1000,domE=1000,incE=1000,incdomE=1000)
        result=[]
        for q,hits in zip(queries,stream,strict=True):
            if float(hits.Z) != Z:
                raise ValueError('Actual auxiliary search space mismatch')
            result.append({'sequence_id':q['sequence_id'],'protein_sha256':hashlib.sha256(q['sequence'].encode()).hexdigest(),
                'Z':float(hits.Z),'domZ':float(hits.domZ),'status':'scanned',
                'hits':[hit_record(h,i,old=old) for i,h in enumerate(hits)]})
        return result

    def apply_records(self, queries, original_records, scans, *, threads):
        # Complete each auxiliary scan before the next, so CPU pools do not overlap.
        old_batches = self._scan(queries,self.inherited_hmms,370,threads,old=True)
        new_batches = self._scan(queries,self.new_hmms,387,threads)
        extra_batches=[self._scan(queries,bank.hmms,bank.card['auxiliary_Z'],threads) for bank in self.extensions]
        manifest_digest = _sha(self.directory/'manifest.json')
        outputs=[]
        for index,(q,original,scan,old,new) in enumerate(zip(queries,original_records,scans,old_batches,new_batches,strict=True)):
            raw = scan['specificity_core_evidence']
            if raw['Z'] != 361:
                raise ValueError('Actual core search space mismatch')
            hits = raw['hits_in_production_order']
            winner = next((h for h in hits if h['profile']==original['best_positive_profile']),None)
            negatives = [h for h in hits if h['intended_use']=='hard_negative']
            best_negative = max(negatives,key=lambda h:h['score']) if negatives else None
            domains = winner['domains'] if winner else []
            domain = max(domains,key=lambda d:d['score']) if domains else None
            core = {**original,'protein_sha256':hashlib.sha256(q['sequence'].encode()).hexdigest(),'protein_length':len(q['sequence']),
                'positive_evalue':winner['evalue'] if winner else None,'positive_reported':winner['reported'] if winner else False,
                'positive_included':winner['included'] if winner else False,'hard_negative_hit_present':bool(negatives),
                'raw_best_negative_profile':best_negative['profile'] if best_negative else None,
                'raw_best_negative_score':best_negative['score'] if best_negative else None,
                'raw_best_negative_evalue':best_negative['evalue'] if best_negative else None,
                'hmm_coverage':(domain['hmm_to']-domain['hmm_from']+1)/domain['hmm_length'] if domain else None,
                'query_coverage':(domain['target_to']-domain['target_from']+1)/domain['target_length'] if domain else None,
                'diagnostic_hit_count':len(hits)}
            inherited=apply_competition(apply_prior_floor(core,self.floor),old['hits'],self.inherited_card)
            result=apply_cumulative(core,inherited,new,self.card,'intersection')
            # The canonical API uses 'no cas' for all rejected records. Reasons
            # live in explicit status/evidence, never an accepted display label.
            result['result']=result['cas_family'] if result['is_cas'] else 'no cas'
            raw_original = copy.deepcopy(original)
            raw_original['evidence']['model_id'] = self.manifest['core_original_model_id']
            result['evidence']['decision_threshold_scope'] = 'original_core_before_specificity_overlay'
            result['evidence']['specificity_overlay_id'] = self.manifest['id']
            result['specificity_evidence']={
                'schema_version':1,'overlay_id':self.manifest['id'],'mode':self.manifest['mode'],
                'manifest_sha256':manifest_digest,
                'core_original_model_id':self.manifest['core_original_model_id'],
                'core_original_prediction':raw_original,
                'core_search':raw,'inherited_search':old,'new_search':new,
                'confidence_probability':None,'live_deployment_approved':False}
            from casandra.specificity_extension import apply_extension
            for bank,batches in zip(self.extensions,extra_batches,strict=True):
                result=apply_extension(core,result,batches[index],bank.card)
                result['additional_specificity_evidence'][-1]['bank_manifest_sha256']=bank.manifest_sha256
            outputs.append(result)
        return outputs
