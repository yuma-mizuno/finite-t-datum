"""Canonical connected-slice loops of the certified higher-rank families."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from collections import defaultdict

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from lower_ranks import slice_data
sys.path.insert(0,str(HERE.parent/'rank4'))
from rank4_slices import equivalence_states


def main():
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);parser.add_argument('--engine',choices=['python','cpp'],default='python');args=parser.parse_args()
    directory=HERE/f'rank{args.rank}';source=directory/'families.jsonl';start=time.monotonic()
    families=[json.loads(x) for x in source.read_text().splitlines()]
    results=[];minima={};groups=defaultdict(list);prepared=[]
    for f in families:
        assert f['coverage_status']=='unsat'
        for index,space in enumerate(f['spaces']):
            cert=space['certificate']
            assert space['scaling_and_shifts'] and cert['positive'] and cert['negative'] and cert['labelled_tropical_seed_period']
            s=slice_data(cert);b=tuple(map(tuple,s['B']));word=tuple(s['mutation_word']);p=tuple(s['relabel_old_to_new'])
            assert len(word)==args.rank
            seen=set();lengths=[]
            for i in range(len(p)):
                if i in seen:continue
                k=i;length=0
                while k not in seen:seen.add(k);length+=1;k=p[k]
                lengths.append(length)
            lengths=tuple(sorted(lengths));groups[lengths].append([f['id'],index])
            prepared.append((f['id'],index,s,b,word,p,lengths))
    checkpoints=directory/'slice_tasks';checkpoints.mkdir(exist_ok=True)
    source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    for cid,index,s,b,word,p,lengths in prepared:
        target=checkpoints/f'{cid}-{index}.json'
        if target.exists():
            saved=json.loads(target.read_text());assert saved['family_source_sha256']==source_hash
            s=saved['slice'];signature=json.dumps(s['canonical_signature'])
        else:
            collision=len(groups[lengths])>1
            if collision and args.engine=='cpp':
                from slice_cpp import orbit
                result=orbit(b,word,p,directory/'slice_cpp',f'{cid}-{index}')
                minimum=result['minimum'];state_count=result['states']
            else:
                states=(equivalence_states(b,word,p)|equivalence_states(tuple(tuple(-x for x in row) for row in b),word,p)) if collision else set()
                minimum=min(states) if collision else None;state_count=len(states) if collision else None
            encoding={'terminal_cycle_lengths':lengths,'orbit_minimum':minimum}
            signature=json.dumps(encoding)
            s.update({'id':cid,'space_index':index,'equivalent_canonical_states_including_sign_exchange':state_count,
                      'canonical_signature':encoding,'signature_sha256':hashlib.sha256(signature.encode()).hexdigest(),
                      'signature_convention':'Catalogue distinctness: terminal cycle lengths; exhaustive cyclic-loop minimum within repeated length classes.',
                      'distinctness':{'terminal_cycle_lengths':lengths,'records_with_same_lengths':groups[lengths],
                                      'method':'exhaustive cyclic-loop orbit' if collision else 'unique invariant terminal cycle lengths'}})
            target.write_text(json.dumps({'family_source_sha256':source_hash,'slice':s},indent=2)+'\n',encoding='utf8',newline='\n')
        assert signature not in minima,(cid,index,minima.get(signature));minima[signature]=[cid,index]
        results.append(s);print(cid,'vertices',s['vertices'],s['distinctness']['method'],flush=True)
    (directory/'slice_signatures.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
    report={'rank':args.rank,'stage':'slice-distinctness','wall_seconds':time.monotonic()-start,'families':len(results),
            'status':'all distinct','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps(report)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
