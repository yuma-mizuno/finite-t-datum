"""Collect complete constant tasks and compare the identity-symmetrizer control."""
import hashlib,itertools,json,math,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main(rank):
    directory=HERE/f'rank{rank}';tasks=[json.loads(p.read_text()) for p in (directory/'constant_tasks').glob('upper-*.json')]
    assert len(tasks)==math.factorial(rank+1) and {x['upper_index'] for x in tasks}==set(range(len(tasks)))
    assert all(x['completed'] and x['rank']==rank for x in tasks)
    def hashes(field):
        result=set()
        for task in tasks:
            value=task.get(field,[]);result.update(value if isinstance(value,list) else [value])
        return sorted(result)
    keys=sorted({tuple(k) for t in tasks for k in t['keys']});records=[]
    for i,k in enumerate(keys,1):
        d=list(k[:rank]);tail=k[rank:]
        records.append({'id':i,'symmetrizer':d,'N_plus_1':[list(tail[j*rank:(j+1)*rank]) for j in range(rank)],
                        'N_minus_1':[list(tail[rank*rank+j*rank:rank*rank+(j+1)*rank]) for j in range(rank)]})
    output={'rank':rank,'count':len(records),'enumeration_complete':True,'candidate_scope':'All primitive symmetrizable necessary constants, after any recorded hereditary exclusions.',
            'task_source_sha256':hashes('source_sha256'),
            'task_binary_sha256':hashes('binary_sha256'),
            'principal_table_sha256':hashes('principal_table_sha256'),
            'weighted_principal_table_sha256':hashes('weighted_principal_table_sha256'),
            'collector_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'candidates':records}
    (directory/'constant_candidates.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8',newline='\n')
    print(rank,len(records),'constant triples;',sum(r['symmetrizer']==[1]*rank for r in records),'with identity symmetrizer',flush=True)
    if rank in (3,4):
        old=json.loads((HERE.parent/f'rank{rank}/constant_candidates.json').read_text())['candidates']
        canonical=lambda p,m:min(tuple(b[i][j] for b in pair for i in perm for j in perm) for pair in ((p,m),(m,p)) for perm in itertools.permutations(range(rank)))
        expected={canonical(c['N_plus_1'],c['N_minus_1']) for c in old if all(sum((2*(i==j)-c['N_plus_1'][i][j])*(2*(i==j)-c['N_minus_1'][i][j]) for j in range(rank))%2==0 for i in range(rank))}
        actual={canonical(c['N_plus_1'],c['N_minus_1']) for c in records if c['symmetrizer']==[1]*rank}
        assert expected==actual,(len(expected),len(actual));print('Exact identity control:',len(actual),flush=True)
if __name__=='__main__':main(int(sys.argv[1]))
