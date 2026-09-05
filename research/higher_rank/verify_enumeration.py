"""Check task coverage and reproduce the rank-four reference set exactly."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent


def task_data(directory):
    # A resumed run may have added tasks beyond the older checkpoint archive.
    paths=list((directory/'constant_tasks').glob('upper-*.json'))
    if paths:return [json.loads(p.read_text()) for p in paths]
    archive=directory/'constant_tasks.zip'
    if archive.exists():
        with zipfile.ZipFile(archive) as z:return [json.loads(z.read(name)) for name in z.namelist()]
    return []


def main():
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);args=parser.parse_args()
    control=[json.loads(p.read_text()) for p in (HERE/'control-rank4-ordered').glob('upper-*.json')]
    assert {x['upper_index'] for x in control}==set(range(120))
    actual={tuple(k) for x in control for k in x['keys']}
    old=json.loads((HERE.parent/'rank4/constant_candidates.json').read_text())['candidates'];expected=set()
    for c in old:
        a,b=([[(2 if i==j else 0)-c[field][i][j] for j in range(4)] for i in range(4)] for field in ('N_plus_1','N_minus_1'))
        if all(sum(x*y for x,y in zip(a[i],b[i]))%2==0 for i in range(4)):
            expected.add(tuple(x for field in ('N_plus_1','N_minus_1') for row in c[field] for x in row))
    assert actual==expected and len(actual)==460
    directory=HERE/f'rank{args.rank}';items=task_data(directory);count=math.factorial(args.rank+1)
    assert len(items)==count and {x['upper_index'] for x in items}==set(range(count))
    assert all(x['rank']==args.rank and x['completed'] and x['triangular_count']==count and x['parity_pruning'] for x in items)
    keys=sorted({tuple(k) for x in items for k in x['keys']})
    cs=json.loads((directory/'constant_candidates.json').read_text())
    assert cs['enumeration_complete'] and cs['count']==len(keys)
    assert keys==[tuple(x for field in ('N_plus_1','N_minus_1') for row in c[field] for x in row) for c in cs['candidates']]
    report={'rank':args.rank,'complete_upper_tasks':count,'candidate_count':len(keys),
            'rank_four_control':'exact equality of all 460 parity-compatible reference candidates',
            'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest(),
            'constant_file_sha256':hashlib.sha256((directory/'constant_candidates.json').read_bytes()).hexdigest()}
    pruned=sum(x.get('hereditary_pruning',False) for x in items)
    if pruned:
        controls=json.loads((HERE/'pruning-controls.json').read_text())
        assert controls['source_sha256']==report['source_sha256']
        report.update({'hereditary_pruning_tasks':pruned,'unpruned_tasks_retained':count-pruned,
                       'candidate_scope':'All retained necessary constants after certified hereditary-prefix exclusions; not the unpruned constant set.',
                       'pruning_controls':controls})
    (directory/'enumeration-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
