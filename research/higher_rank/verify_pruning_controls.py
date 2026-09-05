"""Differential checks against completed ranks; not a substitute for the lemmas."""
import hashlib
import json
from pathlib import Path
import re
import zipfile
from obstructions import key

HERE=Path(__file__).resolve().parent


def read_tasks(name):
    tasks=[json.loads(p.read_text()) for p in (HERE/name).glob('upper-*.json')]
    if not tasks:
        with zipfile.ZipFile(HERE/(name+'.zip')) as z:tasks=[json.loads(z.read(p)) for p in z.namelist()]
    return tasks,{tuple(k) for x in tasks for k in x['keys']}


def main():
    records=json.loads((HERE.parents[1]/'docs/catalogue/catalogue.json').read_text(encoding='utf8'))['records']
    reports=[]
    for n in (4,5):
        tasks,columns=read_tasks(f'control-rank{n}-columns')
        pruned_tasks,pruned=read_tasks(f'control-rank{n}-hereditary')
        if n==4:_,expected=read_tasks('control-rank4-ordered')
        else:
            cs=json.loads((HERE/'rank5/constant_candidates.json').read_text())['candidates']
            expected={tuple(x for field in ('N_plus_1','N_minus_1') for row in c[field] for x in row) for c in cs}
        known={key(*[[[2*(i==j)-r['datum'][a][i][j] for j in range(n)] for i in range(n)] for a in ('A_plus_1','A_minus_1')]) for r in records if r['rank']==n}
        assert columns==expected and known<=pruned and pruned<=expected
        reports.append({'rank':n,'column_bounds_exact_reference_equality':len(expected),
                        'hereditary_retained_candidates':len(pruned),'all_classified_families_preserved':len(known),
                        'completed_tasks':len(tasks),'hereditary_completed_tasks':len(pruned_tasks)})
    report={'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest(),
            'principal_table_sha256':hashlib.sha256((HERE/'principal-constants.txt').read_bytes()).hexdigest(),
            'controls':reports}
    (HERE/'pruning-controls.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print(report)


if __name__=='__main__':main()
