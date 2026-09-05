"""Build higher-rank records from exact family and slice certificates."""
from pathlib import Path
import json
import sys
import hashlib
import subprocess
import sympy as sp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'docs/catalogue'))
from build_catalogue import terms,source,z

def readlines(p):return [json.loads(x) for x in p.read_text().splitlines()]

def main(rank):
    directory=HERE/f'rank{rank}';base=f'research/higher_rank/rank{rank}'
    constants=json.loads((directory/'constant_candidates.json').read_text())
    witnesses={x['id']:x for x in readlines(directory/'lift_feasibility.jsonl')}
    families=readlines(directory/'families.jsonl')
    slices={x['id']:x for x in json.loads((directory/'slice_signatures.json').read_text())}
    queries={x['id']:x for x in readlines(directory/'verification.jsonl')}
    complete=(constants['enumeration_complete'] and len(witnesses)==len(queries)==constants['count']
              and all(x['status'] in ('sat','unsat') for x in witnesses.values())
              and all(x['result']=='unsat' for x in queries.values()))
    records=[]
    for number,family in enumerate(families,1):
        cid=family['id'];w=witnesses[cid];q=queries[cid]
        assert q['result']=='unsat' and family['coverage_status']=='unsat' and len(family['spaces'])==1
        space=family['spaces'][0];cert=space['certificate'];r=cert['delays'];n0=sp.diag(*(1+z**x for x in r))
        ap,am=(sp.Matrix([[sp.sympify(x,locals={'z':z}) for x in row] for row in cert[a]]) for a in ('A_plus','A_minus'))
        assert all(sp.expand(x)==0 for x in ap*am.subs(z,1/z).T-am*ap.subs(z,1/z).T)
        s=slices[cid]
        record={'schema_version':'2.0.0','id':f'r{rank}-c{number:02d}','rank':rank,'class_number':number,'constant_id':cid,
          'scope':{'symmetrizer':'identity','leading_permutation':'identity','indecomposable':True},'index_base':0,
          'datum':{'delays':r,'N_plus':terms(n0-ap),'N_minus':terms(n0-am),
                   'A_plus_1':[[int(x) for x in row] for row in ap.subs(z,1).tolist()],
                   'A_minus_1':[[int(x) for x in row] for row in am.subs(z,1).tolist()]},
          'family':{'parameters':['lambda']+[f's{i+1}' for i in range(rank-1)],'fixed_shift':{'species':rank-1,'value':'0'},
                    'dimension':space['dimension'],'rref':space['rref'],'variable_names':w['variable_names'],
                    'representative_values':space['values'],'coverage':q},
          'periodicity':{'time_coordinate':'displayed representative','h_plus':cert['positive']['h'],'h_minus':cert['negative']['h'],
                         'labelled_period':cert['labelled_tropical_seed_period'],
                         'positive_negative_permutation':cert['positive']['negative_permutation'],
                         'negative_negative_permutation':cert['negative']['negative_permutation']},
          'exchange':{k:cert[k] for k in ('vertices','B','mutation_vertices','relabel_old_to_new')},
          'slice':{k:s[k] for k in ('components','vertices','B','mutation_word','relabel_old_to_new','canonical_signature','signature_sha256','distinctness','signature_convention')},
          'provenance':{'verification_kind':'computer-assisted','source_commit':subprocess.check_output(['git','log','-1','--format=%H','--','research/higher_rank'],cwd=ROOT,text=True).strip(),
             'sources':[source(base+'/'+p,role) for p,role in [('constant_candidates.json','complete parity-compatible constant enumeration'),
                         ('families.jsonl','complete lift spaces and exact mutation certificates'),('verification.jsonl','independent replay'),
                         ('slice_signatures.json','slice distinctness')]],
             'manuscript':'research/higher_rank/methods.html','pdf':None,
             'query_path':base+'/smt_queries.zip','query_member':q['query_file']}}
        records.append(record)
    assert len(records)==sum(w['status']=='sat' for w in witnesses.values())
    (directory/'base-records.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8',newline='\n')
    report={'rank':rank,'complete':complete,'constant_candidates':constants['count'],'families':len(records),
            'unresolved_lifts':[cid for cid,w in witnesses.items() if w['status'] not in ('sat','unsat')],
            'unresolved_replays':[cid for cid,q in queries.items() if q['result']!='unsat'],
            'slice_distinctness':'all distinct by terminal cycle lengths and exhaustive orbits within repeated length classes',
            'scope':'Identity symmetrizer, diagonal leading matrix, indecomposable; rational time rescaling, species shifts, relabeling, sign exchange and slice changes.',
            'computation_seconds':sum(x['wall_seconds'] for x in readlines(directory/'computation-runs.jsonl'))}
    (directory/'classification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print(report)

if __name__=='__main__':main(int(sys.argv[1]))
