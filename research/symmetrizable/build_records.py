"""Prepare new nonidentity records from completed weighted lift certificates."""
import hashlib,json,math,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/catalogue'));from lower_ranks import slice_data
sys.path.insert(0,str(ROOT/'docs/catalogue'));from build_catalogue import terms,source,z,sp,archive_path
def readlines(path):return {r['id']:r for r in map(json.loads,path.read_text().splitlines())}
def main(rank):
    directory=HERE/f'rank{rank}';base=f'research/symmetrizable/rank{rank}'
    constants=json.loads((directory/'constant_candidates.json').read_text());cs={c['id']:c for c in constants['candidates']}
    ws=readlines(directory/'lift_feasibility.jsonl');fs=readlines(directory/'families.jsonl');vs=readlines(directory/'verification.jsonl')
    assert len(ws)==len(vs)==constants['count'] and all(w['status'] in ('sat','unsat') for w in ws.values())
    assert all(v['result']=='unsat' for v in vs.values())
    records=[]
    for cid,f in sorted(fs.items()):
        d=cs[cid]['symmetrizer']
        if d==[1]*rank:continue
        assert f['coverage_status']=='unsat' and len(f['spaces'])==1
        space=f['spaces'][0];cert=space['certificate'];w=ws[cid];q=vs[cid]
        assert space['scaling_and_shifts'] and space['dimension']==rank
        assert cert['positive'] and cert['negative'] and cert['labelled_tropical_seed_period']
        delays=cert['delays'];n0=sp.diag(*(1+z**r for r in delays));ap,am=(sp.Matrix([[sp.sympify(x,locals={'z':z}) for x in row] for row in cert[key]]) for key in ('A_plus','A_minus'))
        s=slice_data(cert);component={0}
        while True:
            more=component|{j for i in component for j,x in enumerate(cert['B'][i]) if x}
            if more==component:break
            component=more
        s['symmetrizer']=[d[cert['vertices'][i][0]] for i in sorted(component)]
        number=len(records)+1
        r={'schema_version':'2.1.0','id':f's{rank}-c{number:02d}','rank':rank,'class_number':number,'constant_id':cid,
           'scope':{'symmetrizer':'positive_diagonal','leading_permutation':'identity','indecomposable':True},'index_base':0,
           'datum':{'delays':delays,'symmetrizer':d,'N_plus':terms(n0-ap),'N_minus':terms(n0-am),
                    'A_plus_1':[[int(x) for x in row] for row in ap.subs(z,1).tolist()],
                    'A_minus_1':[[int(x) for x in row] for row in am.subs(z,1).tolist()]},
           'family':{'parameters':['lambda']+[f's{i+1}' for i in range(rank-1)],'fixed_shift':{'species':rank-1,'value':'0'},
                     'dimension':rank,'rref':space['rref'],'variable_names':w['variable_names'],'representative_values':space['values'],
                     'coefficient_units':[[d[i]//math.gcd(d[i],d[j]) for j in range(rank)] for i in range(rank)],'coverage':q},
           'periodicity':{'time_coordinate':'displayed representative','h_plus':cert['positive']['h'],'h_minus':cert['negative']['h'],
                          'labelled_period':cert['labelled_tropical_seed_period'],'positive_negative_permutation':cert['positive']['negative_permutation'],
                          'negative_negative_permutation':cert['negative']['negative_permutation']},
           'exchange':{**{k:cert[k] for k in ('vertices','B','mutation_vertices','relabel_old_to_new')},'symmetrizer':[d[i] for i,p in cert['vertices']]},
           'slice':s,
           'provenance':{'verification_kind':'computer-assisted','source_commit':subprocess.check_output(['git','log','-1','--format=%H','--','research/symmetrizable'],cwd=ROOT,text=True).strip(),
                         'sources':[source(base+'/'+p,role) for p,role in [('constant_candidates.json','complete primitive weighted constants'),('families.jsonl','complete polynomial lift spaces and periodicity'),('verification.jsonl','exact exclusion and coverage replay')]],
                         'manuscript':'research/symmetrizable/methods.html','pdf':None,'query_path':archive_path(directory,q['query_file']).relative_to(ROOT).as_posix(),'query_member':q['query_file']}}
        records.append(r)
    (directory/'base-records.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8',newline='\n')
    report={'rank':rank,'constant_candidates':constants['count'],'polynomial_families':len(fs),'identity_families':len(fs)-len(records),'new_nonidentity_families':len(records),
            'lift_classification_complete':True,'slice_distinctness_complete':False,'enrichment_complete':False,'complete':False}
    (directory/'classification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n');print(report,flush=True)
if __name__=='__main__':main(int(sys.argv[1]))
