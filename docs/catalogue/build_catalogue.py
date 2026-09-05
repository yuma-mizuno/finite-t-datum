"""Build the offline catalogue and typed records from the checked research data."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
import subprocess

import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
z=sp.Symbol('z')


def read(path):
    p=ROOT/path
    return ([json.loads(line) for line in p.read_text(encoding='utf-8').splitlines()]
            if p.suffix=='.jsonl' else json.loads(p.read_text(encoding='utf-8')))


def source(path,role):
    return {'path':path,'role':role,'sha256':hashlib.sha256((ROOT/path).read_bytes()).hexdigest()}


def terms(a):
    return [[[[int(c),int(e[0])] for e,c in sorted(sp.Poly(a[i,j],z).terms(),key=lambda x:x[0]) if c]
             for j in range(a.cols)] for i in range(a.rows)]


def main():
    records=[]
    commit=subprocess.check_output(['git','log','-1','--format=%H','--','research/rank3','research/rank4'],cwd=ROOT,text=True).strip()
    for rank in (3,4):
        base=f'research/rank{rank}'
        witnesses={x['id']:x for x in read(base+'/lift_feasibility.jsonl') if x['status']=='sat'}
        constants={x['id']:x for x in read(base+'/constant_candidates.json')['candidates']}
        spaces={x['id']:x for x in read(base+('/lift_spaces.json' if rank==3 else '/families.jsonl'))}
        slices={x['id']:x for x in read(base+'/slice_signatures.json')}
        verified=read(base+'/verification.json')
        queries={x['id']:x for x in verified['queries']}
        assert all(x['result']=='unsat' for x in queries.values())
        assert len(witnesses)==(16 if rank==3 else 37)
        class_map=({x['constant_id']:x['class'] for x in read(base+'/classification.json')}
                   if rank==4 else {cid:k for k,cid in enumerate(sorted(witnesses),1)})
        for cid,w in sorted(witnesses.items()):
            space=spaces[cid]['spaces'][0]
            assert spaces[cid]['coverage_status']=='unsat' and len(spaces[cid]['spaces'])==1
            cert=w['certificate'] if rank==3 else space['certificate']
            r=cert['delays'];n0=sp.diag(*(1+z**x for x in r))
            ap,am=(sp.Matrix([[sp.sympify(x,locals={'z':z}) for x in row] for row in cert[field]])
                   for field in ('A_plus','A_minus'))
            np,nm=n0-ap,n0-am
            assert all(sp.expand(x)==0 for x in ap*am.subs(z,1/z).T-am*ap.subs(z,1/z).T)
            matrices={'plus':terms(np),'minus':terms(nm)}
            for sign in matrices.values():
                for i,row in enumerate(sign):
                    for entry in row:
                        assert all(c>0 and 0<p<r[i] for c,p in entry)
            for i in range(rank):
                for j in range(rank):
                    assert not ({p for c,p in matrices['plus'][i][j]} & {p for c,p in matrices['minus'][i][j]})
            assert np.subs(z,1).tolist()==constants[cid]['N_plus_1']
            assert nm.subs(z,1).tolist()==constants[cid]['N_minus_1']
            s=slices[cid]
            record={
                'schema_version':'1.0.0','id':f'r{rank}-c{class_map[cid]:02d}',
                'rank':rank,'class_number':class_map[cid],'constant_id':cid,
                'scope':{'symmetrizer':'identity','leading_permutation':'identity','indecomposable':True},
                'index_base':0,
                'datum':{'delays':r,'N_plus':matrices['plus'],'N_minus':matrices['minus'],
                         'A_plus_1':[[int(x) for x in row] for row in ap.subs(z,1).tolist()],
                         'A_minus_1':[[int(x) for x in row] for row in am.subs(z,1).tolist()]},
                'family':{'parameters':['lambda']+[f's{i+1}' for i in range(rank-1)],
                          'fixed_shift':{'species':rank-1,'value':'0'},
                          'dimension':space['dimension'],'rref':space['rref'],
                          'variable_names':w['variable_names'],
                          'representative_values':w['values'] if rank==3 else space['values'],
                          'coverage':queries[cid]},
                'periodicity':{'time_coordinate':'displayed representative',
                               'h_plus':cert['positive']['h'],'h_minus':cert['negative']['h'],
                               'labelled_period':cert['labelled_tropical_seed_period'],
                               'positive_negative_permutation':cert['positive']['negative_permutation'],
                               'negative_negative_permutation':cert['negative']['negative_permutation']},
                'exchange':{'vertices':cert['vertices'],'B':cert['B'],
                            'mutation_vertices':cert['mutation_vertices'],
                            'relabel_old_to_new':cert['relabel_old_to_new']},
                'slice':{'components':s['components'],'vertices':s['slice_vertices'],'B':s['slice_B'],
                         'mutation_word':s['slice_mutation_word'],
                         'relabel_old_to_new':s['slice_relabel_old_to_new'],
                         'canonical_signature':s['canonical_signature'],
                         'signature_sha256':s['signature_sha256']},
                'provenance':{'source_commit':commit,'verification_kind':'computer-assisted',
                              'sources':[source(base+'/lift_feasibility.jsonl','polynomial witnesses'),
                                         source(base+('/lift_spaces.json' if rank==3 else '/families.jsonl'),'complete lift families'),
                                         source(base+'/constant_candidates.json','independent constant matrices'),
                                         source(base+'/verification.json','exclusion and coverage replay'),
                                         source(base+'/slice_signatures.json','slice equivalence certificates')],
                              'manuscript':base+f'/rank{rank}-classification.tex',
                              'pdf':f'output/pdf/rank{rank}-classification.pdf',
                              'query_path':base+('/smt_queries/'+queries[cid]['file'] if rank==3 else '/smt_queries.zip'),
                              'query_member':queries[cid]['file']}}
            assert record['family']['dimension']==rank
            assert cert['positive'] and cert['negative'] and cert['labelled_tropical_seed_period']
            records.append(record)
    lower=read('research/catalogue/lower-rank-records.json')
    lower_commit=subprocess.check_output(['git','log','-1','--format=%H','--','research/catalogue'],cwd=ROOT,text=True).strip() or commit
    for record in lower:
        record['provenance']['source_commit']=lower_commit
        record['provenance']['sources']=[source('research/catalogue/lower_ranks.py','exact lower-rank representatives and mutation certificates'),
                                       source('research/catalogue/lower-ranks-proof.html','rank-one proof and rank-two literature reference')]
    records.extend(lower)
    spectra=read('research/catalogue/spectral-data.json')
    family_notes=read('research/catalogue/family-notes.json')
    quivers=read('research/catalogue/quiver-data.json')
    verification=read('research/catalogue/enrichment-verification.json')
    assert len(records)==len(spectra)==len(family_notes)==len(quivers)==len(verification)==61
    for record in records:
        record['schema_version']='2.0.0'
        record['matrix_ratios']=spectra[record['id']]['matrix_ratios']
        record['exponents']=spectra[record['id']]['exponents']
        plot=f'research/catalogue/plots/{record["id"]}.svg'
        assert (ROOT/plot).is_file()
        record['exponents']['plot_path']=plot
        record['verification']=verification[record['id']]
        record['notes']={'family':family_notes[record['id']], 'quiver':quivers[record['id']]}
        assert record['notes']['quiver']['status'].startswith('certified-')
        record['provenance']['enrichment_commit']=lower_commit
        record['provenance']['enrichment_sources']=[source('research/catalogue/'+p,role) for p,role in
            [('spectral-data.json','certified exponent multiplicities and exact matrix ratios'),('family-notes.json','named family matches'),
             ('quiver-data.json','quiver mutation-class certificates'),('enrichment-verification.json','independent exact-path replay and interval Jacobian verification'),
             ('sources/mizuno-thesis.pdf','Mizuno thesis, equations (3.2.1), (3.4.1) and Definition 3.4.4')]]
    higher=[]
    for path in sorted((ROOT/'research/higher_rank').glob('rank*/catalogue-records.json')):
        classification=json.loads((path.parent/'classification.json').read_text())
        assert classification['complete'], path
        extra=json.loads(path.read_text(encoding='utf8'))
        assert len(extra)==classification['families']
        records.extend(extra);higher.append(classification)
    records.sort(key=lambda x:(x['rank'],x['class_number']))
    dataset={'schema_version':'2.0.0','title':'Finite T-data catalogue',
             'polynomial_encoding':'Each entry is a list of [coefficient, exponent] pairs, with ascending exponents.',
             'equivalence':['admissible rational time rescaling','species shifts','simultaneous index permutations','sign exchange'],
             'scope':f'Identity symmetrizer and diagonal N0; indecomposable data of ranks 1 through {max(r["rank"] for r in records)}.',
             'higher_rank_classifications':higher,
             'source_commit':commit,'records':records,'proofs':json.loads((HERE/'proofs.json').read_text(encoding='utf-8'))}
    (HERE/'records').mkdir(exist_ok=True)
    for record in records:
        (HERE/'records'/f'{record["id"]}.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    (HERE/'catalogue.json').write_text(json.dumps(dataset,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    html=(HERE/'index.template.html').read_text(encoding='utf-8')
    buttons=''.join(f'<button data-rank="{rank}" aria-pressed="{str(rank==4).lower()}">Rank {rank} <span>{sum(r["rank"]==rank for r in records)}</span></button>' for rank in sorted({r['rank'] for r in records},reverse=True))
    html=re.sub(r'(<div class="rank-picker"[^>]*>).*?(</div>)',lambda m:m[1]+buttons+'<button data-rank="all" aria-pressed="false">All ranks</button>'+m[2],html,count=1)
    for name,path in [('CSS','catalogue.css'),('CORE','core.js'),('APP','app.js')]:
        html=html.replace('/*__'+name+'__*/',(HERE/path).read_text(encoding='utf-8'))
    packed=json.dumps(dataset,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    html=html.replace('/*__DATA__*/',packed)
    (HERE/'index.html').write_text(html,encoding='utf-8',newline='\n')
    print(f'Built {len(records)} checked records and an offline document ({len(html.encode()):,} bytes).')


if __name__=='__main__':main()
