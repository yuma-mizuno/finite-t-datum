"""Certified exponent data and bounded quiver recognition for higher ranks."""
import json
from pathlib import Path
import sys
import time
import signal
from sage.all import matrix,ZZ,QQ,identity_matrix,ClusterQuiver,set_random_seed

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from exponents import compute
from quiver_types import mutate,dynkin_path

def quiver(record):
    b=record['slice']['B'];n=len(b)
    result={'vertices':n,'scope':'One connected component; the mutation class alone does not specify its loop.',
            'status':'unidentified','label':'Mutation class not identified','mutation_finite':None,'cluster_finite':None,
            'reference':'https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/cluster_algebra_quiver/quiver.html'}
    try:
        signal.alarm(4)
        q=ClusterQuiver(matrix(ZZ,b));set_random_seed(91500+record['class_number'])
        raw=q.mutation_type();result['recognition']=str(raw)
        if hasattr(raw,'is_finite') and raw.is_finite():
            result['label']=str(raw);result['status']='recognized';result['cluster_finite']=True;result['mutation_finite']=True
            cert=dynkin_path(b,limit=8000)
            if cert:
                kind,path,target=cert
                result.update({'status':'certified-dynkin','label':kind,'certificate':{'mutation_path':path,'target_B':target,'target_graph':kind}})
        else:
            outcome=q.is_mutation_finite(nr_of_checks=200*n,return_path=True)
            if isinstance(outcome,tuple) and not outcome[0]:
                path=[int(k) for k in outcome[1]];target=tuple(map(tuple,b))
                for k in path:target=mutate(target,k)
                edge=next((i,j,abs(target[i][j])) for i in range(n) for j in range(i+1,n) if abs(target[i][j])>2)
                result.update({'status':'certified-mutation-infinite','label':'Mutation-infinite','mutation_finite':False,'cluster_finite':False,
                    'certificate':{'mutation_path':path,'target_B':[list(row) for row in target],'edge_with_at_least_three_arrows':list(edge)}})
    except TimeoutError:result['search_note']='Recognition stopped after four seconds; no unproved class label is asserted.'
    finally:signal.alarm(0)
    return result

def main():
    rank=int(sys.argv[1]);mode=sys.argv[2];directory=HERE/f'rank{rank}'
    records=json.loads((directory/'base-records.json').read_text());target=directory/(mode+'-data.json')
    results=json.loads(target.read_text()) if target.exists() else {}
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError()))
    for record in records:
        if record['id'] in results:continue
        results[record['id']]=compute(record) if mode=='spectral' else quiver(record)
        if mode=='spectral':
            a,b=(matrix(QQ,record['datum'][key]) for key in ('A_plus_1','A_minus_1'))
            ratios=results[record['id']]['matrix_ratios'];K=matrix(QQ,ratios['A_plus_inverse_A_minus']);Ki=matrix(QQ,ratios['A_minus_inverse_A_plus'])
            assert a*K==b and b*Ki==a and K*Ki==identity_matrix(QQ,rank)
        target.write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
        print(record['id'],mode,results[record['id']].get('label','certified'),flush=True)

if __name__=='__main__':main()
