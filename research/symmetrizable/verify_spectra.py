"""Independent interval logarithmic Jacobian, using the thesis dual fixed point."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from exponents import solution_box,RI,CI
from sage.all import matrix,QQ,ZZ,identity_matrix,diagonal_matrix
def mutate(b,k):
    n=b.nrows();return matrix(ZZ,[[(-b[i,j] if i==k or j==k else b[i,j]+max(b[i,k],0)*max(b[k,j],0)-max(-b[i,k],0)*max(-b[k,j],0)) for j in range(n)] for i in range(n)])
def main(rank):
    directory=HERE/f'rank{rank}';records=json.loads((directory/'base-records.json').read_text());spectra=json.loads((directory/'spectral-data.json').read_text());report={};start=time.monotonic()
    for r in records:
        info=spectra[r['id']];K=matrix(QQ,r['datum']['A_plus_1']).inverse()*matrix(QQ,r['datum']['A_minus_1']);D=diagonal_matrix(QQ,r['datum']['symmetrizer']);f,_=solution_box(D.inverse()*K*D)
        B=matrix(ZZ,r['exchange']['B']);n=B.nrows();J=identity_matrix(RI,n)
        for k in r['exchange']['mutation_vertices']:
            a=r['exchange']['vertices'][k][0];L=identity_matrix(RI,n);L[k,k]=-1
            for i in range(n):
                if i!=k:L[i,k]=max(B[k,i],0)*(1-f[a])+max(-B[k,i],0)*f[a]
            J=L*J;B=mutate(B,k)
        p=r['exchange']['relabel_old_to_new'];inv=[p.index(i) for i in range(n)]
        J=matrix(RI,[J.row(i) for i in inv]);B=matrix(ZZ,[[B[i,j] for j in inv] for i in inv]);assert B==matrix(ZZ,r['exchange']['B'])
        component={0}
        while True:
            more=component|{j for i in component for j in range(n) if B[i,j]}
            if more==component:break
            component=more
        indices=sorted(component);J=(J**r['slice']['components']).matrix_from_rows_and_columns(indices,indices)
        N=info['exponents']['root_order'];root=CI.zeta(N);power=identity_matrix(RI,len(indices));traces=[]
        for k in range(N):traces.append(power.trace());power=power*J
        values=[]
        for m in range(N):
            value=sum(CI(traces[k])*root**(-m*k) for k in range(N))/N
            lo=value.real().lower().ceil();hi=value.real().upper().floor();assert lo==hi and lo>=0 and value.imag().contains_zero()
            values.extend([m]*int(lo))
        assert values==info['exponents']['values'],r['id']
        report[r['id']]={'matrix_ratios':'exact','jacobian_spectrum':'independent interval-certified match','fixed_point_convention':'D^(-1) A_plus(1)^(-1) A_minus(1) D'}
        print(r['id'],'weighted spectrum independently verified',flush=True)
    (directory/'spectral-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps({'rank':rank,'stage':'independent-weighted-spectra','wall_seconds':time.monotonic()-start,'records':len(report)})+'\n')
if __name__=='__main__':main(int(sys.argv[1]))
