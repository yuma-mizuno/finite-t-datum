"""Independent logarithmic mutation-Jacobian verification of higher-rank spectra."""
import json
from pathlib import Path
import sys
from sage.all import matrix,vector,identity_matrix,QQ,ZZ
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from exponents import solution_box,RI,CI
from verify_enrichment import mutation

def main(rank):
    directory=HERE/f'rank{rank}';rs=json.loads((directory/'base-records.json').read_text());spectra=json.loads((directory/'spectral-data.json').read_text())
    target=directory/'spectral-verification.json';out=json.loads(target.read_text()) if target.exists() else {}
    for r in rs:
        if r['id'] in out:continue
        info=spectra[r['id']];K=matrix(QQ,info['matrix_ratios']['A_plus_inverse_A_minus']);f,_=solution_box(K)
        B=matrix(ZZ,r['exchange']['B']);n=B.nrows();J=identity_matrix(RI,n)
        for k in r['exchange']['mutation_vertices']:
            species=r['exchange']['vertices'][k][0];L=identity_matrix(RI,n);L[k,k]=-1
            for i in range(n):
                if i!=k:L[i,k]=max(B[k,i],0)*(1-f[species])+max(-B[k,i],0)*f[species]
            J=L*J;B=mutation(B,k)
        perm=r['exchange']['relabel_old_to_new'];inv=[perm.index(i) for i in range(n)]
        J=matrix(RI,[J.row(i) for i in inv]);B=matrix(ZZ,[[B[i,j] for j in inv] for i in inv]);assert B==matrix(ZZ,r['exchange']['B'])
        component={0}
        while True:
            more=component|{j for i in component for j in range(n) if B[i,j]}
            if more==component:break
            component=more
        indices=sorted(component);J=(J**r['slice']['components']).matrix_from_rows_and_columns(indices,indices)
        N=info['exponents']['root_order'];root=CI.zeta(N);power=identity_matrix(RI,len(indices));traces=[]
        for k in range(N):traces.append(power.trace());power=power*J
        actual=[]
        for m in range(N):
            v=sum(CI(traces[k])*root**(-m*k) for k in range(N))/N
            lo=v.real().lower().ceil();hi=v.real().upper().floor();assert lo==hi and lo>=0 and v.imag().contains_zero()
            actual.extend([m]*int(lo))
        assert actual==info['exponents']['values'],r['id']
        out[r['id']]={'matrix_ratios':'exact','jacobian_spectrum':'independent interval-certified match'}
        target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n');print(r['id'],'Jacobian verified',flush=True)
if __name__=='__main__':main(int(sys.argv[1]))
