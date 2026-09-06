"""Independent rank-two brute force and exact rank-four algorithm controls."""
from fractions import Fraction
from itertools import product,permutations
from math import gcd
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def task_keys(p):return {tuple(k) for f in p.glob('upper-*.json') for k in json.loads(f.read_text())['keys']}
def canonical(d,p,m):
    return min(tuple(d[i] for i in perm)+tuple(b[i][j] for b in pair for i in perm for j in perm) for perm in permutations(range(2)) for pair in ((p,m),(m,p)))
def brute_rank_two():
    result=set();tried=0
    for d in product(range(1,4),repeat=2):
        if gcd(*d)!=1:continue
        for diag in product(range(2),repeat=4):
            for edges in product(range(4),repeat=4):
                p=[[diag[0],edges[0]],[edges[1],diag[1]]];m=[[diag[2],edges[2]],[edges[3],diag[3]]];tried+=1
                if not(max(p[0][1],m[0][1]) and max(p[1][0],m[1][0])):continue
                if any(b[i][j]*d[j]%d[i] or b[i][j]*d[j]//d[i]>3 for b in (p,m) for i in range(2) for j in range(2)):continue
                lower=max(Fraction(b[0][1],2-b[1][1]) for b in (p,m))
                uppers=[Fraction(2-b[0][0],b[1][0]) for b in (p,m) if b[1][0]]
                if uppers and lower>=min(uppers):continue
                a,b=([[(2 if i==j else 0)-c[i][j] for j in range(2)] for i in range(2)] for c in (p,m))
                g=[[sum(a[i][k]*d[k]*b[j][k] for k in range(2)) for j in range(2)] for i in range(2)]
                if g[0][1]!=g[1][0] or any(g[i][i]%(2*d[i]) for i in range(2)):continue
                result.add(canonical(d,p,m))
    assert result==task_keys(HERE/'rank2/constant_tasks') and len(result)==11
    return {'integer_triples_tested':tried,'primitive_diagonal_bound':3,'exact_constant_triples':11}
def main():
    r2=brute_rank_two();expected=task_keys(HERE/'rank4/constant_tasks');assert len(expected)==1040
    controls=[]
    for name in ['control-rank4-ratio','control-rank4-ratio2']:
        assert task_keys(HERE/name)==expected
        controls.append({'directory':name,'exact_equal_constant_triples':1040})
    out={'rank_two_independent_brute_force':r2,'rank_four_exact_controls':controls,
         'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest()}
    (HERE/'constant-controls.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
