"""Check mixed positivity independently and compare complete rank-four/five unions."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib,json,random,sys
from principal import obstruction
from verify_subset_controls import keys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'rank3/_deps'));import z3
import sympy as S
def mixture_positive(p,m):
    n=len(p)
    for signs in product(range(2),repeat=n):
        a=S.Matrix([[2*(i==j)-(m if signs[j] else p)[i][j] for j in range(n)] for i in range(n)])
        if any(a[:k,:k].det()<=0 for k in range(1,n+1)):return False
    return True
def common_positive(p,m):
    n=len(p);v=z3.Reals(' '.join(f'v{i}' for i in range(n)));solver=z3.Solver()
    solver.add(*(x>=1 for x in v));solver.add(*(z3.Sum([(2*(i==j)-b[i][j])*v[i] for i in range(n)])>=1 for b in (p,m) for j in range(n)))
    return solver.check()==z3.sat
def main():
    count=0
    for diag in product(range(2),repeat=4):
        for edges in product(range(4),repeat=4):
            p=[[diag[0],edges[0]],[edges[1],diag[1]]];m=[[diag[2],edges[2]],[edges[3],diag[3]]]
            lo=max(Fraction(b[0][1],2-b[1][1]) for b in (p,m));upper=[Fraction(2-b[0][0],b[1][0]) for b in (p,m) if b[1][0]]
            expected=not upper or lo<min(upper)
            actual=all((2-(m if a else p)[0][0])*(2-(m if b else p)[1][1])-(m if b else p)[0][1]*(m if a else p)[1][0]>0 for a,b in product(range(2),repeat=2))
            assert actual==expected;count+=1
    rng=random.Random(61093);random_count=0
    for n in range(3,7):
        for sample in range(20):
            p,m=[[[rng.randrange(2) if i==j else (rng.randrange(1,4) if rng.random()<.22 else 0) for j in range(n)] for i in range(n)] for _ in range(2)]
            assert mixture_positive(p,m)==common_positive(p,m),(n,p,m);random_count+=1
    comparisons=[]
    for n in (4,5):
        candidates=json.loads((HERE/f'rank{n}/constant_candidates.json').read_text())['candidates']
        expected={tuple(c['symmetrizer'])+tuple(x for a in ('N_plus_1','N_minus_1') for row in c[a] for x in row) for c in candidates if obstruction(c) is None}
        actual=keys(HERE/f'control-rank{n}-mmatrix')
        assert actual==expected,(n,len(expected),len(actual),len(expected-actual),len(actual-expected))
        comparisons.append({'rank':n,'original_candidates':len(candidates),'all_subset_survivors':len(expected),'complete_union_equal':True})
    result={'rank_two_all_4096_pairs':count,'independent_random_linear_solver_controls':random_count,'complete_reference_comparisons':comparisons,
            'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest(),
            'weighted_table_sha256':hashlib.sha256((HERE/'principal-constants-rank6.txt.weights').read_bytes()).hexdigest()}
    (HERE/'mmatrix-controls.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(result,flush=True)
if __name__=='__main__':main()
