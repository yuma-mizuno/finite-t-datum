"""Independently filter exact reference unions by diagonal Laurent reciprocity."""
import hashlib,json,math
from pathlib import Path
from principal import obstruction
from verify_subset_controls import keys
HERE=Path(__file__).resolve().parent

def allowed(key,n):
    d=key[:n];flat=key[n:];p,m=[[[flat[s*n*n+i*n+j] for j in range(n)] for i in range(n)] for s in range(2)]
    for i in range(n):
        if p[i][i] or m[i][i]:continue
        common=[k for k in range(n) if p[i][k]*m[i][k]]
        if len(common)==1 and p[i][common[0]]*m[i][common[0]]%2:return False
        if len(common)==2 and all(p[i][k]*m[i][k]==1 for k in common) and d[common[0]]!=d[common[1]]:return False
    return True

def main():
    results=[]
    for n in (4,5):
        cs=json.loads((HERE/f'rank{n}/constant_candidates.json').read_text())['candidates']
        expected={tuple(c['symmetrizer'])+tuple(x for s in ('N_plus_1','N_minus_1') for row in c[s] for x in row) for c in cs if obstruction(c) is None}
        expected={k for k in expected if allowed(k,n)}
        directory=HERE/f'control-rank{n}-overlap';assert len(list(directory.glob('upper-*.json')))==math.factorial(n+1)
        actual=keys(directory);assert expected==actual,(n,len(expected),len(actual),len(expected-actual),len(actual-expected))
        results.append({'rank':n,'full_upper_tasks':math.factorial(n+1),'filtered_constants':len(actual),'equal_to_independently_filtered_reference':True})
    expected={k for k in keys(HERE/'control-rank6-mmatrix') if allowed(k,6)}
    actual={tuple(k) for k in json.loads((HERE/'control-rank6-overlap/upper-3.json').read_text())['keys']}
    assert actual==expected,(len(actual),len(expected),len(actual-expected),len(expected-actual))
    records=json.loads((HERE.parents[1]/'docs/catalogue/catalogue.json').read_text())['records']
    for r in records:
        n=r['rank'];key=tuple(r['datum'].get('symmetrizer',[1]*n))+tuple(2*(i==j)-r['datum'][s][i][j] for s in ('A_plus_1','A_minus_1') for i in range(n) for j in range(n))
        assert allowed(key,n),r['id']
    result={'comparisons':results,'rank6_upper3_exact_filtered_comparison':len(actual),'all_verified_families_retained':len(records),
            'source_sha256':hashlib.sha256((HERE/'benchmark_overlap_strategy.cpp').read_bytes()).hexdigest()}
    (HERE/'overlap-controls.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(result,flush=True)
if __name__=='__main__':main()
