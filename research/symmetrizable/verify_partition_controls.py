"""Independent signed-mass partition controls for the diagonal Laurent reduction."""
import hashlib,itertools,json,math
from pathlib import Path
from principal import obstruction
from verify_subset_controls import keys
HERE=Path(__file__).resolve().parent

def allowed(key,n):
    d=key[:n];flat=key[n:];p,m=[[[flat[s*n*n+i*n+j] for j in range(n)] for i in range(n)] for s in range(2)]
    for i in range(n):
        terms=[(p[i][k]*m[i][k],d[k]) for k in range(n) if p[i][k]*m[i][k]]
        if not any(c%2 for c,w in terms):continue
        if len(terms)==1:return False
        if math.prod(c+1 for c,w in terms)>512:continue
        possible={0}
        for c,w in terms:possible={a+w*t for a in possible for t in range(-c,c+1,2)}
        if 0 not in possible:return False
    return True

def main():
    results=[]
    for n in (4,5):
        cs=json.loads((HERE/f'rank{n}/constant_candidates.json').read_text())['candidates']
        unfiltered={tuple(c['symmetrizer'])+tuple(x for s in ('N_plus_1','N_minus_1') for row in c[s] for x in row) for c in cs if obstruction(c) is None}
        if n==5:assert keys(HERE/'control-rank5-signs')==unfiltered
        expected={k for k in unfiltered if allowed(k,n)}
        directory=HERE/f'control-rank{n}-partition';assert len(list(directory.glob('upper-*.json')))==math.factorial(n+1)
        actual=keys(directory);assert expected==actual,(n,len(expected),len(actual),len(expected-actual),len(actual-expected))
        results.append({'rank':n,'full_upper_tasks':math.factorial(n+1),'filtered_constants':len(actual),'equal_to_independently_filtered_reference':True})
    records=json.loads((HERE.parents[1]/'docs/catalogue/catalogue.json').read_text())['records']
    for r in records:
        n=r['rank'];d=r['datum'].get('symmetrizer',[1]*n)
        key=tuple(d)+tuple(2*(i==j)-r['datum'][s][i][j] for s in ('A_plus_1','A_minus_1') for i in range(n) for j in range(n))
        assert allowed(key,n),r['id']
        for i in range(n):
            mass=0
            for k in range(n):
                for a,p in r['datum']['N_plus'][i][k]:
                    for b,q in r['datum']['N_minus'][i][k]:
                        assert p!=q
                        mass+=d[k]*a*b*(1 if p>q else -1)
            assert mass==0,(r['id'],i,mass)
    result={'comparisons':results,'all_verified_families_retained':len(records),'signed_positive_negative_mass_checked_on_polynomials':True,
            'source_sha256':hashlib.sha256((HERE/'benchmark_partition_strategy.cpp').read_bytes()).hexdigest()}
    (HERE/'partition-controls.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(result,flush=True)
if __name__=='__main__':main()
