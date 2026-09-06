"""Independent full-subset filtering of reference constants and adaptive controls."""
import hashlib,itertools,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def table(path):
    out={}
    for line in path.read_text().splitlines():
        n,*key=map(int,line.split());out.setdefault(n,set()).add(tuple(key))
    return out
def allowed(key,n,lookup):
    flat=key[n:];p,m=[[[flat[s*n*n+i*n+j] for j in range(n)] for i in range(n)] for s in range(2)]
    for size in range(1,n):
        for S in itertools.combinations(range(n),size):
            if any(any(p[i][k] for i in S) and any(m[i][k] for i in S) for k in range(n) if k not in S):continue
            remaining=set(S)
            while remaining:
                component={min(remaining)}
                while True:
                    more=component|{j for i in component for j in remaining if p[i][j] or p[j][i] or m[i][j] or m[j][i]}
                    if more==component:break
                    component=more
                remaining-=component;T=sorted(component);value=tuple(b[i][j] for b in (p,m) for i in T for j in T)
                if value not in lookup[len(T)]:return False
    return True
def keys(directory):return {tuple(k) for p in directory.glob('upper-*.json') for k in json.loads(p.read_text())['keys']}
def main():
    lookup=table(HERE/'principal-constants.txt');reference=keys(HERE/'rank4/constant_tasks');expected={k for k in reference if allowed(k,4,lookup)}
    assert len(expected)==513
    for name in ('control-rank4-subsets','control-rank4-partial','control-rank4-adaptive'):
        assert len(list((HERE/name).glob('upper-*.json')))==120 and keys(HERE/name)==expected,name
    reference={tuple(k) for k in json.loads((HERE/'rank5/constant_tasks/upper-122.json').read_text())['keys']}
    expected5={k for k in reference if allowed(k,5,lookup)};actual=keys(HERE/'control-rank5-adaptive')
    assert actual==expected5 and len(actual)==283,(len(actual),len(expected5))
    report={'rank4_reference_constants':len(reference),'rank4_all_subset_survivors':513,'rank4_three_algorithms_equal':True,
            'rank5_upper_122_survivors':283,'rank5_exact_filtered_reference_match':True,
            'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest(),
            'table_sha256':hashlib.sha256((HERE/'principal-constants.txt').read_bytes()).hexdigest()}
    report['rank4_reference_constants']=1040
    (HERE/'subset-controls.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n');print(report)
if __name__=='__main__':main()
