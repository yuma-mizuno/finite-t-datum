"""Principal-block exclusions and matching equations, with symmetrizers."""
import itertools,json,re
from functools import lru_cache
from math import gcd
from pathlib import Path
HERE=Path(__file__).resolve().parent

def readlines(p):return {r['id']:r for r in map(json.loads,p.read_text().splitlines())} if p.exists() else {}

def label(d,p,m):
    factor=gcd(*d)
    return tuple(x//factor for x in d)+tuple(x for b in (p,m) for row in b for x in row)

@lru_cache(maxsize=1)
def tables():
    allowed={1:{(1,0,0),(1,0,1),(1,1,0)}};spaces={}
    for n in range(2,6):
        directory=HERE/f'rank{n}';path=directory/'constant_candidates.json'
        if not path.exists():continue
        data=json.loads(path.read_text());assert data['enumeration_complete']
        ws=readlines(directory/'lift_feasibility.jsonl');fs=readlines(directory/'families.jsonl')
        allowed[n]=set();spaces[n]={}
        for c in data['candidates']:
            if ws.get(c['id'],{}).get('status')=='unsat':continue
            p,m=c['N_plus_1'],c['N_minus_1'];d=c['symmetrizer'];family=fs.get(c['id'])
            for perm in itertools.permutations(range(n)):
                for swap in range(2):
                    pair=(p,m) if not swap else (m,p)
                    key=label([d[i] for i in perm],*[[[b[i][j] for j in perm] for i in perm] for b in pair])
                    allowed[n].add(key)
                    if family and family['coverage_status']=='unsat' and len(family['spaces'])==1:
                        spaces[n].setdefault(key,(c,ws[c['id']],family['spaces'][0],perm,swap))
    return allowed,spaces

def principal_blocks(candidate):
    p,m=candidate['N_plus_1'],candidate['N_minus_1'];d=candidate['symmetrizer'];n=len(d)
    for size in range(min(n-1,5),0,-1):
        for S in itertools.combinations(range(n),size):
            exterior=[k for k in range(n) if k not in S]
            if any(any(p[i][k] for i in S) and any(m[i][k] for i in S) for k in exterior):continue
            remaining=set(S)
            while remaining:
                component={min(remaining)}
                while True:
                    more=component|{j for i in component for j in remaining if p[i][j] or m[i][j] or p[j][i] or m[j][i]}
                    if more==component:break
                    component=more
                remaining-=component;T=sorted(component)
                key=label([d[i] for i in T],*[[[b[i][j] for j in T] for i in T] for b in (p,m)])
                yield S,T,key

def obstruction(candidate):
    allowed,_=tables()
    for S,T,key in principal_blocks(candidate):
        if len(T) in allowed and key not in allowed[len(T)]:
            return {'lemma':'Weighted principal T-datum under sign-pure exterior columns','subset':list(S),'component':T,'lower_rank':len(T),'labelled_constant_key':list(key)}
    return None

def add_relations(solver,variables,indices,candidate,z3):
    _,spaces=tables();seen=set()
    for S,T,key in principal_blocks(candidate):
        if tuple(T) in seen:continue
        seen.add(tuple(T));found=spaces.get(len(T),{}).get(key)
        if not found:continue
        c,w,space,perm,swap=found
        old_to_large={old:T[new] for new,old in enumerate(perm)};mapping=[]
        for name in w['variable_names']:
            if name.startswith('r'):mapping.append(old_to_large[int(name[1:])]);continue
            match=re.fullmatch(r'p([01])_(\d)(\d)_(\d+)',name);assert match,name
            s,i,j,k=map(int,match.groups());mapping.append(indices[s^swap,old_to_large[i],old_to_large[j]][k])
        solver.add(*(z3.Sum([z3.RealVal(c)*variables[k] for c,k in zip(row,mapping) if c!='0'])==0 for row in space['rref']))

def write_constant_table(path):
    allowed,_=tables();rows=[]
    for n,keys in allowed.items():
        flat={key[n:] for key in keys}
        rows.extend(str(n)+' '+' '.join(map(str,key)) for key in sorted(flat))
    Path(path).write_text('\n'.join(rows)+'\n',encoding='utf8',newline='\n')
    print({n:len({key[n:] for key in keys}) for n,keys in allowed.items()},flush=True)

def write_weighted_table(path):
    allowed,_=tables();rows=[]
    for n,keys in allowed.items():
        unique={}
        for key in keys:
            assert key[n:] not in unique or unique[key[n:]]==key[:n]
            unique[key[n:]]=key[:n]
        rows.extend(str(n)+' '+' '.join(map(str,unique[key]+key)) for key in sorted(unique))
    Path(path).write_text('\n'.join(rows)+'\n',encoding='utf8',newline='\n')

if __name__=='__main__':write_constant_table(HERE/'principal-constants.txt')
