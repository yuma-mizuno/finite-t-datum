"""Cover all admissible lifts by exact linear spaces, then certify dynamics."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time

import sympy as sp
from rank4_lifts import build_problem, get_matrices, z3, certificate

HERE=Path(__file__).resolve().parent


def matching_space(values,equations):
    rows=[]
    for positive,negative in equations:
        left,right=defaultdict(list),defaultdict(list)
        for e in positive:left[sum(x*y for x,y in zip(e,values))].append(e)
        for e in negative:right[sum(x*y for x,y in zip(e,values))].append(e)
        assert left.keys()==right.keys()
        for v in sorted(left):
            assert len(left[v])==len(right[v])
            rows.extend([[x-y for x,y in zip(a,b)] for a,b in zip(sorted(left[v]),sorted(right[v]))])
    if not rows:return sp.zeros(0,len(values))
    rref,pivots=sp.Matrix(rows).rref()
    return rref[:len(pivots),:]


def reduced(e,rref):
    v=list(map(sp.Rational,e))
    for row in rref.tolist():
        pivot=next(i for i,x in enumerate(row) if x)
        t=v[pivot]
        v=[x-t*y for x,y in zip(v,row)]
    return tuple(v)


def classify(candidate,witness,timeout):
    start=time.monotonic()
    solver,variables,indices,equations=build_problem(candidate,timeout)
    spaces=[]
    values=witness['values']
    status=z3.sat
    while status==z3.sat:
        space=matching_space(values,equations)
        for left,right in equations:
            assert Counter(reduced(e,space) for e in left)==Counter(reduced(e,space) for e in right)
        spaces.append((space,values))
        solver.add(z3.Not(z3.And(*[z3.Sum([z3.RealVal(str(c))*v for c,v in zip(row,variables)])==0
                                  for row in space.tolist()])))
        status=solver.check()
        if status==z3.sat:
            values=[solver.model().eval(v).as_long() for v in variables]
        if len(spaces)>100:raise RuntimeError('More than 100 spaces')
    cid=candidate['id']
    maximal=[]
    for a,values in spaces:
        if any(a.rows>b.rows and b.col_join(a).rank()==a.rows for b,_ in spaces):continue
        if any(a==b for b,_ in maximal):continue
        maximal.append((a,values))
    result={'id':cid,'coverage_status':str(status),
            'raw_space_count':len(spaces),'spaces':[]}
    for a,values in maximal:
        gauge=[]
        for species in range(3):
            v=[0]*len(variables)
            for (s,i,j),ks in indices.items():
                for k in ks:v[k]=int(i==species)-int(j==species)
            gauge.append(v)
        generators=sp.Matrix.hstack(*(sp.Matrix(v) for v in [values]+gauge))
        dimension=len(variables)-a.rows
        pure=(dimension==4 and generators.rank()==4 and a*generators==sp.zeros(a.rows,4))
        ap,am=get_matrices(values,indices)
        cert=certificate(ap,am,bound=150)
        result['spaces'].append({'dimension':dimension,'rref':[[str(c) for c in row] for row in a.tolist()],
                                 'values':values,'scaling_and_shifts':pure,'certificate':cert})
    result['seconds']=time.monotonic()-start
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--timeout',type=int,default=15000)
    args=parser.parse_args()
    candidates=json.loads((HERE/'constant_candidates.json').read_text())['candidates']
    witnesses={x['id']:x for x in map(json.loads,(HERE/'lift_feasibility.jsonl').read_text().splitlines())}
    target=HERE/'families.jsonl'
    previous={}
    if target.exists():previous={x['id']:x for x in map(json.loads,target.read_text().splitlines())}
    todo=[c for c in candidates if witnesses.get(c['id'],{}).get('status')=='sat' and
          (c['id'] not in previous or previous[c['id']]['coverage_status']=='unknown')]
    (HERE/'smt_queries').mkdir(exist_ok=True)
    with target.open('a',encoding='utf-8') as out,ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(classify,c,witnesses[c['id']],args.timeout):c['id'] for c in todo}
        for future in as_completed(futures):
            x=future.result();out.write(json.dumps(x)+'\n');out.flush();previous[x['id']]=x
            print(x['id'],x['coverage_status'],[(s['dimension'],s['scaling_and_shifts'],
                  (s['certificate'].get('positive') or {}).get('h'),
                  (s['certificate'].get('negative') or {}).get('h'),
                  s['certificate'].get('labelled_tropical_seed_period')) for s in x['spaces']],
                  round(x['seconds'],2),flush=True)
    target.write_text(''.join(json.dumps(previous[k])+'\n' for k in sorted(previous)),encoding='utf-8',newline='\n')


if __name__=='__main__':main()
