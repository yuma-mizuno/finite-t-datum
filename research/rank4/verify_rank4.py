"""Independent constant checks and replay of all exclusion/coverage formulas."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import itertools as it
import json
from math import gcd,lcm
from functools import reduce
from pathlib import Path
import platform
import subprocess
import sys

import sympy as sp
from rank4_lifts import build_problem,z3,certificate,get_matrices
from classify_families import reduced

HERE=Path(__file__).resolve().parent


def det(a):
    if not a:return 1
    return sum((-1)**j*a[0][j]*det([r[:j]+r[j+1:] for r in a[1:]]) for j in range(len(a)))


def key(p,m):
    return min(tuple(a[q[i]][q[j]] for a in pair for i in range(4) for j in range(4))
               for q in it.permutations(range(4)) for pair in ((p,m),(m,p)))


def check_constant(c):
    p,m=c['N_plus_1'],c['N_minus_1']
    assert all(len(a)==4 and all(len(row)==4 for row in a) for a in (p,m))
    assert all(0<=a[i][j]<=(1 if i==j else 15) for a in (p,m) for i in range(4) for j in range(4))
    assert key(p,m)==tuple(x for a in (p,m) for row in a for x in row)
    reach=[[i==j or p[i][j]+m[i][j]>0 for j in range(4)] for i in range(4)]
    for k in range(4):
        for i in range(4):
            for j in range(4):reach[i][j]|=reach[i][k] and reach[k][j]
    assert all(all(row) for row in reach)
    ap,am=([[2*(i==j)-a[i][j] for j in range(4)] for i in range(4)] for a in (p,m))
    assert all(sum(ap[i][k]*am[j][k] for k in range(4))==sum(am[i][k]*ap[j][k] for k in range(4))
               for i in range(4) for j in range(4))
    assert all(det([[a[i][j] for j in s] for i in s])>0 for a in (ap,am)
               for k in range(1,5) for s in it.combinations(range(4),k))
    v=z3.Reals('v0 v1 v2 v3');solver=z3.Solver()
    solver.add(*(x>=1 for x in v))
    solver.add(*(z3.Sum([v[i]*a[i][j] for i in range(4)])>=1 for a in (ap,am) for j in range(4)))
    assert solver.check()==z3.sat
    w=[sp.Rational(str(solver.model().eval(x))) for x in v]
    scale=lcm(*(int(x.q) for x in w))
    w=[int(x*scale) for x in w];g=reduce(gcd,w);w=[x//g for x in w]
    assert all(x>0 for x in w)
    assert all(sum(w[i]*a[i][j] for i in range(4))>0 for a in (ap,am) for j in range(4))
    parity=[sum(ap[i][j]*am[i][j] for j in range(4))%2 for i in range(4)]
    return {'id':c['id'],'positive_left_vector':w,'odd_diagonal_rows':[i+1 for i,x in enumerate(parity) if x]}


def replay(c,item,family,timeout):
    cid=c['id']
    encoding=item.get('encoding','multiplicity');arithmetic=item.get('arithmetic','integer')
    solver,variables,indices,equations=build_problem(c,timeout,encoding,arithmetic)
    if item['status']=='sat':
        assert family['coverage_status']=='unsat' and len(family['spaces'])==1
        space=family['spaces'][0]
        a=sp.Matrix([[sp.Rational(x) for x in row] for row in space['rref']])
        assert len(variables)-a.rank()==4 and space['scaling_and_shifts']
        gauge=[]
        for species in range(3):
            v=[0]*len(variables)
            for (s,i,j),ks in indices.items():
                for k in ks:v[k]=int(i==species)-int(j==species)
            gauge.append(v)
        generators=sp.Matrix.hstack(*(sp.Matrix(v) for v in [space['values']]+gauge))
        assert generators.rank()==4 and a*generators==sp.zeros(a.rows,4)
        for left,right in equations:
            assert Counter(reduced(e,a) for e in left)==Counter(reduced(e,a) for e in right)
        ap,am=get_matrices(space['values'],indices)
        z=sp.Symbol('z')
        assert (2*sp.eye(4)-ap.subs(z,1)).tolist()==c['N_plus_1']
        assert (2*sp.eye(4)-am.subs(z,1)).tolist()==c['N_minus_1']
        fresh=certificate(ap,am,bound=150)
        assert json.loads(json.dumps(fresh))==space['certificate']
        assert fresh['positive'] and fresh['negative'] and fresh['labelled_tropical_seed_period']
        solver.add(z3.Not(z3.And(*[z3.Sum([z3.RealVal(str(x))*v for x,v in zip(row,variables)])==0
                                  for row in a.tolist()])))
        suffix='coverage'
    else:
        assert item['status']=='unsat';suffix='exclusion'
    query=solver.to_smt2()
    name=f'{cid:04d}_{suffix}.smt2'
    (HERE/'smt_queries'/name).write_text(query,encoding='utf-8',newline='\n')
    status=solver.check()
    return {'id':cid,'file':name,'result':str(status),'encoding':encoding,'arithmetic':arithmetic,
            'sha256':hashlib.sha256(query.encode()).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--constants-only',action='store_true')
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--timeout',type=int,default=120000)
    args=parser.parse_args()
    cs=json.loads((HERE/'constant_candidates.json').read_text())['candidates']
    constants=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,x in enumerate(pool.map(check_constant,cs,chunksize=20),1):
            constants.append(x)
            if i%1000==0:print('Constant checks',i,flush=True)
    (HERE/'constant_verification.json').write_text(json.dumps(constants,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert len(cs)==4865 and len({key(c['N_plus_1'],c['N_minus_1']) for c in cs})==4865
    odd=sum(bool(x['odd_diagonal_rows']) for x in constants)
    print('Constants verified',len(cs),'parity exclusions',odd,flush=True)
    if args.constants_only:return
    fs={x['id']:x for x in map(json.loads,(HERE/'lift_feasibility.jsonl').read_text().splitlines())}
    families={x['id']:x for x in map(json.loads,(HERE/'families.jsonl').read_text().splitlines())}
    assert len(fs)==4865 and all(x['status'] in ('sat','unsat') for x in fs.values())
    queries=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(replay,c,fs[c['id']],families.get(c['id']),args.timeout):c['id'] for c in cs}
        for i,f in enumerate(as_completed(futures),1):
            x=f.result();queries.append(x)
            if i%500==0 or x['result']!='unsat':print('Replay',i,x['id'],x['result'],flush=True)
    for x in queries:
        if x['result']!='unknown':continue
        completed=subprocess.run([sys.executable,str(HERE/'check_smt.py'),str(HERE/'smt_queries'/x['file']),
                                  '--timeout',str(max(args.timeout,180000))],capture_output=True,text=True,check=True)
        retry=json.loads(completed.stdout)
        assert retry['sha256']==x['sha256']
        x['result']=retry['result'];x['isolated_replay']=retry
        print('Isolated replay',x['id'],x['result'],flush=True)
    report={'python':platform.python_version(),'sympy':sp.__version__,'z3':z3.get_version_string(),
            'constant_pairs':len(cs),'parity_exclusions':odd,
            'lift_status_counts':dict(Counter(x['status'] for x in fs.values())),
            'replay_counts':dict(Counter(x['result'] for x in queries)),
            'queries':sorted(queries,key=lambda x:x['id'])}
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert all(x['result']=='unsat' for x in queries)
    print({k:v for k,v in report.items() if k!='queries'},flush=True)


if __name__=='__main__':main()
