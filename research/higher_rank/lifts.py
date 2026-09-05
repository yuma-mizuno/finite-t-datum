"""Rank-independent lift spaces, exact periodicity and verification.

The full support/symplectic encoding is reused from the checked rank-four
engine. No rank-four-specific chain exclusion applies in higher ranks.
"""
from collections import Counter,defaultdict
import gzip
import hashlib
import itertools as it
import json
from math import gcd,lcm
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'rank4'))
from rank4_lifts import build_problem as base_problem,get_matrices,z3,sp,certificate,polynomial_data
from classify_families import matching_space,reduced
from obstructions import obstruction,add_monomial_relations


def build_problem(candidate,timeout,arithmetic='integer',encoding='multiplicity'):
    solver,variables,indices,equations=base_problem(candidate,timeout,arithmetic=arithmetic,encoding='multiplicity' if encoding in ('flow','principal') else encoding)
    if encoding=='flow':
        def contains_ite(e):return z3.is_app_of(e,z3.Z3_OP_ITE) or any(contains_ite(c) for c in e.children())
        plain=[c for c in solver.assertions() if not contains_ite(c)]
        solver=z3.Then('simplify','propagate-values','solve-eqs','simplify','smt').solver();solver.set(timeout=timeout);solver.add(*plain)
        expression=lambda e:z3.Sum([c*v for c,v in zip(e,variables) if c]) if any(e) else z3.IntVal(0)
        for group,(left,right) in enumerate(equations):
            a,b=Counter(left),Counter(right);aa=list(a);bb=list(b)
            flow=[[z3.Real(f'flow_{group}_{i}_{j}') for j in range(len(bb))] for i in range(len(aa))]
            for i,e in enumerate(aa):
                solver.add(z3.Sum(flow[i])==a[e])
                for j,f in enumerate(bb):solver.add(flow[i][j]>=0,z3.Or(flow[i][j]==0,expression(e)==expression(f)))
            for j,f in enumerate(bb):solver.add(z3.Sum([flow[i][j] for i in range(len(aa))])==b[f])
    add_monomial_relations(solver,variables,indices,candidate)
    if encoding=='principal':
        from principal_relations import add_principal_relations
        add_principal_relations(solver,variables,indices,candidate,z3)
    return solver,variables,indices,equations


def feasibility(candidate,timeout,arithmetic='integer',encoding='multiplicity'):
    start=time.monotonic()
    exclusion=obstruction(candidate)
    if exclusion:return {'id':candidate['id'],'status':'unsat','analytic_exclusion':exclusion,'seconds':time.monotonic()-start}
    solver,variables,indices,equations=build_problem(candidate,timeout,arithmetic=arithmetic,encoding=encoding)
    before=time.monotonic();status=solver.check()
    result={'id':candidate['id'],'status':str(status),'solver_seconds':time.monotonic()-before,
            'delay_bound':None,'encoding':encoding,'arithmetic':arithmetic,'variable_names':[str(v) for v in variables]}
    if status==z3.unknown:result['reason']=solver.reason_unknown()
    if status==z3.sat:
        values=integral_model(solver.model(),variables)
        a=matching_space(values,equations)
        values=small_witness(a,values,indices,len(candidate['N_plus_1']),timeout)
        ap,am=get_matrices(values,indices);polynomial_data(ap,am)
        result.update({'values':values,'A_plus':[[str(x) for x in row] for row in ap.tolist()],
                       'A_minus':[[str(x) for x in row] for row in am.tolist()]})
    result['seconds']=time.monotonic()-start
    return result


def integral_model(model,variables):
    values=[sp.Rational(str(model.eval(v))) for v in variables]
    scale=lcm(*(int(v.q) for v in values))
    return [int(v*scale) for v in values]


def small_witness(a,values,indices,rank,timeout):
    """Minimize inside a proved matching space, using only linear constraints."""
    variables=z3.Ints(' '.join(f'x{i}' for i in range(len(values))))
    opt=z3.Optimize();opt.set(timeout=timeout)
    opt.add(*(v>=1 for v in variables[:rank]))
    for (s,i,j),ks in indices.items():
        opt.add(*(z3.And(variables[k]>0,variables[k]<variables[i]) for k in ks))
        opt.add(*(variables[k]<=variables[l] for k,l in zip(ks,ks[1:])))
        if s==0:opt.add(*(variables[k]!=variables[l] for k in ks for l in indices[1,i,j]))
    opt.add(*(z3.Sum([z3.RealVal(str(c))*v for c,v in zip(row,variables)])==0 for row in a.tolist()))
    opt.minimize(z3.Sum(variables[:rank]));opt.minimize(z3.Sum(variables[rank:]))
    if opt.check()==z3.sat:return [opt.model().eval(v).as_long() for v in variables]
    return values


def spaces(candidate,witness,timeout,arithmetic='integer',encoding='multiplicity'):
    start=time.monotonic();rank=len(candidate['N_plus_1'])
    solver,variables,indices,equations=build_problem(candidate,timeout,arithmetic=arithmetic,encoding=encoding)
    found=[];values=witness['values'];status=z3.sat
    while status==z3.sat and len(found)<100:
        a=matching_space(values,equations)
        for left,right in equations:
            assert Counter(reduced(e,a) for e in left)==Counter(reduced(e,a) for e in right)
        found.append((a,values))
        solver.add(z3.Not(z3.And(*[z3.Sum([z3.RealVal(str(c))*v for c,v in zip(row,variables)])==0 for row in a.tolist()])))
        status=solver.check()
        if status==z3.sat:values=integral_model(solver.model(),variables)
    maximal=[]
    for a,values in found:
        if any(a.rows>b.rows and b.col_join(a).rank()==a.rows for b,_ in found):continue
        if any(a==b for b,_ in maximal):continue
        maximal.append((a,values))
    result={'id':candidate['id'],'coverage_status':str(status),'raw_space_count':len(found),
            'arithmetic':arithmetic,'encoding':encoding,'spaces':[]}
    if status==z3.unknown:result['reason']=solver.reason_unknown()
    for a,values in maximal:
        values=small_witness(a,values,indices,rank,timeout)
        gauge=[]
        for species in range(rank-1):
            v=[0]*len(variables)
            for (s,i,j),ks in indices.items():
                for k in ks:v[k]=int(i==species)-int(j==species)
            gauge.append(v)
        generators=sp.Matrix.hstack(*(sp.Matrix(v) for v in [values]+gauge))
        dimension=len(variables)-a.rows
        pure=dimension==rank and generators.rank()==rank and a*generators==sp.zeros(a.rows,rank)
        ap,am=get_matrices(values,indices)
        cert=certificate(ap,am,bound=300)
        result['spaces'].append({'dimension':dimension,'rref':[[str(c) for c in row] for row in a.tolist()],
                                'values':values,'scaling_and_shifts':pure,'certificate':cert})
    result['seconds']=time.monotonic()-start
    return result


def constant_check(candidate):
    start=time.monotonic();p,m=candidate['N_plus_1'],candidate['N_minus_1'];n=len(p)
    ap,am=(sp.Matrix([[2*(i==j)-a[i][j] for j in range(n)] for i in range(n)]) for a in (p,m))
    assert ap*am.T==am*ap.T
    assert all(0<=a[i][j]<=(1 if i==j else 2**n-1) for a in (p,m) for i in range(n) for j in range(n))
    assert all(int((ap*am.T)[i,i])%2==0 for i in range(n))
    reach=[[i==j or p[i][j]+m[i][j]>0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):reach[i][j]|=reach[i][k] and reach[k][j]
    assert all(all(row) for row in reach)
    for a in (ap,am):
        for k in range(1,n+1):
            for indices in it.combinations(range(n),k):assert a.extract(indices,indices).det()>0
    variables=z3.Reals(' '.join(f'v{i}' for i in range(n)));solver=z3.Solver()
    solver.add(*(v>=1 for v in variables))
    solver.add(*(z3.Sum([int(a[i,j])*variables[i] for i in range(n)])>=1 for a in (ap,am) for j in range(n)))
    assert solver.check()==z3.sat
    v=integral_model(solver.model(),variables);factor=0
    for x in v:factor=gcd(factor,x)
    v=[x//factor for x in v]
    assert all(x>0 for x in v) and all(sum(v[i]*int(a[i,j]) for i in range(n))>0 for a in (ap,am) for j in range(n))
    return {'id':candidate['id'],'status':'verified','positive_left_vector':v,'seconds':time.monotonic()-start}


def verify(candidate,witness,family,timeout,directory):
    start=time.monotonic();n=len(candidate['N_plus_1'])
    info=constant_check(candidate)
    arithmetic=(family or witness).get('arithmetic','integer')
    encoding=(family or witness).get('encoding','multiplicity')
    if witness.get('analytic_exclusion'):
        assert obstruction(candidate)==witness['analytic_exclusion']
        solver=z3.Solver();solver.add(z3.BoolVal(False));suffix='analytic_exclusion'
        info['analytic_exclusion']=witness['analytic_exclusion']
    else:solver,variables,indices,equations=build_problem(candidate,timeout,arithmetic=arithmetic,encoding=encoding)
    if witness['status']=='sat':
        assert family and family['coverage_status']=='unsat'
        membership=[]
        for space in family['spaces']:
            a=sp.Matrix([[sp.Rational(x) for x in row] for row in space['rref']])
            assert len(variables)-a.rank()==n and space['scaling_and_shifts']
            generators=[space['values']]
            for species in range(n-1):
                v=[0]*len(variables)
                for (s,i,j),ks in indices.items():
                    for k in ks:v[k]=int(i==species)-int(j==species)
                generators.append(v)
            generators=sp.Matrix.hstack(*(sp.Matrix(v) for v in generators))
            assert generators.rank()==n and a*generators==sp.zeros(a.rows,n)
            for left,right in equations:
                assert Counter(reduced(e,a) for e in left)==Counter(reduced(e,a) for e in right)
            ap,am=get_matrices(space['values'],indices)
            fresh=certificate(ap,am,bound=300)
            assert json.loads(json.dumps(fresh))==space['certificate']
            assert fresh['positive'] and fresh['negative'] and fresh['labelled_tropical_seed_period']
            membership.append(z3.And(*[z3.Sum([z3.RealVal(str(c))*v for c,v in zip(row,variables)])==0 for row in a.tolist()]))
        solver.add(z3.Not(z3.Or(*membership)));suffix='coverage'
    else:
        assert witness['status']=='unsat'
        if not witness.get('analytic_exclusion'):suffix='exclusion'
    query=solver.to_smt2();target=Path(directory)/'smt_queries';target.mkdir(exist_ok=True)
    filename=f'{candidate["id"]:06d}_{suffix}.smt2.gz'
    (target/filename).write_bytes(gzip.compress(query.encode(),mtime=0))
    status=solver.check()
    info.update({'result':str(status),'arithmetic':arithmetic,'query_file':filename,
                 'sha256':hashlib.sha256(query.encode()).hexdigest(),'seconds':time.monotonic()-start})
    if status==z3.unknown:info['reason']=solver.reason_unknown()
    return info


def task(stage,candidate,witness,family,timeout,arithmetic,directory,encoding='multiplicity'):
    try:
        if stage=='lifts':return feasibility(candidate,timeout,arithmetic,encoding)
        if stage=='families':return spaces(candidate,witness,timeout,arithmetic,encoding)
        if stage=='verify':return verify(candidate,witness,family,timeout,directory)
        raise ValueError(stage)
    except Exception as exc:
        import traceback
        return {'id':candidate['id'],'status':'error','error':str(exc),'traceback':traceback.format_exc()}
