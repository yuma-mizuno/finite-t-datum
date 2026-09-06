"""Exact polynomial lifts with a primitive positive right symmetrizer.

An atom in entry (i,j) has coefficient d_i/gcd(d_i,d_j). Thus dual
integrality is enforced at every exponent, including coincident atoms.
Each worker configures the existing family/coverage engine for one candidate.
"""
from collections import Counter
from math import gcd
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'higher_rank'))
import lifts as engine
from classify_rank3 import step,identity,negative_permutation
from reduced_encoding import encode
from obstructions import obstruction as identity_obstruction,add_monomial_relations
from principal import obstruction as principal_obstruction,add_relations
sp,z3,z=engine.sp,engine.z3,engine.sp.Symbol('z')


class Indices(dict):
    pass


def build_problem(candidate,timeout,arithmetic='integer',encoding='reduced'):
    if encoding not in ('multiplicity','reduced'):raise ValueError('Unsupported weighted Laurent-identity encoding: '+encoding)
    counts=[candidate['N_plus_1'],candidate['N_minus_1']];d=candidate['symmetrizer'];n=len(d)
    names=[f'r{i}' for i in range(n)];indices=Indices();indices.weights={};indices.symmetrizer=d
    for s in range(2):
        for i in range(n):
            for j in range(n):
                weight=d[i]//gcd(d[i],d[j]);assert counts[s][i][j]%weight==0
                indices[s,i,j]=[];indices.weights[s,i,j]=weight
                for k in range(counts[s][i][j]//weight):
                    indices[s,i,j].append(len(names));names.append(f'p{s}_{i}{j}_{k}')
    variables=[(z3.Int if arithmetic=='integer' else z3.Real)(name) for name in names]
    dim=len(variables);zero=(0,)*dim;basis=[tuple(int(i==j) for i in range(dim)) for j in range(dim)]
    a=[[[[(-indices.weights[s,i,j],basis[k]) for k in indices[s,i,j]]+([(1,zero),(1,basis[i])] if i==j else []) for j in range(n)] for i in range(n)] for s in range(2)]
    solver=z3.Solver();solver.set(timeout=timeout)
    solver.add(*(v>=1 if arithmetic=='integer' else v>0 for v in variables[:n]))
    if arithmetic!='integer':solver.add(z3.Sum(variables[:n])==1)
    for (s,i,j),ks in indices.items():
        solver.add(*(z3.And(variables[k]>0,variables[k]<variables[i]) for k in ks))
        solver.add(*(variables[k]<=variables[l] for k,l in zip(ks,ks[1:])))
        if s==0:solver.add(*(variables[k]!=variables[l] for k in ks for l in indices[1,i,j]))
    expression=lambda e:z3.Sum([c*v for c,v in zip(e,variables) if c]) if any(e) else z3.IntVal(0)
    equations=[]
    for i in range(n):
        for j in range(i,n):
            terms=Counter()
            for k in range(n):
                for s in range(2):
                    for c,e in a[s][i][k]:
                        for f,t in a[1-s][j][k]:
                            terms[tuple(x-y for x,y in zip(e,t))]+=(1-2*s)*c*f*d[k]
            divisor=0
            for c in terms.values():divisor=gcd(divisor,abs(c))
            if not divisor:continue
            left=[e for e,c in sorted(terms.items()) for _ in range(max(c//divisor,0))]
            right=[e for e,c in sorted(terms.items()) for _ in range(max(-c//divisor,0))]
            assert len(left)==len(right)
            equations.append((left,right));aa=list(map(expression,left));bb=list(map(expression,right))
            solver.add(z3.Sum(aa)==z3.Sum(bb))
            if encoding=='multiplicity':
                for e in set(left):
                    v=expression(e);solver.add(z3.Sum([z3.If(v==f,1,0) for f in aa])==z3.Sum([z3.If(v==f,1,0) for f in bb]))
    add_relations(solver,variables,indices,candidate,z3)
    if d==[1]*n:add_monomial_relations(solver,variables,indices,candidate)
    if encoding=='reduced':encode(solver,variables,equations,z3,sp)
    return solver,variables,indices,equations


def get_matrices(values,indices):
    n=len(indices.symmetrizer);result=[]
    for s in range(2):
        a=sp.diag(*(1+z**r for r in values[:n]))
        for i in range(n):
            for j in range(n):a[i,j]-=indices.weights[s,i,j]*sum(z**values[k] for k in indices[s,i,j])
        result.append(a)
    return result


def polynomial_data(ap,am,d):
    n=len(d);assert len(d)==ap.rows==ap.cols==am.rows==am.cols
    assert all(isinstance(x,int) and x>0 for x in d) and gcd(*d)==1
    delays=[int(sp.degree(ap[i,i],z)) for i in range(n)];n0=sp.diag(*(1+z**r for r in delays));signs=[]
    for a in (ap,am):
        data={}
        for i in range(n):
            for j in range(n):
                for (p,),c in sp.Poly(sp.expand(n0[i,j]-a[i,j]),z).terms():
                    if c:
                        assert c.is_Integer and c>0 and 0<p<delays[i]
                        assert (int(c)*d[j])%d[i]==0
                        data[i,j,p]=int(c)
        signs.append(data)
    assert not(signs[0].keys()&signs[1].keys())
    D=sp.diag(*d)
    assert all(sp.expand(c)==0 for c in ap*D*am.subs(z,1/z).T-am*D*ap.subs(z,1/z).T)
    return delays,signs


def exchange_matrix(delays,signs,d):
    vertices=[(i,p) for i,r in enumerate(delays) for p in range(r)];plus,minus=signs
    dual=[{(i,j,p):c*d[j]//d[i] for (i,j,p),c in sign.items()} for sign in signs]
    pd,md=dual;b=[]
    for i,p in vertices:
        row=[]
        for j,q in vertices:
            value=-plus.get((i,j,p-q),0)+minus.get((i,j,p-q),0)+pd.get((j,i,q-p),0)-md.get((j,i,q-p),0)
            for k in range(len(d)):
                for v in range(min(p,q)+1):value+=plus.get((i,k,p-v),0)*md.get((j,k,q-v),0)-minus.get((i,k,p-v),0)*pd.get((j,k,q-v),0)
            row.append(value)
        b.append(row)
    weights=[d[i] for i,p in vertices]
    assert all(b[i][j]*weights[j]==-b[j][i]*weights[i] for i in range(len(b)) for j in range(len(b)))
    lookup={v:k for k,v in enumerate(vertices)};permutation=[lookup[i,(p-1)%delays[i]] for i,p in vertices];mutation=[lookup[i,0] for i in range(len(d))]
    assert all(b[i][j]==0 for i in mutation for j in mutation)
    return vertices,b,permutation,mutation


def certificate(ap,am,d,bound=300):
    delays,signs=polynomial_data(ap,am,d);vertices,original,permutation,mutation=exchange_matrix(delays,signs,d);n=len(original)
    result={'delays':delays,'symmetrizer':d,'vertices':vertices,'B':original,'mutation_vertices':mutation,'relabel_old_to_new':permutation,
            'A_plus':[[str(x) for x in row] for row in ap.tolist()],'A_minus':[[str(x) for x in row] for row in am.tolist()]}
    b,c=step(original,identity(n),permutation,mutation,1)
    assert b==original and step(b,c,permutation,mutation,-1)==(original,identity(n))
    for direction,label in ((1,'positive'),(-1,'negative')):
        b,c=original,identity(n)
        for h in range(1,bound+1):
            b,c=step(b,c,permutation,mutation,direction);assert b==original
            if all(x<=0 for row in c for x in row):
                assert negative_permutation(c)
                perm=[row.index(-1) for row in c]
                assert all(d[vertices[i][0]]==d[vertices[perm[i]][0]] for i in range(n))
                result[label]={'h':h,'negative_permutation':perm};break
        else:result[label]=None
    b,c=original,identity(n)
    for period in range(1,2*bound+1):
        b,c=step(b,c,permutation,mutation,1)
        if c==identity(n):result['labelled_tropical_seed_period']=period;break
    else:result['labelled_tropical_seed_period']=None
    return result


def constant_check(candidate):
    import time,itertools
    start=time.monotonic();p,m=candidate['N_plus_1'],candidate['N_minus_1'];d=candidate['symmetrizer'];n=len(d);D=sp.diag(*d)
    ap,am=(2*sp.eye(n)-sp.Matrix(b) for b in (p,m));assert ap*D*am.T==am*D*ap.T
    assert gcd(*d)==1 and min(d)>0
    for i in range(n):
        assert (ap*D*am.T)[i,i]%(2*d[i])==0
        for j in range(n):
            for b in (p,m):
                assert 0<=b[i][j]<=(1 if i==j else 2**n-1)
                assert b[i][j]*d[j]%d[i]==0 and b[i][j]*d[j]//d[i]<=2**n-1
    for a in (ap,am):
        for size in range(1,n+1):
            for S in itertools.combinations(range(n),size):assert a.extract(S,S).det()>0
    vs=z3.Reals(' '.join(f'v{i}' for i in range(n)));solver=z3.Solver();solver.add(*(v>=1 for v in vs))
    solver.add(*(z3.Sum([int(a[i,j])*vs[i] for i in range(n)])>=1 for a in (ap,am) for j in range(n)))
    assert solver.check()==z3.sat;v=engine.integral_model(solver.model(),vs)
    return {'id':candidate['id'],'status':'verified','symmetrizer':d,'positive_left_vector':v,'seconds':time.monotonic()-start}


def diagonal_partition_obstruction(candidate):
    """Necessary equality of positive/negative Laurent masses, in integral dual atoms."""
    from math import prod
    d=candidate['symmetrizer'];p,m=candidate['N_plus_1'],candidate['N_minus_1']
    for i in range(len(d)):
        terms=[]
        for k in range(len(d)):
            g=d[i]//gcd(d[i],d[k]);assert p[i][k]%g==m[i][k]%g==0
            mass=(p[i][k]//g)*(m[i][k]//g)
            if mass:terms.append((k,mass,d[k]*g*g))
        if not any(c%2 for k,c,w in terms):continue
        if len(terms)!=1 and prod(c+1 for k,c,w in terms)>4096:continue
        possible={0}
        for k,c,w in terms:possible={a+w*t for a in possible for t in range(-c,c+1,2)}
        if 0 not in possible:
            return {'lemma':'Diagonal Laurent cross-products have equal positive and negative mass',
                    'row':i,'column_atom_products_and_weights':[list(t) for t in terms],
                    'dependency':'Disjoint opposite signs, strict row support, dual integrality and the diagonal symplectic identity.'}
    return None


def task(stage,candidate,witness,family,timeout,arithmetic,directory,encoding='reduced'):
    d=candidate['symmetrizer']
    engine.build_problem=build_problem;engine.get_matrices=get_matrices
    engine.polynomial_data=lambda ap,am:polynomial_data(ap,am,d)
    engine.certificate=lambda ap,am,bound=300:certificate(ap,am,d,bound)
    engine.constant_check=constant_check
    engine.obstruction=lambda c:(identity_obstruction(c) if d==[1]*len(d) else None) or principal_obstruction(c) or diagonal_partition_obstruction(c)
    # In this weighted entry point the requested replay encoding is explicit.
    # The legacy engine interprets "multiplicity" as "reuse the witness encoding".
    if stage=='verify' and encoding=='multiplicity':
        return engine.verify(candidate,witness,family,timeout,directory,encoding_override='multiplicity')
    return engine.task(stage,candidate,witness,family,timeout,arithmetic,directory,encoding)
