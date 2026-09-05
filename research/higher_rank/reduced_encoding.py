"""Exact linear elimination before matching Laurent-exponent multisets."""
from collections import Counter
from functools import lru_cache
import itertools


def encode(solver,variables,equations,z3,sp):
    dim=len(variables);positions={v.get_id():i for i,v in enumerate(variables)}
    def linear(e):
        if e.get_id() in positions:return {positions[e.get_id()]:sp.Integer(1)}
        if z3.is_rational_value(e) or z3.is_int_value(e):return {dim:sp.Rational(str(e))}
        kind=e.decl().kind();children=e.children()
        if kind==z3.Z3_OP_TO_REAL:return linear(children[0])
        if kind==z3.Z3_OP_UMINUS:return {i:-v for i,v in linear(children[0]).items()}
        if kind in (z3.Z3_OP_ADD,z3.Z3_OP_SUB):
            out={}
            for k,c in enumerate(children):
                sign=-1 if kind==z3.Z3_OP_SUB and k else 1
                for i,v in linear(c).items():out[i]=out.get(i,0)+sign*v
            return {i:v for i,v in out.items() if v}
        if kind==z3.Z3_OP_MUL:
            factor=sp.Integer(1);other=None
            for c in children:
                values=linear(c)
                if all(i==dim for i in values):factor*=values.get(dim,0)
                else:
                    assert other is None,'Nonlinear equality';other=values
            return {i:factor*v for i,v in (other if other is not None else {dim:1}).items()}
        if kind==z3.Z3_OP_DIV:
            denominator=linear(children[1]);assert set(denominator)<={dim}
            return {i:v/denominator[dim] for i,v in linear(children[0]).items()}
        raise ValueError(('Nonlinear equality',str(e)))
    rows=[]
    for c in solver.assertions():
        if z3.is_eq(c):
            values=linear(c.arg(0)-c.arg(1));rows.append([values.get(i,0) for i in range(dim+1)])
    while True:
        matrix,pivots=sp.Matrix(rows).rref() if rows else (sp.zeros(0,dim+1),())
        if dim in pivots:solver.add(z3.BoolVal(False));return
        reduction=[(p,[(j,matrix[k,j]) for j in range(p+1,dim+1) if matrix[k,j]]) for k,p in enumerate(pivots)]
        @lru_cache(maxsize=None)
        def reduce(e):
            out={i:sp.Integer(v) for i,v in enumerate(e) if v}
            for p,row in reduction:
                coefficient=out.pop(p,0)
                if coefficient:
                    for j,c in row:
                        value=out.get(j,0)-coefficient*c
                        if value:out[j]=value
                        else:out.pop(j,None)
            return tuple(sorted(out.items()))
        reduced=[];extra=[]
        for left,right in equations:
            terms=Counter()
            for e in left:terms[reduce(e)]+=1
            for e in right:terms[reduce(e)]-=1
            a={e:c for e,c in terms.items() if c>0};b={e:-c for e,c in terms.items() if c<0}
            assert sum(a.values())==sum(b.values())
            if not a:continue
            if len(a)==1 or len(b)==1:
                one,many=(next(iter(a)),b) if len(a)==1 else (next(iter(b)),a)
                for e in many:
                    row=[0]*(dim+1)
                    for i,c in one:row[i]+=c
                    for i,c in e:row[i]-=c
                    extra.append(row)
            else:reduced.append((a,b))
        if not extra:break
        rows.extend(extra)
        for row in extra:solver.add(z3.Sum([z3.RealVal(str(c))*(variables[i] if i<dim else 1) for i,c in enumerate(row) if c])==0)
    def expression(e):return z3.Sum([z3.RealVal(str(c))*(variables[i] if i<dim else 1) for i,c in e]) if e else z3.RealVal(0)
    for group,(a,b) in enumerate(reduced):
        if len(a)>len(b):a,b=b,a
        left={e:expression(e) for e in a};right={e:expression(e) for e in b}
        if sum(a.values())<=3:
            aa=[left[e] for e,c in a.items() for _ in range(c)];bb=[right[e] for e,c in b.items() for _ in range(c)]
            solver.add(z3.Or(*[z3.And(*[x==y for x,y in zip(aa,permutation)]) for permutation in itertools.permutations(bb)]))
        else:
            for e,v in left.items():
                solver.add(z3.Or(*[v==f for f in right.values()]))
                solver.add(z3.Sum([z3.If(v==f,a[k],0) for k,f in left.items()])==z3.Sum([z3.If(v==f,b[k],0) for k,f in right.items()]))
            low,high=z3.Reals(f'extreme_{group}_min extreme_{group}_max')
            for side in (left,right):
                solver.add(*(low<=v for v in side.values()),*(high>=v for v in side.values()))
                solver.add(z3.Or(*[low==v for v in side.values()]),z3.Or(*[high==v for v in side.values()]))
