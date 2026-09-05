"""Hereditary polynomial-lift obstructions from complete lower-rank results."""
import itertools as it
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def key(p,m):
    n=len(p)
    return min(tuple(b[q[i]][q[j]] for b in pair for i in range(n) for j in range(n))
               for q in it.permutations(range(n)) for pair in ((p,m),(m,p)))


def lower_keys():
    records=json.loads((HERE.parents[1]/'docs/catalogue/catalogue.json').read_text(encoding='utf8'))['records']
    valid={n:set() for n in range(1,6)}
    for r in records:
        n=r['rank']
        if n>5:continue
        p,m=([[(2 if i==j else 0)-r['datum'][a][i][j] for j in range(n)] for i in range(n)] for a in ('A_plus_1','A_minus_1'))
        valid[n].add(key(p,m))
    assert [len(valid[n]) for n in range(1,6)]==[2,6,16,37,55]
    return valid


VALID=lower_keys()


def hereditary(candidate):
    p,m=candidate['N_plus_1'],candidate['N_minus_1'];n=len(p)
    for size in range(min(n-1,5),0,-1):
        for S in it.combinations(range(n),size):
            external=[k for k in range(n) if k not in S]
            if any(any(p[i][k] for i in S) and any(m[i][k] for i in S) for k in external):continue
            remaining=set(S)
            while remaining:
                component={min(remaining)}
                while True:
                    more=component|{j for i in component for j in remaining if p[i][j]+m[i][j]+p[j][i]+m[j][i]}
                    if more==component:break
                    component=more
                remaining-=component;T=sorted(component)
                restricted=key([[p[i][j] for j in T] for i in T],[[m[i][j] for j in T] for i in T])
                if restricted not in VALID[len(T)]:
                    return {'lemma':'sign-pure external columns preserve the principal T-datum',
                            'subset':list(S),'offending_component':T,'lower_rank':len(T),'canonical_constant_key':list(restricted),
                            'external_columns':[{'column':k,'zero_sign':'plus' if not any(p[i][k] for i in S) else 'minus'} for k in external],
                            'dependency':'Complete polynomial-lift classification under simultaneous positivity in ranks at most '+('five.' if n>5 else 'four.')}
    return None


def leaf_extension(candidate):
    """Extremal-degree obstruction for a monomial tree with one attachment."""
    n=len(candidate['N_plus_1'])
    for sign in range(2):
        p,m=(candidate[x] for x in (('N_plus_1','N_minus_1') if sign==0 else ('N_minus_1','N_plus_1')))
        for leaf in range(n):
            top=[i for i in range(n) if i!=leaf]
            if any(p[i][j] for i in top for j in range(n)):continue
            if any(m[i][j] not in (0,1) or m[i][j]!=m[j][i] for i in top for j in top):continue
            edges=[(i,j) for i in top for j in top if i<j and m[i][j]]
            if len(edges)!=len(top)-1:continue
            seen={top[0]}
            while True:
                more=seen|{j for i in seen for j in top if m[i][j]}
                if more==seen:break
                seen=more
            if seen!=set(top):continue
            attachments=[i for i in top if m[i][leaf]]
            if len(attachments)!=1 or m[attachments[0]][leaf]!=1:continue
            a=attachments[0];reverse=m[leaf][a]
            if p[leaf][leaf] and reverse:continue
            if reverse>1:continue
            if reverse==1 and not any(p[leaf][i] or m[leaf][i] for i in top if i!=a) and not p[leaf][a]:continue
            return {'lemma':'monomial-tree leaf-extension extremal-degree obstruction','sign_exchange':bool(sign),
                    'tree_species':top,'extra_species':leaf,'attachment_species':a,'reverse_coefficient':reverse,
                    'tree_edges':[list(e) for e in edges]}
    return None


def obstruction(candidate):
    preliminary=hereditary(candidate) or leaf_extension(candidate) or single_attachment(candidate)
    if preliminary:return preliminary
    from finite_extension import exclusion
    return exclusion(candidate)


def attachments(candidate):
    counts=[candidate['N_plus_1'],candidate['N_minus_1']];n=len(counts[0])
    for s in range(2):
        for ell in range(n):
            top=[i for i in range(n) if i!=ell]
            if any(counts[s][i][ell] for i in top):continue
            aa=[i for i in top if counts[1-s][i][ell]]
            if len(aa)==1 and counts[1-s][aa[0]][ell]==1:yield s,ell,aa[0],top


def single_attachment(candidate):
    counts=[candidate['N_plus_1'],candidate['N_minus_1']]
    for s,ell,a,top in attachments(candidate):
        reverse=counts[1-s][ell][a]
        if not reverse or (reverse==1 and any(counts[t][ell][i] for t in range(2) for i in top if (t,i)!=(1-s,a))):
            return {'lemma':'single monomial attachment extremal-support obstruction','sign_exchange':bool(s),
                    'extra_species':ell,'attachment_species':a,'reverse_coefficient':reverse}
    return None


def add_monomial_relations(solver,variables,indices,candidate):
    """Redundant linear consequences of the exact reciprocal identities."""
    counts=[candidate['N_plus_1'],candidate['N_minus_1']];n=len(counts[0])
    for s,ell,a,top in attachments(candidate):
        b=variables[indices[1-s,a,ell][0]];reverse=indices[1-s,ell,a]
        if not reverse:solver.add(False);continue
        low=variables[a]-b;high=variables[ell]-b
        solver.add(variables[reverse[0]]==low,variables[reverse[-1]]==high)
        if len(reverse)>1:
            solver.add(variables[reverse[0]]<variables[reverse[1]],variables[reverse[-2]]<variables[reverse[-1]])
        for t in range(2):
            for i in top:
                if (t,i)==(1-s,a):continue
                solver.add(*(variables[k]>low for k in indices[t,ell,i]),*(variables[k]<high for k in indices[t,ell,i]))
    # (1+z^R)(T-U)=z^a V, with T,U disjoint and V of mass two.
    # Finite support forces one odd alternating chain. Its length is T(1)+U(1).
    for s in range(2):
        for i in range(n):
            if any(counts[s][i]):continue
            js=[j for j in range(n) if counts[1-s][i][j]]
            if len(js)!=1 or js[0]==i or counts[1-s][i][js[0]]!=1:continue
            j=js[0];a=variables[indices[1-s,i,j][0]];R=variables[i]
            for ell in range(n):
                if ell in (i,j):continue
                mass=counts[1-s][ell][i]
                if counts[s][ell][i]!=mass+1 or counts[s][ell][j]!=2:continue
                P=[variables[k] for k in indices[s,ell,i]]
                M=[variables[k] for k in indices[1-s,ell,i]]
                V=[variables[k] for k in indices[s,ell,j]]
                solver.add(*(x==P[0]+2*k*R for k,x in enumerate(P)))
                solver.add(*(x==P[0]+(2*k+1)*R for k,x in enumerate(M)))
                solver.add(V[0]==P[0]+a-R,V[1]==P[-1]+a)
    for s in range(2):
        for i,j in it.combinations(range(n),2):
            if (counts[1-s][i][j]==counts[1-s][j][i]==1
                and counts[s][i][j]==counts[s][j][i]==0
                and not any(counts[0][i][k]*counts[1][j][k] or counts[1][i][k]*counts[0][j][k] for k in range(n))):
                a=indices[1-s,i,j][0];b=indices[1-s,j,i][0]
                solver.add(variables[i]==variables[j],variables[a]+variables[b]==variables[i])
    for i in range(n):
        if any(counts[0][i][j]*counts[1][i][j] for j in range(n) if j!=i):continue
        for s in range(2):
            if counts[s][i][i]==1 and counts[1-s][i][i]==0:
                solver.add(2*variables[indices[s,i,i][0]]==variables[i])
