"""Exact rank-generic RSG, SG and ADE/tadpole identifications."""
import importlib.util
import itertools as it
import json
from pathlib import Path
import sys
import sympy as S

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('lower_family_notes',HERE.parent/'catalogue/family_notes.py')
lower=importlib.util.module_from_spec(spec);spec.loader.exec_module(lower)
z=lower.z

def constants(pair):return [tuple(map(int,a.subs(z,1))) for a in pair]
def canonical(pair):
    n=pair[0].rows;c=constants(pair)
    return min(tuple(c[s][p[i]*n+p[j]] for s in signs for i in range(n) for j in range(n))
               for p in it.permutations(range(n)) for signs in ((0,1),(1,0)))

def match(source,target):
    n=source[0].rows;sc=constants(source);tc=constants(target);sr=[S.degree(source[0][i,i],z) for i in range(n)]
    for swap in (False,True):
        pair=source[::-1] if swap else source
        for p in it.permutations(range(n)):
            if any(sc[s^int(swap)][p[i]*n+p[j]]!=tc[s][i*n+j] for s in range(2) for i in range(n) for j in range(n)):continue
            arrays=[a.extract(p,p) for a in pair];lam=S.Rational(S.degree(target[0][0,0],z),sr[p[0]])
            if any(lam*sr[p[i]]!=S.degree(target[0][i,i],z) for i in range(n)):continue
            ss=S.symbols('s:'+str(n));eq=[ss[-1]]
            for a,b in zip(arrays,target):
                for i in range(n):
                    for j in range(n):
                        if i!=j and a[i,j]!=0:
                            amin=min(e[0] for e,c in S.Poly(a[i,j],z).terms());bmin=min(e[0] for e,c in S.Poly(b[i,j],z).terms())
                            eq.append(ss[i]-ss[j]-(bmin-lam*amin))
            for shift in S.linsolve(eq,ss):
                if any(s.free_symbols for s in shift):continue
                if all(S.expand(sum(c*z**(lam*e[0]+shift[i]-shift[j]) for e,c in S.Poly(a[i,j],z).terms())-b[i,j])==0
                       for a,b in zip(arrays,target) for i in range(n) for j in range(n)):
                    return {'direction':'named constructor to displayed representative','sign_exchange':swap,'permutation':list(p),'lambda':str(lam),'shifts':list(map(str,shift))}

def adj(label,n):
    a=S.zeros(n)
    edges=([(i,i+1) for i in range(n-3)]+[(n-3,n-2),(n-3,n-1)] if label=='D' else [(i,i+1) for i in range(n-1)])
    if label=='E':edges=[(0,1),(1,2),(2,3),(3,4),(2,5)]+([(4,6)] if n>=7 else [])+([(6,7)] if n==8 else [])
    for i,j in edges:a[i,j]=a[j,i]=1
    if label=='T':a[-1,-1]=1
    return a

def main(rank):
    directory=HERE/f'rank{rank}';named=[]
    for first in range(2,rank+3):
        for tail in lower.compositions(rank+2-first):
            ns=[first]+tail;named.append({'category':'RSG','label':'Reduced sine-Gordon RSG('+', '.join(map(str,ns))+')','parameters':ns,'reference':lower.NS,'pair':lower.rsg(ns)})
    # Use the frozen Sage constructor directly for the SG family.
    sys.path.insert(0,str(HERE.parent/'catalogue/sources'))
    from rsg_constructor_snapshot import SG
    for first in range(2,rank):
        for tail in lower.compositions(rank-1-first):
            ns=[first]+tail;pair=SG(ns).t_datum()
            pair=tuple(S.Matrix([[S.sympify(str(x).replace('^','**'),locals={'z':z}) for x in row] for row in a.rows()]) for a in pair)
            named.append({'category':'SG','label':'Sine-Gordon SG('+', '.join(map(str,ns))+')','parameters':ns,'reference':lower.NS,'pair':pair})
    factors=[(a,k) for a in ('A','T','D','E') for k in range(1,rank+1) if (a not in ('D','E') or (a=='D' and k>=4) or (a=='E' and k in (6,7,8)))]
    for (a,k),(b,l) in it.combinations_with_replacement(factors,2):
        if k*l!=rank or a==b=='T':continue
        pair=((1+z*z)*S.eye(rank)-z*S.kronecker_product(adj(a,k),S.eye(l)),(1+z*z)*S.eye(rank)-z*S.kronecker_product(S.eye(k),adj(b,l)))
        named.append({'category':'Zamolodchikov','label':f'Zamolodchikov ({a}{k}, {b}{l})'+(' · tadpole extension' if 'T' in (a,b) else ''),'parameters':[a+str(k),b+str(l)],'reference':lower.THESIS,'pair':pair})
    buckets={}
    for item in named:
        lower.check_pair(*item['pair']);buckets.setdefault(canonical(item['pair']),[]).append(item)
    out={}
    for record in json.loads((directory/'base-records.json').read_text()):
        target=lower.matrices(record);matches=[]
        for item in buckets.get(canonical(target),[]):
            witness=match(item['pair'],target)
            if witness:matches.append({k:v for k,v in item.items() if k!='pair'}|{'status':'exact-polynomial-match','transformation':witness,
                                      'source':'research/catalogue/sources/rsg_constructor_snapshot.py' if item['category'] in ('RSG','SG') else 'research/higher_rank/family_notes.py'})
        if all(v==2 for v in record['datum']['delays']):
            cp,cm=(2*S.eye(rank)-a.subs(z,1) for a in target)
            assert cp==cp.T and cm==cm.T and cp*cm==cm*cp
            matches.append({'category':'Commuting Cartan','label':'Commuting Cartan / tadpole construction','status':'exact-structural-check',
                            'reference':lower.THESIS,'description':'N±(z) = z C±, with symmetric nonnegative commuting C± and positive definite 2I − C±.'})
        out[record['id']]={'identifications':matches,'status':'identified' if matches else 'not identified in the checked constructors',
          'search_scope':f'Exact polynomial matching against all RSG and SG parameter lists of rank {rank}, ADE/tadpole tensor pairs, and the commuting Cartan construction. An unmatched record is not asserted to be new.',
          'parameter_convention':'RSG and SG parameters follow generation order in the frozen constructor.'}
        print(record['id'],[m['label'] for m in matches],flush=True)
    (directory/'family-notes.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')

if __name__=='__main__':main(int(sys.argv[1]))
