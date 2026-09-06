"""Exact symmetrizable Cartan tensor matches, retaining the diagonal weights."""
import importlib.util,itertools,json,math,sys
from pathlib import Path
import sympy as S
from sage.all import CartanMatrix
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('identity_family_notes',HERE.parent/'higher_rank/family_notes.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
z=old.z
THESIS='https://yuma-mizuno.github.io/thesis.pdf#page=37'

def diagonal(a):
    n=a.rows;values=[None]*n;values[0]=S.Rational(1)
    for _ in range(n):
        for i in range(n):
            if values[i] is None:continue
            for j in range(n):
                if i!=j and a[i,j]:
                    value=values[i]*a[j,i]/a[i,j]
                    assert values[j] is None or values[j]==value
                    values[j]=value
    assert all(v is not None and v>0 for v in values)
    scale=S.ilcm(*[v.q for v in values]) if n>1 else values[0].q
    result=[int(v*scale) for v in values];g=math.gcd(*result)
    return [v//g for v in result]

def factors(rank):
    out=[]
    for n in range(1,rank+1):
        for a in 'ABCDEFGT':
            if a=='T':c=2*S.eye(n)-old.adj(a,n)
            else:
                if (a in 'BC' and n<2) or (a=='D' and n<4) or (a=='E' and n not in (6,7,8)) or (a=='F' and n!=4) or (a=='G' and n!=2):continue
                try:c=S.Matrix(CartanMatrix([a,n]).rows())
                except (ValueError,TypeError):continue
            out.append((a+str(n),c,diagonal(c)))
    return out

def weighted_check(pair,d):
    ap,am=pair;D=S.diag(*d);n=ap.rows;n0=S.diag(*(1+z**int(S.degree(ap[i,i],z)) for i in range(n)))
    ns=[n0-a for a in pair]
    assert all(S.expand(x)==0 for x in ap*D*am.subs(z,1/z).T-am*D*ap.subs(z,1/z).T)
    for i in range(n):
        for j in range(n):
            supports=[]
            for N in ns:
                entries=[(p[0],int(c)) for p,c in S.Poly(N[i,j],z).terms() if c]
                assert all(c>0 and 0<p<S.degree(ap[i,i],z) and c*d[j]%d[i]==0 for p,c in entries)
                supports.append({p for p,c in entries})
            assert not supports[0]&supports[1]

def main(rank):
    named=[]
    for (a,ca,da),(b,cb,db) in itertools.combinations_with_replacement(factors(rank),2):
        if ca.rows*cb.rows!=rank or a[0]==b[0]=='T':continue
        d=[x*y for x in da for y in db]
        if d==[1]*rank:continue
        cp=S.kronecker_product(2*S.eye(ca.rows)-ca,S.eye(cb.rows));cm=S.kronecker_product(S.eye(ca.rows),2*S.eye(cb.rows)-cb)
        pair=tuple((1+z*z)*S.eye(rank)-z*c for c in (cp,cm));weighted_check(pair,d)
        named.append({'category':'Zamolodchikov','label':f'Cartan tensor ({a}, {b})'+(' · tadpole extension' if 'T' in (a[0],b[0]) else ''),
                      'parameters':[a,b],'reference':THESIS,'pair':pair,'symmetrizer':d})
    buckets={}
    for item in named:buckets.setdefault(old.canonical(item['pair']),[]).append(item)
    directory=HERE/f'rank{rank}';out={}
    for r in json.loads((directory/'base-records.json').read_text()):
        target=old.lower.matrices(r);d=r['datum']['symmetrizer'];matches=[];weighted_check(target,d)
        for item in buckets.get(old.canonical(target),[]):
            witness=old.match(item['pair'],target)
            if witness and [item['symmetrizer'][i] for i in witness['permutation']]==d:
                matches.append({k:v for k,v in item.items() if k!='pair'}|{'status':'exact-polynomial-match','transformation':witness,'source':'research/symmetrizable/family_notes.py'})
        if all(v==2 for v in r['datum']['delays']):
            cp,cm=(2*S.eye(rank)-a.subs(z,1) for a in target);D=S.diag(*d)
            assert cp*D==D*cp.T and cm*D==D*cm.T and cp*cm==cm*cp
            assert all((2*S.eye(rank)-c)*D==((2*S.eye(rank)-c)*D).T and all((((2*S.eye(rank)-c)*D)[:i,:i]).det()>0 for i in range(1,rank+1)) for c in (cp,cm))
            matches.append({'category':'Commuting Cartan','label':'Commuting Cartan / tadpole construction','status':'exact-structural-check',
                            'reference':'https://yuma-mizuno.github.io/thesis.pdf#page=36',
                            'description':'N±(z) = z C±. The nonnegative matrices C± have disjoint support, commute, and share the displayed right symmetrizer D; (2I − C±)D is positive definite.'})
        out[r['id']]={'identifications':matches,'status':'identified' if matches else 'not identified in the checked constructors',
                     'search_scope':'Exact weighted polynomial matching against finite Cartan A, B, C, D, E, F, G and tadpole tensor pairs, and the commuting Cartan construction. The RSG and SG constructors previously checked have identity symmetrizer. An unmatched record is not asserted to be new; folded identifications have not been exhausted.',
                     'parameter_convention':'Cartan names and valuations follow Sage CartanMatrix. Source-to-record transformations retain D; Langlands duals are not identified automatically.'}
        print(r['id'],[m['label'] for m in matches],flush=True)
    (directory/'family-notes.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':
    for rank in map(int,sys.argv[1:]):main(rank)
