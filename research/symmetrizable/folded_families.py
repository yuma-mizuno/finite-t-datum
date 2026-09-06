"""Certify equitable polynomial folds of named identity-symmetrizer records."""
import importlib.util,itertools,json,math,sys
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('weighted_family_constructors',HERE/'family_notes.py')
constructors=importlib.util.module_from_spec(spec);spec.loader.exec_module(constructors)
old,weighted_check,z=constructors.old,constructors.weighted_check,constructors.z

def partitions(n):
    def extend(i,blocks):
        if i==n:yield tuple(tuple(b) for b in blocks);return
        for b in blocks:
            b.append(i);yield from extend(i+1,blocks);b.pop()
        blocks.append([i]);yield from extend(i+1,blocks);blocks.pop()
    yield from extend(0,[])

def poly_sum(entries):
    result={}
    for entry in entries:
        for c,p in entry:result[p]=result.get(p,0)+c
    return tuple((c,p) for p,c in sorted(result.items()) if c)

def canonical(d,pair):
    n=len(d);target=sorted(d)
    return min(tuple(target)+tuple(sum(c for c,e in a[i][j]) for a in matrices for i in perm for j in perm)
               for perm in itertools.permutations(range(n)) if [d[i] for i in perm]==target for matrices in (pair,pair[::-1]))

def matrices(delays,pair):
    return tuple(S.diag(*(1+z**p for p in delays))-S.Matrix([[sum(c*z**p for c,p in entry) for entry in row] for row in a]) for a in pair)

def verify_relation(record,note,parents):
    fold=note['folding'];parent=parents[fold['parent_record']];n=parent['rank'];k=record['rank']
    blocks=fold['partition_zero_based'];assert sorted(i for b in blocks for i in b)==list(range(n))
    shift=list(map(S.Rational,fold['parent_centering_shifts']));scale=fold['parent_integer_dilation']
    original=old.lower.matrices(parent)
    centered=tuple(S.Matrix(n,n,lambda i,j:sum(c*z**(scale*(e[0]+shift[i]-shift[j])) for e,c in S.Poly(a[i,j],z).terms())) for a in original)
    weighted_check(centered,[1]*n)
    E=S.zeros(n,k)
    for a,block in enumerate(blocks):
        for i in block:E[i,a]=1
    sizes=E.T*E;assert list(sizes.diagonal())==fold['block_sizes']
    Q=tuple(sizes.inv()*E.T*a*E for a in centered)
    assert all(S.expand(x)==0 for a,q in zip(centered,Q) for x in list(a*E-E*q)+list(a.T*E-E*sizes.inv()*q.T*sizes))
    d=fold['quotient_symmetrizer'];assert S.diag(*d)==math.lcm(*fold['block_sizes'])*sizes.inv()
    source=Q
    if fold['langlands_dual']:
        D=S.diag(*d);source=tuple(D.inv()*a*D for a in Q)
    datum=fold['folded_datum'];assert all(a==b for a,b in zip(source,matrices(datum['delays'],[datum['N_plus'],datum['N_minus']])))
    weighted_check(source,datum['symmetrizer'])
    w=note['transformation'];perm=w['permutation'];lam=S.Rational(w['lambda']);gauge=list(map(S.Rational,w['shifts']))
    if w['sign_exchange']:source=source[::-1]
    target=old.lower.matrices(record)
    for a,b in zip(source,target):
        a=a.extract(perm,perm)
        changed=S.Matrix(k,k,lambda i,j:sum(c*z**(lam*e[0]+gauge[i]-gauge[j]) for e,c in S.Poly(a[i,j],z).terms()))
        assert all(S.expand(x)==0 for x in changed-b)
    assert [datum['symmetrizer'][i] for i in perm]==record['datum']['symmetrizer']
    return True

def main():
    all_records=json.loads((ROOT/'docs/catalogue/catalogue.json').read_text())['records']
    targets={};notes={}
    for path in sorted(HERE.glob('rank*/base-records.json')):
        n=int(path.parent.name[4:]);notes[n]=json.loads((path.parent/'family-notes.json').read_text())
        for r in json.loads(path.read_text()):
            pair=[r['datum'][s] for s in ('N_plus','N_minus')]
            targets.setdefault(canonical(r['datum']['symmetrizer'],pair),[]).append(r)
        for note in notes[n].values():note['identifications']=[m for m in note['identifications'] if m.get('category')!='Fold']
    seen=set();relations=[]
    for parent in all_records:
        names=[m for m in parent['notes']['family']['identifications'] if m['category'] in ('SG','RSG')]
        if parent['scope']['symmetrizer']!='identity' or not names or parent['rank']<3:continue
        n=parent['rank'];shift=list(map(S.Rational,parent['exponents']['cartan_like_centering_shifts']))
        scale=int(S.ilcm(*[x.q for x in shift]));delays=[scale*x for x in parent['datum']['delays']]
        pair=[[[tuple((c,int(scale*(S.Rational(p)+shift[i]-shift[j]))) for c,p in entry) for j,entry in enumerate(row)] for i,row in enumerate(parent['datum'][s])] for s in ('N_plus','N_minus')]
        for blocks in partitions(n):
            k=len(blocks)
            if k<2 or k==n or k not in notes or len({len(b) for b in blocks})==1:continue
            if any(len({delays[i] for i in b})>1 for b in blocks):continue
            quotient=[];valid=True
            for a in pair:
                rows=[]
                for left in blocks:
                    row=[]
                    for right in blocks:
                        sums={poly_sum(a[i][j] for j in right) for i in left}
                        columns={poly_sum(a[i][j] for i in left) for j in right}
                        if len(sums)!=1 or len(columns)!=1:valid=False;break
                        row.append(next(iter(sums)))
                    if not valid:break
                    rows.append(row)
                if not valid:break
                quotient.append(rows)
            if not valid:continue
            multiple=math.lcm(*(len(b) for b in blocks));d=[multiple//len(b) for b in blocks]
            qdelays=[delays[b[0]] for b in blocks]
            for dual in (False,True):
                dd=d if not dual else [math.lcm(*d)//x for x in d]
                aa=quotient if not dual else [[[[ (c*d[j]//d[i],p) for c,p in entry] for j,entry in enumerate(row)] for i,row in enumerate(a)] for a in quotient]
                if dual:assert all(c*d[j]%d[i]==0 for a in quotient for i,row in enumerate(a) for j,entry in enumerate(row) for c,p in entry)
                candidates=targets.get(canonical(dd,aa),[])
                if not candidates:continue
                source=matrices(qdelays,aa)
                try:weighted_check(source,dd)
                except AssertionError:continue # A quotient need not retain opposite-sign disjointness.
                for r in candidates:
                    identity=(r['id'],parent['id'],dual)
                    if identity in seen:continue
                    witness=old.match(source,old.lower.matrices(r))
                    if not witness or [dd[i] for i in witness['permutation']]!=r['datum']['symmetrizer']:continue
                    seen.add(identity)
                    fold={'parent_record':parent['id'],'parent_centering_shifts':list(map(str,shift)),'parent_integer_dilation':scale,
                          'partition_zero_based':[list(b) for b in blocks],'block_sizes':[len(b) for b in blocks],
                          'quotient_symmetrizer':d,'langlands_dual':dual,
                          'folded_datum':{'delays':qdelays,'symmetrizer':dd,'N_plus':aa[0],'N_minus':aa[1]}}
                    for named in names:
                        note={'category':'Fold','label':'Fold of '+named['label']+(' · Langlands dual' if dual else ''),
                              'status':'exact-fold-and-polynomial-match','reference':named['reference'],'source':'research/symmetrizable/folded_families.py',
                              'description':'An exact equitable polynomial fold of '+parent['id']+'. Row and column equitability, dual integrality, strict support, disjoint signs and the weighted symplectic identity were checked. This records a structural relation to the named family; it is not asserted to be a published name for the folded datum.',
                              'transformation':witness,'folding':fold}
                        notes[k][r['id']]['identifications'].append(note)
                    relations.append({'record':r['id'],'parent':parent['id'],'dual':dual,'blocks':blocks})
    for n,table in notes.items():
        for item in table.values():
            item['status']='identified' if item['identifications'] else 'not identified in the checked constructors'
            scope=item['search_scope'].split(' Equitable folds')[0]
            item['search_scope']=scope+' Equitable folds and their Langlands duals of the named SG/RSG identity representatives through rank six were also checked after their recorded Cartan-like centering. One witness per parent and dual choice is retained; this bounded parent library does not exhaust all possible unfoldings.'
        (HERE/f'rank{n}/family-notes.json').write_text(json.dumps(table,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    (HERE/'fold-relations.json').write_text(json.dumps(relations,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Certified',len(relations),'fold relations on',len({r['record'] for r in relations}),'weighted records',flush=True)
if __name__=='__main__':main()
