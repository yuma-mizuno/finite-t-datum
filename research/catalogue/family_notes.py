"""Match named constructors by exact polynomials, including the affine gauge.

The RSG formulas are transcribed from the archived constructor source. Only existing
predecessor families are evaluated, so missing predecessor indices contribute
zero instead of invoking Python negative indexing in the upstream constructor.
All matched instances are checked as T-data, independently of their name.
"""
from pathlib import Path
import itertools as it
import json
import sys
import sympy as S

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];z=S.Symbol('z')
PAPER='https://arxiv.org/html/2301.13239v2#S1.ThmTheorem7'
NS='https://arxiv.org/abs/1212.6853'
THESIS='https://yuma-mizuno.github.io/thesis.pdf#page=36'

def matrices(r):
 d=r['datum'];n=r['rank'];n0=S.diag(*(1+z**p for p in d['delays']))
 return [n0-S.Matrix([[sum(c*z**p for c,p in e) for e in row] for row in d[key]]) for key in ('N_plus','N_minus')]

def rsg(ns):
 F=len(ns);n=lambda a:ns[a-1];nt=lambda a:n(a)-(2 if a==1 else 0)
 p=[None,1]
 for a in range(2,F+1):p.append(n(1) if a==2 else n(a-1)*p[a-1]+p[a-2])
 idx=[(a,m) for a in range(1,F+1) for m in range(1,nt(a)+1)]
 plus=S.zeros(len(idx));minus=S.zeros(len(idx))
 for i,(a,m) in enumerate(idx):
  for j,(b,k) in enumerate(idx):
   if n(1)==2 and (a,m)==(2,1) and (b,k)==(2,1):plus[i,j]+=z**p[a]
   if (a,m)==(2,1) and (b,k)==(1,1):minus[i,j]+=z**p[a]
   target=plus if a%2==0 else minus
   if m==1:
    if a>=3 and (b,k)==(a-2,nt(a-2)):target[i,j]+=z**p[a]
    elif a>=2 and b==a-1:
     h=(nt(a-1)+1-k)*p[a-1];target[i,j]+=z**h+z**(2*p[a]-h)
   other=minus if a%2==0 else plus
   if m==nt(a) and (b,k)==(a+1,1):other[i,j]+=z**p[a]
   if a==b and abs(m-k)==1:target[i,j]+=z**p[a]
 n0=S.diag(*(1+z**(2*p[a]) for a,m in idx))
 return n0-plus,n0-minus

def compositions(n):
 if n==0:yield [];return
 for first in range(1,n+1):
  for tail in compositions(n-first):yield [first]+tail

def adj(label,n):
 a=S.zeros(n)
 edges=([(0,1),(1,2),(1,3)] if label=='D' else [(i,i+1) for i in range(n-1)])
 for i,j in edges:a[i,j]=a[j,i]=1
 if label=='T':a[n-1,n-1]=1
 return a

def match(source,target):
 n=source[0].rows;sr=[S.degree(source[0][i,i],z) for i in range(n)]
 for swap in (False,True):
  pair=source[::-1] if swap else source
  for p in it.permutations(range(n)):
   arrays=[a.extract(p,p) for a in pair]
   if any(a.subs(z,1)!=b.subs(z,1) for a,b in zip(arrays,target)):continue
   lam=S.Rational(S.degree(target[0][0,0],z),sr[p[0]])
   if any(lam*sr[p[i]]!=S.degree(target[0][i,i],z) for i in range(n)):continue
   ss=S.symbols('s:'+str(n));equations=[ss[-1]]
   for a,b in zip(arrays,target):
    for i in range(n):
     for j in range(n):
      if i!=j and a[i,j]!=0:
       amin=min(e[0] for e,c in S.Poly(a[i,j],z).terms());bmin=min(e[0] for e,c in S.Poly(b[i,j],z).terms())
       equations.append(ss[i]-ss[j]-(bmin-lam*amin))
   solutions=S.linsolve(equations,ss)
   for shift in solutions:
    if any(s.free_symbols for s in shift):continue
    good=True
    for a,b in zip(arrays,target):
     for i in range(n):
      for j in range(n):
       changed=sum(c*z**(lam*e[0]+shift[i]-shift[j]) for e,c in S.Poly(a[i,j],z).terms())
       if S.expand(changed-b[i,j])!=0:good=False;break
      if not good:break
     if not good:break
    if good:return {'direction':'named constructor to displayed representative','sign_exchange':swap,'permutation':list(p),'lambda':str(lam),'shifts':[str(s) for s in shift]}
 return None

def check_pair(a,b):
 n=a.rows;r=[S.degree(a[i,i],z) for i in range(n)];zero=S.diag(*(1+z**p for p in r))
 assert all(S.expand(v)==0 for v in a*b.subs(z,1/z).T-b*a.subs(z,1/z).T)
 for i in range(n):
  for j in range(n):
   supports=[]
   for c in (zero-a,zero-b):
    terms={e[0]:v for e,v in S.Poly(c[i,j],z).terms() if v}
    assert all(v>0 and v.is_Integer and 0<p<r[i] for p,v in terms.items());supports.append(set(terms))
   assert not supports[0]&supports[1]

def main():
 rs=json.loads((HERE/'lower-rank-records.json').read_text())+[r for r in json.loads((ROOT/'docs/catalogue/catalogue.json').read_text(encoding='utf8'))['records'] if r['rank']>=3]
 named={n:[] for n in range(1,5)}
 for n in range(1,5):
  for first in range(2,n+3):
   for tail in compositions(n+2-first):
    ns=[first]+tail;pair=rsg(ns);check_pair(*pair)
    named[n].append({'label':'Reduced sine-Gordon RSG('+', '.join(map(str,ns))+')','category':'RSG','parameters':ns,
                     'reference':NS,'source':'research/catalogue/sources/rsg_constructor_snapshot.py','pair':pair})
  factors=[(a,k) for a in ('A','T') for k in range(1,n+1)]+([('D',4)] if n==4 else [])
  for (a,k),(b,l) in it.product(factors,repeat=2):
   if k*l!=n or (a=='T' and b=='T'):continue
   if (a,k)>(b,l):continue
   ap=(1+z*z)*S.eye(n)-z*S.kronecker_product(adj(a,k),S.eye(l))
   am=(1+z*z)*S.eye(n)-z*S.kronecker_product(S.eye(k),adj(b,l));check_pair(ap,am)
   named[n].append({'label':f'Zamolodchikov ({a}{k}, {b}{l})'+(' · tadpole extension' if 'T' in (a,b) else ''),
                    'category':'Zamolodchikov','parameters':[a+str(k),b+str(l)],'reference':THESIS,'pair':(ap,am)})
 # The only SG lists of rank at most four (rank = sum(n_a) + 1).
 for ns in ([2],[3],[2,1]):
  n=sum(ns)+1;ap=(1+z*z)*S.eye(n);am=ap.copy()
  for i,j in [(0,2),(1,2)]+([(2,3)] if ns==[3] else []):am[i,j]=am[j,i]=-z
  if ns==[2,1]:
   ap[3,3]=am[3,3]=1+z**4
   ap[3,0]=ap[3,1]=-z*z;ap[3,2]=-z-z**3;ap[2,3]=-z
  check_pair(ap,am)
  named[n].append({'label':'Sine-Gordon SG('+', '.join(map(str,ns))+')','category':'SG','parameters':ns,
                    'reference':NS,'source':'research/catalogue/sources/rsg_constructor_snapshot.py','pair':(ap,am)})
 sys.path.insert(0,str(ROOT))
 from audit_mizuno_finite_type_clocks import DATA
 for example in DATA:
  if example.a_plus.rows==3:
   named[3].append({'label':'Mizuno 1912.05710, Table 2, row '+example.label.split('-')[-1],
     'category':'Literature example','reference':'https://arxiv.org/abs/1912.05710',
     'pair':(example.a_plus,example.a_minus)})
 notes={}
 remarks={
  'r2-c01':('Zamolodchikov A2','The opposite pair is the Zamolodchikov A2 Y-system.'),
  'r2-c02':('RSG continued fraction 3/4','Remark 1.7 identifies the opposite pair with the reduced sine-Gordon system for 3/4 = [1,3].'),
  'r2-c03':('Mizuno’s nonstandard rank-two family','Remark 1.7 describes this system as apparently new there, and notes its earlier implicit appearance in arXiv:1912.05710, Table 2.'),
  'r2-c04':('Half of the (A2,A2) system','Remark 1.7 identifies this pair as half of the Y-system associated with (A2,A2). It also has the exact tadpole tensor description recorded below.'),
  'r2-c05':('RSG continued fraction 3/5','Remark 1.7 identifies the opposite pair after z ↦ z² with RSG for 3/5 = [1,1,2].'),
  'r2-c06':('Zamolodchikov tadpole T2','The opposite pair is the tadpole T2 Y-system, by Remark 1.7.')}
 for r in rs:
  matches=[];target=matrices(r)
  for constructor in named[r['rank']]:
   witness=match(constructor['pair'],target)
   if witness:matches.append({k:v for k,v in constructor.items() if k!='pair'}|{'status':'exact-polynomial-match','transformation':witness})
  if all(v==2 for v in r['datum']['delays']):
   cp=2*S.eye(r['rank'])-target[0].subs(z,1);cm=2*S.eye(r['rank'])-target[1].subs(z,1)
   assert cp==cp.T and cm==cm.T and cp*cm==cm*cp
   assert all((2*S.eye(r['rank'])-c)[:k,:k].det()>0 for c in (cp,cm) for k in range(1,r['rank']+1))
   matches.append({'label':'Commuting Cartan / tadpole construction','category':'Commuting Cartan',
    'status':'exact-structural-check','reference':THESIS,
    'description':'N±(z) = z C±, where C± are symmetric nonnegative integral matrices and C₊C₋ = C₋C₊. The two matrices 2I − C± have positive principal minors.'})
  note={'identifications':matches,'status':'identified' if matches or r['id'] in remarks else 'not identified in the checked constructors',
        'search_scope':'Exact affine-gauge matching against RSG lists of this rank, the rank-at-most-four SG constructors, ADE/tadpole tensor pairs and the seven rank-three examples in Mizuno’s Table 2; also checked for the commuting Cartan construction. Absence of a match is not a claim of novelty.'}
  if any(m['category']=='RSG' for m in matches):
   note['parameter_convention']='RSG parameters are in generation order in the saved constructor. This reverses the continued-fraction list in Remark 1.7: RSG(3, 1) corresponds to 3/4 = [1,3], and RSG(2, 1, 1) to 3/5 = [1,1,2].'
  if r['id'] in remarks:note['literature_note']={'title':remarks[r['id']][0],'text':remarks[r['id']][1],'reference':PAPER}
  notes[r['id']]=note;print(r['id'],[m['label'] for m in matches],flush=True)
 (HERE/'family-notes.json').write_text(json.dumps(notes,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
