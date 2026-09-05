"""Rank-one elementary classification and the six rows of Mizuno's Table 1.

Run with ordinary Python and SymPy. Existing rank-three mutation routines are
reused, but no earlier classification data are edited.
"""
from pathlib import Path
import hashlib
import itertools as it
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/rank3'))
from classify_rank3 import certificate, z

PAPER='https://arxiv.org/html/2301.13239v2'
M=lambda a:sp.Matrix(a)
ROWS=[
 (1,1,M([[1+z]]),M([[1+z]])),
 (1,2,M([[1+z*z]]),M([[1-z+z*z]])),
 (2,1,M([[1+z*z,-z],[-z,1+z*z]]),sp.diag(1+z*z,1+z*z)),
 (2,2,M([[1+z*z,-z],[-z-z**5,1+z**6]]),M([[1+z*z,0],[-z**3,1+z**6]])),
 (2,3,M([[1+z*z,-z],[-z-z**5-z**9,1+z**10]]),M([[1+z*z,0],[-z**3-z**7,1+z**10]])),
 (2,4,M([[1+z*z,-z],[-z,1+z*z]]),sp.diag(1-z+z*z,1-z+z*z)),
 (2,5,M([[1+z*z,-z],[-z-z*z,1+z**3]]),sp.diag(1-z+z*z,1+z**3)),
 (2,6,M([[1+z*z,-z],[-z,1-z+z*z]]),sp.diag(1+z*z,1+z*z)),
]

def encode(a):
 return [[[[int(c),int(e[0])] for e,c in sorted(sp.Poly(x,z).terms()) if c] for x in row] for row in a.tolist()]

def slice_data(cert):
 b,p,mut=cert['B'],cert['relabel_old_to_new'],cert['mutation_vertices']
 seen={0}
 while True:
  more=seen|{j for i in seen for j,x in enumerate(b[i]) if x}
  if more==seen:break
  seen=more
 component=sorted(seen);index={v:i for i,v in enumerate(component)}
 current=set(component);t=0
 while True:
  t+=1;current={p[i] for i in current}
  if current==seen:break
  assert not current&seen
 assert t*len(component)==len(b)
 forward=list(range(len(b)));word=[]
 for _ in range(t):
  word.extend(index[v] for v in component if forward[v] in mut)
  forward=[p[v] for v in forward]
 sb=[[b[i][j] for j in component] for i in component]
 spm=[index[forward[v]] for v in component]
 signature=[len(component)]+[x for row in sb for x in row]+word+spm
 return {'components':t,'vertices':len(component),'B':sb,'mutation_word':word,
         'relabel_old_to_new':spm,'canonical_signature':signature,
         'signature_sha256':hashlib.sha256(json.dumps(signature,separators=(',',':')).encode()).hexdigest(),
         'signature_convention':'Labelled slice encoding; the literature theorem supplies distinctness.'}

def main():
 records=[]
 for rank,k,ap,am in ROWS:
  cert=certificate(ap,am);r=cert['delays'];n0=sp.diag(*(1+z**d for d in r))
  np,nm=encode(n0-ap),encode(n0-am)
  names=[f'r{i}' for i in range(rank)];values=list(r)
  generators=[list(r)]+[[0]*rank for _ in range(rank-1)]
  for sign,a in enumerate((np,nm)):
   for i,row in enumerate(a):
    for j,entry in enumerate(row):
     unit=0
     for c,p in entry:
      for _ in range(c):
       names.append(f'p{sign}_{i}{j}_{unit}');values.append(p);unit+=1
       generators[0].append(p)
       for l in range(rank-1):generators[l+1].append(int(i==l)-int(j==l))
  equations=M(generators).nullspace()
  rref=M([list(v) for v in equations]).rref()[0].tolist() if equations else []
  p=cert['positive'];m=cert['negative']
  record={'schema_version':'2.0.0','id':f'r{rank}-c{k:02d}','rank':rank,'class_number':k,
   'constant_id':None,'scope':{'symmetrizer':'identity','leading_permutation':'identity','indecomposable':True},'index_base':0,
   'datum':{'delays':r,'N_plus':np,'N_minus':nm,'A_plus_1':[[int(x) for x in row] for row in ap.subs(z,1).tolist()],
            'A_minus_1':[[int(x) for x in row] for row in am.subs(z,1).tolist()]},
   'family':{'parameters':['lambda']+[f's{i+1}' for i in range(rank-1)],'fixed_shift':{'species':rank-1,'value':'0'},
             'dimension':rank,'rref':[[str(x) for x in row] for row in rref], 'variable_names':names,'representative_values':values,
             'coverage':{'kind':'literature' if rank==2 else 'elementary-proof',
                         'reference':PAPER+'#S1.ThmTheorem5' if rank==2 else 'research/catalogue/lower-ranks-proof.html',
                         'statement':'Theorem 1.5 and Lemmas 3.6–3.9' if rank==2 else 'Scalar positivity, parity and reciprocity'}},
   'periodicity':{'time_coordinate':'displayed representative','h_plus':p['h'],'h_minus':m['h'],'labelled_period':cert['labelled_tropical_seed_period'],
                  'positive_negative_permutation':p['negative_permutation'],'negative_negative_permutation':m['negative_permutation']},
   'exchange':{key:cert[key] for key in ['vertices','B','mutation_vertices','relabel_old_to_new']},
   'slice':slice_data(cert),
   'provenance':{'source_commit':None,'verification_kind':'literature and exact mutation check' if rank==2 else 'elementary proof and exact mutation check',
                 'sources':[],'manuscript':'research/catalogue/lower-ranks-proof.html','pdf':None,'query_path':None,'query_member':None,
                 'classification_reference':PAPER+'#S1.T1' if rank==2 else 'research/catalogue/lower-ranks-proof.html'}}
  if rank==2: assert (p['h'],m['h'])==[(3,2),(8,6),(18,10),(3,3),(5,3),(5,2)][k-1]
  records.append(record)
 (HERE/'lower-rank-records.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
 print('Generated eight lower-rank records, with exact two-sided mutation certificates.')

if __name__=='__main__':main()
