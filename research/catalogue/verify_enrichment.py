"""Independently replay quiver certificates and the logarithmic Jacobian.

Run with SageMath. The Jacobian calculation uses thesis (3.4.2), independently
of the rank-by-rank determinant expansion used to produce spectral-data.json.
"""
from pathlib import Path
import json
from sage.all import matrix,vector,identity_matrix,QQ,ZZ,ClusterQuiver
from exponents import records,solution_box,RI,CI

HERE=Path(__file__).resolve().parent

def mutation(b,k):
 n=b.nrows();return matrix(ZZ,[[(-b[i,j] if i==k or j==k else b[i,j]+max(b[i,k],0)*max(b[k,j],0)-max(-b[i,k],0)*max(-b[k,j],0)) for j in range(n)] for i in range(n)])

def main():
 spectral=json.loads((HERE/'spectral-data.json').read_text());quivers=json.loads((HERE/'quiver-data.json').read_text())
 report={}
 for r in records():
  ap=matrix(QQ,r['datum']['A_plus_1']);am=matrix(QQ,r['datum']['A_minus_1']);info=spectral[r['id']]
  K=matrix(QQ,info['matrix_ratios']['A_plus_inverse_A_minus']);Ki=matrix(QQ,info['matrix_ratios']['A_minus_inverse_A_plus'])
  assert ap*K==am and am*Ki==ap and K*Ki==identity_matrix(QQ,r['rank'])
  q=quivers[r['id']];cert=q['certificate'];b=matrix(ZZ,r['slice']['B'])
  for k in cert['mutation_path']:b=mutation(b,k)
  assert b==matrix(ZZ,cert['target_B'])
  if q['status']=='certified-dynkin':
   assert sum(abs(x) for x in b.list())==2*(b.nrows()-1)
   cartan=matrix(ZZ,[[2 if i==j else -abs(b[i,j]) for j in range(b.nrows())] for i in range(b.nrows())])
   assert all(cartan[:i,:i].det()>0 for i in range(1,b.nrows()+1))
  elif q['status']=='certified-mutation-infinite':
   i,j,weight=cert['edge_with_at_least_three_arrows'];assert b.nrows()>=3 and abs(b[i,j])==weight and weight>2
  else:
   if q['status']=='certified-surface':
    surface=cert['surface_triangulation'];edges=[tuple(e) for e in surface['edge_order']];target=matrix(ZZ,len(edges))
    for face in surface['oriented_faces']:
     indices=[edges.index(tuple(sorted((face[i],face[(i+1)%3])))) for i in range(3)]
     for i in range(3):target[indices[i],indices[(i+1)%3]]+=1;target[indices[(i+1)%3],indices[i]]-=1
    assert target==matrix(ZZ,surface['B'])
   else:target=ClusterQuiver(cert['standard_type']).b_matrix()
   p=[cert['isomorphism_to_standard'][str(i)] for i in range(b.nrows())]
   assert sorted(p)==list(range(b.nrows())) and all(b[i,j]==target[p[i],p[j]] for i in range(b.nrows()) for j in range(b.nrows()))
  f,_=solution_box(K);B=matrix(ZZ,r['exchange']['B']);n=B.nrows();J=identity_matrix(RI,n)
  for k in r['exchange']['mutation_vertices']:
   species=r['exchange']['vertices'][k][0];L=identity_matrix(RI,n);L[k,k]=-1
   for i in range(n):
    if i!=k:L[i,k]=max(B[k,i],0)*(1-f[species])+max(-B[k,i],0)*f[species]
   J=L*J;B=mutation(B,k)
  perm=r['exchange']['relabel_old_to_new'];inv=[perm.index(i) for i in range(n)]
  J=matrix(RI,[J.row(i) for i in inv]);B=matrix(ZZ,[[B[i,j] for j in inv] for i in inv]);assert B==matrix(ZZ,r['exchange']['B'])
  component={0}
  while True:
   more=component|{j for i in component for j in range(n) if B[i,j]}
   if more==component:break
   component=more
  indices=sorted(component);J=(J**r['slice']['components']).matrix_from_rows_and_columns(indices,indices)
  N=info['exponents']['root_order'];root=CI.zeta(N);power=identity_matrix(RI,len(indices));traces=[]
  for k in range(N):traces.append(power.trace());power=power*J
  actual=[]
  for m in range(N):
   v=sum(CI(traces[k])*root**(-m*k) for k in range(N))/N
   lo=v.real().lower().ceil();hi=v.real().upper().floor();assert lo==hi and v.imag().contains_zero()
   actual.extend([m]*int(lo))
  assert actual==info['exponents']['values'],r['id']
  report[r['id']]={'matrix_ratios':'exact','mutation_class_path':'exact replay','jacobian_spectrum':'independent interval-certified match'}
  print(r['id'],'verified',flush=True)
 (HERE/'enrichment-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
 print('All 61 enrichment records independently verified.')

if __name__=='__main__':main()
