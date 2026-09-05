"""Recognize mutation classes, preserving replayable integer witnesses."""
from pathlib import Path
import heapq
import itertools
import json
from sage.all import ClusterQuiver, matrix, ZZ, set_random_seed

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
records=json.loads((HERE/'lower-rank-records.json').read_text())+[r for r in json.loads((ROOT/'docs/catalogue/catalogue.json').read_text(encoding='utf-8'))['records'] if r['rank']>=3]

def mutate(b,k):
 return tuple(tuple(-x if i==k or j==k else x+max(b[i][k],0)*max(b[k][j],0)-max(-b[i][k],0)*max(-b[k][j],0) for j,x in enumerate(row)) for i,row in enumerate(b))

def tree_type(b):
 n=len(b);adj=[[j for j,x in enumerate(row) if x] for row in b]
 if any(abs(x)>1 for row in b for x in row) or sum(map(len,adj))!=2*(n-1):return None
 deg=list(map(len,adj))
 if max(deg,default=0)<=2:return 'A'+str(n)
 if sorted(d for d in deg if d>2)!=[3]:return None
 branch=deg.index(3);arms=[]
 for v in adj[branch]:
  prev=branch;length=1
  while len(adj[v])==2:
   nxt=next(w for w in adj[v] if w!=prev);prev,v=v,nxt;length+=1
  arms.append(length)
 arms.sort()
 if arms[:2]==[1,1]:return 'D'+str(n)
 if arms in ([1,2,2],[1,2,3],[1,2,4]):return 'E'+str(n)
 return None

def dynkin_path(initial,limit=12000):
 initial=tuple(map(tuple,initial));queue=[];counter=itertools.count();seen={initial}
 def score(b):return sum(abs(x) for row in b for x in row)//2
 heapq.heappush(queue,(score(initial),0,next(counter),initial,()))
 for _ in range(limit):
  if not queue:break
  _,depth,_,b,word=heapq.heappop(queue);kind=tree_type(b)
  if kind:return kind,list(word),[list(row) for row in b]
  for k in range(len(b)):
   if word and word[-1]==k:continue
   other=mutate(b,k)
   if other in seen or any(abs(x)>2 for row in other for x in row):continue
   seen.add(other);heapq.heappush(queue,(score(other),depth+1,next(counter),other,word+(k,)))
 return None

def main():
 out={}
 for r in records:
  set_random_seed(8100+r['rank']*100+r['class_number']);q=ClusterQuiver(matrix(ZZ,r['slice']['B']))
  raw=q.mutation_type();entry={'software':'SageMath 10.8','recognition':str(raw),'vertices':r['slice']['vertices'],
   'scope':'one connected component of the representative exchange quiver; this class does not specify its mutation loop',
   'reference':'https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/cluster_algebra_quiver/quiver.html#sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver.mutation_type'}
  if hasattr(raw,'is_finite') and raw.is_finite():
   cert=dynkin_path(r['slice']['B']);assert cert,(r['id'],str(raw))
   kind,path,target=cert;entry.update({'status':'certified-dynkin','label':kind,'mutation_finite':True,'cluster_finite':True,
      'certificate':{'mutation_path':path,'target_B':target,'target_graph':kind}})
  else:
   result=q.is_mutation_finite(nr_of_checks=500*q.n(),return_path=True)
   if isinstance(result,tuple) and not result[0]:
    path=[int(v) for v in result[1]];b=tuple(map(tuple,r['slice']['B']))
    for k in path:b=mutate(b,k)
    heavy=next((i,j,abs(b[i][j])) for i in range(len(b)) for j in range(i+1,len(b)) if abs(b[i][j])>2)
    entry.update({'status':'certified-mutation-infinite','label':'Mutation-infinite','mutation_finite':False,'cluster_finite':False,
      'certificate':{'mutation_path':path,'target_B':[list(row) for row in b],'edge_with_at_least_three_arrows':list(heavy)}})
   else:entry.update({'status':'recognized' if not isinstance(raw,str) else 'unidentified','label':str(raw),'mutation_finite':None,'cluster_finite':False})
  out[r['id']]=entry;print(r['id'],entry['label'],entry['status'],flush=True)
  (HERE/'quiver-data.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
