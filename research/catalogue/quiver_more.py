"""Finish exceptional recognitions with explicit mutation paths."""
from pathlib import Path
import json
from sage.all import ClusterQuiver,matrix,ZZ
from quiver_types import records,dynkin_path,mutate
HERE=Path(__file__).resolve().parent

def main():
 out=json.loads((HERE/'quiver-data.json').read_text())
 for r in records:
  entry=out[r['id']]
  if entry['status'] not in ('unidentified','recognized'):continue
  cert=dynkin_path(r['slice']['B'],limit=30000) if entry['status']=='unidentified' else None
  if cert:
   kind,path,target=cert;entry.update({'status':'certified-dynkin','label':kind,'mutation_finite':True,'cluster_finite':True,
    'certificate':{'mutation_path':path,'target_B':target,'target_graph':kind}})
  else:
   q=ClusterQuiver(matrix(ZZ,r['slice']['B']));n=q.n()
   targets={};types={};surfaces={}
   for rank in (6,7,8):
    if n==rank+2:
     target=ClusterQuiver(['E',rank,[1,1]])
     targets[f'E{rank}^(1,1)']=target.digraph()
     types[f'E{rank}^(1,1)']=['E',rank,[1,1]]
   if n in (6,7):
    targets[f'X{n}']=ClusterQuiver(['X',n]).digraph();types[f'X{n}']=['X',n]
   if n>=5:
    targets[f'D{n-1}^(1)']=ClusterQuiver(['D',n-1,1]).digraph();types[f'D{n-1}^(1)']=['D',n-1,1]
   if n in (7,8,9):
    targets[f'E{n-1}^(1)']=ClusterQuiver(['E',n-1,1]).digraph();types[f'E{n-1}^(1)']=['E',n-1,1]
   if n>=6 and n%3==0:
    punctures=n//3+2
    equator=punctures-2
    faces=([(0,1,2),(0,3,1),(0,2,3),(1,3,2)] if n==6 else
           [(0,2+i,2+(i+1)%equator) for i in range(equator)]+
           [(1,2+(i+1)%equator,2+i) for i in range(equator)])
    edges=sorted({tuple(sorted((f[i],f[(i+1)%3]))) for f in faces for i in range(3)})
    sb=matrix(ZZ,n)
    for face in faces:
     triangle=[edges.index(tuple(sorted((face[i],face[(i+1)%3])))) for i in range(3)]
     for i in range(3):sb[triangle[i],triangle[(i+1)%3]]+=1;sb[triangle[(i+1)%3],triangle[i]]-=1
    label=f'Sphere with {punctures} punctures'
    targets[label]=ClusterQuiver(sb).digraph()
    surfaces[label]={'genus':0,'punctures':punctures,'oriented_faces':faces,'edge_order':edges,'B':[[int(v) for v in row] for row in sb.rows()]}
   for i,(other,path) in enumerate(q.mutation_class_iter(return_paths=True)):
    for label,target in targets.items():
     equivalent,iso=other.digraph().is_isomorphic(target,edge_labels=True,certificate=True)
     if equivalent:
      b=tuple(map(tuple,r['slice']['B']))
      for v in path:b=mutate(b,int(v))
      replay=ClusterQuiver(matrix(ZZ,b)).digraph()
      eq,iso=replay.is_isomorphic(target,edge_labels=True,certificate=True);assert eq
      entry.update({'status':'certified-surface' if label in surfaces else 'certified-exceptional','label':label,'mutation_finite':True,'cluster_finite':False,
       'certificate':{'mutation_path':[int(v) for v in path],'target_B':[list(row) for row in b],
                      'standard_type':types.get(label),'surface_triangulation':surfaces.get(label),'isomorphism_to_standard':{str(k):int(v) for k,v in iso.items()}}})
      break
    if entry['status'] in ('certified-exceptional','certified-surface') or i>=30000:break
  print(r['id'],entry['label'],entry['status'],flush=True)
  (HERE/'quiver-data.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
