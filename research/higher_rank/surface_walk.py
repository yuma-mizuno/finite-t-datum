"""Find and certify an oriented sphere triangulation directly from a quiver."""
import itertools as it
import json
from pathlib import Path
import random
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def triangulation(b):
    n=len(b)
    if np.max(abs(b))>1 or any(sum(abs(row))!=4 for row in b):return None
    arrows={(i,j) for i in range(n) for j in range(n) if b[i,j]>0}
    cycles=[]
    for i,j,k in it.combinations(range(n),3):
        for cycle in ((i,j,k),(i,k,j)):
            edges={(cycle[t],cycle[(t+1)%3]) for t in range(3)}
            if edges<=arrows:cycles.append((cycle,edges))
    by_edge={e:[i for i,(_,es) in enumerate(cycles) if e in es] for e in arrows}
    def glue(chosen):
        faces=[cycles[i][0] for i in chosen];parent=list(range(3*len(faces)));occurrences={i:[] for i in range(n)}
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        def merge(i,j):parent[root(i)]=root(j)
        for t,face in enumerate(faces):
            for k,edge in enumerate(face):occurrences[edge].append((3*t+k,3*t+(k+1)%3))
        if any(len(v)!=2 for v in occurrences.values()):return None
        for pair in occurrences.values():
            (a,c),(d,e)=pair;merge(a,e);merge(c,d)
        roots=sorted({root(i) for i in range(len(parent))});labels={v:i for i,v in enumerate(roots)}
        oriented=[[labels[root(3*t+k)] for k in range(3)] for t in range(len(faces))]
        if any(len(set(f))!=3 for f in oriented):return None
        edges=[tuple(sorted((labels[root(v[0][0])],labels[root(v[0][1])]))) for v in occurrences.values()]
        if len(set(edges))!=n or len(roots)-n+len(faces)!=2:return None
        ordered=sorted(edges);permutation={str(i):ordered.index(e) for i,e in enumerate(edges)}
        target=np.zeros((n,n),dtype=np.int16)
        for face in oriented:
            vertices=[ordered.index(tuple(sorted((face[k],face[(k+1)%3])))) for k in range(3)]
            for k in range(3):target[vertices[k],vertices[(k+1)%3]]+=1;target[vertices[(k+1)%3],vertices[k]]-=1
        assert all(b[i,j]==target[permutation[str(i)],permutation[str(j)]] for i in range(n) for j in range(n))
        return {'genus':0,'punctures':len(roots),'oriented_faces':oriented,'edge_order':ordered,'B':target.tolist()},permutation
    def visit(remaining,chosen):
        if not remaining:return glue(chosen)
        edge=min(remaining,key=lambda e:sum(cycles[i][1]<=remaining for i in by_edge[e]))
        for i in by_edge[edge]:
            if cycles[i][1]<=remaining:
                result=visit(remaining-cycles[i][1],chosen+[i])
                if result:return result
        return None
    return visit(arrows,[])


def mutate(b,k):
    out=b+np.outer(np.maximum(b[:,k],0),np.maximum(b[k,:],0))-np.outer(np.maximum(-b[:,k],0),np.maximum(-b[k,:],0))
    out[k,:]=-b[k,:];out[:,k]=-b[:,k];return out


def main(rank,seconds):
    directory=HERE/f'rank{rank}';path=directory/'quiver-data.json';out=json.loads(path.read_text())
    for record in json.loads((directory/'base-records.json').read_text()):
        q=out[record['id']]
        if q['status'].startswith('certified'):continue
        start=time.monotonic();rng=random.Random(record['constant_id']);initial=np.array(record['slice']['B'],dtype=np.int16);n=len(initial);found=None;attempts=0
        while time.monotonic()-start<seconds and found is None:
            b=initial.copy();word=[];previous=-1
            for step in range(100):
                attempts+=1;result=triangulation(b)
                if result:
                    surface,permutation=result
                    # This verifier additionally checks the manifold links and
                    # connectedness; Euler characteristic alone is insufficient.
                    sys.path.insert(0,str(HERE.parent/'catalogue'))
                    from verify_enrichment import verify_sphere
                    verify_sphere(surface)
                    q.update({'status':'certified-surface','label':f'Sphere with {surface["punctures"]} punctures','mutation_finite':True,'cluster_finite':False,
                              'certificate':{'mutation_path':word,'target_B':b.tolist(),'standard_type':None,'surface_triangulation':surface,'isomorphism_to_standard':permutation}})
                    found=True;break
                choices=[]
                for k in range(n):
                    if k==previous:continue
                    other=mutate(b,k)
                    if np.max(abs(other))>2:
                        i,j=next((i,j) for i in range(n) for j in range(i+1,n) if abs(other[i,j])>2)
                        q.update({'status':'certified-mutation-infinite','label':'Mutation-infinite','mutation_finite':False,'cluster_finite':False,
                                  'certificate':{'mutation_path':word+[k],'target_B':other.tolist(),'edge_with_at_least_three_arrows':[i,j,int(abs(other[i,j]))]}})
                        found=True;break
                    score=int(abs(np.abs(other).sum(axis=0)-4).sum())+6*int((abs(other)>1).sum())
                    choices.append((score+rng.random()*8,k,other))
                if found or not choices:break
                _,previous,b=min(choices,key=lambda x:x[0]);word.append(previous)
                if time.monotonic()-start>=seconds:break
        path.write_text(json.dumps(out,indent=2)+'\n')
        print(record['id'],q['status'],'attempts',attempts,'seconds',time.monotonic()-start,flush=True)


if __name__=='__main__':main(int(sys.argv[1]),float(sys.argv[2]) if len(sys.argv)>2 else 120)
