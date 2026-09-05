"""Deterministic heuristic search; success is checked by an exact path replay."""
import json
from pathlib import Path
import random
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from quiver_types import tree_type,mutate

def walk(initial,seed,seconds=35):
    rng=random.Random(seed);start=time.monotonic();n=len(initial);initial=np.array(initial,dtype=np.int16)
    while time.monotonic()-start<seconds:
        b=initial.copy();word=[];previous=-1
        for step in range(600):
            if int(abs(b).sum())==2*(n-1):
                kind=tree_type(b.tolist())
                if kind:return kind,word,b.tolist()
            options=[]
            for k in range(n):
                if k==previous:continue
                other=b+np.outer(np.maximum(b[:,k],0),np.maximum(b[k,:],0))-np.outer(np.maximum(-b[:,k],0),np.maximum(-b[k,:],0))
                other[k,:]=-b[k,:];other[:,k]=-b[:,k]
                if np.max(abs(other))>2:continue
                degrees=np.count_nonzero(other,axis=0)
                score=10*(int(abs(other).sum())//2)+int(np.maximum(degrees-2,0).sum())
                options.append((score+rng.random()*3,k,other))
            if not options:break
            _,previous,b=min(options,key=lambda x:x[0]);word.append(previous)
            if step%100==0 and time.monotonic()-start>=seconds:return None
    return None

def main(rank):
    directory=HERE/f'rank{rank}';target=directory/'quiver-data.json';out=json.loads(target.read_text())
    for r in json.loads((directory/'base-records.json').read_text()):
        q=out[r['id']]
        if q['status']!='recognized':continue
        result=walk(r['slice']['B'],r['class_number'])
        if result:
            kind,path,b=result;current=tuple(map(tuple,r['slice']['B']))
            # Erase closed subwalks while preserving the exact endpoint.
            states=[current];positions={current:0};short=[]
            for k in path:
                current=mutate(current,k)
                if current in positions:
                    pos=positions[current];states=states[:pos+1];short=short[:pos];positions={x:i for i,x in enumerate(states)}
                else:short.append(k);positions[current]=len(states);states.append(current)
            assert current==tuple(map(tuple,b)) and tree_type(b)==kind
            q.update({'status':'certified-dynkin','label':kind,'mutation_finite':True,'cluster_finite':True,
                      'certificate':{'mutation_path':short,'target_B':b,'target_graph':kind}})
        print(r['id'],q['status'],len(q.get('certificate',{}).get('mutation_path',[])),flush=True)
        target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main(int(sys.argv[1]))
