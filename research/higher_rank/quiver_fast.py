"""Exact quiver witnesses with vectorized mutation and bounded best-first search."""
import heapq
import itertools as it
import json
from pathlib import Path
import signal
import sys
import numpy as np
from sage.all import ClusterQuiver,matrix,ZZ
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
from quiver_types import tree_type

def dynkin(initial,limit=12000):
    b=np.array(initial,dtype=np.int16);counter=it.count();queue=[(int(abs(b).sum())//2,0,next(counter),b,())];seen={b.tobytes()}
    for _ in range(limit):
        if not queue:break
        score,depth,_,b,word=heapq.heappop(queue)
        if score==len(b)-1:
            kind=tree_type(b.tolist())
            if kind:return kind,list(word),b.tolist()
        for k in range(len(b)):
            if word and word[-1]==k:continue
            other=b+np.outer(np.maximum(b[:,k],0),np.maximum(b[k,:],0))-np.outer(np.maximum(-b[:,k],0),np.maximum(-b[k,:],0))
            other[k,:]=-b[k,:];other[:,k]=-b[:,k]
            if np.max(abs(other))>2:continue
            key=other.tobytes()
            if key in seen:continue
            seen.add(key);heapq.heappush(queue,(int(abs(other).sum())//2,-len(word)-1,next(counter),other,word+(k,)))
    return None

def main(rank):
    directory=HERE/f'rank{rank}';path=directory/'quiver-data.json';results=json.loads(path.read_text());records=json.loads((directory/'base-records.json').read_text())
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError()))
    for r in records:
        q=results[r['id']]
        if q['status'].startswith('certified'):continue
        try:
            signal.alarm(20);cert=dynkin(r['slice']['B'])
            if cert:
                kind,word,target=cert;q.update({'status':'certified-dynkin','label':kind,'mutation_finite':True,'cluster_finite':True,
                                              'certificate':{'mutation_path':word,'target_B':target,'target_graph':kind}})
        except TimeoutError:pass
        finally:signal.alarm(0)
        path.write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
        print(r['id'],q['status'],q['label'],flush=True)

if __name__=='__main__':main(int(sys.argv[1]))
