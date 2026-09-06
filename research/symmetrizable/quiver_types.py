"""Valued mutation classes with exact integer certificates."""
from functools import lru_cache
import heapq,itertools,json,random,signal,sys,time
from pathlib import Path
from sage.all import CartanMatrix,ClusterQuiver,DiGraph,QuiverMutationType,ZZ,matrix,set_random_seed
HERE=Path(__file__).resolve().parent

def mutate(b,k):
    return tuple(tuple(-x if i==k or j==k else x+max(b[i][k],0)*max(b[k][j],0)-max(-b[i][k],0)*max(-b[k][j],0) for j,x in enumerate(row)) for i,row in enumerate(b))

def graph(b,oriented=False):
    g=DiGraph();g.add_vertices(range(len(b)))
    g.add_edges((i,j,(int(x),int(-b[j][i]))) if oriented else (i,j,abs(int(x)))
                for i,row in enumerate(b) for j,x in enumerate(row) if (x>0 if oriented else x!=0) and i!=j)
    return g

@lru_cache(None)
def templates(n):
    out=[]
    for a in 'ABCDEFG':
        if (a=='D' and n<4) or (a=='E' and n not in (6,7,8)) or (a=='F' and n!=4) or (a=='G' and n!=2) or (a in 'BC' and n<2):continue
        try:c=CartanMatrix([a,n])
        except (ValueError,TypeError):continue
        out.append((a+str(n),graph(c.rows())))
    return out

def tree_type(b):
    n=len(b)
    if sum(x!=0 for row in b for x in row)!=2*(n-1):return None
    g=graph(b)
    for name,t in templates(n):
        if g.is_isomorphic(t,edge_labels=True):return name
    return None

@lru_cache(None)
def standard_templates(n):
    candidates=[[a,k,1] for a in ('BB','BC','BD','CC','CD','D','E','F','G') for k in range(1,n+1)]
    candidates += [[a,k,list(t)] for a,k in (('F',4),('G',2)) for t in ((1,1),(1,2),(2,2),(1,3),(3,3))]
    out={}
    for data in candidates:
        try:
            kind=QuiverMutationType(data)
            if not kind.is_mutation_finite() or kind.is_finite():continue
            q=ClusterQuiver(kind)
            if q.n()!=n:continue
            b=[list(map(int,row)) for row in q.b_matrix().rows()]
            out[str(kind)]=(str(kind),graph(b,True),b)
        except (ValueError,TypeError,NotImplementedError):continue
    return list(out.values())

def heavy(b):
    if len(b)<3:return None
    return next(([i,j,abs(b[i][j]),abs(b[j][i])] for i in range(len(b)) for j in range(i+1,len(b)) if abs(b[i][j]*b[j][i])>4),None)

def certificate(b,path,kind=None):
    h=heavy(b)
    if h:
        return {'status':'certified-mutation-infinite','label':'Mutation-infinite','mutation_finite':False,'cluster_finite':False,
                'certificate':{'mutation_path':list(path),'target_B':list(map(list,b)),'edge_product_greater_than_four':h},
                'reference':'https://arxiv.org/abs/1006.4276'}
    if kind:
        return {'status':'certified-dynkin','label':kind,'mutation_finite':True,'cluster_finite':True,
                'certificate':{'mutation_path':list(path),'target_B':list(map(list,b)),'target_graph':kind}}

def search(initial,standard=None,limit=40000):
    initial=tuple(map(tuple,initial));queue=[];counter=itertools.count();seen={initial}
    def score(b):return sum(abs(b[i][j]*b[j][i]) for i in range(len(b)) for j in range(i+1,len(b)))
    heapq.heappush(queue,(score(initial),0,next(counter),initial,()))
    for _ in range(limit):
        if not queue:break
        _,depth,_,b,word=heapq.heappop(queue)
        cert=certificate(b,word,tree_type(b))
        if cert:return cert
        if standard:
            g=graph(b,True)
            for name,t,target in standard:
                iso=g.is_isomorphic(t,edge_labels=True,certificate=True)
                if iso[0]:
                    return {'status':'certified-standard','label':name,'mutation_finite':True,'cluster_finite':False,
                            'certificate':{'mutation_path':list(word),'target_B':list(map(list,b)),'standard_B':target,
                                           'standard_type':name,'target_to_standard_vertices':[int(iso[1][i]) for i in range(len(b))],
                                           'comparison':'valued directed graph isomorphism'}}
        for k in range(len(b)):
            if word and word[-1]==k:continue
            other=mutate(b,k)
            if other in seen:continue
            h=heavy(other)
            if h:return certificate(other,word+(k,))
            seen.add(other);heapq.heappush(queue,(score(other),depth+1,next(counter),other,word+(k,)))
    return None

def compute(r):
    b=r['slice']['B'];n=len(b)
    base={'vertices':n,'scope':'One connected component, with its edge valuations. The mutation class does not specify the mutation loop.',
          'reference':'https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/cluster_algebra_quiver/quiver.html'}
    initial=tuple(map(tuple,b));rng=random.Random(77100+r['rank']*100+r['class_number'])
    result=search(b,limit=2000)
    if result:return base|result
    for attempt in range(40):
        walk=initial;word=[];previous=-1
        for step in range(300*n):
            k=rng.choice([i for i in range(n) if i!=previous]);previous=k;word.append(k);walk=mutate(walk,k)
            if heavy(walk):return base|certificate(walk,word)
    raw=None;standard=None
    try:
        signal.alarm(10);set_random_seed(91800+r['rank']*100+r['class_number'])
        raw=ClusterQuiver(matrix(ZZ,b)).mutation_type();base['recognition']=str(raw)
        if hasattr(raw,'is_mutation_finite') and raw.is_mutation_finite():
            target=[list(map(int,row)) for row in ClusterQuiver(raw).b_matrix().rows()]
            standard=[(str(raw),graph(target,True),target)]
    except (TimeoutError,ValueError,TypeError):pass
    finally:signal.alarm(0)
    result=search(b,standard or standard_templates(n),limit=60000)
    return base|(result or {'status':'unidentified','label':'Mutation class not identified','mutation_finite':None,'cluster_finite':None,
                           'search_note':'No certificate was found within the recorded search; this is not a finiteness assertion.'})

def verify(r,q):
    b=tuple(map(tuple,r['slice']['B']));d=r['slice']['symmetrizer']
    assert all(b[i][j]*d[j]==-b[j][i]*d[i] for i in range(len(b)) for j in range(len(b)))
    c=q.get('certificate')
    if not c:return False
    for k in c['mutation_path']:b=mutate(b,k)
    assert list(map(list,b))==c['target_B']
    assert all(b[i][j]*d[j]==-b[j][i]*d[i] for i in range(len(b)) for j in range(len(b)))
    if q['status']=='certified-dynkin':assert tree_type(b)==q['label']
    elif q['status']=='certified-mutation-infinite':assert heavy(b)==c['edge_product_greater_than_four']
    elif q['status']=='certified-standard':assert graph(b,True).is_isomorphic(graph(c['standard_B'],True),edge_labels=True)
    elif q['status']=='certified-finite-orbit':
        if c.get('orbit_modulo_relabeling'):
            from unlabelled_quiver_orbits import verify_orbit
            return verify_orbit(r,q)
        import gzip,hashlib
        path=HERE.parents[1]/c['orbit_archive'];assert hashlib.sha256(path.read_bytes()).hexdigest()==c['archive_sha256']
        payload=json.loads(gzip.decompress(path.read_bytes()));states={tuple(map(tuple,a)) for a in payload['states']}
        assert len(states)==c['labelled_orbit_size'] and tuple(map(tuple,r['slice']['B'])) in states and b in states
        assert all(mutate(a,k) in states for a in states for k in range(len(b)))
    else:raise AssertionError(q['status'])
    return True

def main(rank):
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError()))
    directory=HERE/f'rank{rank}';path=directory/'quiver-data.json';out=json.loads(path.read_text()) if path.exists() else {};start=time.monotonic()
    cached_path=directory/'pending-quivers.json';orbit_path=directory/'pending-orbit-certificates.json'
    cached=json.loads(cached_path.read_text()) if cached_path.exists() else {}
    orbits=json.loads(orbit_path.read_text()) if orbit_path.exists() else {}
    for r in json.loads((directory/'base-records.json').read_text()):
        if r['id'] not in out or out[r['id']]['status']=='unidentified':
            import hashlib
            digest=hashlib.sha256(json.dumps({'B':r['slice']['B'],'D':r['slice']['symmetrizer']},separators=(',',':')).encode()).hexdigest()
            item=cached.get(digest)
            if item:
                assert item['initial_B']==r['slice']['B'] and item['symmetrizer']==r['slice']['symmetrizer']
                q=orbits.get(digest,item['result'])
                if q['status']!='unidentified':
                    assert verify(r,q);out[r['id']]=q
        if r['id'] not in out or out[r['id']]['status']=='unidentified':out[r['id']]=compute(r)
        q=out[r['id']];q['certificate_replayed']=verify(r,q)
        path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
        print(r['id'],q['label'],q['status'],flush=True)
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as f:f.write(json.dumps({'stage':'valued-quiver-classes','wall_seconds':time.monotonic()-start,'certified':sum(q['certificate_replayed'] for q in out.values()),'total':len(out)})+'\n')
if __name__=='__main__':
    for rank in map(int,sys.argv[1:]):main(rank)
