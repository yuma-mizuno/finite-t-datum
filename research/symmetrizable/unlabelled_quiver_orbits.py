"""Finite mutation-orbit certificates with explicit relabelling on every transition."""
from collections import Counter,deque
import gzip,hashlib,itertools,json,math,sys,time
from pathlib import Path
from sage.all import DiGraph
from quiver_types import mutate,heavy,tree_type,certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def canonical(b,d):
    n=len(b);g=DiGraph();g.add_vertices(range(n))
    g.add_edges((i,j,(int(b[i][j]),int(-b[j][i]))) for i in range(n) for j in range(n) if b[i][j]>0)
    partition=[[i for i in range(n) if d[i]==value] for value in sorted(set(d))]
    graph,p=g.canonical_label(partition=partition,edge_labels=True,certificate=True)
    p=tuple(int(p[i]) for i in range(n));assert sorted(p)==list(range(n))
    inverse=tuple(p.index(i) for i in range(n));result=tuple(tuple(b[i][j] for j in inverse) for i in inverse)
    weights=tuple(d[i] for i in inverse)
    assert weights==tuple(sorted(d))
    for i,j,label in graph.edges(sort=False):assert label==(result[i][j],-result[j][i])
    return result,p

def replay_transition(b,k,target,p):
    c=mutate(b,k);n=len(b)
    return sorted(p)==list(range(n)) and all(c[i][j]==target[p[i]][p[j]] for i in range(n) for j in range(n))

def verify_orbit(r,q):
    cert=q['certificate'];path=ROOT/cert['orbit_archive'];assert hashlib.sha256(path.read_bytes()).hexdigest()==cert['archive_sha256']
    payload=json.loads(gzip.decompress(path.read_bytes()));states=[tuple(map(tuple,b)) for b in payload['states']];d=payload['symmetrizer'];n=len(d)
    assert len(states)==cert['canonical_orbit_size'] and len(set(states))==len(states)
    p=cert['initial_to_canonical'];initial=r['slice']['B'];root=cert['initial_state']
    assert sorted(p)==list(range(n)) and all(initial[i][j]==states[root][p[i]][p[j]] for i in range(n) for j in range(n))
    assert all(r['slice']['symmetrizer'][i]==d[p[i]] for i in range(n))
    assert all(len(row)==n for row in payload['transitions']) and len(payload['transitions'])==len(states)
    for a,(b,transitions) in enumerate(zip(states,payload['transitions'])):
        assert not heavy(b) and tree_type(b) is None
        for k,(target,p) in enumerate(transitions):
            assert all(d[i]==d[p[i]] for i in range(n))
            assert replay_transition(b,k,states[target],p),(a,k,target)
    reached={root};queue=deque([root])
    while queue:
        for target,p in payload['transitions'][queue.popleft()]:
            if target not in reached:reached.add(target);queue.append(target)
    assert len(reached)==len(states)
    minimum=tuple(d)+min(tuple(x for row in b for x in row) for b in states)
    assert list(minimum)==payload['canonical_minimum']
    digest=hashlib.sha256(json.dumps(minimum,separators=(',',':')).encode()).hexdigest()
    assert digest==cert['canonical_class_sha256']
    assert cert['labelled_orbit_upper_bound']==len(states)*math.prod(math.factorial(c) for c in Counter(d).values())
    return True

def close_orbit(initial,d,archive):
    original=tuple(map(tuple,initial));n=len(d);root,p0=canonical(original,d);weights=tuple(sorted(d))
    states=[root];indices={root:0};queue=deque([0]);transitions=[]
    paths=[()];labels=[tuple(p0.index(i) for i in range(n))];start=time.monotonic();last=start
    while queue:
        index=queue.popleft();b=states[index];path=paths[index];lab=labels[index]
        kind=tree_type(b)
        if heavy(b) or kind:
            actual=original
            for k in path:actual=mutate(actual,k)
            assert all(b[i][j]==actual[lab[i]][lab[j]] for i in range(n) for j in range(n))
            return certificate(actual,path,kind)
        row=[]
        for k in range(n):
            other=mutate(b,k);c,p=canonical(other,weights)
            if c not in indices:
                indices[c]=len(states);states.append(c);queue.append(indices[c])
                paths.append(path+(lab[k],));labels.append(tuple(lab[p.index(i)] for i in range(n)))
            row.append([indices[c],list(p)])
        assert index==len(transitions);transitions.append(row)
        if time.monotonic()-last>20:
            print(n,'vertices:',len(states),'canonical states,',len(queue),'queued;',round(time.monotonic()-start,1),'seconds',flush=True);last=time.monotonic()
    minimum=tuple(weights)+min(tuple(x for row in b for x in row) for b in states)
    digest=hashlib.sha256(json.dumps(minimum,separators=(',',':')).encode()).hexdigest()
    payload={'symmetrizer':weights,'states':states,'transitions':transitions,'canonical_minimum':minimum,
             'convention':'Each transition stores its target index and the explicit old-to-new vertex permutation after one mutation.'}
    archive.parent.mkdir(exist_ok=True);archive.write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
    cert={'mutation_path':[],'target_B':initial,'orbit_archive':archive.relative_to(ROOT).as_posix(),'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
          'orbit_modulo_relabeling':True,'canonical_orbit_size':len(states),'labelled_orbit_upper_bound':len(states)*math.prod(math.factorial(c) for c in Counter(d).values()),
          'initial_state':0,'initial_to_canonical':list(p0),'canonical_class_sha256':digest,'complete_transition_closure':True}
    result={'status':'certified-finite-orbit','label':f'Finite valued mutation class ({len(states)} diagrams)',
            'vertices':n,'mutation_finite':True,'cluster_finite':False,'certificate':cert,
            'reference':'https://arxiv.org/abs/1006.4276',
            'scope':'One connected valued mutation class. All transitions close in the stored finite set after explicitly certified vertex relabellings.',
            'search_note':'Complete mutation orbit up to weight-preserving relabelling, with a permutation certificate on every transition. No finite Dynkin representative occurs.'}
    assert verify_orbit({'slice':{'B':initial,'symmetrizer':d}},result)
    result['certificate_replayed']=True;return result

def main(rank):
    directory=HERE/f'rank{rank}';cache=json.loads((directory/'pending-quivers.json').read_text());file=directory/'pending-orbit-certificates.json'
    out=json.loads(file.read_text()) if file.exists() else {};start=time.monotonic()
    for digest,item in cache.items():
        if item['result']['status']!='unidentified' or digest in out:continue
        b,d=item['initial_B'],item['symmetrizer'];print('Closing',digest[:12],len(b),'vertices',flush=True)
        q=close_orbit(b,d,directory/'quiver_orbits'/('orbit-'+digest[:16]+'.json.gz'));out[digest]=q
        temporary=file.with_suffix('.tmp');temporary.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n');temporary.replace(file)
        print(digest[:12],q['label'],q['status'],flush=True)
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out_file:
        out_file.write(json.dumps({'stage':'complete-quiver-orbits-up-to-relabeling','enumeration_complete':False,'wall_seconds':time.monotonic()-start,'quivers_checked':len(out)})+'\n')
def final_records(rank):
    directory=HERE/f'rank{rank}';path=directory/'quiver-data.json';data=json.loads(path.read_text());start=time.monotonic()
    for r in json.loads((directory/'base-records.json').read_text()):
        if data[r['id']]['status']!='unidentified':continue
        data[r['id']]=close_orbit(r['slice']['B'],r['slice']['symmetrizer'],directory/'quiver_orbits'/(r['id']+'.json.gz'))
        from quiver_types import verify
        data[r['id']]['certificate_replayed']=verify(r,data[r['id']])
        temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n');temporary.replace(path)
        print(r['id'],data[r['id']]['label'],flush=True)
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:
        out.write(json.dumps({'stage':'final-record-quiver-orbits-up-to-relabeling','wall_seconds':time.monotonic()-start,'certified':len(data)})+'\n')
if __name__=='__main__':
    (final_records if '--final' in sys.argv else main)(int(sys.argv[1]))
