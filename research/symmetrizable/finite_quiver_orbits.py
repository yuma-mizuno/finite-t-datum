"""Close the full labelled mutation orbit for the remaining small valued quivers."""
from collections import deque
import gzip,hashlib,itertools,json,sys,time
from pathlib import Path
from quiver_types import mutate,heavy,tree_type
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def main(rank):
    directory=HERE/f'rank{rank}';records=json.loads((directory/'base-records.json').read_text());path=directory/'quiver-data.json';data=json.loads(path.read_text());start=time.monotonic()
    for r in records:
        if data[r['id']]['status']!='unidentified':continue
        initial=tuple(map(tuple,r['slice']['B']));d=r['slice']['symmetrizer'];n=len(initial)
        queue=deque([initial]);words={initial:()};seen=set()
        while queue:
            b=queue.popleft();seen.add(b);assert not heavy(b) and not tree_type(b),r['id']
            for k in range(n):
                other=mutate(b,k)
                if other not in words:words[other]=words[b]+(k,);queue.append(other)
            if len(words)>2000000:raise RuntimeError('Orbit exceeds two million labelled states; no finiteness assertion made.')
        assert len(seen)==len(words)
        permutations=[p for p in itertools.permutations(range(n)) if [d[i] for i in p]==sorted(d)]
        minimum=min(tuple(sorted(d))+tuple(b[i][j] for i in p for j in p) for b in seen for p in permutations)
        identifier=hashlib.sha256(json.dumps(minimum,separators=(',',':')).encode()).hexdigest()
        target=min(seen);archive=directory/'quiver_orbits';archive.mkdir(exist_ok=True);file=archive/(r['id']+'.json.gz')
        payload={'record':r['id'],'symmetrizer':d,'initial_B':r['slice']['B'],'states':[list(map(list,b)) for b in sorted(seen)],'canonical_minimum':minimum}
        file.write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
        # Independently replay the stored archive and check every outgoing transition.
        loaded=json.loads(gzip.decompress(file.read_bytes()));states={tuple(map(tuple,b)) for b in loaded['states']}
        assert initial in states and states==seen
        assert all(mutate(b,k) in states for b in states for k in range(n))
        cert={'mutation_path':list(words[target]),'target_B':list(map(list,target)),
              'orbit_archive':str(file.relative_to(ROOT)).replace('\\','/'),'archive_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),
              'labelled_orbit_size':len(seen),'canonical_class_sha256':identifier,'complete_transition_closure':True}
        data[r['id']].update({'status':'certified-finite-orbit','label':'Finite valued mutation class '+identifier[:10],
                            'mutation_finite':True,'cluster_finite':False,'certificate':cert,'certificate_replayed':True,
                            'search_note':'All labelled mutations were exhausted; the orbit contains no finite Dynkin representative. No standard name is assigned.'})
        path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
        print(r['id'],'closed orbit',len(seen),identifier,flush=True)
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as f:f.write(json.dumps({'stage':'complete-valued-mutation-orbits','wall_seconds':time.monotonic()-start})+'\n')
if __name__=='__main__':
    for rank in map(int,sys.argv[1:]):main(rank)
