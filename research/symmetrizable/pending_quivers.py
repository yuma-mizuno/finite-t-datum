"""Cache exact valued-quiver witnesses from covered families before rank closure."""
import hashlib,json,signal,sys,time
from pathlib import Path
from quiver_types import compute,verify
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'catalogue'));from lower_ranks import slice_data

def key(r):return hashlib.sha256(json.dumps({'B':r['slice']['B'],'D':r['slice']['symmetrizer']},separators=(',',':')).encode()).hexdigest()

def main(rank):
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError()))
    directory=HERE/f'rank{rank}';file=directory/'pending-quivers.json';cache=json.loads(file.read_text()) if file.exists() else {}
    text=(directory/'pending-families.jsonl').read_text();text=text[:text.rfind('\n')+1]
    families={tuple(f['constant_key']):f for f in map(json.loads,text.splitlines())};start=time.monotonic();count=0
    source_sha=hashlib.sha256((HERE/'quiver_types.py').read_bytes()).hexdigest()
    for fullkey,f in families.items():
        if fullkey[:rank]==(1,)*rank or f.get('coverage_status')!='unsat':continue
        assert len(f['spaces'])==1
        cert=f['spaces'][0]['certificate'];s=slice_data(cert);component={0}
        while True:
            more=component|{j for i in component for j,x in enumerate(cert['B'][i]) if x}
            if more==component:break
            component=more
        s['symmetrizer']=[fullkey[cert['vertices'][i][0]] for i in sorted(component)]
        r={'rank':rank,'slice':s};digest=key(r);r['class_number']=int(digest[:8],16)%100000
        if digest in cache:continue
        q=compute(r);q['certificate_replayed']=verify(r,q)
        cache[digest]={'initial_B':s['B'],'symmetrizer':s['symmetrizer'],'result':q,'source_sha256':source_sha}
        temporary=file.with_suffix('.tmp');temporary.write_text(json.dumps(cache,indent=2)+'\n',encoding='utf8',newline='\n');temporary.replace(file)
        print(digest[:12],len(s['B']),q['label'],q['status'],flush=True);count+=1
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:
        out.write(json.dumps({'stage':'quiver-witnesses-before-rank-closure','enumeration_complete':False,'new_quivers':count,'wall_seconds':time.monotonic()-start})+'\n')
if __name__=='__main__':main(int(sys.argv[1]))
