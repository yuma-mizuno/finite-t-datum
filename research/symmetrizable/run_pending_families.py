"""Build reusable full lift spaces while the independent constant search continues."""
import concurrent.futures,hashlib,json,sys,time
from pathlib import Path
from weighted_lifts import task
HERE=Path(__file__).resolve().parent

def read(path):
    if not path.exists():return {}
    text=path.read_text()
    if text and not text.endswith('\n'):text=text.rsplit('\n',1)[0] if '\n' in text else ''
    return {tuple(r['constant_key']):r for r in map(json.loads,text.splitlines())}

def main(rank,adopt=False):
    directory=HERE/f'rank{rank}';output=directory/'pending-families.jsonl';done=read(output)
    if adopt:
        constants=json.loads((directory/'constant_candidates.json').read_text());assert constants['enumeration_complete']
        path=directory/'families.jsonl';results={r['id']:r for r in map(json.loads,path.read_text().splitlines())} if path.exists() else {}
        for c in constants['candidates']:
            key=tuple(c['symmetrizer'])+tuple(x for s in ('N_plus_1','N_minus_1') for row in c[s] for x in row)
            if key in done and c['id'] not in results:
                r=done[key].copy();del r['constant_key'];r['id']=c['id'];results[c['id']]=r
        path.write_text(''.join(json.dumps(results[k])+'\n' for k in sorted(results)),encoding='utf8',newline='\n');print('Adopted',len(results),'family computations',flush=True);return
    witnesses=read(directory/'pending-lifts.jsonl');jobs=[]
    for key,w in witnesses.items():
        if w['status']!='sat' or done.get(key,{}).get('coverage_status')=='unsat':continue
        flat=key[rank:];pair=[[[flat[s*rank*rank+i*rank+j] for j in range(rank)] for i in range(rank)] for s in range(2)]
        c={'id':w['id'],'symmetrizer':list(key[:rank]),'N_plus_1':pair[0],'N_minus_1':pair[1]};jobs.append((key,c,w))
    start=time.monotonic();print('Pending enumeration:',len(jobs),'complete-family jobs',flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool,output.open('a',encoding='utf8',newline='\n') as out:
        futures={pool.submit(task,'families',c,w,None,10000,'integer',str(directory),'reduced'):key for key,c,w in jobs}
        for future in concurrent.futures.as_completed(futures):
            r=future.result();r['constant_key']=list(futures[future]);out.write(json.dumps(r)+'\n');out.flush()
            print(r['id'],r.get('coverage_status',r.get('status')),[(s['dimension'],s['certificate']['labelled_tropical_seed_period']) for s in r.get('spaces',[])],flush=True)
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:
        out.write(json.dumps({'stage':'families-from-incomplete-constant-checkpoints','enumeration_complete':False,'keys_checked':len(jobs),'wall_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256((HERE/'weighted_lifts.py').read_bytes()).hexdigest()})+'\n')
if __name__=='__main__':main(int(sys.argv[1]),len(sys.argv)>2 and sys.argv[2]=='adopt')
