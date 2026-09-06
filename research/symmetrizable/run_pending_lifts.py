"""Compute reusable lifts from finished checkpoints without claiming enumeration complete."""
import argparse,concurrent.futures,hashlib,json,sys,time
from pathlib import Path
from weighted_lifts import task
HERE=Path(__file__).resolve().parent
def main(rank,timeout=2000,retry_unknown=False,encoding='reduced',unknown_only=False,workers=10):
    directory=HERE/f'rank{rank}';path=directory/'pending-lifts.jsonl';done={}
    if path.exists():done={tuple(r['constant_key']):r for r in map(json.loads,path.read_text().splitlines())}
    checkpoints=list((directory/'constant_tasks').glob('upper-*.json'))+list((directory/'constant_parts').glob('upper-*.json'))
    keys=sorted({tuple(k) for p in checkpoints for k in json.loads(p.read_text())['keys']})
    start=time.monotonic();todo=[]
    for number,key in enumerate(keys,1):
        if unknown_only and done.get(key,{}).get('status')!='unknown':continue
        if key in done and (not retry_unknown or done[key]['status']!='unknown'):continue
        d=list(key[:rank]);flat=key[rank:];pair=[[[flat[s*rank*rank+i*rank+j] for j in range(rank)] for i in range(rank)] for s in range(2)]
        c={'id':number,'symmetrizer':d,'N_plus_1':pair[0],'N_minus_1':pair[1]};todo.append((key,c))
    print('Pending enumeration: known constants',len(keys),'lift jobs',len(todo),flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool,path.open('a',encoding='utf8',newline='\n') as out:
        futures={pool.submit(task,'lifts',c,None,None,timeout,'integer',str(directory),encoding):key for key,c in todo}
        last=time.monotonic()
        for future in concurrent.futures.as_completed(futures):
            key=futures[future];result=future.result();result['constant_key']=list(key);done[key]=result
            out.write(json.dumps(result)+'\n');out.flush()
            if time.monotonic()-last>20:
                from collections import Counter
                print(len(done),'/',len(keys),dict(Counter(r['status'] for r in done.values())),round(time.monotonic()-start,1),flush=True);last=time.monotonic()
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as f:f.write(json.dumps({'stage':'lifts-from-incomplete-constant-checkpoints','wall_seconds':time.monotonic()-start,'keys_checked':len(todo),'enumeration_complete':False,'source_sha256':hashlib.sha256((HERE/'weighted_lifts.py').read_bytes()).hexdigest()})+'\n')
def adopt(rank):
    directory=HERE/f'rank{rank}';constants=json.loads((directory/'constant_candidates.json').read_text());assert constants['enumeration_complete']
    old={tuple(r['constant_key']):r for r in map(json.loads,(directory/'pending-lifts.jsonl').read_text().splitlines())}
    path=directory/'lift_feasibility.jsonl';results={r['id']:r for r in map(json.loads,path.read_text().splitlines())} if path.exists() else {}
    for c in constants['candidates']:
        key=tuple(c['symmetrizer'])+tuple(x for a in ('N_plus_1','N_minus_1') for row in c[a] for x in row)
        if key in old and c['id'] not in results:
            r=old[key].copy();del r['constant_key'];r['id']=c['id'];results[c['id']]=r
    path.write_text(''.join(json.dumps(results[k])+'\n' for k in sorted(results)),encoding='utf8',newline='\n');print('Adopted',len(results),'of',constants['count'],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);parser.add_argument('stage',nargs='?',default='lifts',choices=['lifts','adopt'])
    parser.add_argument('--timeout',type=int,default=2000);parser.add_argument('--retry-unknown',action='store_true');parser.add_argument('--unknown-only',action='store_true')
    parser.add_argument('--workers',type=int,default=10)
    parser.add_argument('--encoding',choices=['multiplicity','reduced'],default='reduced');args=parser.parse_args()
    if args.stage=='adopt':adopt(args.rank)
    else:main(args.rank,args.timeout,args.retry_unknown,args.encoding,args.unknown_only,args.workers)
