"""Replay resolved checkpoints now; adopt by full constant key after enumeration closes."""
import argparse,concurrent.futures,hashlib,json,sys,time
from pathlib import Path
from weighted_lifts import task
HERE=Path(__file__).resolve().parent

def keyed(path):
    if not path.exists():return {}
    text=path.read_text()
    if text and not text.endswith('\n'):text=text.rsplit('\n',1)[0] if '\n' in text else ''
    return {tuple(r['constant_key']):r for r in map(json.loads,text.splitlines())}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);parser.add_argument('--workers',type=int,default=8);parser.add_argument('--timeout',type=int,default=5000);parser.add_argument('--adopt',action='store_true');args=parser.parse_args()
    rank=args.rank;directory=HERE/f'rank{rank}';target=directory/'pending-verification.jsonl';done=keyed(target)
    work=directory/'pending-verification';work.mkdir(exist_ok=True)
    if args.adopt:
        constants=json.loads((directory/'constant_candidates.json').read_text());assert constants['enumeration_complete']
        outpath=directory/'verification.jsonl';existing={r['id']:r for r in map(json.loads,outpath.read_text().splitlines())} if outpath.exists() else {}
        queries=directory/'smt_queries';queries.mkdir(exist_ok=True)
        for c in constants['candidates']:
            key=tuple(c['symmetrizer'])+tuple(x for s in ('N_plus_1','N_minus_1') for row in c[s] for x in row)
            if key not in done or c['id'] in existing:continue
            r=done[key].copy();del r['constant_key'];r['id']=c['id']
            if 'query_file' in r:
                payload=(work/'smt_queries'/r['query_file']).read_bytes();name=f'{c["id"]:06d}_'+r['query_file'].split('_',1)[1]
                destination=queries/name
                if destination.exists():assert destination.read_bytes()==payload
                else:destination.write_bytes(payload)
                r['query_file']=name
            existing[c['id']]=r
        outpath.write_text(''.join(json.dumps(existing[k])+'\n' for k in sorted(existing)),encoding='utf8',newline='\n');print('Adopted',len(existing),'exact replays',flush=True);return
    witnesses=keyed(directory/'pending-lifts.jsonl');families=keyed(directory/'pending-families.jsonl')
    ids_path=work/'stable-ids.jsonl';ids={key:r['id'] for key,r in keyed(ids_path).items()};number=max(ids.values(),default=0);jobs=[]
    with ids_path.open('a',encoding='utf8',newline='\n') as out:
        for key,w in witnesses.items():
            if w['status'] not in ('sat','unsat') or done.get(key,{}).get('result')=='unsat':continue
            f=families.get(key)
            if w['status']=='sat' and (not f or f.get('coverage_status')!='unsat'):continue
            if key not in ids:
                number+=1;ids[key]=number;out.write(json.dumps({'id':number,'constant_key':list(key)})+'\n')
            flat=key[rank:];pair=[[[flat[s*rank*rank+i*rank+j] for j in range(rank)] for i in range(rank)] for s in range(2)]
            c={'id':ids[key],'symmetrizer':list(key[:rank]),'N_plus_1':pair[0],'N_minus_1':pair[1]};jobs.append((key,c,w,f))
    source_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ('weighted_lifts.py','run_pending_verification.py')}
    start=time.monotonic();last=start;print('Incomplete enumeration:',len(jobs),'independent replay jobs',flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool,target.open('a',encoding='utf8',newline='\n') as out:
        pending={};iterator=iter(jobs);completed=0
        def submit():
            item=next(iterator,None)
            if item is None:return
            key,c,w,f=item;pending[pool.submit(task,'verify',c,w,f,args.timeout,'integer',str(work),'multiplicity')]=key
        for _ in range(2*args.workers):submit()
        while pending:
            finished,_=concurrent.futures.wait(pending,timeout=5,return_when=concurrent.futures.FIRST_COMPLETED)
            for future in finished:
                r=future.result();r['constant_key']=list(pending.pop(future));out.write(json.dumps(r)+'\n');out.flush();completed+=1;submit()
                if r.get('result')!='unsat':print('Unresolved replay',r['id'],r.get('result',r.get('status')),r.get('error',''),flush=True)
            if time.monotonic()-last>20:
                print(completed,'/',len(jobs),'replayed;',round(time.monotonic()-start,1),'seconds',flush=True);last=time.monotonic()
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:
        out.write(json.dumps({'stage':'replays-from-incomplete-constant-checkpoints','enumeration_complete':False,'keys_checked':len(jobs),'wall_seconds':time.monotonic()-start,
                              'source_hashes':source_hashes})+'\n')
if __name__=='__main__':main()
