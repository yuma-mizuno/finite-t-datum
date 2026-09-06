"""Checkpoint short partner blocks; preserve and combine every exact upper task."""
import argparse,concurrent.futures,hashlib,json,math,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
def run_one(args):
    n,a,lo,hi,directory,table,binary,meta=args
    name=f'upper-{a}-part-{lo}-{hi}';path=directory/'constant_parts'/(name+'.json');log=directory/'constant_logs'/(name+'.txt')
    with log.open('w') as out:
        p=subprocess.run([str(binary),str(n),str(path.parent),'0','1',str(a),str(table),str(lo),str(hi)],stdout=out,stderr=out)
    if p.returncode:raise RuntimeError(f'{name} exited {p.returncode}; see {log}')
    data=json.loads(path.read_text());assert data['completed'];data.update(meta);data['partner_interval']=[lo,hi]
    temporary=path.with_suffix('.meta');temporary.write_text(json.dumps(data)+'\n',encoding='utf8',newline='\n');temporary.replace(path)
    return a
def main():
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);parser.add_argument('--workers',type=int,default=12)
    parser.add_argument('--block',type=int,default=256);parser.add_argument('--output');parser.add_argument('--binary',default='enumerate_constants_cost');args=parser.parse_args()
    n=args.rank;count=math.factorial(n+1);directory=HERE/(args.output or f'rank{n}');directory.mkdir(exist_ok=True)
    for name in ('constant_parts','constant_logs','constant_tasks'):(directory/name).mkdir(exist_ok=True)
    table=HERE/('principal-constants-rank6.txt' if n==6 else 'principal-constants.txt');binary=HERE/'bin'/args.binary
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    meta={'source_sha256':sha(HERE/'enumerate_constants.cpp'),'binary_sha256':sha(binary),'binary_name':binary.name,'principal_table_sha256':sha(table)}
    if Path(str(table)+'.weights').exists():meta['weighted_principal_table_sha256']=sha(Path(str(table)+'.weights'))
    completed={json.loads(p.read_text())['upper_index'] for p in (directory/'constant_tasks').glob('upper-*.json')}
    jobs=[];remaining={};intervals={}
    for a in range(count):
        if a in completed:continue
        intervals[a]=[(lo,min(lo+args.block,count)) for lo in range(a,count,args.block)]
        missing=[]
        for lo,hi in intervals[a]:
            path=directory/'constant_parts'/f'upper-{a}-part-{lo}-{hi}.json'
            if not path.exists():missing.append((n,a,lo,hi,directory,table,binary,meta))
        jobs.extend(missing);remaining[a]=len(missing)
    # Distribute short blocks across upper tasks rather than leave a whole
    # expensive upper triangle to one straggling worker.
    jobs.sort(key=lambda task:(task[2]-task[1],task[1]))
    def aggregate(a):
        parts=[json.loads((directory/'constant_parts'/f'upper-{a}-part-{lo}-{hi}.json').read_text()) for lo,hi in intervals[a]]
        assert all(p['completed'] and p['upper_index']==a for p in parts)
        keys=sorted({tuple(k) for p in parts for k in p['keys']})
        sources=sorted({p['source_sha256'] for p in parts})
        result={'rank':n,'upper_index':a,'triangular_count':count,'completed':True,'keys':keys,**meta,
                'seconds':sum(p['seconds'] for p in parts),'partner_blocks':intervals[a],'full_subset_pruning':True,
                'parts':len(parts),'parts_archive':'constant_parts.zip'}
        result['source_sha256']=sources;result['binary_sha256']=sorted({p['binary_sha256'] for p in parts})
        (directory/'constant_tasks'/f'upper-{a}.json').write_text(json.dumps(result)+'\n',encoding='utf8',newline='\n');completed.add(a)
    for a,number in remaining.items():
        if not number:aggregate(a)
    start=time.monotonic();last=start;done=0
    (directory/'constant-split-running.json').write_text(json.dumps({'started_unix':time.time(),'tasks':len(jobs),'block_size':args.block,**meta})+'\n',encoding='utf8',newline='\n')
    print('Rank',n,'partner blocks',len(jobs),'preserved upper tasks',len(completed),flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_one,job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            a=future.result();remaining[a]-=1;done+=1
            if not remaining[a]:aggregate(a)
            if time.monotonic()-last>20:
                print('Rank',n,'blocks',done,'/',len(jobs),'upper tasks',len(completed),'/',count,'seconds',round(time.monotonic()-start,1),flush=True);last=time.monotonic()
    assert len(completed)==count
    report={'rank':n,'stage':'constant-enumeration-partner-blocks','wall_seconds':time.monotonic()-start,'complete':True,'completed_upper_tasks':count,'partner_blocks':len(jobs),**meta}
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps(report)+'\n')
    print(report,flush=True)
    if not args.output:
        from collect_constants import main as collect
        collect(n)
if __name__=='__main__':main()
