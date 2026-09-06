"""Resumable higher-rank stages sharing one cumulative wall-clock budget."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import hashlib
import json
from pathlib import Path
import time
from weighted_lifts import task
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'higher_rank'))
from budget import remaining_seconds, time_limit

HERE=Path(__file__).resolve().parent


def read_jsonl(path):
    return {x['id']:x for x in map(json.loads,path.read_text().splitlines())} if path.exists() else {}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('rank',type=int)
    parser.add_argument('stage',choices=['lifts','families','verify'])
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--timeout',type=int,default=2000)
    parser.add_argument('--arithmetic',choices=['integer','rational'],default='integer')
    parser.add_argument('--encoding',choices=['multiplicity','reduced'],default='reduced')
    parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--retry',action='store_true')
    parser.add_argument('--partial',action='store_true',help='Verify resolved candidates while explicitly retaining unresolved cases.')
    args=parser.parse_args();directory=HERE/f'rank{args.rank}'
    ledger=directory/'computation-runs.jsonl'
    previous=[json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
    spent=sum(x['wall_seconds'] for x in previous);remaining=remaining_seconds(directory,spent,args.seconds)
    if not remaining:print('Rank computation budget exhausted.',flush=True);return
    constants=json.loads((directory/'constant_candidates.json').read_text())
    assert constants['enumeration_complete']
    witnesses=read_jsonl(directory/'lift_feasibility.jsonl');families=read_jsonl(directory/'families.jsonl')
    filename={'lifts':'lift_feasibility.jsonl','families':'families.jsonl','verify':'verification.jsonl'}[args.stage]
    target=directory/filename;done=read_jsonl(target)
    status_key={'lifts':'status','families':'coverage_status','verify':'result'}[args.stage]
    todo=[c for c in constants['candidates'] if (c['id'] not in done or (args.retry and done[c['id']].get(status_key) not in ('sat','unsat')))
          and (args.stage!='families' or witnesses.get(c['id'],{}).get('status')=='sat')]
    if args.stage=='verify':
        assert len(witnesses)==constants['count']
        unresolved=[c['id'] for c in constants['candidates'] if witnesses[c['id']]['status'] not in ('sat','unsat')]
        assert args.partial or not unresolved, unresolved
        todo=[c for c in todo if c['id'] not in unresolved]
        assert all(families[cid]['coverage_status']=='unsat' for cid,x in witnesses.items() if x['status']=='sat')
    source_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('weighted_lifts.py','run_lifts.py')}
    source_hashes['rank4_lifts.py']=hashlib.sha256((HERE.parent/'rank4/rank4_lifts.py').read_bytes()).hexdigest()
    start=time.monotonic();deadline=start+remaining
    pool=ProcessPoolExecutor(max_workers=args.workers);pending={};it=iter(todo);completed=0;last=0;terminated=False
    def submit():
        c=next(it,None)
        if c is None:return False
        future=pool.submit(task,args.stage,c,witnesses.get(c['id']),families.get(c['id']),args.timeout,args.arithmetic,str(directory),args.encoding)
        pending[future]=c['id'];return True
    print('Starting',args.stage,'rank',args.rank,'tasks',len(todo),'remaining seconds',round(remaining,1),flush=True)
    try:
        with target.open('a',encoding='utf8',newline='\n') as out:
            for _ in range(2*args.workers):submit()
            while pending:
                if time.monotonic()>=deadline:terminated=True;break
                finished,_=wait(pending,timeout=min(5,max(0,deadline-time.monotonic())),return_when=FIRST_COMPLETED)
                for future in finished:
                    cid=pending.pop(future);x=future.result();assert x['id']==cid
                    out.write(json.dumps(x)+'\n');out.flush();done[cid]=x;completed+=1;submit()
                    if x.get(status_key)=='sat' or x.get('status')=='error':
                        print('Record',cid,x.get(status_key),x.get('values',[])[:args.rank],x.get('error',''),flush=True)
                if time.monotonic()-last>=15:
                    print(completed,'/',len(todo),'finished;',dict(Counter(x.get(status_key,'error') for x in done.values())),
                          'seconds',round(time.monotonic()-start,1),flush=True);last=time.monotonic()
    finally:
        if terminated:
            processes=list(pool._processes.values())
            for process in processes:process.terminate()
            for process in processes:process.join()
            pool.shutdown(wait=False,cancel_futures=True)
        else:pool.shutdown()
    target.write_text(''.join(json.dumps(done[k])+'\n' for k in sorted(done)),encoding='utf8',newline='\n')
    report={'rank':args.rank,'stage':args.stage,'wall_seconds':time.monotonic()-start,'previous_wall_seconds':spent,
            'budget_seconds':time_limit(directory,args.seconds),'terminated_at_budget':terminated,'tasks_scheduled':len(todo),'tasks_completed':completed,
            'status_counts':dict(Counter(x.get(status_key,'error') for x in done.values())),
            'source_sha256':source_hashes['weighted_lifts.py'],'source_hashes':source_hashes,
            'unresolved_lift_ids':[cid for cid,x in witnesses.items() if x['status'] not in ('sat','unsat')] if args.stage=='verify' else []}
    with ledger.open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps(report)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
