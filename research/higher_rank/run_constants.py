"""Run independent exhaustive upper-triangle tasks within a wall-clock budget.

Use Linux Python (WSL here), after compiling enumerate_constants.cpp.
Completed task files are atomic checkpoints; a timeout never completes a task.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent


def collect(directory,rank):
    tasks=[json.loads(p.read_text()) for p in (directory/'constant_tasks').glob('upper-*.json')]
    total=math.factorial(rank+1)
    indices=[x['upper_index'] for x in tasks]
    assert len(indices)==len(set(indices))
    assert all(x['rank']==rank and x['triangular_count']==total and x['completed'] for x in tasks)
    keys=sorted({tuple(key) for x in tasks for key in x['keys']})
    complete=set(indices)==set(range(total))
    report={'rank':rank,'stage':'constant-enumeration','complete':complete,'completed_tasks':len(tasks),
            'total_tasks':total,'candidate_count':len(keys),'parity_pruning':True,
            'sum_completed_task_seconds':sum(x['seconds'] for x in tasks)}
    (directory/'constant-progress.json').write_text(json.dumps(report,indent=2)+'\n')
    if complete:
        candidates=[{'id':i,'N_plus_1':[list(k[j*rank:(j+1)*rank]) for j in range(rank)],
                     'N_minus_1':[list(k[rank*rank+j*rank:rank*rank+(j+1)*rank]) for j in range(rank)]}
                    for i,k in enumerate(keys,1)]
        (directory/'constant_candidates.json').write_text(json.dumps({'rank':rank,'count':len(candidates),
            'enumeration_complete':True,'parity_pruning':True,'candidates':candidates},indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('rank',type=int)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--binary',default='enumerate_constants_pruned')
    args=parser.parse_args()
    directory=HERE/f'rank{args.rank}';directory.mkdir(exist_ok=True)
    tasks=directory/'constant_tasks';tasks.mkdir(exist_ok=True)
    binary=HERE/'bin'/args.binary;logs=[];children=[]
    ledger=directory/'computation-runs.jsonl'
    previous=[json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
    spent=sum(x['wall_seconds'] for x in previous)
    remaining=max(0,args.seconds-spent)
    source_hash=hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest()
    binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    start=time.monotonic();deadline=start+remaining;terminated=False
    if not remaining:print('Rank computation budget exhausted.',flush=True);return
    try:
        for shard in range(args.workers):
            log=(directory/f'constants-{shard}.log').open('a');logs.append(log)
            children.append(subprocess.Popen([str(binary),str(args.rank),str(tasks),str(shard),str(args.workers)],stdout=log,stderr=subprocess.STDOUT))
        print('Started',[(p.pid,i) for i,p in enumerate(children)],flush=True)
        while any(p.poll() is None for p in children):
            if time.monotonic()>=deadline:
                terminated=True
                for p in children:
                    if p.poll() is None:p.terminate()
                break
            time.sleep(min(15,max(0,deadline-time.monotonic())))
            progress=collect(directory,args.rank)
            print(round(time.monotonic()-start,1),'seconds;',progress['completed_tasks'],'/',progress['total_tasks'],
                  'tasks;',progress['candidate_count'],'distinct candidates; live',sum(p.poll() is None for p in children),flush=True)
        codes=[p.wait() for p in children]
    finally:
        for p in children:
            if p.poll() is None:p.terminate();p.wait()
        for log in logs:log.close()
    result=collect(directory,args.rank)
    result.update({'wall_seconds':time.monotonic()-start,'workers':args.workers,'budget_seconds':args.seconds,
                   'terminated_at_budget':terminated,'exit_codes':codes,
                   'previous_wall_seconds':spent,'source_sha256':source_hash,'binary_sha256':binary_hash})
    with (directory/'computation-runs.jsonl').open('a') as out:out.write(json.dumps(result)+'\n')
    print(json.dumps(result),flush=True)
    assert terminated or all(c==0 for c in codes),codes


if __name__=='__main__':main()
