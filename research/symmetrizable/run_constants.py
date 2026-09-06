"""Run exact checkpointed upper-triangle tasks, with no time cutoff."""
import argparse,concurrent.futures,hashlib,json,math,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
def run_one(args):
    n,i,table,source_hash,binary_hash,binary=args;directory=HERE/f'rank{n}';logs=directory/'constant_logs';logs.mkdir(exist_ok=True)
    with (logs/f'upper-{i}.txt').open('w') as out:
        result=subprocess.run([str(binary),str(n),str(directory/'constant_tasks'),'0','1',str(i),str(table)],stdout=out,stderr=out)
    if result.returncode:raise RuntimeError(f'Upper {i} exited {result.returncode}')
    path=directory/'constant_tasks'/f'upper-{i}.json';data=json.loads(path.read_text())
    data.update({'source_sha256':source_hash,'binary_sha256':binary_hash,'principal_table_sha256':hashlib.sha256(table.read_bytes()).hexdigest()})
    temporary=path.with_suffix('.json.meta');temporary.write_text(json.dumps(data)+'\n',encoding='utf8',newline='\n');temporary.replace(path)
    return i
def main():
    parser=argparse.ArgumentParser();parser.add_argument('rank',type=int);parser.add_argument('--workers',type=int,default=12)
    parser.add_argument('--binary',default='enumerate_constants_completion');parser.add_argument('--table',default='principal-constants.txt');args=parser.parse_args()
    binary=HERE/'bin'/args.binary
    n=args.rank;directory=HERE/f'rank{n}';directory.mkdir(exist_ok=True);table=HERE/args.table
    todo=[i for i in range(math.factorial(n+1)) if not(directory/'constant_tasks'/f'upper-{i}.json').exists()]
    start=time.monotonic();last=start;done=0
    source_hash=hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest()
    binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    (directory/'constant-running.json').write_text(json.dumps({'started_unix':time.time(),'source_sha256':source_hash,'binary_sha256':binary_hash,'tasks':len(todo)})+'\n',encoding='utf8',newline='\n')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(run_one,(n,i,table,source_hash,binary_hash,binary)) for i in todo]
        for future in concurrent.futures.as_completed(fs):
            future.result();done+=1
            if time.monotonic()-last>=15:
                print(n,'upper tasks',done,'/',len(todo),'seconds',round(time.monotonic()-start,1),flush=True);last=time.monotonic()
    report={'rank':n,'stage':'constant-enumeration','tasks_completed':done,'total_tasks':math.factorial(n+1),'wall_seconds':time.monotonic()-start,
            'source_sha256':source_hash,'binary_sha256':binary_hash,'table_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),'complete':True}
    with (directory/'computation-runs.jsonl').open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps(report)+'\n')
    print(json.dumps(report),flush=True)
    from collect_constants import main as collect
    collect(n)
if __name__=='__main__':main()
