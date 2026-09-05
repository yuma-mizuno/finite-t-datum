"""Resumable process-parallel unbounded lift search; UNKNOWN is preserved."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time

from rank4_lifts import run

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--timeout',type=int,default=2000)
    parser.add_argument('--retry-unknown',action='store_true')
    parser.add_argument('--encoding',choices=('multiplicity','sort'),default='multiplicity')
    parser.add_argument('--arithmetic',choices=('integer','rational'),default='integer')
    args=parser.parse_args()
    candidates=json.loads((HERE/'constant_candidates.json').read_text())['candidates']
    target=HERE/'lift_feasibility.jsonl'
    previous={}
    if target.exists():
        previous={x['id']:x for x in map(json.loads,target.read_text().splitlines())}
    todo=[c for c in candidates if c['id'] not in previous or
          (args.retry_unknown and previous[c['id']]['status']=='unknown')]
    start=time.monotonic()
    with target.open('a',encoding='utf-8') as out, ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run,c,args.timeout,True,False,args.encoding,args.arithmetic):c['id'] for c in todo}
        for i,future in enumerate(as_completed(futures),1):
            result=future.result()
            assert result['variable_names'][:4] == ['r0','r1','r2','r3']
            if result['status']=='sat':
                assert len(result['A_plus'])==len(result['A_minus'])==4
            out.write(json.dumps(result)+'\n');out.flush()
            previous[result['id']]=result
            if i%50==0 or result['status']=='sat':
                print(i,'/',len(todo),dict(Counter(x['status'] for x in previous.values())),
                      'latest',result['id'],result['status'],result.get('values',[])[:4],
                      'seconds',round(time.monotonic()-start,1),flush=True)
    ordered=[previous[k] for k in sorted(previous)]
    target.write_text(''.join(json.dumps(x)+'\n' for x in ordered),encoding='utf-8',newline='\n')
    print('DONE',dict(Counter(x['status'] for x in ordered)),flush=True)


if __name__=='__main__':
    main()
