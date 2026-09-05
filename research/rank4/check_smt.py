"""Replay one saved SMT formula in a fresh process and report exact status."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'rank3'/'_deps'))
import z3

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('path',type=Path)
    parser.add_argument('--timeout',type=int,default=180000)
    parser.add_argument('--seed',type=int,default=0)
    args=parser.parse_args()
    query=args.path.read_bytes()
    solver=z3.Solver();solver.from_string(query.decode('utf-8'))
    solver.set(timeout=args.timeout,random_seed=args.seed)
    start=time.monotonic();status=solver.check()
    print(json.dumps({'result':str(status),'sha256':hashlib.sha256(query).hexdigest(),
                      'timeout_ms':args.timeout,'random_seed':args.seed,'seconds':time.monotonic()-start}))
