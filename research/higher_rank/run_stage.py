"""Charge a subprocess stage to the cumulative one-hour rank budget."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import os
import signal
from budget import remaining_seconds, time_limit

HERE=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('rank',type=int);p.add_argument('stage');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    ledger=HERE/f'rank{a.rank}/computation-runs.jsonl'
    spent=sum(json.loads(x)['wall_seconds'] for x in ledger.read_text().splitlines());remaining=remaining_seconds(ledger.parent,spent)
    if not remaining:print('Rank budget exhausted.');return
    start=time.monotonic();timed_out=False
    child=subprocess.Popen(a.command,start_new_session=(os.name!='nt'))
    try:
        try:code=child.wait(timeout=None if remaining==float('inf') else remaining)
        except subprocess.TimeoutExpired:
            timed_out=True
            if os.name!='nt':os.killpg(child.pid,signal.SIGTERM)
            else:child.terminate()
            code=child.wait()
    finally:
        report={'rank':a.rank,'stage':a.stage,'wall_seconds':time.monotonic()-start,'previous_wall_seconds':spent,
                'budget_seconds':time_limit(ledger.parent),'terminated_at_budget':timed_out,'exit_code':child.returncode,'command':a.command}
        with ledger.open('a',encoding='utf8',newline='\n') as out:out.write(json.dumps(report)+'\n')
        print(json.dumps(report),flush=True)
    if code and not timed_out:raise SystemExit(code)
if __name__=='__main__':main()
