"""Controls of pulse exclusions and normalized delays against known results."""
import hashlib
import json
from pathlib import Path
import time
from fractions import Fraction
from obstructions import attachments
from finite_extension import core,patterns,exclusion

HERE=Path(__file__).resolve().parent


def main():
    start=time.monotonic();report=[]
    for rank in (5,6):
        directory=HERE/f'rank{rank}'
        cs={c['id']:c for c in json.loads((directory/'constant_candidates.json').read_text())['candidates']}
        witnesses={x['id']:x for x in map(json.loads,(directory/'lift_feasibility.jsonl').read_text().splitlines())}
        exclusions=0;positive=[]
        for cid,w in witnesses.items():
            if w['status']=='sat':
                assert exclusion(cs[cid]) is None,cid
                for s,ell,a,top in attachments(cs[cid]):
                    found=core(cs[cid],s,ell,a,top)
                    if found is None:continue
                    record,flip,local_a,target,diagonal,permutation=found
                    result=patterns(record['id'],flip,local_a,target,diagonal)
                    scale=Fraction(w['values'][permutation[0]],record['datum']['delays'][0]);R=Fraction(w['values'][ell])/scale
                    exponent=Fraction(w['values'][w['variable_names'].index(f'p{s}_{ell}{ell}_0')])/scale if diagonal else None
                    assert any(x['delay']==R and x['diagonal_exponent']==exponent for x in result['patterns']),(cid,result,R,exponent)
                    positive.append({'constant_id':cid,'core':record['id'],'normalized_delay':str(R),'matching_pulse_pattern':True})
            elif rank==5 and exclusion(cs[cid]):
                assert w['status']=='unsat',cid;exclusions+=1
        report.append({'rank':rank,'previously_verified_exclusions_reproduced':exclusions,'known_positive_extension_checks':positive})
    result={'source_sha256':hashlib.sha256((HERE/'finite_extension.py').read_bytes()).hexdigest(),'controls':report}
    (HERE/'finite-extension-controls.json').write_text(json.dumps(result,indent=2)+'\n')
    with (HERE/'rank6/computation-runs.jsonl').open('a') as out:out.write(json.dumps({'rank':6,'stage':'finite-pulse-extension-controls','wall_seconds':time.monotonic()-start,'result':'All controls passed'})+'\n')
    print([(x['rank'],x['previously_verified_exclusions_reproduced'],len(x['known_positive_extension_checks'])) for x in report],flush=True)


if __name__=='__main__':main()
