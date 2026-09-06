"""Compare the difficult full rank-six upper task with independently filtered reference keys."""
import hashlib,json
from pathlib import Path
from verify_partition_controls import allowed
HERE=Path(__file__).resolve().parent
def main():
    reference=json.loads((HERE/'control-rank6-mmatrix/upper-3.json').read_text())
    actual=json.loads((HERE/'control-rank6-cost/upper-3.json').read_text())
    expected={tuple(k) for k in reference['keys'] if allowed(tuple(k),6)}
    got={tuple(k) for k in actual['keys']}
    assert expected==got,(len(expected),len(got),len(expected-got),len(got-expected))
    p=HERE/'branch-cost-controls.json';result=json.loads(p.read_text())
    result['rank6_full_upper_control']={'upper_task':3,'reference_constants':len(reference['keys']),'filtered_constants':len(got),'independently_filtered_reference_equal':True,'output_sha256':hashlib.sha256((HERE/'control-rank6-cost/upper-3.json').read_bytes()).hexdigest()}
    p.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(result)
if __name__=='__main__':main()
