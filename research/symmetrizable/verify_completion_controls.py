"""Compare row-completion pruning with the exhaustive row enumerator."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def keys(p):return {tuple(k) for f in p.glob('upper-*.json') for k in json.loads(f.read_text())['keys']}
def main():
    reference=HERE/'control-rank4-hereditary-reference';completion=HERE/'control-rank4-completion'
    for i in range(120):
        a,b=[json.loads((p/f'upper-{i}.json').read_text()) for p in (reference,completion)]
        assert a['completed'] and b['completed'] and sorted(a['keys'])==sorted(b['keys']),i
    a,b=[json.loads(p.read_text()) for p in (HERE/'rank5/constant_tasks/upper-3.json',HERE/'control-rank5-completion/upper-3.json')]
    assert sorted(a['keys'])==sorted(b['keys']) and len(a['keys'])==480
    result={'rank4_all_120_upper_tasks_equal':True,'rank4_retained_constants':len(keys(reference)),
            'rank5_upper_3_equal':True,'rank5_upper_3_constants':480,
            'source_sha256':hashlib.sha256((HERE/'enumerate_constants.cpp').read_bytes()).hexdigest(),
            'reference_source_sha256':'110f9858c686d66fcfbf36931949dddf0c955a96235c06a3bddb8a0bc0dd0216',
            'table_sha256':hashlib.sha256((HERE/'principal-constants.txt').read_bytes()).hexdigest()}
    (HERE/'completion-controls.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(result)
if __name__=='__main__':main()
