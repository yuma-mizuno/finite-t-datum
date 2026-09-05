"""Convert the canonical integer enumeration and check the rank-three control."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
control={tuple(map(int,line.split())) for line in (HERE/'rank3_control.txt').read_text().splitlines()}
old=json.loads((HERE.parent/'rank3'/'constant_candidates.json').read_text())['candidates']
assert control=={tuple(x for field in ('N_plus_1','N_minus_1') for row in c[field] for x in row) for c in old}
keys=[tuple(map(int,line.split())) for line in (HERE/'constant_keys.txt').read_text().splitlines()]
assert keys==sorted(set(keys)) and len(keys)==4865 and all(len(k)==32 for k in keys)
items=[{'id':i,'N_plus_1':[list(k[4*j:4*j+4]) for j in range(4)],
        'N_minus_1':[list(k[16+4*j:20+4*j]) for j in range(4)]} for i,k in enumerate(keys,1)]
(HERE/'constant_candidates.json').write_text(json.dumps({'count':len(items),'candidates':items},indent=2)+'\n',
                                              encoding='utf-8',newline='\n')
print('Rank-three control: exact equality of 180 pairs. Rank-four candidates:',len(items))
