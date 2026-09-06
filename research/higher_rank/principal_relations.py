"""Necessary linear lift relations inherited from classified principal data."""
import itertools as it
import json
from pathlib import Path
import re
from functools import lru_cache


@lru_cache(maxsize=1)
def labelled_records():
    records=json.loads((Path(__file__).resolve().parents[2]/'docs/catalogue/catalogue.json').read_text(encoding='utf8'))['records']
    lookup={n:{} for n in range(3,6)}
    for record in records:
        if record['scope']['symmetrizer']!='identity':continue
        n=record['rank']
        if n not in lookup:continue
        for exchange in range(2):
            for perm in it.permutations(range(n)):
                # perm maps each lower-record species to a position in S.
                inverse=[perm.index(i) for i in range(n)]
                constants=tuple(2*(i==j)-record['datum'][a][inverse[i]][inverse[j]]
                                for a in (('A_plus_1','A_minus_1') if not exchange else ('A_minus_1','A_plus_1'))
                                for i in range(n) for j in range(n))
                lookup[n].setdefault(constants,(record,exchange,perm))
    return lookup


def add_principal_relations(solver,variables,indices,candidate,z3):
    lookup=labelled_records()
    p,m=candidate['N_plus_1'],candidate['N_minus_1'];n=len(p);certificates=[]
    for size in range(min(n-1,5),2,-1):
        for S in it.combinations(range(n),size):
            external=[k for k in range(n) if k not in S]
            if any(any(p[i][k] for i in S) and any(m[i][k] for i in S) for k in external):continue
            constants=tuple(b[i][j] for b in (p,m) for i in S for j in S)
            found=lookup[size].get(constants)
            if not found:continue
            record,exchange,positions=found;perm=tuple(S[i] for i in positions);mapping=[]
            for name in record['family']['variable_names']:
                if name.startswith('r'):mapping.append(perm[int(name[1:])]);continue
                match=re.fullmatch(r'p([01])_(\d)(\d)_(\d+)',name);assert match,name
                s,i,j,k=map(int,match.groups());mapping.append(indices[s^exchange,perm[i],perm[j]][k])
            rows=record['family']['rref']
            solver.add(*(z3.Sum([z3.RealVal(c)*variables[k] for c,k in zip(row,mapping) if c!='0'])==0 for row in rows))
            certificates.append({'subset':list(S),'catalogue_id':record['id'],'sign_exchange':bool(exchange),
                                 'lower_to_original_species':list(perm),'lower_to_original_variables':mapping})
    return certificates
