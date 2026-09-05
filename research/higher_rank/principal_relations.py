"""Necessary linear lift relations inherited from classified principal data."""
import itertools as it
import json
from pathlib import Path
import re


def add_principal_relations(solver,variables,indices,candidate,z3):
    records=json.loads((Path(__file__).resolve().parents[2]/'docs/catalogue/catalogue.json').read_text(encoding='utf8'))['records']
    p,m=candidate['N_plus_1'],candidate['N_minus_1'];n=len(p);certificates=[]
    for size in range(min(n-1,4),2,-1):
        for S in it.combinations(range(n),size):
            external=[k for k in range(n) if k not in S]
            if any(any(p[i][k] for i in S) and any(m[i][k] for i in S) for k in external):continue
            found=None
            for record in records:
                if record['rank']!=size:continue
                for exchange in range(2):
                    counts=[p,m] if not exchange else [m,p]
                    for perm in it.permutations(S):
                        if all(counts[s][perm[i]][perm[j]]==2*(i==j)-record['datum'][a][i][j]
                               for s,a in enumerate(('A_plus_1','A_minus_1')) for i in range(size) for j in range(size)):
                            found=(record,exchange,perm);break
                    if found:break
                if found:break
            if not found:continue
            record,exchange,perm=found;mapping=[]
            for name in record['family']['variable_names']:
                if name.startswith('r'):mapping.append(perm[int(name[1:])]);continue
                match=re.fullmatch(r'p([01])_(\d)(\d)_(\d+)',name);assert match,name
                s,i,j,k=map(int,match.groups());mapping.append(indices[s^exchange,perm[i],perm[j]][k])
            rows=record['family']['rref']
            solver.add(*(z3.Sum([z3.RealVal(c)*variables[k] for c,k in zip(row,mapping) if c!='0'])==0 for row in rows))
            certificates.append({'subset':list(S),'catalogue_id':record['id'],'sign_exchange':bool(exchange),
                                 'lower_to_original_species':list(perm),'lower_to_original_variables':mapping})
    return certificates
