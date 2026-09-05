"""Exhaust the equal-delay-two rank-four subcase with exact certificates."""
import itertools as it
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'rank3'))
from classify_rank3 import canonical, connected, matmul, certificate, z


def det(a):
    if not a:
        return 1
    return sum((-1)**j * a[0][j] * det([row[:j]+row[j+1:] for row in a[1:]])
               for j in range(len(a)))


def positive(a):
    return all(det([[a[i][j] for j in s] for i in s]) > 0
               for r in range(1, len(a)+1) for s in it.combinations(range(len(a)), r))


def main():
    n = 4
    slots = list(it.combinations_with_replacement(range(n), 2))
    reps = set()
    counts = {'signed_symmetric_matrices': 0, 'connected_positive_commuting': 0}
    for values in it.product((-1, 0, 1), repeat=len(slots)):
        counts['signed_symmetric_matrices'] += 1
        a = [[0]*n for _ in range(n)]
        for (i,j), x in zip(slots, values):
            a[i][j] = a[j][i] = x
        if not connected(a):
            continue
        p = [[max(x,0) for x in row] for row in a]
        m = [[max(-x,0) for x in row] for row in a]
        if matmul(p,m) != matmul(m,p):
            continue
        if not all(positive([[2*(i==j)-a[i][j] for j in range(n)] for i in range(n)])
                   for a in (p,m)):
            continue
        counts['connected_positive_commuting'] += 1
        reps.add(canonical(a))
    output = {'scope': 'All delays equal to 2; D=I; diagonal N0', 'counts': counts, 'classes': []}
    for cid,key in enumerate(sorted(reps),1):
        a = [list(key[n*i:n*i+n]) for i in range(n)]
        p = sp.Matrix([[max(x,0) for x in row] for row in a])
        m = sp.Matrix([[max(-x,0) for x in row] for row in a])
        ap,am = ((1+z*z)*sp.eye(n)-z*s for s in (p,m))
        cert = certificate(ap,am,bound=100)
        assert cert['positive'] and cert['negative'] and cert['labelled_tropical_seed_period']
        output['classes'].append({'id':cid, 'N_plus_1':[[int(x) for x in row] for row in p.tolist()],
                                  'N_minus_1':[[int(x) for x in row] for row in m.tolist()],
                                  'certificate':cert})
        print(cid,cert['positive']['h'],cert['negative']['h'],cert['labelled_tropical_seed_period'],flush=True)
    (HERE/'degree_two.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(counts,'orbits',len(reps),flush=True)


if __name__ == '__main__':
    main()
