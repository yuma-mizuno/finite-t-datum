"""Exhaustive constant-matrix necessary conditions, using the cycle bound 7.

This enumerates coefficient sums, not polynomial lifts or finite T-data.
All tests use exact integer arithmetic (NumPy int64 or Python integers).
"""
from __future__ import annotations

import itertools as it
import json
from math import gcd
from functools import reduce
from pathlib import Path

import numpy as np

from classify_rank3 import canonical, matmul


SLOTS = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]


def strongly_connected(x):
    a, b, c, d, e, f = x
    # In rank three, strong connectivity is equivalent to every vertex
    # having at least one incoming and at least one outgoing edge.
    return all((a or b, c or d, e or f, c or e, a or f, b or d))


def pair_key(p, m):
    return min(tuple(a[q[i]][q[j]] for a in pair for i in range(3) for j in range(3))
               for q in it.permutations(range(3)) for pair in ((p, m), (m, p)))


def simultaneous_witness(p, m):
    normals = [[int(i == j) for j in range(3)] for i in range(3)]
    normals += [[2 * (i == j) - a[i][j] for i in range(3)]
                for a in (p, m) for j in range(3)]
    rays = set()
    for a, b in it.combinations(normals, 2):
        cross = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        for s in (-1, 1):
            ray = tuple(s*x for x in cross)
            if any(ray) and all(sum(x*y for x, y in zip(row, ray)) >= 0 for row in normals):
                g = reduce(gcd, map(abs, ray))
                rays.add(tuple(x//g for x in ray))
    v = tuple(sum(ray[i] for ray in rays) for i in range(3))
    if not all(sum(x*y for x, y in zip(row, v)) > 0 for row in normals):
        return None
    g = reduce(gcd, v)
    return [x//g for x in v]


def main():
    diag = np.array(list(it.product((1, 2), repeat=6)), dtype=np.int64)
    candidates = set()
    maxima_count = split_count = labelled_count = 0
    for x in it.product(range(8), repeat=6):
        a, b, c, d, e, f = x
        if a*c >= 4 or b*e >= 4 or d*f >= 4 or a*d*e >= 8 or b*f*c >= 8:
            continue
        if not strongly_connected(x):
            continue
        maxima_count += 1
        splits = [list({(w, t) for t in range(w+1)} | {(t, w) for t in range(w+1)}) for w in x]
        raw = np.array(list(it.product(*splits)), dtype=np.int64)
        split_count += len(raw)
        p, m = raw[:, :, 0], raw[:, :, 1]
        # Off-diagonal equality of (2I-P)(2I-M)^T with its transpose.
        indices = {(i, j): k for k, (i, j) in enumerate(SLOTS)}
        ok = np.ones((len(raw), len(diag)), dtype=bool)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            k = 3-i-j
            ij, ji = indices[i, j], indices[j, i]
            ik, jk = indices[i, k], indices[j, k]
            value = (-m[:, ji, None]*diag[None, :, i]
                     +m[:, ij, None]*diag[None, :, j]
                     -p[:, ij, None]*diag[None, :, 3+j]
                     +p[:, ji, None]*diag[None, :, 3+i]
                     +(p[:, ik]*m[:, jk]-p[:, jk]*m[:, ik])[:, None])
            ok &= value == 0
        for s, offset in ((p, 0), (m, 3)):
            d0, d1, d2 = (diag[None, :, offset+i] for i in range(3))
            ac, be, df = (s[:, i, None]*s[:, j, None] for i, j in ((0, 2), (1, 4), (3, 5)))
            ok &= (d0*d1 > ac) & (d0*d2 > be) & (d1*d2 > df)
            det = d0*d1*d2-d0*df-d1*be-d2*ac-(s[:, 0]*s[:, 3]*s[:, 4]+s[:, 1]*s[:, 2]*s[:, 5])[:, None]
            ok &= det > 0
        rows, ds = np.where(ok)
        labelled_count += len(rows)
        for row, di in zip(rows, ds):
            pp = [[0]*3 for _ in range(3)]
            mm = [[0]*3 for _ in range(3)]
            for i in range(3):
                pp[i][i], mm[i][i] = int(2-diag[di, i]), int(2-diag[di, 3+i])
            for k, (i, j) in enumerate(SLOTS):
                pp[i][j], mm[i][j] = int(p[row, k]), int(m[row, k])
            candidates.add(pair_key(pp, mm))
    result = []
    for key in sorted(candidates):
        p = [list(key[3*i:3*i+3]) for i in range(3)]
        m = [list(key[9+3*i:9+3*i+3]) for i in range(3)]
        v = simultaneous_witness(p, m)
        if v:
            result.append({"id": len(result)+1, "N_plus_1": p, "N_minus_1": m, "positive_left_vector": v})
    counts = {"maxima_patterns": maxima_count, "off_diagonal_splits": split_count,
              "labelled_symplectic_M_matrix_pairs": labelled_count,
              "orbits_before_simultaneous_positivity": len(candidates),
              "orbits_after_simultaneous_positivity": len(result)}
    print(counts, flush=True)
    Path(__file__).with_name("constant_candidates.json").write_text(json.dumps({"counts": counts, "candidates": result}, indent=2)+"\n")


if __name__ == "__main__":
    main()
