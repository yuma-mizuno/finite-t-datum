"""Exact rank-three experiments in the scope of arXiv:2301.13239.

The degree-two enumeration is exhaustive. The arbitrary-delay classification
is assembled by the other scripts; see CLASSIFICATION.md. Run with SymPy.
"""

from __future__ import annotations

import itertools as it
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from audit_mizuno_finite_type_clocks import DATA, z


def determinant(a):
    if len(a) == 1:
        return a[0][0]
    if len(a) == 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def positive_principal_minors(a):
    return all(determinant([[a[i][j] for j in indices] for i in indices]) > 0
               for k in range(1, len(a) + 1)
               for indices in it.combinations(range(len(a)), k))


def connected(a):
    seen = {0}
    while True:
        more = seen | {j for i in seen for j in range(len(a))
                       if a[i][j] or a[j][i]}
        if more == seen:
            return len(seen) == len(a)
        seen = more


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def canonical(a):
    n = len(a)
    return min(tuple(s * a[p[i]][p[j]] for i in range(n) for j in range(n))
               for p in it.permutations(range(n)) for s in (-1, 1))


def degree_two():
    """All 3^6 possible symmetric signed matrices; no numerical eigensolver."""
    representatives = {}
    counts = {"signed_symmetric_matrices": 0, "connected_positive_commuting": 0}
    slots = list(it.combinations_with_replacement(range(3), 2))
    for values in it.product((-1, 0, 1), repeat=6):
        counts["signed_symmetric_matrices"] += 1
        n = [[0] * 3 for _ in range(3)]
        for (i, j), value in zip(slots, values):
            n[i][j] = n[j][i] = value
        if not connected(n):
            continue
        p = [[max(x, 0) for x in row] for row in n]
        m = [[max(-x, 0) for x in row] for row in n]
        if matmul(p, m) != matmul(m, p):
            continue
        if not all(positive_principal_minors(
                [[2 * (i == j) - b[i][j] for j in range(3)] for i in range(3)])
                for b in (p, m)):
            continue
        counts["connected_positive_commuting"] += 1
        key = canonical(n)
        representatives[key] = [list(key[3*i:3*i+3]) for i in range(3)]
    return counts, representatives


def polynomial_data(ap, am):
    n = ap.rows
    delays = [int(sp.degree(ap[i, i], z)) for i in range(n)]
    n0 = sp.diag(*(1 + z**r for r in delays))
    signs = []
    for a in (ap, am):
        data = {}
        for i in range(n):
            for j in range(n):
                for (p,), value in sp.Poly(sp.expand(n0[i, j] - a[i, j]), z).terms():
                    if value:
                        assert value.is_Integer and value > 0
                        assert 0 < p < delays[i]
                        data[i, j, p] = int(value)
        signs.append(data)
    assert not (signs[0].keys() & signs[1].keys())
    assert all(sp.expand(e) == 0 for e in
               ap * am.subs(z, 1/z).T - am * ap.subs(z, 1/z).T)
    return delays, signs


def exchange_matrix(delays, signs):
    vertices = [(i, p) for i, r in enumerate(delays) for p in range(r)]
    plus, minus = signs
    def coefficient(i, j, p):
        return plus.get((i, j, p), 0) - minus.get((i, j, p), 0)
    b = []
    for i, p in vertices:
        row = []
        for j, q in vertices:
            value = -coefficient(i, j, p-q) + coefficient(j, i, q-p)
            for k in range(len(delays)):
                for v in range(min(p, q) + 1):
                    value += (plus.get((i, k, p-v), 0) * minus.get((j, k, q-v), 0)
                              - minus.get((i, k, p-v), 0) * plus.get((j, k, q-v), 0))
            row.append(value)
        b.append(row)
    assert all(b[i][j] == -b[j][i] for i in range(len(b)) for j in range(len(b)))
    lookup = {v: k for k, v in enumerate(vertices)}
    permutation = [lookup[i, (p-1) % delays[i]] for i, p in vertices]
    mutation = [lookup[i, 0] for i in range(len(delays))]
    assert all(b[i][j] == 0 for i in mutation for j in mutation)
    return vertices, b, permutation, mutation


def mutate(b, c, k):
    n = len(b)
    b1 = [[(-b[i][j] if i == k or j == k else
            b[i][j] + max(-b[i][k], 0) * b[k][j] + b[i][k] * max(b[k][j], 0))
           for j in range(n)] for i in range(n)]
    c1 = [[(-c[k][j] if i == k else
            c[i][j] + max(b[k][i], 0) * c[k][j] - b[k][i] * min(c[k][j], 0))
           for j in range(n)] for i in range(n)]
    return b1, c1


def permute(b, c, permutation):
    inv = [permutation.index(i) for i in range(len(b))]
    return ([[b[inv[i]][inv[j]] for j in range(len(b))] for i in range(len(b))],
            [c[inv[i]][:] for i in range(len(b))])


def step(b, c, permutation, mutation, direction):
    if direction == -1:
        inv = [permutation.index(i) for i in range(len(b))]
        b, c = permute(b, c, inv)
    for k in mutation:
        b, c = mutate(b, c, k)
    if direction == 1:
        b, c = permute(b, c, permutation)
    return b, c


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def negative_permutation(c):
    return (all(x in (-1, 0) for row in c for x in row)
            and all(sum(row) == -1 for row in c)
            and all(sum(col) == -1 for col in zip(*c)))


def certificate(ap, am, bound=300):
    delays, signs = polynomial_data(ap, am)
    vertices, original, permutation, mutation = exchange_matrix(delays, signs)
    n = len(original)
    result = {"delays": delays, "vertices": vertices, "B": original,
              "mutation_vertices": mutation, "relabel_old_to_new": permutation,
              "A_plus": [[str(ap[i, j]) for j in range(ap.cols)] for i in range(ap.rows)],
              "A_minus": [[str(am[i, j]) for j in range(am.cols)] for i in range(am.rows)]}
    b, c = step(original, identity(n), permutation, mutation, 1)
    assert b == original, "Polynomial datum does not give a mutation loop"
    assert step(b, c, permutation, mutation, -1) == (original, identity(n))
    for direction, label in ((1, "positive"), (-1, "negative")):
        b, c = original, identity(n)
        for h in range(1, bound + 1):
            b, c = step(b, c, permutation, mutation, direction)
            assert b == original
            if all(x <= 0 for row in c for x in row):
                assert negative_permutation(c)
                result[label] = {"h": h, "negative_permutation": [row.index(-1) for row in c]}
                break
        else:
            result[label] = None
    b, c = original, identity(n)
    for period in range(1, 2 * bound + 1):
        b, c = step(b, c, permutation, mutation, 1)
        if c == identity(n):
            result["labelled_tropical_seed_period"] = period
            break
    else:
        result["labelled_tropical_seed_period"] = None
    return result


def main():
    counts, reps = degree_two()
    output = {"scope": "D=I, diagonal N0; degree-two enumeration and benchmark certificates; full classification in CLASSIFICATION.md",
              "degree_two_counts": counts, "degree_two": {}, "earlier_examples": {}}
    print(counts, "orbits", len(reps), flush=True)
    for k, n in enumerate(reps.values(), 1):
        p = sp.Matrix([[max(x, 0) for x in row] for row in n])
        m = sp.Matrix([[max(-x, 0) for x in row] for row in n])
        ap, am = (1+z*z)*sp.eye(3)-z*p, (1+z*z)*sp.eye(3)-z*m
        cert = certificate(ap, am)
        assert cert["positive"] and cert["negative"]
        output["degree_two"][f"Q{k}"] = {"signed_N": n, **cert}
        print(f"Q{k}", n, cert["positive"], cert["negative"], cert["labelled_tropical_seed_period"], flush=True)
    for datum in DATA:
        cert = certificate(datum.a_plus, datum.a_minus)
        assert cert["positive"] and cert["negative"], datum.label
        output["earlier_examples"][datum.label] = cert
        print(datum.label, cert["positive"], cert["negative"], cert["labelled_tropical_seed_period"], flush=True)
    Path(__file__).with_name("certificates.json").write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
