"""Unbounded integer feasibility for polynomial lifts of constant candidates.

SAT gives one exact witness. UNSAT excludes that constant pair at all delays.
UNKNOWN is retained. A SAT witness does not classify all lifts of its pair.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "_deps"))
import z3
import sympy as sp

from classify_rank3 import certificate, polynomial_data, z


def build_problem(candidate, timeout):
    counts = [candidate["N_plus_1"], candidate["N_minus_1"]]
    names = ["r0", "r1", "r2"]
    indices = {}
    for s in range(2):
        for i in range(3):
            for j in range(3):
                indices[s, i, j] = []
                for k in range(counts[s][i][j]):
                    indices[s, i, j].append(len(names))
                    names.append(f"p{s}_{i}{j}_{k}")
    variables = [z3.Int(name) for name in names]
    dim = len(names)
    zero = (0,) * dim
    basis = [tuple(int(i == j) for i in range(dim)) for j in range(dim)]
    a = [[[[(-1, basis[k]) for k in indices[s, i, j]]
           + ([(1, zero), (1, basis[i])] if i == j else [])
           for j in range(3)] for i in range(3)] for s in range(2)]
    solver = z3.Solver()
    solver.set(timeout=timeout)
    solver.add(*(v >= 1 for v in variables[:3]))
    for s in range(2):
        for i in range(3):
            for j in range(3):
                ks = indices[s, i, j]
                solver.add(*(z3.And(variables[k] > 0, variables[k] < variables[i]) for k in ks))
                solver.add(*(variables[k] <= variables[l] for k, l in zip(ks, ks[1:])))
                if s == 0:
                    solver.add(*(variables[k] != variables[l] for k in ks for l in indices[1, i, j]))
    def expression(v):
        return z3.Sum([a*x for a, x in zip(v, variables) if a]) if any(v) else z3.IntVal(0)
    equations = []
    for i in range(3):
        for j in range(i, 3):
            terms = Counter()
            for k in range(3):
                for s in range(2):
                    for c, e in a[s][i][k]:
                        for d, f in a[1-s][j][k]:
                            terms[tuple(x-y for x, y in zip(e, f))] += (1-2*s)*c*d
            positive = [e for e, c in sorted(terms.items()) for _ in range(max(c, 0))]
            negative = [e for e, c in sorted(terms.items()) for _ in range(max(-c, 0))]
            assert len(positive) == len(negative)
            if not positive:
                continue
            equations.append((positive, negative))
            left, right = list(map(expression, positive)), list(map(expression, negative))
            solver.add(z3.Sum(left) == z3.Sum(right))
            # Equal multiplicities at every positive-side exponent, together
            # with equal total sizes, is exactly equality of the multisets.
            for e in set(positive):
                v = expression(e)
                solver.add(z3.Sum([z3.If(v == f, 1, 0) for f in left])
                           == z3.Sum([z3.If(v == f, 1, 0) for f in right]))
    return solver, variables, indices, equations


def get_matrices(values, indices):
    result = []
    for s in range(2):
        a = sp.diag(*(1 + z**r for r in values[:3]))
        for i in range(3):
            for j in range(3):
                a[i, j] -= sum(z**values[k] for k in indices[s, i, j])
        result.append(a)
    return result


def run(candidate, timeout, optimize, dynamics):
    solver, variables, indices, equations = build_problem(candidate, timeout)
    start = time.monotonic()
    status = solver.check()
    result = {"id": candidate["id"], "status": str(status), "solver_seconds": time.monotonic()-start,
              "delay_bound": None, "variable_names": [str(v) for v in variables]}
    if status == z3.unknown:
        result["reason"] = solver.reason_unknown()
    if status != z3.sat:
        return result
    model = solver.model()
    values = [model.eval(v).as_long() for v in variables]
    if optimize:
        opt = z3.Optimize()
        opt.set(timeout=timeout)
        opt.add(*solver.assertions())
        opt.minimize(z3.Sum(variables[:3]))
        opt.minimize(z3.Sum(variables[3:]))
        if opt.check() == z3.sat:
            model = opt.model()
            values = [model.eval(v).as_long() for v in variables]
    ap, am = get_matrices(values, indices)
    polynomial_data(ap, am)
    result["values"] = values
    result["A_plus"] = [[str(ap[i,j]) for j in range(3)] for i in range(3)]
    result["A_minus"] = [[str(am[i,j]) for j in range(3)] for i in range(3)]
    if dynamics and sum(values[:3]) <= 40:
        cert = certificate(ap, am, bound=100)
        result["certificate"] = cert
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=3000)
    parser.add_argument("--ids", type=int, nargs="*")
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--dynamics", action="store_true")
    parser.add_argument("--output", default="lift_feasibility.jsonl")
    args = parser.parse_args()
    candidates = json.loads((HERE / "constant_candidates.json").read_text())["candidates"]
    with (HERE / args.output).open("w", encoding="utf-8") as out:
        for candidate in candidates:
            if args.ids and candidate["id"] not in args.ids:
                continue
            result = run(candidate, args.timeout, args.optimize, args.dynamics)
            out.write(json.dumps(result) + "\n")
            out.flush()
            print(result["id"], result["status"], round(result["solver_seconds"], 2),
                  result.get("values", [])[:3],
                  {s: (result.get("certificate", {}).get(s) or {}).get("h") for s in ("positive", "negative")},
                  flush=True)


if __name__ == "__main__":
    main()
