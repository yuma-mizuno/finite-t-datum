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
sys.path.insert(0, str(HERE.parent / "rank3" / "_deps"))
sys.path.append(str(HERE.parent / "rank3"))
import z3
import sympy as sp

from classify_rank3 import certificate, polynomial_data, z


def chain_obstruction(candidate):
    """The four-vertex chain leading-term lemma proved in the manuscript."""
    p,m=candidate['N_plus_1'],candidate['N_minus_1']
    if len(p)!=4: return False
    if any(p[i][j] for i in range(4) for j in range(4) if not (i==3 and j<3)):
        return False
    required={(0,1):1,(0,3):1,(1,0):1,(1,2):1,(2,1):1}
    return all(m[i][j]==required.get((i,j),0) for i in range(4) for j in range(4)
               if (i,j)!=(3,2))


def build_problem(candidate, timeout, encoding="multiplicity", arithmetic="integer"):
    counts = [candidate["N_plus_1"], candidate["N_minus_1"]]
    rank = len(counts[0])
    names = [f"r{i}" for i in range(rank)]
    indices = {}
    for s in range(2):
        for i in range(rank):
            for j in range(rank):
                indices[s, i, j] = []
                for k in range(counts[s][i][j]):
                    indices[s, i, j].append(len(names))
                    names.append(f"p{s}_{i}{j}_{k}")
    variables = [(z3.Int if arithmetic=="integer" else z3.Real)(name) for name in names]
    dim = len(names)
    zero = (0,) * dim
    basis = [tuple(int(i == j) for i in range(dim)) for j in range(dim)]
    a = [[[[(-1, basis[k]) for k in indices[s, i, j]]
           + ([(1, zero), (1, basis[i])] if i == j else [])
           for j in range(rank)] for i in range(rank)] for s in range(2)]
    solver = z3.Solver()
    solver.set(timeout=timeout)
    if chain_obstruction(candidate):
        solver.add(z3.BoolVal(False))
    # Diagonal A+ A-* is reciprocal with constant coefficient 2: its value
    # at 1 is even. Disjoint support makes every N+ N-* constant term zero.
    if any(sum((2*(i==j)-counts[0][i][j])*(2*(i==j)-counts[1][i][j])
               for j in range(rank)) % 2 for i in range(rank)):
        solver.add(z3.BoolVal(False))
    if arithmetic=="integer":
        solver.add(*(v >= 1 for v in variables[:rank]))
    else:
        solver.add(*(v > 0 for v in variables[:rank]))
        solver.add(z3.Sum(variables[:rank])==1)
    for s in range(2):
        for i in range(rank):
            for j in range(rank):
                ks = indices[s, i, j]
                solver.add(*(z3.And(variables[k] > 0, variables[k] < variables[i]) for k in ks))
                solver.add(*(variables[k] <= variables[l] for k, l in zip(ks, ks[1:])))
                if s == 0:
                    solver.add(*(variables[k] != variables[l] for k in ks for l in indices[1, i, j]))
    def expression(v):
        return z3.Sum([a*x for a, x in zip(v, variables) if a]) if any(v) else z3.IntVal(0)
    equations = []
    sort_counter = 0
    def sorted_network(values):
        nonlocal sort_counter
        size=1
        while size<len(values): size*=2
        # Every Laurent exponent is between -sum(r_i) and sum(r_i).
        result=values+[z3.Sum(variables[:rank])+1]*(size-len(values))
        width=2
        while width<=size:
            stride=width//2
            while stride:
                for x in range(size):
                    y=x^stride
                    if y<=x: continue
                    a,b=result[x],result[y]
                    lo,hi=(z3.Ints if arithmetic=="integer" else z3.Reals)(f"sort_{sort_counter}_lo sort_{sort_counter}_hi")
                    sort_counter+=1
                    solver.add(lo==z3.If(a<=b,a,b),hi==z3.If(a<=b,b,a))
                    result[x],result[y]=(lo,hi) if (x&width)==0 else (hi,lo)
                stride//=2
            width*=2
        return result
    for i in range(rank):
        for j in range(i, rank):
            terms = Counter()
            for k in range(rank):
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
            if encoding == "sort":
                solver.add(*(a==b for a,b in zip(sorted_network(left),sorted_network(right))))
            else:
                for e in set(positive):
                    v = expression(e)
                    solver.add(z3.Sum([z3.If(v == f, 1, 0) for f in left])
                               == z3.Sum([z3.If(v == f, 1, 0) for f in right]))
    return solver, variables, indices, equations


def get_matrices(values, indices):
    rank = max(i for s,i,j in indices)+1
    result = []
    for s in range(2):
        a = sp.diag(*(1 + z**r for r in values[:rank]))
        for i in range(rank):
            for j in range(rank):
                a[i, j] -= sum(z**values[k] for k in indices[s, i, j])
        result.append(a)
    return result


def run(candidate, timeout, optimize, dynamics, encoding="multiplicity", arithmetic="integer"):
    rank = len(candidate["N_plus_1"])
    solver, variables, indices, equations = build_problem(candidate, timeout, encoding, arithmetic)
    start = time.monotonic()
    status = solver.check()
    result = {"id": candidate["id"], "status": str(status), "solver_seconds": time.monotonic()-start,
              "delay_bound": None, "encoding": encoding, "arithmetic":arithmetic,
              "variable_names": [str(v) for v in variables]}
    if chain_obstruction(candidate):
        result['analytic_exclusion']='four_vertex_chain_leading_term'
    if status == z3.unknown:
        result["reason"] = solver.reason_unknown()
    if status != z3.sat:
        return result
    model = solver.model()
    rational_values = [sp.Rational(str(model.eval(v))) for v in variables]
    denominator=sp.ilcm(*(v.q for v in rational_values))
    values = [int(v*denominator) for v in rational_values]
    if optimize and arithmetic=="integer":
        opt = z3.Optimize()
        opt.set(timeout=timeout)
        opt.add(*solver.assertions())
        opt.minimize(z3.Sum(variables[:rank]))
        opt.minimize(z3.Sum(variables[rank:]))
        if opt.check() == z3.sat:
            model = opt.model()
            values = [model.eval(v).as_long() for v in variables]
    ap, am = get_matrices(values, indices)
    polynomial_data(ap, am)
    result["values"] = values
    result["A_plus"] = [[str(ap[i,j]) for j in range(rank)] for i in range(rank)]
    result["A_minus"] = [[str(am[i,j]) for j in range(rank)] for i in range(rank)]
    if dynamics and sum(values[:rank]) <= 40:
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
            rank = len(candidate["N_plus_1"])
            print(result["id"], result["status"], round(result["solver_seconds"], 2),
                  result.get("values", [])[:rank],
                  {s: (result.get("certificate", {}).get(s) or {}).get("h") for s in ("positive", "negative")},
                  flush=True)


if __name__ == "__main__":
    main()
