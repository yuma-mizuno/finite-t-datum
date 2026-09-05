"""Cover every unbounded polynomial lift by exact rational linear spaces.

For a satisfying delay tuple, pair equal exponents on the two sides of each
Laurent identity. The resulting homogeneous linear equations imply the whole
identity. Block this entire space and repeat; final UNSAT proves coverage.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

import sympy as sp

from solve_lifts import HERE, build_problem, z3


def matching_space(values, equations):
    rows = []
    for positive, negative in equations:
        left, right = defaultdict(list), defaultdict(list)
        for e in positive:
            left[sum(x*y for x, y in zip(e, values))].append(e)
        for e in negative:
            right[sum(x*y for x, y in zip(e, values))].append(e)
        assert left.keys() == right.keys()
        for value in sorted(left):
            assert len(left[value]) == len(right[value])
            for a, b in zip(sorted(left[value]), sorted(right[value])):
                rows.append([x-y for x, y in zip(a, b)])
    rref, pivots = sp.Matrix(rows).rref()
    return rref[:len(pivots), :]


def main():
    candidates = json.loads((HERE / "constant_candidates.json").read_text())["candidates"]
    witnesses = [json.loads(line) for line in (HERE / "lift_feasibility.jsonl").read_text().splitlines()]
    witness_map = {c["id"]: c for c in witnesses if c["status"] == "sat"}
    results = []
    for candidate in candidates:
        cid = candidate["id"]
        if cid not in witness_map:
            continue
        solver, variables, indices, equations = build_problem(candidate, 10000)
        spaces = []
        while True:
            status = solver.check()
            if status != z3.sat:
                break
            model = solver.model()
            values = [model.eval(v).as_long() for v in variables]
            space = matching_space(values, equations)
            spaces.append(space)
            solver.add(z3.Not(z3.And(*[z3.Sum([z3.RealVal(str(c))*v for c, v in zip(row, variables)]) == 0
                                      for row in space.tolist()])))
            if len(spaces) > 100:
                raise RuntimeError(f"Unexpectedly many spaces: {cid}")
        maximal = []
        for a in spaces:
            if any(a.rows > b.rows and b.col_join(a).rank() == a.rows for b in spaces):
                continue
            if a not in maximal:
                maximal.append(a)
        item = {"id": cid, "coverage_status": str(status), "raw_space_count": len(spaces), "spaces": []}
        for a in maximal:
            dim = len(variables)-a.rows
            gauge = []
            for species in (0, 1):
                vector = [0]*len(variables)
                for (s, i, j), ks in indices.items():
                    for k in ks:
                        vector[k] = int(i == species)-int(j == species)
                gauge.append(vector)
            scale = witness_map[cid]["values"]
            gauge_scale = sp.Matrix.hstack(*(sp.Matrix(v) for v in [scale]+gauge))
            exact_gauge_scale = (dim == 3 and gauge_scale.rank() == 3 and a*gauge_scale == sp.zeros(a.rows, 3))
            item["spaces"].append({"dimension": dim, "rref": [[str(c) for c in row] for row in a.tolist()],
                                    "equals_scaling_and_phase_shifts_of_witness": exact_gauge_scale})
        results.append(item)
        (HERE / "lift_spaces.json").write_text(json.dumps(results, indent=2)+"\n")
        print(cid, str(status), "spaces",len(spaces),"maximal",len(maximal),
              [(s["dimension"],s["equals_scaling_and_phase_shifts_of_witness"]) for s in item["spaces"]],flush=True)


if __name__ == "__main__":
    main()
