"""Replay the finite classification reduction and symbolic family certificates.

The verifier trusts Python/SymPy exact arithmetic and Z3's UNSAT answers;
this is a computer-assisted proof, not a Lean formalization.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import platform

import numpy as np
import sympy as sp

from classify_rank3 import DATA, certificate, polynomial_data, z
from enumerate_constants import pair_key, simultaneous_witness
from solve_lifts import HERE, build_problem, get_matrices, z3


def reduced_exponent(e, rref):
    vector = sp.Matrix([list(e)])
    for row in rref.tolist():
        pivot = next(i for i, x in enumerate(row) if x)
        vector -= vector[0, pivot] * sp.Matrix([row])
    return tuple(vector)


def main():
    constants = json.loads((HERE / "constant_candidates.json").read_text())
    feasibility = {x["id"]: x for x in map(json.loads, (HERE / "lift_feasibility.jsonl").read_text().splitlines())}
    spaces = {x["id"]: x for x in json.loads((HERE / "lift_spaces.json").read_text())}
    assert len(constants["candidates"]) == len(feasibility) == 180
    assert Counter(x["status"] for x in feasibility.values()) == {"sat": 16, "unsat": 164}
    queries = HERE / "smt_queries"
    queries.mkdir(exist_ok=True)
    manifest = []
    for candidate in constants["candidates"]:
        cid = candidate["id"]
        p, m = candidate["N_plus_1"], candidate["N_minus_1"]
        assert simultaneous_witness(p, m)
        ap, am = 2*sp.eye(3)-sp.Matrix(p), 2*sp.eye(3)-sp.Matrix(m)
        assert ap*am.T == am*ap.T
        assert all(a.extract(indices, indices).det() > 0 for a in (ap, am)
                   for indices in ((0,), (1,), (2,), (0,1), (0,2), (1,2), (0,1,2)))
        solver, variables, indices, equations = build_problem(candidate, 30000)
        item = feasibility[cid]
        if item["status"] == "sat":
            assert spaces[cid]["coverage_status"] == "unsat"
            assert len(spaces[cid]["spaces"]) == 1
            space = spaces[cid]["spaces"][0]
            rref = sp.Matrix([[sp.Rational(x) for x in row] for row in space["rref"]])
            assert rref.rank() == len(variables)-3
            gauge = []
            for species in (0, 1):
                v = [0]*len(variables)
                for (s, i, j), ks in indices.items():
                    for k in ks:
                        v[k] = int(i == species)-int(j == species)
                gauge.append(v)
            generators = sp.Matrix.hstack(*(sp.Matrix(v) for v in [item["values"]]+gauge))
            assert generators.rank() == 3
            assert rref*generators == sp.zeros(rref.rows, 3)
            # Independently prove the polynomial identity on the entire
            # linear family, rather than merely evaluating sample points.
            for left, right in equations:
                assert Counter(reduced_exponent(e, rref) for e in left) == Counter(reduced_exponent(e, rref) for e in right)
            ap, am = get_matrices(item["values"], indices)
            assert (2*sp.eye(3)-ap.subs(z, 1)).tolist() == p
            assert (2*sp.eye(3)-am.subs(z, 1)).tolist() == m
            fresh = certificate(ap, am, bound=100)
            assert json.loads(json.dumps(fresh)) == item["certificate"]
            assert fresh["positive"] and fresh["negative"]
            solver.add(z3.Not(z3.And(*[
                z3.Sum([z3.RealVal(str(c))*v for c, v in zip(row, variables)]) == 0
                for row in rref.tolist()])))
            query_name = f"{cid:03d}_no_lift_outside_family.smt2"
        else:
            query_name = f"{cid:03d}_no_lift.smt2"
        query = solver.to_smt2()
        (queries / query_name).write_bytes(query.encode("utf-8"))
        status = solver.check()
        assert status == z3.unsat, (cid, status)
        manifest.append({"id": cid, "file": query_name, "result": "unsat",
                         "sha256": hashlib.sha256(query.encode()).hexdigest()})
    # All seven historical examples must map into the enumerated universe.
    mapping = {}
    by_key = {pair_key(c["N_plus_1"], c["N_minus_1"]): c["id"] for c in constants["candidates"]}
    for datum in DATA:
        if datum.a_plus.rows != 3:
            continue
        n0 = sp.diag(*(1+z**sp.degree(datum.a_plus[i,i], z) for i in range(3)))
        p = (n0-datum.a_plus).subs(z,1).tolist()
        m = (n0-datum.a_minus).subs(z,1).tolist()
        cid = by_key[pair_key(p,m)]
        assert feasibility[cid]["status"] == "sat"
        mapping[datum.label] = cid
    # Sixth rank-two control absent from the earlier five-row table.
    t2p = sp.Matrix([[1+z*z,-z],[-z,1-z+z*z]])
    t2m = (1+z*z)*sp.eye(2)
    control = certificate(t2p,t2m)
    assert (control["positive"]["h"], control["negative"]["h"]) == (5,2)
    report = {"all_180_unbounded_queries": "unsat", "all_16_symbolic_families": "verified",
              "all_16_two_sided_reddening_certificates": "verified", "historical_rank3_mapping": mapping,
              "python": platform.python_version(), "sympy": sp.__version__, "numpy": np.__version__,
              "z3": z3.get_version_string(), "queries": manifest}
    (HERE / "verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print({k:v for k,v in report.items() if k != "queries"})


if __name__ == "__main__":
    main()
