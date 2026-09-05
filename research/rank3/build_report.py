"""Render the proved reduction and exact computed representatives as Markdown."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from classify_rank3 import z
from solve_lifts import HERE


INTRO = r"""# Rank-three finite T-data: a computer-assisted classification

The classification below has **16 indecomposable classes** in the precise
scope of [Mizuno, arXiv:2301.13239v2](https://arxiv.org/html/2301.13239v2):
the symmetrizer is the identity, and the leading matrix is diagonal.
The proof combines an elementary, delay-independent reduction with exhaustive
integer computations. The supporting programs and exact certificates are included with the classification.

## 1. Scope and statement

Let
\[
 N_0(z)=\operatorname{diag}(1+z^{r_1},1+z^{r_2},1+z^{r_3}),\qquad
 A_\pm(z)=N_0(z)-N_\pm(z),
\]
where \(r_i\) are positive integers and \(N_\pm\) have nonnegative integer
coefficients. A coefficient in entry \((i,j)\) is supported at an integer
\(p\) satisfying \(0<p<r_i\); the two signs have disjoint support in each
entry. Require
\[
 A_+(z)A_-(z^{-1})^{\mathsf T}
 =A_-(z)A_+(z^{-1})^{\mathsf T}.
\]
Finite type means periodicity of the associated universal Y-system, equivalently
the finite-type condition in Definition 1.4 of the cited paper.

**Classification.** An indecomposable pair in this scope is of finite type if
and only if, after a simultaneous permutation of the three indices and possibly
exchanging the two signs, it has the form
\[
 A_\pm(z)=E(z)B^{(k)}_\pm(z^\lambda)E(z)^{-1},\qquad
 E(z)=\operatorname{diag}(z^{s_1},z^{s_2},1),                 \tag{1}
\]
for exactly one row \(k\) of the table below. Here \(\lambda\in\mathbb Q_{>0}\)
and \(s_1,s_2\in\mathbb Q\), and the parameters must make all the displayed
exponents integral and satisfy the support inequalities above. The expression
is interpreted first in the group algebra of rational powers of \(z\).

More explicitly, if a representative has delay \(r_i^{(k)}\) and a term
\(z^p\) in entry \((i,j)\), its transformed delay and exponent are
\[
 r_i=\lambda r_i^{(k)},\qquad
 p'=\lambda p+s_i-s_j,\qquad s_3=0.
\]
Every such admissible transformation gives a finite-type datum. There are no
additional polynomial lift families, including at arbitrarily large delays.

**Corollary.** For an indecomposable rank-three datum in this scope, finite
type is equivalent to simultaneous positivity: there exists a row vector
\(v>0\) such that \(vA_+(1)>0\) and \(vA_-(1)>0\). Necessity is the cited
positivity theorem; sufficiency follows from the exhaustive reduction and
the periodicity certificates below.

These 16 families also give 16 classes under the paper's change-of-slices
equivalence; the argument and the distinctness check are in Section 6.
Decomposable rank-three data are direct sums with rank partitions
\(3=2+1\) or \(3=1+1+1\), using the classified rank-two and rank-one data.
The theorem here does not cover nontrivial symmetrizers or non-diagonal \(N_0\),
which belong to the more general definition of T-data.

## 2. Explicit representatives

In each row,
\[
 B^{(k)}_\pm(z)=\operatorname{diag}(1+z^{r_1},1+z^{r_2},1+z^{r_3})-N_\pm(z).
\]
The notation \(ij:f\) in a table cell means that entry \((i,j)\) of
\(N_\pm\) is \(f\). Unlisted entries vanish. The indices are one-based.
The integer ID refers to the exhaustive constant-candidate enumeration.

"""

PROOF = r"""
## 3. Why the coefficient search is exhaustive

Write \(P=N_+(1)\) and \(M=N_-(1)\).
The simultaneous positivity theorem, recalled as Lemma 3.3 in the rank-two
paper and valid in arbitrary rank, gives a row vector \(v>0\) satisfying
\[
 v(2I-P)>0,\qquad v(2I-M)>0.                              \tag{2}
\]
Consequently both diagonal entries \(2-P_{ii}\) and \(2-M_{ii}\) are
positive integers, so \(P_{ii},M_{ii}\in\{0,1\}\).
Both Z-matrices \(2I-P\), \(2I-M\) are nonsingular M-matrices; in particular,
all their principal minors are positive. This follows directly by diagonal
scaling using \(v\), or by the usual M-matrix criterion.

**Strong connectivity lemma.** The directed graph with an arrow \(i\to j\)
when \(i\ne j\) and \(P_{ij}+M_{ij}>0\) is strongly connected for every
indecomposable datum in this scope, even without finite type.

To prove this, suppose a cut makes both polynomial matrices block lower
triangular. Their constant terms are identity matrices, so
\(K=A_+^{-1}A_-\) exists over \(\mathbb Q(z)\) and is regular at zero.
The symplectic identity gives
\[
 K(z)=K(z^{-1})^{\mathsf T}.
\]
Thus this block lower triangular matrix is block diagonal. If there were
nonzero entries crossing the cut, let \(h>0\) be their smallest polynomial
degree. Since
\[
 K=I+A_+^{-1}(N_+-N_-),
\]
the coefficient of degree \(h\) in the crossing block is the corresponding
coefficient of \(N_+-N_-\). Contributions involving a nonconstant factor of
\(A_+^{-1}\) have larger degree: a crossing factor already has degree at
least \(h\). Disjoint support makes the degree-\(h\) crossing coefficient
nonzero, a contradiction. A reducible directed graph therefore has no edges
between its strongly connected components, and the datum is decomposable.

**Cycle bound.** Put \(w_{ij}=\max(P_{ij},M_{ij})\) for \(i\ne j\).
From (2), every nonzero \(w_{ij}\) satisfies
\[
 w_{ij}v_i<2v_j.
\]
Every arrow in a strongly connected three-vertex graph belongs to a simple
directed cycle of length two or three. Multiplying along that cycle gives
\[
 w_{ij}w_{ji}<4,\qquad
 w_{12}w_{23}w_{31}<8,\qquad w_{13}w_{32}w_{21}<8.          \tag{3}
\]
In particular \(0\le P_{ij},M_{ij}\le7\) for \(i\ne j\).
This bound concerns coefficient sums, not a delay cutoff.

The program `enumerate_constants.py` exhausts the resulting finite set.
It first enumerates the six off-diagonal maxima \(w_{ij}\), then every
ordered split \((P_{ij},M_{ij})\) with that maximum, and all 64 diagonal
choices. It imposes
\[
 (2I-P)(2I-M)^{\mathsf T}=(2I-M)(2I-P)^{\mathsf T},
\]
all principal-minor inequalities, and exact simultaneous positivity.
All arithmetic in this enumeration is integer arithmetic. The counts are:

| Stage | Count |
|---|---:|
| Strongly connected off-diagonal maxima satisfying (3) | 1,134 |
| Ordered off-diagonal splits | 1,474,290 |
| Labelled pairs after symplectic and principal-minor tests | 1,835 |
| Orbits under index permutations and sign exchange | 180 |
| Orbits also satisfying simultaneous positivity | 180 |

For the last test, the program constructs an explicit positive integer vector
\(v\) for each retained pair. Feasibility is checked by enumerating the
extreme rays of the cone given by the nine homogeneous weak inequalities
\(v\ge0\), \(v(2I-P)\ge0\), \(v(2I-M)\ge0\), and testing their sum for
strict positivity. In dimension three each extreme ray is obtained from the
cross product of two constraint normals. The cone is pointed because it lies
in the nonnegative orthant. This is an exact feasibility test.

## 4. Unbounded classification of polynomial lifts

For each of the 180 constant pairs, introduce an integer variable for every
delay and for every unit of coefficient mass in every entry. Several units
of the same sign may have equal exponents; thus coefficients greater than
one are included. Exponents of opposite signs in the same entry must differ.
Order the exponents within each sign and entry to remove duplicate labellings.

Expand each entry of
\(A_+A_-^*-A_-A_+^*\), where \(F^*(z)=F(z^{-1})^{\mathsf T}\).
The resulting expression is a difference of two finite multisets of Laurent
monomials, and each exponent is an integral linear form in the delay variables.
Polynomial equality is exactly equality of the two exponent multisets.
The solver encodes the equality of multiplicities at every exponent appearing
on the positive side; equality of total multiplicities is already imposed by
the constant-matrix identity. Thus the encoding is exact even when monomials
coalesce.

There is **no upper bound on any delay or exponent** in these formulas.
They are formulas in linear integer arithmetic with Boolean combinations of
equalities. The results are:

* 164 constant pairs have no polynomial lift: the unbounded formula is UNSAT.
* The 16 constant pairs in Section 2 have polynomial lifts.

For each survivor, `classify_lift_spaces.py` pairs equal exponents in one
satisfying solution. The matching gives a homogeneous rational linear
subspace on which the entire Laurent identity holds. The program then asks
whether any admissible solution lies outside that subspace. For all 16
survivors the answer is UNSAT, and the subspace has dimension exactly three.
Its three independent generators are the listed representative's exponent
vector and the two species shifts \(s_1,s_2\) in (1).

`verify_classification.py` separately reduces every exponent modulo the
row-reduced linear equations and checks equality of the resulting monomial
multisets. Hence the identities hold symbolically throughout each family.
The verifier also replays all 164 exclusion queries and all 16
no-solution-outside-the-family queries. Their complete SMT-LIB2 formulas and
SHA-256 hashes are saved, with every result equal to UNSAT.

The exhaustive necessity argument is now complete: finite type implies one
of 180 constant pairs, and an admissible symplectic polynomial lift of such
a pair belongs to exactly one of the 16 families (1).

## 5. Sufficiency and periods

The programs reconstruct the exchange matrix from equation (2.2) of the
rank-two paper. For a delay tuple \(\boldsymbol r\), the vertices are
\(R=\{(i,p):0\le p<r_i\}\). One time step mutates the three vertices
\((i,0)\), then relabels \((i,p)\mapsto(i,p-1\bmod r_i)\).
The implementation verifies skew-symmetry, commutativity of the three
mutations, return of the exchange matrix after every step, and inverse-step
consistency.

Starting from the identity principal C-matrix, exact integer mutation gives
a negative permutation matrix in both time directions for every representative.
The lengths are \(h_+,h_-\) in the table. Proposition 2.7 of the cited paper
therefore proves finite type for every representative in every semifield.
The independently detected return to the identity C-matrix gives the labelled
seed periods in the table; these are not periods inferred from numerical
specializations. The five historical rank-two examples and the sixth
rank-two tadpole case reproduce the rank-two reddening lengths as controls.

Finite type is preserved throughout (1). For integral shifts, the change is
the bijection of solutions \(Y'_i(u)=Y_i(u-s_i)\). Replacing \(z\) by
\(z^m\), for a positive integer \(m\), splits the system into \(m\)
independent residue classes of time, so preserves finite type in both
directions. Given rational \(\lambda,s_1,s_2\), choose a positive integer
\(L\) clearing their denominators. Then the transformed datum at \(z^L\)
is an integral species shift of the representative at \(z^{L\lambda}\).
The integral argument applies to both sides. This proves sufficiency for
every admissible member of every family, with no computation at a delay bound.

## 6. Equivalences and distinctness

Formula (1) is already a literal classification of polynomial pairs; it does
not require an informal interpretation of equivalence. The 16 constant pairs
are distinct under permutations and sign exchange, and specialization at
\(z=1\) is invariant under (1). Thus this literal parameterization has no
overlap between its 16 rows.

For comparison with Section 3.1 of the paper, integral time dilations repeat
the same slice sequence. Species shifts change the placement of the mutation
events and the initial cut in that sequence. They can be connected by a
straight path of admissible rational shifts: the support conditions are
strict linear inequalities, hence define a convex domain. Along this path,
the order of events changes only when their times coincide. After clearing
denominators, such a coincidence is a simultaneous mutation in a valid
T-datum and the relevant vertices have zero exchange entry. Its two orders
therefore commute. Moving the initial cut rotates the cyclic sequence.
Consequently admissible shifts and rational rescalings give change of slices
in the paper's sense.

Distinctness under that equivalence is checked directly, without treating
\(N_\pm(1)\) as an assumed invariant of arbitrary slice changes.
For each representative, `slice_equivalence.py` follows one connected
component through the cyclic permutation of components. Its returning mutation
loop has exactly three mutation events and a terminal vertex permutation.
Each event roots one orbit of the terminal permutation, and these three orbits
cover the vertices. Ordering each orbit from its root gives a canonical
labelling of every vertex. The loop is therefore encoded completely by its
three orbit lengths and the exchange matrix in this labelling.

The program exhausts cyclic rotations and adjacent commuting mutations until
closure. If a loop is represented by \((B,(a,b,c),\pi)\), rotation gives
\((\mu_a B,(b,c,\pi^{-1}a),\pi)\); adjacent events can be exchanged precisely
when their exchange entry at that point vanishes. Sign exchange is handled
by \(B\mapsto-B\) with the mutation word and permutation unchanged.
The sixteen sets of canonical encodings are disjoint. The full minimal
encoding, not merely its hash, is saved for every row.

## 7. Reproduction and proof status

The programs and certificates are in `research/rank3/`.
The checked runtime is Python 3.10.6, SymPy 1.14.0, NumPy 2.2.6,
and Z3 5.1.0. Z3 is installed locally under the ignored `_deps/` directory.

From the repository root, run:

```powershell
python -m pip install --target research\rank3\_deps z3-solver==5.1.0.0
python research\rank3\classify_rank3.py
python research\rank3\enumerate_constants.py
python research\rank3\solve_lifts.py --timeout 2000 --optimize --dynamics
python research\rank3\classify_lift_spaces.py
python research\rank3\verify_classification.py
python research\rank3\slice_equivalence.py
python research\rank3\build_report.py
```

A solver timeout is reported as UNKNOWN and cannot pass the verifier as an
exclusion. The successful run had no UNKNOWN results. The mathematical
classification depends on the exact-arithmetic programs and the correctness
of Z3's UNSAT results. It is not a proof-assistant formalization. The outputs are exact certificates for the finite enumeration and the
unbounded feasibility and coverage checks.

Files:

* [Constant candidates](constant_candidates.json): all 180 pairs and positive vectors.
* [Polynomial witnesses and periodicity certificates](lift_feasibility.jsonl): all 180 solver outcomes and the 16 explicit witnesses.
* [Complete lift spaces](lift_spaces.json): exact row-reduced equations for all 16 families.
* [Verification and query hashes](verification.json): the 180 replayed UNSAT answers and runtime versions.
* [SMT-LIB2 queries](smt_queries/001_no_lift_outside_family.smt2): the directory contains all 180 unbounded formulas.
* [Slice signatures](slice_signatures.json): complete encodings proving the 16 classes distinct.
* [Earlier paper](https://arxiv.org/abs/1912.05710): Table 3 supplies seven examples, all recovered here.

The seven earlier examples correspond, in their original order, to constant
IDs 1, 75, 47, 24, 2, 4, and 21. The other nine rows are additional to that
particular table; this comparison is not a claim that those nine systems
are new to the literature.
"""


def entries(a, delays):
    n0 = sp.diag(*(1+z**r for r in delays))
    matrix = n0-sp.Matrix([[sp.sympify(x, locals={"z":z}) for x in row] for row in a])
    result = [rf"{i+1}{j+1}: {sp.latex(matrix[i,j])}"
              for i in range(3) for j in range(3) if matrix[i,j] != 0]
    return r"\(" + (r";\quad ".join(result) if result else "0") + r"\)"


def main():
    data = [x for x in map(json.loads, (HERE / "lift_feasibility.jsonl").read_text().splitlines()) if x["status"] == "sat"]
    slices = {x["id"]:x for x in json.loads((HERE / "slice_signatures.json").read_text())}
    lines = [INTRO, "| Class | ID | Delays | Nonzero entries of \\(N_+\\) | Nonzero entries of \\(N_-\\) |\n|---|---:|---|---|---|\n"]
    for k, item in enumerate(data, 1):
        delays = item["values"][:3]
        lines.append(f"| {k} | {item['id']} | \\({tuple(delays)}\\) | {entries(item['A_plus'],delays)} | {entries(item['A_minus'],delays)} |\n")
    lines.append("\nThe following are exact labelled seed periods for these representatives.\nThe reddening lengths are measured in their displayed time coordinates.\n\n| Class | \\(h_+\\) | \\(h_-\\) | Labelled period | Vertices per slice |\n|---|---:|---:|---:|---:|\n")
    for k, item in enumerate(data, 1):
        cert = item["certificate"]
        lines.append(f"| {k} | {cert['positive']['h']} | {cert['negative']['h']} | {cert['labelled_tropical_seed_period']} | {slices[item['id']]['slice_vertices']} |\n")
    lines.append(PROOF)
    (HERE / "CLASSIFICATION.md").write_text("".join(lines),encoding="utf-8")
    print("Wrote", HERE / "CLASSIFICATION.md")


if __name__ == "__main__":
    main()
