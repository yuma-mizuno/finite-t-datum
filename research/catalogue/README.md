# Catalogue enrichment: ranks 1–4

The readable entry point is [the interactive catalogue](../../docs/catalogue/index.html).
[Methods and conventions](methods.html) gives a rendered account of the spectral
calculation and its certificates. The package contains 61 indecomposable records:
2 in rank one, 6 in rank two, 16 in rank three and 37 in rank four. The scope is
identity symmetrizer and diagonal leading matrix, modulo admissible rescaling,
species shifts, index permutations and exchange of signs.

## Sources and mathematical status

- `lower_ranks.py` supplies the scalar representatives and the six rows of
  [Mizuno's Table 1](https://arxiv.org/html/2301.13239v2#S1.T1). It reuses the
  exact principal-coefficient mutation engine in `research/rank3` to recompute
  their reddening and period certificates. Rank-one completeness has an
  [elementary proof](lower-ranks-proof.html); rank-two completeness is the
  published Theorem 1.5 and its lift analysis in Lemmas 3.6–3.9.
- The rank-three and rank-four classification certificates remain in their
  original research directories. This enrichment does not change them.
- `family_notes.py` checks polynomial identities, including the recorded sign
  exchange, permutation, rational dilation and species shifts. It searches the
  RSG constructors of the relevant rank, SG constructors of rank at most four,
  finite ADE/tadpole tensor pairs and the seven rank-three examples in Mizuno's
  arXiv:1912.05710 Table 2. Degree-two pairs are also checked for the commuting
  Cartan construction. A missing name means only that these checks found no
  identification. Literature notes retain the author's sign and slice conventions.
- `quiver_types.py` and `quiver_more.py` save mutation paths to Dynkin trees,
  standard exceptional representatives, triangulated punctured spheres, or
  a quiver containing at least three parallel arrows. All 61 labels have an
  exact witness. Surface certificates include oriented faces, edge order and
  the vertex bijection. The three sphere cases are `r3-c05`, `r4-c32` (four
  punctures), and `r4-c36` (five punctures). Mutation-infinite labels use the
  three-arrow criterion for connected skew-symmetric quivers on at least
  three vertices. These labels concern quivers, while finite T-type concerns
  a particular periodic mutation loop.
- `exponents.py` implements Mizuno's thesis (3.2.1), (3.4.1), Definition 3.4.4
  and Theorem 3.4.5, with the Cartan-like hypothesis checked after centering.
  It certifies the fixed point with interval arithmetic and determines exact
  integer multiplicities on the period grid. `verify_enrichment.py` independently
  derives the spectrum from the logarithmic mutation Jacobian (3.4.2), replays
  every quiver witness, and verifies both rational matrix products.

The resulting `enrichment-verification.json` records successful independent
checks for all 61 entries. These are reproducible computer-assisted certificates,
not proof-assistant formalizations.

## Conventions

All JSON indices are zero-based. The reader displays species and vertices from
one. Rational matrix entries and gauge parameters are strings, such as `3/2`.
The matrices recorded are exactly `A_plus(1)^(-1) A_minus(1)` and its inverse.

RSG parameter lists use the generation order of the archived constructor.
Thus `RSG(3,1)` corresponds to the continued fraction `3/4 = [1,3]` in Remark
1.7, and `RSG(2,1,1)` to `3/5 = [1,1,2]`. Do not read a constructor list as
the paper's continued-fraction list in the same order.

`sources/rsg_constructor_snapshot.py` contains the archived continued-fraction constructors
used for the family identifications. Its SHA-256 is
`d8068387034e7ceba2509983f6b9a034aa4d7c761ae1f898168e30ca9d9d0689`.
The standalone transcription in `family_notes.py` treats absent predecessor
families as zero; it guards the archived source's negative predecessor index
at the first generation. Every generated candidate is checked for the full
symplectic identity and every recorded match is an exact polynomial identity.
The mathematical source is [Nakanishi–Stella, Wonder of sine-Gordon
Y-systems](https://arxiv.org/abs/1212.6853).

`sources/mizuno-thesis.pdf` is the March 2021 thesis *Difference equations
arising from cluster algebras*, downloaded from
[the author's site](https://yuma-mizuno.github.io/thesis.pdf).
Printed pages 55–58
are PDF pages 56–59 because of the cover. The catalogue records the PDF's
SHA-256 as source provenance.

## Reproduce

Run from the repository root. Python with SymPy suffices for the first two
commands; the subsequent calculations use SageMath 10.8, mpmath and Matplotlib.
The pre-existing rank-three and rank-four records in `docs/catalogue/catalogue.json`
provide the polynomial and exchange inputs; enrichment fields are ignored.

```text
python research/catalogue/lower_ranks.py
python research/catalogue/family_notes.py
sage -python research/catalogue/exponents.py
sage -python research/catalogue/quiver_types.py
sage -python research/catalogue/quiver_more.py
sage -python research/catalogue/verify_enrichment.py
sage -python research/catalogue/plot_exponents.py
python docs/catalogue/build_catalogue.py
python docs/catalogue/test_records.py
node docs/catalogue/test_core.js
node docs/catalogue/test_browser.js
```

The quiver search scripts seek witnesses; verification replays the saved paths
without relying on the search or on Sage's type-recognition verdict. The SVG
plots are standalone Matplotlib figures suitable for export. The interactive
plots use the same exact exponent data and can be inspected without a server.

Primary references for quiver notes are
[Felikson–Shapiro–Tumarkin](https://arxiv.org/abs/0811.1703),
[Fomin–Shapiro–Thurston](https://arxiv.org/abs/math/0608367), and
[SageMath's quiver documentation](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/cluster_algebra_quiver/quiver.html).
