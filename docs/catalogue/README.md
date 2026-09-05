# Finite T-data: structured, interactive document

Open **[index.html](index.html)** in a browser. It is a self-contained offline
document: no server, network request, package installation or external math
renderer is needed. A fragment such as `#r4-c19/mutation` addresses a particular
record and view. Keep it in this repository if you want the relative links to
research sources and the optional print editions to work.

The catalogue contains 2 rank-one, 6 rank-two, 16 rank-three and 37 rank-four
indecomposable families, with identity symmetrizer and diagonal leading
matrix. Rank two follows the six rows of Mizuno's published Table 1; the
rank-one argument is included. The higher ranks retain their original class
numbers and constant IDs.

- **Matrices:** exact polynomial representatives, switching between A and N,
  both A+(1) and A−(1), and the exact rational matrices A+(1)⁻¹A−(1) and
  A−(1)⁻¹A+(1).
- **Notes:** RSG, SG, Zamolodchikov and other identifications, with the exact
  changes of variables and source conventions. Quiver notes have exact
  mutation paths to Dynkin, exceptional or surface representatives, or a
  witness proving mutation-infinite type.
- **Exponents:** Mizuno's thesis determinant formula, certified root-of-unity
  multiplicities, an interactive unit-disk plot, an exact table, a standalone
  SVG figure and an exportable interval certificate. A separate Jacobian
  calculation independently checks all 61 spectra.
- **Admissible lifts:** rational time rescaling and species shifts, checked
  with reduced BigInt fractions; export the transformed datum.
- **Mutation loop:** step through the stored connected-slice mutation word
  and its terminal relabelling, or replay a mutation-class certificate;
  inspect vertices and the exchange matrix.
- **Certificates:** family equations, negative permutations, canonical slice
  encodings, source links, query metadata and SHA-256 hashes.
- **Classification:** expandable statements, reductions, computational steps
  and proof dependencies. The full manuscripts remain linked.

## Data as the durable record

`catalogue.json` holds the entire document's records and structured proof
outline. `records/` provides each datum separately, using stable IDs such as
`r4-c19`. `record.schema.json` specifies the base record format.

`datum.N_plus[i][j]` and `datum.N_minus[i][j]` store an entry as a list of
`[coefficient, exponent]` pairs, in increasing exponent order. Zero is `[]`.
For example `[[1,1],[2,3]]` means z + 2z³. Delays and these two matrices
determine A± completely. The two specializations at 1 are also stored and
independently checked against the constant enumeration.

All indices in JSON are **zero-based**. The reader displays species, vertex
and permutation indices from one. Time offsets remain zero-based. In the
negative-permutation arrays, entry i is the column containing −1 in row i
of the terminal C-matrix. An `old_to_new` permutation sends each old vertex
index to its new index.

The `family` field records the homogeneous equations R x = 0, their original
variable order, a representative exponent vector and the coverage query.
The browser applies r′ᵢ = λrᵢ and p′ = λp + sᵢ − sⱼ with the last shift fixed
at zero, and verifies integrality and the strict row support inequalities.
The theorem identifies these admissible transforms with the entire family.

Lift exports have `kind: "admissible-lift"` and refer to their parent record.
They use decimal strings for all polynomial integers to avoid loss beyond
JavaScript's safe integer range. This derived format is distinct from the
base-record schema. It does not inherit the representative's time period or
mutation certificate: a new period is not computed by the reader.

The stored verification status describes the source's computer-assisted
argument. The UI does not rerun the exhaustive enumeration or the SMT solver.
The connected-slice animation does recompute every displayed mutation and
relabel operation from its exact integer exchange matrix.

## Regeneration and checks

The source of truth is the committed JSON/JSONL under `research/rank3`,
`research/rank4` and `research/catalogue`. Do not edit the generated records
by hand. [Methods and conventions](../../research/catalogue/methods.html)
explains the new certificates; the [research package](../../research/catalogue/README.md)
gives their reproduction commands. The build uses
Python 3.10+ and SymPy; the core tests use Node.js 18+.

```text
python docs/catalogue/build_catalogue.py
python docs/catalogue/test_records.py
node docs/catalogue/test_core.js
```

The builder verifies support, sign disjointness, the full symplectic identity
and the constants at 1 using symbolic arithmetic. The tests check package
consistency, source and query hashes, all 61 slice-loop returns, the RREF
witnesses, permutations, rational lift examples and arithmetic edge cases.
Source provenance records the latest commit affecting the original research
directories, with a separate enrichment commit and source hashes. A
presentation-only commit does not change those references.

The presentation sources are `index.template.html`, `catalogue.css`, `app.js`,
`core.js` and `proofs.json`. `index.html`, `catalogue.json` and `records/*.json`
are generated and committed. No external runtime is needed to read them.

For browser regression checks, install Playwright (or point `PLAYWRIGHT_MODULE`
to an existing installation) and optionally set `BROWSER_EXECUTABLE` to Chrome
or Edge. Run `node docs/catalogue/test_browser.js`. The test opens the document
from disk with networking disabled and saves screenshots in ignored `.qa/`.
