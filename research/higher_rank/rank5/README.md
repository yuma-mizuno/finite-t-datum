# Rank-five finite T-data

The classification contains 55 indecomposable families with identity symmetrizer
and diagonal leading matrix. The [proof reductions](../methods.html) explain
the finite enumeration, unbounded lift checks and distinctness argument.
The [interactive catalogue](../../../docs/catalogue/index.html#r5-c01/matrices)
shows each representative and its certificates.

## Stored evidence

| File | Content |
| --- | --- |
| `constant_tasks.zip` | All 720 completed upper-triangle tasks |
| `constant_candidates.json` | All 4,697 parity-compatible constant pairs |
| `enumeration-verification.json` | Task coverage and rank-four control |
| `lift_feasibility.jsonl` | 55 feasible pairs and 4,642 excluded pairs |
| `families.jsonl` | The 55 five-dimensional lift spaces and periodicity certificates |
| `verification.jsonl` | Successful UNSAT replays and hashes for all 4,697 queries |
| `smt_queries.zip` | Gzip-compressed SMT-LIB2 formulas, one archive member per query |
| `slice_signatures.json` | Terminal cycle lengths and exhaustive comparisons within repeated-length groups |
| `slice_tasks/` | Individual connected-slice certificates |
| `spectral-data.json`, `spectral-verification.json` | Root multiplicities and an independent spectral check |
| `family-notes.json`, `quiver-data.json` | Exact named-family matches and mutation witnesses |
| `enrichment-verification.json` | Independent matrix, quiver and Jacobian checks |
| `base-records.json`, `catalogue-records.json` | Structured representative records |

`computation-runs.jsonl` preserves runtimes, retries and outcomes. The recorded
Sage executable paths are normalized to `sage`; this does not change the
computation results. Earlier failures remain in the ledger. Completeness is
determined by the final verification files, not by every attempt succeeding.

## Check the distributed catalogue

From the repository root, with Python 3.10+, SymPy and Node.js 18+:

```sh
python research/higher_rank/verify_enumeration.py 5
python docs/catalogue/test_records.py
node docs/catalogue/test_core.js
```

The enumeration check compares the retained 120 rank-four control tasks with
all 460 parity-compatible reference candidates. It then checks all 720 rank-five
tasks and the exact combined candidate set. The catalogue checks cover all
116 records. Browser checks additionally require Playwright and Chromium:

```sh
node docs/catalogue/test_browser.js
node docs/catalogue/test_rank5_browser.js
```

## Rebuild from stored classification data

```sh
python research/higher_rank/build_records.py 5
sage -python research/higher_rank/finish_enrichment.py 5
python docs/catalogue/build_catalogue.py
```

The first command checks representative polynomial identities. The second
replays the existing quiver and spectral certificates and renders the SVG
plots using SageMath, mpmath and Matplotlib. It uses the distributed named-family,
quiver and spectral data; it does not search for new witnesses.

## Rerun computations

The search programs are in `research/higher_rank/`. The constant enumerator
requires a C++17 compiler and runs with Linux Python (including WSL):

```sh
mkdir -p research/higher_rank/bin
g++ -O3 -std=c++17 research/higher_rank/enumerate_constants.cpp -o research/higher_rank/bin/enumerate_constants_pruned
python research/higher_rank/run_constants.py 5 --workers 8 --seconds 7200
```

The lift driver requires SymPy and Z3; the saved run used Z3 5.1.0. Its stages
are `lifts`, `families` and `verify`, with `--retry` for unresolved results.
`--arithmetic` selects integer or rational constraints; `--encoding` selects
multiplicity, sorting, transport or principal-subdatum constraints.

```sh
python -m pip install --target research/rank3/_deps z3-solver==5.1.0.0
python research/higher_rank/run_lifts.py 5 lifts --workers 8 --timeout 60000 --seconds 7200 --retry
python research/higher_rank/run_lifts.py 5 families --workers 8 --timeout 60000 --seconds 7200 --retry
python research/higher_rank/run_lifts.py 5 verify --workers 8 --timeout 60000 --seconds 7200 --retry
```

The time limit is cumulative across the recorded stages. Increase it when
needed; a timeout does not establish an exclusion. These commands resume saved
results. A full fresh search requires a separate checkout with fresh task,
lift, family and verification outputs and a fresh computation ledger. Keep the
distributed certificates in the original checkout for comparison. The commands
above are not a promise that every query finishes within the stated timeout.
