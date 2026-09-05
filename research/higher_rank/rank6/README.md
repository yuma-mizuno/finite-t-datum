# Rank-six finite T-data

The completed classification contains 108 indecomposable families with identity
symmetrizer and diagonal leading matrix. The [proof reductions](../methods.html)
explain the finite constant search, hereditary exclusions, finite pulse bound,
unbounded lift equations, periodicity and slice comparisons.
[Open rank six in the catalogue](../../../docs/catalogue/index.html#r6-c01/matrices).

All 5,040 upper-triangle tasks are complete. The 32,087 retained necessary
constant pairs include earlier unpruned checkpoints; this count is **after
necessary hereditary exclusions**, not the size of the unpruned constant set.
There are 31,979 excluded pairs and 108 complete six-dimensional lift spaces.
All exclusion and coverage checks have successful saved replays. The 45 finite
pulse exclusions have an additional independent C++ verification of 2,887,440
pulse positions.

## Files

| File | Content |
| --- | --- |
| `constant_tasks.zip` | Complete upper-triangle task outputs |
| `constant_candidates.json` | Retained necessary constants |
| `enumeration-verification.json` | Task coverage and lower-rank pruning controls |
| `lift_feasibility.jsonl`, `families.jsonl` | Lift outcomes and full matching spaces |
| `verification.jsonl` | All 32,087 replay outcomes, formula hashes and analytic witnesses |
| `smt_queries-01.zip`, `smt_queries-02.zip`, `smt_queries-03.zip` | Independently readable query archives |
| `smt_queries_manifest.json` | Complete member index and archive hashes |
| `pulse-verification.json`, `pulse_replays/` | Independent pulse inputs and results |
| `slice_signatures.json`, `slice_tasks/`, `slice_cpp/` | Periodic connected slices and exact orbit comparisons |
| `family-notes.json`, `quiver-data.json` | Named-family identities and mutation paths |
| `spectral-data.json`, `spectral-verification.json` | Root multiplicities and independent checks |
| `enrichment-verification.json` | Matrix, quiver and Jacobian verification |
| `base-records.json`, `catalogue-records.json`, `plots/` | Records and standalone exponent plots |

Each ZIP contains `.smt2.gz` members. A record's `provenance.query_path` points
to the ZIP containing its coverage formula. Every original compressed member
is preserved byte for byte. The original single ZIP digest is retained as
`original_smt_queries_zip_sha256`; the distributed files have their own hashes
in `smt_query_archives` and the manifest. The single source ZIP is not required.

The computation ledger preserves attempts and runtimes, with machine-specific
executable paths normalized to portable command names. `computation-policy.json`
sets no wall-clock cutoff for rank six. Earlier timed attempts remain recorded.

## Verify the distributed data

Run from the repository root with Python 3.10+, SymPy and Node.js 18+:

```sh
python research/higher_rank/verify_pruning_controls.py
python research/higher_rank/verify_enumeration.py 6
python research/higher_rank/query_archives.py 6 --verify
python docs/catalogue/test_records.py
node docs/catalogue/test_core.js
```

The archive verifier checks every query hash and the complete member set. It
does not solve the queries again. An analytic exclusion may have a trivial
false formula; its mathematical proof is the lemma and witness recorded in the
verification file and the linked argument.

Browser checks require Playwright and Chromium:

```sh
node docs/catalogue/test_browser.js
node docs/catalogue/test_rank6_browser.js
```

## Reproduce the independent pulse check

The saved inputs are plain integer files. On Linux or WSL, compile the verifier
and replay all 45 cases:

```sh
mkdir -p research/higher_rank/bin
g++ -O3 -std=c++17 research/higher_rank/verify_pulses.cpp -o research/higher_rank/bin/verify_pulses
python3 research/higher_rank/verify_pulse_exclusions.py
```

The script rebuilds the inputs from the classified cores and the recorded
analytic witnesses, then checks every permitted pulse position. The C++ source
and the input files are distributed; local binaries are not.

## Rebuild the catalogue

```sh
python research/higher_rank/build_records.py 6
sage -python research/higher_rank/finish_enrichment.py 6
python docs/catalogue/build_catalogue.py
```

The Sage step verifies stored quiver and Jacobian certificates and renders
exponent plots. It requires SageMath, mpmath and Matplotlib. The query manifest
is used automatically to select the correct archive for each record.
For recomputing the constant enumeration and lifts, see the commands and exact
hypotheses in the [methods](../methods.html#pulse-reduction). Use a separate
checkout for fresh searches and retain the distributed certificates for
comparison. Timeouts remain unresolved outcomes.

After regenerating individual queries in `rank6/smt_queries/`, run
`python research/higher_rank/package_checkpoints.py 6 queries` to verify and
package them into size-limited ZIP files again.
