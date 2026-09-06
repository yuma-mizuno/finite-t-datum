# Reproducing the rank-four classification

The readable report is `docs/proofs/rank4-classification.pdf`, relative to the
repository root. Its theorem is
restricted to identity symmetrizer and diagonal leading matrix.

The final result is 37 indecomposable families. There are 4,865 constant pairs,
4,405 parity exclusions, two further chain exclusions, 421 remaining unbounded
exclusions, and 37 complete four-dimensional lift spaces. Every representative
has exact two-sided reddening and labelled-period certificates.

Run commands from the repository root. Python 3.10.6, SymPy 1.14.0, NumPy 2.2.6,
and Z3 5.1.0 were used. The existing rank-three folder supplies the mutation
routines. Install the Z3 package as shown below.

```powershell
python -m pip install --target research/rank3/_deps z3-solver==5.1.0.0
wsl.exe --exec sh -c 'mkdir -p research/rank4/bin && g++ -O3 -std=c++17 research/rank4/enumerate_constants.cpp -o research/rank4/bin/enumerate_constants'
wsl.exe --exec research/rank4/bin/enumerate_constants 3 0 1 research/rank4/rank3_control.txt
wsl.exe --exec research/rank4/bin/enumerate_constants 4 0 1 research/rank4/constant_keys.txt
python research/rank4/prepare_constants.py
python research/rank4/run_lifts.py --workers 6 --timeout 2000
python research/rank4/run_lifts.py --workers 6 --timeout 60000 --retry-unknown
python research/rank4/run_lifts.py --workers 2 --timeout 180000 --retry-unknown --arithmetic rational
python research/rank4/classify_families.py --workers 6
python research/rank4/degree_two.py
python research/rank4/rank4_slices.py
python research/rank4/verify_rank4.py --workers 6 --timeout 180000
python research/rank4/package_certificates.py
python research/rank4/make_pdf_tables.py
python research/rank4/audit_delivery.py
```

The searches resume from existing JSONL files. For a completely fresh run, use
a separate checkout and remove only the generated `lift_feasibility.jsonl`
and `families.jsonl` there before running the searches. UNKNOWN results are
retained; a timeout never counts as an exclusion. Solver runtimes depend on
the machine and solver state. The verifier retries unresolved saved formulas
in a fresh process and refuses to complete unless every result is UNSAT.

The archived SMT formulas can be replayed without rerunning the enumeration:

```powershell
python -m zipfile -e research/rank4/smt_queries.zip research/rank4/smt_queries
python research/rank4/check_smt.py research/rank4/smt_queries/1117_exclusion.smt2 --timeout 180000
```

`verification.json` records the exact SHA-256 of each uncompressed formula.
`package_certificates.py` checks all hashes both before and after archiving.
Formulas for elementary exclusions contain `false`, justified by the parity
and chain lemmas in the PDF. The other exclusions and every coverage query
encode the full support and symplectic constraints. The archive is the tracked
copy; the extracted `smt_queries/` folder is ignored by Git.

Build the PDF from `research/rank4/`:

```powershell
latexmk -lualatex -synctex=-1 -interaction=nonstopmode -halt-on-error '-outdir=../../docs/proofs' rank4-classification.tex
```
