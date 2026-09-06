# Symmetrizable classification

The [classification document](methods.html) states the mathematical scope,
reductions and computer-assisted proof obligations. The complete catalogue
contains 385 families in ranks one through six, including the 224 families
with identity symmetrizer.

Each `rankN/` directory contains the constant candidates, lift spaces, query
replays, slice and spectral certificates, and records for the nonidentity
symmetrizers. Shared programs are kept here; the identity-symmetrizer engine
is reused from `../higher_rank/` and `../rank4/`.

Run these commands from the repository root:

```sh
python tools/prepare_verification.py
python research/symmetrizable/audit_completion.py
python research/symmetrizable/audit_constant_archives.py 6
```

The first command restores compressed JSON sources and archived control
inputs. Restored files are ignored by Git. `distribution.json` records both
the distributed and original byte hashes. The large rank-six SMT archive
is distributed as independent ZIP parts with a complete member manifest;
the audit reads them directly.

`source-manifest.json` and `certificate-sources.zip` preserve the distributed
verification sources and the original constant-enumerator revisions needed
by the saved tasks. The imported revision identifiers refer to the research
repository; its Git history is separate from this public repository.

The audits check saved finite certificates. They do not rerun the entire
constant search or every SMT solver call. Full reproduction commands and
additional SageMath, SymPy, Z3 and C++ dependencies are in the classification
document.
