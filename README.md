# Finite T-data catalogue

An interactive catalogue of 385 indecomposable finite T-datum families in ranks
one through six, with exact matrices, admissible lifts, mutation loops,
exponents and reproducible classification certificates.

The scope is primitive positive diagonal symmetrizers and diagonal leading
matrix. There are 2, 10, 24, 66, 95 and 188 families in ranks one through six.
The 224 identity-symmetrizer families retain their IDs; 161 nonidentity
families have `s` IDs. The identity-symmetrizer rank-two completeness result
is published; the remaining classifications are computer-assisted and depend
on the included exact-arithmetic programs and solver certificates. A missing family name means that no identification
was found among the constructors and literature examples checked.

## Read the catalogue

**[Open the published catalogue](https://yuma-mizuno.github.io/finite-t-datum-catalog/).**

Download the repository and open [the catalogue](docs/catalogue/index.html)
in a browser. The HTML is self-contained and works offline. Keep the repository
layout to follow its links to proofs, data and the optional PDF editions.
The HTML source view on GitHub does not run the interactive document.

- [Catalogue documentation](docs/catalogue/README.md): controls, data format and checks.
- [Symmetrizable classification](research/symmetrizable/methods.html): all 385 families and proof obligations.
- [Methods and conventions](research/catalogue/methods.html): spectra and witnesses.
- [Ranks one and two](research/catalogue/lower-ranks-proof.html): proof and references.
- [Rank-three classification](docs/proofs/rank3-classification.pdf).
- [Rank-four classification](docs/proofs/rank4-classification.pdf).
- [Rank-five and rank-six classifications](research/higher_rank/methods.html).

## Reproduce and verify

Python 3.10+ with SymPy and Node.js 18+ run the catalogue checks:

```sh
python -m pip install sympy jsonschema
python docs/catalogue/test_records.py
node docs/catalogue/test_core.js
```

To rebuild the HTML and JSON records from the research data:

```sh
python tools/prepare_verification.py
python docs/catalogue/build_catalogue.py
```

Full computation instructions are in the
[rank-three report](research/rank3/CLASSIFICATION.md),
[rank-four reproduction guide](research/rank4/REPRODUCE.md),
[rank-five reproduction guide](research/higher_rank/rank5/README.md),
[rank-six reproduction guide](research/higher_rank/rank6/README.md), and
[family and spectral calculations](research/catalogue/README.md).
The [symmetrizable reproduction guide](research/symmetrizable/README.md) covers
the weighted classification, compressed inputs and query archive parts.
These guides describe the additional Z3, C++ and SageMath dependencies and distinguish
verification of stored certificates from rerunning the searches.

For browser checks, install Playwright and Chromium, then run
`node docs/catalogue/test_browser.js`. The document is tested with networking
disabled. Browser screenshots and downloaded test artifacts are ignored by Git.

## Website

This repository is prepared for static hosting. The root `index.html` opens the
catalogue, and `.nojekyll` allows plain static publishing. On GitHub Pages,
publish the root of `main` so that links to `research/` and `docs/proofs/` resolve.

## Sources and provenance

The records contain SHA-256 hashes of their source files and SMT queries.
The initial imported records also identify revisions of the original research
repository. Those revision identifiers document their origin; the original
Git history is not included here. Subsequent rebuilds record revisions of
this repository. See [source references](REFERENCES.md).

## Repository structure

- `docs/catalogue/`: the offline reader, presentation sources, generated records and UI tests.
- `docs/proofs/`: PDF editions of the classification documents.
- `research/catalogue/`: shared family, spectrum and quiver calculations.
- `research/rank3/` and `research/rank4/`: rank-specific engines and proof sources.
- `research/higher_rank/`: the shared identity-symmetrizer engine and rank-five/six certificates.
- `research/symmetrizable/`: weighted classifications, shared programs and certificates for ranks two through six.
- `tools/`: lossless data restoration and source-archive maintenance.

The `research/` packages are the reproducibility evidence for this catalogue.
Separate research notes, exploratory work and local working files are maintained
in a different repository. Presentation changes belong here. Preserve the
mathematical package paths when rebuilding archived certificates.

Large text inputs are distributed with gzip; SMT queries use independently
readable ZIP parts. `distribution.json` and the query manifests retain the
original and distributed hashes. Run `python tools/prepare_verification.py`
before executing computational programs that read unpacked inputs. Restored
files are ignored by Git. To refresh the distributed source archive after
editing verification programs, run `python tools/archive_sources.py`.
