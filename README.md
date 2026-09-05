# Finite T-data

An interactive catalogue of 224 indecomposable finite T-datum families in ranks
one through six, with exact matrices, admissible lifts, mutation loops,
exponents and reproducible classification certificates.

The scope is identity symmetrizer and diagonal leading matrix. The catalogue
contains 2 rank-one, 6 rank-two, 16 rank-three, 37 rank-four, 55 rank-five and 108 rank-six families.
Rank-two completeness is a published result; the higher-rank classifications
are computer-assisted and depend on the included exact-arithmetic programs
and solver certificates. A missing family name means that no identification
was found among the constructors and literature examples checked.

## Read the catalogue

Download the repository and open [the catalogue](docs/catalogue/index.html)
in a browser. The HTML is self-contained and works offline. Keep the repository
layout to follow its links to proofs, data and the optional PDF editions.
The HTML source view on GitHub does not run the interactive document.

- [Catalogue documentation](docs/catalogue/README.md): controls, data format and checks.
- [Methods and conventions](research/catalogue/methods.html): spectra and witnesses.
- [Ranks one and two](research/catalogue/lower-ranks-proof.html): proof and references.
- [Rank-three classification](output/pdf/rank3-classification.pdf).
- [Rank-four classification](output/pdf/rank4-classification.pdf).
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
python docs/catalogue/build_catalogue.py
```

Full computation instructions are in the
[rank-three report](research/rank3/CLASSIFICATION.md),
[rank-four reproduction guide](research/rank4/REPRODUCE.md),
[rank-five reproduction guide](research/higher_rank/rank5/README.md),
[rank-six reproduction guide](research/higher_rank/rank6/README.md), and
[family and spectral calculations](research/catalogue/README.md).
They describe the additional Z3, C++ and SageMath dependencies and distinguish
verification of stored certificates from rerunning the searches.

For browser checks, install Playwright and Chromium, then run
`node docs/catalogue/test_browser.js`. The document is tested with networking
disabled. Browser screenshots and downloaded test artifacts are ignored by Git.

## Website

This repository is prepared for static hosting. The root `index.html` opens the
catalogue, and `.nojekyll` allows plain static publishing. On GitHub Pages,
publish the root of `main` so that links to `research/` and `output/pdf/` resolve.
Hosting is configured separately; pushing this repository does not enable Pages.

## Sources and provenance

The records contain SHA-256 hashes of their source files and SMT queries.
The initial imported records also identify revisions of the original research
repository. Those revision identifiers document their origin; the original
Git history is not included here. Subsequent rebuilds record revisions of
this repository. See [source references](REFERENCES.md).
