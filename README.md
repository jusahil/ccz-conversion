# Exact CCZ Conversion under Conserved Binary Charges

Research software and exact data for the article of the same title. This package evaluates Boolean cubic polarization over GF(2), finite concentration certificates, incidence-network and cohomological examples, and decoder stability. The analytical theorems and their physical hypotheses are in the article; finite computations do not prove arbitrary-size statements. In particular, 121/128 is a certified upper bound on squared target fidelity, not an attained optimum.

## Quick start (Linux)

Reference environment: Python 3.13.5, GNU g++ 14.2.0 (C++17), Poppler (`pdfinfo`, `pdftoppm`). Install the compiler and Poppler via your system package manager, then from the extracted repository root run:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip check
python verification/check_manifest.py
python -m unittest discover -s tests -v
python -m examples.quick_start
python verification/reproduce.py --profile full --out reproduction/full
```

A fresh non-existing output directory is required. Use `--profile smoke` for the three signed-branch programs and `--profile core` for the 16 retained baseline entry points and table/plot verification. The full profile additionally runs seven independent mathematical checks; it currently comprises 23 scientific entry points, the table-semantic verifier, and rendering of two quantitative figures. `--jobs 2` limits simultaneous computations. All coordinates in binary masks are zero-indexed (bit 0 corresponds to the article's first coordinate). Do not disable assertions. The manifest detects missing, changed and unexpected files, not authenticity. The reference stdout comparisons deliberately use the documented environment; exact CSV and signed-branch comparisons supply additional numerical/semantic checks.

## Contents

| Path | Scientific purpose |
|---|---|
| `src/` | Cubic algebra, finite concentration and coherent physical checks; reader-facing interface and exact spectra |
| `certificates/` | Independent Gaussian-integer projected-amplitude and Pauli-product branch solvers; secondary Walsh Pauli check |
| `checks/` | Independent finite regressions, network, cohomology, determinant, robustness and table checks |
| `data/` | 12 exact reference tables plus the complete signed Pauli branch table |
| `scripts/` | Typed data extraction and quantitative figure generation |
| `verification/` | Strict integrity, execution, reference outputs and signed-record validation |
| `tests/` | Fast mathematical and negative tests; not a substitute for full execution |
| `figures/` | Two derived quantitative figure families, not additional experimental data |
| `examples/` | Minimal reproducible examples |
| `docs/RESULTS.md` | Article-to-code correspondence and interpretation |
| `.github/workflows/` | Linux continuous-reproduction configuration |
| `SHA256SUMS` | Snapshot integrity inventory |

## Reproducibility and scope

The two primary certificate solvers emit each `(a,b,sign)` outcome with nullity, denominator, and, where nullity is six, exact moment and rearrangement bound. A third Pauli-based routine independently checks the aggregate bins. Their common key space has 32,766 outcomes. A full run writes logs, structured JSON, exact CSV and regenerated figures without overwriting reference data. The program-level tests of graph, tensor and physical constructions have explicitly bounded domains. Seeded tests are not exhaustive. The finite graph count, finite cohomology checks and phase-noise formula checks are not hardware simulations.

The quantitative graphs are generated from the exact stored values or regenerated equivalents. The article's conceptual diagrams are in its LaTeX source. Plot renderability is tested; identical PDF metadata or pixel equality across platforms is not asserted. See `docs/RESULTS.md` for the meaning of each file.

## Availability and citation

Author: Sahil Imtiyaz. The research software is distributed under the MIT License (see `LICENSE`). The accompanying scientific article is a separate work. Cite the article and the versioned software release together; machine-readable software citation metadata is in `CITATION.cff`.

Repository: https://github.com/jusahil/ccz-conversion

The initial GitHub Actions reproduction run on commit `dacc2d92537966df3510853b97be4d2f77959257` completed successfully, including integrity verification, nine tests, full scientific reproduction and evidence upload. See [the workflow run](https://github.com/jusahil/ccz-conversion/actions/runs/44565234749). The planned first public release is v1.0.0; no Zenodo DOI is asserted until it is issued. Consult the article for operational assumptions and complete mathematical proofs.

After any intentional, reviewed repository changes, regenerate the manifest with `python scripts/update_manifest.py --confirm-reviewed-change` and rerun the tests. Do not regenerate it to conceal an unexplained integrity failure.
