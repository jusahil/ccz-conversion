# Scientific result index

The appendix names correspond to the integrated article. Assertions are enabled. `C` and `N` identifiers are stable program identifiers, not theorem numbers. Different routines intentionally use independent mathematical realizations. All computations are finite; expected values that are analytical inputs are distinguished from computed certificates.

| Main result / appendix | Entry points and datasets | Evidence type |
|---|---|---|
| Compression, witness and activation (A–C, J.2) | `src/algebra.py`, `src/extended.py`, `checks/foundational.py`, `checks/protocols.py`, `checks/independent_cubics.py`, `checks/relative_pencils.py`, `checks/relative_examples.py` | Cubic polarization, kernels, reference orientation, finite decoder checks; analytical completeness is in article |
| Exact frontier and concentration (D, E, J.5) | `src/new_certificate.py`, `checks/independent_concentration.py`, `checks/concentration_families.py`, `copy_ranks.csv`, `decoder_checks.csv`, `reference_nullities.csv`, `monomial_moments.csv`, `batch_yield.csv` | Contractions and decoder enumeration; yield table evaluates a proved formula |
| Seven-qubit separation (E.2, J.3) | `certificates/projected_amplitudes.cpp`, `certificates/pauli_products.cpp`, `certificates/pauli_walsh.cpp`, `signed_pauli_branches.csv`, `fano_histograms.csv` | Exhaustive finite certificate after analytical first-nullity-drop reduction; upper bound not achievable optimum |
| Incidence networks (F, H.3) | `checks/graph_census.py`, `checks/weighted_graphs.py`, `checks/tensor_graphs.py`, `checks/series_filtration.py`, `checks/network_examples.py`, `graph_counts.csv` | Graph enumeration, weighted seeded tests, support filtration and named examples; fixed primitive library |
| Cup and determinant resources (G, H, J.4) | `checks/cohomology.py`, `src/physical.py`, `src/new_certificate.py`, `checks/determinant_css.py`, `checks/relative_examples.py`, `mesh_counts.csv`, `determinant_planes.csv` | Finite determinant planes, exact cochain identities, full branch maps, logical parity checks; not a native-gate simulator |
| Robustness and resource accounting (I) | `checks/finite_controls.py`, `checks/robustness.py`, `cycle_noise.csv`, `reference_price.csv` | Exact formula checks in declared noise and abstract cost models; no fault-tolerance threshold |
| Fidelity ceilings (E, G, J) | `src/spectra.py`, `scripts/export_data.py`, `checks/table_semantics.py`, `fidelity_ceilings.csv` | Computed direct Walsh spectra, rearrangement and signed-branch bounds; A's attained value follows the analytic decoder |

## File and schema details

- `signed_pauli_branches.csv`: all nonidentity signed Paulis in lexicographic key order. `denominator` equals 128 plus the signed integer Pauli numerator, and success probability for that measurement outcome is denominator/256. Moment and bound columns are populated only for nullity six. All fractions are reduced.
- `fano_histograms.csv`: separate marginal histograms of nullity, fourth Pauli moment and sorted-spectrum bound. They must not be interpreted as a joint distribution.
- `determinant_planes.csv`: six-dimensional restricted phase, distinct three-planes, exact contraction histograms grouped by nullity.
- `copy_ranks.csv`: displacement ranks in each fully charged copy experiment; `decoder_checks.csv` counts different objects (outcome/coordinate assignments).
- `mesh_counts.csv`: explicit mesh cardinalities, first cohomology and cup/shift checks; not a general manifold classifier.
- `reference_nullities.csv`, `monomial_moments.csv`: exact finite diagnostics in the indicated reference/copy families.
- `batch_yield.csv`: reference dimension 1 represents both 1 and 2; `cycle_noise.csv`: exact rational probabilities and squared decoder fidelity; `reference_price.csv`: normalized abstract cost lines and minimum.
- `fidelity_ceilings.csv`: distinguish attained analytical value of A, computed bounds for B, the trace phase and determinant. Different target counts are explicit.
- `verification/programs.json`: program identity, source checksum and reference-comparison contract. `verification/reference_stdout/` retains baseline scientific output (not a substitute for execution).
- Programs under `checks/` and the computational scripts in `src/` are intended to execute as standalone entry points through `verification/reproduce.py`; only the small `src.ccz_research` interface is advertised as a reusable API.
