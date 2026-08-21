# Changelog

## coda.base 1.0.7 (Release date: 2026-05-12)

- [`dist_coda()`](https://mcomas.net/coda.base/reference/dist_coda.md)
  now accepts `method = "L2"` as an alias for the Euclidean distance.
- The partial principal-balance search validates grouped inputs more
  thoroughly and improves the handling of tabu-search neighbourhoods.
- Additional tests cover distance aliases, basis construction and
  partial tabu-search edge cases.
- The pkgdown site now uses Bootstrap 5 and includes a redesigned logo,
  favicons and a revised Get Started vignette.

## coda.base 1.0.6 (Release date: 2026-05-08)

CRAN release: 2026-05-08

- Added the standard compositional operations
  [`closure()`](https://mcomas.net/coda.base/reference/closure.md),
  [`perturbation()`](https://mcomas.net/coda.base/reference/perturbation.md)
  and
  [`powering()`](https://mcomas.net/coda.base/reference/powering.md).
- Added
  [`dist_coda()`](https://mcomas.net/coda.base/reference/dist_coda.md)
  for Aitchison, angular, L1 and Euclidean distances between
  compositions.
- Added exact, constrained and tabu-search methods for partial principal
  balances:
  [`partial_pb_exact()`](https://mcomas.net/coda.base/reference/partial_pb_exact.md),
  [`partial_pb_constrained()`](https://mcomas.net/coda.base/reference/partial_pb_constrained.md)
  and
  [`partial_pb_tabu_search()`](https://mcomas.net/coda.base/reference/partial_pb_tabu_search.md).
- Added
  [`pb_tabu_search()`](https://mcomas.net/coda.base/reference/pb_tabu_search.md)
  to construct a sequential binary partition by tabu search, together
  with
  [`random_composition_with_fixed_pb()`](https://mcomas.net/coda.base/reference/random_composition_with_fixed_pb.md)
  for simulation.
- Extended
  [`pb_basis()`](https://mcomas.net/coda.base/reference/pb_basis.md)
  with the tabu-search method and improved the C++ implementation of
  principal-balance searches.
- Expanded validation and tests for compositional operations, distances,
  zero replacement and tabu-search neighbourhoods.
- Corrected the `statistitian_time` data set name to
  `statistician_time`.

## coda.base 1.0.5 (Release date: 2026-03-03)

CRAN release: 2026-03-04

- Reworked the basis and coordinate code, with stronger checks for
  dimensions, names and invalid compositions.
- Added
  [`cdp_basis()`](https://mcomas.net/coda.base/reference/cdp_basis.md)
  and
  [`conditional_obasis()`](https://mcomas.net/coda.base/reference/conditional_obasis.md)
  to the public API.
- Added
  [`gen_coda_with_zeros_and_missings()`](https://mcomas.net/coda.base/reference/gen_coda_with_zeros_and_missings.md)
  for generating data sets used to assess zero and missing-value
  replacement methods.
- Improved
  [`coda_replacement()`](https://mcomas.net/coda.base/reference/coda_replacement.md)
  and the conditional orthonormal-basis calculations used for zero and
  missing-value imputation.
- Improved
  [`sbp_basis()`](https://mcomas.net/coda.base/reference/sbp_basis.md),
  constrained principal balances and the C++ principal-balance
  implementation.
- [`plot_balance()`](https://mcomas.net/coda.base/reference/plot_balance.md)
  can display statistics associated with balances.
- Optimised and reorganised the C++ code and expanded the basis tests.

## coda.base 1.0.4 (Release date: 2025-09-20)

- Changed the sign convention of CoDaPack bases to agree with CoDaPack.

## coda.base 1.0.3 (Release date: 2025-07-02)

CRAN release: 2025-07-02

- Improved completion of sequential binary partitions in
  [`sbp_basis()`](https://mcomas.net/coda.base/reference/sbp_basis.md).
- Improved constrained principal-balance calculations and their handling
  of subcompositions.

## coda.base 1.0.2 (Release date: 2025-05-07)

- Added the `kilauea_iki` compositional data set and its documentation.
- Improved documentation for zero replacement and coordinate
  transformations.
- Improved dimension checks for pairwise bases and composition
  reconstruction.

## coda.base 1.0.0 (Release date: 2025-05-02)

CRAN release: 2025-05-02

- Added
  [`coda_replacement()`](https://mcomas.net/coda.base/reference/coda_replacement.md)
  for imputing rounded zeros and missing values in compositional data.
- Added conditional orthonormal-basis calculations used by the
  imputation procedure.
- Reorganised basis, coordinate and constrained principal-balance code,
  with more consistent argument validation and naming.
- Updated the package documentation and test suite for the 1.0 release.

## coda.base 0.5.5 (Release date: 2023-11-25)

CRAN release: 2023-11-25

- Added direct coordinate helpers
  [`alr_c()`](https://mcomas.net/coda.base/reference/coordinates.md),
  [`clr_c()`](https://mcomas.net/coda.base/reference/coordinates.md),
  [`ilr_c()`](https://mcomas.net/coda.base/reference/coordinates.md) and
  [`olr_c()`](https://mcomas.net/coda.base/reference/coordinates.md).
- Added
  [`olr_basis()`](https://mcomas.net/coda.base/reference/ilr_basis.md)
  and support for orthonormal log-ratio coordinates in
  [`coordinates()`](https://mcomas.net/coda.base/reference/coordinates.md)
  and
  [`composition()`](https://mcomas.net/coda.base/reference/composition.md).
- Extended
  [`variation_array()`](https://mcomas.net/coda.base/reference/variation_array.md)
  with optional log-ratio means and changed its defaults to return the
  variation matrix directly.
- Improved printing of coordinate objects and preservation of basis
  metadata.

## coda.base 0.5.4 (Release date: 2023-03-01)

- Added a collection of documented compositional data sets, including
  `alimentation`, `blood_mn`, `bmi_activity`, `eurostat_employment`,
  `foraminiferals`, `house_expend`, `mammals_milk`, `milk_cows`,
  `montana`, `petrafm`, `pollen`, `serprot`, `statistitian_time` and
  `weibo_hotels`.
- [`plot_balance()`](https://mcomas.net/coda.base/reference/plot_balance.md)
  now uses balance names in dendrogram labels.
- Fixed completion of incomplete sequential binary partitions in
  [`sbp_basis()`](https://mcomas.net/coda.base/reference/sbp_basis.md).
- Removed
  [`cc_basis()`](https://mcomas.net/coda.base/reference/cc_basis.md)
  because its results did not match the intended method; it was
  redesigned and restored in a later release.

## coda.base 0.5.3 (Release date: 2023-02-08)

- [`sbp_basis()`](https://mcomas.net/coda.base/reference/sbp_basis.md)
  can complete an incomplete sequential binary partition to obtain a
  full basis.

## coda.base 0.5.2 (Release date: 2022-07-14)

CRAN release: 2022-07-18

- Basis constructors consistently return standard R matrices; sparse
  matrices are used internally where they improve performance.
- Updated the minimum `RcppArmadillo` requirement for the compiled code.

## coda.base 0.5.1 (Release date: 2022-07-07)

- Added the Wilks statistic as a criterion for canonical-correlation
  balances.

## coda.base 0.4.3

- center() and gmean() functions incorporated.

## coda.base 0.4.2

- Pairwise default name changed to lr instead of alr.
- Sparse matrices are used for logcontrast matrices of balances.
- Problem with constrained principal balances detected and solved.

## coda.base 0.3.2 (Release date: 2022-03-13)

- [`coord()`](https://mcomas.net/coda.base/reference/coordinates.md) and
  [`comp()`](https://mcomas.net/coda.base/reference/composition.md) are
  included as a shorter names of
  [`coordinates()`](https://mcomas.net/coda.base/reference/coordinates.md)
  and
  [`composition()`](https://mcomas.net/coda.base/reference/composition.md).
- Basis names taken from composition.
- Coordinates based on pairwise log-ratios available with
  [`pairwise_basis()`](https://mcomas.net/coda.base/reference/pairwise_basis.md)
  or basis = ‘pw’ in function coordinates().
- Parameter label removed in function coordinates().

## coda.base 0.3.1

CRAN release: 2020-05-14

- PB constrained approximation available. Adaptation of code provided by
  J.A.Martín-Fernández.
- PB exact algorithm rewritten to improve performance.
- PB local search algorithm removed.

## coda.base 0.2.2

CRAN release: 2020-04-16

- MASS dependence removed
- Code reorganised to avoid repetition.

## coda.base 0.2.1

CRAN release: 2019-07-02

- ?\_basis functions always return a matrix
- consistencies between methods improved

## coda.base 0.1.12 (Release date: 2018-04-15)

CRAN release: 2019-04-15

- coordinates function keeps row.names attribute
- Principal components calculations are done using svd function which is
  faster than princomp
- sbp_basis function does not need package to be loaded, can run as
  coda.base::sbp_basis()

## coda.base 0.1.11 (Release date: 2018-12-23)

CRAN release: 2019-01-15

- Principal components coordinates uncentered.
- Sequential binary partition can be defined using a matrix with -1,0
  and 1’s.
- CoDaPack’s basis included.
- A vignette explaining how to build coordinates has been included.

## coda.base 0.1.10 (Release date: 2018-08-03)

CRAN release: 2018-08-06

- Principal Balances included as direct method in coordinates function.
- Coordinates function filters non-valid compositions.

## coda.base 0.1.9 (Release date: 2018-05-03)

CRAN release: 2018-05-09

- Modifications included to rcran policies.

## coda.base 0.1.8 (Release date: 2017-12-03)

CRAN release: 2017-12-03

- Exact algorithm for Principal Balances included.
- Bug in lsearch base in principal components solved.

## coda.base 0.1.7 (Release date: 2017-08-25)

CRAN release: 2017-09-26

- Problem with function composition using ‘clr’ basis solved.
- ‘coda.base.basis’ option included to hide the basis when printing the
  coordinates.
- New method which uses the principal components included in pb_basis.
- Documentaion improved.
