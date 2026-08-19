# Generate a random composition with a prescribed first principal balance

Simulates a random composition whose coordinate system is constructed
from a sequential binary partition induced by a given first balance. The
supplied balance is completed to a full orthonormal basis using
[`sbp_basis`](https://mcomas.net/coda.base/reference/sbp_basis.md) with
`fill = TRUE`.

## Usage

``` r
random_composition_with_fixed_pb(principal_balance, n = 100, sd1 = 5)
```

## Arguments

- principal_balance:

  An integer or numeric vector in \\\\-1,0,1\\\\ defining the first
  balance of the SBP.

- n:

  Integer. Number of observations to generate.

- sd1:

  Numeric value used to scale the first latent coordinate before
  rotating the simulated coordinates.

## Value

A composition matrix with `n` rows and `length(principal_balance)`
columns.

## Details

Standard normal latent coordinates are first generated in dimension
\\D - 1\\, where \\D\\ is the number of parts. Their sample covariance
matrix is then diagonalized, and the associated eigenvectors are used to
rotate the latent coordinates before mapping them back to the simplex
using the basis induced by `principal_balance`.

This function is mainly intended for examples, simulation studies, and
experiments where a specific first balance structure is desired.

## See also

[`sbp_basis`](https://mcomas.net/coda.base/reference/sbp_basis.md),
[`composition`](https://mcomas.net/coda.base/reference/composition.md)
