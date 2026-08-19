# Compositions from coordinates with respect to a basis

Reconstruct a composition from coordinates with respect to a given
basis.

## Usage

``` r
composition(H, basis = "ilr")

comp(H, basis = "ilr")
```

## Arguments

- H:

  Coordinates of a composition. It can be a numeric matrix, a data
  frame, or a numeric vector.

- basis:

  Basis used to interpret the coordinates. Either a character string
  naming a predefined basis or a matrix.

## Value

A composition corresponding to the given coordinates.

## See also

[`coordinates`](https://mcomas.net/coda.base/reference/coordinates.md),
[`ilr_basis`](https://mcomas.net/coda.base/reference/ilr_basis.md),
[`alr_basis`](https://mcomas.net/coda.base/reference/alr_basis.md),
[`clr_basis`](https://mcomas.net/coda.base/reference/clr_basis.md),
[`sbp_basis`](https://mcomas.net/coda.base/reference/sbp_basis.md)
