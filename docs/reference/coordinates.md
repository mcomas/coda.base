# Coordinates of compositions with respect to a basis

Compute coordinates of a composition or a compositional data set with
respect to a given log-ratio basis.

The \`basis\` argument can be either:

- a character string identifying a predefined coordinate system, or

- a matrix whose columns define a system of log-contrasts.

The predefined options are:

- \`"ilr"\`: isometric log-ratio coordinates,

- \`"olr"\`: orthonormal log-ratio coordinates,

- \`"clr"\`: centered log-ratio coordinates,

- \`"alr"\`: additive log-ratio coordinates,

- \`"pw"\`: pairwise log-ratios,

- \`"pc"\`: principal component log-ratio coordinates,

- \`"pb"\`: principal balance coordinates,

- \`"cdp"\`: CoDaPack default balances.

## Usage

``` r
coordinates(X, basis = "ilr")

coord(..., basis = "ilr")

alr_c(X)

clr_c(X)

ilr_c(X)

olr_c(X)
```

## Arguments

- X:

  A compositional data set. It can be a numeric matrix, a data frame, or
  a numeric vector.

- basis:

  Basis used to compute the coordinates. Either a character string
  naming a predefined basis or a matrix with log-ratio basis vectors in
  columns.

- ...:

  components of the composition

## Value

Coordinates of \`X\` with respect to the given \`basis\`. The returned
object has the same general type as the input when possible.

## See also

[`ilr_basis`](https://mcomas.net/coda.base/reference/ilr_basis.md),
[`alr_basis`](https://mcomas.net/coda.base/reference/alr_basis.md),
[`clr_basis`](https://mcomas.net/coda.base/reference/clr_basis.md),
[`sbp_basis`](https://mcomas.net/coda.base/reference/sbp_basis.md),
[`composition`](https://mcomas.net/coda.base/reference/composition.md)

## Examples

``` r
coordinates(1:5)
#>       ilr1       ilr2       ilr3       ilr4 
#> -0.4901291 -0.6140370 -0.6833297 -0.7288906 

B <- ilr_basis(5)
coordinates(1:5, B)
#>       ilr1       ilr2       ilr3       ilr4 
#> -0.4901291 -0.6140370 -0.6833297 -0.7288906 

X <- rbind(1:5, 2:6)
coordinates(X, "clr")
#>            clr1       clr2       clr3      clr4      clr5
#> [1,] -0.9574983 -0.2643512 0.14111394 0.4287960 0.6519396
#> [2,] -0.6227031 -0.2172380 0.07044412 0.2935877 0.4759092
```
