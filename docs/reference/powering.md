# Powering of compositional data

The powering operation raises each part of a composition to a scalar
exponent and then applies closure to re-normalize the result as a
composition.

## Usage

``` r
powering(X, alpha)
```

## Arguments

- X:

  A numeric vector, matrix or data.frame containing compositions.

- alpha:

  A numeric scalar or vector. If `X` is a matrix or data.frame, `alpha`
  may have length 1, length equal to the number of rows of `X`, or
  length equal to the number of parts of `X`. If `X` is a vector and
  `alpha` has length greater than 1, one powered composition is returned
  for each value of `alpha`.

## Value

An object with the same format as `X` containing the powered
compositions, except that vector `X` with vector `alpha` returns a
matrix.

## Details

Powering is the analogue of scalar multiplication in the simplex. Each
part is raised to `alpha`, and the result is closed with
[`closure`](https://mcomas.net/coda.base/reference/closure.md). When
`alpha` has one value per row, each composition is powered by its
corresponding value. When it has one value per part, each part receives
its corresponding exponent. For vector `X` and vector `alpha`, each row
of the result is `X` powered by the corresponding element of `alpha`.

## Examples

``` r
x <- c(a = 1, b = 2, c = 3)
powering(x, 2)
#> [1] 0.07142857 0.28571429 0.64285714
powering(x, c(1, 2))
#>               a         b         c
#> [1,] 0.16666667 0.3333333 0.5000000
#> [2,] 0.07142857 0.2857143 0.6428571

X <- rbind(
  c(1, 2, 3),
  c(4, 5, 6)
)
powering(X, c(1, 2))
#>           [,1]      [,2]      [,3]
#> [1,] 0.1666667 0.3333333 0.5000000
#> [2,] 0.2077922 0.3246753 0.4675325
```
