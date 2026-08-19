# Additive log-ratio basis

Construct the transformation matrix associated with additive log-ratio
(alr) coordinates.

## Usage

``` r
alr_basis(dim, denominator = NULL, numerator = NULL)
```

## Arguments

- dim:

  Number of parts. It can be a single integer, a matrix or data frame,
  or a character vector of part names.

- denominator:

  Part used as denominator. By default, the last part is used.

- numerator:

  Parts used as numerators. By default, all parts except the denominator
  are used, preserving their original order.

## Value

A matrix defining the alr coordinate system.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman & Hall, London.

## Examples

``` r
alr_basis(5)
#>    alr1 alr2 alr3 alr4
#> c1    1    0    0    0
#> c2    0    1    0    0
#> c3    0    0    1    0
#> c4    0    0    0    1
#> c5   -1   -1   -1   -1
alr_basis(5, 3)
#>    alr1 alr2 alr3 alr4
#> c1    1    0    0    0
#> c2    0    1    0    0
#> c3   -1   -1   -1   -1
#> c4    0    0    1    0
#> c5    0    0    0    1
alr_basis(5, 3, c(1, 5, 2, 4))
#>    alr1 alr2 alr3 alr4
#> c1    1    0    0    0
#> c2    0    0    1    0
#> c3   -1   -1   -1   -1
#> c4    0    0    0    1
#> c5    0    1    0    0
```
