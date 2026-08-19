# Centered log-ratio basis

Construct the transformation matrix associated with centered log-ratio
(clr) coordinates.

## Usage

``` r
clr_basis(dim)
```

## Arguments

- dim:

  Number of parts. It can be a single integer, a matrix or data frame,
  or a character vector of part names.

## Value

A square matrix defining the clr coordinate system.

## Details

CLR coordinates are linearly dependent and lie in the \\D - 1\\
dimensional clr-plane.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman & Hall, London.

## Examples

``` r
B <- clr_basis(5)
clr_coordinates <- coordinates(c(1, 2, 3, 4, 5), B)
sum(clr_coordinates) < 1e-15
#> [1] TRUE
```
