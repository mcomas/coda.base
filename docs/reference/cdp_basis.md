# CoDaPack default ilr basis

Construct the default isometric log-ratio basis used in CoDaPack.

## Usage

``` r
cdp_basis(dim)
```

## Arguments

- dim:

  Number of parts. It can be a single integer, a matrix or data frame,
  or a character vector of part names.

## Value

A matrix with \\D\\ rows and \\D - 1\\ columns containing the CoDaPack
default ilr basis.

## Examples

``` r
cdp_basis(5)
#>          ilr1       ilr2       ilr3       ilr4
#> c1  0.3651484  0.4082483  0.7071068  0.0000000
#> c2  0.3651484  0.4082483 -0.7071068  0.0000000
#> c3  0.3651484 -0.8164966  0.0000000  0.0000000
#> c4 -0.5477226  0.0000000  0.0000000  0.7071068
#> c5 -0.5477226  0.0000000  0.0000000 -0.7071068
cdp_basis(c("a", "b", "c", "d"))
#>   ilr1       ilr2       ilr3
#> a  0.5  0.7071068  0.0000000
#> b  0.5 -0.7071068  0.0000000
#> c -0.5  0.0000000  0.7071068
#> d -0.5  0.0000000 -0.7071068
```
