# CoDaPack's default binary partition

Compute the default binary partition used in CoDaPack's software

## Usage

``` r
cdp_partition(ncomp)
```

## Arguments

- ncomp:

  number of parts

## Value

matrix

## Examples

``` r
cdp_partition(4)
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
#> [2,]    1   -1    0
#> [3,]   -1    0    1
#> [4,]   -1    0   -1
```
