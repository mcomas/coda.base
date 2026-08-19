# Dataset center

Generic function to calculate the center of a compositional dataset

## Usage

``` r
center(X, zero.rm = FALSE, na.rm = FALSE)
```

## Arguments

- X:

  compositional dataset

- zero.rm:

  a logical value indicating whether zero values should be stripped
  before the computation proceeds.

- na.rm:

  a logical value indicating whether NA values should be stripped before
  the computation proceeds.

## Examples

``` r
X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
g = rep(c('a','b','c','d'), 25)
center(X)
#> [1] 0.2059718 0.2152515 0.1862152 0.2317245 0.1608370
(by_g <- by(X, g, center))
#> INDICES: a
#>        V1        V2        V3        V4        V5 
#> 0.2318552 0.1643589 0.1621590 0.2570751 0.1845518 
#> ----------------------------------------------------- 
#> INDICES: b
#>        V1        V2        V3        V4        V5 
#> 0.2201313 0.2056451 0.2210786 0.2197650 0.1333801 
#> ----------------------------------------------------- 
#> INDICES: c
#>        V1        V2        V3        V4        V5 
#> 0.1683561 0.2157972 0.2166346 0.1995019 0.1997101 
#> ----------------------------------------------------- 
#> INDICES: d
#>        V1        V2        V3        V4        V5 
#> 0.1993826 0.2801630 0.1473765 0.2435038 0.1295741 
center(t(simplify2array(by_g)))
#>        V1        V2        V3        V4        V5 
#> 0.2059718 0.2152515 0.1862152 0.2317245 0.1608370 
```
