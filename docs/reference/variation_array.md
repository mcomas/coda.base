# Variation array is returned.

Variation array is returned.

## Usage

``` r
variation_array(X, include_means = FALSE, ml_covariance = FALSE)
```

## Arguments

- X:

  Compositional dataset

- include_means:

  if TRUE logratio means are included in the lower-left triangle

- ml_covariance:

  if TRUE Maximum-likelihood estimation of the covariance for the
  multivariate normal distribution is used (dividing the scatter matrix
  by n instead of n-1)

## Value

variation array matrix

## Examples

``` r
set.seed(1)
X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
variation_array(X)
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 0.000000 1.726005 1.842314 1.880032 1.865787
#> [2,] 1.726005 0.000000 2.085389 2.011943 1.982605
#> [3,] 1.842314 2.085389 0.000000 1.825880 2.365846
#> [4,] 1.880032 2.011943 1.825880 0.000000 2.295994
#> [5,] 1.865787 1.982605 2.365846 2.295994 0.000000
variation_array(X, include_means = TRUE)
#>             [,1]         [,2]        [,3]       [,4]     [,5]
#> [1,]  0.00000000  1.726005360  1.84231397  1.8800325 1.865787
#> [2,] -0.14669544  0.000000000  2.08538862  2.0119435 1.982605
#> [3,] -0.07921383  0.067481612  0.00000000  1.8258798 2.365846
#> [4,] -0.05728551  0.089409935  0.02192832  0.0000000 2.295994
#> [5,] -0.14802161 -0.001326164 -0.06880778 -0.0907361 0.000000
```
