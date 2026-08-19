# Replacement of missing values and below-detection zeros in compositional data

Performs imputation of missing values and/or values below the detection
limit in compositional data using an EM algorithm assuming normality on
the simplex.

## Usage

``` r
coda_replacement(
  X,
  DL = NULL,
  dl_prop = 0.65,
  eps = 1e-04,
  parameters = FALSE,
  debug = FALSE,
  maxit = 500
)
```

## Arguments

- X:

  A compositional data set: numeric matrix or data frame where rows
  represent observations and columns represent parts.

- DL:

  An optional matrix or vector of detection limits. If \`NULL\`, the
  minimum non-zero value in each column of \`X\` is used.

- dl_prop:

  A numeric value between 0 and 1 used for initialization in the EM
  algorithm.

- eps:

  Convergence tolerance.

- parameters:

  Logical; if \`TRUE\`, return additional estimated parameters.

- debug:

  Logical; if \`TRUE\`, print the log-likelihood at each iteration.

- maxit:

  Maximum number of iterations

## Value

If \`parameters = FALSE\`, the imputed object with the same format as
\`X\` (\`matrix\` or \`data.frame\`, preserving data-frame subclasses
when possible) and preserving original names. If \`parameters = TRUE\`,
a list with the estimated clr mean, clr covariance, and imputed clr
coordinates.

## Examples

``` r
X <- matrix(c(
  0.00, 0.30, 0.70,
  0.20,   NA, 0.80,
  0.40, 0.60, 0.00,
  0.25, 0.25, 0.50,
  0.10, 0.30, 0.60
), ncol = 3, byrow = TRUE)
colnames(X) <- c("sand", "silt", "clay")

DL <- c(0.05, 0.05, 0.05)

X_imp <- coda_replacement(X, DL = DL, maxit = 20)
X_imp
#>           sand      silt      clay
#> [1,] 0.0314770 0.2905569 0.6779661
#> [2,] 0.1246784 0.3766081 0.4987135
#> [3,] 0.3874092 0.5811138 0.0314770
#> [4,] 0.2500000 0.2500000 0.5000000
#> [5,] 0.1000000 0.3000000 0.6000000

set.seed(10)
X <- composition(matrix(rnorm(3*10), ncol = 3))
X[sample(c(TRUE, FALSE), 4*10, replace = TRUE, c(1, 3))] <- NA
params <- coda_replacement(X, parameters = TRUE, debug = TRUE)
#> Iteration: 1, LogLik: -21.509
#> Iteration: 2, LogLik: -21.4807
#> Iteration: 3, LogLik: -21.4702
#> Iteration: 4, LogLik: -21.4663
#> Iteration: 5, LogLik: -21.4649
names(params)
#> [1] "clr_mu"    "clr_sigma" "clr_h"    
params$clr_mu
#>            [,1]
#> [1,] -0.5855132
#> [2,]  0.3387310
#> [3,] -0.6110079
#> [4,]  0.8577901
composition(params$clr_h)
#>               c1        c2         c3        c4         c5
#>  [1,] 0.19717945 0.1818282 0.18233840 0.2869558 0.15169813
#>  [2,] 0.08937430 0.1695037 0.15571250 0.5624073 0.02300224
#>  [3,] 0.09229793 0.5497201 0.09823144 0.1708750 0.08887548
#>  [4,] 0.08137958 0.1989978 0.09778718 0.5992803 0.02255514
#>  [5,] 0.11811059 0.2404526 0.18020017 0.4118082 0.04942845
#>  [6,] 0.14389646 0.2562383 0.26652872 0.2408355 0.09250100
#>  [7,] 0.09847978 0.4169787 0.08725439 0.3375058 0.05978130
#>  [8,] 0.12951236 0.3877787 0.12353546 0.2536799 0.10549366
#>  [9,] 0.12750912 0.3962536 0.03585914 0.3018078 0.13857035
#> [10,] 0.19223702 0.2085054 0.13775352 0.3032841 0.15821987
```
