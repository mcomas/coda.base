# \`coda.base\` features' enumeration

``` r
library(coda.base)
#> 
#> Attaching package: 'coda.base'
#> The following object is masked from 'package:stats':
#> 
#>     dist
```

## Main functions

- Simulated data set:

``` r
nP = 6
nO = 8
set.seed(1)
X = matrix(rlnorm(nP * nO), ncol = nP, nrow = nO)
colnames(X) = paste0('P', 1:nP)
rownames(X) = paste0('O', 1:nO)
X
#>           P1        P2        P3        P4        P5        P6
#> O1 0.5344838 1.7785196 0.9839401 1.8586041 1.4735458 0.8482977
#> O2 1.2015872 0.7368371 2.5698209 0.9454174 0.9476168 0.7761871
#> O3 0.4336018 4.5348008 2.2732743 0.8557342 0.2523194 2.0076470
#> O4 4.9297132 1.4767493 1.8110401 0.2297526 0.6603439 1.7448406
#> O5 1.3902836 0.5372775 2.5067256 0.6199292 0.6741586 0.5022006
#> O6 0.4402254 0.1091863 2.1861375 1.5188319 0.9424114 0.4928772
#> O7 1.6281250 3.0800041 1.0774154 3.8910520 3.0042422 1.4399119
#> O8 2.0924271 0.9560610 0.1367841 0.9023185 2.1450776 2.1566000
```

##### Descriptive statistics

- Center:

``` r
center(X)
#>        P1        P2        P3        P4        P5        P6 
#> 0.1732294 0.1599997 0.1982046 0.1544340 0.1509712 0.1631611
```

- Variation array:

``` r
variation_array(X)
#>           P1        P2       P3        P4        P5        P6
#> P1 0.0000000 1.7807354 2.154229 2.1457756 0.9638581 0.6955516
#> P2 1.7807354 0.0000000 2.614290 1.9008474 2.0593876 0.6793622
#> P3 2.1542290 2.6142898 0.000000 1.9166209 2.4966568 1.9518749
#> P4 2.1457756 1.9008474 1.916621 0.0000000 0.5296424 1.2179905
#> P5 0.9638581 2.0593876 2.496657 0.5296424 0.0000000 0.9414218
#> P6 0.6955516 0.6793622 1.951875 1.2179905 0.9414218 0.0000000
```

##### Aitchison distances between compositions

[`stats::dist()`](https://rdrr.io/r/stats/dist.html) is rewritten to
include the Aitchison distance between compositions:

``` r
dist(X, method = 'aitchison')
#> Warning: Use dist_coda(x, method = "aitchison") instead. The Aitchison
#> extension in coda.base::dist() will be removed in a future version.
#>           O1        O2        O3        O4        O5        O6
#> O2 1.7312698                                                  
#> O3 2.4651454 2.6482563                                        
#> O4 3.2924718 2.2982518 3.1363020                              
#> O5 2.1967501 0.5310050 2.9124924 2.1147115                    
#> O6 2.6637001 1.8840274 4.0189064 3.9336894 2.0813815          
#> O7 0.7482808 1.9679763 2.9204289 3.2563331 2.3587517 2.8831297
#> O8 2.7579161 3.2673090 4.1643719 3.2700093 3.5174190 4.0863918
#>           O7
#> O2          
#> O3          
#> O4          
#> O5          
#> O6          
#> O7          
#> O8 2.2159481
```

##### Log-ratio coordinates. `coordinates()`

``` r
coordinates(X)
#>          ilr1        ilr2       ilr3       ilr4       ilr5
#> O1 -0.8501086 -0.00746765 -0.5560864 -0.2230977  0.3219244
#> O2  0.3457976 -0.82034125  0.2859262  0.2193991  0.3613081
#> O3 -1.6598694 -0.39448617  0.5671774  1.5316655 -0.6427156
#> O4  0.8523731  0.32550510  2.0182094  0.6189998 -0.3815869
#> O5  0.6722806 -0.86944130  0.5951600  0.3860022  0.5839787
#> O6  0.9858706 -1.87771387 -1.0123423 -0.3572872  0.2999825
#> O7 -0.4507819  0.59736115 -0.6896777 -0.3028747  0.4240687
#> O8  0.5538473  1.90737459 -0.2850948 -0.9953748 -0.8176105
```

By default, `coda.base` uses the isometric log-ratio coordinates defined
in Egozcue et al. 2003.

- `coordinates(X, 'ilr')`: isometric log-ratio coordinates (Egozcue et
  al. 2003, defaults)
- `coordinates(X, 'olr')`: orthonormal log-ratio coordinates (equivalent
  to `ilr`)
- `coordinates(X, 'alr')`: additive log-ratio coordinates
- `coordinates(X, 'clr')`: centered log-ratio coordinates
- `coordinates(X, 'pw')`: pairwise log-ratio coordinates
- `coordinates(X, 'pc')`: principal component coordinates
- `coordinates(X, 'pb')`: principal balance coordinates
- `coordinates(X, 'cdp')`: balanced isometric log-ratio coordinates

To reduce typing,
[`alr_c()`](https://mcomas.net/coda.base/reference/coordinates.md),
[`clr_c()`](https://mcomas.net/coda.base/reference/coordinates.md),
[`ilr_c()`](https://mcomas.net/coda.base/reference/coordinates.md) and
[`olr_c()`](https://mcomas.net/coda.base/reference/coordinates.md) are
functions that call
[`coordinates()`](https://mcomas.net/coda.base/reference/coordinates.md)
function with the option given by their name.

`coordinates(X, B)` accepts a log-contrast matrix $`B`$ to build the
log-ratio coordinates. Different log-contrast matrices $`B`$ can be
constructed (following section).

##### Functions to build log-contrast matrices

``` r
ilr_basis(nP)
#>          ilr1       ilr2       ilr3       ilr4       ilr5
#> c1  0.7071068  0.4082483  0.2886751  0.2236068  0.1825742
#> c2 -0.7071068  0.4082483  0.2886751  0.2236068  0.1825742
#> c3  0.0000000 -0.8164966  0.2886751  0.2236068  0.1825742
#> c4  0.0000000  0.0000000 -0.8660254  0.2236068  0.1825742
#> c5  0.0000000  0.0000000  0.0000000 -0.8944272  0.1825742
#> c6  0.0000000  0.0000000  0.0000000  0.0000000 -0.9128709
```

``` r
all.equal(as.numeric(coordinates(X, 'ilr')),
          as.numeric(log(X) %*% ilr_basis(nP)))
#> [1] TRUE
```

Log-ratio matrix transformations:

- `ilr_basis(nP)` or `olr_basis(nP)` (Egozcue et al. 2003, defaults)
- `ilr_basis(nP, type = 'pivot')` or `olr_basis(nP, type = 'pivot')` to
  pivot log-ratio coordinates
- `ilr_basis(nP, type = 'cdp')` or `olr_basis(nP, type = 'cdp')` to
  balanced log-ratio coordinates (CoDaPack’s default)
- `alr_basis(nP)` to additive-log ratio coordinates. Numerator order and
  denominator can be modified. For example,
  `alr_basis(nP, denominator = 1, numerator = nP:2)`.
- `clr_basis(nP)` to centered log-ratio coordinates.
- `pc_basis(X)` to principal components log-ratio coordinates.
- `cc_basis(X, X2)` to canonical correlations log-ratio coordinates.
- `pb_basis(X, method = 'exact')` to principal balances using the exact
  algorithm.
- `pb_basis(X, method = 'constrained')` to principal balances using pca
  constrained algorithm.
- `pb_basis(X, method = 'cluster')` to principal balances obtained using
  parts clustering algorithm.
- `pairwise_basis(nP)` to pairwise log-ratio coordinates.
- `sbp_basis(b0 = b1~b2, b1 = P1~P2+P3, b2 = P4~P5+P6, data=X)`
