# Generate compositional data with zeros and missing values

Simulate compositional data and optionally introduce structural zeros
(interpreted as values below a detection limit) and missing values.

The function first generates a compositional data set \`X0\`, then
creates a modified version \`X\` by:

- replacing values below \`dl_par\` by zero, if \`zeros = TRUE\`,

- introducing missing values at random, if \`missings = TRUE\`.

A matrix of detection limits \`DL\` is also returned. It contains
\`dl_par\` in the positions that were censored to zero, and \`0\`
elsewhere.

## Usage

``` r
gen_coda_with_zeros_and_missings(
  n,
  d,
  missings = TRUE,
  zeros = TRUE,
  dl_par = 0.05,
  na_p = 0.15
)
```

## Arguments

- n:

  Number of observations.

- d:

  Dimension of the latent coordinate space used to generate the
  compositions.

- missings:

  Logical; if \`TRUE\`, introduce missing values at random.

- zeros:

  Logical; if \`TRUE\`, replace values below \`dl_par\` by zero.

- dl_par:

  Detection-limit threshold used to generate zeros.

- na_p:

  Probability that any entry is replaced by \`NA\` when \`missings =
  TRUE\`.

## Value

A list with three components:

- X:

  The generated compositional data set with simulated zeros and/or
  missing values.

- DL:

  A matrix of detection limits, with \`dl_par\` in censored positions
  and \`0\` elsewhere.

- X0:

  The original simulated compositional data set before introducing zeros
  or missing values.

## Details

Compositions are generated from multivariate normal coordinates and
mapped to the simplex through \`composition()\`. The eigenvector
rotation is included to induce a non-trivial covariance structure in the
generated coordinates.

Missing values are introduced completely at random, independently for
each cell, with probability \`na_p\`.

## Examples

``` r
set.seed(123)
sim <- gen_coda_with_zeros_and_missings(100, 4)

str(sim)
#> List of 3
#>  $ X : num [1:100, 1:5] 0 0.0591 0.1205 0.0656 0.1226 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:5] "c1" "c2" "c3" "c4" ...
#>  $ DL: num [1:100, 1:5] 0.05 0 0 0 0 0 0.05 0 0 0 ...
#>  $ X0: num [1:100, 1:5] 0.0233 0.0591 0.1205 0.0656 0.1226 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:5] "c1" "c2" "c3" "c4" ...
summary(sim$X0)
#>        c1                c2                 c3         
#>  Min.   :0.01082   Min.   :0.007262   Min.   :0.01432  
#>  1st Qu.:0.07792   1st Qu.:0.088177   1st Qu.:0.08084  
#>  Median :0.15044   Median :0.172085   Median :0.15281  
#>  Mean   :0.19450   Mean   :0.220267   Mean   :0.19687  
#>  3rd Qu.:0.27525   3rd Qu.:0.330935   3rd Qu.:0.27919  
#>  Max.   :0.76669   Max.   :0.668187   Max.   :0.83365  
#>        c4                c5          
#>  Min.   :0.01230   Min.   :0.009548  
#>  1st Qu.:0.09951   1st Qu.:0.073729  
#>  Median :0.17773   Median :0.131945  
#>  Mean   :0.21763   Mean   :0.170740  
#>  3rd Qu.:0.30501   3rd Qu.:0.233579  
#>  Max.   :0.78509   Max.   :0.724166  
summary(sim$X)
#>        c1                c2                c3         
#>  Min.   :0.00000   Min.   :0.00000   Min.   :0.00000  
#>  1st Qu.:0.06834   1st Qu.:0.09224   1st Qu.:0.07759  
#>  Median :0.14290   Median :0.18188   Median :0.14802  
#>  Mean   :0.18021   Mean   :0.22393   Mean   :0.19039  
#>  3rd Qu.:0.26717   3rd Qu.:0.33782   3rd Qu.:0.27905  
#>  Max.   :0.60122   Max.   :0.66819   Max.   :0.83365  
#>  NAs    :12        NAs    :17        NAs    :7        
#>        c4               c5         
#>  Min.   :0.0000   Min.   :0.00000  
#>  1st Qu.:0.1049   1st Qu.:0.07564  
#>  Median :0.1777   Median :0.13214  
#>  Mean   :0.2170   Mean   :0.17136  
#>  3rd Qu.:0.3002   3rd Qu.:0.23383  
#>  Max.   :0.7851   Max.   :0.72417  
#>  NAs    :14       NAs    :15       
table(sim$X == 0, useNA = "ifany")
#> 
#> FALSE  TRUE  <NA> 
#>   379    56    65 
```
