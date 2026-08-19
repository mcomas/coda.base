# Basis from a sequential binary partition

Construct a balance basis from a sequential binary partition (SBP) or
from a more general collection of balances.

## Usage

``` r
sbp_basis(sbp, data = NULL, fill = FALSE, silent = FALSE)
```

## Arguments

- sbp:

  A list of formulas or a matrix describing balances.

- data:

  Optional compositional data set used to extract part names when
  \`sbp\` is given as a list of formulas. If \`data = NULL\`, part names
  are inferred from the formulas themselves.

- fill:

  Logical; if \`TRUE\`, complete the supplied balances to obtain a full
  basis.

- silent:

  Logical; if \`FALSE\`, report whether the resulting balances form a
  basis, and whether they are orthogonal or orthonormal.

## Value

A matrix whose columns are balances.

## Details

The argument \`sbp\` can be specified in two ways:

- as a list of formulas, where each formula defines the numerator and
  the denominator groups of a balance,

- as a matrix with one column per balance and one row per part. Positive
  entries indicate parts in the numerator, negative entries indicate
  parts in the denominator, and zeros indicate unused parts.

## Examples

``` r
X <- data.frame(
  a = 1:2, b = 2:3, c = 4:5,
  d = 5:6, e = 10:11, f = 100:101, g = 1:2
)

# Sequential SBP construction
sbp_basis(list(
  b1 = a ~ b + c + d + e + f + g,
  b2 = b ~ c + d + e + f + g,
  b3 = c ~ d + e + f + g,
  b4 = d ~ e + f + g,
  b5 = e ~ f + g,
  b6 = f ~ g
), data = X)
#>           b1         b2         b3         b4         b5         b6
#> a  0.9258201  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> b -0.1543033  0.9128709  0.0000000  0.0000000  0.0000000  0.0000000
#> c -0.1543033 -0.1825742  0.8944272  0.0000000  0.0000000  0.0000000
#> d -0.1543033 -0.1825742 -0.2236068  0.8660254  0.0000000  0.0000000
#> e -0.1543033 -0.1825742 -0.2236068 -0.2886751  0.8164966  0.0000000
#> f -0.1543033 -0.1825742 -0.2236068 -0.2886751 -0.4082483  0.7071068
#> g -0.1543033 -0.1825742 -0.2236068 -0.2886751 -0.4082483 -0.7071068

# Chain construction
sbp_basis(list(
  b1 = a ~ b,
  b2 = b1 ~ c,
  b3 = b2 ~ d,
  b4 = b3 ~ e,
  b5 = b4 ~ f,
  b6 = b5 ~ g
), data = X)
#>           b1         b2         b3         b4         b5         b6
#> a  0.7071068  0.4082483  0.2886751  0.2236068  0.1825742  0.1543033
#> b -0.7071068  0.4082483  0.2886751  0.2236068  0.1825742  0.1543033
#> c  0.0000000 -0.8164966  0.2886751  0.2236068  0.1825742  0.1543033
#> d  0.0000000  0.0000000 -0.8660254  0.2236068  0.1825742  0.1543033
#> e  0.0000000  0.0000000  0.0000000 -0.8944272  0.1825742  0.1543033
#> f  0.0000000  0.0000000  0.0000000  0.0000000 -0.9128709  0.1543033
#> g  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000 -0.9258201

# Formula parts can be inferred when data is not supplied
sbp_basis(list(
  b1 = a ~ b,
  b2 = b1 ~ c
))
#>           b1         b2
#> a  0.7071068  0.4082483
#> b -0.7071068  0.4082483
#> c  0.0000000 -0.8164966

# Non-orthogonal system of balances
sbp_basis(list(
  b1 = a + b + c ~ e + f + g,
  b2 = d ~ a + b + c,
  b3 = d ~ e + g,
  b4 = a ~ e + b,
  b5 = b ~ f,
  b6 = c ~ g
), data = X)
#> Warning: Given basis is not orthogonal
#>           b1         b2         b3         b4         b5         b6
#> a  0.4082483 -0.2886751  0.0000000  0.8164966  0.0000000  0.0000000
#> b  0.4082483 -0.2886751  0.0000000 -0.4082483  0.7071068  0.0000000
#> c  0.4082483 -0.2886751  0.0000000  0.0000000  0.0000000  0.7071068
#> d  0.0000000  0.8660254  0.8164966  0.0000000  0.0000000  0.0000000
#> e -0.4082483  0.0000000 -0.4082483 -0.4082483  0.0000000  0.0000000
#> f -0.4082483  0.0000000  0.0000000  0.0000000 -0.7071068  0.0000000
#> g -0.4082483  0.0000000 -0.4082483  0.0000000  0.0000000 -0.7071068

# Direct construction from a contrast matrix
sbp_basis(cbind(
  c( 1,  1, -1, -1),
  c( 1, -1,  1, -1),
  c( 1, -1, -1,  1)
))
#>      b1   b2   b3
#> c1  0.5  0.5  0.5
#> c2  0.5 -0.5 -0.5
#> c3 -0.5  0.5 -0.5
#> c4 -0.5 -0.5  0.5
```
