# Basis from a sequential binary partition

Construct a balance basis from a sequential binary partition (SBP) or
from a more general collection of balances.

## Usage

``` r
sbp_basis(sbp, data = NULL, fill = FALSE, silent = FALSE)
```

## Arguments

- sbp:

  A non-empty list of formulas or a numeric matrix with finite entries
  describing balances. At least two parts are required. A matrix with
  zero columns can be supplied to construct a full basis with \`fill =
  TRUE\`.

- data:

  Optional compositional data set with column names, or a character
  vector of part names, used when \`sbp\` is a list of formulas. If
  \`data = NULL\`, part names are inferred from the formulas themselves.

- fill:

  Logical; if \`TRUE\`, complete the supplied balances to span the
  clr-plane, using hierarchical balances when compatible and an
  orthogonal complement obtained by QR otherwise.

- silent:

  Logical; if \`FALSE\`, warn when the result is not a basis, or when
  the basis is not orthogonal or orthonormal. Does not suppress errors
  for invalid input.

## Value

A numeric matrix whose columns are normalized balances and, when needed
for completion, general log-contrasts. Unnamed supplied balances receive
names \`b1\`, \`b2\`, etc. Added hierarchical balances have names \`.b\`
followed by their column number; QR complement columns use \`.q\`
instead. Names are made unique without changing supplied balance names.
With \`fill = TRUE\`, the matrix has \`D - 1\` columns for independent
inputs, or more columns if the supplied balances are linearly dependent.

## Details

The argument \`sbp\` can be specified in two ways:

- as a list of formulas, where each formula defines the numerator and
  the denominator groups of a balance,

- as a matrix with one column per balance and one row per part. Positive
  entries indicate parts in the numerator, negative entries indicate
  parts in the denominator, and zeros indicate unused parts.

Formula groups contain part names joined by \`+\`, optionally enclosed
in parentheses. Repeated parts within a group are counted once;
numerator and denominator groups must be non-empty and disjoint. A named
balance can be referenced in another formula to include all its parts,
regardless of their signs. Forward references are allowed, but circular
references and ambiguous names shared by balances and parts are
rejected.

With \`fill = TRUE\`, supplied balances are preserved in their original
order. If they are compatible with a hierarchical SBP, the hierarchy is
completed by adding balances. Compatibility means that each pair
involves disjoint sets of parts, or one balance lies wholly within one
group of the other. Otherwise, QR decomposition supplies an orthonormal
basis for the orthogonal complement of the supplied balances within the
clr-plane. These added columns sum to zero and are orthogonal to all
supplied columns, but need not be balances.

Completion adds \`D - 1 - r\` columns, where \`D\` is the number of
parts and \`r\` is the rank of the supplied balances. The result is
orthonormal if the supplied balances are mutually orthogonal.
Non-orthogonal inputs are kept and produce a warning unless \`silent =
TRUE\`. Linearly dependent inputs are also kept: completion then returns
a redundant generating system rather than a basis, and warns
accordingly. Without completion, an incomplete collection also produces
a warning unless \`silent = TRUE\`.

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

# Complete a partial hierarchy with further balances
B_sbp <- sbp_basis(list(ab = a ~ b), data = letters[1:4], fill = TRUE)
B_sbp
#>           ab        .b2        .b3
#> a  0.7071068  0.4082483  0.2886751
#> b -0.7071068  0.4082483  0.2886751
#> c  0.0000000 -0.8164966  0.2886751
#> d  0.0000000  0.0000000 -0.8660254
crossprod(B_sbp)
#>                ab          .b2           .b3
#> ab   1.000000e+00 1.844191e-17 -1.695163e-18
#> .b2  1.844191e-17 1.000000e+00  1.076573e-17
#> .b3 -1.695163e-18 1.076573e-17  1.000000e+00

# Orthogonal balances with crossing groups: complete by QR
# The added .q5 column is a general log-contrast, not a balance
B_qr <- sbp_basis(list(
  ef = e ~ f,
  cd = c ~ d,
  aef = a + e + f ~ b + c + d,
  ab = a + b ~ c + d + e + f
), data = letters[1:6], fill = TRUE)
round(B_qr, 3)
#>       ef     cd    aef     ab    .q5
#> a  0.000  0.000  0.408  0.577 -0.577
#> b  0.000  0.000 -0.408  0.577  0.577
#> c  0.000  0.707 -0.408 -0.289 -0.289
#> d  0.000 -0.707 -0.408 -0.289 -0.289
#> e  0.707  0.000  0.408 -0.289  0.289
#> f -0.707  0.000  0.408 -0.289  0.289
crossprod(B_qr)
#>               ef            cd           aef           ab
#> ef  1.000000e+00  0.000000e+00  1.844191e-17 1.695163e-18
#> cd  0.000000e+00  1.000000e+00 -1.844191e-17 1.695163e-18
#> aef 1.844191e-17 -1.844191e-17  1.000000e+00 5.382864e-18
#> ab  1.695163e-18  1.695163e-18  5.382864e-18 1.000000e+00
#> .q5 1.045650e-16  1.695163e-18 -6.458532e-17 1.619798e-16
#>               .q5
#> ef   1.045650e-16
#> cd   1.695163e-18
#> aef -6.458532e-17
#> ab   1.619798e-16
#> .q5  1.000000e+00

# Non-orthogonal inputs are preserved; completion warns
B_nonorthogonal <- sbp_basis(list(ab = a ~ b, bc = b ~ c),
                            data = letters[1:4], fill = TRUE)
#> Warning: Given basis is not orthogonal
crossprod(B_nonorthogonal)
#>                ab            bc          .q3
#> ab   1.000000e+00 -5.000000e-01 1.695163e-18
#> bc  -5.000000e-01  1.000000e+00 1.695163e-18
#> .q3  1.695163e-18  1.695163e-18 1.000000e+00

# Dependent inputs are also preserved, producing a redundant system
B_redundant <- sbp_basis(list(ab = a ~ b, ab_again = a ~ b),
                        data = letters[1:4], fill = TRUE)
#> Warning: Given partition is not a basis
dim(B_redundant)
#> [1] 4 4

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
