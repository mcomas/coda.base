# Isometric and orthonormal log-ratio bases

Construct an isometric log-ratio (ilr) basis for a composition with
\\D\\ parts. The ilr basis is an orthonormal basis of the clr-plane and
provides \\D - 1\\ coordinates. The same basis is sometimes referred to
as an orthonormal log-ratio (olr) basis.

## Usage

``` r
ilr_basis(dim, type = "default")

olr_basis(dim, type = "default")
```

## Arguments

- dim:

  Number of parts. It can be:

  - a single integer,

  - a matrix or data frame, in which case the number of columns is used,

  - a character vector of part names, in which case its length is used.

- type:

  Type of ilr basis to construct. Available options are:

  - \`"default"\`: standard Helmert-type ilr basis,

  - \`"pivot"\`: pivot balance basis,

  - \`"cdp"\`: CoDaPack default basis.

## Value

A matrix with \\D\\ rows and \\D - 1\\ columns representing an
orthonormal log-ratio basis.

## Details

For \`type = "default"\`, the function returns the standard Helmert-type
ilr basis. Alternative constructions are available through \`type =
"pivot"\` and \`type = "cdp"\`.

The default basis vectors are: \$\$ h_i = \sqrt{\frac{i}{i+1}} \log
\frac{\sqrt\[i\]{\prod\_{j=1}^i x_j}}{x\_{i+1}}, \qquad i = 1, \ldots,
D - 1 \$\$

## References

Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
Barceló-Vidal, C. (2003). *Isometric logratio transformations for
compositional data analysis*. Mathematical Geology, **35**(3), 279–300.

## Examples

``` r
ilr_basis(5)
#>          ilr1       ilr2       ilr3       ilr4
#> c1  0.7071068  0.4082483  0.2886751  0.2236068
#> c2 -0.7071068  0.4082483  0.2886751  0.2236068
#> c3  0.0000000 -0.8164966  0.2886751  0.2236068
#> c4  0.0000000  0.0000000 -0.8660254  0.2236068
#> c5  0.0000000  0.0000000  0.0000000 -0.8944272
ilr_basis(alimentation[, 1:9])
#>          ilr1       ilr2       ilr3       ilr4       ilr5       ilr6
#> RM  0.7071068  0.4082483  0.2886751  0.2236068  0.1825742  0.1543033
#> WM -0.7071068  0.4082483  0.2886751  0.2236068  0.1825742  0.1543033
#> E   0.0000000 -0.8164966  0.2886751  0.2236068  0.1825742  0.1543033
#> M   0.0000000  0.0000000 -0.8660254  0.2236068  0.1825742  0.1543033
#> F   0.0000000  0.0000000  0.0000000 -0.8944272  0.1825742  0.1543033
#> C   0.0000000  0.0000000  0.0000000  0.0000000 -0.9128709  0.1543033
#> S   0.0000000  0.0000000  0.0000000  0.0000000  0.0000000 -0.9258201
#> N   0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> FV  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#>          ilr7       ilr8
#> RM  0.1336306  0.1178511
#> WM  0.1336306  0.1178511
#> E   0.1336306  0.1178511
#> M   0.1336306  0.1178511
#> F   0.1336306  0.1178511
#> C   0.1336306  0.1178511
#> S   0.1336306  0.1178511
#> N  -0.9354143  0.1178511
#> FV  0.0000000 -0.9428090
ilr_basis(c("a", "b", "c", "d"), type = "pivot")
#>         ilr1       ilr2       ilr3
#> a  0.8660254  0.0000000  0.0000000
#> b -0.2886751  0.8164966  0.0000000
#> c -0.2886751 -0.4082483  0.7071068
#> d -0.2886751 -0.4082483 -0.7071068
```
