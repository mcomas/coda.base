# Tabu search approximation to a sequential binary partition

Builds a sequential binary partition (SBP) by repeatedly applying
grouped tabu search to select balances over the current sets of parts.
At each step, the best candidate split is retained and the remaining
candidate subproblems are explored until an SBP with \\D - 1\\ balances
is obtained.

## Usage

``` r
pb_tabu_search(X, iter = 100, debug = FALSE)
```

## Arguments

- X:

  A numeric matrix with strictly positive finite entries. Rows are
  observations and columns are compositional parts.

- iter:

  Integer. Maximum number of tabu search iterations used in each partial
  search.

- debug:

  Logical. If `TRUE`, progress information is printed during the partial
  searches.

## Value

An integer matrix representing a sequential binary partition. Rows
correspond to the original parts of `X` and columns correspond to
balances. Entries are in \\\\-1,0,1\\\\. The returned matrix has
attribute `"max_steps"`, giving the largest iteration index at which a
best partial solution was found among all partial searches performed.

## Details

This function provides a heuristic approximation to a principal balance
basis. The first balance is searched on the full set of parts, and the
subsequent balances are obtained by recursively refining the best
currently available split.

All partial searches are initialized with the constrained principal
balance of the corresponding grouped composition.

The procedure starts from the trivial grouping where each part forms its
own singleton group. After each partial tabu search, up to three
candidate subproblems may be generated from the selected solution:

- the split between active and inactive groups,

- the left active branch,

- the right active branch.

All generated candidates are stored, and at each stage the candidate
with the largest variance criterion is selected for inclusion in the SBP
and for further refinement.

This is a heuristic search strategy and does not guarantee a globally
optimal SBP.

## See also

[`partial_pb_tabu_search`](https://mcomas.net/coda.base/reference/partial_pb_tabu_search.md),
[`sbp_basis`](https://mcomas.net/coda.base/reference/sbp_basis.md)

## Examples

``` r
set.seed(1)
X <- matrix(rexp(500), ncol = 5)

SBP <- pb_tabu_search(X, iter = 30)
SBP
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    1    1
#> [2,]    1    0    1   -1
#> [3,]   -1    0    0    0
#> [4,]    1   -1   -1    0
#> [5,]    1    1   -1    0
#> attr(,"max_steps")
#> [1] 2
attr(SBP, "max_steps")
#> [1] 2
```
