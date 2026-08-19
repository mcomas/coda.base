# Constrained search for a partial principal balance on grouped parts

Builds a single grouped constrained principal balance from the first
principal component of the grouped composition.

## Usage

``` r
partial_pb_constrained(X, lI = NULL, constrained.criterion = "variance")
```

## Arguments

- X:

  A numeric matrix with strictly positive finite entries. Rows are
  observations and columns are compositional parts.

- lI:

  A list defining a partition of a subset of the columns of `X`. If
  `NULL`, each column of `X` is used as a singleton group.

- constrained.criterion:

  Criterion used to choose the constrained balance. Either `"variance"`
  (default) or `"angle"`.

## Value

A list with the following elements:

- `dim`:

  Dimension of the grouped problem, equal to `length(lI) - 1`.

- `lI`:

  The input grouping structure.

- `variance`:

  Variance criterion of the selected grouped balance.

- `balance_raw`:

  Integer vector in \\\\-1,0,1\\\\ describing the selected grouped
  split.

- `balance`:

  The corresponding one-column balance basis.

- `constrained.criterion`:

  Criterion used to construct the balance.
