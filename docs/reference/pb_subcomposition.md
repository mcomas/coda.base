# Recursive constrained principal balances on subcompositions

Recursively construct balances on selected subcompositions, optionally
enforcing groups of variables to remain together through constraints.

## Usage

``` r
pb_subcomposition(
  X,
  variables = seq_len(ncol(X)),
  constraints = NULL,
  angle = FALSE
)
```

## Arguments

- X:

  Compositional data set.

- variables:

  Indices of the variables currently considered.

- constraints:

  Optional list of groups of variables to be constrained together during
  the recursive search.

- angle:

  Logical; if \`TRUE\`, use the angle criterion instead of the variance
  criterion when computing constrained balances.

## Value

A list of balance vectors.
