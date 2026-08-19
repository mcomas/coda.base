# Exact search for a partial principal balance on grouped parts

Finds the grouped balance with maximum variance among all assignments
whose number of active groups is between `min_parts` and `max_parts`.

## Usage

``` r
partial_pb_exact(
  X,
  lI = NULL,
  min_parts = 2,
  max_parts = NULL,
  method = "restricted"
)
```

## Arguments

- X:

  A numeric matrix with strictly positive finite entries. Rows are
  observations and columns are compositional parts.

- lI:

  A list defining a partition of a subset of the columns of `X`. If
  `NULL`, each column of `X` is used as a singleton group.

- min_parts:

  Integer. Minimum number of active groups.

- max_parts:

  Integer or `NULL`. Maximum number of active groups. If `NULL`, all
  groups may be active.

- method:

  Exhaustive search method. Currently only `"restricted"` is
  implemented; it enumerates only supports whose sizes are inside the
  requested range and assigns signs in binary Gray-code order.

## Value

A list with the following elements:

- `dim`:

  Dimension of the grouped problem, equal to `length(lI) - 1`.

- `lI`:

  The input grouping structure.

- `variance`:

  Variance criterion of the best grouped balance.

- `balance_raw`:

  Integer vector in \\\\-1,0,1\\\\ describing the best grouped split.

- `balance`:

  The corresponding one-column balance basis.

- `min_parts`:

  Minimum number of active groups.

- `max_parts`:

  Maximum number of active groups.

## Details

The search enumerates only supports whose size is between `min_parts`
and `max_parts`. For each support, signs are generated in binary
Gray-code order, fixing the first active group on the left side to avoid
evaluating both a balance and its sign reversal.
